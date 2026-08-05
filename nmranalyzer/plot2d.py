"""Contour viewer for 2D spectra (COSY, HSQC, HMBC, ...).

Opens as its own window so the 1D canvas stays uncluttered.  Contours are
recomputed for whatever region is on screen, at the resolution of the screen,
so zooming in reveals detail instead of magnifying a coarse grid.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import messagebox, ttk

from . import contour
from . import nmr2d

PROJ = 42          # width/height of the projection strips
AXIS_W = 48        # room for the F1 tick labels
AXIS_H = 30        # room for the F2 tick labels
PAD = 10

POSITIVE_COLOURS = ["#1f4e9c", "#2f6fc0", "#4a8bd4", "#6aa5e0"]
NEGATIVE_COLOURS = ["#a52020", "#c14545", "#d46a6a"]


class Plot2DWindow(tk.Toplevel):
    def __init__(self, master, spec):
        super().__init__(master)
        self.spec = spec
        self.title("%s  -  2D" % spec.name)
        self.geometry("1000x760")
        self.minsize(680, 520)

        f2_hi, f2_lo = spec.f2.limits
        f1_hi, f1_lo = spec.f1.limits
        self.view = [f2_hi, f2_lo, f1_hi, f1_lo]     # f2 hi, f2 lo, f1 hi, f1 lo

        self.base_pct = tk.DoubleVar(value=2.0)      # % of the strongest point
        self.n_levels = tk.IntVar(value=10)
        self.factor = tk.DoubleVar(value=1.5)
        self.show_negative = tk.BooleanVar(value=True)
        self.show_diagonal = tk.BooleanVar(value=spec.is_homonuclear())
        self.show_projections = tk.BooleanVar(value=True)
        self.sensitivity = tk.DoubleVar(value=12.0)
        self.status = tk.StringVar(value="")

        self._drag = None
        self._cache_key = None
        self._cache = None

        self._build()
        self.after(60, self.redraw)

    # -- construction -------------------------------------------------------

    def _build(self):
        bar = ttk.Frame(self, padding=(6, 4))
        bar.pack(side="top", fill="x")

        ttk.Label(bar, text="Base level %").pack(side="left")
        ttk.Spinbox(bar, from_=0.05, to=50.0, increment=0.25, width=6,
                    textvariable=self.base_pct,
                    command=self.redraw).pack(side="left", padx=(2, 8))
        ttk.Label(bar, text="Levels").pack(side="left")
        ttk.Spinbox(bar, from_=1, to=30, width=4, textvariable=self.n_levels,
                    command=self.redraw).pack(side="left", padx=(2, 8))
        ttk.Label(bar, text="x").pack(side="left")
        ttk.Spinbox(bar, from_=1.05, to=4.0, increment=0.05, width=5,
                    textvariable=self.factor,
                    command=self.redraw).pack(side="left", padx=(2, 8))

        ttk.Button(bar, text="Lower", width=6,
                   command=lambda: self._scale_base(1 / 1.6)).pack(side="left")
        ttk.Button(bar, text="Raise", width=6,
                   command=lambda: self._scale_base(1.6)).pack(side="left", padx=(2, 8))

        ttk.Checkbutton(bar, text="Negative", variable=self.show_negative,
                        command=self.redraw).pack(side="left")
        ttk.Checkbutton(bar, text="Diagonal", variable=self.show_diagonal,
                        command=self.redraw).pack(side="left")
        ttk.Checkbutton(bar, text="Projections",
                        variable=self.show_projections,
                        command=self.redraw).pack(side="left")

        ttk.Separator(bar, orient="vertical").pack(side="left", fill="y", padx=8)
        ttk.Label(bar, text="Sensitivity").pack(side="left")
        ttk.Spinbox(bar, from_=2.0, to=500.0, increment=2.0, width=6,
                    textvariable=self.sensitivity,
                    command=self.pick_peaks).pack(side="left", padx=(2, 6))
        ttk.Button(bar, text="Pick Peaks",
                   command=self.pick_peaks).pack(side="left", padx=2)
        ttk.Button(bar, text="Reset",
                   command=self.reset_view).pack(side="left", padx=2)

        panes = ttk.PanedWindow(self, orient="vertical")
        panes.pack(side="top", fill="both", expand=True)

        holder = ttk.Frame(panes)
        panes.add(holder, weight=4)
        self.canvas = tk.Canvas(holder, background="white", highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)

        table = ttk.Frame(panes)
        panes.add(table, weight=1)
        cols = ("f2", "f1", "intensity")
        self.tree = ttk.Treeview(table, columns=cols, show="headings", height=6)
        for key, title, width in (("f2", "F2 (ppm)", 130),
                                  ("f1", "F1 (ppm)", 130),
                                  ("intensity", "Intensity", 160)):
            self.tree.heading(key, text=title)
            self.tree.column(key, width=width, anchor="center")
        self.tree.pack(side="left", fill="both", expand=True)
        ttk.Scrollbar(table, orient="vertical", command=self.tree.yview
                      ).pack(side="left", fill="y")
        self.tree.bind("<<TreeviewSelect>>", self._peak_selected)

        ttk.Label(self, textvariable=self.status, padding=(8, 2)).pack(
            side="bottom", fill="x")

        self.canvas.bind("<Configure>", lambda e: self.redraw())
        self.canvas.bind("<ButtonPress-1>", self._press)
        self.canvas.bind("<B1-Motion>", self._drag_move)
        self.canvas.bind("<ButtonRelease-1>", self._release)
        self.canvas.bind("<Double-Button-1>", lambda e: self.reset_view())
        self.canvas.bind("<ButtonPress-3>", self._pan_press)
        self.canvas.bind("<B3-Motion>", self._pan_move)
        self.canvas.bind("<Motion>", self._motion)
        self.canvas.bind("<Leave>", lambda e: self.canvas.delete("cursor"))
        self.canvas.bind("<MouseWheel>", self._wheel)
        self.canvas.bind("<Button-4>", lambda e: self._wheel(e, 1))
        self.canvas.bind("<Button-5>", lambda e: self._wheel(e, -1))

    # -- geometry -----------------------------------------------------------

    @property
    def plot_left(self):
        return AXIS_W + (PROJ if self.show_projections.get() else 0)

    @property
    def plot_top(self):
        return PAD + (PROJ if self.show_projections.get() else 0)

    @property
    def plot_right(self):
        return max(self.plot_left + 20, self.canvas.winfo_width() - PAD)

    @property
    def plot_bottom(self):
        return max(self.plot_top + 20, self.canvas.winfo_height() - AXIS_H)

    def x_of(self, ppm):
        hi, lo = self.view[0], self.view[1]
        if hi == lo:
            return self.plot_left
        return self.plot_left + (hi - ppm) / (hi - lo) * (self.plot_right - self.plot_left)

    def y_of(self, ppm):
        hi, lo = self.view[2], self.view[3]
        if hi == lo:
            return self.plot_top
        return self.plot_top + (hi - ppm) / (hi - lo) * (self.plot_bottom - self.plot_top)

    def ppm_x(self, x):
        hi, lo = self.view[0], self.view[1]
        width = self.plot_right - self.plot_left
        return hi - (x - self.plot_left) / width * (hi - lo) if width else hi

    def ppm_y(self, y):
        hi, lo = self.view[2], self.view[3]
        height = self.plot_bottom - self.plot_top
        return hi - (y - self.plot_top) / height * (hi - lo) if height else hi

    # -- drawing ------------------------------------------------------------

    def _scale_base(self, factor):
        self.base_pct.set(max(0.01, min(90.0, self.base_pct.get() * factor)))
        self.redraw()

    def reset_view(self):
        f2_hi, f2_lo = self.spec.f2.limits
        f1_hi, f1_lo = self.spec.f1.limits
        self.view = [f2_hi, f2_lo, f1_hi, f1_lo]
        self.redraw()

    def _visible_block(self):
        """Rows/cols of the matrix currently on screen."""
        spec = self.spec
        c0 = spec.f2.clamp(spec.f2.index(self.view[0]))
        c1 = spec.f2.clamp(spec.f2.index(self.view[1]))
        r0 = spec.f1.clamp(spec.f1.index(self.view[2]))
        r1 = spec.f1.clamp(spec.f1.index(self.view[3]))
        if c0 > c1:
            c0, c1 = c1, c0
        if r0 > r1:
            r0, r1 = r1, r0
        return r0, r1, c0, c1

    def _contours(self):
        """Contour segments for the visible region, cached per view and level."""
        r0, r1, c0, c1 = self._visible_block()
        width = max(1, int(self.plot_right - self.plot_left))
        height = max(1, int(self.plot_bottom - self.plot_top))
        key = (r0, r1, c0, c1, width, height, self.base_pct.get(),
               self.n_levels.get(), self.factor.get(), self.show_negative.get())
        if key == self._cache_key:
            return self._cache

        block = [row[c0:c1 + 1] for row in self.spec.data[r0:r1 + 1]]
        rows, cols = len(block), len(block[0]) if block else 0
        if rows < 2 or cols < 2:
            self._cache_key, self._cache = key, ([], 1, 1, r0, c0)
            return self._cache

        grid, row_step, col_step = contour.downsample(block, rows, cols,
                                                      height, width)
        base = self.spec.max_intensity() * self.base_pct.get() / 100.0
        positive = contour.levels(base, self.n_levels.get(), self.factor.get())
        wanted = list(positive)
        if self.show_negative.get():
            wanted += contour.levels(base, self.n_levels.get(),
                                     self.factor.get(), negative=True)
        segments = contour.segments(grid, wanted)

        ordered = []
        for i, level in enumerate(positive):
            ordered.append((level, segments.get(level, []),
                            POSITIVE_COLOURS[min(i, len(POSITIVE_COLOURS) - 1)]))
        if self.show_negative.get():
            for i, level in enumerate(contour.levels(base, self.n_levels.get(),
                                                     self.factor.get(),
                                                     negative=True)):
                ordered.append((level, segments.get(level, []),
                                NEGATIVE_COLOURS[min(i, len(NEGATIVE_COLOURS) - 1)]))

        self._cache_key = key
        self._cache = (ordered, row_step, col_step, r0, c0)
        return self._cache

    def redraw(self):
        canvas = self.canvas
        canvas.delete("all")
        if self.spec.rows < 2 or self.spec.cols < 2:
            return

        self._draw_axes()
        ordered, row_step, col_step, r0, c0 = self._contours()

        spec = self.spec
        total = 0
        for _level, segs, colour in ordered:
            if not segs:
                continue
            coords = []
            for (gc0, gr0), (gc1, gr1) in segs:
                x0 = self.x_of(spec.f2.ppm(c0 + (gc0 + 0.5) * col_step - 0.5))
                y0 = self.y_of(spec.f1.ppm(r0 + (gr0 + 0.5) * row_step - 0.5))
                x1 = self.x_of(spec.f2.ppm(c0 + (gc1 + 0.5) * col_step - 0.5))
                y1 = self.y_of(spec.f1.ppm(r0 + (gr1 + 0.5) * row_step - 0.5))
                coords.append((x0, y0, x1, y1))
            total += len(coords)
            for x0, y0, x1, y1 in coords:
                canvas.create_line(x0, y0, x1, y1, fill=colour)

        if self.show_diagonal.get() and spec.is_homonuclear():
            lo = max(min(self.view[1], self.view[3]), -1e6)
            hi = min(max(self.view[0], self.view[2]), 1e6)
            canvas.create_line(self.x_of(hi), self.y_of(hi),
                               self.x_of(lo), self.y_of(lo),
                               fill="#bbbbbb", dash=(4, 4))

        if self.show_projections.get():
            self._draw_projections()
        self._draw_peaks()
        self.status.set("%d contour segments   |   base level %.2f%% of maximum"
                        % (total, self.base_pct.get()))

    def _draw_projections(self):
        spec = self.spec
        r0, r1, c0, c1 = self._visible_block()
        canvas = self.canvas

        top = [max(spec.data[r][c] for r in range(r0, r1 + 1))
               for c in range(c0, c1 + 1)]
        peak = max(top) or 1.0
        coords = []
        for i, v in enumerate(top):
            x = self.x_of(spec.f2.ppm(c0 + i))
            y = self.plot_top - 4 - (v / peak) * (PROJ - 8)
            coords.extend((x, max(PAD, y)))
        if len(coords) >= 4:
            canvas.create_line(*coords, fill="#666666")

        left = [max(row[c0:c1 + 1]) for row in spec.data[r0:r1 + 1]]
        peak = max(left) or 1.0
        coords = []
        for i, v in enumerate(left):
            y = self.y_of(spec.f1.ppm(r0 + i))
            x = self.plot_left - 4 - (v / peak) * (PROJ - 8)
            coords.extend((max(AXIS_W, x), y))
        if len(coords) >= 4:
            canvas.create_line(*coords, fill="#666666")

    def _draw_peaks(self):
        spec = self.spec
        for peak in spec.peaks:
            if not (self.view[1] <= peak.f2_ppm <= self.view[0]):
                continue
            if not (self.view[3] <= peak.f1_ppm <= self.view[2]):
                continue
            x, y = self.x_of(peak.f2_ppm), self.y_of(peak.f1_ppm)
            self.canvas.create_line(x - 5, y, x + 5, y, fill="#cc3333")
            self.canvas.create_line(x, y - 5, x, y + 5, fill="#cc3333")

    def _draw_axes(self):
        canvas = self.canvas
        left, right = self.plot_left, self.plot_right
        top, bottom = self.plot_top, self.plot_bottom
        canvas.create_rectangle(left, top, right, bottom, outline="#333333")

        step = _nice_step(self.view[0] - self.view[1])
        value = int(self.view[1] / step) * step
        while value <= self.view[0] + step:
            if self.view[1] <= value <= self.view[0]:
                x = self.x_of(value)
                canvas.create_line(x, bottom, x, bottom + 5, fill="#333333")
                canvas.create_text(x, bottom + 7, text="%g" % round(value, 4),
                                   anchor="n", font=("TkDefaultFont", 9))
            value += step
        canvas.create_text(right, bottom + 18,
                           text="F2  %s (ppm)" % (self.spec.f2.nucleus or ""),
                           anchor="e", fill="#666666", font=("TkDefaultFont", 9))

        step = _nice_step(self.view[2] - self.view[3])
        value = int(self.view[3] / step) * step
        while value <= self.view[2] + step:
            if self.view[3] <= value <= self.view[2]:
                y = self.y_of(value)
                canvas.create_line(left - 5, y, left, y, fill="#333333")
                canvas.create_text(left - 7, y, text="%g" % round(value, 4),
                                   anchor="e", font=("TkDefaultFont", 9))
            value += step
        canvas.create_text(4, top - 6,
                           text="F1  %s (ppm)" % (self.spec.f1.nucleus or ""),
                           anchor="w", fill="#666666", font=("TkDefaultFont", 9))

    # -- interaction --------------------------------------------------------

    def _press(self, event):
        self._drag = (event.x, event.y)

    def _drag_move(self, event):
        if not self._drag:
            return
        self.canvas.delete("rubber")
        self.canvas.delete("cursor")
        self.canvas.create_rectangle(self._drag[0], self._drag[1],
                                     event.x, event.y, outline="#3366cc",
                                     dash=(4, 3), tags="rubber")

    def _release(self, event):
        if not self._drag:
            return
        x0, y0 = self._drag
        self._drag = None
        self.canvas.delete("rubber")
        if abs(event.x - x0) < 5 or abs(event.y - y0) < 5:
            return
        f2 = sorted((self.ppm_x(x0), self.ppm_x(event.x)), reverse=True)
        f1 = sorted((self.ppm_y(y0), self.ppm_y(event.y)), reverse=True)
        self.view = [f2[0], f2[1], f1[0], f1[1]]
        self.redraw()

    def _pan_press(self, event):
        self._pan = (event.x, event.y, list(self.view))

    def _pan_move(self, event):
        anchor = getattr(self, "_pan", None)
        if not anchor:
            return
        x0, y0, view = anchor
        width = self.plot_right - self.plot_left
        height = self.plot_bottom - self.plot_top
        if width <= 0 or height <= 0:
            return
        dx = (event.x - x0) / width * (view[0] - view[1])
        dy = (event.y - y0) / height * (view[2] - view[3])
        self.view = [view[0] + dx, view[1] + dx, view[2] + dy, view[3] + dy]
        self.redraw()

    def _wheel(self, event, direction=None):
        if direction is None:
            direction = 1 if event.delta > 0 else -1
        factor = 0.8 if direction > 0 else 1.25
        cx, cy = self.ppm_x(event.x), self.ppm_y(event.y)
        self.view = [cx + (self.view[0] - cx) * factor,
                     cx + (self.view[1] - cx) * factor,
                     cy + (self.view[2] - cy) * factor,
                     cy + (self.view[3] - cy) * factor]
        self.redraw()

    def _motion(self, event):
        if event.x < self.plot_left or event.x > self.plot_right:
            return
        if event.y < self.plot_top or event.y > self.plot_bottom:
            return
        f2, f1 = self.ppm_x(event.x), self.ppm_y(event.y)
        self.canvas.delete("cursor")
        self.canvas.create_line(self.plot_left, event.y, self.plot_right, event.y,
                                fill="#aaaacc", dash=(2, 4), tags="cursor")
        self.canvas.create_line(event.x, self.plot_top, event.x, self.plot_bottom,
                                fill="#aaaacc", dash=(2, 4), tags="cursor")
        self.status.set("F2 %.4f ppm    F1 %.4f ppm" % (f2, f1))

    # -- analysis -----------------------------------------------------------

    def pick_peaks(self):
        spec = self.spec
        window = (self.view[1], self.view[0], self.view[3], self.view[2])
        skip = 0.05 * (spec.f2.sw_ppm) if spec.is_homonuclear() else 0.0
        self.config(cursor="watch")
        self.update_idletasks()
        try:
            spec.peaks = nmr2d.pick_peaks_2d(
                spec, sensitivity=max(2.0, float(self.sensitivity.get())),
                window=window, skip_diagonal_ppm=skip)
        finally:
            self.config(cursor="")

        self.tree.delete(*self.tree.get_children())
        for i, peak in enumerate(spec.peaks):
            self.tree.insert("", "end", iid=str(i),
                             values=("%.4f" % peak.f2_ppm,
                                     "%.4f" % peak.f1_ppm,
                                     "%.4g" % peak.intensity))
        self.redraw()
        note = " (diagonal excluded)" if skip else ""
        hint = ("   raise Sensitivity to drop the weak ones"
                if len(spec.peaks) > 30 else "")
        self.status.set("Found %d cross peaks%s at %g x noise%s"
                        % (len(spec.peaks), note, self.sensitivity.get(), hint))

    def _peak_selected(self, _event):
        sel = self.tree.selection()
        if not sel:
            return
        peak = self.spec.peaks[int(sel[0])]
        self.canvas.delete("cursor")
        x, y = self.x_of(peak.f2_ppm), self.y_of(peak.f1_ppm)
        self.canvas.create_oval(x - 8, y - 8, x + 8, y + 8, outline="#cc3333",
                                width=2, tags="cursor")


def _nice_step(span):
    import math
    if span <= 0:
        return 1.0
    raw = span / 8.0
    mag = 10.0 ** math.floor(math.log10(raw))
    for mult in (1.0, 2.0, 5.0, 10.0):
        if raw <= mult * mag:
            return mult * mag
    return 10.0 * mag
