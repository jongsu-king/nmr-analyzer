#!/usr/bin/env python3
"""NMR Analyzer - a desktop viewer and analyser for 1D NMR spectra.

Reads Bruker experiment folders and zip archives, ACD/Labs ``.esp`` files and
JCAMP-DX files; processes raw FIDs; and supports overlay comparison, peak
picking, multiplet analysis and integration with molar-ratio readout.

Standard library only.  Run with ``python3 nmr_analyzer.py``.
"""

from __future__ import annotations

import csv
import os
import sys
import tkinter as tk
from tkinter import colorchooser, filedialog, messagebox, ttk

import analysis
import dsp
import export
import fitting
import nmr2d
import nmrio
import plot2d
import prefs
import solvents
import structwin

APP_NAME = "NMR Analyzer"
APP_VERSION = "1.1"

IS_MAC = sys.platform == "darwin"
ACCEL_MOD = "Command" if IS_MAC else "Control"
ACCEL_TEXT = "Cmd+" if IS_MAC else "Ctrl+"

DATA_FILETYPES = [
    ("All supported", "*.esp *.zip *.jdx *.dx *.jcamp"),
    ("ACD spectra", "*.esp"),
    ("Zipped Bruker data", "*.zip"),
    ("JCAMP-DX", "*.jdx *.dx *.jcamp"),
    ("All files", "*.*"),
]
SESSION_FILETYPES = [("NMR Analyzer session", "*.nmrs"), ("JSON", "*.json")]

PALETTE = [
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
    "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f",
]

MODE_ZOOM = "zoom"
MODE_INTEGRATE = "integrate"
MODE_PEAK = "peak"
MODE_REFERENCE = "reference"

MARGIN_L = 12
MARGIN_R = 12
MARGIN_T = 14
AXIS_H = 34


# ---------------------------------------------------------------------------
# Plot canvas
# ---------------------------------------------------------------------------


class PlotCanvas(tk.Canvas):
    """Draws the spectra and owns all mouse interaction with the plot."""

    def __init__(self, master, app, **kw):
        super().__init__(master, background="white", highlightthickness=0, **kw)
        self.app = app
        self.view_left = 12.0
        self.view_right = -1.0
        self.y_scale = 1.0
        self.stack = False
        self.normalise_each = True
        self.show_peaks = True
        self.show_integrals = True
        self.show_grid = False
        self.show_cursor = True

        self._drag_start = None
        self._drag_kind = None
        self._marker = None

        self.bind("<Configure>", lambda e: self.redraw())
        self.bind("<ButtonPress-1>", self._on_press)
        self.bind("<B1-Motion>", self._on_drag)
        self.bind("<ButtonRelease-1>", self._on_release)
        self.bind("<Double-Button-1>", lambda e: self.app.reset_zoom())
        self.bind("<ButtonPress-3>", self._on_pan_press)
        self.bind("<B3-Motion>", self._on_pan_drag)
        self.bind("<Motion>", self._on_motion)
        self.bind("<Leave>", self._on_leave)
        self.bind("<MouseWheel>", self._on_wheel)
        self.bind("<Button-4>", lambda e: self._on_wheel(e, 1))
        self.bind("<Button-5>", lambda e: self._on_wheel(e, -1))

    # -- geometry -----------------------------------------------------------

    @property
    def plot_left(self):
        return MARGIN_L

    @property
    def plot_right(self):
        return max(MARGIN_L + 10, self.winfo_width() - MARGIN_R)

    @property
    def plot_top(self):
        return MARGIN_T

    @property
    def plot_bottom(self):
        return max(MARGIN_T + 10, self.winfo_height() - AXIS_H)

    def x_of(self, ppm):
        span = self.view_left - self.view_right
        if span == 0:
            return self.plot_left
        frac = (self.view_left - ppm) / span
        return self.plot_left + frac * (self.plot_right - self.plot_left)

    def ppm_of(self, x):
        width = self.plot_right - self.plot_left
        if width <= 0:
            return self.view_left
        frac = (x - self.plot_left) / width
        return self.view_left - frac * (self.view_left - self.view_right)

    # -- view control -------------------------------------------------------

    def data_limits(self):
        """(highest ppm, lowest ppm) across the visible spectra."""
        specs = [s for s in self.app.spectra if s.visible]
        if not specs:
            return None
        return (max(s.limits[0] for s in specs),
                min(s.limits[1] for s in specs))

    def set_view(self, left, right):
        if left < right:
            left, right = right, left
        if left - right < 1e-4:
            return

        # Keep the window on the data: without this, a few scroll-wheel
        # zoom-outs walk the view thousands of ppm away and the plot goes blank.
        limits = self.data_limits()
        if limits:
            hi, lo = limits
            margin = max((hi - lo) * 0.02, 1e-3)
            width = left - right
            if width > (hi - lo) + 2 * margin:
                left, right = hi + margin, lo - margin
            elif left > hi + margin:
                left = hi + margin
                right = left - width
            elif right < lo - margin:
                right = lo - margin
                left = right + width

        self.view_left, self.view_right = left, right
        self.redraw()

    def zoom_at(self, ppm, factor):
        left = ppm + (self.view_left - ppm) * factor
        right = ppm + (self.view_right - ppm) * factor
        self.set_view(left, right)

    # -- drawing ------------------------------------------------------------

    def redraw(self):
        self.delete("all")
        specs = [s for s in self.app.spectra if s.visible]
        if not specs:
            self.create_text(self.winfo_width() // 2, self.winfo_height() // 2,
                             text="Open a Bruker folder/zip, an ACD .esp or a "
                                  "JCAMP-DX file to begin",
                             fill="#999999", font=("TkDefaultFont", 12))
            return

        self._draw_axis()
        active = self.app.active_spectrum()

        height = self.plot_bottom - self.plot_top
        lanes = len(specs) if self.stack else 1
        lane_h = height / lanes
        global_max = self._reference_max(specs)

        for i, spec in enumerate(specs):
            base_y = self.plot_bottom - (i * lane_h if self.stack else 0)
            usable = lane_h * 0.92 if self.stack else height * 0.95
            ref = self._trace_max(spec) if self.normalise_each else global_max
            if ref <= 0:
                ref = 1.0
            scale = usable * self.y_scale * spec.scale / ref
            width = 2 if spec is active else 1
            self._draw_trace(spec, base_y, scale, width)
            if self.show_integrals:
                self._draw_regions(spec, base_y, scale)
            if self.show_peaks and spec is active:
                self._draw_peaks(spec, base_y, scale)
            if self.stack:
                self.create_text(self.plot_right - 6, base_y - usable - 2,
                                 text=spec.name, anchor="ne",
                                 fill=spec.color, font=("TkDefaultFont", 9))

        if not self.stack and len(specs) > 1:
            self._draw_legend(specs)

        if self._marker is not None:
            x = self.x_of(self._marker)
            self.create_line(x, self.plot_top, x, self.plot_bottom,
                             fill="#cc0000", dash=(3, 3))

    def _visible_slice(self, spec):
        a = spec.clamp(spec.index(self.view_left))
        b = spec.clamp(spec.index(self.view_right))
        if a > b:
            a, b = b, a
        return a, b

    def _trace_max(self, spec):
        """Tallest signal in view, ignoring one- or two-point spikes.

        A real line spans many points, so taking the maximum of a short moving
        average stops a single glitch from flattening the whole trace.
        """
        a, b = self._visible_slice(spec)
        data = spec.real
        if b - a < 4:
            chunk = data[a:b + 1]
            return max(chunk) if chunk else 1.0
        best = 0.0
        for i in range(a + 1, b):
            smoothed = (data[i - 1] + data[i] + data[i + 1]) / 3.0
            if smoothed > best:
                best = smoothed
        return best or 1.0

    def _reference_max(self, specs):
        return max((self._trace_max(s) for s in specs), default=1.0)

    def _draw_trace(self, spec, base_y, scale, width):
        a, b = self._visible_slice(spec)
        if b <= a:
            return
        px_left, px_right = self.plot_left, self.plot_right
        columns = max(1, int(px_right - px_left))
        count = b - a + 1
        data = spec.real
        top = self.plot_top

        coords = []
        if count <= columns * 2:
            # Enough screen space for every point.
            for i in range(a, b + 1):
                x = self.x_of(spec.ppm(i))
                y = base_y - data[i] * scale
                coords.append(x)
                coords.append(max(top - 200, min(base_y + 200, y)))
        else:
            # Min/max envelope: one vertical segment per pixel column.
            per = count / columns
            idx = a
            for col in range(columns):
                stop = a + int((col + 1) * per)
                stop = min(stop, b + 1)
                if stop <= idx:
                    stop = idx + 1
                chunk = data[idx:stop]
                if not chunk:
                    break
                lo, hi = min(chunk), max(chunk)
                x = px_left + col
                y_hi = max(top - 200, min(base_y + 200, base_y - hi * scale))
                y_lo = max(top - 200, min(base_y + 200, base_y - lo * scale))
                coords.extend((x, y_hi, x, y_lo))
                idx = stop
        if len(coords) >= 4:
            self.create_line(*coords, fill=spec.color, width=width)

    def _draw_regions(self, spec, base_y, scale):
        if not spec.regions:
            return
        norm = analysis.normalise(spec.regions)
        for region, value in zip(spec.regions, norm):
            x1 = self.x_of(region.hi)
            x2 = self.x_of(region.lo)
            if x2 < self.plot_left or x1 > self.plot_right:
                continue
            shade = self.create_rectangle(x1, self.plot_top, x2, base_y,
                                          fill="#eef4fb", outline="", width=0)
            self.tag_lower(shade)
            self._draw_integral_curve(spec, region, base_y, scale, x1, x2)
            if region.fit is not None:
                self._draw_fit(spec, region, base_y, scale)
            label = "%.2f" % value
            if region.protons:
                label = "%gH" % region.protons
            # Above the shaded band, where it cannot collide with the ppm axis.
            self.create_text((x1 + x2) / 2, self.plot_top + 1, text=label,
                             anchor="n", fill="#2266aa",
                             font=("TkDefaultFont", 9))
            for x in (x1, x2):
                self.create_line(x, base_y - 6, x, base_y + 6, fill="#2266aa")

    def _draw_integral_curve(self, spec, region, base_y, scale, x1, x2):
        a = spec.clamp(spec.index(region.hi))
        b = spec.clamp(spec.index(region.lo))
        if a > b:
            a, b = b, a
        if b - a < 2:
            return
        total = sum(spec.real[a:b + 1])
        if total == 0:
            return
        height = (base_y - self.plot_top) * 0.35
        steps = min(400, b - a)
        stride = max(1, (b - a) // steps)
        run = 0.0
        coords = []
        for i in range(a, b + 1, stride):
            run += sum(spec.real[i:i + stride])
            x = self.x_of(spec.ppm(i))
            y = base_y - height * (run / total) - 4
            coords.extend((x, y))
        if len(coords) >= 4:
            self.create_line(*coords, fill="#2266aa", width=1)

    def _draw_fit(self, spec, region, base_y, scale):
        """Overlay the fitted envelope and its individual components."""
        a = spec.clamp(spec.index(region.hi))
        b = spec.clamp(spec.index(region.lo))
        if a > b:
            a, b = b, a
        if b - a < 3:
            return
        stride = max(1, (b - a) // 300)
        ppm_values = [spec.ppm(i) for i in range(a, b + 1, stride)]
        xs = [self.x_of(p) for p in ppm_values]

        for index in range(len(region.fit.peaks)):
            comp = region.fit.component(index, ppm_values)
            coords = []
            for x, v in zip(xs, comp):
                coords.extend((x, base_y - v * scale))
            if len(coords) >= 4:
                self.create_line(*coords, fill="#88bb55", width=1)

        total = region.fit.evaluate(ppm_values)
        coords = []
        for x, v in zip(xs, total):
            coords.extend((x, base_y - v * scale))
        if len(coords) >= 4:
            self.create_line(*coords, fill="#cc7700", width=1, dash=(4, 2))

    def _draw_peaks(self, spec, base_y, scale):
        shown = 0
        last_label_x = None
        for peak in spec.peaks:
            if not (self.view_right <= peak.ppm <= self.view_left):
                continue
            if shown > 80:
                break
            x = self.x_of(peak.ppm)
            y = base_y - peak.height * scale
            self.create_line(x, y - 4, x, y - 12, fill="#aa3333")
            # Only label a peak when it will not collide with the previous one.
            if last_label_x is None or abs(x - last_label_x) > 34:
                self.create_text(x, y - 14, text="%.3f" % peak.ppm, anchor="s",
                                 fill="#aa3333", font=("TkDefaultFont", 8))
                last_label_x = x
            shown += 1

    def _draw_axis(self):
        y = self.plot_bottom
        self.create_line(self.plot_left, y, self.plot_right, y, fill="#333333")
        span = self.view_left - self.view_right
        step = _nice_step(span)
        minor = step / 5.0

        # Minor ticks first, so the labelled ones draw over them.
        value = int(self.view_right / minor) * minor
        while value <= self.view_left + minor:
            if self.view_right <= value <= self.view_left:
                x = self.x_of(value)
                self.create_line(x, y, x, y + 3, fill="#999999")
            value += minor

        value = int(self.view_right / step) * step
        last_label_x = None
        while value <= self.view_left + step:
            if self.view_right - 1e-9 <= value <= self.view_left + 1e-9:
                x = self.x_of(value)
                if self.show_grid:
                    self.create_line(x, self.plot_top, x, y, fill="#eeeeee")
                self.create_line(x, y, x, y + 6, fill="#333333")
                # The ppm axis runs right to left, so x *decreases* as the
                # value increases: the spacing test needs the absolute gap.
                if x < self.plot_right - 26 and (last_label_x is None
                                                 or abs(x - last_label_x) > 28):
                    self.create_text(x, y + 8, text="%g" % round(value, 4),
                                     anchor="n", font=("TkDefaultFont", 9))
                    last_label_x = x
            value += step
        self.create_text(self.plot_right, y + 8, text="ppm", anchor="ne",
                         font=("TkDefaultFont", 9), fill="#666666")

    def _draw_legend(self, specs):
        """Name each trace in overlay mode; stacked mode labels its own lanes.

        Drawn on an opaque strip because integral and peak labels also live at
        the top of the plot and would otherwise show through.
        """
        y = self.plot_top + 2
        for spec in specs:
            item = self.create_text(self.plot_right - 5, y, text=spec.name,
                                    anchor="ne", fill=spec.color,
                                    font=("TkDefaultFont", 9))
            x1, y1, x2, y2 = self.bbox(item)
            backing = self.create_rectangle(x1 - 3, y1 - 1, x2 + 2, y2 + 1,
                                            fill="white", outline="")
            self.tag_lower(backing, item)
            y += 14

    # -- interaction --------------------------------------------------------

    def _on_press(self, event):
        self._drag_start = event.x
        self._drag_kind = self.app.mode.get()
        if self._drag_kind == MODE_REFERENCE:
            self.app.set_reference(self.ppm_of(event.x))
            self._drag_start = None
        elif self._drag_kind == MODE_PEAK:
            self._drag_start = None
            self.app.pick_peaks_here(self.ppm_of(event.x))

    def _on_drag(self, event):
        if self._drag_start is None:
            return
        self.delete("rubber")
        self.delete("cursor")
        colour = "#cc0000" if self._drag_kind == MODE_INTEGRATE else "#3366cc"
        self.create_rectangle(self._drag_start, self.plot_top, event.x,
                              self.plot_bottom, outline=colour,
                              dash=(4, 3), tags="rubber")
        # Live readout of what the selection covers.
        p1, p2 = self.ppm_of(self._drag_start), self.ppm_of(event.x)
        lo, hi = min(p1, p2), max(p1, p2)
        spec = self.app.active_spectrum()
        text = "%.3f - %.3f ppm  (%.3f ppm" % (hi, lo, hi - lo)
        text += ", %.1f Hz)" % ((hi - lo) * spec.sf) if spec else ")"
        self.create_text((self._drag_start + event.x) / 2, self.plot_top + 2,
                         text=text, anchor="n", fill=colour,
                         font=("TkDefaultFont", 9), tags="rubber")

    def _on_release(self, event):
        if self._drag_start is None:
            return
        self.delete("rubber")
        x1, x2 = self._drag_start, event.x
        self._drag_start = None
        if abs(x2 - x1) < 4:
            return
        p1, p2 = self.ppm_of(x1), self.ppm_of(x2)
        if self._drag_kind == MODE_INTEGRATE:
            self.app.add_region(min(p1, p2), max(p1, p2))
        else:
            self.set_view(max(p1, p2), min(p1, p2))

    def _on_pan_press(self, event):
        self._pan_anchor = (event.x, self.view_left, self.view_right)

    def _on_pan_drag(self, event):
        anchor = getattr(self, "_pan_anchor", None)
        if not anchor:
            return
        x0, left0, right0 = anchor
        width = self.plot_right - self.plot_left
        if width <= 0:
            return
        shift = (event.x - x0) / width * (left0 - right0)
        self.set_view(left0 + shift, right0 + shift)

    def _on_motion(self, event):
        ppm = self.ppm_of(event.x)
        spec = self.app.active_spectrum()
        hz = ppm * spec.sf if spec else 0.0
        self.app.status.set("%.4f ppm    %.2f Hz" % (ppm, hz))
        self._draw_cursor(event.x, ppm, spec)

    def _draw_cursor(self, x, ppm, spec):
        """A crosshair with a live ppm readout, redrawn without a full repaint."""
        self.delete("cursor")
        if not self.show_cursor or x < self.plot_left or x > self.plot_right:
            return
        self.create_line(x, self.plot_top, x, self.plot_bottom,
                         fill="#9999bb", dash=(2, 4), tags="cursor")
        label = "%.3f" % ppm
        anchor = "nw" if x < self.plot_right - 60 else "ne"
        offset = 4 if anchor == "nw" else -4
        self.create_text(x + offset, self.plot_top + 2, text=label,
                         anchor=anchor, fill="#5555aa",
                         font=("TkDefaultFont", 9), tags="cursor")

    def _on_leave(self, _event):
        self.delete("cursor")

    def _on_wheel(self, event, direction=None):
        if direction is None:
            direction = 1 if event.delta > 0 else -1
        state = getattr(event, "state", 0)
        if state & 0x0001 or state & 0x0004:   # Shift or Control: intensity
            self.y_scale *= 1.25 if direction > 0 else 0.8
            self.redraw()
        else:
            self.zoom_at(self.ppm_of(event.x), 0.8 if direction > 0 else 1.25)


def _nice_step(span):
    """Pick a human-friendly axis tick spacing for the given ppm span."""
    import math
    if span <= 0:
        return 1.0
    raw = span / 10.0
    mag = 10.0 ** math.floor(math.log10(raw))
    for mult in (1.0, 2.0, 5.0, 10.0):
        if raw <= mult * mag:
            return mult * mag
    return 10.0 * mag


# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.prefs = prefs.Preferences()

        self.minsize(940, 620)
        self._restore_geometry()

        self.spectra = []
        self.spectra_2d = []
        self.mode = tk.StringVar(value=MODE_ZOOM)
        self.status = tk.StringVar(value="Ready")

        # Shared view options: created before the menus, which bind to them.
        self.var_stack = tk.BooleanVar(value=False)
        self.var_norm = tk.BooleanVar(value=True)
        self.var_show_peaks = tk.BooleanVar(value=True)
        self.var_show_integrals = tk.BooleanVar(value=True)
        self.var_grid = tk.BooleanVar(value=False)

        # Document state
        self.session_path = None        # None until saved under a name
        self.dirty = False
        self._undo_stack = []
        self._redo_stack = []

        self._build_menu()
        self._build_layout()
        self._apply_preferences()
        self.reset_zoom()
        self._update_title()

        self.protocol("WM_DELETE_WINDOW", self.quit_app)
        if IS_MAC:
            # Route the macOS application menu through the same guards.
            self.createcommand("::tk::mac::Quit", self.quit_app)
            self.createcommand("tk::mac::ShowPreferences", self.show_preferences)

    # -- preferences --------------------------------------------------------

    def _restore_geometry(self):
        saved = self.prefs.get("geometry", "")
        if saved:
            try:
                size = saved.split("+")[0]
                w, h = (int(v) for v in size.split("x"))
                # Never restore a window bigger than the current display.
                if w <= self.winfo_screenwidth() and h <= self.winfo_screenheight():
                    self.geometry(saved)
                    return
            except (ValueError, IndexError):
                pass
        self.geometry("1280x800")

    def _apply_preferences(self):
        self.var_sens.set(float(self.prefs.get("sensitivity", 8.0)))
        self.var_lb.set(float(self.prefs.get("line_broadening", 0.3)))
        self.var_si.set(str(self.prefs.get("zero_fill", "65536")))
        self.var_bl_method.set(self.prefs.get("baseline_method", "spline"))
        self.var_stack.set(bool(self.prefs.get("stack", False)))
        self.var_norm.set(bool(self.prefs.get("normalise_each", True)))
        self.var_show_peaks.set(bool(self.prefs.get("show_peaks", True)))
        self.var_show_integrals.set(bool(self.prefs.get("show_integrals", True)))
        self.var_grid.set(bool(self.prefs.get("show_grid", False)))
        self._sync_view_options()
        self._refresh_undo_labels()

    def _store_preferences(self):
        self.prefs.set("geometry", self.winfo_geometry())
        self.prefs.set("sensitivity", float(self.var_sens.get()))
        self.prefs.set("line_broadening", float(self.var_lb.get()))
        self.prefs.set("zero_fill", self.var_si.get())
        self.prefs.set("baseline_method", self.var_bl_method.get())
        self.prefs.set("stack", bool(self.var_stack.get()))
        self.prefs.set("normalise_each", bool(self.var_norm.get()))
        self.prefs.set("show_peaks", bool(self.var_show_peaks.get()))
        self.prefs.set("show_integrals", bool(self.var_show_integrals.get()))
        self.prefs.set("show_grid", bool(self.var_grid.get()))
        self.prefs.save()

    def show_preferences(self):
        messagebox.showinfo(
            APP_NAME,
            "Settings are remembered automatically from the panels you use "
            "and are stored in\n\n%s" % prefs.config_path())

    # -- document state -----------------------------------------------------

    def _update_title(self):
        name = (os.path.basename(self.session_path) if self.session_path
                else "Untitled")
        self.title("%s%s - %s" % (name, " *" if self.dirty else "", APP_NAME))

    def mark_dirty(self, dirty=True):
        if self.dirty != dirty:
            self.dirty = dirty
            self._update_title()

    def confirm_discard(self):
        """Give the user a chance to save. False means "cancel the action"."""
        if not self.dirty or not self.spectra:
            return True
        answer = messagebox.askyesnocancel(
            APP_NAME, "Save changes to %s before continuing?"
            % (os.path.basename(self.session_path) if self.session_path
               else "this session"))
        if answer is None:
            return False
        if answer:
            return self.save_session()
        return True

    # -- undo ---------------------------------------------------------------

    def _capture(self):
        """A cheap snapshot of the analysis state (not the spectral data)."""
        return [(spec,
                 [(r.lo, r.hi, r.protons, r.label) for r in spec.regions],
                 list(spec.peaks),
                 spec.ref_shift)
                for spec in self.spectra]

    def push_undo(self, label):
        self._undo_stack.append((label, self._capture()))
        del self._undo_stack[:-40]
        self._redo_stack.clear()
        self.mark_dirty()
        self._refresh_undo_labels()

    def _restore(self, snapshot):
        for spec, regions, peaks, ref_shift in snapshot:
            spec.ref_shift = ref_shift
            spec.peaks = list(peaks)
            spec.regions = []
            for lo, hi, protons, label in regions:
                region = analysis.Region(lo, hi, label=label)
                region.protons = protons
                analysis.integrate_region(spec, region)
                spec.regions.append(region)
        self.refresh_tables()
        self.plot.redraw()

    def undo(self):
        if not self._undo_stack:
            return
        label, snapshot = self._undo_stack.pop()
        self._redo_stack.append((label, self._capture()))
        self._restore(snapshot)
        self.mark_dirty()
        self._refresh_undo_labels()
        self.status.set("Undid %s" % label)

    def redo(self):
        if not self._redo_stack:
            return
        label, snapshot = self._redo_stack.pop()
        self._undo_stack.append((label, self._capture()))
        self._restore(snapshot)
        self.mark_dirty()
        self._refresh_undo_labels()
        self.status.set("Redid %s" % label)

    def _refresh_undo_labels(self):
        undo_label = ("Undo %s" % self._undo_stack[-1][0]
                      if self._undo_stack else "Undo")
        redo_label = ("Redo %s" % self._redo_stack[-1][0]
                      if self._redo_stack else "Redo")
        self.edit_menu.entryconfig(0, label=undo_label,
                                   state="normal" if self._undo_stack else "disabled")
        self.edit_menu.entryconfig(1, label=redo_label,
                                   state="normal" if self._redo_stack else "disabled")
        if hasattr(self, "btn_undo"):
            self.btn_undo.state(("!disabled",) if self._undo_stack else ("disabled",))
            self.btn_redo.state(("!disabled",) if self._redo_stack else ("disabled",))

    # -- construction -------------------------------------------------------

    def _accel(self, key, shift=False):
        """Menu accelerator text for the current platform."""
        return "%s%s%s" % (ACCEL_TEXT, "Shift+" if shift else "", key.upper())

    def _bind_accel(self, key, handler, shift=False):
        sequence = "<%s-%s%s>" % (ACCEL_MOD, "Shift-" if shift else "",
                                  key.upper() if shift else key.lower())
        self.bind_all(sequence, lambda e: (handler(), "break")[1])

    def _build_menu(self):
        menu = tk.Menu(self)

        # -- File
        file_menu = tk.Menu(menu, tearoff=0)
        file_menu.add_command(label="New", accelerator=self._accel("n"),
                              command=self.new_session)
        file_menu.add_separator()
        file_menu.add_command(label="Open Data Files...",
                              accelerator=self._accel("o"),
                              command=self.open_files)
        file_menu.add_command(label="Open Bruker Folder...",
                              command=self.open_folder)
        file_menu.add_command(label="Open Session...",
                              accelerator=self._accel("o", shift=True),
                              command=self.open_session)
        self.recent_menu = tk.Menu(file_menu, tearoff=0)
        file_menu.add_cascade(label="Open Recent", menu=self.recent_menu)
        file_menu.add_separator()
        file_menu.add_command(label="Save", accelerator=self._accel("s"),
                              command=self.save_session)
        file_menu.add_command(label="Save As...",
                              accelerator=self._accel("s", shift=True),
                              command=self.save_session_as)
        file_menu.add_separator()

        export_menu = tk.Menu(file_menu, tearoff=0)
        export_menu.add_command(label="Plot as SVG...", command=self.export_svg)
        export_menu.add_command(label="Plot as PostScript...",
                                command=self.export_plot)
        export_menu.add_separator()
        export_menu.add_command(label="Peak List (CSV)...",
                                command=self.export_peaks)
        export_menu.add_command(label="Integrals (CSV)...",
                                command=self.export_integrals)
        export_menu.add_command(label="Report (Text)...",
                                command=self.export_report)
        file_menu.add_cascade(label="Export", menu=export_menu)
        file_menu.add_separator()
        file_menu.add_command(label="Close Spectrum", command=self.remove_spectrum)
        file_menu.add_command(label="Quit", accelerator=self._accel("q"),
                              command=self.quit_app)
        menu.add_cascade(label="File", menu=file_menu)

        # -- Edit
        self.edit_menu = tk.Menu(menu, tearoff=0)
        self.edit_menu.add_command(label="Undo", accelerator=self._accel("z"),
                                   command=self.undo, state="disabled")
        self.edit_menu.add_command(label="Redo",
                                   accelerator=self._accel("z", shift=True),
                                   command=self.redo, state="disabled")
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Copy Report",
                                   accelerator=self._accel("c"),
                                   command=self.copy_report)
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Delete Selected Region",
                                   command=self.delete_region)
        self.edit_menu.add_command(label="Clear All Regions",
                                   command=self.clear_regions)
        self.edit_menu.add_command(label="Clear Peak List",
                                   command=self.clear_peaks)
        menu.add_cascade(label="Edit", menu=self.edit_menu)

        # -- View
        view_menu = tk.Menu(menu, tearoff=0)
        view_menu.add_command(label="Reset Zoom", accelerator=self._accel("r"),
                              command=self.reset_zoom)
        view_menu.add_command(label="Zoom Aromatic (6-9 ppm)",
                              command=lambda: self.plot.set_view(9.0, 6.0))
        view_menu.add_command(label="Zoom Aliphatic (0-5 ppm)",
                              command=lambda: self.plot.set_view(5.0, 0.0))
        view_menu.add_separator()
        view_menu.add_command(label="Taller", accelerator=self._accel("+"),
                              command=lambda: self.scale_intensity(1.4))
        view_menu.add_command(label="Shorter", accelerator=self._accel("-"),
                              command=lambda: self.scale_intensity(1 / 1.4))
        view_menu.add_separator()
        view_menu.add_checkbutton(label="Stack Spectra", variable=self.var_stack,
                                  command=self._sync_view_options)
        view_menu.add_checkbutton(label="Normalise Each",
                                  variable=self.var_norm,
                                  command=self._sync_view_options)
        view_menu.add_checkbutton(label="Show Peak Labels",
                                  variable=self.var_show_peaks,
                                  command=self._sync_view_options)
        view_menu.add_checkbutton(label="Show Integrals",
                                  variable=self.var_show_integrals,
                                  command=self._sync_view_options)
        menu.add_cascade(label="View", menu=view_menu)

        # -- Tools
        tools_menu = tk.Menu(menu, tearoff=0)
        tools_menu.add_command(label="Auto Phase", accelerator=self._accel("p"),
                               command=self.auto_phase)
        tools_menu.add_command(label="Reprocess Active Spectrum",
                               command=self.apply_processing)
        tools_menu.add_command(label="Reprocess All Spectra",
                               command=self.reprocess_all)
        tools_menu.add_separator()
        tools_menu.add_command(label="Pick Peaks", accelerator=self._accel("k"),
                               command=self.pick_peaks_all)
        tools_menu.add_command(label="Auto Integrate",
                               accelerator=self._accel("i"),
                               command=self.auto_integrate)
        tools_menu.add_command(label="Fit Lines in Selected Region",
                               accelerator=self._accel("f"),
                               command=self.fit_region)
        tools_menu.add_separator()
        tools_menu.add_command(label="Calibrate to Solvent Peak",
                               accelerator=self._accel("l"),
                               command=self.calibrate_to_solvent)
        tools_menu.add_command(label="Clear Calibration",
                               command=self.clear_reference)
        tools_menu.add_separator()
        self.twod_menu = tk.Menu(tools_menu, tearoff=0)
        tools_menu.add_cascade(label="2D Spectra", menu=self.twod_menu)
        tools_menu.add_command(label="Structure...",
                               accelerator=self._accel("d"),
                               command=self.open_structure)
        tools_menu.add_separator()
        mode_menu = tk.Menu(tools_menu, tearoff=0)
        for label, value, key in (("Zoom", MODE_ZOOM, "1"),
                                  ("Integrate", MODE_INTEGRATE, "2"),
                                  ("Peak", MODE_PEAK, "3"),
                                  ("Reference", MODE_REFERENCE, "4")):
            mode_menu.add_radiobutton(label=label, value=value,
                                      variable=self.mode, accelerator=key)
            self.bind_all("<Key-%s>" % key,
                          lambda e, v=value: self.mode.set(v))
        tools_menu.add_cascade(label="Mouse Mode", menu=mode_menu)
        menu.add_cascade(label="Tools", menu=tools_menu)

        # -- Help
        help_menu = tk.Menu(menu, tearoff=0, name="help")
        help_menu.add_command(label="Keyboard Shortcuts",
                              command=self.show_shortcuts)
        help_menu.add_command(label="About %s" % APP_NAME, command=self.show_about)
        menu.add_cascade(label="Help", menu=help_menu)

        self.config(menu=menu)
        self.refresh_recent_menu()
        self._refresh_2d_menu()

        self._bind_accel("n", self.new_session)
        self._bind_accel("o", self.open_files)
        self._bind_accel("O", self.open_session, shift=True)
        self._bind_accel("s", self.save_session)
        self._bind_accel("S", self.save_session_as, shift=True)
        self._bind_accel("q", self.quit_app)
        self._bind_accel("z", self.undo)
        self._bind_accel("Z", self.redo, shift=True)
        self._bind_accel("c", self.copy_report)
        self._bind_accel("r", self.reset_zoom)
        self._bind_accel("p", self.auto_phase)
        self._bind_accel("k", self.pick_peaks_all)
        self._bind_accel("i", self.auto_integrate)
        self._bind_accel("f", self.fit_region)
        self._bind_accel("l", self.calibrate_to_solvent)
        self._bind_accel("d", self.open_structure)
        self.bind_all("<%s-plus>" % ACCEL_MOD,
                      lambda e: self.scale_intensity(1.4))
        self.bind_all("<%s-equal>" % ACCEL_MOD,
                      lambda e: self.scale_intensity(1.4))
        self.bind_all("<%s-minus>" % ACCEL_MOD,
                      lambda e: self.scale_intensity(1 / 1.4))

    def refresh_recent_menu(self):
        self.recent_menu.delete(0, "end")
        entries = self.prefs.recent
        if not entries:
            self.recent_menu.add_command(label="(nothing yet)", state="disabled")
            return
        for path in entries:
            self.recent_menu.add_command(
                label=os.path.basename(path.rstrip(os.sep)) or path,
                command=lambda p=path: self.open_path(p))
        self.recent_menu.add_separator()
        self.recent_menu.add_command(label="Clear Menu",
                                     command=self.clear_recent)

    def clear_recent(self):
        self.prefs.clear_recent()
        self.prefs.save()
        self.refresh_recent_menu()

    def _build_layout(self):
        self._build_toolbar()

        outer = ttk.PanedWindow(self, orient="horizontal")
        outer.pack(side="top", fill="both", expand=True)

        left = ttk.Frame(outer, padding=4)
        outer.add(left, weight=0)
        self._build_left_panel(left)

        right = ttk.PanedWindow(outer, orient="vertical")
        outer.add(right, weight=1)

        plot_frame = ttk.Frame(right)
        right.add(plot_frame, weight=3)
        self.plot = PlotCanvas(plot_frame, self)
        self.plot.pack(fill="both", expand=True)

        tabs = ttk.Notebook(right)
        right.add(tabs, weight=1)
        self._build_tabs(tabs)

        # Weights alone do not place the initial sash, so set it explicitly.
        self.after(80, lambda: right.sashpos(0, int(right.winfo_height() * 0.62)))

        bar = ttk.Frame(self, padding=(6, 2))
        bar.pack(side="bottom", fill="x")
        ttk.Label(bar, textvariable=self.status).pack(side="left")
        ttk.Label(bar, foreground="#777777",
                  text="drag = zoom/integrate · right-drag = pan · "
                       "wheel = zoom · shift+wheel = intensity · "
                       "double-click = reset").pack(side="right")

    def _build_toolbar(self):
        """Two rows: document and mouse mode on top, analysis and view below."""
        bar = ttk.Frame(self, padding=(6, 3))
        bar.pack(side="top", fill="x")

        def group(parent, first=False):
            if not first:
                ttk.Separator(parent, orient="vertical").pack(
                    side="left", fill="y", padx=7, pady=1)
            frame = ttk.Frame(parent)
            frame.pack(side="left")
            return frame

        def button(parent, text, command, tip):
            widget = ttk.Button(parent, text=text, command=command, width=len(text) + 2)
            widget.pack(side="left", padx=1)
            self._add_tip(widget, tip)
            return widget

        # -- row 1: document + mouse mode
        row = ttk.Frame(bar)
        row.pack(side="top", fill="x")

        g = group(row, first=True)
        button(g, "New", self.new_session, "Start a new session (%sN)" % ACCEL_TEXT)
        button(g, "Open", self.open_files, "Open data files (%sO)" % ACCEL_TEXT)
        button(g, "Save", self.save_session, "Save the session (%sS)" % ACCEL_TEXT)

        g = group(row)
        self.btn_undo = button(g, "Undo", self.undo, "Undo (%sZ)" % ACCEL_TEXT)
        self.btn_redo = button(g, "Redo", self.redo, "Redo (%sShift+Z)" % ACCEL_TEXT)

        g = group(row)
        ttk.Label(g, text="Mouse:").pack(side="left", padx=(0, 4))
        for label, value, key in (("Zoom", MODE_ZOOM, "1"),
                                  ("Integrate", MODE_INTEGRATE, "2"),
                                  ("Peak", MODE_PEAK, "3"),
                                  ("Reference", MODE_REFERENCE, "4")):
            widget = ttk.Radiobutton(g, text=label, value=value,
                                     variable=self.mode)
            widget.pack(side="left", padx=2)
            self._add_tip(widget, "%s mode (key %s)" % (label, key))

        # -- row 2: analysis + view
        row = ttk.Frame(bar)
        row.pack(side="top", fill="x", pady=(3, 0))

        g = group(row, first=True)
        button(g, "Auto Phase", self.auto_phase,
               "Fit PH0/PH1 by minimising negative signal (%sP)" % ACCEL_TEXT)
        button(g, "Baseline", self.toggle_baseline,
               "Toggle baseline correction on the active spectrum")
        button(g, "Calibrate", self.calibrate_to_solvent,
               "Move the residual solvent line to its book value (%sL)" % ACCEL_TEXT)

        g = group(row)
        button(g, "Pick Peaks", self.pick_peaks_all,
               "Pick peaks in the displayed region (%sK)" % ACCEL_TEXT)
        button(g, "Integrate", self.auto_integrate,
               "Group the picked peaks into integration regions (%sI)" % ACCEL_TEXT)
        button(g, "Fit", self.fit_region,
               "Deconvolute the selected region (%sF)" % ACCEL_TEXT)

        g = group(row)
        button(g, "Reset", self.reset_zoom, "Show the whole spectrum (%sR)" % ACCEL_TEXT)
        button(g, "Arom.", lambda: self.plot.set_view(9.0, 6.0),
               "Zoom to the aromatic region, 6-9 ppm")
        button(g, "Aliph.", lambda: self.plot.set_view(5.0, 0.0),
               "Zoom to the aliphatic region, 0-5 ppm")
        button(g, "+", lambda: self.scale_intensity(1.4), "Taller")
        button(g, "-", lambda: self.scale_intensity(1 / 1.4), "Shorter")

        g = group(row)
        for text, var in (("Stack", self.var_stack),
                          ("Normalise", self.var_norm)):
            ttk.Checkbutton(g, text=text, variable=var,
                            command=self._sync_view_options).pack(side="left")

    def _add_tip(self, widget, text):
        """Minimal hover tooltip; ttk has none built in."""
        def enter(_event):
            self._hide_tip()
            x = widget.winfo_rootx() + 12
            y = widget.winfo_rooty() + widget.winfo_height() + 4
            tip = tk.Toplevel(widget)
            tip.wm_overrideredirect(True)
            tip.wm_geometry("+%d+%d" % (x, y))
            tk.Label(tip, text=text, background="#ffffe0", relief="solid",
                     borderwidth=1, font=("TkDefaultFont", 9),
                     padx=5, pady=2).pack()
            self._tip = tip

        widget.bind("<Enter>", enter)
        widget.bind("<Leave>", lambda e: self._hide_tip())
        widget.bind("<ButtonPress>", lambda e: self._hide_tip(), add="+")

    def _hide_tip(self):
        tip = getattr(self, "_tip", None)
        if tip is not None:
            tip.destroy()
            self._tip = None

    def toggle_baseline(self):
        self.var_baseline.set(not self.var_baseline.get())
        self.apply_processing()
        self.status.set("Baseline correction %s"
                        % ("on" if self.var_baseline.get() else "off"))

    def _build_left_panel(self, parent):
        ttk.Label(parent, text="Spectra", font=("TkDefaultFont", 10, "bold")
                  ).pack(anchor="w")
        self.tree = ttk.Treeview(parent, columns=("show",), show="tree headings",
                                 height=8, selectmode="browse")
        self.tree.heading("#0", text="Name")
        self.tree.heading("show", text="On")
        self.tree.column("#0", width=210, stretch=True)
        self.tree.column("show", width=34, anchor="center", stretch=False)
        self.tree.pack(fill="both", expand=False, pady=(2, 4))
        self.tree.bind("<<TreeviewSelect>>", lambda e: self.on_select_spectrum())
        self.tree.bind("<Button-1>", self._tree_click)

        btns = ttk.Frame(parent)
        btns.pack(fill="x", pady=(0, 6))
        ttk.Button(btns, text="Colour", width=8,
                   command=self.choose_colour).pack(side="left")
        ttk.Button(btns, text="Remove", width=8,
                   command=self.remove_spectrum).pack(side="left", padx=4)

        proc = ttk.LabelFrame(parent, text="Processing", padding=6)
        proc.pack(fill="x")

        row = 0
        ttk.Label(proc, text="Line broadening (Hz)").grid(row=row, column=0,
                                                          sticky="w", columnspan=2)
        row += 1
        self.var_lb = tk.DoubleVar(value=0.3)
        ttk.Spinbox(proc, from_=0.0, to=20.0, increment=0.1, width=8,
                    textvariable=self.var_lb,
                    command=self.apply_processing).grid(row=row, column=0, sticky="w")
        ttk.Label(proc, text="Zero fill").grid(row=row, column=1, sticky="e")
        row += 1
        self.var_si = tk.StringVar(value="65536")
        ttk.Combobox(proc, textvariable=self.var_si, width=10, state="readonly",
                     values=("16384", "32768", "65536", "131072")
                     ).grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 6))

        row += 1
        ttk.Label(proc, text="Phase PH0").grid(row=row, column=0, sticky="w")
        self.lbl_p0 = ttk.Label(proc, text="0.0")
        self.lbl_p0.grid(row=row, column=1, sticky="e")
        row += 1
        self.var_p0 = tk.DoubleVar(value=0.0)
        tk.Scale(proc, from_=-180, to=180, resolution=0.5, orient="horizontal",
                 variable=self.var_p0, showvalue=False, length=200,
                 command=lambda v: self._phase_changed()
                 ).grid(row=row, column=0, columnspan=2, sticky="ew")
        row += 1
        ttk.Label(proc, text="Phase PH1").grid(row=row, column=0, sticky="w")
        self.lbl_p1 = ttk.Label(proc, text="0.0")
        self.lbl_p1.grid(row=row, column=1, sticky="e")
        row += 1
        self.var_p1 = tk.DoubleVar(value=0.0)
        tk.Scale(proc, from_=-360, to=360, resolution=1.0, orient="horizontal",
                 variable=self.var_p1, showvalue=False, length=200,
                 command=lambda v: self._phase_changed()
                 ).grid(row=row, column=0, columnspan=2, sticky="ew")

        row += 1
        self.var_baseline = tk.BooleanVar(value=False)
        ttk.Checkbutton(proc, text="Baseline correction",
                        variable=self.var_baseline,
                        command=self.apply_processing
                        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(6, 0))
        row += 1
        self.var_bl_method = tk.StringVar(value="spline")
        ttk.Combobox(proc, textvariable=self.var_bl_method, width=10,
                     state="readonly", values=("spline", "poly"),
                     ).grid(row=row, column=0, columnspan=2, sticky="ew")
        row += 1
        self.var_apply_all = tk.BooleanVar(value=False)
        ttk.Checkbutton(proc, text="Apply to all spectra",
                        variable=self.var_apply_all
                        ).grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1
        ttk.Button(proc, text="Apply / Reprocess", command=self.apply_processing
                   ).grid(row=row, column=0, columnspan=2, sticky="ew", pady=(6, 0))
        proc.columnconfigure(0, weight=1)

        ref = ttk.LabelFrame(parent, text="Reference", padding=6)
        ref.pack(fill="x", pady=(8, 0))
        ttk.Button(ref, text="Calibrate to solvent peak",
                   command=self.calibrate_to_solvent).pack(fill="x")
        self.lbl_solvent = ttk.Label(ref, text="", foreground="#666666")
        self.lbl_solvent.pack(anchor="w", pady=(2, 6))
        ttk.Label(ref, text="Or click a peak in Reference mode,\n"
                            "then set its true shift:").pack(anchor="w")
        rowf = ttk.Frame(ref)
        rowf.pack(fill="x", pady=(4, 0))
        self.var_ref_target = tk.StringVar(value="0.00")
        ttk.Entry(rowf, textvariable=self.var_ref_target, width=8).pack(side="left")
        ttk.Label(rowf, text="ppm").pack(side="left", padx=4)
        ttk.Button(rowf, text="Apply", command=self.apply_reference).pack(side="right")
        ttk.Button(ref, text="Clear calibration", command=self.clear_reference
                   ).pack(fill="x", pady=(4, 0))

    def _build_tabs(self, tabs):
        # Integrals
        frame = ttk.Frame(tabs)
        tabs.add(frame, text="Integrals")
        cols = ("range", "center", "raw", "norm", "protons", "multiplet", "j")
        self.int_tree = ttk.Treeview(frame, columns=cols, show="headings")
        for key, title, width in (
            ("range", "Range (ppm)", 130), ("center", "Centre (ppm)", 95),
            ("raw", "Integral", 110), ("norm", "Relative", 80),
            ("protons", "H", 50), ("multiplet", "Multiplet", 90),
            ("j", "J (Hz)", 150),
        ):
            self.int_tree.heading(key, text=title)
            self.int_tree.column(key, width=width, anchor="center")
        self.int_tree.pack(side="left", fill="both", expand=True)
        self.int_tree.bind("<<TreeviewSelect>>", self._region_selected)
        ttk.Scrollbar(frame, orient="vertical", command=self.int_tree.yview
                      ).pack(side="left", fill="y")
        side = ttk.Frame(frame, padding=4)
        side.pack(side="left", fill="y")
        ttk.Label(side, text="Set H =").pack(anchor="w")
        self.var_protons = tk.StringVar(value="1")
        ttk.Entry(side, textvariable=self.var_protons, width=6).pack(anchor="w")
        ttk.Button(side, text="Assign", command=self.assign_protons
                   ).pack(fill="x", pady=2)
        ttk.Button(side, text="Fit lines", command=self.fit_region
                   ).pack(fill="x", pady=2)
        ttk.Button(side, text="Delete", command=self.delete_region
                   ).pack(fill="x", pady=2)
        ttk.Button(side, text="Clear all", command=self.clear_regions
                   ).pack(fill="x", pady=2)

        # Peaks
        frame = ttk.Frame(tabs)
        tabs.add(frame, text="Peaks")
        cols = ("ppm", "hz", "height", "rel")
        self.peak_tree = ttk.Treeview(frame, columns=cols, show="headings")
        for key, title, width in (("ppm", "Shift (ppm)", 120),
                                  ("hz", "Shift (Hz)", 120),
                                  ("height", "Height", 140),
                                  ("rel", "Relative (%)", 110)):
            self.peak_tree.heading(key, text=title)
            self.peak_tree.column(key, width=width, anchor="center")
        self.peak_tree.pack(side="left", fill="both", expand=True)
        ttk.Scrollbar(frame, orient="vertical", command=self.peak_tree.yview
                      ).pack(side="left", fill="y")
        side = ttk.Frame(frame, padding=4)
        side.pack(side="left", fill="y")
        ttk.Label(side, text="Sensitivity\n(x noise)").pack(anchor="w")
        self.var_sens = tk.DoubleVar(value=8.0)
        ttk.Spinbox(side, from_=2.0, to=100.0, increment=1.0, width=7,
                    textvariable=self.var_sens).pack(anchor="w")
        ttk.Button(side, text="Pick peaks", command=self.pick_peaks_all
                   ).pack(fill="x", pady=4)
        self.peak_tree.bind("<<TreeviewSelect>>", self._peak_selected)

        # Fit
        frame = ttk.Frame(tabs)
        tabs.add(frame, text="Fit")
        self.lbl_fit = ttk.Label(frame, text="Select a region in the Integrals "
                                             "tab and press \"Fit lines\".",
                                 foreground="#666666")
        self.lbl_fit.pack(anchor="w", padx=6, pady=(4, 2))
        holder = ttk.Frame(frame)
        holder.pack(fill="both", expand=True)
        cols = ("ppm", "fwhm", "area", "rel")
        self.fit_tree = ttk.Treeview(holder, columns=cols, show="headings")
        for key, title, width in (("ppm", "Shift (ppm)", 130),
                                  ("fwhm", "FWHM (Hz)", 120),
                                  ("area", "Area", 150),
                                  ("rel", "Relative", 110)):
            self.fit_tree.heading(key, text=title)
            self.fit_tree.column(key, width=width, anchor="center")
        self.fit_tree.pack(side="left", fill="both", expand=True)
        ttk.Scrollbar(holder, orient="vertical", command=self.fit_tree.yview
                      ).pack(side="left", fill="y")

        # Report
        frame = ttk.Frame(tabs)
        tabs.add(frame, text="Report")
        self.report = tk.Text(frame, wrap="word", height=8,
                              font=("TkFixedFont", 11))
        self.report.pack(side="left", fill="both", expand=True)
        ttk.Scrollbar(frame, orient="vertical", command=self.report.yview
                      ).pack(side="left", fill="y")

        # Info
        frame = ttk.Frame(tabs)
        tabs.add(frame, text="Info")
        self.info = tk.Text(frame, wrap="word", height=8,
                            font=("TkFixedFont", 11))
        self.info.pack(fill="both", expand=True)

    # -- spectrum management ------------------------------------------------

    def new_session(self):
        """Start over, after offering to save whatever is open."""
        if not self.confirm_discard():
            return
        self.spectra = []
        self.spectra_2d = []
        self._refresh_2d_menu()
        self.session_path = None
        self._undo_stack.clear()
        self._redo_stack.clear()
        self._refresh_undo_labels()
        self.refresh_tree()
        self.refresh_tables()
        self.plot._marker = None
        self.reset_zoom()
        self.dirty = False
        self._update_title()
        self.status.set("New session")

    def open_files(self):
        paths = filedialog.askopenfilenames(
            title="Open NMR data", filetypes=DATA_FILETYPES,
            initialdir=self.prefs.get("last_dir") or None)
        for path in paths:
            self._load(path)

    def _register_2d(self, path, spectra):
        """Remember any 2D experiments and offer to open the viewer."""
        for spec in spectra:
            if all(spec.source != s.source or spec.name != s.name
                   for s in self.spectra_2d):
                self.spectra_2d.append(spec)
        self._refresh_2d_menu()
        self.prefs.add_recent(path)
        self.prefs.save()
        self.refresh_recent_menu()
        names = ", ".join(s.name for s in spectra)
        self.status.set("Found %d 2D spectrum/spectra (%s) - Tools > 2D Spectra"
                        % (len(spectra), names))
        if len(self.spectra_2d) == len(spectra):
            self.open_2d(spectra[0])

    def _refresh_2d_menu(self):
        self.twod_menu.delete(0, "end")
        if not self.spectra_2d:
            self.twod_menu.add_command(label="(none loaded)", state="disabled")
            return
        for spec in self.spectra_2d:
            self.twod_menu.add_command(
                label="%s  [%s]" % (spec.name, spec.meta.get("Pulse program", "")),
                command=lambda s=spec: self.open_2d(s))

    def open_structure(self):
        """Open (or raise) the structure window."""
        existing = getattr(self, "_structure_window", None)
        if existing is not None and existing.winfo_exists():
            existing.deiconify()
            existing.lift()
            existing.focus_set()
            return existing
        window = structwin.StructureWindow(self, self)
        self._structure_window = window
        window.focus_set()
        return window

    def open_2d(self, spec):
        window = plot2d.Plot2DWindow(self, spec)
        window.focus_set()
        return window

    def open_folder(self):
        path = filedialog.askdirectory(
            title="Open a Bruker experiment folder",
            initialdir=self.prefs.get("last_dir") or None)
        if path:
            self._load(path)

    def open_path(self, path):
        """Open a path from the recent menu, picking the right handler."""
        if not os.path.exists(path):
            messagebox.showerror(APP_NAME, "No longer there:\n\n%s" % path)
            self.refresh_recent_menu()
            return
        if path.lower().endswith((".nmrs", ".json")):
            self._open_session_path(path)
        else:
            self._load(path)

    def _load(self, path):
        two_d = nmr2d.load_2d(path)
        try:
            loaded = nmrio.load(path)
        except Exception as exc:                       # surfaced to the user
            if two_d:
                # A 2D-only dataset has no 1D spectrum to fall back on.
                self._register_2d(path, two_d)
                return
            messagebox.showerror(APP_NAME, "Could not read %s\n\n%s"
                                 % (os.path.basename(path), exc))
            return
        if two_d:
            self._register_2d(path, two_d)
        for spec in loaded:
            spec.color = PALETTE[len(self.spectra) % len(PALETTE)]
            self.spectra.append(spec)
        self.refresh_tree()
        if loaded:
            self.tree.selection_set(str(self.spectra.index(loaded[-1])))
            self.on_select_spectrum()
            self.reset_zoom()
            self.mark_dirty()
        self.prefs.add_recent(path)
        self.prefs.save()
        self.refresh_recent_menu()
        self.status.set("Loaded %d spectrum/spectra from %s"
                        % (len(loaded), os.path.basename(path)))

    def refresh_tree(self):
        selected = self.tree.selection()
        self.tree.delete(*self.tree.get_children())
        for i, spec in enumerate(self.spectra):
            self.tree.insert("", "end", iid=str(i), text=spec.name,
                             values=("Y" if spec.visible else "",),
                             tags=("s%d" % i,))
            self.tree.tag_configure("s%d" % i, foreground=spec.color)
        if selected and selected[0] in self.tree.get_children():
            self.tree.selection_set(selected)

    def _tree_click(self, event):
        if self.tree.identify_column(event.x) != "#1":
            return
        item = self.tree.identify_row(event.y)
        if not item:
            return
        spec = self.spectra[int(item)]
        spec.visible = not spec.visible
        self.refresh_tree()
        self.plot.redraw()
        return "break"

    def active_spectrum(self):
        sel = self.tree.selection()
        if not sel:
            return self.spectra[0] if self.spectra else None
        return self.spectra[int(sel[0])]

    def on_select_spectrum(self):
        spec = self.active_spectrum()
        if not spec:
            return
        self.var_lb.set(spec.lb)
        self.var_si.set(str(spec.si))
        self.var_p0.set(spec.p0)
        self.var_p1.set(spec.p1)
        self.var_baseline.set(spec.baseline_on)
        self.var_bl_method.set(spec.baseline_method)
        self.lbl_p0.config(text="%.1f" % spec.p0)
        self.lbl_p1.config(text="%.1f" % spec.p1)
        solvent = solvents.identify(spec.meta.get("Solvent", ""))
        self.lbl_solvent.config(
            text="%s, residual line at %.2f ppm" % (solvent.label, solvent.primary)
            if solvent else "Solvent not recognised (%r)"
                            % spec.meta.get("Solvent", ""))
        self.refresh_tables()
        self.plot.redraw()

    def choose_colour(self):
        spec = self.active_spectrum()
        if not spec:
            return
        _, hexval = colorchooser.askcolor(color=spec.color, title="Trace colour")
        if hexval:
            spec.color = hexval
            self.refresh_tree()
            self.plot.redraw()

    def remove_spectrum(self):
        spec = self.active_spectrum()
        if not spec:
            return
        self.spectra.remove(spec)
        self.mark_dirty()
        self.refresh_tree()
        self.refresh_tables()
        self.plot.redraw()

    # -- view ---------------------------------------------------------------

    def _sync_view_options(self):
        self.plot.stack = self.var_stack.get()
        self.plot.normalise_each = self.var_norm.get()
        self.plot.show_peaks = self.var_show_peaks.get()
        self.plot.show_integrals = self.var_show_integrals.get()
        self.plot.show_grid = self.var_grid.get()
        self.plot.redraw()

    def scale_intensity(self, factor):
        self.plot.y_scale *= factor
        self.plot.redraw()
        self.status.set("Intensity x%.2f" % self.plot.y_scale)

    def reprocess_all(self):
        if not self.spectra:
            return
        self.config(cursor="watch")
        self.update_idletasks()
        try:
            for spec in self.spectra:
                spec.reprocess()
                self.reintegrate(spec)
        finally:
            self.config(cursor="")
        self.refresh_tables()
        self.plot.redraw()
        self.status.set("Reprocessed %d spectra" % len(self.spectra))

    def clear_peaks(self):
        spec = self.active_spectrum()
        if not spec or not spec.peaks:
            return
        self.push_undo("Clear Peak List")
        spec.peaks = []
        self.refresh_tables()
        self.plot.redraw()

    def copy_report(self):
        text = self.report.get("1.0", "end").strip()
        if not text:
            self.status.set("Nothing to copy.")
            return
        self.clipboard_clear()
        self.clipboard_append(text)
        self.status.set("Report copied to the clipboard")

    def export_report(self):
        text = self.report.get("1.0", "end").strip()
        if not text:
            messagebox.showinfo(APP_NAME, "Build a report first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".txt", filetypes=[("Text", "*.txt")],
            initialdir=self.prefs.get("last_dir") or None)
        if not path:
            return
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text + "\n")
        self.status.set("Wrote %s" % path)

    # -- help ---------------------------------------------------------------

    def show_about(self):
        messagebox.showinfo(
            "About %s" % APP_NAME,
            "%s %s\n\n"
            "A viewer and analyser for 1D NMR spectra.\n"
            "Reads Bruker folders and zips, ACD/Labs .esp and JCAMP-DX.\n\n"
            "Python %s, Tk %s\nStandard library only.\n\nMIT licence."
            % (APP_NAME, APP_VERSION,
               "%d.%d.%d" % sys.version_info[:3],
               self.tk.call("info", "patchlevel")))

    def show_shortcuts(self):
        m = ACCEL_TEXT
        rows = [
            ("File", ""),
            ("New / Open data / Open session", "%sN  %sO  %sShift+O" % (m, m, m)),
            ("Save / Save As", "%sS  %sShift+S" % (m, m)),
            ("Quit", "%sQ" % m),
            ("", ""),
            ("Edit", ""),
            ("Undo / Redo", "%sZ  %sShift+Z" % (m, m)),
            ("Copy report", "%sC" % m),
            ("", ""),
            ("View", ""),
            ("Reset zoom", "%sR" % m),
            ("Taller / shorter", "%s+  %s-" % (m, m)),
            ("", ""),
            ("Tools", ""),
            ("Auto phase", "%sP" % m),
            ("Pick peaks", "%sK" % m),
            ("Auto integrate", "%sI" % m),
            ("Fit lines", "%sF" % m),
            ("Calibrate to solvent", "%sL" % m),
            ("Mouse mode zoom/integrate/peak/reference", "1  2  3  4"),
            ("", ""),
            ("Mouse", ""),
            ("Zoom to selection / create region", "drag"),
            ("Pan", "right-drag"),
            ("Zoom in and out", "wheel"),
            ("Intensity", "shift+wheel"),
            ("Reset zoom", "double-click"),
        ]
        lines = []
        for left, right in rows:
            if not left:
                lines.append("")
            elif not right:
                lines.append(left)
            else:
                lines.append("   %-42s %s" % (left, right))
        messagebox.showinfo("Keyboard Shortcuts", "\n".join(lines))

    def reset_zoom(self):
        visible = [s for s in self.spectra if s.visible]
        if not visible:
            self.plot.set_view(12.0, -1.0)
            return
        left = max(s.limits[0] for s in visible)
        right = min(s.limits[1] for s in visible)
        self.plot.y_scale = 1.0
        self.plot.set_view(left, right)

    # -- processing ---------------------------------------------------------

    def _phase_changed(self):
        spec = self.active_spectrum()
        if not spec:
            return
        self.lbl_p0.config(text="%.1f" % self.var_p0.get())
        self.lbl_p1.config(text="%.1f" % self.var_p1.get())
        if getattr(self, "_phase_job", None):
            self.after_cancel(self._phase_job)
        self._phase_job = self.after(120, self.apply_processing)

    def apply_processing(self):
        spec = self.active_spectrum()
        if not spec:
            return
        self._phase_job = None
        targets = self.spectra if self.var_apply_all.get() else [spec]
        for target in targets:
            target.lb = float(self.var_lb.get())
            target.si = int(self.var_si.get())
            target.p0 = float(self.var_p0.get())
            target.p1 = float(self.var_p1.get())
            target.baseline_on = bool(self.var_baseline.get())
            target.baseline_method = self.var_bl_method.get()
        if not spec.has_fid() and (spec.lb or spec.si != spec.npoints):
            self.status.set("Line broadening and zero filling need the raw FID; "
                            "phase and baseline still apply.")
        self.config(cursor="watch")
        self.update_idletasks()
        try:
            for target in targets:
                target.reprocess()
                self.reintegrate(target)
        finally:
            self.config(cursor="")
        self.mark_dirty()
        self.refresh_tables()
        self.plot.redraw()

    def auto_phase(self):
        spec = self.active_spectrum()
        if not spec:
            return
        self.config(cursor="watch")
        self.update_idletasks()
        try:
            base = spec.complex_spectrum()
            p0, p1 = dsp.autophase(base)
        finally:
            self.config(cursor="")
        self.var_p0.set(round(spec.p0 + p0, 1))
        self.var_p1.set(round(spec.p1 + p1, 1))
        self.apply_processing()
        self.status.set("Auto phase: PH0 = %.1f, PH1 = %.1f"
                        % (spec.p0, spec.p1))

    # -- reference ----------------------------------------------------------

    def set_reference(self, ppm):
        spec = self.active_spectrum()
        if not spec:
            return
        # Snap to the tallest point within +/- 0.05 ppm of the click.
        lo, hi = spec.index(ppm + 0.05), spec.index(ppm - 0.05)
        lo, hi = spec.clamp(min(lo, hi)), spec.clamp(max(lo, hi))
        if hi > lo:
            best = max(range(lo, hi + 1), key=lambda i: spec.real[i])
            ppm = spec.ppm(best)
        self._ref_click = ppm
        self.plot._marker = ppm
        self.plot.redraw()
        self.status.set("Reference peak at %.4f ppm - enter its true shift and "
                        "press Apply" % ppm)

    def apply_reference(self):
        spec = self.active_spectrum()
        clicked = getattr(self, "_ref_click", None)
        if not spec or clicked is None:
            messagebox.showinfo(APP_NAME, "Click a peak in Reference mode first.")
            return
        try:
            target = float(self.var_ref_target.get())
        except ValueError:
            messagebox.showerror(APP_NAME, "Enter a numeric ppm value.")
            return
        self._recalibrate(spec, target - clicked)
        self._ref_click = None
        self.plot._marker = None
        self.status.set("Calibrated: offset %+.4f ppm" % spec.ref_shift)

    def clear_reference(self):
        spec = self.active_spectrum()
        if spec:
            self._recalibrate(spec, -spec.ref_shift)
            self.plot._marker = None
            self.status.set("Calibration cleared")

    def _recalibrate(self, spec, delta):
        """Move the ppm axis by ``delta`` and carry the annotations with it."""
        if not delta:
            return
        self.push_undo("Calibrate")
        spec.ref_shift += delta
        for peak in spec.peaks:
            peak.ppm += delta
        for region in spec.regions:
            region.lo += delta
            region.hi += delta
        self.reintegrate(spec)
        self.refresh_tables()
        self.plot.redraw()

    # -- analysis -----------------------------------------------------------

    def view_window(self):
        """The ppm range currently on screen, low first."""
        return self.plot.view_right, self.plot.view_left

    def pick_peaks_all(self):
        """Pick peaks across the displayed region (reset zoom for the lot)."""
        spec = self.active_spectrum()
        if not spec:
            return
        self.push_undo("Pick Peaks")
        lo, hi = self.view_window()
        spec.peaks = analysis.pick_peaks(spec, lo, hi,
                                         sensitivity=self.var_sens.get())
        self.refresh_tables()
        self.plot.redraw()
        self.status.set("Found %d peaks between %.2f and %.2f ppm"
                        % (len(spec.peaks), hi, lo))

    def pick_peaks_here(self, ppm):
        """Pick peaks in a narrow window around the click, and add them."""
        spec = self.active_spectrum()
        if not spec:
            return
        span = abs(self.plot.view_left - self.plot.view_right) * 0.05
        found = analysis.pick_peaks(spec, ppm - span, ppm + span,
                                    sensitivity=self.var_sens.get())
        self.push_undo("Pick Peak")
        known = {round(p.ppm, 4) for p in spec.peaks}
        for peak in found:
            if round(peak.ppm, 4) not in known:
                spec.peaks.append(peak)
        spec.peaks.sort(key=lambda p: -p.ppm)
        self.refresh_tables()
        self.plot.redraw()

    def add_region(self, lo, hi):
        spec = self.active_spectrum()
        if not spec:
            return
        self.push_undo("Add Region")
        region = analysis.Region(lo, hi)
        analysis.integrate_region(spec, region)
        spec.regions.append(region)
        spec.regions.sort(key=lambda r: -r.center)
        self.refresh_tables()
        self.plot.redraw()

    def auto_integrate(self, max_regions=40):
        """Build integration regions from the peaks in the displayed region."""
        spec = self.active_spectrum()
        if not spec:
            return
        lo_view, hi_view = self.view_window()
        peaks = [p for p in spec.peaks if lo_view <= p.ppm <= hi_view]
        if not peaks:
            self.pick_peaks_all()
            peaks = [p for p in spec.peaks if lo_view <= p.ppm <= hi_view]
        if not peaks:
            self.status.set("No peaks found in the displayed region.")
            return

        # Group peaks that sit within ~20 Hz of each other into one multiplet.
        gap_ppm = 20.0 / spec.sf
        peaks = sorted(peaks, key=lambda p: -p.ppm)
        groups = [[peaks[0]]]
        for peak in peaks[1:]:
            if abs(groups[-1][-1].ppm - peak.ppm) <= gap_ppm:
                groups[-1].append(peak)
            else:
                groups.append([peak])

        # Keep the most intense groups if the region is crowded.
        truncated = len(groups) > max_regions
        if truncated:
            groups.sort(key=lambda g: -sum(p.height for p in g))
            groups = groups[:max_regions]

        self.push_undo("Auto Integrate")
        pad = 6.0 / spec.sf
        spec.regions = []
        for group in groups:
            hi = max(p.ppm for p in group) + pad
            lo = min(p.ppm for p in group) - pad
            region = analysis.Region(lo, hi)
            analysis.integrate_region(spec, region)
            spec.regions.append(region)
        spec.regions.sort(key=lambda r: -r.center)
        self.refresh_tables()
        self.plot.redraw()
        note = " (strongest %d kept)" % max_regions if truncated else ""
        self.status.set("Created %d integration regions%s"
                        % (len(spec.regions), note))

    def reintegrate(self, spec):
        analysis.refresh_peaks(spec)
        for region in spec.regions:
            analysis.integrate_region(spec, region)
            region.fit = None       # the data changed; any old fit is stale

    def fit_region(self):
        """Deconvolute the selected region into individual lines."""
        spec = self.active_spectrum()
        sel = self.int_tree.selection()
        if not spec or not sel:
            messagebox.showinfo(APP_NAME, "Select a region in the Integrals tab "
                                          "first.")
            return
        region = spec.regions[int(sel[0])]
        self.config(cursor="watch")
        self.update_idletasks()
        try:
            result = fitting.fit_region(spec, region.lo, region.hi,
                                        seed_peaks=region.peaks)
        finally:
            self.config(cursor="")
        if result is None:
            self.status.set("Region too narrow to fit.")
            return
        region.fit = result
        self.show_fit(region)
        self.plot.redraw()

    def show_fit(self, region):
        result = region.fit
        self.fit_tree.delete(*self.fit_tree.get_children())
        if result is None:
            self.lbl_fit.config(text="No fit for this region.")
            return
        smallest = min((p.area for p in result.peaks if p.area > 0), default=1.0)
        for i, peak in enumerate(result.peaks):
            self.fit_tree.insert("", "end", iid=str(i),
                                 values=("%.4f" % peak.ppm,
                                         "%.2f" % peak.fwhm_hz,
                                         "%.4g" % peak.area,
                                         "%.2f" % (peak.area / smallest)))
        quality = "converged" if result.converged else \
                  "stopped after %d iterations" % result.iterations
        self.lbl_fit.config(
            text="%d lines, %s, residual %.1f%% of the tallest line, "
                 "lineshape %.0f%% Lorentzian.  Fitted areas exclude the "
                 "baseline, so they run below the plain integral."
                 % (len(result.peaks), quality, 100.0 * result.rel_rms,
                    100.0 * result.eta))
        self.status.set("Fitted %d lines in %.3f-%.3f ppm"
                        % (len(result.peaks), region.hi, region.lo))

    def calibrate_to_solvent(self):
        """Move the axis so the residual solvent line sits at its book value."""
        spec = self.active_spectrum()
        if not spec:
            return
        found = solvents.calibrate(spec)
        if found is None:
            messagebox.showinfo(
                APP_NAME,
                "Could not identify the solvent for this spectrum.\n\n"
                "Metadata says: %r\n\nUse Reference mode to calibrate by hand."
                % spec.meta.get("Solvent", ""))
            return
        delta, at_ppm, solvent = found
        self._recalibrate(spec, delta)
        self.status.set("%s line found at %.4f ppm, moved to %.2f (%+.4f ppm)"
                        % (solvent.label, at_ppm, solvent.primary, delta))

    def assign_protons(self):
        spec = self.active_spectrum()
        sel = self.int_tree.selection()
        if not spec or not sel:
            return
        try:
            value = float(self.var_protons.get())
        except ValueError:
            messagebox.showerror(APP_NAME, "Enter a numeric proton count.")
            return
        self.push_undo("Assign Protons")
        for item in sel:
            spec.regions[int(item)].protons = value
        self.refresh_tables()
        self.plot.redraw()

    def delete_region(self):
        spec = self.active_spectrum()
        sel = self.int_tree.selection()
        if not spec or not sel:
            return
        self.push_undo("Delete Region")
        for item in sorted((int(i) for i in sel), reverse=True):
            del spec.regions[item]
        self.refresh_tables()
        self.plot.redraw()

    def clear_regions(self):
        spec = self.active_spectrum()
        if spec and spec.regions:
            self.push_undo("Clear Regions")
            spec.regions = []
            self.refresh_tables()
            self.plot.redraw()

    def _region_selected(self, _event):
        """Show whichever fit belongs to the region the user just selected."""
        spec = self.active_spectrum()
        sel = self.int_tree.selection()
        if not spec or not sel:
            return
        index = int(sel[0])
        if index < len(spec.regions):
            self.show_fit(spec.regions[index])

    def _peak_selected(self, _event):
        spec = self.active_spectrum()
        sel = self.peak_tree.selection()
        if not spec or not sel:
            return
        peak = spec.peaks[int(sel[0])]
        self.plot._marker = peak.ppm
        self.plot.redraw()

    # -- tables -------------------------------------------------------------

    def refresh_tables(self):
        spec = self.active_spectrum()
        self.int_tree.delete(*self.int_tree.get_children())
        self.peak_tree.delete(*self.peak_tree.get_children())
        self.fit_tree.delete(*self.fit_tree.get_children())
        self.lbl_fit.config(text="Select a region in the Integrals tab and "
                                 "press \"Fit lines\".")
        self.report.delete("1.0", "end")
        self.info.delete("1.0", "end")
        if not spec:
            return

        norm = analysis.normalise(spec.regions)
        for i, (region, value) in enumerate(zip(spec.regions, norm)):
            m = region.multiplet
            self.int_tree.insert(
                "", "end", iid=str(i),
                values=("%.3f - %.3f" % (region.hi, region.lo),
                        "%.4f" % (m.center_ppm if m else region.center),
                        "%.4g" % region.value,
                        "%.2f" % value,
                        "%g" % region.protons if region.protons else "",
                        m.pattern if m else "",
                        ", ".join("%.1f" % j for j in m.couplings) if m else ""))

        top = max((p.height for p in spec.peaks), default=1.0) or 1.0
        for i, peak in enumerate(spec.peaks):
            self.peak_tree.insert(
                "", "end", iid=str(i),
                values=("%.4f" % peak.ppm,
                        "%.2f" % (peak.ppm * spec.sf),
                        "%.4g" % peak.height,
                        "%.1f" % (100.0 * peak.height / top)))

        self.report.insert("1.0", self._build_report(spec))
        self.info.insert("1.0", self._build_info(spec))

    def _build_report(self, spec):
        if not spec.regions:
            return ("Drag across a signal in Integrate mode, or press "
                    "Auto Integrate, to build a report.")

        lines = [analysis.format_report(spec, spec.regions), "",
                 "Molar ratio (relative integrals)"]
        norm = analysis.normalise(spec.regions)
        for region, value in zip(spec.regions, norm):
            lines.append("  %7.3f ppm   %8.3f%s"
                         % (region.center, value,
                            "   (assigned %gH)" % region.protons
                            if region.protons else ""))

        # Composition needs at least two regions with a known proton count.
        assigned = [r for r in spec.regions if r.protons]
        if len(assigned) >= 2:
            components = [analysis.Component("%.2f ppm" % r.center, r, r.protons)
                          for r in assigned]
            lines.append("")
            lines.append("Composition")
            lines.append(analysis.format_composition(components))
            lines.append("")
            lines.append("(The first assigned region, highest ppm, is taken as "
                         "starting material.)")
        else:
            lines.append("")
            lines.append("Assign a proton count to two or more regions to get "
                         "mol %% and conversion.")
        return "\n".join(lines)

    def _build_info(self, spec):
        left, right = spec.limits
        lines = ["Name:     %s" % spec.name,
                 "Source:   %s" % spec.source,
                 "Points:   %d" % spec.npoints,
                 "Range:    %.3f to %.3f ppm" % (left, right),
                 "Digital resolution: %.3f Hz/point" % spec.hz_per_point(),
                 "Raw FID:  %s" % ("yes (%d complex points)" % len(spec.fid)
                                   if spec.fid else "no"),
                 "Calibration offset: %+.4f ppm" % spec.ref_shift,
                 ""]
        for key, value in spec.meta.items():
            lines.append("%-20s %s" % (key + ":", value))
        return "\n".join(lines)

    # -- export -------------------------------------------------------------

    def export_peaks(self):
        spec = self.active_spectrum()
        if not spec or not spec.peaks:
            messagebox.showinfo(APP_NAME, "Pick some peaks first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")])
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["ppm", "Hz", "height"])
            for peak in spec.peaks:
                writer.writerow(["%.4f" % peak.ppm,
                                 "%.2f" % (peak.ppm * spec.sf),
                                 "%.6g" % peak.height])
        self.status.set("Wrote %s" % path)

    def export_integrals(self):
        spec = self.active_spectrum()
        if not spec or not spec.regions:
            messagebox.showinfo(APP_NAME, "Create some integration regions first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")])
        if not path:
            return
        norm = analysis.normalise(spec.regions)
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["from_ppm", "to_ppm", "centre_ppm", "integral",
                             "relative", "protons", "multiplet", "J_Hz"])
            for region, value in zip(spec.regions, norm):
                m = region.multiplet
                writer.writerow([
                    "%.4f" % region.hi, "%.4f" % region.lo,
                    "%.4f" % (m.center_ppm if m else region.center),
                    "%.6g" % region.value, "%.3f" % value,
                    region.protons or "", m.pattern if m else "",
                    " ".join("%.2f" % j for j in m.couplings) if m else ""])
        self.status.set("Wrote %s" % path)

    def export_plot(self):
        path = filedialog.asksaveasfilename(defaultextension=".ps",
                                            filetypes=[("PostScript", "*.ps")])
        if not path:
            return
        self.plot.postscript(file=path, colormode="color")
        self.status.set("Wrote %s" % path)

    def export_svg(self):
        if not [s for s in self.spectra if s.visible]:
            messagebox.showinfo(APP_NAME, "Nothing to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".svg",
                                            filetypes=[("SVG", "*.svg")])
        if not path:
            return
        try:
            export.write_svg(path, self.spectra,
                             self.plot.view_left, self.plot.view_right,
                             stack=self.plot.stack,
                             normalise_each=self.plot.normalise_each,
                             show_peaks=self.plot.show_peaks,
                             show_integrals=self.plot.show_integrals,
                             active=self.active_spectrum())
        except Exception as exc:
            messagebox.showerror(APP_NAME, "Could not write the SVG:\n\n%s" % exc)
            return
        self.status.set("Wrote %s" % path)

    # -- sessions -----------------------------------------------------------

    def save_session(self):
        """Save to the current file, asking for a name the first time."""
        if not self.spectra:
            messagebox.showinfo(APP_NAME, "Open some data first.")
            return False
        if not self.session_path:
            return self.save_session_as()
        return self._write_session(self.session_path)

    def save_session_as(self):
        if not self.spectra:
            messagebox.showinfo(APP_NAME, "Open some data first.")
            return False
        path = filedialog.asksaveasfilename(
            title="Save session as", defaultextension=".nmrs",
            filetypes=SESSION_FILETYPES,
            initialdir=self.prefs.get("last_dir") or None,
            initialfile=os.path.basename(self.session_path)
            if self.session_path else "session.nmrs")
        if not path:
            return False
        return self._write_session(path)

    def _write_session(self, path):
        try:
            export.save_session(
                path, self.spectra,
                (self.plot.view_left, self.plot.view_right),
                {"stack": self.plot.stack,
                 "normalise_each": self.plot.normalise_each,
                 "show_peaks": self.plot.show_peaks,
                 "show_integrals": self.plot.show_integrals,
                 "sensitivity": float(self.var_sens.get())})
        except Exception as exc:
            messagebox.showerror(APP_NAME, "Could not save:\n\n%s" % exc)
            return False
        self.session_path = path
        self.mark_dirty(False)
        self.prefs.add_recent(path)
        self.prefs.save()
        self.refresh_recent_menu()
        self._update_title()
        self.status.set("Saved to %s" % path)
        return True

    def open_session(self):
        if not self.confirm_discard():
            return
        path = filedialog.askopenfilename(
            title="Open session",
            filetypes=SESSION_FILETYPES + [("All files", "*.*")],
            initialdir=self.prefs.get("last_dir") or None)
        if path:
            self._open_session_path(path)

    def _open_session_path(self, path):
        self.config(cursor="watch")
        self.update_idletasks()
        try:
            spectra, view, options, warnings = export.load_session(path)
        except Exception as exc:
            messagebox.showerror(APP_NAME, "Could not open the session:\n\n%s" % exc)
            return
        finally:
            self.config(cursor="")

        if not spectra:
            messagebox.showerror(APP_NAME, "The session referred to no data that "
                                           "could be opened.\n\n%s"
                                 % "\n".join(warnings))
            return

        self.spectra = spectra
        self.session_path = path
        self._undo_stack.clear()
        self._redo_stack.clear()
        self._refresh_undo_labels()
        self.var_stack.set(bool(options.get("stack", False)))
        self.var_norm.set(bool(options.get("normalise_each", True)))
        self.var_show_peaks.set(bool(options.get("show_peaks", True)))
        self.var_show_integrals.set(bool(options.get("show_integrals", True)))
        self.var_sens.set(float(options.get("sensitivity", 8.0)))
        self._sync_view_options()
        self.refresh_tree()
        self.tree.selection_set("0")
        self.on_select_spectrum()
        self.plot.set_view(*view)
        self.mark_dirty(False)
        self.prefs.add_recent(path)
        self.prefs.save()
        self.refresh_recent_menu()
        if warnings:
            messagebox.showwarning(APP_NAME, "Some sources could not be "
                                             "reopened:\n\n%s" % "\n".join(warnings))
        self.status.set("Restored %d spectrum/spectra from %s"
                        % (len(spectra), os.path.basename(path)))

    def quit_app(self):
        if not self.confirm_discard():
            return
        self._store_preferences()
        self.destroy()


def main():
    import sys
    app = App()
    for path in sys.argv[1:]:
        app._load(path)
    app.mainloop()


if __name__ == "__main__":
    main()
