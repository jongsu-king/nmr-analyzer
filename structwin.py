"""Structure window: draw a SMILES and check it against the integrals.

The question this answers is the everyday one -- "is what I made what I think
I made?".  Enter the structure you were aiming for, and the integrals of the
spectrum are normalised so that their total equals the proton count of that
formula.  Each region then reads directly as a number of protons, which either
matches an expected environment or does not.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import messagebox, ttk

import analysis
import depict
import smiles as smiles_mod


class StructureWindow(tk.Toplevel):
    def __init__(self, master, app):
        super().__init__(master)
        self.app = app
        self.molecule = None
        self.title("Structure")
        self.geometry("980x620")
        self.minsize(760, 500)

        self.smiles_text = tk.StringVar(value="")
        self.status = tk.StringVar(value="Enter a SMILES string and press Draw.")
        self.exchangeable = tk.BooleanVar(value=False)

        self._build()

    # -- construction -------------------------------------------------------

    def _build(self):
        top = ttk.Frame(self, padding=(8, 6))
        top.pack(side="top", fill="x")
        ttk.Label(top, text="SMILES").pack(side="left")
        entry = ttk.Entry(top, textvariable=self.smiles_text)
        entry.pack(side="left", fill="x", expand=True, padx=6)
        entry.bind("<Return>", lambda e: self.draw())
        ttk.Button(top, text="Draw", command=self.draw).pack(side="left")
        ttk.Button(top, text="Export SVG...",
                   command=self.export_svg).pack(side="left", padx=(6, 0))

        examples = ttk.Frame(self, padding=(8, 0))
        examples.pack(side="top", fill="x")
        ttk.Label(examples, text="Examples:",
                  foreground="#666666").pack(side="left")
        for label, text in (("carbazole", "c1ccc2[nH]c3ccccc3c2c1"),
                            ("3-iodocarbazole", "Ic1ccc2[nH]c3ccccc3c2c1"),
                            ("toluene", "Cc1ccccc1"),
                            ("ethyl acetate", "CCOC(C)=O")):
            ttk.Button(examples, text=label, width=len(label) + 2,
                       command=lambda t=text: (self.smiles_text.set(t), self.draw())
                       ).pack(side="left", padx=2)

        panes = ttk.PanedWindow(self, orient="horizontal")
        panes.pack(side="top", fill="both", expand=True)

        left = ttk.Frame(panes)
        panes.add(left, weight=3)
        self.canvas = tk.Canvas(left, background="white", highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Configure>", lambda e: self._render())

        right = ttk.Frame(panes, padding=6)
        panes.add(right, weight=2)

        self.summary = tk.Text(right, height=6, wrap="word",
                               font=("TkFixedFont", 11))
        self.summary.pack(fill="x")

        ttk.Label(right, text="Expected proton environments",
                  font=("TkDefaultFont", 10, "bold")).pack(anchor="w", pady=(8, 2))
        self.env_tree = ttk.Treeview(right, columns=("h", "kind", "sites"),
                                     show="headings", height=6)
        for key, title, width in (("h", "H", 40), ("kind", "Environment", 150),
                                  ("sites", "Sites", 50)):
            self.env_tree.heading(key, text=title)
            self.env_tree.column(key, width=width, anchor="center")
        self.env_tree.pack(fill="both", expand=True)

        ttk.Label(right, text="Measured integrals, scaled to this formula",
                  font=("TkDefaultFont", 10, "bold")).pack(anchor="w", pady=(8, 2))
        self.fit_tree = ttk.Treeview(right, columns=("range", "h", "near"),
                                     show="headings", height=6)
        for key, title, width in (("range", "Range (ppm)", 120),
                                  ("h", "Protons", 70),
                                  ("near", "Matches", 150)):
            self.fit_tree.heading(key, text=title)
            self.fit_tree.column(key, width=width, anchor="center")
        self.fit_tree.pack(fill="both", expand=True)

        row = ttk.Frame(right)
        row.pack(fill="x", pady=(6, 0))
        ttk.Checkbutton(row, text="Ignore OH/NH (exchangeable)",
                        variable=self.exchangeable,
                        command=self.check_against_spectrum).pack(side="left")
        ttk.Button(row, text="Check against spectrum",
                   command=self.check_against_spectrum).pack(side="right")

        ttk.Label(self, textvariable=self.status, padding=(8, 4),
                  wraplength=940, justify="left").pack(side="bottom", fill="x")

    # -- structure ----------------------------------------------------------

    def draw(self):
        text = self.smiles_text.get().strip()
        if not text:
            return
        try:
            self.molecule = smiles_mod.parse(text)
        except smiles_mod.SmilesError as exc:
            self.molecule = None
            self.canvas.delete("all")
            messagebox.showerror("Structure", "Could not read that SMILES:\n\n%s"
                                 % exc)
            self.status.set("SMILES error: %s" % exc)
            return
        self._render()
        self._fill_summary()
        self.check_against_spectrum()

    def _render(self):
        canvas = self.canvas
        canvas.delete("all")
        if self.molecule is None:
            return
        width = max(canvas.winfo_width(), 60)
        height = max(canvas.winfo_height(), 60)

        coords = depict.layout(self.molecule)
        screen, scale = depict.transform(coords, width, height)
        font_size = int(max(9, min(16, scale * 0.34)))

        for x1, y1, x2, y2 in depict.bond_segments(self.molecule, screen, scale):
            canvas.create_line(x1, y1, x2, y2, fill="#222222", width=2)

        for atom in self.molecule.atoms:
            label = depict.visible_label(atom)
            if not label:
                continue
            x, y = screen[atom.index]
            radius = font_size * 0.8
            canvas.create_oval(x - radius, y - radius, x + radius, y + radius,
                               fill="white", outline="")
            canvas.create_text(x, y, text=label, fill="#111111",
                               font=("TkDefaultFont", font_size))

    def _fill_summary(self):
        mol = self.molecule
        self.summary.delete("1.0", "end")
        counts = mol.formula_counts()
        lines = [
            "Formula        %s" % mol.formula(),
            "Molecular mass %.2f" % mol.molecular_weight(),
            "DBE            %g   (rings + pi bonds)" % mol.dbe(),
            "Rings          %d" % mol.ring_count(),
            "Protons        %d in %d environment(s)"
            % (counts.get("H", 0), len(mol.proton_environments())),
        ]
        self.summary.insert("1.0", "\n".join(lines))

        self.env_tree.delete(*self.env_tree.get_children())
        for i, env in enumerate(mol.proton_environments()):
            self.env_tree.insert("", "end", iid=str(i),
                                 values=(env.count, env.label, len(env.atoms)))

    # -- comparison with the spectrum ---------------------------------------

    def _expected_protons(self):
        """Protons the formula predicts, optionally minus exchangeable ones."""
        mol = self.molecule
        total = mol.formula_counts().get("H", 0)
        if not self.exchangeable.get():
            return total, 0
        skipped = 0
        for env in mol.proton_environments():
            if env.carrier.symbol in ("O", "N", "S"):
                skipped += env.count
        return total - skipped, skipped

    def check_against_spectrum(self):
        """Scale the integrals to the formula and match them to environments.

        Proximity to a whole number is a weak test on its own -- with a free
        scale factor almost any structure can be made to land near integers.
        What is diagnostic is whether the *pattern* of integrals matches the
        set of proton environments the structure predicts.
        """
        self.fit_tree.delete(*self.fit_tree.get_children())
        if self.molecule is None:
            return
        spec = self.app.active_spectrum()
        if spec is None or not spec.regions:
            self.status.set("Integrate some regions in the main window first, "
                            "then press Check.")
            return

        expected, skipped = self._expected_protons()
        if expected <= 0:
            self.status.set("This structure has no protons to compare.")
            return

        usable = [r for r in spec.regions if r.value > 0]
        dropped = len(spec.regions) - len(usable)
        total = sum(r.value for r in usable)
        if not usable or total <= 0:
            self.status.set("No region has a positive integral - check the "
                            "phase and the baseline.")
            return

        scale = expected / total
        environments = list(self.molecule.proton_environments())
        if self.exchangeable.get():
            environments = [e for e in environments
                            if e.carrier.symbol not in ("O", "N", "S")]
        unmatched = list(environments)

        rows = []
        for region in sorted(usable, key=lambda r: -r.center):
            protons = region.value * scale
            best = None
            if unmatched:
                best = min(unmatched, key=lambda e: abs(e.count - protons))
                tolerance = max(0.3, 0.15 * best.count)
                if abs(best.count - protons) <= tolerance:
                    unmatched.remove(best)
                else:
                    best = None
            rows.append((region, protons, best))

        for i, (region, protons, env) in enumerate(rows):
            self.fit_tree.insert(
                "", "end", iid=str(i),
                values=("%.3f - %.3f" % (region.hi, region.lo),
                        "%.2f" % protons,
                        "%dH %s" % (env.count, env.label) if env else "-"))

        matched = sum(1 for _r, _p, env in rows if env)
        problems = []
        if dropped:
            problems.append("%d region(s) with a non-positive integral ignored"
                            % dropped)
        if matched < len(rows):
            problems.append("%d region(s) match no predicted environment"
                            % (len(rows) - matched))
        if unmatched:
            problems.append("%d predicted environment(s) unaccounted for (%s)"
                            % (len(unmatched),
                               ", ".join("%dH %s" % (e.count, e.label)
                                         for e in unmatched)))
        # An excluded exchangeable proton is a choice the user made, not a
        # defect in the match, so it is reported separately.
        note = (" %d exchangeable proton(s) were excluded." % skipped
                if skipped else "")

        head = ("Integrals scaled so the total is %d H for %s. "
                % (expected, self.molecule.formula()))
        if not problems:
            self.status.set(head + "All %d region(s) match a predicted proton "
                                   "environment - consistent with this "
                                   "structure.%s" % (len(rows), note))
        else:
            self.status.set(head + "Not a clean match: " + "; ".join(problems)
                            + "." + note
                            + "  This is a consistency check, not proof - "
                              "overlapping signals and missed regions look the "
                              "same as a wrong structure.")

    def export_svg(self):
        if self.molecule is None:
            messagebox.showinfo("Structure", "Draw a structure first.")
            return
        from tkinter import filedialog
        path = filedialog.asksaveasfilename(defaultextension=".svg",
                                            filetypes=[("SVG", "*.svg")])
        if not path:
            return
        svg = depict.render_svg(self.molecule, 460, 360,
                                title="%s   %s" % (self.smiles_text.get().strip(),
                                                   self.molecule.formula()))
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(svg)
        self.status.set("Wrote %s" % path)
