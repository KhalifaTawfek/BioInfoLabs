"""DNA patterns from gene promoter sequences.

For every sliding window (SW) of the sequence we compute:
  * (C+G)%                - percentage of C and G letters in the window
  * information content   - sigma(S), self-sequence alignment (Kappa IC):

        sigma(S) = 100 - [ sum_{u=1}^{|S|-1} ( sum_{i=1}^{|S|-1} f(x_i, x_{u+i}) / (|S|-u) * 100 ) ] / (|S|-1)

        f(x_i, x_{u+i}) = 1 if x_i == x_{u+i}, else 0

Each window gives one point (C+G%, IC); all the points form the pattern.
The center of weight of a pattern is the mean of its points.
"""
import sys
import tkinter as tk
from tkinter import filedialog, font as tkfont

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

S_TEST = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTA"
WINDOW = 30


# ---------------------------------------------------------------- algorithm
def f(a, b):
    return 1 if a == b else 0


def cg_content(s):
    """(C+G)% of a sequence."""
    return 100 * (s.count("C") + s.count("G")) / len(s)


def information_content(s):
    """sigma(S) from the board."""
    n = len(s)
    total = 0
    for u in range(1, n):                        # u = 1 .. |S|-1
        m = 0
        for i in range(n - 1):                   # i = 1 .. |S|-1 (0-based here)
            if i + u < n:                        # x_{u+i} exists
                m += f(s[i], s[i + u])
        total += m / (n - u) * 100
    return 100 - total / (n - 1)


def pattern(seq, w=WINDOW):
    """One point (C+G%, IC) for every sliding window, step 1."""
    return [(cg_content(seq[k:k + w]), information_content(seq[k:k + w]))
            for k in range(len(seq) - w + 1)]


def center_of_weight(points):
    n = len(points)
    return sum(p[0] for p in points) / n, sum(p[1] for p in points) / n


def clean(text):
    """Keep only A, C, G, T (drops FASTA headers, spaces, numbers)."""
    body = "".join(l for l in text.splitlines() if not l.startswith(">"))
    return "".join(c for c in body.upper() if c in "ACGT")


# ---------------------------------------------------------------- theme
BG, CARD, CARD2, LINE = "#0e1117", "#161b26", "#1f2636", "#2a3245"
TEXT, MUTED, ACCENT, ACCENT_HI = "#e6e9ef", "#8b93a7", "#7c5cff", "#9a80ff"
CYAN, PINK, BAD = "#22d3ee", "#f472b6", "#f87171"
PALETTE = ["#22d3ee", "#f472b6", "#34d399", "#fbbf24", "#a78bfa", "#fb923c", "#60a5fa", "#f87171"]
UI, MONO = "TkDefaultFont", "TkFixedFont"


def init_fonts(root):
    global UI, MONO
    fams = set(tkfont.families(root))
    UI = next((x for x in ("Segoe UI", "Inter", "DejaVu Sans") if x in fams), UI)
    MONO = next((x for x in ("Cascadia Mono", "Consolas", "DejaVu Sans Mono") if x in fams), MONO)


def dark_titlebar(win):
    if sys.platform != "win32":
        return
    try:
        import ctypes
        win.update_idletasks()
        hwnd = ctypes.windll.user32.GetParent(win.winfo_id())
        v = ctypes.c_int(1)
        ctypes.windll.dwmapi.DwmSetWindowAttribute(hwnd, 20, ctypes.byref(v), ctypes.sizeof(v))
    except Exception:
        pass


class Pill(tk.Label):
    def __init__(self, master, text, command, primary=False):
        self.base, self.hover = (ACCENT, ACCENT_HI) if primary else (CARD2, LINE)
        super().__init__(master, text=text, bg=self.base, fg=TEXT, font=(UI, 10, "bold"),
                         padx=14, pady=6, cursor="hand2")
        self.bind("<Enter>", lambda e: self.config(bg=self.hover))
        self.bind("<Leave>", lambda e: self.config(bg=self.base))
        self.bind("<Button-1>", lambda e: command())


def card(master, **kw):
    return tk.Frame(master, bg=CARD, highlightthickness=1, highlightbackground=LINE, **kw)


def caption(master, text):
    return tk.Label(master, text=text.upper(), bg=CARD, fg=MUTED, font=(UI, 8, "bold"))


# ---------------------------------------------------------------- app
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        init_fonts(self)
        self.title("DNA Patterns")
        self.configure(bg=BG, padx=18, pady=16)
        dark_titlebar(self)
        self.patterns = []          # (name, points, center, color)

        # header
        tk.Label(self, text="DNA Patterns", bg=BG, fg=TEXT, font=(UI, 18, "bold")).pack(anchor="w")
        tk.Label(self, text="(C+G)% and information content from every sliding window of a gene promoter",
                 bg=BG, fg=MUTED, font=(UI, 10)).pack(anchor="w", pady=(0, 12))

        # input card
        ic = card(self, padx=14, pady=12)
        ic.pack(fill="x")
        caption(ic, "DNA sequence").grid(row=0, column=0, sticky="w")
        self.len_lbl = tk.Label(ic, text="", bg=CARD, fg=MUTED, font=(MONO, 9))
        self.len_lbl.grid(row=0, column=1, sticky="e")
        self.seq = tk.Text(ic, height=4, width=110, font=(MONO, 10), wrap="char", bg=CARD2, fg=TEXT,
                           insertbackground=TEXT, relief="flat", highlightthickness=1,
                           highlightbackground=LINE, highlightcolor=ACCENT, padx=8, pady=6)
        self.seq.grid(row=1, column=0, columnspan=2, sticky="we", pady=(6, 10))
        self.seq.insert("1.0", S_TEST)
        self.seq.bind("<KeyRelease>", lambda e: self._update_len())
        ic.columnconfigure(0, weight=1)

        row = tk.Frame(ic, bg=CARD)
        row.grid(row=2, column=0, columnspan=2, sticky="we")
        caption(row, "Name").pack(side="left")
        self.name = tk.StringVar(value="Test S")
        self._entry(row, self.name, 18).pack(side="left", padx=(6, 16), ipady=4)
        caption(row, "Window").pack(side="left")
        self.win = tk.StringVar(value=str(WINDOW))
        self._entry(row, self.win, 5).pack(side="left", padx=(6, 16), ipady=4)
        Pill(row, "Generate pattern", self.run, primary=True).pack(side="right")
        Pill(row, "Clear", self.clear).pack(side="right", padx=6)
        Pill(row, "Open FASTA", self.open_fasta).pack(side="right")
        Pill(row, "Test sequence S", self.load_test).pack(side="right", padx=6)
        self.err = tk.Label(ic, text="", bg=CARD, fg=BAD, font=(UI, 9))
        self.err.grid(row=3, column=0, columnspan=2, sticky="w")

        # stat cards
        sc = tk.Frame(self, bg=BG)
        sc.pack(fill="x", pady=10)
        self.stats = {}
        for k, name in enumerate(("Length", "Windows", "(C+G)% of S", "IC of S",
                                  "Center (C+G)%", "Center IC")):
            c = card(sc, padx=12, pady=8)
            c.grid(row=0, column=k, sticky="nsew", padx=(0 if k == 0 else 8, 0))
            sc.columnconfigure(k, weight=1)
            caption(c, name).pack(anchor="w")
            v = tk.Label(c, text="-", bg=CARD, fg=TEXT, font=(UI, 16, "bold"))
            v.pack(anchor="w")
            self.stats[name] = v

        # charts
        chc = card(self, padx=8, pady=8)
        chc.pack(fill="both", expand=True)
        self.fig = Figure(figsize=(12, 5), dpi=90, facecolor=CARD)
        self.ax1 = self.fig.add_subplot(121)
        self.ax2 = self.fig.add_subplot(122)
        self.canvas = FigureCanvasTkAgg(self.fig, master=chc)
        self.canvas.get_tk_widget().configure(bg=CARD, highlightthickness=0)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        tk.Label(self, text="IC:  σ(S) = 100 − [ Σu=1..|S|−1 ( Σi=1..|S|−1 "
                            "f(xᵢ, xᵤ₊ᵢ) / (|S|−u) × 100 ) ] / (|S|−1)"
                            "      f = 1 if xᵢ = xᵤ₊ᵢ, else 0",
                 bg=BG, fg=MUTED, font=(MONO, 9)).pack(anchor="w", pady=(8, 0))
        self.run()

    def _entry(self, master, var, width):
        return tk.Entry(master, textvariable=var, width=width, font=(MONO, 11), bg=CARD2, fg=TEXT,
                        insertbackground=TEXT, relief="flat", highlightthickness=1,
                        highlightbackground=LINE, highlightcolor=ACCENT)

    def _update_len(self):
        self.len_lbl.config(text=f"{len(clean(self.seq.get('1.0', 'end')))} b")

    # ---------- actions
    def load_test(self):
        self.seq.delete("1.0", "end")
        self.seq.insert("1.0", S_TEST)
        self.name.set("Test S")
        self._update_len()

    def open_fasta(self):
        path = filedialog.askopenfilename(filetypes=[("FASTA / text", "*.fa *.fasta *.fna *.txt"),
                                                     ("All files", "*.*")])
        if not path:
            return
        with open(path) as fh:
            text = fh.read()
        self.seq.delete("1.0", "end")
        self.seq.insert("1.0", clean(text))
        first = text.splitlines()[0] if text.strip() else ""
        self.name.set(first[1:25].strip() if first.startswith(">")
                      else path.replace("\\", "/").split("/")[-1])
        self._update_len()

    def clear(self):
        self.patterns = []
        for v in self.stats.values():
            v.config(text="-")
        self.draw()

    def run(self):
        self._update_len()
        s = clean(self.seq.get("1.0", "end"))
        try:
            w = int(self.win.get())
            assert w >= 2
        except (ValueError, AssertionError):
            self.err.config(text="Window must be a whole number ≥ 2.")
            return
        if len(s) < w:
            self.err.config(text=f"The sequence ({len(s)} b) is shorter than the window ({w}).")
            return
        self.err.config(text="")

        pts = pattern(s, w)
        center = center_of_weight(pts)
        name = self.name.get().strip() or f"Pattern {len(self.patterns) + 1}"
        color = PALETTE[len(self.patterns) % len(PALETTE)]
        self.patterns.append((name, pts, center, color))

        vals = {"Length": f"{len(s)} b", "Windows": len(pts),
                "(C+G)% of S": f"{cg_content(s):.2f}", "IC of S": f"{information_content(s):.2f}",
                "Center (C+G)%": f"{center[0]:.2f}", "Center IC": f"{center[1]:.2f}"}
        for k, v in vals.items():
            self.stats[k].config(text=str(v))
        self.draw()

    # ---------- charts
    def _style(self, ax, title):
        ax.set_facecolor(CARD2)
        ax.set_title(title, color=TEXT, fontsize=11, loc="left", pad=10)
        ax.set_xlabel("(C+G) %", color=MUTED)
        ax.set_ylabel("Information content (IC)", color=MUTED)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.grid(True, color=LINE, linewidth=0.6)
        for side in ax.spines.values():
            side.set_color(LINE)

    def draw(self):
        a1, a2 = self.ax1, self.ax2
        a1.clear()
        a2.clear()
        if self.patterns:
            name, pts, (cx, cy), color = self.patterns[-1]
            a1.scatter([p[0] for p in pts], [p[1] for p in pts], s=14, color=color, alpha=0.75,
                       edgecolors="none", label="sliding windows")
            a1.scatter([cx], [cy], s=160, color="white", marker="X", edgecolors=BG, linewidths=1,
                       label=f"center ({cx:.2f}, {cy:.2f})", zorder=5)
            self._style(a1, f"Pattern — {name}")
            leg = a1.legend(loc="lower right", fontsize=8, facecolor=CARD, edgecolor=LINE)
            for t in leg.get_texts():
                t.set_color(TEXT)
        else:
            self._style(a1, "Pattern")

        for name, _, (cx, cy), color in self.patterns:
            a2.scatter([cx], [cy], s=90, color=color, edgecolors="white", linewidths=0.8, zorder=3)
            a2.annotate(name, (cx, cy), textcoords="offset points", xytext=(7, 6),
                        fontsize=8, color=TEXT)
        self._style(a2, "Centers of weight (one per pattern)")
        self.fig.tight_layout(pad=1.5)
        self.canvas.draw()


if __name__ == "__main__":
    App().mainloop()
