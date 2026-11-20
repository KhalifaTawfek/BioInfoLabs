"""DNA patterns from promoter sequences: C+G% vs information content over a sliding window."""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

S_TEST = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTCTAGCGATCTA"


# ---------------------------------------------------------------- algorithm
def cg_content(s):
    """C+G percentage of a window."""
    return 100 * (s.count("C") + s.count("G")) / len(s)


def information_content(s):
    """Information content of a window (self-sequence alignment, Kappa IC based).
    The window is compared with itself shifted by u = 1..n-1; for each shift the
    percentage of matching positions is taken; IC = 100 - average of these."""
    n = len(s)
    t = 0
    for u in range(1, n):
        matches = sum(1 for i in range(n - u) if s[i] == s[i + u])
        t += matches / (n - u) * 100
    return 100 - t / (n - 1)


def pattern(seq, window=30):
    """One (C+G%, IC) point for every window position (step 1)."""
    return [(cg_content(seq[i:i + window]), information_content(seq[i:i + window]))
            for i in range(len(seq) - window + 1)]


def center_of_weight(points):
    n = len(points)
    return sum(x for x, _ in points) / n, sum(y for _, y in points) / n


def clean(text):
    """Keep only A, C, G, T (drops FASTA headers, spaces, numbers)."""
    lines = [l for l in text.splitlines() if not l.startswith(">")]
    return "".join(c for c in "".join(lines).upper() if c in "ACGT")


# ---------------------------------------------------------------- app
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DNA Patterns - C+G% vs Information Content")
        self.centers = []  # (name, cg, ic)

        top = ttk.Frame(self, padding=8)
        top.pack(fill="x")
        ttk.Label(top, text="DNA sequence (paste a promoter, or open a FASTA file):").grid(
            row=0, column=0, columnspan=6, sticky="w")
        self.seq = tk.Text(top, height=4, width=100, font=("Consolas", 10), wrap="char")
        self.seq.grid(row=1, column=0, columnspan=6, sticky="we", pady=4)
        self.seq.insert("1.0", S_TEST)

        ttk.Label(top, text="Name:").grid(row=2, column=0, sticky="e")
        self.name = tk.StringVar(value="Test S")
        ttk.Entry(top, textvariable=self.name, width=18).grid(row=2, column=1, sticky="w")
        ttk.Label(top, text="Window:").grid(row=2, column=2, sticky="e")
        self.win = tk.StringVar(value="30")
        ttk.Spinbox(top, from_=5, to=500, width=6, textvariable=self.win).grid(row=2, column=3, sticky="w")
        bar = ttk.Frame(top)
        bar.grid(row=2, column=4, columnspan=2, sticky="e")
        ttk.Button(bar, text="Test sequence S", command=self.load_test).pack(side="left", padx=2)
        ttk.Button(bar, text="Open FASTA...", command=self.open_fasta).pack(side="left", padx=2)
        ttk.Button(bar, text="Generate pattern", command=self.run).pack(side="left", padx=2)
        ttk.Button(bar, text="Clear centers", command=self.clear).pack(side="left", padx=2)

        self.info = ttk.Label(self, padding=(8, 0), font=("Consolas", 10), justify="left")
        self.info.pack(fill="x")

        self.fig = Figure(figsize=(11, 4.8), dpi=90)
        self.ax1 = self.fig.add_subplot(121)
        self.ax2 = self.fig.add_subplot(122)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.run()

    def load_test(self):
        self.seq.delete("1.0", "end")
        self.seq.insert("1.0", S_TEST)
        self.name.set("Test S")

    def open_fasta(self):
        path = filedialog.askopenfilename(filetypes=[("FASTA / text", "*.fa *.fasta *.txt"), ("All", "*.*")])
        if not path:
            return
        with open(path) as f:
            text = f.read()
        self.seq.delete("1.0", "end")
        self.seq.insert("1.0", clean(text))
        first = text.splitlines()[0] if text else ""
        self.name.set(first[1:30].strip() if first.startswith(">") else path.split("/")[-1])

    def clear(self):
        self.centers = []
        self.draw_centers()
        self.canvas.draw()

    def run(self):
        s = clean(self.seq.get("1.0", "end"))
        try:
            w = int(self.win.get())
        except ValueError:
            messagebox.showerror("Error", "Window must be an integer.")
            return
        if len(s) < w:
            messagebox.showerror("Error", f"Sequence ({len(s)} b) is shorter than the window ({w}).")
            return

        pts = pattern(s, w)
        cx, cy = center_of_weight(pts)
        name = self.name.get().strip() or f"Pattern {len(self.centers) + 1}"
        self.centers.append((name, cx, cy))

        self.info.config(text=(
            f"Length = {len(s)} b    windows = {len(pts)}    "
            f"first window: CG = {pts[0][0]:.2f} %   IC = {pts[0][1]:.2f}    "
            f"center of weight: CG = {cx:.2f} %   IC = {cy:.2f}"))

        a = self.ax1
        a.clear()
        a.scatter([p[0] for p in pts], [p[1] for p in pts], s=12, c=range(len(pts)),
                  cmap="viridis", alpha=0.8, label="windows")
        a.scatter([cx], [cy], s=120, c="red", marker="X", label="center of weight")
        a.set(title=f"Pattern: {name}", xlabel="C+G %", ylabel="Information content (IC)",
              xlim=(0, 100), ylim=(0, 100))
        a.grid(alpha=0.3)
        a.legend(loc="upper left", fontsize=8)
        self.draw_centers()
        self.fig.tight_layout()
        self.canvas.draw()

    def draw_centers(self):
        a = self.ax2
        a.clear()
        for name, x, y in self.centers:
            a.scatter([x], [y], s=60)
            a.annotate(name, (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8)
        a.set(title="Centers of weight", xlabel="C+G %", ylabel="Information content (IC)",
              xlim=(0, 100), ylim=(0, 100))
        a.grid(alpha=0.3)


if __name__ == "__main__":
    App().mainloop()
