"""
bioinfolab3.py - Genome Codon Analyzer (desktop window app)

Part 1: Reads a genome (FASTA or plain text), slides a window of length 3
        with step 1 across it, counts every overlapping codon in a
        dictionary, and shows the top 10 codons with counts and frequencies.
Part 2: Maps the top codons to amino acids (standard genetic code), asks
        the question about foods low in those amino acids, and answers it
        from built-in general nutritional knowledge.

Run:
    python bioinfolab3.py                 # opens the window, pick a file
    python bioinfolab3.py genome.fna      # opens the window and analyzes it
Requires only the Python standard library (tkinter ships with Python).
"""

import os
import sys
from collections import Counter
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

VALID_BASES = set("ACGT")
TOP_N = 10

# ===========================================================================
# Standard genetic code (DNA codon -> amino acid), built from the TCAG table
# ===========================================================================
_BASES = "TCAG"
_AA_ONE_LETTER = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_ONE_TO_NAME = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic acid",
    "C": "Cysteine", "E": "Glutamic acid", "Q": "Glutamine", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
    "*": "Stop",
}
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
GENETIC_CODE = {cod: _ONE_TO_NAME[aa] for cod, aa in zip(_CODONS, _AA_ONE_LETTER)}

# ===========================================================================
# Nutrition knowledge: food groups that are relatively LOW in each amino acid
# (general nutritional knowledge, approximate)
# ===========================================================================
LOW_FOODS = {
    "Lysine":        ["Cereal grains (wheat, white bread, rice, corn, oats)", "Nuts (almonds)"],
    "Threonine":     ["Cereal grains (wheat, white bread, rice, corn, oats)", "Most fruits"],
    "Tryptophan":    ["Corn (maize) and cornmeal", "Gelatin", "Most fruits"],
    "Methionine":    ["Legumes (beans, lentils, peas)", "Potatoes", "Most fruits"],
    "Cysteine":      ["Legumes (beans, lentils, peas)", "Potatoes", "Most fruits"],
    "Isoleucine":    ["Corn (maize) and cornmeal", "Gelatin", "Most fruits"],
    "Leucine":       ["Gelatin", "Leafy greens", "Most fruits"],
    "Valine":        ["Gelatin", "Most fruits"],
    "Phenylalanine": ["Most fruits"],
    "Tyrosine":      ["Most fruits"],
    "Histidine":     ["Most fruits", "Lettuce and cucumbers"],
    "Arginine":      ["Most fruits", "Milk and dairy (relatively low)"],
    "Glycine":       ["Most fruits", "Leafy greens"],
    "Proline":       ["Most fruits", "Leafy greens"],
}
_DEFAULT_LOW = ["Most fruits"]
# Foods that are genuinely almost protein-free, so low in every amino acid.
_PROTEIN_FREE = "Sugar, sweets, soft drinks, butter and cooking oils"


# ===========================================================================
# Part 1 - analysis logic
# ===========================================================================
def read_genome(path):
    """Read FASTA or plain-text DNA; return one uppercase string.
    Header/comment lines ('>' or ';') and all whitespace are ignored."""
    parts = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in ">;":
                continue
            parts.append("".join(line.split()).upper())
    return "".join(parts)


def count_codons(sequence):
    """Sliding window of length 3, step 1. Windows with any non-ACGT
    character are skipped. Returns (codon -> count dict, skipped count)."""
    counts = Counter()
    skipped = 0
    for i in range(len(sequence) - 2):
        codon = sequence[i:i + 3]
        if VALID_BASES.issuperset(codon):
            counts[codon] += 1
        else:
            skipped += 1
    return dict(counts), skipped


def top_codons(codon_counts, n=TOP_N):
    """Top n codons as [(codon, count, percent)]; ties sorted alphabetically."""
    total = sum(codon_counts.values())
    ranked = sorted(codon_counts.items(), key=lambda kv: (-kv[1], kv[0]))[:n]
    return [(c, k, 100.0 * k / total) for c, k in ranked]


def analyze(path):
    """Run Part 1 on a file. Raises ValueError with a readable message."""
    seq = read_genome(path)
    if len(seq) < 3:
        raise ValueError(f"The sequence has only {len(seq)} base(s); at least 3 are needed.")
    counts, skipped = count_codons(seq)
    if not counts:
        raise ValueError("No valid A/C/G/T codons were found in this file.")
    return {
        "length": len(seq),
        "total": sum(counts.values()),
        "skipped": skipped,
        "distinct": len(counts),
        "top": top_codons(counts),
        "seq": seq,
        "counts": counts,
    }


# ===========================================================================
# Part 2 - amino acids, question and answer
# ===========================================================================
def unique_amino_acids(top):
    """Distinct amino acids in rank order (stop codons excluded)."""
    result = []
    for codon, _, _ in top:
        aa = GENETIC_CODE[codon]
        if aa != "Stop" and aa not in result:
            result.append(aa)
    return result


def build_question(amino_acids):
    return ("Given these amino acids that appear most frequently according to "
            f"the codon counts: {', '.join(amino_acids)}, which common foods "
            "are relatively low in or lack these amino acids?")


def answer_question(amino_acids):
    """Answer the food question: group foods and list which of the
    given amino acids each one is low in."""
    food_to_aas = {}
    for aa in amino_acids:
        for food in LOW_FOODS.get(aa, _DEFAULT_LOW):
            food_to_aas.setdefault(food, []).append(aa)

    # Foods that are low in the most amino acids come first.
    ranked = sorted(food_to_aas.items(), key=lambda kv: -len(kv[1]))
    lines = ["Common foods that are relatively low in these amino acids:", ""]
    for food, aas in ranked:
        lines.append(f"  • {food}: low in {', '.join(aas)}")
    lines += [
        f"  • {_PROTEIN_FREE}: contain almost no protein, so they lack "
        "all of these amino acids",
        "",
        "In short: cereal grains (especially corn) are the most notable foods "
        "low in lysine and threonine, while fruits, fats and sugars are low in "
        "nearly all amino acids because they contain very little protein.",
        "Meat, fish, eggs, dairy, soy and quinoa are rich in all of them, and "
        "combining grains with legumes (e.g. rice and beans) fills most gaps.",
    ]
    return "\n".join(lines)


# ===========================================================================
# GUI
# ===========================================================================
class CodonApp(tk.Tk):
    BAR_COLOR = "#3b7dd8"

    def __init__(self, initial_path=None):
        super().__init__()
        self.title("Genome Codon Analyzer")
        self.geometry("940x900")
        self.minsize(760, 700)
        self._seq, self._counts, self._playing = "", {}, False
        self._build_widgets()
        if initial_path:
            self.after(100, lambda: self.run_analysis(initial_path))

    # ---- layout -----------------------------------------------------------
    def _build_widgets(self):
        style = ttk.Style(self)
        style.configure("Treeview", rowheight=24)
        style.configure("Header.TLabel", font=("Segoe UI", 11, "bold"))

        top = ttk.Frame(self, padding=10)
        top.pack(fill="x")
        ttk.Button(top, text="Open genome file…", command=self.choose_file).pack(side="left")
        self.file_label = ttk.Label(top, text="No file loaded", foreground="#666")
        self.file_label.pack(side="left", padx=10)

        self.stats_label = ttk.Label(self, text="", padding=(10, 0))
        self.stats_label.pack(fill="x")

        # Sliding window viewer: slider over genome positions
        win = ttk.Frame(self, padding=(10, 8, 10, 0))
        win.pack(fill="x")
        ttk.Label(win, text="Sliding window (length 3, step 1)",
                  style="Header.TLabel").pack(anchor="w")
        self.seq_view = tk.Text(win, height=1, font=("Consolas", 13), relief="flat",
                                bg="#f4f6fa", padx=6, pady=6, wrap="none")
        self.seq_view.tag_configure("window", background="#ffd54f",
                                    font=("Consolas", 13, "bold"))
        self.seq_view.tag_configure("center", justify="center")
        self.seq_view.pack(fill="x", pady=(4, 2))
        self.seq_view.configure(state="disabled")

        self.pos_var = tk.IntVar(value=0)
        self.slider = tk.Scale(win, from_=0, to=0, orient="horizontal",
                               variable=self.pos_var, showvalue=False,
                               command=lambda v: self.update_window(), state="disabled")
        self.slider.pack(fill="x")

        ctrl = ttk.Frame(win)
        ctrl.pack(fill="x", pady=(2, 0))
        ttk.Button(ctrl, text="\u25c0 Step", width=8,
                   command=lambda: self.step(-1)).pack(side="left")
        self.play_btn = ttk.Button(ctrl, text="\u25b6 Play", width=8, command=self.toggle_play)
        self.play_btn.pack(side="left", padx=4)
        ttk.Button(ctrl, text="Step \u25b6", width=8,
                   command=lambda: self.step(1)).pack(side="left")
        ttk.Label(ctrl, text="Speed").pack(side="left", padx=(14, 4))
        self.speed = tk.Scale(ctrl, from_=1, to=50, orient="horizontal", length=120,
                              showvalue=False)
        self.speed.set(10)
        self.speed.pack(side="left")
        self.window_label = ttk.Label(ctrl, text="", font=("Segoe UI", 10))
        self.window_label.pack(side="left", padx=14)

        # Part 1: table + bar chart side by side
        mid = ttk.Frame(self, padding=10)
        mid.pack(fill="both", expand=False)
        ttk.Label(mid, text="Top 10 codons (sliding window, step = 1)",
                  style="Header.TLabel").pack(anchor="w")
        row = ttk.Frame(mid)
        row.pack(fill="x", pady=5)

        cols = ("rank", "codon", "count", "freq", "aa")
        self.table = ttk.Treeview(row, columns=cols, show="headings", height=10)
        for col, title, width, anchor in [
            ("rank", "#", 40, "center"), ("codon", "Codon", 70, "center"),
            ("count", "Count", 80, "e"), ("freq", "Freq %", 80, "e"),
            ("aa", "Amino acid", 130, "w"),
        ]:
            self.table.heading(col, text=title)
            self.table.column(col, width=width, anchor=anchor, stretch=False)
        self.table.pack(side="left")

        self.chart = tk.Canvas(row, height=260, bg="white", highlightthickness=1,
                               highlightbackground="#ccc")
        self.chart.pack(side="left", fill="both", expand=True, padx=(10, 0))
        self.chart.bind("<Configure>", lambda e: self.draw_chart())
        self._top = []

        # Part 2: question + answer
        bottom = ttk.Frame(self, padding=10)
        bottom.pack(fill="both", expand=True)
        ttk.Label(bottom, text="Question", style="Header.TLabel").pack(anchor="w")
        self.question = tk.Text(bottom, height=3, wrap="word", font=("Segoe UI", 10),
                                bg="#f4f6fa", relief="flat", padx=8, pady=6)
        self.question.pack(fill="x", pady=(2, 8))
        ttk.Label(bottom, text="Answer", style="Header.TLabel").pack(anchor="w")
        self.answer = tk.Text(bottom, wrap="word", font=("Segoe UI", 10),
                              relief="flat", padx=8, pady=6)
        sb = ttk.Scrollbar(bottom, command=self.answer.yview)
        sb.pack(side="right", fill="y")
        self.answer.configure(yscrollcommand=sb.set)
        self.answer.pack(fill="both", expand=True, pady=(2, 0))
        for box in (self.question, self.answer):
            box.configure(state="disabled")

    # ---- actions ----------------------------------------------------------
    def choose_file(self):
        path = filedialog.askopenfilename(
            title="Select genome file",
            filetypes=[("FASTA / DNA files", "*.fasta *.fa *.fna *.ffn *.txt *.seq"),
                       ("All files", "*.*")])
        if path:
            self.run_analysis(path)

    def run_analysis(self, path):
        if not os.path.isfile(path):
            messagebox.showerror("File not found", path)
            return
        self.config(cursor="watch")
        self.update_idletasks()
        try:
            result = analyze(path)
        except (ValueError, OSError) as exc:
            messagebox.showerror("Cannot analyze file", str(exc))
            return
        finally:
            self.config(cursor="")
        self.show_results(path, result)

    def show_results(self, path, r):
        self.file_label.config(text=os.path.basename(path), foreground="black")
        self.stats_label.config(
            text=f"Genome length: {r['length']:,} bases   |   Codons counted: "
                 f"{r['total']:,}   |   Distinct codons: {r['distinct']}   |   "
                 f"Skipped (non-ACGT): {r['skipped']:,}")

        self.table.delete(*self.table.get_children())
        for i, (codon, count, pct) in enumerate(r["top"], 1):
            self.table.insert("", "end", values=(i, codon, f"{count:,}",
                                                 f"{pct:.3f}", GENETIC_CODE[codon]))
        self._top = r["top"]
        self.draw_chart()

        self._seq, self._counts = r["seq"], r["counts"]
        self.slider.configure(state="normal", to=len(self._seq) - 3)
        self.pos_var.set(0)
        self.update_window()

        aas = unique_amino_acids(r["top"])
        self._set_text(self.question, build_question(aas))
        self._set_text(self.answer, answer_question(aas))

    # ---- sliding window ---------------------------------------------------
    CONTEXT = 30  # bases shown on each side of the window

    def update_window(self):
        if not self._seq:
            return
        i = self.pos_var.get()
        lo = max(0, i - self.CONTEXT)
        hi = min(len(self._seq), i + 3 + self.CONTEXT)
        before, codon, after = self._seq[lo:i], self._seq[i:i + 3], self._seq[i + 3:hi]
        pad = " " * (self.CONTEXT - len(before))  # keeps the window centred

        v = self.seq_view
        v.configure(state="normal")
        v.delete("1.0", "end")
        v.insert("end", pad + before, "center")
        v.insert("end", codon, ("center", "window"))
        v.insert("end", after, "center")
        v.configure(state="disabled")

        if set(codon) <= VALID_BASES:
            info = (f"{codon} \u2192 {GENETIC_CODE[codon]}   "
                    f"(found {self._counts.get(codon, 0):,} times in genome)")
        else:
            info = f"{codon} contains non-ACGT bases \u2192 skipped"
        self.window_label.config(
            text=f"Positions {i + 1}-{i + 3} of {len(self._seq):,}:   {info}")

    def step(self, delta):
        if self._seq:
            new = min(max(0, self.pos_var.get() + delta), len(self._seq) - 3)
            self.pos_var.set(new)
            self.update_window()

    def toggle_play(self):
        if not self._seq:
            return
        self._playing = not self._playing
        self.play_btn.config(text="\u23f8 Pause" if self._playing else "\u25b6 Play")
        if self._playing:
            self._tick()

    def _tick(self):
        if not self._playing:
            return
        if self.pos_var.get() >= len(self._seq) - 3:
            self.toggle_play()  # reached the end
            return
        self.step(1)
        self.after(max(10, 1000 // self.speed.get()), self._tick)

    def draw_chart(self):
        c = self.chart
        c.delete("all")
        if not self._top:
            c.create_text(c.winfo_width() / 2, c.winfo_height() / 2,
                          text="Open a genome file to see results", fill="#888")
            return
        w, h = c.winfo_width(), c.winfo_height()
        left, right, top_pad, bottom_pad = 50, 15, 15, 30
        n = len(self._top)
        max_pct = max(p for _, _, p in self._top)
        slot = (w - left - right) / n
        c.create_line(left, h - bottom_pad, w - right, h - bottom_pad, fill="#999")
        for i, (codon, _, pct) in enumerate(self._top):
            bar_h = (h - top_pad - bottom_pad) * pct / max_pct
            x0 = left + i * slot + slot * 0.15
            x1 = left + (i + 1) * slot - slot * 0.15
            y0 = h - bottom_pad - bar_h
            c.create_rectangle(x0, y0, x1, h - bottom_pad, fill=self.BAR_COLOR, width=0)
            c.create_text((x0 + x1) / 2, y0 - 8, text=f"{pct:.1f}", font=("Segoe UI", 8))
            c.create_text((x0 + x1) / 2, h - bottom_pad + 12, text=codon,
                          font=("Consolas", 9))
        c.create_text(12, h / 2, text="%", angle=90, fill="#666")

    @staticmethod
    def _set_text(box, text):
        box.configure(state="normal")
        box.delete("1.0", "end")
        box.insert("1.0", text)
        box.configure(state="disabled")


if __name__ == "__main__":
    CodonApp(sys.argv[1] if len(sys.argv) > 1 else None).mainloop()