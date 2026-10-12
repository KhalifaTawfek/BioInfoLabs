import random
import tkinter as tk
from tkinter import ttk

CODON_TABLE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

# 3-letter code -> (1-letter code, full name)
AMINO_INFO = {
    "Phe": ("F", "Phenylalanine"), "Leu": ("L", "Leucine"),
    "Ile": ("I", "Isoleucine"), "Met": ("M", "Methionine"),
    "Val": ("V", "Valine"), "Ser": ("S", "Serine"),
    "Pro": ("P", "Proline"), "Thr": ("T", "Threonine"),
    "Ala": ("A", "Alanine"), "Tyr": ("Y", "Tyrosine"),
    "His": ("H", "Histidine"), "Gln": ("Q", "Glutamine"),
    "Asn": ("N", "Asparagine"), "Lys": ("K", "Lysine"),
    "Asp": ("D", "Aspartic acid"), "Glu": ("E", "Glutamic acid"),
    "Cys": ("C", "Cysteine"), "Trp": ("W", "Tryptophan"),
    "Arg": ("R", "Arginine"), "Gly": ("G", "Glycine"),
}

START_CODON = "AUG"
STOP_CODONS = ("UAA", "UAG", "UGA")


def transcribe(dna_gene: str) -> str:
    """Transcription: copy the entire gene into RNA form (T -> U)."""
    return dna_gene.upper().replace("T", "U")


def random_rna(length: int = 100) -> str:
    """Random RNA of exactly `length` letters (A, C, G, U)."""
    return "".join(random.choice("ACGU") for _ in range(length))


def translate(rna: str):
    """
    Translation: find the first AUG (start, Met). From there take 3 letters
    at a time, look each group up in the table, until a stop codon
    (UAA / UAG / UGA) is reached.
    Returns (rows, protein, status); rows = [(codon, amino_acid), ...].
    """
    start = rna.find(START_CODON)
    if start == -1:
        return [], "", "no start codon (AUG) found"

    rows = []
    protein = ""
    for i in range(start, len(rna) - 2, 3):
        codon = rna[i:i + 3]
        aa = CODON_TABLE[codon]
        if aa == "STOP":
            rows.append((codon, "STOP"))
            return rows, protein, f"stop codon {codon} reached"
        rows.append((codon, aa))
        protein += "-" + aa
    return rows, protein, "end of sequence reached without a stop codon"


def split_triplets(seq: str) -> str:
    """Show a sequence in groups of 3 separated by '-', e.g. GUU-GCA-UAA."""
    return "-".join(seq[i:i + 3] for i in range(0, len(seq), 3))


def random_rna_with_gene(length: int = 100) -> str:
    """Random RNA (exactly `length` letters) that contains AUG ... in-frame stop,
    so the demo always produces a result. Retries until it does."""
    while True:
        rna = random_rna(length)
        rows, _, status = translate(rna)
        if rows and rows[-1][1] == "STOP":
            return rna


def build_output(rna: str) -> str:
    """Text shown when the user presses 'Show output'."""
    rows, protein, status = translate(rna)
    if not rows:
        return status

    start = rna.find(START_CODON)
    lines = [
        "starting from AUG:",
        (rna[:start] + " | " if start else "") + split_triplets(rna[start:]),
        "",
        f"Translation starts at position {start} (AUG)",
        "",
        "Codons (AUG to stop):",
        "-".join(codon for codon, _ in rows),
        "",
        f"Status : {status}",
    ]
    return "\n".join(lines)


def build_protein_text(rna: str) -> str:
    """Text for the 'Name proteins' window: amino acid name for every codon."""
    rows, protein, status = translate(rna)
    if not rows:
        return status

    lines = ["Codon  Amino acid"]
    for codon, aa in rows:
        lines.append(f"{codon}    {'-' + aa if aa != 'STOP' else 'STOP'}")
    lines += ["", f"Protein: {protein}", f"Status : {status}"]

    # Conversion summary: letter sequence -> amino acid names
    amino = [aa for _, aa in rows if aa != "STOP"]
    codons = [c for c, aa in rows if aa != "STOP"]
    lines += [
        "",
        "Letters sequence to amino acids:",
        "-".join(codons),
        "-> " + "-".join(AMINO_INFO[a][1] for a in amino),
        "-> " + "".join(AMINO_INFO[a][0] for a in amino) + "  (one-letter codes)",
    ]
    return "\n".join(lines)


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("RNA Translation")
        self.geometry("720x480")
        self.minsize(520, 360)

        self.rna = ""

        buttons = ttk.Frame(self, padding=10)
        buttons.pack(fill="x")
        ttk.Button(buttons, text="New sequence", command=self.new_sequence).pack(side="left")
        ttk.Button(buttons, text="Show output", command=self.show_output).pack(side="left", padx=8)
        ttk.Button(buttons, text="Name proteins", command=self.name_proteins).pack(side="left")

        ttk.Label(self, text="RNA sequence", padding=(10, 0)).pack(anchor="w")
        self.seq_box = self._make_text(height=4)

        ttk.Label(self, text="Output", padding=(10, 8, 10, 0)).pack(anchor="w")
        self.out_box = self._make_text(height=12, expand=True)

        self.new_sequence()

    def _make_text(self, height, expand=False):
        box = tk.Text(self, height=height, wrap="word", font=("Courier New", 11),
                      state="disabled", padx=6, pady=6)
        box.pack(fill="both" if expand else "x", expand=expand, padx=10, pady=(2, 6))
        return box

    @staticmethod
    def _set(box, text):
        box.config(state="normal")
        box.delete("1.0", "end")
        box.insert("1.0", text)
        box.config(state="disabled")

    def new_sequence(self):
        self.rna = random_rna_with_gene(100)
        self._set(self.seq_box, f"{self.rna}\n({len(self.rna)} letters)")
        self._set(self.out_box, "")

    def show_output(self):
        self._set(self.out_box, build_output(self.rna))

    def name_proteins(self):
        """Open a new window with the amino acid name for every codon."""
        win = tk.Toplevel(self)
        win.title("Protein names")
        win.geometry("460x520")
        win.minsize(360, 300)

        box = tk.Text(win, wrap="word", font=("Courier New", 11), padx=8, pady=8)
        box.pack(fill="both", expand=True, padx=10, pady=10)
        box.insert("1.0", build_protein_text(self.rna))
        box.config(state="disabled")

        ttk.Button(win, text="Close", command=win.destroy).pack(pady=(0, 10))


if __name__ == "__main__":
    App().mainloop()