#!/usr/bin/env python3
"""
Greedy de novo assembly – live visualiser (tkinter, stdlib only)

1. Fetches a real 1–3 kb DNA sequence from NCBI (default NC_005148.1,
   Porcine circovirus 2 complete genome, 1,768 bp). Falls back to an embedded
   copy of the same record if offline.
2. Samples 2,000 random reads of 100–120 bp, stored in a list.
3. Reassembles with a greedy overlap-layout-consensus assembler
   (containment removal -> suffix/prefix overlap graph -> greedy longest-first merges).
4. Plays the assembly back live (speed slider + draggable timeline) and writes
   an analysis of why reconstruction fails.

Run:  python assembly_viz.py
"""
import random
import threading
import urllib.request
import webbrowser
import tkinter as tk
from tkinter import ttk
from collections import defaultdict

# ---------------------------------------------------------------- parameters
ACCESSION = "NC_005148.1"
N_READS, RMIN, RMAX = 2000, 100, 120
MIN_OV = 20                      # minimum suffix/prefix overlap accepted
EMBEDDED_HEADER = "NC_005148.1 Porcine circovirus 2, complete genome"
EMBEDDED_SEQ = "".join("""
ACCAGCGCACTTCGGCAGCGGCAGCACCTCGGCAGCACCTCAGCAGCAACATGCCCAGCAAGAAGAATGGAAGAAGCGGA
CCCCAACCACACAAAAGGTGGGTGTTCACGCTGAATAATCCTTCCGAAGACGAGCGCAAGAAAATACGGGAGCTTCCAAT
CTCCCTTTTTGATTATTTTATTGTTGGCGAGGAGGGTAATGAGGAAGGACGAACACCCCACCTCCAGGGGTTCGCTAATT
TTGTGAAGAAGCAAACATTTAATAAAGTGAAATGGTATTTCGGTGCCCGCTGCCACATCGAGAAAGCGAAAGGAACTGAT
CAGCAGAATAAAGAATACTGCAGTAAAGAAGGCAACTTACTGATGGAATGTGGAGCTCCTAGATCTCAAGGACAACGGAG
TGACCTGTCTACTGCTGTGAGTACCTTGTTGGAGAGCGGGAGTCTGGTGACCGTTGCAGAGCAGCACCCTGTAACGTTTG
TCAGAAATTTCCGCGGGCTGGCTGAACTTTTGAAAGTGAGCGGGAAAATGCAGAAGCGTGATTGGAAGACCAATGTACAC
GTCATTGTGGGGCCACCTGGGTGTGGTAAAAGCAAATGGGCTGCTAATTTTGCAGACCCGGAAACCACATACTGGAAACC
ACCTAGAAACAAGTGGTGGGACGGTTACCATGGTGAAGAAGTGGTTGTTATTGATGACTTTTATGGCTGGCTGCCGTGGG
ATGATCTACTGAGATTGTGTGATCGATATCCATTGACTGTAGAGACTAAAGGTGGAACTGTACCTTTTTTGGCCCGCAGT
ATTCTGATTACCAGCAATCAGACCCCGTTGGAATGGTACTCCTCAACTGCTGTCCCAGCTGTAGAAGCTCTCTATCGGAG
GATTACTTCCTTGGTATTTTGGAAGAATGCTACAGAACAATCCACGGAGGAAGGGGGCCAGTTCGTCACCCTTTCCCCCC
CATGCCCTGAATTTCCATATGAAATAAATTACTGAGTCTTTTTTATCACTTCGTAATGGTTTTTATTATTCATTTAGGGT
TTAAGTGGGGGGTCTTTAAGATTAAATTCTCTGAATTGTACATACATGGTTACACGGATATTGTAGTCCTGGTCGTATTT
ACTGTTTTCGAACGCAGTGCCGAGGCCTACGTGGTCCACATTTCTACTGGTTTGTAGTCTCATCCACAGCTGATTCCTTT
TGTTATTTGGTTGGAAGTAATCAATAGTGGAATCAAGAACAGGTTTGGGTGTGAAGTAACGGGAGTGGTAGGAGAAGGGT
TGGGGGATTGTATGGCGGGAGGAGTAGTTTACATAGGGGTCATAGGTTTGGGCTGTGGCCTTTATTACAAAGTTGTCATC
TAGAATAATAGCACTGGATCCAACTCCCCTGTCACCCTGGGTGATCGGGGAGCAGGGCCAGAATTCAACCTTAACCTTTC
TTATTCTGTAGTATTCAAAGGGTATAGAGATTTTGTTGGTCCCCCCTCCCGGGGGAACAAAGTCGTCAAGATTAAATCTC
AGCATGTCCACCGCCCAGGAGGGCGTGCTGACTGTGGTAGCCTTGACAGTATATCCGAAGGTGCGGGAGAGGCGGGCGTT
GAAGATGCCATTTTTCCTTCTCCAGCGGTAACGGTGGCGGGGGTGGACGAGCCAGGGGCGGCGGCGGAGGATCTGGCCAA
GATGGCTGCGGGGGCGGTGTCTTCTCCTTCGGTAACGCCTCCTTGGATACGTCATACCTGAAAACGAAAGAAGTGCGCTG
TAAGTATT
""".split())

# ---------------------------------------------------------------- theme
BG, PANEL, PANEL2 = "#0e1320", "#151c2c", "#1d2638"
FG, MUTED, GRID = "#e6e9f0", "#8a93a8", "#2a3449"
ACC, GREEN, PURPLE = "#4cc9f0", "#43d17a", "#b48cff"
ORANGE, RED, YEL = "#ffb347", "#ff5c7a", "#ffe066"
UI, UIB, MONO = ("Segoe UI", 10), ("Segoe UI Semibold", 10), ("Consolas", 10)
KIND_COL = {"fwd": GREEN, "rev": PURPLE, "multi": ORANGE, "chimeric": RED}
KIND_TXT = {"fwd": "placed", "rev": "placed (rev-comp)",
            "multi": "ambiguous placement (repeat)", "chimeric": "CHIMERIC"}

BASE_COL = {"A": "#43d17a", "C": "#4cc9f0", "G": "#ffe066", "T": "#ff5c7a"}
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


# ================================================================ core logic
def fetch_ncbi(acc, timeout=10):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={acc}&rettype=fasta&retmode=text")
    req = urllib.request.Request(url, headers={"User-Agent": "assembly-demo/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        txt = r.read().decode()
    lines = txt.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("NCBI did not return FASTA")
    return lines[0][1:].strip(), "".join(lines[1:]).upper()


def simulate_reads(ref, rng, both_strands=False):
    """Return reads (plain list of strings) + hidden truth (start, len, strand)."""
    reads, truth = [], []
    for _ in range(N_READS):
        ln = rng.randint(RMIN, RMAX)
        st = rng.randint(0, len(ref) - ln)
        r, strand = ref[st:st + ln], "+"
        if both_strands and rng.random() < 0.5:
            r, strand = revcomp(r), "-"
        reads.append(r)
        truth.append((st, ln, strand))
    return reads, truth


def coverage(L, truth):
    d = [0] * (L + 1)
    for s, ln, _ in truth:
        d[s] += 1
        d[s + ln] -= 1
    cov, c = [], 0
    for i in range(L):
        c += d[i]
        cov.append(c)
    return cov


def remove_contained(reads):
    uniq = sorted(set(reads), key=lambda r: (-len(r), r))
    kept, big = [], ""
    for r in uniq:
        if r not in big:
            kept.append(r)
            big += "|" + r
    return kept


def find_overlaps(reads, k=MIN_OV):
    idx = defaultdict(list)
    for j, b in enumerate(reads):
        idx[b[:k]].append(j)
    edges = []
    for i, a in enumerate(reads):
        for p in range(1, len(a) - k + 1):
            suf = a[p:]
            for j in idx.get(a[p:p + k], ()):
                if j != i and reads[j].startswith(suf):
                    edges.append((len(suf), i, j))
    return edges


def place(s, ref):
    p = ref.find(s)
    if p >= 0:
        return (p, p + len(s), "multi" if ref.find(s, p + 1) >= 0 else "fwd")
    rc = revcomp(s)
    p = ref.find(rc)
    if p >= 0:
        return (p, p + len(s), "multi" if ref.find(rc, p + 1) >= 0 else "rev")
    return (None, len(s), "chimeric")


def assemble(kept, ref):
    """Greedy: take overlaps longest-first, merge if tail->head and no cycle."""
    edges = find_overlaps(kept)
    out = defaultdict(dict)
    for ov, i, j in edges:
        out[i][j] = max(out[i].get(j, 0), ov)
    ambiguous = 0
    for d in out.values():
        m = max(d.values())
        if sum(v == m for v in d.values()) > 1:
            ambiguous += 1

    edges.sort(key=lambda e: (-e[0], e[1], e[2]))
    n = len(kept)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    contigs = {i: kept[i] for i in range(n)}
    places = {i: place(kept[i], ref) for i in range(n)}
    has_out, has_in = [False] * n, [False] * n
    snaps, events = [(dict(places), None)], []
    for ov, i, j in edges:
        if has_out[i] or has_in[j]:
            continue
        ra, rb = find(i), find(j)
        if ra == rb:
            continue
        has_out[i] = has_in[j] = True
        s = contigs.pop(ra) + contigs.pop(rb)[ov:]
        parent[rb] = ra
        contigs[ra] = s
        del places[rb]
        places[ra] = place(s, ref)
        events.append((ov, len(s), places[ra][2]))
        snaps.append((dict(places), ra))
    return dict(n_edges=len(edges), ambiguous=ambiguous,
                contigs=list(contigs.values()), snaps=snaps, events=events)


def snap_stats(pl, L):
    lens, iv, chim = [], [], 0
    for s, e, k in pl.values():
        if s is None:
            lens.append(e)
            chim += 1
        else:
            lens.append(e - s)
            iv.append((s, e))
    lens.sort(reverse=True)
    tot, acc, n50 = sum(lens), 0, 0
    for x in lens:
        acc += x
        if acc * 2 >= tot:
            n50 = x
            break
    iv.sort()
    covered, cs, ce = 0, None, None
    for s, e in iv:
        if cs is None or s > ce:
            if cs is not None:
                covered += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if cs is not None:
        covered += ce - cs
    return dict(n=len(lens), largest=lens[0] if lens else 0, n50=n50,
                frac=covered / L, chim=chim)


def longest_repeat(ref):
    def has_dup(m):
        seen = set()
        for i in range(len(ref) - m + 1):
            w = ref[i:i + m]
            if w in seen:
                return True
            seen.add(w)
        return False
    lo, hi = 0, len(ref)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if has_dup(mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


def zero_ranges(cov):
    out, start = [], None
    for i, c in enumerate(cov + [1]):
        if c == 0 and start is None:
            start = i
        elif c != 0 and start is not None:
            out.append((start, i - 1))
            start = None
    return out


def run_pipeline(base, seed, both, repeat):
    rng = random.Random(seed)
    ref, repeats = base, []
    if repeat:
        seg_s, seg_e = 200, 500
        pos = len(base) * 3 // 4
        ref = base[:pos] + base[seg_s:seg_e] + base[pos:]
        repeats = [(seg_s, seg_e), (pos, pos + seg_e - seg_s)]
    reads, truth = simulate_reads(ref, rng, both)
    cov = coverage(len(ref), truth)
    kept = remove_contained(reads)
    asm = assemble(kept, ref)
    final = asm["snaps"][-1][0]
    st = snap_stats(final, len(ref))
    exact = any(c == ref or c == revcomp(ref) for c in asm["contigs"])
    return dict(ref=ref, repeats=repeats, reads=reads, truth=truth, cov=cov,
                kept=kept, asm=asm, final=st, exact=exact,
                dup=len(reads) - len(set(reads)),
                contained=len(set(reads)) - len(kept),
                rep_len=longest_repeat(ref), zeros=zero_ranges(cov),
                both=both, repeat=repeat, seed=seed)


# ================================================================ report
def build_report(r, header, source):
    L, cov, f = len(r["ref"]), r["cov"], r["final"]
    mean = sum(cov) / L
    zeros = r["zeros"]
    nz = sum(b - a + 1 for a, b in zeros)
    z_txt = ", ".join(f"{a + 1}-{b + 1}" if a != b else f"{a + 1}" for a, b in zeros) or "none"
    R = r["rep_len"]
    out = []
    H1 = lambda t: out.append((t + "\n", "h1"))
    H2 = lambda t: out.append(("\n" + t + "\n", "h2"))
    P = lambda t, tag="body": out.append((t, tag))

    H1("Assembly report")
    P("Original sequence:  ", "dim"); P(header + "\n", "num")
    P("Accession / link:   ", "dim")
    P(f"{ACCESSION}  https://www.ncbi.nlm.nih.gov/nuccore/{ACCESSION}\n", "link")
    P("Source:             ", "dim"); P(source + "\n")
    P("Length:             ", "dim"); P(f"{len(r['ref']):,} bp", "num")
    if r["repeat"]:
        P(f"  (original {L - 300:,} bp + 300 bp synthetic repeat inserted)", "warn")
    P("\n")
    P("Reads:              ", "dim")
    P(f"{N_READS:,} random substrings, {RMIN}-{RMAX} bp"
      f"{', random strand' if r['both'] else ', forward strand only'}, seed {r['seed']}\n")
    P("Coverage:           ", "dim")
    P(f"mean {mean:.0f}x, min {min(cov)}x, uncovered bases: {nz} ({z_txt})\n", "num")
    P("Pre-processing:     ", "dim")
    P(f"{r['dup']:,} exact duplicates, {r['contained']:,} contained reads removed -> "
      f"{len(r['kept']):,} reads in overlap graph, {r['asm']['n_edges']:,} overlap edges (>= {MIN_OV} bp)\n")

    H2("Assembly result")
    P(f"{f['n']} contig(s), largest {f['largest']:,} bp, N50 {f['n50']:,} bp, "
      f"{f['frac'] * 100:.1f}% of the reference recovered, "
      f"{f['chim']} chimeric contig(s), {r['asm']['ambiguous']} ambiguous overlap decision(s).\n", "num")
    if r["exact"]:
        P("Exact reconstruction: YES\n", "good")
    else:
        P("Exact reconstruction: NO", "warn")
        P(" - the reads alone did not determine the sequence uniquely.\n")

    H2("1. Insufficient coverage (sampling)")
    P(f"Mean coverage is about {mean:.0f}x, yet coverage is not uniform. Base i is only covered by "
      "reads that start within the previous ~110 bp, so the first and last bases of a linear sequence "
      "are covered only by reads that start exactly at the edge. With uniform random starts, a few "
      f"end bases are routinely missed (here: {z_txt}). Any base no read covers is simply absent "
      "from the data, and no algorithm can recover it. Lander-Waterman statistics predict gaps "
      "wherever coverage drops, and in real data the drops come from GC bias and hard-to-sequence "
      "regions, not only from the ends.\n")

    H2("2. Ambiguous overlaps")
    P(f"The assembler only sees string matches. A {MIN_OV} bp threshold is a trade-off: shorter "
      "thresholds create spurious random overlaps (about 4^-k chance per position pair), longer ones "
      "drop true joins in low-coverage spots. Each read can have many candidate successors. "
      f"In this run, {r['asm']['ambiguous']} reads had two or more different successors with an equal best overlap. "
      "Greedy chooses one by tie-break order, not by evidence.\n")

    H2("3. Repeats")
    P(f"The longest exact repeat in this reference is {R} bp. ")
    if R < RMIN - MIN_OV:
        P("That is shorter than a read, so reads span every repeat copy plus unique flanks, "
          "and the repeat is resolvable. That is why the plain genome assembles cleanly. ")
    else:
        P("That is LONGER than a read, so reads inside the repeat are identical whichever copy "
          "they came from. ", "warn")
    P("When a repeat exceeds the read length, the reads are consistent with two or more different "
      "genomes (the copies can be swapped, collapsed, or the intervening region excised). This is "
      "an information limit, not an engineering one: the input does not contain the answer. "
      "Toggle 'Insert 300 bp repeat' to watch contigs break at the repeat and turn orange "
      "(ambiguous placement) or red (chimeric).\n")

    H2("4. Chimeric joins")
    P("Greedy merges are locally optimal and irreversible. One wrong merge at a repeat boundary "
      "glues two distant regions into a contig that exists nowhere in the genome, and every later "
      "merge builds on that error. "
      f"This run produced {f['chim']} chimeric contig(s). With noisy real reads, sequencing errors "
      "add false overlaps and hide true ones, so chimeras become more common.\n")

    H2("5. Loss of order and orientation")
    P("A set of reads is an unordered multiset. Positions, strand, and the fact that the genome "
      "is circular (PCV2) are all discarded. The assembler has to infer order from overlaps alone. "
      "With 'Random strand' on, half the reads are reverse complements. This naive assembler "
      "only checks forward overlaps, so it builds two parallel, partly duplicated sets of contigs "
      "(purple = reverse orientation). Handling both strands doubles the graph and adds "
      "palindrome ambiguity. The circular junction is never closed, because no read crosses "
      "the arbitrary linear cut.\n")

    H2("Core algorithmic limitation")
    P("Reconstruction from reads is a shortest-common-superstring / overlap-graph path problem. "
      "It is NP-hard, and greedy is only an approximation. More fundamentally, the most compressed "
      "superstring is not the true genome: repeats get collapsed. The sequence is uniquely "
      "determined ONLY IF every base is covered, every adjacent pair of reads overlaps above the "
      "threshold, and every repeat is spanned by at least one read. Whenever coverage gaps or "
      "repeats longer than the reads exist, many different sequences explain the same reads "
      "equally well. Perfect reconstruction is then impossible for any algorithm. Only longer "
      "reads, paired ends, or extra positional information can fix this.\n", "core")
    return out


# ================================================================ GUI
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Greedy DNA assembly - live")
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        w, h = min(1420, sw - 40), min(900, sh - 90)
        self.geometry(f"{w}x{h}+{max(0, (sw - w) // 2)}+10")
        self.minsize(900, 560)
        self.bind("<space>", lambda e: self.toggle())
        self.configure(bg=BG)
        self.header, self.base, self.source = EMBEDDED_HEADER, EMBEDDED_SEQ, "embedded copy"
        self.res, self.step, self.playing = None, 0, False
        self.speed = tk.DoubleVar(value=25)
        self.both, self.rep = tk.BooleanVar(), tk.BooleanVar()
        self.seed = tk.IntVar(value=42)
        self.acc = tk.StringVar(value=ACCESSION)
        self._redraw_job = None
        self._style()
        self._build()
        self.after(60, self.fetch)

    # ---------------------------------------------------------- styling
    def _style(self):
        s = ttk.Style(self)
        s.theme_use("clam")
        s.configure(".", background=BG, foreground=FG, font=UI)
        s.configure("TFrame", background=BG)
        s.configure("Panel.TFrame", background=PANEL)
        s.configure("TLabel", background=BG, foreground=FG)
        s.configure("Panel.TLabel", background=PANEL, foreground=FG)
        s.configure("Dim.TLabel", background=PANEL, foreground=MUTED)
        s.configure("Big.TLabel", background=PANEL, foreground=ACC, font=("Segoe UI Semibold", 15))
        s.configure("Title.TLabel", background=BG, foreground=FG, font=("Segoe UI Semibold", 16))
        s.configure("Sub.TLabel", background=BG, foreground=MUTED)
        s.configure("TButton", background=PANEL2, foreground=FG, borderwidth=0,
                    focusthickness=0, padding=(12, 6), font=UIB)
        s.map("TButton", background=[("active", GRID), ("pressed", ACC)])
        s.configure("Accent.TButton", background=ACC, foreground="#06121d")
        s.map("Accent.TButton", background=[("active", "#7bdaf5")])
        s.configure("TCheckbutton", background=PANEL, foreground=FG)
        s.map("TCheckbutton", background=[("active", PANEL)])
        s.configure("TEntry", fieldbackground=PANEL2, foreground=FG, insertcolor=FG, borderwidth=0)
        s.configure("TSpinbox", fieldbackground=PANEL2, foreground=FG, arrowcolor=FG,
                    background=PANEL2, borderwidth=0)
        s.configure("TNotebook", background=BG, borderwidth=0)
        s.configure("TNotebook.Tab", background=PANEL, foreground=MUTED, padding=(16, 6), font=UIB)
        s.map("TNotebook.Tab", background=[("selected", PANEL2)], foreground=[("selected", FG)])

    def _scale(self, parent, **kw):
        return tk.Scale(parent, orient="horizontal", bg=ACC, fg=FG, troughcolor=PANEL2,
                        activebackground="#7bdaf5", highlightthickness=0, bd=0, sliderrelief="flat",
                        sliderlength=18, width=12, font=UI, **kw)

    # ---------------------------------------------------------- layout
    def _build(self):
        top = ttk.Frame(self, padding=(16, 12, 16, 6))
        top.pack(fill="x")
        left = ttk.Frame(top)
        left.pack(side="left", fill="x", expand=True)
        ttk.Label(left, text="De novo assembly: greedy overlap-layout-consensus",
                  style="Title.TLabel").pack(anchor="w")
        self.info = ttk.Label(left, text="Fetching sequence from NCBI...", style="Sub.TLabel")
        self.info.pack(anchor="w")
        self.link = tk.Label(left, text="", bg=BG, fg=ACC, cursor="hand2", font=UI)
        self.link.pack(anchor="w")
        self.link.bind("<Button-1>", lambda e: webbrowser.open(
            f"https://www.ncbi.nlm.nih.gov/nuccore/{self.acc.get().strip()}"))
        right = ttk.Frame(top)
        right.pack(side="right")
        ttk.Label(right, text="Accession", style="Sub.TLabel").pack(side="left", padx=(0, 6))
        ttk.Entry(right, textvariable=self.acc, width=14).pack(side="left", ipady=4)
        ttk.Button(right, text="Fetch", command=self.fetch).pack(side="left", padx=6)

        # controls are packed first at the bottom so they are never pushed off-screen
        self.ctrl = ttk.Frame(self)
        self.ctrl.pack(side="bottom", fill="x")
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=16)
        live = ttk.Frame(nb)
        nb.add(live, text="Live assembly")
        rep = ttk.Frame(nb)
        nb.add(rep, text="Analysis")

        self.cv = tk.Canvas(live, bg=PANEL, highlightthickness=0)
        self.cv.pack(side="left", fill="both", expand=True, pady=(8, 0))
        self.cv.bind("<Configure>", self._on_resize)

        side = ttk.Frame(live, style="Panel.TFrame", padding=14, width=330)
        side.pack(side="right", fill="y", padx=(10, 0), pady=(8, 0))
        side.pack_propagate(False)
        self.stat = {}
        grid = ttk.Frame(side, style="Panel.TFrame")
        grid.pack(fill="x")
        for n, (key, label) in enumerate([("step", "Merge step"), ("n", "Contigs"),
                                          ("largest", "Largest contig"), ("n50", "N50"),
                                          ("frac", "Recovered"), ("chim", "Chimeric")]):
            cell = ttk.Frame(grid, style="Panel.TFrame")
            cell.grid(row=n // 2, column=n % 2, sticky="w", padx=(0, 18), pady=(0, 6))
            ttk.Label(cell, text=label.upper(), style="Dim.TLabel").pack(anchor="w")
            v = ttk.Label(cell, text="-", style="Big.TLabel")
            v.pack(anchor="w")
            self.stat[key] = v
        ttk.Label(side, text="RECENT MERGES", style="Dim.TLabel").pack(anchor="w", pady=(6, 2))
        self.log = tk.Text(side, bg=PANEL2, fg=FG, font=("Consolas", 9), bd=0, height=14,
                           highlightthickness=0, padx=8, pady=6, wrap="none")
        self.log.pack(fill="both", expand=True)
        for k, c in KIND_COL.items():
            self.log.tag_configure(k, foreground=c)

        self.report_txt = tk.Text(rep, bg=PANEL, fg=FG, font=("Segoe UI", 11), bd=0, wrap="word",
                           highlightthickness=0, padx=28, pady=20, spacing1=2, spacing3=4)
        sb = ttk.Scrollbar(rep, command=self.report_txt.yview)
        self.report_txt.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y", pady=(8, 0))
        self.report_txt.pack(fill="both", expand=True, pady=(8, 0))
        for tag, kw in {"h1": dict(font=("Segoe UI Semibold", 18), foreground=ACC, spacing3=10),
                        "h2": dict(font=("Segoe UI Semibold", 13), foreground=FG, spacing1=10),
                        "dim": dict(foreground=MUTED, font=MONO), "num": dict(foreground=YEL),
                        "link": dict(foreground=ACC, underline=True), "warn": dict(foreground=RED),
                        "good": dict(foreground=GREEN), "body": dict(foreground="#c9cfdc"),
                        "core": dict(foreground=FG, background=PANEL2, lmargin1=12, lmargin2=12,
                                     rmargin=12, font=("Segoe UI", 11))}.items():
            self.report_txt.tag_configure(tag, **kw)

        # ---- Sequence tab: original vs reassembled, letter by letter
        seqtab = ttk.Frame(nb)
        nb.add(seqtab, text="Sequence")
        tools = ttk.Frame(seqtab)
        tools.pack(fill="x", pady=(8, 4))
        ttk.Button(tools, text="Copy original (FASTA)", command=lambda: self._copy("ref")).pack(side="left")
        ttk.Button(tools, text="Copy assembly (FASTA)", command=lambda: self._copy("asm")).pack(side="left", padx=6)
        self.copy_lbl = ttk.Label(tools, text="", style="Sub.TLabel")
        self.copy_lbl.pack(side="left", padx=10)
        self.seq_txt = self._mono_text(seqtab)

        # ---- Reads tab: the Python list itself
        readtab = ttk.Frame(nb)
        nb.add(readtab, text="Reads (list)")
        self.reads_txt = self._mono_text(readtab)

        # controls
        bar = ttk.Frame(self.ctrl, style="Panel.TFrame", padding=(14, 10))
        bar.pack(fill="x", padx=16, pady=12)
        self.bar_top = bar
        self.play_btn = ttk.Button(bar, text="▶  Play", style="Accent.TButton",
                                   command=self.toggle, width=9)
        self.play_btn.pack(side="left")
        ttk.Button(bar, text="Step ›", command=lambda: self.set_step(self.step + 1)).pack(side="left", padx=4)
        ttk.Button(bar, text="⏮", width=3, command=lambda: self.set_step(0)).pack(side="left")
        ttk.Button(bar, text="⏭ End", command=self.to_end).pack(side="left", padx=4)

        ttk.Label(bar, text="Speed", style="Dim.TLabel").pack(side="left", padx=(16, 4))
        self._scale(bar, from_=1, to=300, variable=self.speed, length=140,
                    showvalue=0).pack(side="left")
        self.speed_lbl = ttk.Label(bar, text="", style="Panel.TLabel", width=11)
        self.speed_lbl.pack(side="left")
        self.speed.trace_add("write", lambda *a: self.speed_lbl.config(
            text=f"{self.speed.get():.0f} merges/s"))
        self.speed.set(25)

        ttk.Label(bar, text="Timeline", style="Dim.TLabel").pack(side="left", padx=(12, 4))
        self.tl = self._scale(bar, from_=0, to=1, resolution=1, showvalue=0, command=self.on_scrub)
        self.tl.pack(side="left", fill="x", expand=True)
        self.tl_lbl = ttk.Label(bar, text="0 / 0", style="Panel.TLabel", width=10)
        self.tl_lbl.pack(side="left", padx=(6, 0))

        bar = ttk.Frame(self.ctrl, style="Panel.TFrame", padding=(14, 0, 14, 10))
        bar.pack(fill="x", padx=16, pady=(0, 12))
        self.bar_top.pack_configure(pady=(12, 0))
        ttk.Label(bar, text="Scenario", style="Dim.TLabel").pack(side="left", padx=(0, 8))
        ttk.Checkbutton(bar, text="Random strand", variable=self.both,
                        command=self.rerun).pack(side="left", padx=4)
        ttk.Checkbutton(bar, text="Insert 300 bp repeat", variable=self.rep,
                        command=self.rerun).pack(side="left", padx=4)
        ttk.Label(bar, text="Seed", style="Dim.TLabel").pack(side="left", padx=(10, 4))
        sp = ttk.Spinbox(bar, from_=0, to=99999, textvariable=self.seed, width=6, command=self.rerun)
        sp.pack(side="left")
        sp.bind("<Return>", lambda e: self.rerun())
        sp.bind("<FocusOut>", lambda e: self.rerun() if self.res and self.res["seed"] != self._seed_val() else None)
        ttk.Button(bar, text="↻ New reads", command=self.new_seed).pack(side="left", padx=(6, 0))

    def _mono_text(self, parent):
        wrap = ttk.Frame(parent)
        wrap.pack(fill="both", expand=True)
        t = tk.Text(wrap, bg=PANEL, fg=FG, font=("Consolas", 11), bd=0, wrap="none",
                    highlightthickness=0, padx=20, pady=14, insertbackground=FG)
        vs = ttk.Scrollbar(wrap, command=t.yview)
        hs = ttk.Scrollbar(wrap, orient="horizontal", command=t.xview)
        t.configure(yscrollcommand=vs.set, xscrollcommand=hs.set)
        vs.pack(side="right", fill="y")
        hs.pack(side="bottom", fill="x")
        t.pack(fill="both", expand=True)
        for b, c in BASE_COL.items():
            t.tag_configure(b, foreground=c)
        t.tag_configure("pos", foreground=MUTED)
        t.tag_configure("dim", foreground=MUTED)
        t.tag_configure("code", foreground="#c9cfdc")
        t.tag_configure("kw", foreground=PURPLE)
        t.tag_configure("str", foreground=GREEN)
        t.tag_configure("num", foreground=YEL)
        t.tag_configure("h", foreground=ACC, font=("Segoe UI Semibold", 13))
        t.tag_configure("miss", foreground=RED, background="#3a1622")
        t.tag_configure("mism", foreground="#ffffff", background=RED)
        t.tag_configure("chim", foreground=RED)
        return t

    def _fill(self, widget, parts):
        widget.config(state="normal")
        widget.delete("1.0", "end")
        buf, cur = [], None
        for txt, tag in parts:          # merge runs of the same tag -> fast insert
            if tag != cur and buf:
                widget.insert("end", "".join(buf), cur or ())
                buf = []
            cur = tag
            buf.append(txt)
        if buf:
            widget.insert("end", "".join(buf), cur or ())
        widget.config(state="disabled")

    def _write_sequence(self):
        r = self.res
        ref, L = r["ref"], len(r["ref"])
        asm = [None] * L
        contigs = sorted(r["asm"]["contigs"], key=len, reverse=True)
        chim = []
        for c in contigs:
            st, _, k = place(c, ref)
            if st is None:
                chim.append(c)
                continue
            seq = c if ref.startswith(c, st) else revcomp(c)
            for i, ch in enumerate(seq):
                if asm[st + i] is None:
                    asm[st + i] = ch
        n_miss = asm.count(None)
        n_mm = sum(1 for a, b in zip(asm, ref) if a is not None and a != b)
        p = []
        p.append(("Original vs reassembled sequence\n", "h"))
        p.append((f"{self.header}  ·  {L:,} bp\n", "dim"))
        p.append((f"Assembled into {len(contigs)} contig(s)  ·  missing bases: {n_miss}  ·  "
                  f"wrong bases: {n_mm}  ·  chimeric contigs: {len(chim)}\n", "num"))
        p.append(("ref = original from NCBI    asm = what the assembler rebuilt    ", "dim"))
        p.append(("-", "miss")); p.append((" = not recovered   ", "dim"))
        p.append(("X", "mism")); p.append((" = wrong base\n\n", "dim"))
        for start in range(0, L, 60):
            end = min(L, start + 60)
            p.append((f"ref {start + 1:>6}  ", "pos"))
            for i in range(start, end):
                p.append((ref[i], ref[i] if ref[i] in BASE_COL else "dim"))
                if (i - start) % 10 == 9:
                    p.append((" ", "dim"))
            p.append((f" {end}\n", "pos"))
            p.append(("asm         ", "pos"))
            for i in range(start, end):
                a = asm[i]
                if a is None:
                    p.append(("-", "miss"))
                elif a != ref[i]:
                    p.append((a, "mism"))
                else:
                    p.append((a, a if a in BASE_COL else "dim"))
                if (i - start) % 10 == 9:
                    p.append((" ", "dim"))
            p.append(("\n\n", "dim"))
        if chim:
            p.append(("Chimeric contigs (exist nowhere in the original)\n", "h"))
            for n, c in enumerate(chim, 1):
                p.append((f">chimeric_{n}  {len(c):,} bp\n", "chim"))
                for i in range(0, len(c), 60):
                    p.append((c[i:i + 60] + "\n", "code"))
        self._fill(self.seq_txt, p)

    def _write_reads(self):
        r = self.res
        reads, truth = r["reads"], r["truth"]
        lens = [len(x) for x in reads]
        p = [("The reads are stored in an ordinary Python list\n", "h"),
             ("reads", "code"), (" = simulate_reads(ref, rng)", "dim"),
             ("        # list of str\n\n", "dim"),
             (">>> ", "dim"), ("type", "kw"), ("(reads)\n", "code"),
             (f"{type(reads)}\n", "str"),
             (">>> ", "dim"), ("len", "kw"), ("(reads)\n", "code"), (f"{len(reads)}\n", "num"),
             (">>> ", "dim"), ("min", "kw"), ("(map(len, reads)), ", "code"),
             ("max", "kw"), ("(map(len, reads))\n", "code"),
             (f"({min(lens)}, {max(lens)})\n", "num"),
             (">>> ", "dim"), ("reads[0]\n", "code"), (f"{reads[0]!r}\n", "str"),
             (">>> ", "dim"), ("reads[:3]\n", "code"),
             ("[" + ",\n ".join(repr(x) for x in reads[:3]) + ", ...]\n\n", "str"),
             (f"All {len(reads):,} entries   (start/strand shown as a comment = the hidden truth, "
              "the assembler never sees it)\n\n", "dim"),
             ("reads = [\n", "code")]
        for i, (x, (st, ln, strand)) in enumerate(zip(reads, truth)):
            p.append(("    ", "pos"))
            p.append((f"{x!r}", "str"))
            p.append((f",   # reads[{i}]  len={ln}  true start={st + 1}  strand={strand}\n", "dim"))
        p.append(("]\n", "code"))
        self._fill(self.reads_txt, p)

    def _copy(self, which):
        r = self.res
        if not r:
            return
        wrap = lambda s: "\n".join(s[i:i + 70] for i in range(0, len(s), 70))
        if which == "ref":
            txt = f">{self.header}\n{wrap(r['ref'])}\n"
        else:
            cs = sorted(r["asm"]["contigs"], key=len, reverse=True)
            txt = "".join(f">contig_{n} length={len(c)}\n{wrap(c)}\n" for n, c in enumerate(cs, 1))
        self.clipboard_clear()
        self.clipboard_append(txt)
        self.copy_lbl.config(text="Copied to clipboard")
        self.after(2000, lambda: self.copy_lbl.config(text=""))

    # ---------------------------------------------------------- data
    def fetch(self):
        acc = self.acc.get().strip()
        self.info.config(text=f"Fetching {acc} from NCBI…")
        self._fetched = None

        def work():
            try:
                h, s = fetch_ncbi(acc)
                if not 1000 <= len(s) <= 3000:
                    raise ValueError(f"{acc} is {len(s):,} bp (need 1,000–3,000)")
                self._fetched = ("ok", h, s)
            except Exception as e:
                self._fetched = ("err", str(e), None)
        threading.Thread(target=work, daemon=True).start()
        self._poll_fetch(acc)

    def _poll_fetch(self, acc):
        if self._fetched is None:
            self.after(100, self._poll_fetch, acc)
            return
        st, a, b = self._fetched
        if st == "ok":
            self.header, self.base, self.source = a, b, "NCBI E-utilities (live fetch)"
        elif acc == ACCESSION:
            self.header, self.base = EMBEDDED_HEADER, EMBEDDED_SEQ
            self.source = "embedded copy (NCBI unreachable)"
        else:
            self.info.config(text=f"Fetch failed: {a}")
            return
        self.rerun()

    def _seed_val(self):
        try:
            return int(self.seed.get())
        except (tk.TclError, ValueError):
            return 42

    def new_seed(self):
        self.seed.set(random.randint(0, 99999))
        self.rerun()

    def rerun(self):
        self.playing = False
        self.play_btn.config(text="▶  Play")
        self.config(cursor="watch")
        self.update_idletasks()
        seed = self._seed_val()
        self.res = run_pipeline(self.base, seed, self.both.get(), self.rep.get())
        self.config(cursor="")
        r = self.res
        self.info.config(text=f"{self.header}  ·  {len(self.base):,} bp  ·  {self.source}")
        self.link.config(text=f"https://www.ncbi.nlm.nih.gov/nuccore/{self.acc.get().strip()}")
        self.tl.config(to=len(r["asm"]["events"]))
        self.step = -1
        self._write_report()
        self._write_sequence()
        self._write_reads()
        print(f"[run] reads is a {type(r['reads']).__name__} of {len(r['reads'])} strings; "
              f"reads[0] = {r['reads'][0]!r}")
        self.draw_static()
        self.set_step(0)

    def _write_report(self):
        self.report_txt.config(state="normal")
        self.report_txt.delete("1.0", "end")
        for txt, tag in build_report(self.res, self.header, self.source):
            self.report_txt.insert("end", txt, tag)
        self.report_txt.config(state="disabled")

    # ---------------------------------------------------------- playback
    @property
    def n_steps(self):
        return len(self.res["asm"]["events"]) if self.res else 0

    def set_step(self, k):
        if not self.res:
            return
        k = max(0, min(self.n_steps, int(k)))
        if k == self.step:
            return
        self.step = k
        if int(float(self.tl.get())) != k:
            self.tl.set(k)
        self.draw_dynamic()

    def on_scrub(self, v):
        self.set_step(int(float(v)))

    def to_end(self):
        self.playing = False
        self.play_btn.config(text="▶  Play")
        self.set_step(self.n_steps)

    def toggle(self):
        if not self.res:
            return
        self.playing = not self.playing
        if self.playing and self.step >= self.n_steps:
            self.set_step(0)
        self.play_btn.config(text="❚❚  Pause" if self.playing else "▶  Play")
        if self.playing:
            self.tick()

    def tick(self):
        if not self.playing:
            return
        sp = max(1.0, self.speed.get())
        delay = max(16, int(1000 / sp))
        self.set_step(self.step + max(1, round(sp * delay / 1000)))
        if self.step >= self.n_steps:
            self.playing = False
            self.play_btn.config(text="▶  Play")
        else:
            self.after(delay, self.tick)

    # ---------------------------------------------------------- drawing
    def _on_resize(self, _):
        if self._redraw_job:
            self.after_cancel(self._redraw_job)
        self._redraw_job = self.after(80, self._full_redraw)

    def _full_redraw(self):
        if self.res:
            self.draw_static()
            self.draw_dynamic()

    def _geom(self):
        W, H = max(self.cv.winfo_width(), 600), max(self.cv.winfo_height(), 420)
        X0, X1 = 70, W - 24
        L = len(self.res["ref"])
        return W, H, X0, X1, L, (lambda p: X0 + p / L * (X1 - X0))

    def draw_static(self):
        cv, r = self.cv, self.res
        cv.delete("all")
        W, H, X0, X1, L, X = self._geom()
        # title + ruler
        cv.create_text(X0, 18, anchor="w", fill=FG, font=UIB,
                       text=f"Reference  {self.acc.get().strip()}  ·  {L:,} bp")
        cv.create_line(X0, 44, X1, 44, fill=MUTED, width=2)
        step = 250 if L < 2500 else 500
        for p in range(0, L + 1, step):
            cv.create_line(X(p), 40, X(p), 48, fill=MUTED)
            cv.create_text(X(p), 32, text=f"{p:,}", fill=MUTED, font=("Segoe UI", 8))
        for a, b in r["repeats"]:
            cv.create_rectangle(X(a), 40, X(b), 48, fill=ORANGE, outline="")
        # coverage
        y0, y1 = 64, 144
        cov, mx = r["cov"], max(r["cov"])
        cv.create_text(12, (y0 + y1) / 2, anchor="w", fill=MUTED, font=("Segoe UI", 8),
                       text="coverage")
        cv.create_text(X0 - 6, y0, anchor="e", fill=MUTED, font=("Segoe UI", 8), text=f"{mx}x")
        cv.create_text(X0 - 6, y1, anchor="e", fill=MUTED, font=("Segoe UI", 8), text="0")
        pts = [X0, y1]
        cols = int(X1 - X0)
        for c in range(cols + 1):
            p = min(L - 1, int(c / cols * L))
            pts += [X0 + c, y1 - cov[p] / mx * (y1 - y0)]
        pts += [X1, y1]
        cv.create_polygon(pts, fill="#17465c", outline=ACC)
        for a, b in r["zeros"]:
            cv.create_rectangle(X(a) - 1, y0, X(b + 1) + 1, y1 + 4, fill=RED, outline="")
            cv.create_text(X(a), y1 + 12, text="gap", fill=RED, font=("Segoe UI", 8))
        # reads
        ry0 = 166
        band = 100 if H >= 640 else 40
        cv.create_text(12, ry0 + band / 2, anchor="w", fill=MUTED, font=("Segoe UI", 8),
                       text=f"{N_READS:,}\nreads")
        for i, (s, ln, strand) in enumerate(r["truth"]):
            y = ry0 + (i * 7919 % (band // 2)) * 2
            cv.create_line(X(s), y, X(s + ln), y, fill="#2f4466" if strand == "+" else "#4b3670")
        cv.create_text(X1, ry0 - 8, anchor="e", fill=MUTED, font=("Segoe UI", 8),
                       text="true read positions (hidden from the assembler)")
        cv.create_line(X0, ry0 + band + 16, X1, ry0 + band + 16, fill=GRID)
        self._ctop = ry0 + band + 46

    def draw_dynamic(self):
        cv, r = self.cv, self.res
        cv.delete("dyn")
        W, H, X0, X1, L, X = self._geom()
        pl, hi = r["asm"]["snaps"][self.step]
        if self.step == self.n_steps:
            hi = None
        top = getattr(self, "_ctop", 312)
        bot = H - 78
        ly = top - 14
        # legend
        cv.create_text(X0, ly, anchor="w", fill=FG, font=UIB, tags="dyn",
                       text=f"Contigs at step {self.step:,} / {self.n_steps:,}")
        lx = X1
        for k in ("chimeric", "multi", "rev", "fwd"):
            t = cv.create_text(lx, ly, anchor="e", fill=MUTED, font=("Segoe UI", 8),
                               text=KIND_TXT[k], tags="dyn")
            bx = cv.bbox(t)[0]
            cv.create_rectangle(bx - 14, ly - 5, bx - 4, ly + 5, fill=KIND_COL[k], outline="", tags="dyn")
            lx = bx - 24
        # pack placed contigs into rows
        placed = sorted(((s, e, k, cid) for cid, (s, e, k) in pl.items() if s is not None))
        rows_end, rows = [], []
        for s, e, k, cid in placed:
            xs, xe = X(s), X(e)
            for ri, re_ in enumerate(rows_end):
                if xs > re_ + 2:
                    rows_end[ri] = xe
                    rows[ri].append((xs, xe, k, cid))
                    break
            else:
                rows_end.append(xe)
                rows.append([(xs, xe, k, cid)])
        n_rows = max(1, len(rows))
        rh = min(14.0, (bot - top) / n_rows)
        bh = max(1.0, rh - (3 if rh > 6 else 1 if rh > 2 else 0))
        for ri, row in enumerate(rows):
            y = top + ri * rh
            for xs, xe, k, cid in row:
                is_hi = cid == hi
                cv.create_rectangle(xs, y, max(xe, xs + 1), y + bh, fill=YEL if is_hi else KIND_COL[k],
                                    outline=YEL if is_hi else "", width=2 if is_hi else 0, tags="dyn")
        # chimeric contigs (no place in the genome)
        chim = [(e, cid) for cid, (s, e, k) in pl.items() if s is None]
        yc = H - 56
        cv.create_text(X0, yc - 12, anchor="w", fill=RED, font=("Segoe UI", 8), tags="dyn",
                       text=f"unplaceable / chimeric contigs: {len(chim)}")
        x = X0
        for ln, cid in sorted(chim, reverse=True):
            w = ln / L * (X1 - X0)
            cv.create_rectangle(x, yc, min(X1, x + w), yc + 10,
                                fill=YEL if cid == hi else RED, outline="", tags="dyn")
            x += w + 4
            if x > X1:
                break
        if self.step == self.n_steps:
            f = r["final"]
            msg = ("✔ EXACT reconstruction" if r["exact"] else "✘ not reconstructed exactly") + \
                  f"  ·  {f['n']} contig(s)  ·  largest {f['largest']:,} bp  ·  {f['frac'] * 100:.1f}% recovered"
            cv.create_text(X1, H - 20, anchor="e", fill=GREEN if r["exact"] else RED,
                           font=UIB, text=msg, tags="dyn")
        self._update_side(pl)

    def _update_side(self, pl):
        st = snap_stats(pl, len(self.res["ref"]))
        self.stat["step"].config(text=f"{self.step:,} / {self.n_steps:,}")
        self.stat["n"].config(text=f"{st['n']:,}")
        self.stat["largest"].config(text=f"{st['largest']:,} bp")
        self.stat["n50"].config(text=f"{st['n50']:,} bp")
        self.stat["frac"].config(text=f"{st['frac'] * 100:.1f} %")
        self.stat["chim"].config(text=str(st["chim"]),
                                 foreground=RED if st["chim"] else ACC)
        self.tl_lbl.config(text=f"{self.step} / {self.n_steps}")
        self.log.config(state="normal")
        self.log.delete("1.0", "end")
        ev = self.res["asm"]["events"]
        for idx in range(self.step, max(0, self.step - 40), -1):
            ov, ln, k = ev[idx - 1]
            self.log.insert("end", f"#{idx:<4} ov {ov:>3}  → {ln:>5,} bp  ", "")
            self.log.insert("end", f"{KIND_TXT[k] if k != 'fwd' else 'ok'}\n", k)
        if self.step == 0:
            self.log.insert("end", "Unordered reads; no merges yet.\nPress Play.\n")
        self.log.config(state="disabled")


if __name__ == "__main__":
    App().mainloop()
