"""Needleman-Wunsch global alignment of two DNA sequences - live visual Tkinter app."""
import sys
import tkinter as tk
from tkinter import font as tkfont

S1_DEFAULT = "ACCGTGAAGCCAATAC"
S2_DEFAULT = "AGCGTGCAGCCAATAC"


# ---------------------------------------------------------------- algorithm
def score_matrix(s1, s2, match=1, mismatch=-1, gap=-1):
    n, m = len(s1), len(s2)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = F[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            F[i][j] = max(diag, F[i - 1][j] + gap, F[i][j - 1] + gap)
    return F


def traceback(F, s1, s2, match=1, mismatch=-1, gap=-1):
    """Return (aligned1, aligned2, path) where path goes from (n, m) to (0, 0)."""
    i, j = len(s1), len(s2)
    a1, a2, path = [], [], [(i, j)]
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + (
                match if s1[i - 1] == s2[j - 1] else mismatch):
            a1.append(s1[i - 1]); a2.append(s2[j - 1]); i -= 1; j -= 1
        elif i > 0 and F[i][j] == F[i - 1][j] + gap:
            a1.append(s1[i - 1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(s2[j - 1]); j -= 1
        path.append((i, j))
    return "".join(reversed(a1)), "".join(reversed(a2)), path


def markers(a1, a2):
    return "".join("|" if x == y else (" " if "-" in (x, y) else ".")
                   for x, y in zip(a1, a2))


def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-1):
    """Return (score_matrix, aligned1, markers, aligned2, score)."""
    F = score_matrix(s1, s2, match, mismatch, gap)
    a1, a2, _ = traceback(F, s1, s2, match, mismatch, gap)
    return F, a1, markers(a1, a2), a2, F[len(s1)][len(s2)]


# ---------------------------------------------------------------- theme
BG, CARD, CARD2, LINE = "#0e1117", "#161b26", "#1f2636", "#2a3245"
TEXT, MUTED = "#e6e9ef", "#8b93a7"
ACCENT, ACCENT_HI, CYAN, PINK = "#7c5cff", "#9a80ff", "#22d3ee", "#f472b6"
GOOD, WARN, BAD = "#34d399", "#fbbf24", "#f87171"
EMPTY = "#1a2030"
BASE = {"A": "#34d399", "C": "#60a5fa", "G": "#fbbf24", "T": "#f87171", "-": "#4b5563"}
HEAT = [(32, 36, 60), (96, 72, 210), (34, 211, 238)]  # low -> mid -> high
UI, MONO = "TkDefaultFont", "TkFixedFont"


def init_fonts(root):
    global UI, MONO
    fams = set(tkfont.families(root))
    UI = next((f for f in ("Segoe UI", "Inter", "Helvetica Neue", "DejaVu Sans") if f in fams), UI)
    MONO = next((f for f in ("Cascadia Mono", "JetBrains Mono", "Consolas", "DejaVu Sans Mono")
                 if f in fams), MONO)


def heat(t):
    t = min(1.0, max(0.0, t))
    seg = t * (len(HEAT) - 1)
    k = min(int(seg), len(HEAT) - 2)
    u = seg - k
    a, b = HEAT[k], HEAT[k + 1]
    rgb = [int(a[c] + (b[c] - a[c]) * u) for c in range(3)]
    lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    return "#%02x%02x%02x" % tuple(rgb), ("#0e1117" if lum > 150 else TEXT)


def dark_titlebar(win):
    if sys.platform != "win32":
        return
    try:
        import ctypes
        win.update_idletasks()
        hwnd = ctypes.windll.user32.GetParent(win.winfo_id())
        val = ctypes.c_int(1)
        ctypes.windll.dwmapi.DwmSetWindowAttribute(hwnd, 20, ctypes.byref(val), ctypes.sizeof(val))
    except Exception:
        pass


def round_rect(cv, x0, y0, x1, y1, r, **kw):
    pts = [x0 + r, y0, x1 - r, y0, x1, y0, x1, y0 + r, x1, y1 - r, x1, y1,
           x1 - r, y1, x0 + r, y1, x0, y1, x0, y1 - r, x0, y0 + r, x0, y0]
    return cv.create_polygon(pts, smooth=True, **kw)


class Pill(tk.Label):
    """Flat rounded-looking button."""

    def __init__(self, master, text, command, primary=False, small=False):
        self.base = ACCENT if primary else CARD2
        self.hover = ACCENT_HI if primary else LINE
        super().__init__(master, text=text, bg=self.base, fg=TEXT,
                         font=(UI, 9 if small else 10, "bold"),
                         padx=8 if small else 14, pady=2 if small else 6, cursor="hand2")
        self.command, self.enabled = command, True
        self.bind("<Enter>", lambda e: self.enabled and self.config(bg=self.hover))
        self.bind("<Leave>", lambda e: self.config(bg=self.base))
        self.bind("<Button-1>", lambda e: self.enabled and self.command())

    def set_enabled(self, on):
        self.enabled = on
        self.config(fg=TEXT if on else MUTED, cursor="hand2" if on else "arrow")


def card(master, **kw):
    return tk.Frame(master, bg=CARD, highlightthickness=1, highlightbackground=LINE, **kw)


# ---------------------------------------------------------------- alignment view
class AlignView(tk.Frame):
    """Heatmap matrix + traceback + tile alignment + stat cards.
    animate=True  -> fills cell by cell with playback controls.
    animate=False -> draws the finished result instantly (live mode)."""

    def __init__(self, master, s1, s2, match, mismatch, gap, animate=True, max_px=520):
        super().__init__(master, bg=BG)
        self.animate = animate
        self.s1, self.s2 = s1, s2
        self.match, self.mismatch, self.gap = match, mismatch, gap
        self.F = score_matrix(s1, s2, match, mismatch, gap)
        self.a1, self.a2, self.path = traceback(self.F, s1, s2, match, mismatch, gap)
        self.mk = markers(self.a1, self.a2)
        vals = [v for row in self.F for v in row]
        self.lo, self.hi = min(vals), max(vals)
        n, m = len(s1), len(s2)
        self.cs = max(8, min(42, max_px // (max(n, m) + 2)))
        self.base_cs = self.cs
        self.min_cs, self.max_cs = max(6, self.cs // 2), 140

        self.events = [("fill", i, j) for i in range(n + 1) for j in range(m + 1)]
        self.events += [("trace", k) for k in range(len(self.path))]
        self.events.append(("done",))
        self.pos, self.running, self.job = 0, True, None

        self._build()
        self._draw_grid()
        self._draw_alignment(0)
        if animate:
            self._tick()
        else:
            self.skip()

    # ---------- layout
    def _build(self):
        n, m = len(self.s1), len(self.s2)
        self.W, self.H = (m + 2) * self.cs, (n + 2) * self.cs

        if self.animate:
            bar = tk.Frame(self, bg=BG)
            bar.pack(fill="x", pady=(0, 8))
            self.play_btn = Pill(bar, "Pause", self.toggle, primary=True)
            self.play_btn.pack(side="left")
            for txt, cmd in (("Step", self.step), ("Skip to end", self.skip), ("Restart", self.restart)):
                Pill(bar, txt, cmd).pack(side="left", padx=(6, 0))
            self.speed = tk.DoubleVar(value=40)
            tk.Label(bar, text="speed", bg=BG, fg=MUTED, font=(UI, 9)).pack(side="left", padx=(16, 4))
            tk.Scale(bar, from_=1, to=100, orient="horizontal", variable=self.speed, showvalue=0,
                     bg=ACCENT, troughcolor=CARD2, activebackground=ACCENT_HI, highlightthickness=0,
                     bd=0, sliderrelief="flat", width=8, sliderlength=18, length=150).pack(side="left")

        mc = card(self, padx=12, pady=12)
        mc.pack(fill="both", expand=True)
        hdr = tk.Frame(mc, bg=CARD)
        hdr.pack(fill="x")
        tk.Label(hdr, text="SCORE MATRIX", bg=CARD, fg=MUTED, font=(UI, 8, "bold")).pack(side="left")
        tk.Label(hdr, text="   Ctrl+scroll to zoom \u00b7 drag to pan \u00b7 double-click to fit",
                 bg=CARD, fg="#5b6478", font=(UI, 8)).pack(side="left")
        Pill(hdr, "Fit", self.zoom_fit, small=True).pack(side="right")
        Pill(hdr, "+", lambda: self.zoom(1.25), small=True).pack(side="right", padx=(4, 4))
        self.zoom_lbl = tk.Label(hdr, text="100%", width=5, bg=CARD, fg=TEXT, font=(MONO, 9))
        self.zoom_lbl.pack(side="right")
        Pill(hdr, "\u2212", lambda: self.zoom(0.8), small=True).pack(side="right", padx=(0, 4))
        if self.animate:
            self.prog = tk.Canvas(mc, height=3, bg=CARD2, highlightthickness=0)
            self.prog.pack(fill="x", pady=(6, 0))
        self.cv = tk.Canvas(mc, width=self.W, height=self.H, bg=CARD, highlightthickness=0,
                            xscrollincrement=1, yscrollincrement=1, cursor="fleur")
        self.cv.pack(pady=(8, 0), fill="both", expand=True)
        cv = self.cv
        cv.bind("<Control-MouseWheel>", lambda e: self._wheel_zoom(e, e.delta > 0))
        cv.bind("<Control-Button-4>", lambda e: self._wheel_zoom(e, True))
        cv.bind("<Control-Button-5>", lambda e: self._wheel_zoom(e, False))
        cv.bind("<MouseWheel>", lambda e: (cv.yview_scroll(-e.delta // 3, "units"), "break")[1])
        cv.bind("<Shift-MouseWheel>", lambda e: (cv.xview_scroll(-e.delta // 3, "units"), "break")[1])
        cv.bind("<Button-4>", lambda e: (cv.yview_scroll(-40, "units"), "break")[1])
        cv.bind("<Button-5>", lambda e: (cv.yview_scroll(40, "units"), "break")[1])
        cv.bind("<ButtonPress-1>", lambda e: cv.scan_mark(e.x, e.y))
        cv.bind("<B1-Motion>", lambda e: cv.scan_dragto(e.x, e.y, gain=1))
        cv.bind("<Double-Button-1>", lambda e: self.zoom_fit())

        self.status = tk.Label(self, text="", width=1, bg=BG, fg=MUTED, font=(MONO, 10), anchor="w")
        self.status.pack(fill="x", pady=(8, 8))

        ac = card(self, padx=12, pady=10)
        ac.pack(fill="x")
        tk.Label(ac, text="ALIGNMENT", bg=CARD, fg=MUTED, font=(UI, 8, "bold")).pack(anchor="w")
        L = len(self.a1)
        self.aw = max(self.W, 480)
        self.tw = max(8, min(26, (self.aw - 50) // max(L, 1)))
        self.al = tk.Canvas(ac, width=self.aw, height=self.tw * 3 + 14, bg=CARD, highlightthickness=0)
        self.al.pack(pady=(6, 0), anchor="w")

        sc = tk.Frame(self, bg=BG)
        sc.pack(fill="x", pady=(8, 0))
        self.stat = {}
        for k, name in enumerate(("Score", "Length", "Matches", "Mismatches", "Gaps", "Similarity")):
            c = card(sc, padx=10, pady=6)
            c.grid(row=0, column=k, sticky="nsew", padx=(0 if k == 0 else 6, 0))
            sc.columnconfigure(k, weight=1)
            tk.Label(c, text=name.upper(), bg=CARD, fg=MUTED, font=(UI, 7, "bold")).pack(anchor="w")
            v = tk.Label(c, text="-", bg=CARD, fg=TEXT, font=(UI, 15, "bold"))
            v.pack(anchor="w")
            self.stat[name] = v
        self.simbar = tk.Canvas(self.stat["Similarity"].master, width=90, height=4, bg=CARD2,
                                highlightthickness=0)
        self.simbar.pack(anchor="w", pady=(2, 0))

    def _cell_box(self, i, j, pad=1.5):
        cs = self.cs
        x, y = (j + 1) * cs, (i + 1) * cs
        return x + pad, y + pad, x + cs - pad, y + cs - pad

    def _center(self, i, j):
        return (j + 1.5) * self.cs, (i + 1.5) * self.cs

    def _draw_grid(self):
        cv, cs = self.cv, self.cs
        n, m = len(self.s1), len(self.s2)
        self.show_nums = cs >= 20
        self.r = max(2, cs // 6)
        cv.configure(scrollregion=(0, 0, (m + 2) * cs, (n + 2) * cs))
        hf = (MONO, max(7, cs // 3 + 1), "bold")
        self.nf = (MONO, max(7, cs // 3 - 1))
        for j in range(m + 1):
            ch = "·" if j == 0 else self.s2[j - 1]
            cv.create_text((j + 1.5) * cs, cs / 2, text=ch, font=hf, fill=BASE.get(ch, MUTED),
                           tags="hdr")
        for i in range(n + 1):
            ch = "·" if i == 0 else self.s1[i - 1]
            cv.create_text(cs / 2, (i + 1.5) * cs, text=ch, font=hf, fill=BASE.get(ch, MUTED),
                           tags="hdr")
        self.cells, self.texts = {}, {}
        for i in range(n + 1):
            for j in range(m + 1):
                self.cells[i, j] = round_rect(cv, *self._cell_box(i, j), self.r, fill=EMPTY, outline="")
                cx, cy = self._center(i, j)
                self.texts[i, j] = cv.create_text(cx, cy, font=self.nf, fill=TEXT, tags="txt")

    # ---------- drawing helpers
    def _ring(self, i, j, color, width, tag):
        round_rect(self.cv, *self._cell_box(i, j, pad=1), self.r, fill="", outline=color,
                   width=width, tags=tag)

    def _link(self, a, b, color, width, tag, arrow="last"):
        (x1, y1), (x2, y2) = self._center(*a), self._center(*b)
        dx, dy = x2 - x1, y2 - y1
        d = (dx * dx + dy * dy) ** 0.5 or 1
        f = self.cs * 0.32 / d
        self.cv.create_line(x1 + dx * f, y1 + dy * f, x2 - dx * f, y2 - dy * f, fill=color,
                            width=width, arrow=arrow, capstyle="round", tags=tag)

    def _header_hl(self, i=None, j=None):
        self.cv.delete("hdrhl")
        if i is None:
            return
        cs = self.cs
        round_rect(self.cv, (j + 1) * cs + 2, 2, (j + 2) * cs - 2, cs - 2, self.r,
                   fill=CARD2, outline=ACCENT, tags="hdrhl")
        round_rect(self.cv, 2, (i + 1) * cs + 2, cs - 2, (i + 2) * cs - 2, self.r,
                   fill=CARD2, outline=ACCENT, tags="hdrhl")
        self.cv.tag_raise("hdr")

    # ---------- animation steps
    def _do(self, ev):
        if ev[0] == "fill":
            self._fill(ev[1], ev[2])
        elif ev[0] == "trace":
            self._trace(ev[1])
        else:
            self._finish()

    def _fill(self, i, j, quiet=False):
        F, g, v = self.F, self.gap, self.F[i][j]
        t = 0.5 if self.hi == self.lo else (v - self.lo) / (self.hi - self.lo)
        bg, fg = heat(t)
        self.cv.itemconfig(self.cells[i, j], fill=bg)
        if self.show_nums:
            self.cv.itemconfig(self.texts[i, j], text=str(v), fill=fg)
        if not self.animate or quiet:
            return
        self.cv.delete("hl")
        self._header_hl(i, j)
        self._ring(i, j, "white", 2, "hl")
        if i == 0 and j == 0:
            msg = "start   F[0][0] = 0"
        elif i == 0 or j == 0:
            k = j if i == 0 else i
            msg = f"border  F[{i}][{j}] = {k} × gap({g}) = {v}"
        else:
            a, b = self.s1[i - 1], self.s2[j - 1]
            s = self.match if a == b else self.mismatch
            cand = {"↘ diag": (F[i - 1][j - 1] + s, (i - 1, j - 1)),
                    "↓ up": (F[i - 1][j] + g, (i - 1, j)),
                    "→ left": (F[i][j - 1] + g, (i, j - 1))}
            for val, src in cand.values():
                if val == v:
                    self._ring(*src, CYAN, 2, "hl")
                    self._link(src, (i, j), CYAN, 2, "hl")
            parts = "   ".join(f"{k} {val}" for k, (val, _) in cand.items())
            msg = (f"F[{i}][{j}]  {a}/{b} {'match' if a == b else 'mismatch'} {s:+d}   "
                   f"{parts}   → max {v}")
        self.status.config(text=msg)

    def _trace(self, k, quiet=False):
        if k == 0:
            self.cv.delete("hl")
            self._header_hl()
        i, j = self.path[k]
        self._ring(i, j, PINK, 2, "path")
        if k == 0:
            msg = f"traceback from bottom-right   F[{i}][{j}] = {self.F[i][j]}"
        else:
            pi, pj = self.path[k - 1]
            self._link((pi, pj), (i, j), PINK, 3, "path", arrow="none")
            if pi - i == 1 and pj - j == 1:
                a, b = self.s1[pi - 1], self.s2[pj - 1]
                what = f"↖ diagonal  {a}/{b} {'match' if a == b else 'mismatch'}"
            elif pi - i == 1:
                what = f"↑ up  {self.s1[pi - 1]} vs gap"
            else:
                what = f"← left  gap vs {self.s2[pj - 1]}"
            msg = f"traceback  F[{pi}][{pj}] → F[{i}][{j}]   {what}"
        if self.animate and not quiet:
            self.status.config(text=msg)
            self._draw_alignment(k)

    def _draw_alignment(self, k):
        al, tw = self.al, self.tw
        al.delete("all")
        L = len(self.a1)
        x0, y1, y2 = 40, 4, 4 + tw * 2 + 6
        f = (MONO, max(7, tw // 2 + 1), "bold")
        al.create_text(18, y1 + tw / 2, text="S1", fill=MUTED, font=(UI, 8, "bold"))
        al.create_text(18, y2 + tw / 2, text="S2", fill=MUTED, font=(UI, 8, "bold"))
        for c in range(L):
            x = x0 + c * tw
            if c < L - k:
                for y in (y1, y2):
                    round_rect(al, x + 1, y, x + tw - 1, y + tw, max(2, tw // 5), fill=EMPTY, outline="")
                continue
            p, q = self.a1[c], self.a2[c]
            for y, ch in ((y1, p), (y2, q)):
                round_rect(al, x + 1, y, x + tw - 1, y + tw, max(2, tw // 5), fill=CARD2, outline="")
                al.create_text(x + tw / 2, y + tw / 2, text=ch, fill=BASE[ch], font=f)
            cx, top, bot = x + tw / 2, y1 + tw + 1, y2 - 1
            if p == q:
                al.create_line(cx, top, cx, bot, fill=GOOD, width=2)
            elif "-" not in (p, q):
                al.create_oval(cx - 2, (top + bot) / 2 - 2, cx + 2, (top + bot) / 2 + 2,
                               fill=WARN, outline="")

    def _finish(self):
        self._draw_alignment(len(self.a1))
        self.cv.delete("hl")
        self._header_hl()
        L = len(self.a1)
        matches = self.mk.count("|")
        sim = matches / L
        vals = {"Score": self.F[-1][-1], "Length": L, "Matches": matches,
                "Mismatches": self.mk.count("."), "Gaps": self.a1.count("-") + self.a2.count("-"),
                "Similarity": f"{100 * sim:.1f}%"}
        for k, v in vals.items():
            self.stat[k].config(text=str(v))
        self.stat["Score"].config(fg=CYAN)
        self.simbar.delete("all")
        self.simbar.create_rectangle(0, 0, 90 * sim, 4, fill=GOOD, outline="")
        self.status.config(text=f"done — optimal global score {self.F[-1][-1]}")
        self.running = False
        if self.animate:
            self.play_btn.config(text="Play")
            self.play_btn.set_enabled(False)

    # ---------- zoom
    def _render(self):
        """Redraw the matrix at the current zoom, replaying everything shown so far."""
        self.cv.delete("all")
        self._draw_grid()
        last = self.pos - 1
        for idx in range(self.pos):
            ev = self.events[idx]
            if ev[0] == "fill":
                self._fill(ev[1], ev[2], quiet=idx != last)
            elif ev[0] == "trace":
                self._trace(ev[1], quiet=idx != last)
        if self.pos and self.events[last][0] == "done":
            self.cv.delete("hl")
        self.zoom_lbl.config(text=f"{round(100 * self.cs / self.base_cs)}%")

    def zoom(self, factor, x=None, y=None):
        cv = self.cv
        if x is None:
            x, y = cv.winfo_width() / 2, cv.winfo_height() / 2
        new = max(self.min_cs, min(self.max_cs, round(self.cs * factor)))
        if new == self.cs:
            return
        cx, cy = cv.canvasx(x), cv.canvasy(y)
        f = new / self.cs
        self.cs = new
        self._render()
        W = (len(self.s2) + 2) * new
        H = (len(self.s1) + 2) * new
        cv.xview_moveto(max(0.0, (cx * f - x) / W))
        cv.yview_moveto(max(0.0, (cy * f - y) / H))

    def _wheel_zoom(self, e, zoom_in):
        self.zoom(1.15 if zoom_in else 1 / 1.15, e.x, e.y)
        return "break"

    def zoom_fit(self):
        self.cs = self.base_cs
        self._render()
        self.cv.xview_moveto(0)
        self.cv.yview_moveto(0)

    # ---------- controls
    def _progress(self):
        if self.animate:
            self.prog.delete("all")
            w = self.prog.winfo_width() * self.pos / len(self.events)
            self.prog.create_rectangle(0, 0, w, 3, fill=ACCENT, outline="")

    def _delay(self):
        d = int(600 / max(1.0, self.speed.get()))
        return d if self.events[min(self.pos, len(self.events) - 1)][0] == "fill" else d * 3

    def _tick(self):
        self.job = None
        if not self.running or self.pos >= len(self.events):
            return
        self._do(self.events[self.pos])
        self.pos += 1
        self._progress()
        if self.running:
            self.job = self.after(self._delay(), self._tick)

    def _cancel(self):
        if self.job:
            self.after_cancel(self.job)
            self.job = None

    def toggle(self):
        if self.pos >= len(self.events):
            return
        self.running = not self.running
        self.play_btn.config(text="Pause" if self.running else "Play")
        if self.running:
            self._tick()
        else:
            self._cancel()

    def step(self):
        if self.running:
            self.toggle()
        if self.pos < len(self.events):
            self._do(self.events[self.pos])
            self.pos += 1
            self._progress()

    def skip(self):
        self._cancel()
        self.running = False
        while self.pos < len(self.events):
            self._do(self.events[self.pos])
            self.pos += 1
        self._progress()

    def restart(self):
        self._cancel()
        self.pos, self.running = 0, True
        self._render()
        self._draw_alignment(0)
        self.status.config(text="")
        for v in self.stat.values():
            v.config(text="-", fg=TEXT)
        self.simbar.delete("all")
        self.play_btn.config(text="Pause")
        self.play_btn.set_enabled(True)
        self._tick()

    def destroy(self):
        self._cancel()
        super().destroy()


class LiveViewer(tk.Toplevel):
    """Popup window with the step-by-step animation."""

    def __init__(self, master, s1, s2, match, mismatch, gap):
        super().__init__(master, bg=BG, padx=16, pady=16)
        self.title("Needleman–Wunsch · step by step")
        dark_titlebar(self)
        AlignView(self, s1, s2, match, mismatch, gap, animate=True).pack(fill="both", expand=True)
        self.update_idletasks()
        w, h = self.winfo_reqwidth(), self.winfo_reqheight()
        self.geometry(f"{w}x{h}")
        self.minsize(min(w, 520), min(h, 480))


# ---------------------------------------------------------------- main window
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        init_fonts(self)
        self.title("DNA Sequence Aligner")
        self.configure(bg=BG, padx=18, pady=16)
        dark_titlebar(self)

        head = tk.Frame(self, bg=BG)
        head.pack(fill="x", pady=(0, 12))
        tk.Label(head, text="DNA Sequence Aligner", bg=BG, fg=TEXT, font=(UI, 18, "bold")).pack(anchor="w")
        tk.Label(head, text="Needleman–Wunsch global alignment  ·  updates live as you type",
                 bg=BG, fg=MUTED, font=(UI, 10)).pack(anchor="w")
        legend = tk.Frame(head, bg=BG)
        legend.pack(anchor="w", pady=(6, 0))
        for b in "ACGT":
            tk.Label(legend, text="●", bg=BG, fg=BASE[b], font=(UI, 10)).pack(side="left")
            tk.Label(legend, text=b + "  ", bg=BG, fg=MUTED, font=(MONO, 9, "bold")).pack(side="left")

        ic = card(self, padx=14, pady=12)
        ic.pack(fill="x")
        self.vars, self.lens = {}, {}
        for row, (name, val) in enumerate([("Sequence 1", S1_DEFAULT), ("Sequence 2", S2_DEFAULT)]):
            tk.Label(ic, text=name.upper(), bg=CARD, fg=MUTED, font=(UI, 8, "bold")).grid(
                row=row, column=0, sticky="w", padx=(0, 12))
            v = tk.StringVar(value=val)
            tk.Entry(ic, textvariable=v, width=48, font=(MONO, 12), bg=CARD2, fg=TEXT,
                     insertbackground=TEXT, relief="flat", highlightthickness=1,
                     highlightbackground=LINE, highlightcolor=ACCENT).grid(
                row=row, column=1, columnspan=6, sticky="we", ipady=6, pady=3)
            self.lens[name] = tk.Label(ic, text="", bg=CARD, fg=MUTED, font=(MONO, 9))
            self.lens[name].grid(row=row, column=7, sticky="w", padx=(10, 0))
            self.vars[name] = v

        sr = tk.Frame(ic, bg=CARD)
        sr.grid(row=2, column=0, columnspan=8, sticky="we", pady=(10, 0))
        for name, val, col in (("Match", 1, GOOD), ("Mismatch", -1, WARN), ("Gap", -1, BAD)):
            tk.Label(sr, text="●", bg=CARD, fg=col, font=(UI, 9)).pack(side="left")
            tk.Label(sr, text=name.upper(), bg=CARD, fg=MUTED, font=(UI, 8, "bold")).pack(
                side="left", padx=(3, 6))
            v = tk.StringVar(value=str(val))
            tk.Spinbox(sr, from_=-10, to=10, width=4, textvariable=v, font=(MONO, 11), bg=CARD2,
                       fg=TEXT, buttonbackground=CARD2, insertbackground=TEXT, relief="flat",
                       highlightthickness=1, highlightbackground=LINE, highlightcolor=ACCENT,
                       justify="center").pack(side="left", padx=(0, 18), ipady=3)
            self.vars[name] = v
        Pill(sr, "▶  Animate step by step", self.align, primary=True).pack(side="right")

        self.err = tk.Label(self, text="", bg=BG, fg=BAD, font=(UI, 9))
        self.err.pack(anchor="w", pady=(6, 0))
        self.holder = tk.Frame(self, bg=BG)
        self.holder.pack(fill="both", expand=True)

        self.view, self.pending = None, None
        for v in self.vars.values():
            v.trace_add("write", lambda *_: self._schedule())
        self.bind("<Return>", lambda e: self.align())
        self.refresh()

    def read_inputs(self):
        s1 = self.vars["Sequence 1"].get().strip().upper()
        s2 = self.vars["Sequence 2"].get().strip().upper()
        if not s1 or not s2 or set(s1 + s2) - set("ACGT"):
            raise ValueError("Sequences must be non-empty and contain only A, C, G, T.")
        try:
            mt, mm, gp = (int(self.vars[k].get()) for k in ("Match", "Mismatch", "Gap"))
        except ValueError:
            raise ValueError("Match, mismatch and gap must be integers.")
        return s1, s2, mt, mm, gp

    def _schedule(self):
        if self.pending:
            self.after_cancel(self.pending)
        self.pending = self.after(150, self.refresh)

    def refresh(self):
        self.pending = None
        for name in ("Sequence 1", "Sequence 2"):
            self.lens[name].config(text=f"{len(self.vars[name].get().strip())} bp")
        try:
            args = self.read_inputs()
        except ValueError as e:
            self.err.config(text=str(e))
            return
        self.err.config(text="")
        if self.view:
            self.view.destroy()
        self.view = AlignView(self.holder, *args, animate=False, max_px=470)
        self.view.pack(fill="both", expand=True)

    def align(self):
        try:
            args = self.read_inputs()
        except ValueError as e:
            self.err.config(text=str(e))
            return
        LiveViewer(self, *args)


if __name__ == "__main__":
    App().mainloop()
