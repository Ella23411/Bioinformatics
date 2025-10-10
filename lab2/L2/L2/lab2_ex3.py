import tkinter as tk
from tkinter import filedialog, messagebox
from typing import List, Dict
import os

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from collections import Counter, deque
import string
import math

def read_fasta(path: str) -> str:
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts).upper()

def detect_alphabet(seq: str, restrict_acgt: bool) -> List[str]:
    seq_set = set(seq)
    if restrict_acgt:
        core = {'A', 'C', 'G', 'T'}
        return sorted(core.intersection(seq_set))
    else:
        return sorted(ch for ch in seq_set if ch in set(string.ascii_uppercase))

def sliding_window_relative_freq(seq: str, k: int, alphabet: List[str]) -> Dict[str, List[float]]:
    n = len(seq)
    res = {sym: [] for sym in alphabet}
    if n < k or k <= 0:
        return res

    window = seq[0:k]
    cnt = Counter(ch for ch in window if ch in alphabet)

    def append_freqs():
        for sym in alphabet:
            res[sym].append(cnt.get(sym, 0) / k * 100.0) 

    append_freqs()

    for i in range(1, n - k + 1):
        out_ch = seq[i - 1]
        in_ch = seq[i + k - 1]
        if out_ch in alphabet:
            cnt[out_ch] -= 1
            if cnt[out_ch] == 0:
                del cnt[out_ch]
        if in_ch in alphabet:
            cnt[in_ch] = cnt.get(in_ch, 0) + 1
        append_freqs()

    return res

def smooth_series(values: List[float], span: int) -> List[float]:
    n = len(values)
    if n == 0 or span <= 1 or span > n:
        return values[:]
    if span % 2 == 0:
        span += 1
    half = span // 2
    padded = [values[0]] * half + values + [values[-1]] * half
    win_sum = sum(padded[:span])
    out = [win_sum / span]
    dq = deque(padded[:span], maxlen=span)
    for x in padded[span:]:
        left = dq[0]
        dq.append(x)
        win_sum += x - left
        out.append(win_sum / span)
    start = 0
    end = start + n
    return out[start:end]

# ---------- GUI ----------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Sliding Window Relative Frequencies (%)")
        self.geometry("1000x700")
        self.seq = ""
        self.file_path = ""
        self.restrict_var = tk.BooleanVar(value=True)
        self.smooth_var = tk.BooleanVar(value=False)  
        self.smooth_span = tk.StringVar(value="9")

        top = tk.Frame(self)
        top.pack(fill="x", padx=10, pady=8)

        self.path_var = tk.StringVar(value="No file selected")
        tk.Button(top, text="Open FASTA…", command=self.open_fasta).pack(side="left")
        tk.Label(top, textvariable=self.path_var).pack(side="left", padx=10)

        tk.Label(top, text="Window (k):").pack(side="left", padx=(20, 4))
        self.k_entry = tk.Entry(top, width=6)
        self.k_entry.insert(0, "30")
        self.k_entry.pack(side="left")

        tk.Checkbutton(top, text="Restrict to A/C/G/T", variable=self.restrict_var).pack(side="left", padx=10)

  
        tk.Checkbutton(top, text="Smooth", variable=self.smooth_var).pack(side="left", padx=(20, 4))
        tk.Label(top, text="Span:").pack(side="left")
        self.span_entry = tk.Entry(top, width=6, textvariable=self.smooth_span)
        self.span_entry.pack(side="left")

        tk.Button(top, text="Analyze & Plot", command=self.analyze_and_plot).pack(side="right")

        stats = tk.Frame(self)
        stats.pack(fill="x", padx=10, pady=(0, 6))
        self.stats_var = tk.StringVar(value="Sequence length: –  |  Windows: –  |  Alphabet: –")
        tk.Label(stats, textvariable=self.stats_var).pack(anchor="w")

        self.fig = Figure(figsize=(9.6, 5.8), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("Sliding Window Relative Frequencies (%)")
        self.ax.set_xlabel("Window start index (0-based)")
        self.ax.set_ylabel("Relative frequency (%)")
        self.ax.set_ylim(0, 100)
        self.ax.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

    def open_fasta(self):
        path = filedialog.askopenfilename(
            title="Choose a FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            seq = read_fasta(path)
            if not seq:
                messagebox.showwarning("Empty", "The selected FASTA contains no sequence data.")
                return
            self.seq = seq
            self.file_path = path
            self.path_var.set(os.path.basename(path))
            self.stats_var.set(f"Sequence length: {len(seq)}  |  Windows: –  |  Alphabet: –")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read FASTA:\n{e}")

    def analyze_and_plot(self):
        if not self.seq:
            messagebox.showwarning("No sequence", "Please open a FASTA file first.")
            return
        try:
            k = int(self.k_entry.get())
            if k <= 0:
                raise ValueError
        except Exception:
            messagebox.showwarning("Invalid window", "Window size (k) must be a positive integer.")
            return

  
        try:
            span = int(self.span_entry.get())
            if span < 1:
                span = 1
        except Exception:
            span = 1

        alphabet = detect_alphabet(self.seq, self.restrict_var.get())
        if not alphabet:
            messagebox.showwarning("No alphabet", "No valid symbols found.")
            return

        freqs = sliding_window_relative_freq(self.seq, k, alphabet)
        windows = max(len(self.seq) - k + 1, 0)
        x = list(range(windows))

        self.ax.clear()
        self.ax.set_title("Sliding Window Relative Frequencies (%)")
        self.ax.set_xlabel("Window start index (0-based)")
        self.ax.set_ylabel("Relative frequency (%)")
        self.ax.set_ylim(0, 100)
        self.ax.grid(True)

        do_smooth = self.smooth_var.get() and span > 1
        for sym in alphabet:
            y = freqs[sym]
            if not y:
                continue
            y_plot = smooth_series(y, span) if do_smooth else y
            self.ax.plot(x, y_plot, label=sym)

        self.ax.legend(title="Symbol", ncol=min(len(alphabet), 6))
        self.canvas.draw_idle()

        self.stats_var.set(
            f"Sequence length: {len(self.seq)}  |  Windows: {windows}  |  Alphabet: {','.join(alphabet)}"
        )

if __name__ == "__main__":
    app = App()
    app.mainloop()
