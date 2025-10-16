import math
import tkinter as tk
from tkinter import filedialog, messagebox
from typing import List, Tuple

import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt


def read_fasta(path: str) -> str:
    sequence_parts: List[str] = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            sequence_parts.append(line.upper())
    if not sequence_parts:
        raise ValueError(f"No sequence data found in {path}")
    return "".join(sequence_parts)


def compute_melting_temperatures(
    sequence: str, window_size: int, na_conc: float
) -> Tuple[List[int], List[float], List[float]]:
    if window_size <= 0:
        raise ValueError("Window size must be positive")
    if na_conc <= 0:
        raise ValueError("Sodium concentration must be positive")
    seq_len = len(sequence)
    if seq_len < window_size:
        raise ValueError(
            f"Sequence length ({seq_len}) is shorter than the window size ({window_size})"
        )

    positions: List[int] = []
    tm1_list: List[float] = []
    tm2_list: List[float] = []
    half_window = window_size // 2

    for start in range(seq_len - window_size + 1):
        window_seq = sequence[start : start + window_size]
        window_seq_up = window_seq.upper()
        a_count = window_seq_up.count("A")
        t_count = window_seq_up.count("T")
        g_count = window_seq_up.count("G")
        c_count = window_seq_up.count("C")

        tm1 = 4 * (g_count + c_count) + 2 * (a_count + t_count)
        length = window_size
        gc_fraction = (g_count + c_count) / float(length)
        gc_percent = gc_fraction * 100.0
        tm2 = (
            81.5
            + 16.6 * math.log10(na_conc)
            + 0.41 * gc_percent
            - 600.0 / length
        ) * -1

        centre_position = start + half_window + 1
        positions.append(centre_position)
        tm1_list.append(tm1)
        tm2_list.append(tm2)

    return positions, tm1_list, tm2_list


def segments_above_threshold(positions: List[int], values: List[float], thr: float):
    segs = []
    in_seg = False
    seg_start = None
    last_pos = None
    for pos, val in zip(positions, values):
        if val >= thr:
            if not in_seg:
                in_seg = True
                seg_start = pos
            last_pos = pos
        else:
            if in_seg:
                width = last_pos - seg_start + 1
                segs.append((seg_start, width))
                in_seg = False
    if in_seg:
        width = last_pos - seg_start + 1
        segs.append((seg_start, width))
    return segs


class TmViewer(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("DNA Melting Temperature Viewer")
        self.geometry("1000x700")
        self.sequence = None
        self.positions = []
        self.tm1_vals = []
        self.tm2_vals = []
        self._create_widgets()

    def _create_widgets(self) -> None:
        ctrl_frame = tk.Frame(self)
        ctrl_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=8)

        open_btn = tk.Button(
            ctrl_frame, text="Open FASTA File", command=self._open_fasta
        )
        open_btn.pack(side=tk.LEFT, padx=(0, 10))

        tk.Label(ctrl_frame, text="Window size:").pack(side=tk.LEFT)
        self.window_size_var = tk.IntVar(value=9)
        tk.Entry(ctrl_frame, textvariable=self.window_size_var, width=5).pack(
            side=tk.LEFT, padx=(4, 12)
        )

        tk.Label(ctrl_frame, text="[Na⁺] (M):").pack(side=tk.LEFT)
        self.na_conc_var = tk.DoubleVar(value=0.001)
        tk.Entry(ctrl_frame, textvariable=self.na_conc_var, width=8).pack(
            side=tk.LEFT, padx=(4, 12)
        )

        tk.Label(ctrl_frame, text="Threshold (°C):").pack(side=tk.LEFT)
        self.threshold_var = tk.DoubleVar(value=0.0)
        tk.Entry(ctrl_frame, textvariable=self.threshold_var, width=10).pack(
            side=tk.LEFT, padx=(4, 8)
        )

        tk.Button(ctrl_frame, text="Update Plot", command=self._update_plots).pack(
            side=tk.LEFT, padx=(6, 0)
        )

        self.seq_label = tk.Label(ctrl_frame, text="No file loaded")
        self.seq_label.pack(side=tk.LEFT, padx=(16, 0))

        stats_frame = tk.Frame(self)
        stats_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=(0, 8))
        self.stats_label = tk.Label(
            stats_frame,
            text="Min/Max — F1: –, – | F2: –, –",
            anchor="w",
            justify="left",
        )
        self.stats_label.pack(side=tk.LEFT)

        self.figure, (self.ax_main, self.ax_mask) = plt.subplots(
            2, 1, figsize=(9.5, 5.6), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
        )
        self.figure.subplots_adjust(hspace=0.25)

        self.ax_main.set_xlabel("Position (centre of window)")
        self.ax_main.set_ylabel("Melting Temperature (°C)")
        self.ax_main.set_title("Melting Temperature along Sequence")
        self.ax_main.grid(True)

        self.ax_mask.set_xlabel("Position (centre of window)")
        self.ax_mask.set_yticks([0.5, 1.5])
        self.ax_mask.set_yticklabels(
            ["Above thr: Formula 1", "Above thr: Formula 2"]
        )
        self.ax_mask.set_ylim(0, 2)
        self.ax_mask.grid(True, axis="x")

        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _open_fasta(self) -> None:
        file_path = filedialog.askopenfilename(
            filetypes=[
                ("FASTA files", "*.fasta *.fa *.fna *.txt"),
                ("All files", "*.*"),
            ]
        )
        if not file_path:
            return
        try:
            self.sequence = read_fasta(file_path)
            self.seq_label.config(text=f"Sequence length: {len(self.sequence)} bases")
            self._recompute_and_plot()
        except Exception as exc:
            messagebox.showerror("Error", str(exc))

    def _recompute_and_plot(self):
        if not self.sequence:
            return

        window_size = self.window_size_var.get()
        na_conc = self.na_conc_var.get()

        self.positions, self.tm1_vals, self.tm2_vals = compute_melting_temperatures(
            self.sequence, window_size, na_conc
        )
        self._update_plots()

    def _update_plots(self):
        if not self.positions:
            return

        thr = self.threshold_var.get()

        self.ax_main.clear()
        self.ax_main.plot(
            self.positions,
            self.tm1_vals,
            label="Formula 1: 4×(G+C) + 2×(A+T)",
            color="tab:blue",
        )
        self.ax_main.plot(
            self.positions,
            self.tm2_vals,
            label="Formula 2: 81.5 + 16.6·log₁₀([Na⁺]) + 0.41·%GC – 600/length",
            color="tab:orange",
        )
        self.ax_main.axhline(thr, linestyle="--", linewidth=1, label=f"Threshold = {thr:+.4f} °C")
        self.ax_main.set_xlabel("Position (centre of window)")
        self.ax_main.set_ylabel("Melting Temperature (°C)")
        self.ax_main.set_title("Melting Temperature along Sequence")
        self.ax_main.legend(loc="best")
        self.ax_main.grid(True)

        f1_min, f1_max = min(self.tm1_vals), max(self.tm1_vals)
        f2_min, f2_max = min(self.tm2_vals), max(self.tm2_vals)
        self.stats_label.config(
            text=f"Min/Max — F1: {f1_min:+.4f} °C, {f1_max:+.4f} °C | F2: {f2_min:+.4f} °C, {f2_max:+.4f} °C"
        )

        self.ax_mask.clear()
        self.ax_mask.set_ylim(0, 2)
        self.ax_mask.set_yticks([0.5, 1.5])
        self.ax_mask.set_yticklabels(["Above thr: Formula 1", "Above thr: Formula 2"])
        self.ax_mask.set_xlabel("Position (centre of window)")
        self.ax_mask.grid(True, axis="x")

        segs_f1 = segments_above_threshold(self.positions, self.tm1_vals, thr)
        segs_f2 = segments_above_threshold(self.positions, self.tm2_vals, thr)

        for (start, width) in segs_f1:
            self.ax_mask.broken_barh([(start, width)], (0.1, 0.8))
        for (start, width) in segs_f2:
            self.ax_mask.broken_barh([(start, width)], (1.1, 0.8))

        x_min, x_max = self.positions[0], self.positions[-1]
        self.ax_mask.set_xlim(x_min, x_max)

        self.canvas.draw()


def main() -> None:
    matplotlib.use("TkAgg")
    app = TmViewer()
    app.mainloop()


if __name__ == "__main__":
    main()
