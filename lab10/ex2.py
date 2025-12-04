import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import matplotlib as mpl

mpl.rcParams.update({
    "font.size": 8,
})


def clean_seq(seq: str) -> str:
    return ''.join(b for b in seq.upper() if b in "ACGT")


def cg_total(seq: str) -> float:
    s = clean_seq(seq)
    if not s:
        return 0.0
    cg = sum(1 for b in s if b in "CG")
    return 100.0 * cg / len(s)


def cpg_total(seq: str) -> int:
    s = clean_seq(seq)
    return sum(1 for i in range(len(s) - 1) if s[i:i+2] == "CG")


def cg_window_relative(window: str, total_len: int) -> float:
    w = clean_seq(window)
    if not w or total_len == 0:
        return 0.0
    cg = sum(1 for b in w if b in "CG")
    return 100.0 * cg / total_len


def cpg_window(window: str) -> int:
    w = clean_seq(window)
    return sum(1 for i in range(len(w) - 1) if w[i:i+2] == "CG")


def kappa_ic_window(window: str) -> float:
    A = clean_seq(window)
    L = len(A)
    if L < 2:
        return 0.0

    N = L - 1
    T = 0.0
    for u in range(1, N + 1):
        lengthB = L - u
        C = 0
        for i in range(lengthB):
            if A[i] == A[i + u]:
                C += 1
        T += (C / lengthB) * 100.0
    return T / N


def sliding_metrics(seq: str, window_len: int, step: int = 1):
    s = clean_seq(seq)
    if len(s) < window_len:
        raise ValueError("Sequence shorter than window length.")

    total_len = len(s)

    positions, cg_list, ic_list, cpg_list = [], [], [], []

    for start in range(0, len(s) - window_len + 1, step):
        w = s[start:start + window_len]
        center_pos = start + window_len // 2 + 1

        positions.append(center_pos)
        cg_list.append(cg_window_relative(w, total_len))
        ic_list.append(kappa_ic_window(w))
        cpg_list.append(cpg_window(w))

    return positions, cg_list, ic_list, cpg_list


def center_of_weight(xs, ys):
    if not xs:
        return (0.0, 0.0)
    x_c = sum(xs) / len(xs)
    y_c = sum(ys) / len(ys)
    return (x_c, y_c)


class TwoFastasGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PromKappa – Influenza vs COVID comparison")
        self.geometry("1400x800")

        self.fasta_paths = ["", ""]
        self.fasta_seqs = [[], []]

        self.create_widgets()
        self.create_plots()

    def create_widgets(self):
        left = ttk.Frame(self, padding=5)
        left.pack(side=tk.LEFT, fill=tk.Y)

        params = ttk.LabelFrame(left, text="Parameters", padding=5)
        params.pack(side=tk.TOP, fill=tk.X, pady=5)

        ttk.Label(params, text="Window length:").grid(row=0, column=0, sticky="w")
        self.window_len_var = tk.StringVar(value="30")
        ttk.Entry(params, textvariable=self.window_len_var, width=6).grid(
            row=0, column=1, sticky="w"
        )

        ttk.Label(params, text="Window step:").grid(row=1, column=0, sticky="w")
        self.window_step_var = tk.StringVar(value="1")
        ttk.Entry(params, textvariable=self.window_step_var, width=6).grid(
            row=1, column=1, sticky="w"
        )

        fasta_frame = ttk.LabelFrame(left, text="FASTA files", padding=5)
        fasta_frame.pack(side=tk.TOP, fill=tk.X, pady=5)

        row = 0
        ttk.Label(fasta_frame, text="Influenza FASTA:").grid(row=row, column=0, sticky="w")
        self.flu_label = ttk.Label(fasta_frame, text="(none)")
        self.flu_label.grid(row=row, column=1, sticky="w")
        ttk.Button(
            fasta_frame, text="Load...", command=lambda: self.load_fasta(0)
        ).grid(row=row, column=2, padx=2)

        row += 1
        ttk.Label(fasta_frame, text="COVID FASTA:").grid(row=row, column=0, sticky="w")
        self.covid_label = ttk.Label(fasta_frame, text="(none)")
        self.covid_label.grid(row=row, column=1, sticky="w")
        ttk.Button(
            fasta_frame, text="Load...", command=lambda: self.load_fasta(1)
        ).grid(row=row, column=2, padx=2)

        ttk.Button(left, text="Analyze FASTA files", command=self.analyze_fastas).pack(
            side=tk.TOP, pady=10, fill=tk.X
        )

        res = ttk.LabelFrame(left, text="Info", padding=5)
        res.pack(side=tk.TOP, fill=tk.X, pady=5)
        self.info_label = ttk.Label(
            res, text="Load Influenza and COVID FASTA files and press Analyze."
        )
        self.info_label.pack(anchor="w")

    def create_plots(self):
        mid = ttk.Frame(self, padding=5)
        mid.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.fig = Figure(figsize=(10, 7), dpi=100)

        gs = GridSpec(5, 2, figure=self.fig)

        self.fig.subplots_adjust(
            left=0.06,
            right=0.99,
            top=0.97,
            bottom=0.06,
            wspace=0.35,
            hspace=0.75,
        )

        self.axes = [[None for _ in range(2)] for _ in range(5)]

        for col in range(2):
            self.axes[0][col] = self.fig.add_subplot(gs[0, col])
            self.axes[1][col] = self.fig.add_subplot(gs[1, col])
            self.axes[2][col] = self.fig.add_subplot(gs[2, col])
            self.axes[3][col] = self.fig.add_subplot(gs[3, col])
            self.axes[4][col] = self.fig.add_subplot(gs[4, col])

        for row in range(5):
            for col in range(2):
                self.axes[row][col].tick_params(labelsize=7)

        self.canvas = FigureCanvasTkAgg(self.fig, master=mid)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.draw()

    def load_fasta(self, idx: int):
        label_name = "Influenza" if idx == 0 else "COVID"
        path = filedialog.askopenfilename(
            title=f"Open {label_name} FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.txt"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "r") as f:
                lines = f.readlines()
        except Exception as e:
            messagebox.showerror("Error", f"Could not read file:\n{e}")
            return

        seqs = []
        current = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                    current = []
            else:
                current.append(line)
        if current:
            seqs.append("".join(current))

        if not seqs:
            messagebox.showwarning("Warning", "No sequences found in FASTA.")
            return

        self.fasta_paths[idx] = path
        self.fasta_seqs[idx] = seqs

        label = self.flu_label if idx == 0 else self.covid_label
        label.config(text=f"{len(seqs)} sequences loaded")

    def analyze_fastas(self):
        try:
            win_len = int(self.window_len_var.get())
            step = int(self.window_step_var.get())
            if win_len <= 1 or step <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Error", "Window length must be >1 and step must be >0 (integers)."
            )
            return

        if not self.fasta_seqs[0] or not self.fasta_seqs[1]:
            messagebox.showwarning("Warning", "Load both FASTA files first.")
            return

        for col in range(2):
            self.update_plots_for_fasta(col, win_len, step)

        self.canvas.draw()

    def update_plots_for_fasta(self, col: int, win_len: int, step: int):
        seqs = self.fasta_seqs[col]

        sliding_results = []
        pattern_cg_x, pattern_cg_y = [], []
        pattern_cpg_x, pattern_cpg_y = [], []
        centers_cg = []
        centers_cpg = []

        for seq in seqs:
            try:
                positions, cg_list, ic_list, cpg_list = sliding_metrics(
                    seq, win_len, step
                )
            except ValueError:
                continue

            sliding_results.append((positions, cg_list, ic_list, cpg_list))

            pattern_cg_x.extend(cg_list)
            pattern_cg_y.extend(ic_list)
            pattern_cpg_x.extend(cpg_list)
            pattern_cpg_y.extend(ic_list)

            centers_cg.append(center_of_weight(cg_list, ic_list))
            centers_cpg.append(center_of_weight(cpg_list, ic_list))

        ax_top = self.axes[0][col]
        ax_patt_cg = self.axes[1][col]
        ax_patt_cpg = self.axes[2][col]
        ax_cent_cg = self.axes[3][col]
        ax_cent_cpg = self.axes[4][col]

        ax_top.clear()
        first = True
        for positions, cg_list, ic_list, cpg_list in sliding_results:
            ax_top.plot(
                positions,
                ic_list,
                color="tab:blue",
                alpha=0.5,
                label="Kappa IC" if first else "_nolegend_",
            )
            ax_top.plot(
                positions,
                cg_list,
                color="tab:orange",
                alpha=0.5,
                label="(C+G)% (rel.)" if first else "_nolegend_",
            )
            ax_top.plot(
                positions,
                cpg_list,
                color="tab:green",
                alpha=0.5,
                label="CG count" if first else "_nolegend_",
            )
            first = False

        ax_top.set_xlabel("Position (bp)", fontsize=8)
        ax_top.set_ylabel("Value", fontsize=8)
        label_name = "Influenza" if col == 0 else "COVID"
        title = f"Sliding-window metrics – {label_name}"
        ax_top.set_title(title, fontsize=9)
        ax_top.grid(True)
        ax_top.legend(fontsize=7, loc="upper right")
        ax_top.tick_params(labelsize=7)

        ax_patt_cg.clear()
        if pattern_cg_x:
            cx, cy = center_of_weight(pattern_cg_x, pattern_cg_y)
            ax_patt_cg.scatter(pattern_cg_x, pattern_cg_y, s=6, color="tab:blue")
            ax_patt_cg.axvline(cx, linestyle="--", linewidth=0.7)
            ax_patt_cg.axhline(cy, linestyle="--", linewidth=0.7)
        ax_patt_cg.set_xlabel("C+G% (relative)", fontsize=8)
        ax_patt_cg.set_ylabel("Kappa IC", fontsize=8)
        ax_patt_cg.set_title(f"Pattern IC vs C+G% – {label_name}", fontsize=9)
        ax_patt_cg.grid(True)
        ax_patt_cg.tick_params(labelsize=7)

        ax_patt_cpg.clear()
        if pattern_cpg_x:
            cx2, cy2 = center_of_weight(pattern_cpg_x, pattern_cpg_y)
            ax_patt_cpg.scatter(pattern_cpg_x, pattern_cpg_y, s=6, color="tab:red")
            ax_patt_cpg.axvline(cx2, linestyle="--", linewidth=0.7)
            ax_patt_cpg.axhline(cy2, linestyle="--", linewidth=0.7)
        ax_patt_cpg.set_xlabel("CG (count)", fontsize=8)
        ax_patt_cpg.set_ylabel("Kappa IC", fontsize=8)
        ax_patt_cpg.set_title(f"Pattern IC vs CG count – {label_name}", fontsize=9)
        ax_patt_cpg.grid(True)
        ax_patt_cpg.tick_params(labelsize=7)

        ax_cent_cg.clear()
        if centers_cg:
            xs = [c[0] for c in centers_cg]
            ys = [c[1] for c in centers_cg]
            ax_cent_cg.scatter(xs, ys, s=25, marker="o")
        ax_cent_cg.set_xlabel("C+G% (center per genome)", fontsize=8)
        ax_cent_cg.set_ylabel("Kappa IC", fontsize=8)
        ax_cent_cg.set_title(f"Centers (C+G%, IC) – {label_name}", fontsize=9)
        ax_cent_cg.grid(True)
        ax_cent_cg.tick_params(labelsize=7)

        ax_cent_cpg.clear()
        if centers_cpg:
            xs2 = [c[0] for c in centers_cpg]
            ys2 = [c[1] for c in centers_cpg]
            ax_cent_cpg.scatter(xs2, ys2, s=25, marker="o")
        ax_cent_cpg.set_xlabel("CG (center per genome)", fontsize=8)
        ax_cent_cpg.set_ylabel("Kappa IC", fontsize=8)
        ax_cent_cpg.set_title(f"Centers (CG, IC) – {label_name}", fontsize=9)
        ax_cent_cpg.grid(True)
        ax_cent_cpg.tick_params(labelsize=7)

        info = []
        if self.fasta_seqs[0]:
            info.append(f"Influenza: {len(self.fasta_seqs[0])} sequences")
        if self.fasta_seqs[1]:
            info.append(f"COVID: {len(self.fasta_seqs[1])} sequences")
        self.info_label.config(text="; ".join(info))


if __name__ == "__main__":
    app = TwoFastasGUI()
    app.mainloop()