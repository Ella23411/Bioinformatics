from dataclasses import dataclass
from typing import List, Optional
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


WINDOW_SIZE = 3000
OVERLAP = 500
SCORE_THRESHOLD = 80


@dataclass
class LocalAlignmentResult:
    score: int
    s1_start: int
    s1_end: int
    s2_start: int
    s2_end: int
    aligned_s1: str
    aligned_s2: str


@dataclass
class Hit:
    score: int
    s1_start: int
    s1_end: int
    s2_start: int
    s2_end: int
    aligned_s1: str
    aligned_s2: str
    window_idx: int


def read_single_fasta(path: str) -> str:
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    return "".join(seq)


def smith_waterman(
    s1: str,
    s2: str,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> LocalAlignmentResult:
    n, m = len(s1), len(s2)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = H[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            delete = H[i - 1][j] + gap
            insert = H[i][j - 1] + gap
            val = max(0, diag, delete, insert)
            H[i][j] = val
            if val > max_score:
                max_score = val
                max_pos = (i, j)

    i, j = max_pos
    aligned1 = []
    aligned2 = []

    while i > 0 and j > 0 and H[i][j] > 0:
        score_current = H[i][j]
        score_diag = H[i - 1][j - 1]
        score_up = H[i - 1][j]
        score_left = H[i][j - 1]

        if score_current == score_diag + (match if s1[i - 1] == s2[j - 1] else mismatch):
            aligned1.append(s1[i - 1])
            aligned2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif score_current == score_up + gap:
            aligned1.append(s1[i - 1])
            aligned2.append("-")
            i -= 1
        elif score_current == score_left + gap:
            aligned1.append("-")
            aligned2.append(s2[j - 1])
            j -= 1
        else:
            break

    aligned1.reverse()
    aligned2.reverse()

    s1_start = i
    s1_end = max_pos[0]
    s2_start = j
    s2_end = max_pos[1]

    return LocalAlignmentResult(
        score=max_score,
        s1_start=s1_start,
        s1_end=s1_end,
        s2_start=s2_start,
        s2_end=s2_end,
        aligned_s1="".join(aligned1),
        aligned_s2="".join(aligned2),
    )


def windowed_local_align(
    s1: str,
    s2: str,
    window: int,
    overlap: int,
    score_threshold: int,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> List[Hit]:
    hits: List[Hit] = []

    if window <= 0:
        raise ValueError("window must be > 0")
    if overlap < 0:
        raise ValueError("overlap must be >= 0")
    step = window - overlap
    if step <= 0:
        raise ValueError("overlap must be smaller than window")

    w_idx = 0
    for start in range(0, len(s1), step):
        end = min(start + window, len(s1))
        if end <= start:
            break
        subseq = s1[start:end]

        res = smith_waterman(subseq, s2, match=match, mismatch=mismatch, gap=gap)

        if res.score >= score_threshold and res.s1_end > res.s1_start:
            global_s1_start = start + res.s1_start
            global_s1_end = start + res.s1_end
            hits.append(
                Hit(
                    score=res.score,
                    s1_start=global_s1_start,
                    s1_end=global_s1_end,
                    s2_start=res.s2_start,
                    s2_end=res.s2_end,
                    aligned_s1=res.aligned_s1,
                    aligned_s2=res.aligned_s2,
                    window_idx=w_idx,
                )
            )
        w_idx += 1

    return hits


def chain_hits(
    hits: List[Hit],
    min_distance: int = -100,
    max_distance: int = 5000,
) -> List[Hit]:
    if not hits:
        return []

    hits_sorted = sorted(hits, key=lambda h: h.s1_start)
    chain: List[Hit] = [hits_sorted[0]]
    last = hits_sorted[0]

    for h in hits_sorted[1:]:
        if (
            h.s1_start >= last.s1_end + min_distance
            and h.s2_start >= last.s2_end + min_distance
            and h.s1_start - last.s1_end <= max_distance
            and h.s2_start - last.s2_end <= max_distance
        ):
            chain.append(h)
            last = h

    return chain


def compute_similarity_metrics(
    aligned_s1: str,
    aligned_s2: str,
    match_score: int,
    mismatch_score: int,
    gap_score: int,
):
    M = X = G = 0

    for a, b in zip(aligned_s1, aligned_s2):
        if a == "-" or b == "-":
            G += 1
        elif a == b:
            M += 1
        else:
            X += 1

    identity_no_gaps = M / (M + X) if (M + X) > 0 else 0.0
    L = M + X + G
    identity_with_gaps = M / L if L > 0 else 0.0

    score_raw = M * match_score + X * mismatch_score + G * gap_score
    score_max = L * match_score if match_score > 0 else 1
    normalized_score = score_raw / score_max if score_max != 0 else 0.0

    return {
        "M": M,
        "X": X,
        "G": G,
        "L": L,
        "identity_no_gaps": identity_no_gaps,
        "identity_with_gaps": identity_with_gaps,
        "normalized_score": normalized_score,
    }


class AlignmentApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("COVID-19 vs Influenza Local Alignment (Windowed)")
        self.geometry("1200x800")

        self.covid_seq = ""
        self.flu_seq = ""

        self.chained_hits: List[Hit] = []

        self.current_match = 2
        self.current_mismatch = -1
        self.current_gap = -2

        self.fig = Figure(figsize=(6, 2.5), dpi=100)
        self.ax_overview = self.fig.add_subplot(1, 1, 1)

        self._build_ui()
        self.update_overview_plot(None)

    def _build_ui(self):
        file_frame = ttk.LabelFrame(self, text="Genome FASTA files")
        file_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(file_frame, text="COVID-19 FASTA:").grid(
            row=0, column=0, sticky="w", padx=5, pady=3
        )
        self.covid_path_var = tk.StringVar()
        covid_entry = ttk.Entry(file_frame, textvariable=self.covid_path_var, width=70)
        covid_entry.grid(row=0, column=1, padx=5, pady=3, sticky="we")
        ttk.Button(file_frame, text="Browse...", command=self.load_covid).grid(
            row=0, column=2, padx=5, pady=3
        )

        ttk.Label(file_frame, text="Influenza FASTA:").grid(
            row=1, column=0, sticky="w", padx=5, pady=3
        )
        self.flu_path_var = tk.StringVar()
        flu_entry = ttk.Entry(file_frame, textvariable=self.flu_path_var, width=70)
        flu_entry.grid(row=1, column=1, padx=5, pady=3, sticky="we")
        ttk.Button(file_frame, text="Browse...", command=self.load_flu).grid(
            row=1, column=2, padx=5, pady=3
        )

        file_frame.columnconfigure(1, weight=1)

        self.length_label = ttk.Label(
            file_frame,
            text="COVID length: -    Influenza length: -",
        )
        self.length_label.grid(row=2, column=0, columnspan=3, sticky="w", padx=5, pady=3)

        param_frame = ttk.LabelFrame(self, text="Alignment scoring parameters")
        param_frame.pack(fill="x", padx=10, pady=5)

        self.match_var = tk.StringVar(value="2")
        self.mismatch_var = tk.StringVar(value="-1")
        self.gap_var = tk.StringVar(value="-2")

        row = 0
        ttk.Label(param_frame, text="Match:").grid(
            row=row, column=0, sticky="e", padx=5, pady=3
        )
        ttk.Entry(param_frame, textvariable=self.match_var, width=8).grid(
            row=row, column=1, sticky="w", padx=5, pady=3
        )

        ttk.Label(param_frame, text="Mismatch:").grid(
            row=row, column=2, sticky="e", padx=5, pady=3
        )
        ttk.Entry(param_frame, textvariable=self.mismatch_var, width=8).grid(
            row=row, column=3, sticky="w", padx=5, pady=3
        )

        ttk.Label(param_frame, text="Gap penalty:").grid(
            row=row, column=4, sticky="e", padx=5, pady=3
        )
        ttk.Entry(param_frame, textvariable=self.gap_var, width=8).grid(
            row=row, column=5, sticky="w", padx=5, pady=3
        )

        run_frame = ttk.Frame(self)
        run_frame.pack(fill="x", padx=10, pady=5)
        self.run_button = ttk.Button(
            run_frame,
            text="Compare genomes (windowed local alignment)",
            command=self.run_alignment,
        )
        self.run_button.pack(side="left")

        results_pane = ttk.PanedWindow(self, orient="vertical")
        results_pane.pack(fill="both", expand=True, padx=10, pady=5)

        align_frame = ttk.LabelFrame(
            results_pane,
            text="Selected hit alignment (Green=Match, Red=Mismatch, Gray=Gap)",
        )
        self.align_text = tk.Text(
            align_frame,
            height=14,
            wrap="none",
            font=("Courier", 9),
        )
        self.align_text.pack(side="left", fill="both", expand=True)

        self.align_text.tag_config("match", background="#ccffcc", foreground="green")
        self.align_text.tag_config("mismatch", background="#ffcccc", foreground="red")
        self.align_text.tag_config("gap", background="#dddddd", foreground="black")

        align_scroll_y = ttk.Scrollbar(
            align_frame, orient="vertical", command=self.align_text.yview
        )
        align_scroll_y.pack(side="right", fill="y")

        align_scroll_x = ttk.Scrollbar(
            align_frame, orient="horizontal", command=self.align_text.xview
        )
        align_scroll_x.pack(side="bottom", fill="x")

        self.align_text.configure(
            yscrollcommand=align_scroll_y.set,
            xscrollcommand=align_scroll_x.set,
        )

        results_pane.add(align_frame, weight=3)

        plot_frame = ttk.LabelFrame(results_pane, text="Selected hit on COVID genome")
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        results_pane.add(plot_frame, weight=2)

    def load_covid(self):
        path = filedialog.askopenfilename(
            title="Select COVID-19 FASTA file",
            filetypes=(("FASTA files", "*.fa *.fasta *.fna *.txt"), ("All files", "*.*")),
        )
        if not path:
            return
        try:
            seq = read_single_fasta(path)
            if not seq:
                raise ValueError("No sequence found in file.")
            self.covid_seq = seq
            self.covid_path_var.set(path)
            self._update_lengths()
            self.update_overview_plot(None)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load COVID sequence:\n{e}")

    def load_flu(self):
        path = filedialog.askopenfilename(
            title="Select Influenza FASTA file",
            filetypes=(("FASTA files", "*.fa *.fasta *.fna *.txt"), ("All files", "*.*")),
        )
        if not path:
            return
        try:
            seq = read_single_fasta(path)
            if not seq:
                raise ValueError("No sequence found in file.")
            self.flu_seq = seq
            self.flu_path_var.set(path)
            self._update_lengths()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load Influenza sequence:\n{e}")

    def _update_lengths(self):
        covid_len = len(self.covid_seq) if self.covid_seq else 0
        flu_len = len(self.flu_seq) if self.flu_seq else 0
        self.length_label.config(
            text=f"COVID length: {covid_len} bp    Influenza length: {flu_len} bp"
        )

    def run_alignment(self):
        if not self.covid_seq or not self.flu_seq:
            messagebox.showwarning("Missing data", "Please load both FASTA files first.")
            return

        try:
            match = int(self.match_var.get())
            mismatch = int(self.mismatch_var.get())
            gap = int(self.gap_var.get())
        except ValueError:
            messagebox.showerror(
                "Invalid parameters", "Match/mismatch/gap must be integers."
            )
            return

        self.current_match = match
        self.current_mismatch = mismatch
        self.current_gap = gap

        try:
            hits_all = windowed_local_align(
                self.covid_seq,
                self.flu_seq,
                window=WINDOW_SIZE,
                overlap=OVERLAP,
                score_threshold=SCORE_THRESHOLD,
                match=match,
                mismatch=mismatch,
                gap=gap,
            )
        except Exception as e:
            messagebox.showerror("Error during alignment", str(e))
            return

        chained = chain_hits(hits_all)
        self.chained_hits = chained

        self.align_text.delete("1.0", "end")
        if chained:
            best_hit = max(chained, key=lambda h: h.score)
            self.show_hit(best_hit)
        else:
            self.update_overview_plot(None)
            messagebox.showinfo("Result", "No hits above threshold were found.")

    def update_overview_plot(self, hit: Optional[Hit]):
        self.ax_overview.clear()

        covid_len = len(self.covid_seq) if self.covid_seq else 0
        if covid_len == 0:
            self.canvas.draw()
            return

        self.ax_overview.plot([0, covid_len], [0, 0], linewidth=2)
        self.ax_overview.set_xlim(0, covid_len)
        self.ax_overview.set_ylim(-0.5, 1.5)
        self.ax_overview.set_yticks([0, 1])
        self.ax_overview.set_yticklabels(["COVID-19", "Influenza"])

        if hit is not None:
            runs = []
            cur_pos = hit.s1_start
            in_run = False
            run_start = None

            for a, b in zip(hit.aligned_s1, hit.aligned_s2):
                if a != "-":
                    covid_coord = cur_pos
                    if a == b and a != "-":
                        if not in_run:
                            in_run = True
                            run_start = covid_coord
                    else:
                        if in_run:
                            runs.append((run_start, covid_coord))
                            in_run = False
                    cur_pos += 1
                else:
                    if in_run:
                        runs.append((run_start, cur_pos - 1))
                        in_run = False

            if in_run:
                runs.append((run_start, cur_pos - 1))

            for start, end in runs:
                if start is None or end is None:
                    continue
                if end < start:
                    continue
                self.ax_overview.plot([start, end], [1, 1], linewidth=3)

        self.fig.tight_layout()
        self.canvas.draw()

    def show_hit(self, h: Hit):
        self.align_text.delete("1.0", "end")

        line1 = []
        line2 = []
        line3 = []

        for a, b in zip(h.aligned_s1, h.aligned_s2):
            line1.append(a)
            line2.append(b)
            if a == b and a != "-":
                line3.append("|")
            else:
                line3.append(" ")

        s1_line = "".join(line1)
        s2_line = "".join(line2)
        match_line = "".join(line3)

        metrics = compute_similarity_metrics(
            h.aligned_s1,
            h.aligned_s2,
            self.current_match,
            self.current_mismatch,
            self.current_gap,
        )

        id_no_gaps = 100.0 * metrics["identity_no_gaps"]
        id_with_gaps = 100.0 * metrics["identity_with_gaps"]
        norm_score = 100.0 * metrics["normalized_score"]
        M = metrics["M"]
        X = metrics["X"]
        G = metrics["G"]
        L = metrics["L"]

        header = (
            f"Score (Smith-Waterman): {h.score}  (window {h.window_idx})\n"
            f"M = {M}, X = {X}, G = {G}, L = {L}\n"
            "Similarity equations:\n"
            "  S1 = M / (M + X)                 (identity without gaps)\n"
            "  S2 = M / (M + X + G)             (identity with gaps)\n"
            "  S3 = (M*match + X*mismatch + G*gap) / (L*match)\n"
            f"Values:\n"
            f"  S1 = {id_no_gaps:.2f}%\n"
            f"  S2 = {id_with_gaps:.2f}%\n"
            f"  S3 = {norm_score:.2f}%\n\n"
        )
        self.align_text.insert("end", header)

        start_line = int(self.align_text.index("end").split(".")[0])

        self.align_text.insert("end", s1_line + "\n")
        self.align_text.insert("end", match_line + "\n")
        self.align_text.insert("end", s2_line + "\n")

        row_s1 = start_line
        row_mid = start_line + 1
        row_s2 = start_line + 2

        for i, (a, b) in enumerate(zip(h.aligned_s1, h.aligned_s2)):
            if a == "-" or b == "-":
                tag = "gap"
            elif a == b:
                tag = "match"
            else:
                tag = "mismatch"

            col_start = i
            col_end = i + 1

            self.align_text.tag_add(tag, f"{row_s1}.{col_start}", f"{row_s1}.{col_end}")
            self.align_text.tag_add(tag, f"{row_mid}.{col_start}", f"{row_mid}.{col_end}")
            self.align_text.tag_add(tag, f"{row_s2}.{col_start}", f"{row_s2}.{col_end}")

        self.update_overview_plot(h)


if __name__ == "__main__":
    app = AlignmentApp()
    app.mainloop()