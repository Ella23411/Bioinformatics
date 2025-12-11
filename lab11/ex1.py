import tkinter as tk
from tkinter import ttk, messagebox

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle


def needleman_wunsch(s1, s2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    s1 = s1.upper()
    s2 = s2.upper()
    n = len(s1)
    m = len(s2)

    # scoring matrix
    score = np.zeros((n + 1, m + 1), dtype=int)
    # traceback matrix: 0=diag, 1=up, 2=left
    traceback = np.zeros((n + 1, m + 1), dtype=int)

    # init first row/col
    for i in range(1, n + 1):
        score[i, 0] = score[i - 1, 0] + gap_penalty
        traceback[i, 0] = 1  # up
    for j in range(1, m + 1):
        score[0, j] = score[0, j - 1] + gap_penalty
        traceback[0, j] = 2  # left

    # fill matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if s1[i - 1] == s2[j - 1]:
                diag = score[i - 1, j - 1] + match_score
            else:
                diag = score[i - 1, j - 1] + mismatch_score

            up = score[i - 1, j] + gap_penalty
            left = score[i, j - 1] + gap_penalty

            best = max(diag, up, left)
            score[i, j] = best

            if best == diag:
                traceback[i, j] = 0
            elif best == up:
                traceback[i, j] = 1
            else:
                traceback[i, j] = 2

    # traceback path and alignment strings
    i, j = n, m
    aligned1 = []
    aligned2 = []
    path_cells = []

    while i > 0 or j > 0:
        path_cells.append((i, j))
        tb = traceback[i, j]
        if i > 0 and j > 0 and tb == 0:
            aligned1.append(s1[i - 1])
            aligned2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or tb == 1):
            aligned1.append(s1[i - 1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(s2[j - 1])
            j -= 1

    path_cells.append((0, 0))
    path_cells = path_cells[::-1]

    aligned1 = ''.join(aligned1[::-1])
    aligned2 = ''.join(aligned2[::-1])

    return score, traceback, aligned1, aligned2, path_cells


class NWGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Needleman–Wunsch DNA Aligner")

        # --- left panel: inputs ---
        left = ttk.Frame(self, padding=5)
        left.grid(row=0, column=0, sticky="nsw")

        ttk.Label(left, text="Sequences").grid(row=0, column=0, sticky="w")

        ttk.Label(left, text="Sq1:").grid(row=1, column=0, sticky="w")
        self.seq1_var = tk.StringVar(value="ACCGTGAAGCCAATAC")
        e1 = ttk.Entry(left, textvariable=self.seq1_var, width=25)
        e1.grid(row=2, column=0, sticky="we", pady=(0, 5))

        ttk.Label(left, text="Sq2:").grid(row=3, column=0, sticky="w")
        self.seq2_var = tk.StringVar(value="AGCGTGCAGCCAATAC")
        e2 = ttk.Entry(left, textvariable=self.seq2_var, width=25)
        e2.grid(row=4, column=0, sticky="we", pady=(0, 10))

        # parameters
        ttk.Label(left, text="Parameters").grid(row=5, column=0, sticky="w")

        params = ttk.Frame(left)
        params.grid(row=6, column=0, sticky="we", pady=(2, 10))

        ttk.Label(params, text="Gap:").grid(row=0, column=0, sticky="e")
        self.gap_var = tk.IntVar(value=-1)
        ttk.Entry(params, textvariable=self.gap_var, width=5).grid(
            row=0, column=1, padx=2
        )

        ttk.Label(params, text="Match:").grid(row=1, column=0, sticky="e")
        self.match_var = tk.IntVar(value=1)
        ttk.Entry(params, textvariable=self.match_var, width=5).grid(
            row=1, column=1, padx=2
        )

        ttk.Label(params, text="Mismatch:").grid(row=2, column=0, sticky="e")
        self.mismatch_var = tk.IntVar(value=-1)
        ttk.Entry(params, textvariable=self.mismatch_var, width=5).grid(
            row=2, column=1, padx=2
        )

        ttk.Button(left, text="Align", command=self.run_alignment).grid(
            row=7, column=0, sticky="we", pady=(0, 10)
        )

        # --- right panel: plots + text ---
        right = ttk.Frame(self, padding=5)
        right.grid(row=0, column=1, sticky="nsew")
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        fig = Figure(figsize=(6, 3), dpi=100)
        self.ax_matrix = fig.add_subplot(1, 2, 1)
        self.ax_path = fig.add_subplot(1, 2, 2)

        self.canvas = FigureCanvasTkAgg(fig, master=right)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        right.rowconfigure(0, weight=1)
        right.columnconfigure(0, weight=1)

        self.text = tk.Text(right, height=10, wrap="word")
        self.text.grid(row=1, column=0, sticky="nsew", pady=(5, 0))
        right.rowconfigure(1, weight=0)

    def run_alignment(self):
        s1 = self.seq1_var.get().strip()
        s2 = self.seq2_var.get().strip()

        if not s1 or not s2:
            messagebox.showerror("Error", "Both sequences must be non-empty.")
            return

        try:
            gap = int(self.gap_var.get())
            match = int(self.match_var.get())
            mismatch = int(self.mismatch_var.get())
        except ValueError:
            messagebox.showerror("Error", "Parameters must be integers.")
            return

        score, tb, a1, a2, path = needleman_wunsch(
            s1, s2, match_score=match, mismatch_score=mismatch, gap_penalty=gap
        )

        # --- update matrix heatmap ---
        self.ax_matrix.clear()
        self.ax_matrix.imshow(score, origin="upper")
        self.ax_matrix.set_title("Score matrix")
        self.ax_matrix.set_xlabel("Seq2")
        self.ax_matrix.set_ylabel("Seq1")
        self.ax_matrix.set_xticks(range(len(s2) + 1))
        self.ax_matrix.set_yticks(range(len(s1) + 1))

        # --- update traceback grid ---
        self.ax_path.clear()
        n, m = score.shape
        grid = np.zeros_like(score, dtype=float)
        self.ax_path.imshow(grid, origin="upper")

        # draw grid lines
        for i in range(n + 1):
            self.ax_path.axhline(i - 0.5, linewidth=0.5, color="black")
        for j in range(m + 1):
            self.ax_path.axvline(j - 0.5, linewidth=0.5, color="black")

        # color traceback path cells
        for (i, j) in path:
            rect = Rectangle((j - 0.5, i - 0.5), 1, 1,
                             facecolor="red", alpha=0.7)
            self.ax_path.add_patch(rect)

        self.ax_path.set_xlim(-0.5, m - 0.5)
        self.ax_path.set_ylim(n - 0.5, -0.5)
        self.ax_path.set_title("Traceback path")

        self.canvas.draw()

        # --- update text output ---
        self.text.delete("1.0", tk.END)
        self.show_alignment_text(a1, a2, score)

    def show_alignment_text(self, a1, a2, score):
        # alignment statistics
        matches = 0
        length = len(a1)
        match_line = []

        for c1, c2 in zip(a1, a2):
            if c1 == c2 and c1 != "-":
                matches += 1
                match_line.append("|")
            else:
                match_line.append(" ")
        match_line = "".join(match_line)

        similarity = 100.0 * matches / length if length > 0 else 0.0

        self.text.insert(tk.END, "Show Alignment:\n")
        self.text.insert(tk.END, a1 + "\n")
        self.text.insert(tk.END, match_line + "\n")
        self.text.insert(tk.END, a2 + "\n\n")
        self.text.insert(tk.END, f"Matches = {matches}\n")
        self.text.insert(tk.END, f"Length = {length}\n")
        self.text.insert(tk.END, f"Similarity = {similarity:.1f} %\n")
        self.text.insert(tk.END, f"Final score = {score[-1, -1]}\n")


if __name__ == "__main__":
    app = NWGui()
    app.mainloop()