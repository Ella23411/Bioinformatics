from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import pandas as pd
import matplotlib.pyplot as plt


ALPHABET = ["A", "C", "G", "T"]

MOTIFS: List[str] = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT",
]


@dataclass
class MotifModel:
    count: pd.DataFrame
    weight: pd.DataFrame
    rel_freq: pd.DataFrame
    log_likelihood: pd.DataFrame
    background: Dict[str, float]
    pseudocount: float


def build_model(
    motifs: List[str],
    pseudocount: float = 1.0,
    background: Optional[Dict[str, float]] = None,
) -> MotifModel:
    L = len(motifs[0])
    if any(len(m) != L for m in motifs):
        raise ValueError("All motifs must have the same length")

    if background is None:
        background = {b: 0.25 for b in ALPHABET}

    # 1) Count matrix
    count = pd.DataFrame(0, index=ALPHABET, columns=list(range(1, L + 1)), dtype=float)
    for m in motifs:
        for pos, base in enumerate(m.upper(), start=1):
            count.loc[base, pos] += 1.0

    # 2) Weight matrix
    weight = count + pseudocount

    # 3) Relative frequencies
    rel = weight / weight.sum(axis=0)

    # 4) Log-likelihood matrix
    ll = pd.DataFrame(index=ALPHABET, columns=rel.columns)
    for b in ALPHABET:
        ll.loc[b, :] = (rel.loc[b, :] / background[b]).apply(math.log)

    return MotifModel(count, weight, rel, ll, background, pseudocount)


def score_window(window: str, ll: pd.DataFrame) -> float:
    score = 0.0
    for pos, base in enumerate(window.upper(), start=1):
        if base not in ALPHABET:
            return float("nan")
        score += ll.loc[base, pos]
    return score


def scan_scores(seq: str, ll: pd.DataFrame) -> Tuple[List[int], List[float]]:
    L = ll.shape[1]
    positions, scores = [], []
    for i in range(len(seq) - L + 1):
        positions.append(i + 1)
        scores.append(score_window(seq[i:i+L], ll))
    return positions, scores


def motif_scores(motifs: List[str], ll: pd.DataFrame) -> List[float]:
    return [score_window(m, ll) for m in motifs]


def read_fasta(path: str) -> List[Tuple[str, str]]:
    records = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            records.append((header, "".join(seq)))
    return records


def plot_all_genomes_signals(
    results: List[Tuple[str, List[int], List[float]]],
    threshold: float,
) -> None:
    n = len(results)
    ncols = 2
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(16, 2.2 * nrows), 
        sharey=True,
    )

    axes = axes.flatten()

    for i, (ax, (_, positions, scores)) in enumerate(zip(axes, results), start=1):
        ax.plot(positions, scores, linewidth=1.0)
        ax.axhline(threshold, linestyle="--", linewidth=1.0)
        ax.set_title(f"Influenza {i}", fontsize=11)
        ax.grid(True)
        ax.set_ylabel("")

        ax.tick_params(axis="both", labelsize=9)

    for j in range(len(results), len(axes)):
        axes[j].axis("off")

    fig.text(
        0.02,
        0.5,
        "Log-likelihood",
        va="center",
        rotation=90,
        fontsize=13,
    )

    plt.subplots_adjust(
        wspace=0.25,
        hspace=0.9,   
    )

    plt.tight_layout(rect=[0.06, 0.03, 0.98, 0.97])
    plt.show()



def main():
    fasta_path = "influenza.fasta"

    model = build_model(MOTIFS)

    train_scores = motif_scores(MOTIFS, model.log_likelihood)
    threshold = min(train_scores)

    records = read_fasta(fasta_path)

    results_for_plot = []

    for name, seq in records:
        positions, scores = scan_scores(seq, model.log_likelihood)
        results_for_plot.append((name, positions, scores))

    plot_all_genomes_signals(results_for_plot, threshold)


if __name__ == "__main__":
    main()
