from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import pandas as pd


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

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"


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
    background: Dict[str, float] | None = None,
) -> MotifModel:
    if not motifs:
        raise ValueError("motifs must not be empty")

    L = len(motifs[0])
    if any(len(m) != L for m in motifs):
        raise ValueError("All motifs must have the same length")

    if background is None:
        background = {b: 0.25 for b in ALPHABET}  # uniform null model

    # 1) Count matrix
    count = pd.DataFrame(0, index=ALPHABET, columns=list(range(1, L + 1)), dtype=float)
    for m in motifs:
        m = m.upper().strip()
        for pos, base in enumerate(m, start=1):
            if base not in ALPHABET:
                raise ValueError(f"Invalid base '{base}' in motif '{m}'")
            count.loc[base, pos] += 1.0

    # 2) Weight matrix
    weight = count + float(pseudocount)

    # 3) Relative frequencies matrix 
    rel = weight / weight.sum(axis=0)

    # 4) Log-likelihood matrix: ln(P(base|motif) / P(base|null))
    ll = pd.DataFrame(index=ALPHABET, columns=rel.columns, dtype=float)
    for b in ALPHABET:
        p0 = float(background[b])
        if p0 <= 0.0:
            raise ValueError("Background probabilities must be > 0")
        ll.loc[b, :] = (rel.loc[b, :] / p0).apply(lambda x: math.log(x))

    return MotifModel(
        count=count,
        weight=weight,
        rel_freq=rel,
        log_likelihood=ll,
        background=background,
        pseudocount=float(pseudocount),
    )


def score_window(window: str, ll: pd.DataFrame) -> float:
    """Sum per-position log-likelihoods for a window."""
    window = window.upper()
    if len(window) != ll.shape[1]:
        raise ValueError("Window length must match motif length")
    score = 0.0
    for pos, base in enumerate(window, start=1):
        if base not in ALPHABET:
            # If you have Ns, you can choose to skip, penalize, or average.
            raise ValueError(f"Invalid base '{base}' in window '{window}'")
        score += float(ll.loc[base, pos])
    return score


def scan_sequence(seq: str, ll: pd.DataFrame) -> List[Tuple[int, str, float]]:
    """Return (1-based start, window, score) for each sliding window."""
    seq = seq.upper()
    L = ll.shape[1]
    out: List[Tuple[int, str, float]] = []
    for i in range(0, len(seq) - L + 1):
        w = seq[i : i + L]
        out.append((i + 1, w, score_window(w, ll)))
    return out


def motif_scores(motifs: List[str], ll: pd.DataFrame) -> List[float]:
    return [score_window(m, ll) for m in motifs]


def main() -> None:
    model = build_model(MOTIFS, pseudocount=1.0, background={b: 0.25 for b in ALPHABET})

    print("\n1) COUNT MATRIX")
    print(model.count.to_string())

    print("\n2) WEIGHT MATRIX (count + pseudocount, pseudocount = %.3g)" % model.pseudocount)
    print(model.weight.to_string())

    print("\n3) RELATIVE FREQUENCIES MATRIX (column-normalized)")
    print(model.rel_freq.round(4).to_string())

    print("\n4) LOG-LIKELIHOODS MATRIX  ln(P(base|motif)/P(base|null)), null=0.25")
    print(model.log_likelihood.round(6).to_string())

    # 5) Scan S
    print("\n5) SLIDING WINDOW SCORES FOR S")
    results = scan_sequence(S, model.log_likelihood)

    train = motif_scores(MOTIFS, model.log_likelihood)
    min_train = min(train)
    mean_train = sum(train) / len(train)

    print(f"\nTraining motif scores: min={min_train:.6f}, mean={mean_train:.6f}, max={max(train):.6f}")
    print("Heuristic hit rule: score >= min(training motif score)")

    hits = []
    for start, window, score in results:
        flag = " <== HIT" if score >= min_train else ""
        if flag:
            hits.append((start, window, score))
        print(f"{start:>2d}  {window}  {score:> .6f}{flag}")

    print("\nTop 5 windows by score:")
    for start, window, score in sorted(results, key=lambda x: x[2], reverse=True)[:5]:
        print(f"  start={start:>2d}  {window}  score={score: .6f}")

    if hits:
        print("\nConclusion: YES — there are strong signals of an exon–intron boundary-like motif in S.")
        print("Hit start positions (1-based): " + ", ".join(str(h[0]) for h in hits))
    else:
        print("\nConclusion: No strong signals found in S with this model/threshold.")


if __name__ == "__main__":
    main()
