import json
from collections import Counter
import numpy as np

DNA = ["A", "C", "G", "T"]

def order0_probs(seq: str):
    counts = Counter(seq)
    total = len(seq)
    probs = {b: counts.get(b, 0) / total for b in DNA}
    return probs, counts, total

def order1_transition_matrix(seq: str):
    idx = {b: i for i, b in enumerate(DNA)}
    counts = np.zeros((4, 4), dtype=int)
    row_totals = np.zeros(4, dtype=int)

    for a, b in zip(seq[:-1], seq[1:]):
        i = idx[a]
        j = idx[b]
        counts[i, j] += 1
        row_totals[i] += 1

    probs = np.zeros((4, 4), dtype=float)
    for i in range(4):
        if row_totals[i] > 0:
            probs[i, :] = counts[i, :] / row_totals[i]

    return probs, counts, row_totals

def save_transition_json(out_path: str, seq: str, p0: dict, p1: np.ndarray):
    data = {
        "alphabet": DNA,
        "sequence_length": len(seq),
        "order0": p0,
        "order1": {
            DNA[i]: {DNA[j]: float(p1[i, j]) for j in range(4)}
            for i in range(4)
        },
        "order1_matrix": p1.tolist(),
        "matrix_convention": "rows=current, cols=next (A,C,G,T)"
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

if __name__ == "__main__":

    dna_seq = "TGCCATAGGCGTTGAACGCTACACGGACGATACGAATTTACGTATAGAGC"
    dna_seq = "".join(b for b in dna_seq.upper() if b in DNA)

    out_json = "transition_matrix.json"

    p0, c0, n = order0_probs(dna_seq)
    p1, c1, row_totals = order1_transition_matrix(dna_seq)

    save_transition_json(out_json, dna_seq, p0, p1)

    print("Sequence:", dna_seq)
    print("\nOrder-0 probabilities:")
    for b in DNA:
        print(f"P({b}) = {p0[b]:.4f}")

    print("\nOrder-1 transition matrix (rows=current, cols=next) [A C G T]:")
    print(p1)

    print(f"\nSaved JSON to: {out_json}")