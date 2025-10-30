#!/usr/bin/env python3
import random
from collections import Counter
from typing import List, Tuple
import matplotlib.pyplot as plt


def sample_reads(sequence, num_reads=2000, min_read_len=100, max_read_len=150):
    reads = []
    N = len(sequence)

    # random sampling
    for _ in range(num_reads):
        L = random.randint(min_read_len, max_read_len)
        start = random.randint(0, N - L)
        reads.append((sequence[start:start + L], start))

    # ensure first and last regions are covered
    L0 = random.randint(min_read_len, max_read_len)
    reads.append((sequence[0:L0], 0))

    L1 = random.randint(min_read_len, max_read_len)
    start_last = N - L1
    reads.append((sequence[start_last:start_last + L1], start_last))

    return reads


def reconstruct_sequence(reads: List[Tuple[str, int]], original_length: int) -> str:
    coverage: List[List[str]] = [[] for _ in range(original_length)]

    for read, start in reads:
        for offset, base in enumerate(read):
            idx = start + offset
            if 0 <= idx < original_length:
                coverage[idx].append(base)

    reconstructed = []
    for bases in coverage:
        if bases:
            counts = Counter(bases)
            consensus_base = counts.most_common(1)[0][0]
        else:
            consensus_base = 'N'
        reconstructed.append(consensus_base)

    return ''.join(reconstructed)


def plot_coverage(reads: List[Tuple[str, int]], seq_length: int):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title("Sample Coverage over DNA Sequence")
    ax.set_xlabel("Sample along sequence ")
    ax.set_ylabel("Sample")


    for i, (_, start) in enumerate(reads[:300]):  # limit to 500
        L = len(reads[i][0])
        ax.hlines(y=i, xmin=start, xmax=start + L, color='tab:blue', alpha=0.5)

    ax.set_xlim(0, seq_length)
    ax.set_ylim(-5, min(len(reads[:300]) + 5, 305))
    ax.set_yticks([])
    ax.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.show()


def main() -> None:
    dna_sequence = input("Enter DNA sequence (A/C/G/T only): ").strip().upper()
    orig_len = len(dna_sequence)

    # 1. Sample random reads
    reads = sample_reads(dna_sequence, num_reads=2000, min_read_len=100, max_read_len=150)

    # 2. Plot read coverage
    plot_coverage(reads, orig_len)

    # 3. Reconstruct the sequence
    reconstructed_seq = reconstruct_sequence(reads, orig_len)

    # 4. Report results
    match = reconstructed_seq == dna_sequence
    print(f"Original sequence length: {orig_len}")
    print(f"Reconstructed sequence length: {len(reconstructed_seq)}")
    print(f"Reconstruction successful: {match}")

    preview = 100
    print("\nOriginal (first 100 bases):", dna_sequence[:preview])
    print("Original (last 100 bases):", dna_sequence[-preview:])
    print("\nReconstructed (first 100 bases):", reconstructed_seq[:preview])
    print("Reconstructed (last 100 bases):", reconstructed_seq[-preview:])


if __name__ == "__main__":
    main()
