from itertools import product

S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA".upper()
ALPHABET = ("A", "C", "G", "T")

def all_kmers(k: int):
    return ["".join(p) for p in product(ALPHABET, repeat=k)]

def kmer_counts(sequence: str, k: int):
    seq = sequence.upper().strip().replace(" ", "").replace("\n", "")
    total = max(0, len(seq) - k + 1)
    counts = {kmer: 0 for kmer in all_kmers(k)}
    for i in range(total):
        kmer = seq[i:i+k]
        if len(kmer) == k and set(kmer).issubset(ALPHABET):
            counts[kmer] += 1
    return counts, total

def print_report(sequence: str, k: int):
    counts, total = kmer_counts(sequence, k)
    print(f"\n{k}-mer report")
    print(f"Sequence length: {len(sequence)} | total {k}-mer windows: {total}")
    print("Combination\tCount\tPercentage")
    for kmer in sorted(counts):
        pct = (counts[kmer] / total * 100.0) if total > 0 else 0.0
        print(f"{kmer}\t\t{counts[kmer]}\t{pct:.4f}%")

if __name__ == "__main__":
    print_report(S, 2)
    print_report(S, 3)
