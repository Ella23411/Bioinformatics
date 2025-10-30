import random
import time
from typing import List, Tuple
import matplotlib.pyplot as plt


USER_SEQUENCES: List[str] = [
    # "TCTGGCCCAGTTCCTAGGTAGTGACGAATTCGTGGTGGTGA",
    # "CAGACGGGCGATTTTGTTAAAGCCACCTTCTTTAGTCAAATTCTCA",
    # "AGAGTTGAAAGGCACATTTGGTTGTTTATGGGTAGAATTCGATCTG",
    # "TGCAGCATTGTTAGCAGGAT",
    # "TACGGGCACCCCGGGCCGGCGGAAGGGTCCCGCTCT",
    # "CAGCCTGCGTAGACGGTC",
    # "TGCGCGTTGATGTGAGG",
    # "GCAGTTGCAGTGCTGCTG",
    # "TGTTCCTCCCTCGTAGGTT",
    # "TGCAGGGACTCAGCCAGCGGTAGGGGATGGGCTTAGGC"
]

DNA_ALPH = "ACGT"

def clean_dna(s: str) -> str:
    return "".join([c for c in s.upper() if c in DNA_ALPH])

def random_dna(n: int) -> str:
    return "".join(random.choice(DNA_ALPH) for _ in range(n))

def make_10_random(min_len=1000, max_len=3000) -> List[str]:
    return [random_dna(random.randint(min_len, max_len)) for _ in range(10)]

def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)

def sample_reads(seq: str, n_reads: int = 2000, read_len_range: Tuple[int,int]=(100,150)) -> List[str]:
    reads = []
    L = len(seq)
    for _ in range(n_reads):
        rlen = random.randint(*read_len_range)
        if rlen > L: rlen = L
        start = random.randint(0, L - rlen)
        reads.append(seq[start:start+rlen])
    return reads

def overlap(a: str, b: str, min_ov: int) -> int:
    start = 0
    while True:
        start = a.find(b[:min_ov], start)
        if start == -1:
            return 0
        ov = len(a) - start
        if b.startswith(a[start:]):
            return ov
        start += 1

def greedy_scs(reads: List[str], min_overlap: int = 20) -> str:
    reads = list(dict.fromkeys(reads))  
    reads.sort(key=len, reverse=True)
    pruned = []
    for r in reads:
        if not any(r in q for q in pruned[:300]): 
            pruned.append(r)
    reads = pruned

    while len(reads) > 1:
        best_i, best_j, best_ov = -1, -1, -1
        for i in range(len(reads)):
            ai = reads[i]
            for j in range(len(reads)):
                if i == j: 
                    continue
                ov = overlap(ai, reads[j], min_overlap)
                if ov > best_ov:
                    best_i, best_j, best_ov = i, j, ov
        if best_ov < min_overlap:
            reads[0] = reads[0] + reads[-1]
            reads.pop()
        else:
            i, j = best_i, best_j
            merged = reads[i] + reads[j][best_ov:]
            new_reads = []
            for k, r in enumerate(reads):
                if k not in (i, j):
                    new_reads.append(r)
            new_reads.append(merged)
            reads = new_reads
    return reads[0]

def assemble_with_cap(reads: List[str], min_overlap: int = 20, max_reads: int = 400) -> str:
    sorted_reads = sorted(reads, key=len, reverse=True)
    core = sorted_reads[:max_reads]
    if len(sorted_reads) > max_reads:
        pick = min(max_reads // 4, len(sorted_reads) - max_reads)
        if pick > 0:
            core += random.sample(sorted_reads[max_reads:], k=pick)
    return greedy_scs(core, min_overlap=min_overlap)

def assemble_benchmark(
    sequences: List[str] = None,
    n_reads: int = 2000,
    read_len_range: Tuple[int,int] = (100,150),
    min_overlap: int = 20,
    max_reads_for_assembly: int = 400,
    seed: int = 42
):
    random.seed(seed)

    if not sequences or len(sequences) == 0:
        sequences = make_10_random()
    else:
        sequences = [clean_dna(s) for s in sequences]
        if len(sequences) != 10:
            raise ValueError(f"Please provide exactly 10 sequences (got {len(sequences)}).")

    gc_list = []
    time_ms_list = []

    for idx, seq in enumerate(sequences, start=1):
        gc = gc_percent(seq)
        reads = sample_reads(seq, n_reads=n_reads, read_len_range=read_len_range)
        t0 = time.perf_counter()
        _assembled = assemble_with_cap(reads, min_overlap=min_overlap, max_reads=max_reads_for_assembly)
        t1 = time.perf_counter()
        elapsed_ms = (t1 - t0) * 1000.0

        print(f"Virus {idx:02d}: length={len(seq)} | GC={gc:.2f}% | assembly_time={elapsed_ms:.2f} ms")
        gc_list.append(gc)
        time_ms_list.append(elapsed_ms)


    plt.figure()
    plt.scatter(gc_list, time_ms_list)
    for i, (x, y) in enumerate(zip(gc_list, time_ms_list), start=1):
        plt.annotate(str(i), (x, y), xytext=(3, 3), textcoords="offset points")
    plt.xlabel("Overall GC% (per sequence)")
    plt.ylabel("Assembly time (ms)")
    plt.title("Assembly Time vs GC% (10 sequences)")
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.savefig("assembly_vs_gc.png", dpi=150)
    print("\nSaved plot: assembly_vs_gc.png")

    return {
        "gc_percentages": gc_list,
        "assembly_times_ms": time_ms_list
    }

if __name__ == "__main__":
    assemble_benchmark(USER_SEQUENCES)
