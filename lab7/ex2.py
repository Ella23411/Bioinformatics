import matplotlib.pyplot as plt

def read_fasta_multi(path):
    sequences = []
    header = None
    seq_lines = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None and seq_lines:
                    sequences.append((header, "".join(seq_lines).upper()))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)

    # last one
    if header is not None and seq_lines:
        sequences.append((header, "".join(seq_lines).upper()))

    return sequences


def find_tandem_repeats(seq, min_len=6, max_len=10, min_repeats=2):
    seq = seq.replace("\n", "").replace(" ", "").upper()
    n = len(seq)
    results = {}

    for k in range(min_len, max_len + 1):
        repeats_for_k = []
        if 2 * k > n:
            continue

        i = 0
        while i <= n - 2 * k:
            pattern = seq[i:i + k]
            count = 1
            j = i + k

            while j + k <= n and seq[j:j + k] == pattern:
                count += 1
                j += k

            if count >= min_repeats:
                repeats_for_k.append({
                    "pattern": pattern,
                    "start": i,
                    "repeat_count": count
                })
                i = j
            else:
                i += 1

        if repeats_for_k:
            results[k] = repeats_for_k

    return results


def most_frequent_repeat(repeats):
    best = None
    for k, repeat_list in repeats.items():
        for r in repeat_list:
            rc = r["repeat_count"]
            if best is None:
                best = {
                    "pattern": r["pattern"],
                    "pattern_len": k,
                    "start": r["start"],
                    "repeat_count": rc
                }
            else:
                if (rc > best["repeat_count"] or
                    (rc == best["repeat_count"] and k > best["pattern_len"])):

                    best = {
                        "pattern": r["pattern"],
                        "pattern_len": k,
                        "start": r["start"],
                        "repeat_count": rc
                    }
    return best


if __name__ == "__main__":
    fasta_file = "sequences.fasta"

    seqs = read_fasta_multi(fasta_file)
    print(f"Loaded {len(seqs)} sequences from {fasta_file}")

    labels = []
    max_repeats_values = []
    best_info_per_seq = []

    for idx, (header, seq) in enumerate(seqs, start=1):
        print(f"\nSequence {idx} ({header}): length={len(seq)}")
        repeats = find_tandem_repeats(seq, min_len=3, max_len=6, min_repeats=2)
        best = most_frequent_repeat(repeats)

        if best is None:
            print("  No tandem repeats found.")
            max_rep = 0
        else:
            print(f"  Most frequent repeat:")
            print(f"    pattern      = {best['pattern']}")
            print(f"    pattern_len  = {best['pattern_len']}")
            print(f"    repeat_count = {best['repeat_count']}")
            print(f"    start        = {best['start']}")
            max_rep = best["repeat_count"]

        labels.append(f"Seq{idx}") 
        max_repeats_values.append(max_rep)
        best_info_per_seq.append(best)

    x = range(len(labels))

    plt.figure()
    bars = plt.bar(x, max_repeats_values)
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.ylabel("Max repeat count (most frequent repeat)")
    plt.xlabel("Sequence")
    plt.title("Most frequent repetition per DNA sequence")


    for i, best in enumerate(best_info_per_seq):
        if best is None:
            text = "None"
        else:
            text = best["pattern"]

        height = max_repeats_values[i]

        plt.text(
            i,                       
            height + 0.01,            
            text,                    
            ha='center',
            va='bottom',
            fontsize=8,
            rotation=0              #
        )

    plt.tight_layout()
    plt.show()