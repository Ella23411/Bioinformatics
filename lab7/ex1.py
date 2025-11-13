def read_fasta(path):
    seq_lines = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq_lines.append(line.strip())
    return "".join(seq_lines).upper()


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


if __name__ == "__main__":
    fasta_file = "file.fasta"

    seq = read_fasta(fasta_file)
    print("Sequence length:", len(seq))

    repeats = find_tandem_repeats(seq, min_len=3, max_len=6, min_repeats=2)

    total_patterns = sum(len(v) for v in repeats.values())
    print("\nTotal repeat patterns found:", total_patterns)

    for k in sorted(repeats.keys()):
        print(f"\nRepeats with length {k}:")
        for r in repeats[k].items() if isinstance(repeats[k], dict) else repeats[k]:
            if isinstance(r, dict):
                pattern = r["pattern"]
                start = r["start"]
                repeat_count = r["repeat_count"]
            else:
                pattern = r.get("pattern")
                start = r.get("start")
                repeat_count = r.get("repeat_count")

            print(f"  pattern={pattern}, start={start}, repeats={repeat_count}")