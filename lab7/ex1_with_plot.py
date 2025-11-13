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


def build_repeat_track(seq_len, repeats, value="pattern_length"):
    track = [0] * seq_len

    for k, repeat_list in repeats.items():
        for r in repeat_list:
            start = r["start"]
            pattern_len = k
            total_len = pattern_len * r["repeat_count"]
            end = start + total_len 

            if value == "pattern_length":
                v = pattern_len
            elif value == "repeat_count":
                v = r["repeat_count"]
            else:
                v = 1

            for i in range(start, min(end, seq_len)):
                track[i] = max(track[i], v)

    return track


import matplotlib.pyplot as plt

def plot_repeat_track(track):
    x = list(range(len(track)))
    plt.figure()
    plt.plot(x, track)
    plt.xlabel("Position in sequence")
    plt.ylabel("Repeat length")
    plt.title("Repeats along the sequence")
    plt.tight_layout()
    plt.show()


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

    track = build_repeat_track(len(seq), repeats, value="pattern_length")
    plot_repeat_track(track)
