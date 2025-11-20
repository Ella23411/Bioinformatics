import random


NUCLEOTIDES = "ACGT"

TRANSPOSONS = {
    "T1": "ATGGCATACGTTAGC",
    "T2": "TACGATGCATAG",
    "T3": "CGTCTAGAGGCTA",
    "T4": "GATTCGATTCGA",
}


def random_dna(min_len=200, max_len=400):
    length = random.randint(min_len, max_len)
    
    seq = ""
    for _ in range(length):
        seq += random.choice(NUCLEOTIDES)
    
    return seq


# Insert several transposons at random positions.
def insert_transposons(host_seq, transposons_to_use, allow_overlap=True):
  
    s = list(host_seq)
    placements = {}

    for name, t_seq in transposons_to_use.items():
        t_len = len(t_seq)
        max_start = len(s) - t_len
        pos0 = random.randint(0, max_start)  

        for i, base in enumerate(t_seq):
            s[pos0 + i] = base

        placements[name] = (pos0 + 1, pos0 + t_len)

    return ''.join(s), placements



def find_transposons(seq, transposon_dict):
    
    hits = []
    n = len(seq)

    for name, pattern in transposon_dict.items():
        m = len(pattern)
        for i in range(n - m + 1):
            if seq[i:i+m] == pattern:
                hits.append((name, i + 1, i + m)) 
    return hits



def build_annotation_line(length, hits):
    ann = ['.'] * length
    for idx, (_, start, end) in enumerate(hits, start=1):
        ch = str(idx % 10)
        for i in range(start - 1, end):
            if ann[i] == '.':
                ann[i] = ch
            else:
                ann[i] = '*' 
    return ''.join(ann)


def visualize_sequence(seq, ground_truth_hits, detected_hits, width=80):
    L = len(seq)
    gt_ann = build_annotation_line(L, ground_truth_hits)
    dt_ann = build_annotation_line(L, detected_hits)

    print("\nVISUALIZATION")
    print("Legend: digits = different transposons, * = overlap, . = none")
    print()

    for i in range(0, L, width):
        chunk_seq = seq[i:i+width]
        chunk_gt = gt_ann[i:i+width]
        chunk_dt = dt_ann[i:i+width]
        start_pos = i + 1

        print(f"{start_pos:6d}  {chunk_seq}")
        print("       GT:", chunk_gt)
        print("       DT:", chunk_dt)
        print()


# ==========================
# MAIN EXECUTION
# ==========================

if __name__ == "__main__":
    random.seed(0)

    host = random_dna()

    chosen_names = random.sample(list(TRANSPOSONS.keys()), k=4)
    chosen = {name: TRANSPOSONS[name] for name in chosen_names}

    dna_seq, true_positions_dict = insert_transposons(
        host, chosen, allow_overlap=True
    )

    true_positions = [
        (name, start, end) for name, (start, end) in true_positions_dict.items()
    ]

    detected_positions = find_transposons(dna_seq, chosen)

 
    print("\nDNA SEQUENCE")
    print("Length:", len(dna_seq), "bp")

    print("\nTRUE TRANSPOSON POSITIONS")
    for name, (start, end) in true_positions_dict.items():
        print(f"{name}: {start}–{end}")

    print("\nDETECTED TRANSPOSONS")
    if not detected_positions:
        print("No complete transposons detected.")
    else:
        for name, start, end in detected_positions:
            print(f"{name}: {start}–{end}")

    visualize_sequence(dna_seq, true_positions, detected_positions)