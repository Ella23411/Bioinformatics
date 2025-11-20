import sys
from dataclasses import dataclass
from typing import List, Dict, Tuple


MIN_IR_LEN = 10
MAX_IR_LEN = 40
MAX_IR_MISMATCHES = 2

MIN_TE_LEN = 1000
MAX_TE_LEN = 20000

MIN_SPACER = 0
MAX_SPACER = 20000



def read_fasta(path: str):
    header = None
    seq_chunks = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        yield header, "".join(seq_chunks)



def rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        max_len = max(len(a), len(b))
        min_len = min(len(a), len(b))
        diff = sum(1 for x, y in zip(a[:min_len], b[:min_len]) if x != y)
        return diff + (max_len - min_len)
    return sum(1 for x, y in zip(a, b) if x != y)


@dataclass
class TransposonHit:
    genome_id: str
    start: int
    end: int
    length: int
    left_ir_start: int
    left_ir_end: int
    right_ir_start: int
    right_ir_end: int
    ir_length: int


def detect_transposons_in_sequence(
    genome_id: str,
    seq: str,
    min_ir_len: int = MIN_IR_LEN,
    max_ir_len: int = MAX_IR_LEN,
    max_ir_mismatches: int = MAX_IR_MISMATCHES,
    min_te_len: int = MIN_TE_LEN,
    max_te_len: int = MAX_TE_LEN,
    min_spacer: int = MIN_SPACER,
    max_spacer: int = MAX_SPACER
) -> List[TransposonHit]:

    seq = seq.upper()
    n = len(seq)
    hits = []

    for left_start in range(0, n - 2 * min_ir_len):
        for ir_len in range(min_ir_len, max_ir_len + 1):

            left_end = left_start + ir_len
            if left_end >= n:
                break

            left_ir = seq[left_start:left_end]
            expected = rev_comp(left_ir)

            min_spacer_here = max(min_spacer, min_te_len - 2 * ir_len)
            max_spacer_here = min(max_spacer, max_te_len - 2 * ir_len)

            if min_spacer_here > max_spacer_here:
                continue

            for spacer in range(min_spacer_here, max_spacer_here + 1):

                right_start = left_end + spacer
                right_end = right_start + ir_len
                if right_end > n:
                    break

                right_ir = seq[right_start:right_end]
                mism = hamming_distance(right_ir, expected)

                if mism <= max_ir_mismatches:
                    hits.append(
                        TransposonHit(
                            genome_id=genome_id,
                            start=left_start,
                            end=right_end,
                            length=right_end - left_start,
                            left_ir_start=left_start,
                            left_ir_end=left_end,
                            right_ir_start=right_start,
                            right_ir_end=right_end,
                            ir_length=ir_len
                        )
                    )
    return hits



def classify_relations(
    hits: List[TransposonHit]
):
    hits_sorted = sorted(hits, key=lambda h: (h.genome_id, h.start, h.end))
    n = len(hits_sorted)

    overlaps = {i: [] for i in range(n)}
    nested_in = {i: [] for i in range(n)}

    for i in range(n):
        hi = hits_sorted[i]
        for j in range(i + 1, n):
            hj = hits_sorted[j]

            if hj.genome_id != hi.genome_id:
                continue

            a1, a2 = hi.start, hi.end
            b1, b2 = hj.start, hj.end

            if b1 >= a2:
                break

            if max(a1, b1) < min(a2, b2):
                if a1 <= b1 and b2 <= a2:
                    nested_in[j].append(i)
                elif b1 <= a1 and a2 <= b2:
                    nested_in[i].append(j)
                else:
                    overlaps[i].append(j)
                    overlaps[j].append(i)

    return hits_sorted, overlaps, nested_in


def write_hits(
    hits,
    overlaps,
    nested_in,
    out_path="results.tsv"
):
    with open(out_path, "w") as out:
        out.write(
            "index\tgenome_id\tstart\tend\tlength\tleft_ir_start\tleft_ir_end\t"
            "right_ir_start\tright_ir_end\tir_length\toverlaps_with\tnested_in\n"
        )

        for i, h in enumerate(hits):
            o = ",".join(map(str, overlaps[i])) if overlaps[i] else "-"
            n = ",".join(map(str, nested_in[i])) if nested_in[i] else "-"

            out.write(
                f"{i}\t{h.genome_id}\t{h.start+1}\t{h.end}\t{h.length}\t"
                f"{h.left_ir_start+1}\t{h.left_ir_end}\t"
                f"{h.right_ir_start+1}\t{h.right_ir_end}\t"
                f"{h.ir_length}\t{o}\t{n}\n"
            )



def main():

    fasta_file = "file.fasta"  
    all_hits = []

    sys.stderr.write(f"Loading genomes from {fasta_file}...\n")

    for header, seq in read_fasta(fasta_file):
        sys.stderr.write(f"Processing genome: {header} (length={len(seq)})...\n")
        hits = detect_transposons_in_sequence(header, seq)
        sys.stderr.write(f"  Found {len(hits)} transposons.\n")
        all_hits.extend(hits)

    hits_sorted, overlaps, nested_in = classify_relations(all_hits)

    sys.stderr.write("Writing results to results.tsv...\n")
    write_hits(hits_sorted, overlaps, nested_in)

    sys.stderr.write("Done!\n")


if __name__ == "__main__":
    main()