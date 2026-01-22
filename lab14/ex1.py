import math

seq = "CAGGTTGGAAACGTAA"


S1 = "ATCGATTCGATATCATACACGTAT"       
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"      

BASES = ["A", "C", "G", "T"]
IDX = {b: i for i, b in enumerate(BASES)}


def clean_dna(s: str) -> str:
    s = s.upper().strip()
    for ch in s:
        if ch not in IDX:
            raise ValueError(f"Invalid DNA base '{ch}' in sequence: {s}")
    if len(s) < 2:
        raise ValueError("Sequence must have length >= 2 to compute transitions.")
    return s


def count_transitions(s: str):
    s = clean_dna(s)
    counts = [[0 for _ in range(4)] for _ in range(4)]
    for a, b in zip(s[:-1], s[1:]):
        counts[IDX[a]][IDX[b]] += 1
    return counts


def normalize_rows(counts):
    probs = [[0.0 for _ in range(4)] for _ in range(4)]
    for i in range(4):
        row_sum = sum(counts[i])
        if row_sum == 0:
            continue
        for j in range(4):
            probs[i][j] = counts[i][j] / row_sum
    return probs


def llr_matrix(p_plus, p_minus):
    beta = [[0.0 for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            a = p_plus[i][j]
            b = p_minus[i][j]
            if a == 0.0 or b == 0.0:
                beta[i][j] = 0.0
            else:
                beta[i][j] = math.log(a / b, 2)
    return beta


def score_sequence(s: str, beta):
    s = clean_dna(s)
    total = 0.0
    for a, b in zip(s[:-1], s[1:]):
        total += beta[IDX[a]][IDX[b]]
    return total


def print_matrix(name, M, fmt="{:>8.3f}"):
    print(f"\n{name}")
    print("      " + "".join([f"{b:>8}" for b in BASES]))
    for i, r in enumerate(M):
        print(f"{BASES[i]:>4}  " + "".join(fmt.format(x) for x in r))


def main():
    c_plus = count_transitions(S1)
    c_minus = count_transitions(S2)

    p_plus = normalize_rows(c_plus)
    p_minus = normalize_rows(c_minus)

    beta = llr_matrix(p_plus, p_minus)

    llr_score = score_sequence(seq, beta)
    decision = "CpG ISLAND (CpG+)" if llr_score > 0 else "NON-ISLAND (CpG-)"

    print("=== TRAINING SEQUENCES ===")
    print("S1 (CpG+):", S1)
    print("S2 (CpG-):", S2)

    print_matrix("Transition probabilities P(+):", p_plus)
    print_matrix("Transition probabilities P(-):", p_minus)
    print_matrix("Log-likelihood ratio beta = log2(P(+)/P(-)):", beta, fmt="{:>8.3f}")

    print("\n=== TEST SEQUENCE ===")
    print("seq =", seq)
    print(f"LLR score (sum of beta over transitions) = {llr_score:.3f}")
    print("Decision:", decision)


if __name__ == "__main__":
    main()