def counts_existing(seq: str, k: int) -> dict:
    seq = seq.upper()
    out = {}
    n = len(seq)
    if n < k:
        return out
    for i in range(n - k + 1):
        aux = seq[i:i+k]
        out[aux] = out.get(aux, 0) + 1
    return out

def percentages_existing(seq: str, k: int) -> dict:
    cnt = counts_existing(seq, k)
    total = max(len(seq) - k + 1, 0)
    return {aux: (c / total) * 100.0 for aux, c in cnt.items()} if total else {}

def main():
    S = "ABAA"
    for k in (2, 3):
        cnt = counts_existing(S, k)
        pct = percentages_existing(S, k)
        print(f"S = {S!r} | k = {k}")
        print("-" * 30)
        for km in sorted(cnt): 
            print(f"{km}\tcount={cnt[km]}   percent={pct[km]:6.2f}%")
        print()

  
    S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
    for k in (2, 3):
        cnt = counts_existing(S, k)
        pct = percentages_existing(S, k)
        print(f"S = {S!r}| k = {k}")
        print("-" * 30)
        for km in sorted(cnt): 
            print(f"{km}\tcount={cnt[km]}   percent={pct[km]:6.2f}%")
        print()

if __name__ == "__main__":
    main()
