from collections import Counter
import matplotlib.pyplot as plt
import re

GENETIC_CODE = {
    "TTT":"Phe","TTC":"Phe","TTA":"Leu","TTG":"Leu",
    "CTT":"Leu","CTC":"Leu","CTA":"Leu","CTG":"Leu",
    "ATT":"Ile","ATC":"Ile","ATA":"Ile","ATG":"Met",
    "GTT":"Val","GTC":"Val","GTA":"Val","GTG":"Val",
    "TCT":"Ser","TCC":"Ser","TCA":"Ser","TCG":"Ser",
    "CCT":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "ACT":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "GCT":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "TAT":"Tyr","TAC":"Tyr","TAA":"Stop","TAG":"Stop",
    "CAT":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "AAT":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
    "GAT":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "TGT":"Cys","TGC":"Cys","TGA":"Stop","TGG":"Trp",
    "CGT":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "AGT":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
    "GGT":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly"
}


def read_fasta(path):
    seq = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    seq = seq.replace("U","T")
    seq = re.sub(r"[^ACGT]","",seq)
    return seq

def get_codons(seq):
    codons = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            codons.append(codon)
    return codons

def plot_top10(counter, title, filename):
    top10 = counter.most_common(10)
    codons = [c for c,_ in top10]
    counts = [n for _,n in top10]
    plt.bar(codons, counts)
    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()
    plt.close()

def codon_to_aa_counts(codon_counts):
    aa = Counter()
    for c, v in codon_counts.items():
        name = GENETIC_CODE.get(c)
        if name and name != "Stop":
            aa[name] += v
    return aa

if __name__ == "__main__":
    f1 = r"D:\bioinformatics\lab4\covid19.fasta"
    f2 = r"D:\bioinformatics\lab4\influenza.fasta"

    s1 = read_fasta(f1)
    s2 = read_fasta(f2)

    c1 = Counter(get_codons(s1))
    c2 = Counter(get_codons(s2))

    plot_top10(c1, "Top 10 codons (COVID19)", "COVID19.png")
    plot_top10(c2, "Top 10 codons (INFLUENZA)", "INFLUENZA.png")

    print("\nTop 10 codons for COVID19:")
    for c,v in c1.most_common(10):
        print(f"{c}: {v}")
    print("\nTop 10 codons for INFLUENZA:")
    for c,v in c2.most_common(10):
        print(f"{c}: {v}")

    set1 = {c for c,_ in c1.most_common(10)}
    set2 = {c for c,_ in c2.most_common(10)}
    common = set1 & set2
    print("\nCommon codons in both Top 10 lists:", ", ".join(common) if common else "None")

    aa1 = codon_to_aa_counts(c1)
    aa2 = codon_to_aa_counts(c2)

    print("\nTop 3 amino acids (COVID19):")
    for a,v in aa1.most_common(3):
        print(f"{a}: {v}")
    print("\nTop 3 amino acids (INFLUENZA):")
    for a,v in aa2.most_common(3):
        print(f"{a}: {v}")

