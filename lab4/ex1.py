GENETIC_CODE = {
    "UUU":"Phe","UUC":"Phe","UUA":"Leu","UUG":"Leu",
    "CUU":"Leu","CUC":"Leu","CUA":"Leu","CUG":"Leu",
    "AUU":"Ile","AUC":"Ile","AUA":"Ile","AUG":"Met",
    "GUU":"Val","GUC":"Val","GUA":"Val","GUG":"Val",
    "UCU":"Ser","UCC":"Ser","UCA":"Ser","UCG":"Ser",
    "CCU":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "ACU":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "GCU":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "UAU":"Tyr","UAC":"Tyr","UAA":"Stop","UAG":"Stop",
    "CAU":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "AAU":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
    "GAU":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "UGU":"Cys","UGC":"Cys","UGA":"Stop","UGG":"Trp",
    "CGU":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "AGU":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
    "GGU":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly"
}

def translate_sequence(seq: str) -> tuple[str, str]:
   
    clean_seq = ""
    for ch in seq.upper():
        if ch in "ACGTU":
            if ch == "T":
                clean_seq += "U"
            else:
                clean_seq += ch
    seq = clean_seq


    start = seq.find("AUG")
    if start == -1:
        return seq, "No start codon found."

    protein = []
    for i in range(start, len(seq) - 2, 3):
        codon = seq[i:i+3]
        amino_acid = GENETIC_CODE.get(codon, "")
        if amino_acid == "Stop":
            break
        if amino_acid == "":
            break
        protein.append(amino_acid)

    if not protein:
        return seq, "No protein could be translated."

    return seq, "-".join(protein)


# --- MAIN ---
sequence = input("Enter your DNA or RNA sequence: ")

rna, protein = translate_sequence(sequence)

print("\nRNA sequence:")
print(rna)

print("\nAmino acid sequence:")
print(protein)
