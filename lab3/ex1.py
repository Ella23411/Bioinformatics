import math

# Tm = 4(G + C) + 2(A + T) 
def calculate_tm_1(dna_sequence: str) -> float:
  
    dna = dna_sequence.upper().strip()
    
    A = dna.count('A')
    T = dna.count('T')
    G = dna.count('G')
    C = dna.count('C')

    tm = 4 * (G + C) + 2 * (A + T)

    return tm

# Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length
def calculate_tm_2(dna_sequence: str, na_concentration: float) -> float:

    dna = dna_sequence.upper().strip()
    
    length = len(dna)
    G = dna.count('G')
    C = dna.count('C')
    gc_percent = ((G + C) / length) * 100

    tm = -1 * (81.5 + 16.6 * math.log10(na_concentration) + 0.41 * gc_percent - (600 / length))

    return tm



if __name__ == "__main__":
    # TGAGCT
    dna_input = input("Enter DNA sequence: ").strip()
    tm_value_1 = calculate_tm_1(dna_input)
    tm_value_2 = calculate_tm_2(dna_input, 0.001)
    print(f"Melting Temperature using formula 1: {tm_value_1:.2f} ")
    print(f"Melting Temperature using formula 2: {tm_value_2:.2f} ")

