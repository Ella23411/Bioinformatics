from collections import Counter

seq = "ATTGCCCCGAAT"
alph = list(dict.fromkeys(seq))         
counts = Counter(seq)
freq = [counts[ch] for ch in alph]     

print(alph)
n = len(seq)
for j in freq:
    print(j / n)
