import json
import re
from collections import Counter, defaultdict

def tokenize_words(txt: str):
    return re.findall(r"[A-Za-z]+(?:'[A-Za-z]+)?", txt.lower())

def build_word_symbol_map(words):
    uniq = sorted(set(words))
    word_to_sym = {w: f"w{i}" for i, w in enumerate(uniq)}
    sym_to_word = {v: k for k, v in word_to_sym.items()}
    return word_to_sym, sym_to_word

def order1_transition_probs(symbols):
    trans_counts = defaultdict(Counter)
    row_totals = Counter()

    for a, b in zip(symbols[:-1], symbols[1:]):
        trans_counts[a][b] += 1
        row_totals[a] += 1

    trans_probs = {}
    for a, counter in trans_counts.items():
        total = row_totals[a]
        trans_probs[a] = {b: counter[b] / total for b in counter}

    return trans_probs, trans_counts, row_totals

def save_transition_json(out_path, txt, words, word_to_sym, sym_to_word, trans_probs):
    data = {
        "text_length_chars": len(txt),
        "num_words": len(words),
        "unique_words": len(word_to_sym),
        "word_to_symbol": word_to_sym,
        "symbol_to_word": sym_to_word,
        "matrix_convention": "P(next_word | current_word) using symbols",
        "order1": trans_probs
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

if __name__ == "__main__":
    txt = (
        "Yesterday, I walked to the river and watched the light fold across the water. "
        "A dog barked twice; someone laughed, then the street went quiet. "
        "I wrote a note to myself: keep going, even when plans change."
    )

    words = tokenize_words(txt)
    if len(words) < 2:
        raise ValueError("Need at least 2 words to compute transitions.")

    word_to_sym, sym_to_word = build_word_symbol_map(words)
    symbols = [word_to_sym[w] for w in words]

    trans_probs, trans_counts, row_totals = order1_transition_probs(symbols)

    out_json = "word_transition_matrix.json"
    save_transition_json(out_json, txt, words, word_to_sym, sym_to_word, trans_probs)

    print("Text:")
    print(txt)
    print("\nFirst 20 tokens:")
    print(words[:20])
    print(f"\nSaved JSON to: {out_json}")