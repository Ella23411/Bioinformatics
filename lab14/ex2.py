import re
import math
from collections import Counter

EMINESCU_TEXT = """
Somnoroase păsărele
Pe la cuiburi se adună,
Se ascund în rămurele —
Noapte bună!

Doar izvoarele suspină,
Pe când codrul negru tace;
Dorm și florile-n grădină —
Dormi în pace!
"""

STANESCU_TEXT = """
Leoaică tânără, iubirea
mi-a sărit în faţă.
Mă pândise-n încordare
mai demult.

Colţii albi mi i-a înfipt în faţă,
m-a muşcat, leoaica, azi, de faţă.
Şi deodată-n jurul meu, natura
se făcu un cerc, de-a dura.
"""

ACCUSED_TEXT = """
Somnoroase păsărele pe la cuiburi se adună.
Leoaică tânără, iubirea mi-a sărit în faţă.
În jurul meu natura se făcu un cerc de-a dura.
Noapte bună, dormi în pace.
Iubirea pândise mai demult.
"""

WINDOW_WORDS = 18      
STEP_WORDS = 6       
ALPHA = 0.2            
THRESH = 0.35          


def normalize_text(text: str) -> str:
    text = text.lower()
    text = re.sub(r"[^a-zăâîșşțţ\- \n]", " ", text, flags=re.IGNORECASE)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def tokenize_words(text: str):
    text = normalize_text(text)
    return [w for w in text.split(" ") if w]


class BigramModel:
    def __init__(self, alpha=0.2):
        self.alpha = float(alpha)
        self.unigrams = Counter()
        self.bigrams = Counter()
        self.vocab = set()

    def fit(self, tokens):
        self.unigrams.update(tokens)
        self.vocab.update(tokens)
        for a, b in zip(tokens[:-1], tokens[1:]):
            self.bigrams[(a, b)] += 1

    def prob(self, prev_word, next_word, vocab_union):
        V = max(1, len(vocab_union))
        num = self.bigrams[(prev_word, next_word)] + self.alpha
        den = self.unigrams[prev_word] + self.alpha * V
        return num / den


def llr_beta(prev_word, next_word, model_E, model_S, vocab_union):
    pE = model_E.prob(prev_word, next_word, vocab_union)
    pS = model_S.prob(prev_word, next_word, vocab_union)
    return math.log(pE / pS, 2)


def score_window(tokens, start, window_words, model_E, model_S, vocab_union):
    end = min(len(tokens), start + window_words)
    w = tokens[start:end]
    if len(w) < 2:
        return 0.0, start, end

    s = 0.0
    n = 0
    for a, b in zip(w[:-1], w[1:]):
        s += llr_beta(a, b, model_E, model_S, vocab_union)
        n += 1

    return (s / n), start, end  


def label_from_score(score, thresh):
    if score > thresh:
        return "EMINESCU"
    if score < -thresh:
        return "STANESCU"
    return "NEITHER"


def sliding_attribution(text, model_E, model_S, vocab_union, title):
    tokens = tokenize_words(text)
    if len(tokens) < 2:
        print(f"\n=== {title} ===")
        print("Text too short.")
        return

    results = []
    i = 0
    while i < len(tokens) - 1:
        score, s, e = score_window(tokens, i, WINDOW_WORDS, model_E, model_S, vocab_union)
        lab = label_from_score(score, THRESH)
        results.append((s, e, score, lab))
        i += STEP_WORDS

    counts = Counter([lab for _, _, _, lab in results])
    total = sum(counts.values()) if results else 1
    def pct(x): return 100.0 * x / total

    print(f"\n=== {title} ===")
    print(f"Tokens: {len(tokens)} | windows: {len(results)}")
    print(f"Threshold: {THRESH} (avg LLR per bigram)")
    print("Window label distribution:")
    for lab in ["EMINESCU", "STANESCU", "NEITHER"]:
        print(f"  {lab:<9}: {counts.get(lab, 0):>3}  ({pct(counts.get(lab, 0)):.1f}%)")

    print("\nFirst windows (preview):")
    for (s, e, score, lab) in results[:min(6, len(results))]:
        snippet = " ".join(tokens[s:e])
        print(f"  [{s:>4}:{e:<4}] score={score:+.3f} -> {lab}")
        print(f"    {snippet}")

    if title.lower().startswith("sanity: only eminescu"):
        print("\nSanity expectation: mostly EMINESCU.")
    if title.lower().startswith("sanity: only stanescu"):
        print("\nSanity expectation: mostly STANESCU.")


def main():
    tok_E = tokenize_words(EMINESCU_TEXT)
    tok_S = tokenize_words(STANESCU_TEXT)

    model_E = BigramModel(alpha=ALPHA)
    model_S = BigramModel(alpha=ALPHA)
    model_E.fit(tok_E)
    model_S.fit(tok_S)

    vocab_union = set(model_E.vocab) | set(model_S.vocab)

    print("=== Training summary ===")
    print(f"Eminescu tokens: {len(tok_E)} | vocab: {len(model_E.vocab)}")
    print(f"Stănescu tokens: {len(tok_S)} | vocab: {len(model_S.vocab)}")
    print(f"Union vocab: {len(vocab_union)}")

    sliding_attribution(EMINESCU_TEXT, model_E, model_S, vocab_union, "Sanity: ONLY Eminescu text (should be Eminescu)")
    sliding_attribution(STANESCU_TEXT, model_E, model_S, vocab_union, "Sanity: ONLY Stanescu text (should be Stanescu)")

    sliding_attribution(ACCUSED_TEXT, model_E, model_S, vocab_union, "Accused text (mixed / unknown)")


if __name__ == "__main__":
    main()