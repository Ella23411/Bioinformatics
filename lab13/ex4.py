from __future__ import annotations

import csv
import json
import random
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union


@dataclass
class Jump:
    step: int
    from_state: str
    to_state: str
    prob: float
    rng: float


def _ensure_model_dict(loaded: Any) -> Dict[str, Any]:
    if isinstance(loaded, dict):
        return loaded
    if isinstance(loaded, list):
        if not loaded or not isinstance(loaded[0], dict):
            raise ValueError("JSON root is a list but does not contain an object model at index 0.")
        return loaded[0]
    raise ValueError(f"Unsupported JSON root type: {type(loaded).__name__}")


def _normalize_dist(dist: Dict[str, float], alphabet: Sequence[str]) -> List[Tuple[str, float]]:
    weights = [max(0.0, float(dist.get(s, 0.0))) for s in alphabet]
    total = sum(weights)
    if total <= 0:
        raise ValueError("Distribution has zero total weight.")
    return [(alphabet[i], weights[i] / total) for i in range(len(alphabet))]


def _normalize_row(row: Sequence[float], alphabet: Sequence[str]) -> List[Tuple[str, float]]:
    weights = [max(0.0, float(x)) for x in row]
    total = sum(weights)
    if total <= 0:
        return []
    return [(alphabet[i], weights[i] / total) for i in range(len(alphabet))]


def _sample(rng: random.Random, probs: List[Tuple[str, float]]) -> Tuple[str, float, float]:
    r = rng.random()
    cum = 0.0
    for sym, p in probs:
        cum += p
        if r <= cum:
            return sym, p, r
    sym, p = probs[-1]
    return sym, p, r


def _validate_square_matrix(m: List[List[float]], n: int) -> None:
    if len(m) != n:
        raise ValueError(f"Transition matrix must have {n} rows; got {len(m)}.")
    for i, row in enumerate(m):
        if len(row) != n:
            raise ValueError(f"Transition matrix row {i} must have {n} cols; got {len(row)}.")


def _write_jumps_csv(path: Path, jumps: List[Jump]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["step", "from_state", "to_state", "prob", "rng"])
        w.writeheader()
        for j in jumps:
            w.writerow(asdict(j))


class MarkovEngine:
    def __init__(
        self,
        alphabet: List[str],
        order0: Dict[str, float],
        order1_matrix: List[List[float]],
        seed: Optional[int] = None,
    ):
        self.alphabet = alphabet
        self.index = {s: i for i, s in enumerate(alphabet)}
        self.rng = random.Random(seed)
        self.p0 = _normalize_dist(order0, alphabet)
        _validate_square_matrix(order1_matrix, len(alphabet))
        self.matrix = order1_matrix

    def pick_start(self) -> str:
        s, _, _ = _sample(self.rng, self.p0)
        return s

    def jump(self, current: str) -> Tuple[Optional[str], Optional[Jump]]:
        if current not in self.index:
            raise ValueError(f"Unknown state '{current}'")
        row = self.matrix[self.index[current]]
        probs = _normalize_row(row, self.alphabet)
        if not probs:
            return None, None
        nxt, p, r = _sample(self.rng, probs)
        return nxt, Jump(step=-1, from_state=current, to_state=nxt, prob=p, rng=r)

    def generate_tokens(self, length: int, start: Optional[str] = None) -> Tuple[List[str], List[Jump]]:
        if length < 1:
            raise ValueError("length must be >= 1")
        current = start if start is not None else self.pick_start()
        if current not in self.index:
            raise ValueError(f"Start state '{current}' not in alphabet")
        tokens = [current]
        jumps: List[Jump] = []
        for step in range(length - 1):
            nxt, j = self.jump(current)
            if nxt is None or j is None:
                break
            j.step = step
            jumps.append(j)
            tokens.append(nxt)
            current = nxt
        return tokens, jumps

    def synthesize(self, length: int, join: str = "", start: Optional[str] = None) -> Tuple[str, List[Jump]]:
        tokens, jumps = self.generate_tokens(length=length, start=start)
        return join.join(tokens), jumps


def load_engine_from_json(path: Union[str, Path], seed: Optional[int] = None) -> Tuple[MarkovEngine, Dict[str, Any]]:
    p = Path(path)
    loaded = json.loads(p.read_text(encoding="utf-8"))
    data = _ensure_model_dict(loaded)
    alphabet = data["alphabet"]
    order0 = data["order0"]
    order1_matrix = data["order1_matrix"]
    engine = MarkovEngine(alphabet=alphabet, order0=order0, order1_matrix=order1_matrix, seed=seed)
    return engine, data


engine, data = load_engine_from_json("transition_matrix.json", seed=42)

sequence_length = int(data.get("sequence_length", 50))
mode = data.get("mode", "dna")

if mode == "text":
    joiner = data.get("joiner", " ")
    seq, jumps = engine.synthesize(sequence_length, join=joiner)
else:
    seq, jumps = engine.synthesize(sequence_length, join="")

Path("jumps.json").write_text(json.dumps([asdict(j) for j in jumps], indent=2), encoding="utf-8")
_write_jumps_csv(Path("jumps.csv"), jumps)