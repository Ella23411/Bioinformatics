from pathlib import Path
from typing import List, Tuple
import random
import math
import matplotlib.pyplot as plt

FASTA_MULTI = Path(r"C:\Users\Ella\Desktop\Bioinformatics\lab6\sequences.fasta")
RANDOM_SEED = None
N_FRAGMENTS = 10
MIN_BP = 100
MAX_BP = 3000
CAL_BP_MIN = 100
CAL_BP_MAX = 10000
LADDER_BANDS = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]
LADDER_LABELS = {3000: "3000 bp -", 1500: "1500 bp -", 500: "500 bp -"}
BACKGROUND = "black"
FIGSIZE_ALL = (8, 7)
FIGSIZE_SINGLE = (4, 7)
LANE_WIDTH = 0.08
WELL_TOP_BP = 11000
WELL_HEIGHT_BP = 1000
OUT_ALL = Path("screenshot_2.png")
OUT_BEST = Path("screenshot_3.png")
VISIBLE_SEP = 0.02


def read_multi_fasta(path: Path) -> List[str]:
    try:
        text = path.read_text(encoding="utf-8")
    except:
        return []
    sequences = []
    current = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                seq = "".join(current).upper()
                seq = "".join(ch for ch in seq if ch in "ACGTN")
                sequences.append(seq)
                current = []
        else:
            current.append(line)
    if current:
        seq = "".join(current).upper()
        seq = "".join(ch for ch in seq if ch in "ACGTN")
        sequences.append(seq)
    return sequences


def sample_fragments(seq: str, n=N_FRAGMENTS, min_bp=MIN_BP, max_bp=MAX_BP, seed=None) -> Tuple[List[str], List[int]]:
    rng = random.Random(seed)
    if not seq:
        lengths = sorted([rng.randint(min_bp, max_bp) for _ in range(n)], reverse=True)
        return [], lengths
    max_len = min(max_bp, len(seq))
    min_len = min(min_bp, max_len)
    frags = []
    for _ in range(n):
        L = rng.randint(min_len, max_len)
        start = 0 if L >= len(seq) else rng.randint(0, len(seq)-L)
        frags.append(seq[start:start+L])
    lengths = sorted((len(f) for f in frags), reverse=True)
    return frags, lengths


def norm_distance_log(bp, bp_min=CAL_BP_MIN, bp_max=CAL_BP_MAX):
    lmin, lmax = math.log10(max(bp_min, 50)), math.log10(bp_max)
    x = (math.log10(bp) - lmin) / (lmax - lmin + 1e-9)
    return 1.0 - max(0.0, min(1.0, x))


def pos(bp):
    return norm_distance_log(bp)


def distinct_visible_count(lengths: List[int], sep=VISIBLE_SEP) -> int:
    if not lengths:
        return 0
    ys = sorted((pos(bp) for bp in lengths))
    count = 0
    last = None
    for y in ys:
        if last is None or (y - last) >= sep:
            count += 1
            last = y
    return count


def setup_axes(ax):
    ax.set_facecolor(BACKGROUND)
    y_min, y_max = 50, CAL_BP_MAX * 1.1
    ax.set_xlim(0, 1)
    ax.set_ylim(y_min, y_max)
    ax.set_yscale("log")
    ax.invert_yaxis()
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])


def draw_well(ax, cx):
    r = plt.Rectangle((cx - LANE_WIDTH/2, WELL_TOP_BP), LANE_WIDTH, WELL_HEIGHT_BP, facecolor="#555555", edgecolor="#888888")
    ax.add_patch(r)


def draw_lane_outline(ax, cx):
    ax.vlines([cx - LANE_WIDTH/2, cx + LANE_WIDTH/2], ymin=60, ymax=CAL_BP_MAX, colors="#666666", linewidth=0.8)


def draw_bands(ax, cx, bands, color="white", heavy=None, lw=2.5):
    heavy = heavy or set()
    x0, x1 = cx - LANE_WIDTH/2, cx + LANE_WIDTH/2
    for bp in bands:
        ax.plot([x0, x1], [bp, bp], color=color, linewidth=(4 if bp in heavy else lw), solid_capstyle="butt")


def draw_labels(ax, cx, labels):
    for bp, text in labels.items():
        ax.text(cx - LANE_WIDTH*0.7, bp, text, color="white", fontsize=11, ha="right", va="center")


def render_gel_all(all_lengths: List[List[int]]):
    fig, ax = plt.subplots(figsize=FIGSIZE_ALL)
    fig.patch.set_facecolor(BACKGROUND)
    setup_axes(ax)
    ladder_x = 0.08
    lane_centers = [0.18 + i * (0.78 / max(1, len(all_lengths)-1)) for i in range(len(all_lengths))]
    draw_lane_outline(ax, ladder_x)
    draw_well(ax, ladder_x)
    draw_bands(ax, ladder_x, LADDER_BANDS, heavy={3000,1000})
    draw_labels(ax, ladder_x, LADDER_LABELS)
    cmap = plt.get_cmap("tab10")
    for i, (cx, lens) in enumerate(zip(lane_centers, all_lengths)):
        c = cmap(i % 10)
        draw_lane_outline(ax, cx)
        draw_well(ax, cx)
        draw_bands(ax, cx, lens, color=c, lw=3)
    ax.set_title("All genomes", color="white", pad=8)
    fig.tight_layout()
    fig.savefig(OUT_ALL, dpi=160)
    plt.close(fig)


def render_gel_single(lengths: List[int]):
    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    fig.patch.set_facecolor(BACKGROUND)
    setup_axes(ax)
    ladder_x = 0.30
    sample_x = 0.70
    draw_lane_outline(ax, ladder_x)
    draw_well(ax, ladder_x)
    draw_bands(ax, ladder_x, LADDER_BANDS, heavy={3000,1000})
    draw_labels(ax, ladder_x, LADDER_LABELS)
    draw_lane_outline(ax, sample_x)
    draw_well(ax, sample_x)
    draw_bands(ax, sample_x, lengths, lw=3)
    ax.set_title("Best genome (most DNA segments)", color="white", pad=8)
    fig.tight_layout()
    fig.savefig(OUT_BEST, dpi=170)
    plt.close(fig)


def main():
    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)

    seqs = read_multi_fasta(FASTA_MULTI)[:10]
    all_lengths = []
    for seq in seqs:
        _, lens = sample_fragments(seq, seed=RANDOM_SEED)
        all_lengths.append(lens)

    render_gel_all(all_lengths)

    visible_counts = [distinct_visible_count(lens) for lens in all_lengths]
    best_idx = max(range(len(visible_counts)), key=lambda i: visible_counts[i])

    render_gel_single(all_lengths[best_idx])


if __name__ == "__main__":
    main()