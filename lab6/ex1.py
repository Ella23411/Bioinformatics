from pathlib import Path
from typing import List, Tuple
import random
import matplotlib.pyplot as plt


FASTA_PATH = Path(r"C:\Users\Ella\Desktop\Bioinformatics\lab6\file.fasta")
RANDOM_SEED = None
N_FRAGMENTS = 10
MIN_BP = 100
MAX_BP = 3000
LADDER_BANDS = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]
LADDER_LABELS = {3000: "3000 bp -", 1500: "1500 bp -", 500: "500 bp -"}
FIGSIZE = (4, 7)
BACKGROUND = "black"
BAND_COLOR = "white"
LANE_WIDTH = 0.28
LADDER_X = 0.30
SAMPLE_X = 0.70   
WELL_TOP_BP = 11000
WELL_HEIGHT_BP = 1000
OUT_PNG = Path("screenshot_1.png")



def read_fasta(path: Path) -> str:
    try:
        text = path.read_text(encoding="utf-8")
    except Exception as e:
        print(f"Error reading file: {e}. Using an empty sequence.")
        return ""

    seq_lines = [ln.strip() for ln in text.splitlines() if not ln.startswith(">")]
    seq = "".join(seq_lines).upper()
    allowed = set("ACGTN")
    seq = "".join(ch for ch in seq if ch in allowed)
    return seq


def sample_fragments(seq: str,
                     n: int = N_FRAGMENTS,
                     min_bp: int = MIN_BP,
                     max_bp: int = MAX_BP,
                     seed=None) -> Tuple[List[str], List[int]]:
    rng = random.Random(seed)
    if not seq:
        lengths = [rng.randint(min_bp, max_bp) for _ in range(n)]
        return [], sorted(lengths, reverse=True)

    max_len = min(max_bp, len(seq))
    min_len = min(min_bp, max_len)

    fragments = []
    for _ in range(n):
        L = rng.randint(min_len, max_len)
        start = 0 if L >= len(seq) else rng.randint(0, len(seq) - L)
        fragments.append(seq[start:start + L])

    lengths = [len(f) for f in fragments]
    lengths.sort(reverse=True)
    return fragments, lengths


def setup_axes(ax):
    ax.set_facecolor(BACKGROUND)
    ax.set_xlim(0, 1)
    y_min, y_max = 50, max(LADDER_BANDS) * 1.1
    ax.set_ylim(y_min, y_max)
    ax.set_yscale("log")
    ax.invert_yaxis()

    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])


def draw_well(ax, center_x: float):
    rect = plt.Rectangle(
        (center_x - LANE_WIDTH / 2, WELL_TOP_BP),
        LANE_WIDTH,
        WELL_HEIGHT_BP,
        facecolor="#555555",
        edgecolor="#888888"
    )
    ax.add_patch(rect)


def draw_lane_outline(ax, center_x: float):
    ax.vlines([center_x - LANE_WIDTH / 2, center_x + LANE_WIDTH / 2],
              ymin=60, ymax=max(LADDER_BANDS), colors="#666666", linewidth=0.8)


def draw_bands(ax, center_x: float, bands: List[int], heavy: set = None):
    heavy = heavy or set()
    x0, x1 = center_x - LANE_WIDTH / 2, center_x + LANE_WIDTH / 2
    for bp in bands:
        lw = 4 if bp in heavy else 2.5
        ax.plot([x0, x1], [bp, bp], color=BAND_COLOR, linewidth=lw, solid_capstyle="butt")


def draw_labels(ax, center_x: float, labels: dict):
    for bp, text in labels.items():
        ax.text(center_x - LANE_WIDTH / 2 - 0.05, bp, text,
                color="white", fontsize=12, ha="right", va="center")


def render_gel(fragment_lengths: List[int]):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    fig.patch.set_facecolor(BACKGROUND)
    setup_axes(ax)

    for cx in (LADDER_X, SAMPLE_X):
        draw_lane_outline(ax, cx)
        draw_well(ax, cx)

   
    draw_bands(ax, LADDER_X, LADDER_BANDS, heavy={3000, 1000})
    draw_labels(ax, LADDER_X, LADDER_LABELS)

  
    draw_bands(ax, SAMPLE_X, fragment_lengths)

    ax.set_title("Simulated Gel Electrophoresis", color="white", pad=8)

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)
    plt.close(fig)


def main():
    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)

    seq = read_fasta(FASTA_PATH)
    genome_len = len(seq)
    if genome_len == 0:
        print("Warning: genome sequence empty or unreadable. Proceeding with random lengths only.")
        genome_len = 30000
    print(f"Parsed genome. Total length: {genome_len} bp")

    fragments, fragment_lengths = sample_fragments(seq, n=N_FRAGMENTS,
                                                  min_bp=MIN_BP, max_bp=MAX_BP,
                                                  seed=RANDOM_SEED)
    print(f"Generated {len(fragment_lengths)} fragment lengths (bp): {fragment_lengths}")

    render_gel(fragment_lengths)
    print(f"Generated image: {OUT_PNG}")


if __name__ == "__main__":
    main()