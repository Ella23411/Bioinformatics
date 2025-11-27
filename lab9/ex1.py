import re
import textwrap
from collections import OrderedDict
import math

import matplotlib.pyplot as plt
import matplotlib


ENZYME_DEFINITIONS = OrderedDict({
    "EcoRI":  ("GAATTC", 1),
    "BamHI":  ("GGATCC", 1),
    "HindIII":("AAGCTT", 1),
    "TaqI":   ("TCGA", 1),
    "HaeIII": ("GGCC", 2),
})


def find_cut_positions(sequence, pattern, cut_index):
    seq = re.sub(r"\s+", "", sequence).upper()
    pat = pattern.upper()
    positions = []
    start = 0
    while True:
        idx = seq.find(pat, start)
        if idx == -1:
            break
        cut_pos = idx + cut_index
        positions.append(cut_pos)
        start = idx + 1  
    return positions


def compute_fragments(seq_length, cut_positions):
    if not cut_positions:
        return [seq_length]
    cut_positions = sorted(cut_positions)
    fragments = []
    prev = 0
    for cut in cut_positions:
        frag_len = cut - prev
        if frag_len > 0:
            fragments.append(frag_len)
        prev = cut
    if seq_length > prev:
        fragments.append(seq_length - prev)
    return fragments

def simulate_gel(fragment_map, output_file="gel_simulation.png"):
    if not fragment_map:
        print("No fragments to simulate.")
        return

    all_frags = [f for frags in fragment_map.values() for f in frags if f > 0]
    if not all_frags:
        print("Invalid fragment lengths for simulation.")
        return

    max_len = max(all_frags)
    min_len = min(all_frags)

    if max_len == min_len:
        max_len *= 1.1
        min_len *= 0.9

    log_max = math.log10(max_len)
    log_min = math.log10(min_len)

    def frag_to_y(frag_len):
        lf = math.log10(frag_len)
        y = (log_max - lf) / (log_max - log_min)
        return 0.05 + 0.9 * y

    num_lanes = len(fragment_map)

    fig_height = 8
    fig_width = 2.0 * num_lanes
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    ax.set_facecolor("white")

    for lane_idx, (enzyme, fragments) in enumerate(fragment_map.items()):
        x_center = lane_idx   
        lane_half_width = 0.3

        well_y = 0.98
        ax.hlines(
            well_y,
            x_center - lane_half_width,
            x_center + lane_half_width,
            linewidth=4,
            color="grey"
        )

        for frag_len in fragments:
            y_center = frag_to_y(frag_len)
            ax.hlines(
                y_center,
                x_center - lane_half_width,
                x_center + lane_half_width,
                linewidth=2,
                color="black"
            )

        ax.text(
            x_center,
            1.02,
            enzyme,
            ha="center",
            va="bottom",
            fontsize=9
        )

    label_frags = sorted(set(all_frags), reverse=True)
    x_label = -0.8
    min_gap = 0.03  

    last_y = None
    for frag_len in label_frags:
        y_center = frag_to_y(frag_len)
        if last_y is not None and abs(y_center - last_y) < min_gap:
            continue
        ax.text(
            x_label,
            y_center,
            f"{frag_len} bp",
            ha="right",
            va="center",
            fontsize=9
        )
        last_y = y_center

    ax.set_xlim(-1.0, num_lanes - 0.5)
    ax.set_ylim(0, 1.05)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, facecolor=fig.get_facecolor())
    print(f"Gel simulation saved to '{output_file}'.")


def report_digest(enzyme_name, fragment_lengths, cut_positions):
    num_cuts = len(cut_positions)
    print(f"\n=== {enzyme_name} Digest ===")
    if num_cuts == 0:
        print("No cleavage sites detected in the provided DNA sequence.")
        print(f"Fragment size: {fragment_lengths[0]} bp (undigested)")
    else:
        print(f"Number of cleavage sites: {num_cuts}")
        pos_str = ', '.join(str(p + 1) for p in cut_positions)
        print(f"Cut positions (1-based): {pos_str}")
        frag_str = ', '.join(str(fl) for fl in fragment_lengths)
        print(f"Fragment lengths (bp): {frag_str}")


def get_dna_sequence_from_fasta(filepath="file.fasta"):
    try:
        with open(filepath, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise ValueError(f"FASTA file '{filepath}' not found.")

    seq_raw = "".join(
        line.strip()
        for line in lines
        if not line.strip().startswith(">")
    )

    seq = re.sub(r"[^ACGTacgt]", "", seq_raw)

    if not seq:
        raise ValueError("No valid DNA sequence found in FASTA file.")

    return seq


def get_selected_enzymes_from_user():
    available = list(ENZYME_DEFINITIONS.keys())
    print("Available enzymes:")
    for idx, name in enumerate(available, 1):
        pattern, cut_idx = ENZYME_DEFINITIONS[name]
        print(f"  {idx}. {name}\t(Recognition: {pattern}, cut at index {cut_idx + 1} of the pattern)")

    selection = input(
        "Enter the enzyme names separated by commas (or type 'all' to use all enzymes):\n"
    ).strip()
    if not selection:
        raise ValueError("No enzyme selection provided.")

    if selection.lower() == 'all':
        return available

    chosen = [s.strip() for s in selection.split(',') if s.strip()]

    invalid = [c for c in chosen if c not in ENZYME_DEFINITIONS]
    if invalid:
        raise ValueError(f"Invalid enzyme name(s): {', '.join(invalid)}")
    return chosen


def main():
    print("Restriction Enzyme Digest Simulator")
    print("This tool uses EcoRI, BamHI, HindIII, TaqI, and HaeIII by default.")
    print("DNA sequence will be loaded from 'file.fasta' in the current directory.\n")

    try:
        seq = get_dna_sequence_from_fasta("file.fasta")
        print("Successfully loaded DNA sequence from 'file.fasta'.")
    except ValueError as ve:
        print(ve)
        return

    try:
        selected_enzymes = get_selected_enzymes_from_user()
    except ValueError as ve:
        print(ve)
        return

    seq_len = len(seq)
    print(f"\nSequence length: {seq_len} bp")

    digest_results = OrderedDict()
    for enzyme in selected_enzymes:
        pattern, cut_index = ENZYME_DEFINITIONS[enzyme]
        cut_positions = find_cut_positions(seq, pattern, cut_index)
        fragments = compute_fragments(seq_len, cut_positions)
        digest_results[enzyme] = {
            'cuts': cut_positions,
            'fragments': fragments,
        }
        report_digest(enzyme, fragments, cut_positions)

    do_gel = input("\nWould you like to simulate an electrophoresis gel? (y/n): ").strip().lower()
    if do_gel.startswith('y'):
        fragment_map = {ename: info['fragments'] for ename, info in digest_results.items()}
        simulate_gel(fragment_map, output_file="gel_simulation.png")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\nSimulation aborted by user.")