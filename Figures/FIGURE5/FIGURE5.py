#!/usr/bin/env python3
"""
Influenza Mutation Subway Map (PB2 to NS1, Alternating Labels)
"""

import numpy as np
import matplotlib.pyplot as plt


# ---------- DATA ----------

mutations_data = {
    'PB2': [
        ('T58A', 58), ('V109I', 109), ('V139I', 139),
        ('E362G', 362), ('D441N', 441), ('V495I', 495),
        ('M631L', 631), ('V649I', 649),
        ('K670R', 670), ('T676A', 676)
    ],
    'PB1': [
        ('E75D', 75), ('M171V', 171),
        ('M179I', 179), ('A587P', 587)
    ],
    'PA': [
        ('K113R', 113), ('L219I', 219), ('S277P', 277),
        ('V432I', 432), ('K497R', 497), ('S558L', 558)
    ],
    'HA': [
        ('D104G', 104), ('L131Q', 131),
        ('T211I', 211), ('S336N', 336)
    ],
    'NP': [
        ('Y52H', 52), ('I119V', 119), ('S482N', 482)
    ],
    'NA': [
        ('V67I', 67), ('N71S', 71), ('L269M', 269),
        ('V321I', 321), ('S339P', 339)
    ],
    'M1': [
        ('N82S', 82), ('N85S', 85), ('N87T', 87), ('A227T', 227)
    ],
    'M2': [
        ('D88N', 88)
    ],
    'NS1': [
        ('S7L', 7), ('R67G', 67), ('C116S', 116), ('A223E', 223)
    ]
}

segment_aa_lengths = {
    'PB2': 759,
    'PB1': 757,
    'PA': 716,
    'HA': 566,
    'NP': 498,
    'NA': 469,
    'M1': 252,
    'M2': 97,
    'NS1': 230
}

line_colors = {
    'PB2': '#003049',
    'PB1': '#264653',
    'PA':  '#006D5B',
    'HA':  '#BB3E03',
    'NP':  '#8FB339',
    'NA':  '#BC6C25',
    'M1':  '#6A4C93',
    'M2':  '#6A4C93',
    'NS1': '#006D5B'
}


# ---------- STYLE KNOBS ----------

MM_TO_INCH = 1 / 25.4

FIG_WIDTH_MM = 177.8
FIG_HEIGHT_MM = 210

FIG_WIDTH_INCH = FIG_WIDTH_MM * MM_TO_INCH
FIG_HEIGHT_INCH = FIG_HEIGHT_MM * MM_TO_INCH

LABEL_FONTSIZE = 8
PROTEIN_FONTSIZE = 8
LENGTH_FONTSIZE = 8
PANEL_FONTSIZE = 8

BASE_GAP = 1.2
MIN_INTERVAL_GAP = 0.55
LINE_SPACING = 4.8
LINE_WIDTH = 0.8
STATION_RADIUS = 0.15

MIN_X_SEPARATION = 0.42


# ---------- HELPERS ----------

def estimate_label_width(label: str) -> float:
    return 0.33 * len(label) + 1.45


def assign_tracks_interval(xs, labels, min_gap=0.4):
    if len(xs) == 0:
        return []

    intervals = []
    for i, (x, txt) in enumerate(zip(xs, labels)):
        w = estimate_label_width(txt)
        half_w = w / 2.0
        intervals.append((i, x, x - half_w, x + half_w))

    intervals_sorted = sorted(intervals, key=lambda t: t[1])

    tracks_right_edge = []
    track_for_sorted = [None] * len(intervals_sorted)

    for s_idx, (orig_idx, x, left, right) in enumerate(intervals_sorted):
        placed = False
        for t_idx, t_right in enumerate(tracks_right_edge):
            if left - t_right >= min_gap:
                tracks_right_edge[t_idx] = right
                track_for_sorted[s_idx] = t_idx
                placed = True
                break
        if not placed:
            tracks_right_edge.append(right)
            track_for_sorted[s_idx] = len(tracks_right_edge) - 1

    track_for_original = [None] * len(xs)
    for s_idx, (orig_idx, *_rest) in enumerate(intervals_sorted):
        track_for_original[orig_idx] = track_for_sorted[s_idx]

    return track_for_original


def vertical_offset_alternating(track_idx: int, base_gap: float) -> float:
    level = (track_idx // 2) + 1
    sign = 1 if track_idx % 2 == 0 else -1
    return sign * level * base_gap


def spread_close_positions(xs, min_sep=0.35):
    if len(xs) <= 1:
        return np.array(xs, dtype=float)
    xs = np.array(xs, dtype=float)
    adjusted = xs.copy()
    for i in range(1, len(adjusted)):
        if adjusted[i] - adjusted[i - 1] < min_sep:
            adjusted[i] = adjusted[i - 1] + min_sep
    return adjusted


# ---------- MAIN ----------

def main():
    import matplotlib.font_manager as fm
    fm.fontManager.addfont('/usr/share/fonts/truetype/Arial.ttf')
    fm.fontManager.addfont('/usr/share/fonts/truetype/Arial_Bold.ttf')
    available = [f.name for f in fm.fontManager.ttflist if 'Arial' in f.name]
    if not available:
        raise RuntimeError("Arial not found. Please install Arial before running.")

    plt.rcParams.update({
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 600,
        "savefig.dpi": 600,
        "font.size": 8,
        "font.weight": "bold"
    })

    fig, ax = plt.subplots(figsize=(FIG_WIDTH_INCH, FIG_HEIGHT_INCH))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    line_start_x = 3.0
    max_line_width = 40.0
    start_y = len(mutations_data) * LINE_SPACING
    pb2_length = segment_aa_lengths['PB2']

    top_label_y = start_y
    bottom_label_y = 0

    for i, (protein, mutations) in enumerate(mutations_data.items()):
        y = start_y - (i * LINE_SPACING)
        protein_length = segment_aa_lengths[protein]

        line_width = (protein_length / pb2_length) * max_line_width
        line_end_x = line_start_x + line_width

        ax.plot(
            [line_start_x, line_end_x],
            [y, y],
            color=line_colors[protein],
            linewidth=LINE_WIDTH,
            solid_capstyle='round',
            zorder=1
        )

        ax.text(
            line_start_x - 0.35,
            y,
            protein,
            fontsize=PROTEIN_FONTSIZE,
            fontweight='bold',
            ha='right',
            va='center',
            color='black',
            zorder=6
        )

        ax.text(
            line_end_x + 0.35,
            y,
            f"{protein_length} aa",
            fontsize=LENGTH_FONTSIZE,
            fontweight='bold',
            ha='left',
            va='center',
            color='black',
            zorder=6
        )

        if mutations:
            xs_raw = [
                line_start_x + (pos / protein_length) * line_width
                for (_, pos) in mutations
            ]
            labels = [mut for (mut, _) in mutations]

            xs = spread_close_positions(xs_raw, min_sep=MIN_X_SEPARATION)
            tracks = assign_tracks_interval(xs, labels, min_gap=MIN_INTERVAL_GAP)

            for (mut, pos), x, track_idx in zip(mutations, xs, tracks):
                dy = vertical_offset_alternating(track_idx, BASE_GAP)

                if protein == 'NA':
                    if mut == 'N71S':
                        dy = BASE_GAP * 1.0
                    elif mut == 'V67I':
                        dy = BASE_GAP * 2.0

                if protein == 'NP':
                    if mut == 'I119V':
                        dy = BASE_GAP * 1.0

                station = plt.Circle(
                    (x, y),
                    STATION_RADIUS,
                    facecolor=line_colors[protein],
                    edgecolor='white',
                    linewidth=0.6,
                    zorder=5
                )
                ax.add_patch(station)

                label_y = y + dy
                sign = np.sign(dy) if dy != 0 else 1

                ax.plot(
                    [x, x],
                    [y + sign * (STATION_RADIUS + 0.10),
                     label_y - sign * 0.20],
                    color=line_colors[protein],
                    linewidth=0.6,
                    zorder=3
                )

                ax.text(
                    x,
                    label_y,
                    mut,
                    fontsize=LABEL_FONTSIZE,
                    fontweight='bold',
                    ha='center',
                    va='bottom' if dy > 0 else 'top',
                    color='black',
                    bbox=dict(
                        boxstyle='round,pad=0.22',
                        facecolor='white',
                        edgecolor=line_colors[protein],
                        linewidth=0.6
                    ),
                    zorder=4
                )

                top_label_y = max(top_label_y, label_y)
                bottom_label_y = min(bottom_label_y, label_y)

    # Panel label: uppercase A
    ax.text(
        -0.02,
        1.02,
        "A",
        transform=ax.transAxes,
        fontsize=PANEL_FONTSIZE,
        fontweight="bold",
        va="top",
        ha="left",
        color="black"
    )

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xlim(0, line_start_x + max_line_width + 4.5)
    ax.set_ylim(bottom_label_y - 1.2, top_label_y + 1.2)

    plt.tight_layout(pad=0.4)

    plt.savefig(
        "/content/FIGURE5.pdf",
        
        facecolor="white",
        edgecolor="none"
    )

    plt.savefig(
        "/content/FIGURE5.png",
        dpi=600,
        
        facecolor="white",
        edgecolor="none"
    )

    print("Saved PDF and PNG.")


if __name__ == "__main__":
    main()
