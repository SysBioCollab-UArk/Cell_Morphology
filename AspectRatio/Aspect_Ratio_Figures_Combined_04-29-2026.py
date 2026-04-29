import os
import matplotlib

# Force standalone window for Mac
matplotlib.use('module://backend_interagg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# ═══════════════════════════════════════════════════════════════
#  CONFIGURATION
# ═══════════════════════════════════════════════════════════════
CSV_PATH = (
    '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio'
    '/2026_04_09_10_38_48--BioRep2_Double_Cells_04_09_2026_Final/BioRep2_RAW_ALL_DATA.csv')
    #'/2026_04_09_10_38_48--BioRep2_Double_Cells_04_09_2026_Final/BioRep2_HighPrecision_ALL_DATA.csv')
    #'/BioRep1_Double_03_24_2026/BioRep1_RAW_ALL_DATA.csv')
    #'/BioRep1_Double_03_24_2026/BioRep1_HighPrecision_ALL_DATA.csv')


COND_ORDER = ['TC', 'Random', 'Aligned']

COLOR_MAP = {
    'Bone Clone': '#1A466A',
    'Parental':   '#2E8B57'
}


# ═══════════════════════════════════════════════════════════════
#  DATA LOADING
# ═══════════════════════════════════════════════════════════════
def load_real_data(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing: {path}")
    df = pd.read_csv(path)
    df = df.dropna(subset=['AspectRatio'])
    if 'Status' in df.columns:
        df = df[df['Status'] == 'Viable'].copy()
        print("Status column found — using only 'Viable' cells.")
    else:
        print("No Status column found — using all cells in CSV.")
    df = df[df['Condition'].isin(COND_ORDER)].copy()
    df['Condition'] = pd.Categorical(df['Condition'],
                                     categories=COND_ORDER, ordered=True)
    return df


# ═══════════════════════════════════════════════════════════════
#  HELPERS
# ═══════════════════════════════════════════════════════════════
def stars(p):
    if p < 0.0001: return '****'
    if p < 0.001:  return '***'
    if p < 0.01:   return '**'
    if p < 0.05:   return '*'
    return 'ns'


def sig_bracket(ax, x1, x2, y, p, bracket_color, arm, gap):
    ax.plot([x1, x1, x2, x2],
            [y, y + arm, y + arm, y],
            lw=1.2, color=bracket_color, clip_on=False)
    ax.text((x1 + x2) / 2, y + arm + gap, stars(p),
            ha='center', va='bottom', fontsize=11,
            fontweight='bold', color='black')


def whisker_top(series):
    q1, q3 = np.percentile(series, [25, 75])
    return min(q3 + 1.5 * (q3 - q1), series.max())


def add_brackets(ax, plot_data, color, anchor_top, ylim_top):
    """Stack three pairwise brackets; arm/gap scale with the shared y-range."""
    pairs       = [(0, 1), (0, 2), (1, 2)]
    data_range  = ylim_top                        # use shared top for consistent spacing
    bracket_gap = data_range * 0.10
    arm         = bracket_gap * 0.30
    gap         = bracket_gap * 0.08
    base        = anchor_top + bracket_gap * 0.35

    for rank, (i, j) in enumerate(pairs):
        s1, s2 = plot_data[i], plot_data[j]
        p_val  = (stats.mannwhitneyu(s1, s2, alternative='two-sided')[1]
                  if len(s1) > 1 and len(s2) > 1 else 1.0)
        sig_bracket(ax, i, j, base + rank * bracket_gap,
                    p_val, color, arm=arm, gap=gap)


def draw_key(ax):
    ax.axis('off')
    ax.text(0.5, 0.80, "Significance:",
            transform=ax.transAxes, ha='center',
            fontsize=9, fontweight='bold', va='top')
    ax.text(0.5, 0.42,
            "**** p≤5e-05  |  *** p≤5e-04  |  ** p≤5e-03  |  * p≤5e-02  |  ns p>0.05",
            transform=ax.transAxes, ha='center', fontsize=8, va='top')
    ax.text(0.5, 0.08, "Error bars = Mean ± SEM",
            transform=ax.transAxes, ha='center',
            fontsize=8, color='#444444', va='top')


# ═══════════════════════════════════════════════════════════════
#  GLOBAL Y-LIMIT CALCULATION
#  Call once before plotting so both figures share identical axes.
# ═══════════════════════════════════════════════════════════════
def compute_shared_ylims(df, cell_lines):
    """
    Returns (box_ylim, bar_ylim) — the shared upper y-limit for
    box plots and bar graphs respectively, computed across all cell lines.
    """
    box_anchors = []
    bar_anchors = []

    for cl in cell_lines:
        df_cl = df[df['CellLine'] == cl]
        for cond in COND_ORDER:
            series = df_cl[df_cl['Condition'] == cond]['AspectRatio'].dropna()
            if len(series) > 1:
                box_anchors.append(whisker_top(series))
                bar_anchors.append(series.mean() + stats.sem(series))

    # Bracket headroom: 3 brackets × 10 % of range each, plus a small top pad
    def with_brackets(anchor_max):
        gap      = anchor_max * 0.10
        base     = anchor_max + gap * 0.35
        top      = base + 3 * gap + gap * 0.60   # 3 pairs
        return top * 1.05                         # 5 % breathing room

    box_ylim = with_brackets(max(box_anchors))
    bar_ylim = with_brackets(max(bar_anchors))
    return box_ylim, bar_ylim


# ═══════════════════════════════════════════════════════════════
#  COMBINED FIGURE
# ═══════════════════════════════════════════════════════════════
def plot_combined(df, cell_line, box_ylim, bar_ylim):
    df_cl = df[df['CellLine'] == cell_line]
    if df_cl.empty:
        print(f"Skipping {cell_line}: no data.")
        return

    color     = COLOR_MAP.get(cell_line, '#808080')
    plot_data = [df_cl[df_cl['Condition'] == c]['AspectRatio'].dropna()
                 for c in COND_ORDER]
    means     = np.array([s.mean() if len(s) else 0 for s in plot_data])
    sems      = np.array([stats.sem(s) if len(s) > 1 else 0 for s in plot_data])
    xs        = np.arange(len(COND_ORDER))

    # ── Layout ──────────────────────────────────────────────────
    fig = plt.figure(figsize=(13, 8.5))
    fig.canvas.manager.set_window_title(f'{cell_line}')

    # Taller chart rows, slim key row, generous top margin for suptitle
    gs = fig.add_gridspec(2, 2,
                          height_ratios=[5, 0.45],
                          width_ratios=[1, 1],
                          hspace=0.10,    # tight vertical gap between chart and key
                          wspace=0.30)

    ax_box = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_key = fig.add_subplot(gs[1, :])

    fig.suptitle(f"Aspect Ratio (BioRep 2): {cell_line}",
                 fontweight='bold', fontsize=17,
                 y=0.96)                  # sits clearly above subplot titles

    # ── BOX PLOT ────────────────────────────────────────────────
    bp = ax_box.boxplot(plot_data, positions=xs,
                        patch_artist=True, widths=0.45, zorder=3)

    for box in bp['boxes']:
        box.set(facecolor=color, edgecolor='black',
                linewidth=1.2, alpha=0.85)
    for med in bp['medians']:
        med.set(color='black', linewidth=1.5)
    for el in bp['whiskers'] + bp['caps']:
        el.set(color='black', linewidth=1.2)
    for flier in bp['fliers']:
        flier.set(marker='o', color=color, alpha=0.3,
                  markersize=4, markeredgecolor='none')

    for i, series in enumerate(plot_data):
        if len(series) < 2:
            continue
        ax_box.errorbar(i, means[i], yerr=sems[i],
                        fmt='D', color='white',
                        markeredgecolor='black', markeredgewidth=1.2,
                        markersize=7, capsize=5, capthick=1.4,
                        elinewidth=1.4, ecolor='black', zorder=5,
                        label='Mean ± SEM' if i == 0 else '_nolegend_')

    wt = max(whisker_top(s) for s in plot_data if len(s) > 0)
    add_brackets(ax_box, plot_data, color,
                 anchor_top=wt, ylim_top=box_ylim)

    ax_box.set_ylim(0, box_ylim)
    ax_box.set_title("Box Plot", fontsize=13, fontweight='bold', pad=8)
    ax_box.set_xticks(xs)
    ax_box.set_xticklabels(COND_ORDER, fontweight='bold', fontsize=11)
    ax_box.set_ylabel('Aspect Ratio', fontweight='bold', fontsize=11)
    ax_box.spines['top'].set_visible(False)
    ax_box.spines['right'].set_visible(False)
    ax_box.grid(axis='y', linestyle='--', alpha=0.4, zorder=0)
    ax_box.legend(loc='upper left', fontsize=8, framealpha=0.7)

    # ── BAR GRAPH ───────────────────────────────────────────────
    ax_bar.bar(xs, means, yerr=sems,
               color=color, edgecolor='black', linewidth=1.2,
               alpha=0.85, width=0.5, capsize=6,
               error_kw=dict(elinewidth=1.4, ecolor='black', capthick=1.4),
               zorder=3)

    bar_anchor = (means + sems).max()
    add_brackets(ax_bar, plot_data, color,
                 anchor_top=bar_anchor, ylim_top=bar_ylim)

    ax_bar.set_ylim(0, bar_ylim)
    ax_bar.set_title("Bar Graph (Mean ± SEM)", fontsize=13,
                     fontweight='bold', pad=8)
    ax_bar.set_xticks(xs)
    ax_bar.set_xticklabels(COND_ORDER, fontweight='bold', fontsize=11)
    ax_bar.set_ylabel('Aspect Ratio', fontweight='bold', fontsize=11)
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_bar.grid(axis='y', linestyle='--', alpha=0.4, zorder=0)

    # ── SHARED KEY ──────────────────────────────────────────────
    draw_key(ax_key)

    fig.subplots_adjust(left=0.08, right=0.96,
                        bottom=0.10, top=0.90)


# ═══════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    try:
        data = load_real_data(CSV_PATH)

        print(f"Viable cells loaded: {len(data)}")
        print(data.groupby(['CellLine', 'Condition']).size(), "\n")

        cell_lines = ['Parental', 'Bone Clone']

        # Compute shared limits once across both cell lines
        box_ylim, bar_ylim = compute_shared_ylims(data, cell_lines)

        for cl in cell_lines:
            plot_combined(data, cl, box_ylim, bar_ylim)

        plt.show(block=True)

    except Exception as e:
        print(f"Error: {e}")