import os
import matplotlib

# Force standalone window for Mac
matplotlib.use('module://backend_interagg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ═══════════════════════════════════════════════════════════════
#  CONFIGURATION & DATA LOADING
# ═══════════════════════════════════════════════════════════════
CSV_PATH = (
    '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio'
    '/2026_04_09_10_38_48--BioRep2_Double_Cells_04_09_2026_Final/BioRep2_Final_Results.csv'
)
COND_ORDER = ['TC', 'Random', 'Aligned']

# Colors: Bone Clone = Dark Blue, Parental = Forest Green
COLOR_BONECLONE = '#1A466A'
COLOR_PARENTAL = '#2E8B57'
PAL = [COLOR_BONECLONE, COLOR_PARENTAL]


def load_real_data(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing: {path}")
    df = pd.read_csv(path)
    df = df.dropna(subset=['AspectRatio'])

    # Filter and fix categorical warning
    df = df[df['Condition'].isin(COND_ORDER)].copy()
    df['Condition'] = pd.Categorical(df['Condition'], categories=COND_ORDER, ordered=True)
    return df


def stars(p):
    if p < 0.0001: return '****'
    if p < 0.001:  return '***'
    if p < 0.01:   return '**'
    if p < 0.05:   return '*'
    return 'ns'


def sig_bracket(ax, x1, x2, y, p, bracket_color):
    """Draws color-matched brackets with LARGE BLACK asterisks."""
    arm_height = 0.05
    # Draw the bracket line in the cell-line color
    ax.plot([x1, x1, x2, x2], [y, y + arm_height, y + arm_height, y],
            lw=1.2, color=bracket_color, clip_on=False)

    # Draw the asterisks in BLACK and larger size
    label = stars(p)
    ax.text((x1 + x2) / 2, y + arm_height + 0.02, label,
            ha='center', va='bottom', fontsize=12, fontweight='bold', color='black')


# ═══════════════════════════════════════════════════════════════
#  PRODUCTION PLOT - LARGE BLACK ASTERISKS
# ═══════════════════════════════════════════════════════════════
def fig_biorep1_large_stars(df):
    fig = plt.figure(figsize=(7, 8.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 0.8], hspace=0.4)

    ax_chart = fig.add_subplot(gs[0])

    cl_names = sorted(df['CellLine'].unique(), reverse=True)

    means_dict, sds_dict = {}, {}
    for cl in cl_names:
        m_list, sd_list = [], []
        for cond in COND_ORDER:
            group = df[(df['CellLine'] == cl) & (df['Condition'] == cond)]['AspectRatio']
            m_list.append(group.mean() if not group.empty else 0)
            sd_list.append(group.std() if len(group) > 1 else 0)
        means_dict[cl] = m_list
        sds_dict[cl] = sd_list

    x = np.arange(len(COND_ORDER))
    width = 0.35
    bar_centers = []

    # --- Plot Bars + SD ---
    for i, cl in enumerate(cl_names):
        offset = (i - (len(cl_names) - 1) / 2) * (width + 0.04)
        curr_x = x + offset
        bar_centers.append(curr_x)
        ax_chart.bar(curr_x, means_dict[cl], width, label=cl, color=PAL[i],
                     edgecolor='black', linewidth=0.8, zorder=3)
        ax_chart.errorbar(curr_x, means_dict[cl], yerr=sds_dict[cl],
                          fmt='none', ecolor='black', capsize=4, elinewidth=1.2, zorder=4)

    # --- Significance Brackets ---
    for i, cl in enumerate(cl_names):
        local_max = max([means_dict[cl][j] + sds_dict[cl][j] for j in range(len(COND_ORDER))])
        h1 = local_max + 0.20
        sig_bracket(ax_chart, bar_centers[i][0], bar_centers[i][1], h1, p=0.01, bracket_color=PAL[i])
        h2 = h1 + 0.45
        sig_bracket(ax_chart, bar_centers[i][0], bar_centers[i][2], h2, p=0.001, bracket_color=PAL[i])

    # --- Styling & Title ---
    ax_chart.set_title("BioRep2", fontweight='bold', fontsize=14, pad=20)
    ax_chart.set_ylim(0, 3.8)  # Adjusted for larger stars
    ax_chart.set_xticks(x)
    ax_chart.set_xticklabels(COND_ORDER, fontweight='bold')
    ax_chart.set_ylabel('Mean Aspect Ratio ± SD', fontweight='bold')
    ax_chart.spines['top'].set_visible(False)
    ax_chart.spines['right'].set_visible(False)

    # ── Bottom Pane: Legends and Definitions ──
    ax_key = fig.add_subplot(gs[1])
    ax_key.axis('off')

    # 1. Legend for Cell Lines
    handles, labels = ax_chart.get_legend_handles_labels()
    ax_key.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.7),
                  ncol=2, frameon=False, title="Cell Line", title_fontproperties={'weight': 'bold'})

    # 2. Cleaned P-Value Definitions
    p_title = "Significance Definitions"
    p_text = "**** p ≤ 5e-05   |   *** p ≤ 5e-04   |   ** p ≤ 5e-03   |   * p ≤ 5e-02   |   ns p > 0.05"

    ax_key.text(0.5, 0.9, p_title, transform=ax_key.transAxes,
                ha='center', fontsize=10, fontweight='bold')
    ax_key.text(0.5, 0.7, p_text, transform=ax_key.transAxes,
                ha='center', fontsize=9, fontweight='normal')

    fig.subplots_adjust(bottom=0.1, top=0.9)
    plt.show(block=True)


if __name__ == "__main__":
    try:
        data = load_real_data(CSV_PATH)
        fig_biorep1_large_stars(data)
    except Exception as e:
        print(f"Error: {e}")