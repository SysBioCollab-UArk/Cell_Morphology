import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import warnings
import os
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

# =========================================================
# 1. Data Setup & Outlier Filtering
# =========================================================
file_path = '/Users/alexandravega/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026/BioRep1_CellData.xlsx'
df = pd.read_excel(file_path)

df['LogAspectRatio'] = np.log10(df['AspectRatio'])
df['CellLine'] = df['Group'].apply(lambda x: x.split('-')[0])
df['Condition'] = df['Group'].apply(lambda x: x.split('-')[1])

env_order = ['TC', 'Random', 'Aligned']

condition_labels = {
    'TC': 'Tissue Culture Plastic',
    'Random': 'Random Collagen',
    'Aligned': 'Aligned Collagen'
}

condition_palette = {
    'TC': '#003f5c',
    'Random': '#2f70b7',
    'Aligned': '#a6cee3'
}

output_folder = '/Users/alexandravega/PycharmProjects/Cell_Morphology/AspectRatio/Publication_Figures'
os.makedirs(output_folder, exist_ok=True)


def filter_outliers_and_track(data):
    clean_list = []
    outlier_list = []

    for group_name, group_df in data.groupby('Group'):
        Q1 = group_df['LogAspectRatio'].quantile(0.25)
        Q3 = group_df['LogAspectRatio'].quantile(0.75)
        IQR = Q3 - Q1

        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        is_outlier = (
            (group_df['LogAspectRatio'] < lower_bound) |
            (group_df['LogAspectRatio'] > upper_bound)
        )

        clean_list.append(group_df[~is_outlier])
        outlier_list.append(group_df[is_outlier])

    return pd.concat(clean_list), pd.concat(outlier_list)


df_clean, df_outliers = filter_outliers_and_track(df)

outlier_export_path = os.path.join(output_folder, 'Excluded_Outlier_Cells.csv')
df_outliers.to_csv(outlier_export_path, index=False)

print(f"Successfully removed {len(df_outliers)} outlier cells.")
print(f"List of excluded cells saved to: {outlier_export_path}")


# =========================================================
# 2. Statistical & Plotting Helpers
# =========================================================
def p_label(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def format_p_value(p):
    if p == 0:
        return "0"

    superscript_map = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    base, exponent = f"{p:.2e}".split("e")
    exponent = int(exponent)
    exponent_super = str(exponent).translate(superscript_map)

    return f"{base} × 10{exponent_super}"


def get_stats(data1_log, data2_log, data1_original, data2_original, label1, label2):
    n1 = len(data1_log)
    n2 = len(data2_log)

    mean1_original = np.mean(data1_original)
    sd1_log = np.std(data1_log, ddof=1)

    mean2_original = np.mean(data2_original)
    sd2_log = np.std(data2_log, ddof=1)

    _, p_welch = stats.ttest_ind(data1_log, data2_log, equal_var=False)
    _, p_lev = stats.levene(data1_log, data2_log)

    row = [
        f"{label1} vs {label2}",
        f"{n1:,}",
        f"{mean1_original:.3f}",
        f"{sd1_log:.3f}",
        f"{n2:,}",
        f"{mean2_original:.3f}",
        f"{sd2_log:.3f}",
        format_p_value(p_welch),
        format_p_value(p_lev)
    ]

    return row, p_label(p_welch)


def draw_bracket(ax, x1, x2, y, h, text):
    ax.plot(
        [x1, x1, x2, x2],
        [y, y + h, y + h, y],
        lw=1.5,
        c='black'
    )

    ax.text(
        (x1 + x2) * 0.5,
        y + h + 0.008,
        text,
        ha='center',
        va='bottom',
        color='black',
        fontweight='bold',
        fontsize=12
    )


def add_mean_red_dots(ax, data, x_order):
    for i, condition in enumerate(x_order):
        subset = data[data['Condition'] == condition]

        mean_original = subset['AspectRatio'].mean()
        mean_log10 = np.log10(mean_original)

        ax.scatter(
            i,
            mean_log10,
            color='red',
            edgecolor='black',
            s=95,
            zorder=12
        )


def add_sd_side_brackets(ax, data, x_order):
    for i, condition in enumerate(x_order):
        subset = data[data['Condition'] == condition]

        mean_original = subset['AspectRatio'].mean()
        mean_log10 = np.log10(mean_original)

        sd_log = subset['LogAspectRatio'].std(ddof=1)

        y_low = mean_log10 - sd_log
        y_high = mean_log10 + sd_log

        x = i + 0.50
        cap = 0.025

        ax.plot([x, x], [y_low, y_high], color='black', lw=1.6, zorder=11)
        ax.plot([x - cap, x + cap], [y_low, y_low], color='black', lw=1.6, zorder=11)
        ax.plot([x - cap, x + cap], [y_high, y_high], color='black', lw=1.6, zorder=11)


def add_thick_median_lines(ax, data, x_order):
    for i, condition in enumerate(x_order):
        subset = data[data['Condition'] == condition]
        median_value = subset['LogAspectRatio'].median()

        ax.hlines(
            median_value,
            i - 0.28,
            i + 0.28,
            colors='black',
            linewidth=3.5,
            zorder=9
        )


def format_axis(ax, data):
    ax.set_xlabel(
        "Condition",
        fontsize=14,
        fontweight='bold'
    )

    ax.set_ylabel(
        "log10(Aspect Ratio)",
        fontsize=14,
        fontweight='bold'
    )

    x_labels = [
        condition_labels[condition]
        for condition in env_order
    ]

    ax.set_xticklabels(x_labels, fontsize=11)

    for tick in ax.get_xticklabels():
        tick.set_fontweight('bold')

    ax.tick_params(axis='both', labelsize=12)
    sns.despine(ax=ax)


def add_legend(ax):
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker='o',
            color='red',
            label='Mean',
            markerfacecolor='red',
            markeredgecolor='black',
            markersize=8,
            linestyle='None'
        ),
        Line2D(
            [0],
            [0],
            color='black',
            lw=3.5,
            label='Median'
        ),
        Line2D(
            [0],
            [0],
            color='black',
            lw=1.6,
            label='±SD'
        )
    ]

    ax.legend(
        handles=legend_elements,
        frameon=False,
        loc='upper left',
        fontsize=10
    )


headers = [
    "Comparison",
    "n1",
    "Mean 1",
    "SD 1",
    "n2",
    "Mean 2",
    "SD 2",
    "Welch p",
    "Levene p"
]


# =========================================================
# 3. Function for Single Cell-Line Plot
# =========================================================
def plot_single_cell_line(data, cell_line_name, title, file_prefix):
    df_line = data[data['CellLine'] == cell_line_name]

    tc_log = df_line[df_line['Condition'] == 'TC']['LogAspectRatio'].dropna()
    rd_log = df_line[df_line['Condition'] == 'Random']['LogAspectRatio'].dropna()
    al_log = df_line[df_line['Condition'] == 'Aligned']['LogAspectRatio'].dropna()

    tc_original = df_line[df_line['Condition'] == 'TC']['AspectRatio'].dropna()
    rd_original = df_line[df_line['Condition'] == 'Random']['AspectRatio'].dropna()
    al_original = df_line[df_line['Condition'] == 'Aligned']['AspectRatio'].dropna()

    stats_table = []

    row, s1 = get_stats(tc_log, rd_log, tc_original, rd_original, "TC", "Random")
    stats_table.append(row)

    row, s2 = get_stats(rd_log, al_log, rd_original, al_original, "Random", "Aligned")
    stats_table.append(row)

    row, s3 = get_stats(tc_log, al_log, tc_original, al_original, "TC", "Aligned")
    stats_table.append(row)

    fig, ax = plt.subplots(figsize=(10, 6.5))

    sns.boxenplot(
        data=df_line,
        x='Condition',
        y='LogAspectRatio',
        order=env_order,
        hue='Condition',
        palette=condition_palette,
        legend=False,
        ax=ax
    )

    add_thick_median_lines(ax, df_line, env_order)
    add_mean_red_dots(ax, df_line, env_order)
    add_sd_side_brackets(ax, df_line, env_order)

    data_max = df_line['LogAspectRatio'].max()
    y1 = data_max + 0.18
    y2 = data_max + 0.30
    y3 = data_max + 0.42
    bracket_h = 0.035
    y_top = data_max + 0.55

    draw_bracket(ax, 0, 1, y1, bracket_h, s1)
    draw_bracket(ax, 1, 2, y2, bracket_h, s2)
    draw_bracket(ax, 0, 2, y3, bracket_h, s3)

    ax.set_ylim(bottom=-0.2, top=y_top)

    ax.set_title(
        title,
        fontweight='bold',
        fontsize=16
    )

    format_axis(ax, df_line)
    add_legend(ax)

    fig.tight_layout()

    fig_path = os.path.join(output_folder, f"{file_prefix}_BoxenPlot.png")
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

    fig_table, ax_table = plt.subplots(figsize=(13, 3))
    ax_table.axis('off')

    table = ax_table.table(
        cellText=stats_table,
        colLabels=headers,
        loc='center',
        cellLoc='center',
        colWidths=[
            0.20,
            0.08,
            0.10,
            0.10,
            0.08,
            0.10,
            0.10,
            0.12,
            0.12
        ]
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    for (row_idx, col_idx), cell in table.get_celld().items():
        if row_idx == 0:
            cell.set_text_props(weight='bold', ha='center')
            cell.set_facecolor('#eeeeee')
        else:
            cell.set_text_props(ha='center')

    ax_table.set_title(
        f"{title} Statistical Summary",
        fontweight='bold',
        fontsize=13
    )

    fig_table.tight_layout()

    table_path = os.path.join(output_folder, f"{file_prefix}_StatsTable.png")
    fig_table.savefig(table_path, dpi=300, bbox_inches='tight')

    table_csv_path = os.path.join(output_folder, f"{file_prefix}_StatsTable.csv")
    pd.DataFrame(stats_table, columns=headers).to_csv(table_csv_path, index=False)

    print(f"Saved plot: {fig_path}")
    print(f"Saved table image: {table_path}")
    print(f"Saved table CSV: {table_csv_path}")

    return stats_table


# =========================================================
# 4. Function for Cross-Line Comparison Plot
# =========================================================
def plot_cross_line_comparison(data, title, file_prefix):
    data = data.copy()

    cell_line_order = ['Parental', 'Bone Clone']

    data['CellLine'] = pd.Categorical(
        data['CellLine'],
        categories=cell_line_order,
        ordered=True
    )

    df_par = data[data['CellLine'] == 'Parental']
    df_bc = data[data['CellLine'] == 'Bone Clone']

    tc_par_log = df_par[df_par['Condition'] == 'TC']['LogAspectRatio'].dropna()
    rd_par_log = df_par[df_par['Condition'] == 'Random']['LogAspectRatio'].dropna()
    al_par_log = df_par[df_par['Condition'] == 'Aligned']['LogAspectRatio'].dropna()

    tc_bc_log = df_bc[df_bc['Condition'] == 'TC']['LogAspectRatio'].dropna()
    rd_bc_log = df_bc[df_bc['Condition'] == 'Random']['LogAspectRatio'].dropna()
    al_bc_log = df_bc[df_bc['Condition'] == 'Aligned']['LogAspectRatio'].dropna()

    tc_par_original = df_par[df_par['Condition'] == 'TC']['AspectRatio'].dropna()
    rd_par_original = df_par[df_par['Condition'] == 'Random']['AspectRatio'].dropna()
    al_par_original = df_par[df_par['Condition'] == 'Aligned']['AspectRatio'].dropna()

    tc_bc_original = df_bc[df_bc['Condition'] == 'TC']['AspectRatio'].dropna()
    rd_bc_original = df_bc[df_bc['Condition'] == 'Random']['AspectRatio'].dropna()
    al_bc_original = df_bc[df_bc['Condition'] == 'Aligned']['AspectRatio'].dropna()

    stats_table = []

    row, s_tc = get_stats(
        tc_par_log, tc_bc_log,
        tc_par_original, tc_bc_original,
        "Parental TC", "Bone Clone TC"
    )
    stats_table.append(row)

    row, s_rd = get_stats(
        rd_par_log, rd_bc_log,
        rd_par_original, rd_bc_original,
        "Parental Random", "Bone Clone Random"
    )
    stats_table.append(row)

    row, s_al = get_stats(
        al_par_log, al_bc_log,
        al_par_original, al_bc_original,
        "Parental Aligned", "Bone Clone Aligned"
    )
    stats_table.append(row)

    fig, ax = plt.subplots(figsize=(11, 6.5))

    sns.boxenplot(
        data=data,
        x='Condition',
        y='LogAspectRatio',
        hue='CellLine',
        order=env_order,
        hue_order=cell_line_order,
        palette='viridis',
        width=0.85,
        ax=ax
    )

    data_max = data['LogAspectRatio'].max()
    y1 = data_max + 0.12
    bracket_h = 0.035
    y_top = data_max + 0.30

    draw_bracket(ax, -0.2, 0.2, y1, bracket_h, s_tc)
    draw_bracket(ax, 0.8, 1.2, y1, bracket_h, s_rd)
    draw_bracket(ax, 1.8, 2.2, y1, bracket_h, s_al)

    ax.set_ylim(bottom=-0.2, top=y_top)

    ax.set_title(
        title,
        fontweight='bold',
        fontsize=16
    )

    ax.set_xlabel(
        "Condition",
        fontsize=14,
        fontweight='bold'
    )

    ax.set_ylabel(
        "log10(Aspect Ratio)",
        fontsize=14,
        fontweight='bold'
    )

    x_labels = [
        condition_labels[condition]
        for condition in env_order
    ]

    ax.set_xticklabels(x_labels, fontsize=11)

    for tick in ax.get_xticklabels():
        tick.set_fontweight('bold')

    ax.tick_params(axis='both', labelsize=12)

    ax.legend(
        title='Cell Line',
        frameon=False,
        loc='upper left',
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
        fontsize=11,
        title_fontsize=11
    )

    sns.despine(ax=ax)

    fig.text(
        0.5,
        -0.02,
        "n represents total cells per scaffold condition across both cell lines after outlier removal. "
        "Exact n values for each comparison are shown in the statistical summary table.",
        ha='center',
        va='top',
        fontsize=10
    )

    fig.tight_layout()

    fig_path = os.path.join(output_folder, f"{file_prefix}_BoxenPlot.png")
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')

    fig_table, ax_table = plt.subplots(figsize=(13, 3))
    ax_table.axis('off')

    table = ax_table.table(
        cellText=stats_table,
        colLabels=headers,
        loc='center',
        cellLoc='center',
        colWidths=[
            0.20,
            0.08,
            0.10,
            0.10,
            0.08,
            0.10,
            0.10,
            0.12,
            0.12
        ]
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    for (row_idx, col_idx), cell in table.get_celld().items():
        if row_idx == 0:
            cell.set_text_props(weight='bold', ha='center')
            cell.set_facecolor('#eeeeee')
        else:
            cell.set_text_props(ha='center')

    ax_table.set_title(
        f"{title} Statistical Summary",
        fontweight='bold',
        fontsize=13
    )

    fig_table.tight_layout()

    table_path = os.path.join(output_folder, f"{file_prefix}_StatsTable.png")
    fig_table.savefig(table_path, dpi=300, bbox_inches='tight')

    table_csv_path = os.path.join(output_folder, f"{file_prefix}_StatsTable.csv")
    pd.DataFrame(stats_table, columns=headers).to_csv(table_csv_path, index=False)

    print(f"Saved plot: {fig_path}")
    print(f"Saved table image: {table_path}")
    print(f"Saved table CSV: {table_csv_path}")

    return stats_table


# =========================================================
# 5. Image 1: Parental Line Only
# =========================================================
stats_par = plot_single_cell_line(
    data=df_clean,
    cell_line_name='Parental',
    title='Parental Cells',
    file_prefix='Parental'
)


# =========================================================
# 6. Image 2: Bone Clone Line Only
# =========================================================
stats_bc = plot_single_cell_line(
    data=df_clean,
    cell_line_name='Bone Clone',
    title='Bone Clone Cells',
    file_prefix='BoneClone'
)


# =========================================================
# 7. Image 3: Parental vs Bone Clone Together
# =========================================================
stats_cross = plot_cross_line_comparison(
    data=df_clean,
    title='Parental vs Bone Clone by Scaffold',
    file_prefix='Parental_vs_BoneClone'
)

plt.show()