import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import warnings

warnings.filterwarnings("ignore")

# =========================================================
# 1. Data Setup & Outlier Filtering
# =========================================================
file_path = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026/BioRep1_CellData.xlsx'
df = pd.read_excel(file_path)

df['LogAspectRatio'] = np.log(df['AspectRatio'])
df['CellLine'] = df['Group'].apply(lambda x: x.split('-')[0])
df['Condition'] = df['Group'].apply(lambda x: x.split('-')[1])

env_order = ['TC', 'Random', 'Aligned']


def filter_outliers_and_track(data):
    """Filters outliers using the IQR method per group and tracks dropped cells."""
    clean_list = []
    outlier_list = []

    # Calculate IQR bounds independently for EVERY group
    for group_name, group_df in data.groupby('Group'):
        Q1 = group_df['LogAspectRatio'].quantile(0.25)
        Q3 = group_df['LogAspectRatio'].quantile(0.75)
        IQR = Q3 - Q1

        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        # Mask for outliers
        is_outlier = (group_df['LogAspectRatio'] < lower_bound) | (group_df['LogAspectRatio'] > upper_bound)

        clean_list.append(group_df[~is_outlier])
        outlier_list.append(group_df[is_outlier])

    df_clean = pd.concat(clean_list)
    df_outliers = pd.concat(outlier_list)
    return df_clean, df_outliers


# Apply the filter
df_clean, df_outliers = filter_outliers_and_track(df)

# Export the dropped cells so you can inspect them later
outlier_export_path = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/Excluded_Outlier_Cells.csv'
df_outliers.to_csv(outlier_export_path, index=False)
print(f"Successfully removed {len(df_outliers)} outlier cells.")
print(f"List of excluded cells saved to: {outlier_export_path}")


# =========================================================
# 2. Statistical & Annotation Helpers (Using df_clean)
# =========================================================
def get_stats(data1, data2, label1, label2):
    m1, v1 = np.mean(data1), np.var(data1, ddof=1)
    m2, v2 = np.mean(data2), np.var(data2, ddof=1)

    _, p_welch = stats.ttest_ind(data1, data2, equal_var=False)
    _, p_lev = stats.levene(data1, data2)
    _, p_ks = stats.ks_2samp(data1, data2)

    stars = "***" if p_welch < 0.001 else "**" if p_welch < 0.01 else "*" if p_welch < 0.05 else "ns"

    row = [f"{label1} vs {label2}", f"{m1:.3f}", f"{v1:.3f}", f"{m2:.3f}", f"{v2:.3f}",
           f"{p_welch:.2e}", f"{p_lev:.2e}", f"{p_ks:.3f}"]
    return row, stars


def draw_bracket(ax, x1, x2, y, h, text):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c='black')
    ax.text((x1 + x2) * 0.5, y + h + 0.05, text, ha='center', va='bottom', color='black', fontweight='bold')


headers = ["Comparison", "Mean 1", "Var 1", "Mean 2", "Var 2", "Welch (p)", "Levene (p)", "KS (p)"]

# Set plot ceiling based on the CLEAN data
y_max_global = df_clean['LogAspectRatio'].max()
bracket_h = 0.08

# =========================================================
# Image 1: Parental Line Only
# =========================================================
fig1 = plt.figure(figsize=(14, 12))
gs1 = fig1.add_gridspec(2, 1, height_ratios=[2, 1])

df_par = df_clean[df_clean['CellLine'] == 'Parental']
tc_par = df_par[df_par['Condition'] == 'TC']['LogAspectRatio'].dropna()
rd_par = df_par[df_par['Condition'] == 'Random']['LogAspectRatio'].dropna()
al_par = df_par[df_par['Condition'] == 'Aligned']['LogAspectRatio'].dropna()

stats_par = []
row, s1 = get_stats(tc_par, rd_par, "TC", "Random")
stats_par.append(row)
row, s2 = get_stats(rd_par, al_par, "Random", "Aligned")
stats_par.append(row)
row, s3 = get_stats(tc_par, al_par, "TC", "Aligned")
stats_par.append(row)

ax1 = fig1.add_subplot(gs1[0])
sns.boxenplot(data=df_par, x='Condition', y='LogAspectRatio', order=env_order,
              hue='Condition', palette='crest', legend=False, ax=ax1)

draw_bracket(ax1, 0, 1, y_max_global + 0.1, bracket_h, s1)
draw_bracket(ax1, 1, 2, y_max_global + 0.35, bracket_h, s2)
draw_bracket(ax1, 0, 2, y_max_global + 0.6, bracket_h, s3)
ax1.set_ylim(bottom=-0.2, top=y_max_global + 0.9)
ax1.set_title('Parental Line: Bulk Population Response (Outliers Removed)', fontweight='bold', fontsize=14)

ax1_tab = fig1.add_subplot(gs1[1])
ax1_tab.axis('off')
table1 = ax1_tab.table(cellText=stats_par, colLabels=headers, loc='center', cellLoc='center')
table1.scale(1, 2)
fig1.tight_layout()

# =========================================================
# Image 2: Bone Clone Line Only
# =========================================================
fig2 = plt.figure(figsize=(14, 12))
gs2 = fig2.add_gridspec(2, 1, height_ratios=[2, 1])

df_bc = df_clean[df_clean['CellLine'] == 'Bone Clone']
tc_bc = df_bc[df_bc['Condition'] == 'TC']['LogAspectRatio'].dropna()
rd_bc = df_bc[df_bc['Condition'] == 'Random']['LogAspectRatio'].dropna()
al_bc = df_bc[df_bc['Condition'] == 'Aligned']['LogAspectRatio'].dropna()

stats_bc = []
row, s1 = get_stats(tc_bc, rd_bc, "TC", "Random")
stats_bc.append(row)
row, s2 = get_stats(rd_bc, al_bc, "Random", "Aligned")
stats_bc.append(row)
row, s3 = get_stats(tc_bc, al_bc, "TC", "Aligned")
stats_bc.append(row)

ax2 = fig2.add_subplot(gs2[0])
sns.boxenplot(data=df_bc, x='Condition', y='LogAspectRatio', order=env_order,
              hue='Condition', palette='flare', legend=False, ax=ax2)

draw_bracket(ax2, 0, 1, y_max_global + 0.1, bracket_h, s1)
draw_bracket(ax2, 1, 2, y_max_global + 0.35, bracket_h, s2)
draw_bracket(ax2, 0, 2, y_max_global + 0.6, bracket_h, s3)
ax2.set_ylim(bottom=-0.2, top=y_max_global + 0.9)
ax2.set_title('Bone Clone Line: Bulk Population Response (Outliers Removed)', fontweight='bold', fontsize=14)

ax2_tab = fig2.add_subplot(gs2[1])
ax2_tab.axis('off')
table2 = ax2_tab.table(cellText=stats_bc, colLabels=headers, loc='center', cellLoc='center')
table2.scale(1, 2)
fig2.tight_layout()

# =========================================================
# Image 3: Cross-Line Comparison (Same Scaffold Types)
# =========================================================
fig3 = plt.figure(figsize=(14, 12))
gs3 = fig3.add_gridspec(2, 1, height_ratios=[2, 1])

stats_cross = []
row, s_tc = get_stats(tc_par, tc_bc, "Par|TC", "BC|TC")
stats_cross.append(row)
row, s_rd = get_stats(rd_par, rd_bc, "Par|Random", "BC|Random")
stats_cross.append(row)
row, s_al = get_stats(al_par, al_bc, "Par|Aligned", "BC|Aligned")
stats_cross.append(row)

ax3 = fig3.add_subplot(gs3[0])
sns.boxenplot(data=df_clean, x='Condition', y='LogAspectRatio', hue='CellLine',
              order=env_order, palette='viridis', ax=ax3)

draw_bracket(ax3, -0.2, 0.2, y_max_global + 0.1, bracket_h, s_tc)
draw_bracket(ax3, 0.8, 1.2, y_max_global + 0.1, bracket_h, s_rd)
draw_bracket(ax3, 1.8, 2.2, y_max_global + 0.1, bracket_h, s_al)
ax3.set_ylim(bottom=-0.2, top=y_max_global + 0.4)
ax3.set_title('Direct Comparison: Parental vs Bone Clone by Scaffold (Outliers Removed)', fontweight='bold',
              fontsize=14)
ax3.legend(loc='upper left')

ax3_tab = fig3.add_subplot(gs3[1])
ax3_tab.axis('off')
table3 = ax3_tab.table(cellText=stats_cross, colLabels=headers, loc='center', cellLoc='center')
table3.scale(1, 2)
fig3.tight_layout()

plt.show()