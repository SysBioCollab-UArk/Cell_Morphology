import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import warnings

# Suppress technical warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

# =========================================================
# 1. Data Loading & Engineering
# =========================================================
file_path = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026/BioRep1_CellData.xlsx'

# Load the dataframe first so it is defined for the rest of the script
df = pd.read_excel(file_path)

# Perform engineering
df['LogAspectRatio'] = np.log(df['AspectRatio'])
df['CellLine'] = df['Group'].apply(lambda x: x.split('-')[0])
df['Condition'] = df['Group'].apply(lambda x: x.split('-')[1])

condition_order = ['TC', 'Random', 'Aligned']
cell_line_order = ['Parental', 'Bone Clone']
palette = 'Set2'

# =========================================================
# 2. Statistical Functions (Defined BEFORE use)
# =========================================================
def cohens_d(x, y):
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    pooled_var = ((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / dof
    return (np.mean(x) - np.mean(y)) / np.sqrt(pooled_var)

def variance_ratio(x, y):
    return np.var(x, ddof=1) / np.var(y, ddof=1)

def get_comprehensive_stats(label, c1, l1, c2, l2):
    a = df[(df['Condition'] == c1) & (df['CellLine'] == l1)]['LogAspectRatio'].dropna()
    b = df[(df['Condition'] == c2) & (df['CellLine'] == l2)]['LogAspectRatio'].dropna()

    if len(a) < 2 or len(b) < 2:
        return [label] + ["NaN"] * 8

    # Location Tests
    _, p_mwu = stats.mannwhitneyu(a, b)
    _, p_welch = stats.ttest_ind(a, b, equal_var=False)
    d_val = cohens_d(a, b)

    # Spread Tests
    var_rat = variance_ratio(a, b)
    _, p_lev = stats.levene(a, b, center='mean')
    _, p_bf = stats.levene(a, b, center='median')

    # Shape Tests
    _, p_ks = stats.ks_2samp(a, b)
    try:
        _, _, p_ad = stats.anderson_ksamp([a, b])
    except:
        p_ad = np.nan

    return [label, f"{p_mwu:.2e}", f"{p_welch:.2e}", f"{d_val:.2f}",
            f"{var_rat:.2f}", f"{p_lev:.2e}", f"{p_bf:.2e}", f"{p_ks:.2e}", f"{p_ad:.3f}"]

# =========================================================
# 3. Execution & Visualization
# =========================================================

# Generate Statistics Rows
rows = [
    get_comprehensive_stats("TC: Par vs BC", "TC", "Parental", "TC", "Bone Clone"),
    get_comprehensive_stats("Random: Par vs BC", "Random", "Parental", "Random", "Bone Clone"),
    get_comprehensive_stats("Aligned: Par vs BC", "Aligned", "Parental", "Aligned", "Bone Clone")
]

# --- Figure 1: Distribution Panels ---
fig_dist, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

sns.boxplot(data=df, x='Condition', y='LogAspectRatio', hue='CellLine',
            order=condition_order, hue_order=cell_line_order, palette=palette, ax=ax1)
ax1.set_title('B. Morphological Spread (Box Plot)', fontweight='bold')

sns.violinplot(data=df, x='Condition', y='LogAspectRatio', hue='CellLine',
               order=condition_order, hue_order=cell_line_order, palette=palette,
               split=True, inner='quartile', ax=ax2)
ax2.set_title('C. Density Distribution (Violin Plot)', fontweight='bold')
plt.show()

# --- Figure 2: Stats Table ---
fig_tbl, ax_tbl = plt.subplots(figsize=(16, 4))
ax_tbl.axis('off')
headers = ["Comparison", "MWU (p)", "Welch (p)", "Cohen's d", "Var. Ratio", "Levene (p)", "BF (p)", "KS (p)", "AD (p)"]
tbl = ax_tbl.table(cellText=rows, colLabels=headers, loc='center', cellLoc='center')
tbl.scale(1, 3.5)
plt.title('Comprehensive Statistical Analysis with TC Baseline', fontweight='bold', pad=20)
plt.show()