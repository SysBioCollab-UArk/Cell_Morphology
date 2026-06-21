import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# =========================================================
# 1. Data Setup & Global Definitions
# =========================================================
file_path = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026/BioRep1_CellData.xlsx'
df = pd.read_excel(file_path)

df['LogAspectRatio'] = np.log(df['AspectRatio'])
df['CellLine'] = df['Group'].apply(lambda x: x.split('-')[0])
df['Condition'] = df['Group'].apply(lambda x: x.split('-')[1])

env_order = ['TC', 'Random', 'Aligned']
line_order = ['Parental', 'Bone Clone']

# =========================================================
# 2. Automated Statistical Engine
# =========================================================
def get_env_stats(cell_line, env1, env2):
    a = df[(df['CellLine'] == cell_line) & (df['Condition'] == env1)]['LogAspectRatio'].dropna()
    b = df[(df['CellLine'] == cell_line) & (df['Condition'] == env2)]['LogAspectRatio'].dropna()

    label = f"{cell_line}: {env1} vs {env2}"
    if len(a) < 2 or len(b) < 2: return [label] + ["NaN"] * 4

    _, p_welch = stats.ttest_ind(a, b, equal_var=False)
    nx, ny = len(a), len(b)
    d_val = (np.mean(a) - np.mean(b)) / np.sqrt(((nx - 1) * np.var(a) + (ny - 1) * np.var(b)) / (nx + ny - 2))
    _, p_lev = stats.levene(a, b)
    _, _, p_ad = stats.anderson_ksamp([a, b])

    return [label, f"{p_welch:.2e}", f"{d_val:.2f}", f"{p_lev:.2e}", f"{p_ad:.3f}"]

# Automatically generate all 6 environmental comparisons
stats_rows = []
comparisons = [("TC", "Random"), ("Random", "Aligned"), ("TC", "Aligned")]

for line in line_order:
    for env1, env2 in comparisons:
        stats_rows.append(get_env_stats(line, env1, env2))

# =========================================================
# 3. Two-Way ANOVA & Tukey HSD
# =========================================================
# Proves if the environment affects Bone Clone differently than Parental
model = ols('LogAspectRatio ~ C(Condition) * C(CellLine)', data=df).fit()
anova_results = sm.stats.anova_lm(model, typ=2)

# Tukey HSD for all groups
df['Combined'] = df['CellLine'] + "_" + df['Condition']
tukey = pairwise_tukeyhsd(endog=df['LogAspectRatio'], groups=df['Combined'], alpha=0.05)

# =========================================================
# 4. Sensitivity Gain Calculation
# =========================================================
def calculate_d(line, e1, e2):
    g1 = df[(df['CellLine'] == line) & (df['Condition'] == e1)]['LogAspectRatio']
    g2 = df[(df['CellLine'] == line) & (df['Condition'] == e2)]['LogAspectRatio']
    return abs((np.mean(g1) - np.mean(g2)) / np.sqrt((np.var(g1) + np.var(g2)) / 2))

d_bc = calculate_d('Bone Clone', 'TC', 'Aligned')
d_par = calculate_d('Parental', 'TC', 'Aligned')
sensitivity_gain = d_bc - d_par

# ================= ========================================
# 5. Visualization Suite
# =========================================================
fig = plt.figure(figsize=(20, 14))
gs = fig.add_gridspec(3, 2)

# A. Boxenplot
ax1 = fig.add_subplot(gs[0, 0])
sns.boxenplot(data=df, x='Condition', y='LogAspectRatio', hue='CellLine',
              order=env_order, palette='viridis', ax=ax1)
ax1.set_title('A. Environmental Response Spread', fontweight='bold')

# B. Sensitivity Gain Chart
ax2 = fig.add_subplot(gs[0, 1])
sns.barplot(x=['Parental', 'Bone Clone'], y=[d_par, d_bc], palette='viridis', ax=ax2)
ax2.set_title(f'B. Topographic Sensitivity (Gain: {sensitivity_gain:.2f})', fontweight='bold')
ax2.set_ylabel("Cohen's d (TC to Aligned)")

# C. Stats Table (Now with 6 rows)
ax3 = fig.add_subplot(gs[1, :])
ax3.axis('off')
headers = ["Comparison", "Welch (p)", "Cohen's d", "Levene (p)", "Anderson-Darling (p)"]
table = ax3.table(cellText=stats_rows, colLabels=headers, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.8)

# D. ANOVA Results Output
ax4 = fig.add_subplot(gs[2, :])
ax4.axis('off')
anova_text = f"Two-Way ANOVA Interaction P-Value: {anova_results.loc['C(Condition):C(CellLine)', 'PR(>F)']:.2e}"
ax4.text(0.5, 0.5, anova_text, size=14, ha='center', weight='bold', bbox=dict(facecolor='none', edgecolor='black'))

plt.suptitle('Mechanobiological Response: Parental vs. Bone Clone', fontsize=18, fontweight='bold')
plt.tight_layout()
plt.show()


