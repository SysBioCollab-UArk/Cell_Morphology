import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from bokeh.transform import dodge
from scipy.stats import mannwhitneyu, skew
from itertools import combinations

# Define the main directory
base_dir = "output"

# Define the subfolders to process
folders = ["BONE_TC", "BONE_RC", "BONE_AC", "PAR_TC", "PAR_RC", "PAR_AC"]

# Create an empty list to hold all data
data = []

# Loop through each folder
for folder in folders:
    folder_path = os.path.join(base_dir, folder)

    # Loop through all CSV files in the folder
    for file in os.listdir(folder_path):
        if file.endswith(".csv"):
            file_path = os.path.join(folder_path, file)

            # Read CSV
            df = pd.read_csv(file_path)

            # Make sure the expected columns exist
            if "Class" in df.columns and "AspectRatio" in df.columns:
                # Filter for alive cells only (case-insensitive)
                alive_df = df[df["Class"].astype(str).str.lower() == "alive"].copy()

                # Add a column for the folder label
                alive_df.loc[:, "condition"] = folder

                # Keep only aspect ratio and condition columns
                data.append(alive_df[["AspectRatio", "condition"]])

# Combine all data
if data:
    all_data = pd.concat(data, ignore_index=True)
else:
    raise ValueError("No data found. Please check your folder structure or column names.")

# --- Violin Plot (Seaborn ≥ 0.14 compatible) ---
plt.figure(figsize=(10, 6))
sns.violinplot(
    data=all_data,
    x="condition",
    y="AspectRatio",
    hue="condition",       # explicitly assign hue
    dodge=False,           # overlays each hue on its own x category
    legend=False,          # hides duplicate legend
    palette="Set2",
    inner="box"
)

plt.title("Aspect Ratio Distributions by Condition")
plt.xlabel("Condition")
plt.ylabel("Aspect Ratio")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
# -------------------------
# Summary Statistics (Extended)
# -------------------------
print("\n===== Summary Statistics =====")
summary = (
    all_data.groupby("condition")["AspectRatio"]
    .agg(["mean", "median", "std", "count"])
    .reset_index()
)

# Add skewness for each condition
skews = all_data.groupby("condition")["AspectRatio"].apply(skew).reset_index(name="skew")
summary = pd.merge(summary, skews, on="condition")

print(summary.to_string(index=False))

# Save summary statistics
summary.to_csv("summary_statistics.csv", index=False)
print("\nSummary statistics saved as: summary_statistics.csv")

# -------------------------
#  Mann–Whitney U Tests
# -------------------------
# Define Parental and Bone groups
par_groups = ["PAR_TC", "PAR_AC", "PAR_RC"]
bone_groups = ["BONE_TC", "BONE_AC", "BONE_RC"]

# Function to perform pairwise Mann–Whitney U tests
def pairwise_mannwhitney(groups_list, data, label):
    print(f"\n===== Pairwise Mann–Whitney U Tests ({label}) =====")
    results = []
    for cond1, cond2 in combinations(groups_list, 2):
        data1 = data.loc[data["condition"] == cond1, "AspectRatio"]
        data2 = data.loc[data["condition"] == cond2, "AspectRatio"]
        stat, p = mannwhitneyu(data1, data2, alternative="two-sided")
        results.append((cond1, cond2, stat, p))
        print(f"{cond1} vs {cond2}: U={stat:.2f}, p={p:.4e}")
    return pd.DataFrame(results, columns=["Group1", "Group2", "U-statistic", "p-value"])

# Run pairwise tests
par_results = pairwise_mannwhitney(par_groups, all_data, "Parental Groups")
bone_results = pairwise_mannwhitney(bone_groups, all_data, "Bone Clone Groups")

# Save statistical results
par_results.to_csv("mannwhitney_par_results.csv", index=False)
bone_results.to_csv("mannwhitney_bone_results.csv", index=False)

print("\nStatistical results saved as:")
print(" - mannwhitney_par_results.csv")
print(" - mannwhitney_bone_results.csv")

# =====================================================
#  NEW SECTION: BAR PLOTS FOR SUMMARY STATISTICS
# =====================================================

sns.set(style="whitegrid")

# Mean & Median Plot with Errors Bars (±SD on Mean)
plt.figure(figsize=(10, 6))

# Mean barplot with SD error bars
ax = sns.barplot(
    data=summary,
    x="condition",
    y="mean",
    color="skyblue",
    label="Mean(±SD)",
    errorbar=("sd"), #It plots the  SD error bars
    capsize=0.2
)

# Overlay median values as orange points
plt.scatter(
    x=range(len(summary)),
    y=summary["median"],
    color="orange",
    s=100,
    zorder=3,
    label="Median"
)

# Add text labels for mean ± SD (Optional)
for i, row in summary.iterrows():
    plt.text(i, row["mean"] + row["std"] + 0.05,f"{row['mean']:.2f}]± {row['std']:.2f}",ha='center', va='bottom', fontsize=8, color="black")

# Aesthetics
plt.title("Mean (±SD)and Median Aspect Ratio per Condition", fontsize=14, weight="bold")
plt.xlabel("Condition", fontsize=12)
plt.ylabel("Aspect Ratio", fontsize=12)
plt.xticks(rotation=45)
plt.legend(title="Statistic")
plt.tight_layout()
plt.savefig("aspect_ratio_mean_median_barplot.png", dpi=300)
plt.show()

# 2️Standard Deviation Bar Plot
plt.figure(figsize=(10, 6))
sns.barplot(
    data=summary,
    x="condition",
    y="std",
    hue="condition",
    dodge=False,
    legend=False,
    palette="muted"
)

plt.title("Standard Deviation of Aspect Ratio per Condition", fontsize=14, weight="bold")
plt.xlabel("Condition", fontsize=12)
plt.ylabel("Standard Deviation", fontsize=12)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("aspect_ratio_std_barplot.png", dpi=300)
plt.show()

# Skewness Plot
plt.figure(figsize=(10,6))
sns.barplot(
    data=summary,
    x="condition",
    y="skew",
    hue="condition",
    dodge=False,
    legend=False,
    palette="muted"
)
plt.title("Skewness of Aspect Ratio per Condition", fontsize=14, weight="bold")
plt.xlabel("Condition", fontsize=12)
plt.ylabel("Skewness", fontsize=12)
plt.axhline(0, color="black", linewidth=1)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("aspect_ratio_skew_barplot.png", dpi=300)
plt.show()
