import seaborn as sns
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import os

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "pdf.fonttype": 42,
})

out_path="plot/OverlapCoefficient/"
os.makedirs(out_path,exist_ok=True)

df = pd.read_csv("data/OverlapCoefficient/config.csv")
df["TCN"] = df["TCN"].astype(str).str.strip()
tcn_levels = [f"TCN{i}" for i in range(1, 11)]
df["TCN"] = pd.Categorical(df["TCN"], categories=tcn_levels, ordered=True)
df = df.dropna(subset=["TCN"])

fig, ax = plt.subplots(figsize=(8, 6))

sns.violinplot(
    data=df,
    x="TCN",
    y="OverlapCoefficient",
    order=tcn_levels,
    inner=None,
    cut=0,
    color="white",
    linewidth=1,
    ax=ax
)

for poly in ax.collections:
    poly.set_edgecolor("black")
    poly.set_facecolor("white")

means = df.groupby("TCN", observed=False)["OverlapCoefficient"].mean()

for i, tcn in enumerate(tcn_levels):
    if pd.notna(means[tcn]):
        ax.plot([i - 0.25, i + 0.25], [means[tcn]] * 2, color="red", lw=1.5)

ax.set_xticks(range(len(tcn_levels)))
ax.set_xticklabels([x.replace("TCN", "CN") for x in tcn_levels])
ax.set_xlabel("")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_facecolor("white")
fig.patch.set_facecolor("white")

plt.tight_layout()
plt.savefig(out_path+ "OC.pdf", dpi=1200)
plt.close()