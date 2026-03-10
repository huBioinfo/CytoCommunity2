import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "pdf.fonttype": 42,
})


csv_path = "data/Cohernece/config.csv"
out_dir  = "plot/Cohernece"
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(csv_path)
df = df[~df["Sample"].astype(str).str.strip().str.lower().eq("average")].copy()

for c in ["CHAOS", "PAS"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna(subset=["CHAOS", "PAS"])

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["grid.alpha"] = 0.25

def boxplot(y_col, title, fname, AX_MAX=0.5):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=1200)

    sns.boxplot(
        data=df,
        y=y_col,
        width=0.45,
        showfliers=False,
        linewidth=1.4,
        ax=ax
    )

    sns.stripplot(
        data=df,
        y=y_col,
        jitter=0.22,
        size=5,
        alpha=0.75,
        ax=ax
    )

    ax.set_title(title, pad=12, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")  
    
    AX_MIN = 0.0
    ax.set_ylim(AX_MIN , AX_MAX)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    out_path = os.path.join(out_dir, fname)
    fig.savefig(out_path, dpi=1200, bbox_inches="tight")
    plt.close(fig)

    print(f"[SAVE] {out_path}")

boxplot(
    y_col="CHAOS",
    title=None,
    fname="CHAOS.pdf",
    AX_MAX=0.0175
)

boxplot(
    y_col="PAS",
    title=None,
    fname="PAS.pdf"
)
