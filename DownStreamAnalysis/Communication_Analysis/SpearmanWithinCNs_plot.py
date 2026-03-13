import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "pdf.fonttype": 42,
})

SCORES_LONG_CSV = "data/Communication/config/EnrichScoreMatrix_long.csv"
SIG_PAIRS_CSV   = "data/Communication/config/SpearmanWithinCNs.csv"
OUT_DIR         = "plot/Communication/SpearmanWithinCNs_plot"
os.makedirs(OUT_DIR, exist_ok=True)

RHO_CUTOFF = 0.4

def plot_one_pair(wide, cn, ct1, ct2, rho, p, out_pdf):
    xy = wide[[ct1, ct2]].dropna()
    x = xy[ct1].to_numpy()
    y = xy[ct2].to_numpy()

    x_lo, x_hi = float(np.min(x)), float(np.max(x))

    fig, ax = plt.subplots(figsize=(6, 6))

    # 线性拟合线，仅作视觉引导
    try:
        m, b = np.polyfit(x, y, 1)
        xs = np.linspace(x_lo, x_hi, 200)
        ys = m * xs + b
        ax.plot(xs, ys, linewidth=1.2, color="black", zorder=1)
    except Exception:
        pass

    ax.scatter(x, y, s=30, edgecolor="black", alpha=0.9, zorder=2)

    ax.set_xlabel(f"{ct1} enrichment score in CN{cn}", fontsize=10)
    ax.set_ylabel(f"{ct2} enrichment score in CN{cn}", fontsize=10)

    AX_MIN = 0.0
    AX_MAX = 20.0
    AX_PAD = 0.5
    TICK_STEP = 2.5

    ax.set_xlim(AX_MIN - AX_PAD, AX_MAX + AX_PAD)
    ax.set_ylim(AX_MIN - AX_PAD, AX_MAX + AX_PAD)

    ticks = np.arange(AX_MIN, AX_MAX + TICK_STEP, TICK_STEP)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.1f"))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1f"))
    ax.tick_params(axis="x", labelrotation=30, labelsize=10)
    ax.tick_params(axis="y", labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.text(
        0.80, 0.98,
        f"ρ = {rho:.2f}\np = {p:.2g}",
        transform=ax.transAxes,
        ha="left", va="top", fontsize=10
    )

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight", dpi=1200)
    plt.close(fig)
    print(f"[OK] save: {out_pdf}")
    return True

scores_long = pd.read_csv(SCORES_LONG_CSV)
sig_pairs = pd.read_csv(SIG_PAIRS_CSV)

sig_pairs = sig_pairs[sig_pairs["rho"].abs() > RHO_CUTOFF]

wide_by_cn = {
    str(cn): sub.pivot(index="Sample", columns="CellType", values="Score")
    for cn, sub in scores_long.groupby("CN")
}

for _, row in sig_pairs.iterrows():
    cn = str(row["CN"])
    ct1 = str(row["CellType1"])
    ct2 = str(row["CellType2"])
    rho = float(row["rho"])
    p = float(row["p"])

    cn_dir = os.path.join(OUT_DIR, f"CN{cn}")
    os.makedirs(cn_dir, exist_ok=True)
    out_path = os.path.join(cn_dir, f"{ct1} vs {ct2} CN{cn}.pdf")
    plot_one_pair(wide_by_cn[cn], cn, ct1, ct2, rho, p, out_path)
        
