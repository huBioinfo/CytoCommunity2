# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
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

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

SCORES_LONG_CSV = "data/Communication/config/EnrichScoreMatrix_long.csv"
SPEARMAN_CSV    = "data/Communication/config/Spearman_TopCCA_Pairs.csv"
OUT_DIR         = "plot/Communication/SpearmanBetweenCNs_plot/Sperman_plots"
os.makedirs(OUT_DIR, exist_ok=True)


def plot_cross_cn_pair(scores_long, cnA, ctA, cnB, ctB, rho_sp, p_sp, out_pdf):
    subA = scores_long.loc[
        (scores_long["CN"] == cnA) & (scores_long["CellType"] == ctA),
        ["Sample", "Score"]
    ].rename(columns={"Score": "Score_A"})

    subB = scores_long.loc[
        (scores_long["CN"] == cnB) & (scores_long["CellType"] == ctB),
        ["Sample", "Score"]
    ].rename(columns={"Score": "Score_B"})

    merged = pd.merge(subA, subB, on="Sample", how="inner")

    x = merged["Score_A"].to_numpy()
    y = merged["Score_B"].to_numpy()

    x_lo, x_hi = float(np.min(x)), float(np.max(x))

    fig, ax = plt.subplots(figsize=(6, 6))

    try:
        m, b = np.polyfit(x, y, 1)
        xs = np.linspace(x_lo, x_hi, 200)
        ys = m * xs + b
        ax.plot(xs, ys, linewidth=1.2, color="black", zorder=1)
    except Exception:
        pass

    ax.scatter(x, y, s=30, edgecolor="black", alpha=0.9, zorder=2)

    ax.set_xlabel(f"{ctA} enrichment score in CN{cnA}", fontsize=10)
    ax.set_ylabel(f"{ctB} enrichment score in CN{cnB}", fontsize=10)

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
        f"ρ = {rho_sp:.2f}\np = {p_sp:.2g}",
        transform=ax.transAxes,
        ha="left", va="top", fontsize=10
    )

    fig.tight_layout()
    fig.savefig(out_pdf, dpi=1200)
    plt.close(fig)

    print(f"[OK] save: {out_pdf}")
    return True


scores_long = pd.read_csv(SCORES_LONG_CSV)
scores_long = scores_long.dropna(subset=["Sample", "CN", "CellType"])
scores_long["CN"] = scores_long["CN"].astype(str).str.strip()
scores_long["CellType"] = scores_long["CellType"].astype(str).str.strip()

spearman_df = pd.read_csv(SPEARMAN_CSV)
spearman_df["CN_A"] = spearman_df["CN_A"].astype(str).str.strip()
spearman_df["CN_B"] = spearman_df["CN_B"].astype(str).str.strip()


for _, r in spearman_df.iterrows():
    cnA = r["CN_A"]
    cnB = r["CN_B"]
    ctA = str(r["Top_CellType_A"]).strip()
    ctB = str(r["Top_CellType_B"]).strip()

    out_png = os.path.join(OUT_DIR, f"{ctA} vs {ctB}.pdf")

    ok = plot_cross_cn_pair(
        scores_long,
        cnA, ctA,
        cnB, ctB,
        r.get("Spearman_rho", np.nan),
        r.get("Spearman_p", np.nan),
        out_png
    )
