# -*- coding: utf-8 -*-
import os
import glob
import numpy as np
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

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# =========================
# Config
# =========================
ROOT = "data"
ES_DIR = os.path.join(ROOT, "EnrichScoreMatrix")
CELLTYPE_FILE = os.path.join(ES_DIR, "CellTypeVec_List.csv")

TTEST_BASE = os.path.join(ES_DIR, "t_test", "bias_result")
PAIRS = ["_t0-t1", "_t0-t2", "_t1-t2"]

OUT_DIR = os.path.join("plot/", "DotPlot_CN_Dynamics")
os.makedirs(OUT_DIR, exist_ok=True)

CAP_MLOG10 = 5
SIZE_MIN = 40.0
SIZE_MAX = 1200.0

CMAP_NAME = "coolwarm"



# =========================
# Helpers
# =========================
def read_celltypes(path):
    df = pd.read_csv(path)
    return df[df.columns[0]].astype(str).tolist()


def read_matrix_csv(path):
    return pd.read_csv(path, header=None).values.astype(float)


def infer_condition_from_filename(fname):
    # 0_xxx_EnrichScoreMatrix.csv
    base = os.path.basename(fname)
    suffix = "_EnrichScoreMatrix.csv"
    base = base[:-len(suffix)]
    return base.split("_")[0]   # "0" / "1" / "2"


def load_condition_means(es_dir, celltypes):
    files = glob.glob(os.path.join(es_dir, "*_EnrichScoreMatrix.csv"))
    cond_to_mats = {}

    for f in files:
        cond = infer_condition_from_filename(f)
        mat = read_matrix_csv(f)
        if mat.shape[0] != len(celltypes):
            raise ValueError("Cell type row mismatch in " + f)
        cond_to_mats.setdefault(cond, []).append(mat)

    cond_mean = {}
    for cond, mats in cond_to_mats.items():
        stack = np.stack(mats, axis=0)
        cond_mean[cond] = np.nanmean(stack, axis=0)
        print(f"[OK] condition={cond}, samples={len(mats)}")

    return cond_mean


def load_pair_pvals(pair_dir, n_cn, celltypes):
    pairname = os.path.basename(pair_dir).lstrip("_")  # t0-t1
    pmat = np.ones((len(celltypes), n_cn))

    ct_to_i = {ct: i for i, ct in enumerate(celltypes)}

    for j in range(1, n_cn + 1):
        fn = os.path.join(pair_dir, f"TCN{j}_{pairname}.csv")
        if not os.path.exists(fn):
            continue

        df = pd.read_csv(fn)
        p = np.minimum(
            df["p_value_left"].astype(float).values,
            df["p_value_right"].astype(float).values
        )
        for ct, pv in zip(df["CellType"], p):
            if ct in ct_to_i:
                pmat[ct_to_i[ct], j - 1] = pv

    return np.clip(pmat, 1e-300, 1.0)


def compute_log2fc(a_mat, b_mat):
    log2fc = np.full_like(a_mat, np.nan, dtype=float)
    mask = (a_mat > 0) & (b_mat > 0)
    log2fc[mask] = np.log2(b_mat[mask] / a_mat[mask])
    return log2fc


def compute_global_vmax(cond_mean, pair_list, percentile=98, hard_cap=4.0):
    """
    One global vmax across all pairs -> consistent colorbar.
    """
    all_abs = []

    for pair_folder in pair_list:
        pair = pair_folder.lstrip("_")  # t0-t1
        a_raw, b_raw = pair.split("-")
        a = a_raw.replace("t", "")
        b = b_raw.replace("t", "")

        if a not in cond_mean or b not in cond_mean:
            continue

        log2fc = compute_log2fc(cond_mean[a], cond_mean[b])
        finite = log2fc[np.isfinite(log2fc)]
        if finite.size > 0:
            all_abs.append(np.abs(finite))

    if len(all_abs) == 0:
        vmax = hard_cap if hard_cap is not None else 4.0
    else:
        merged = np.concatenate(all_abs, axis=0)
        vmax = np.nanpercentile(merged, percentile)
        vmax = max(vmax, 1e-6)

    if hard_cap is not None:
        vmax = min(vmax, hard_cap)

    return vmax


# =========================
# Plotting (PDF only)
# =========================
def plot_dot_only_pdf(log2fc, pvals, celltypes, out_pdf, norm):
    """
    Save dotplot ONLY (no colorbar, no size legend) as PDF.
    Return map_size (size mapping) so we can generate a shared size legend.
    """
    n_ct, n_cn = log2fc.shape

    # size mapping
    pvals = np.clip(pvals, 1e-300, 1.0)
    mlog10 = -np.log10(pvals)
    mlog10 = np.minimum(mlog10, CAP_MLOG10)
    mlog10 = np.nan_to_num(mlog10, nan=0.0, posinf=CAP_MLOG10, neginf=0.0)

    mmin, mmax = np.nanmin(mlog10), np.nanmax(mlog10)

    def map_size(m):
        if mmax == mmin:
            return (SIZE_MIN + SIZE_MAX) / 2
        return SIZE_MIN + (m - mmin) / (mmax - mmin) * (SIZE_MAX - SIZE_MIN)

    # scatter data
    xs, ys, cs, ss = [], [], [], []
    for i in range(n_ct):
        for j in range(n_cn):
            xs.append(i)
            ys.append(j)
            cs.append(log2fc[i, j])
            ss.append(map_size(mlog10[i, j]))

    xs = np.array(xs)
    ys = np.array(ys)
    cs = np.array(cs, dtype=float)
    cs[~np.isfinite(cs)] = np.nan
    ss = np.array(ss, dtype=float)

    # figure size
    fig, ax = plt.subplots(figsize=(10,8))

    ax.set_facecolor("white")
    ax.set_axisbelow(True)
    ax.grid(False)
    
    # dots
    ax.scatter(
        xs, ys,
        s=ss,
        c=cs,
        cmap=CMAP_NAME,
        norm=norm,
        edgecolors="white",
        linewidths=0.35,
        alpha=0.95
    )

    # stars: -log10(p) >= 1.3
    sig_mask = (mlog10.flatten() >= 1.3)
    if np.any(sig_mask):
        ax.scatter(
            xs[sig_mask], ys[sig_mask],
            marker="*",
            s=85,
            c="black",
            alpha=0.65,
            linewidths=0.0,
            zorder=5
        )

    # axes
    ax.set_xticks(range(n_ct))
    ax.set_xticklabels(celltypes, rotation=90,  ha="center", va="top")
    ax.set_yticks(range(n_cn))
    ax.set_yticklabels([f"CN{k}" for k in range(1, n_cn + 1)])

    ax.tick_params(axis="both", which="both", direction="out", length=4, width=0.8)
    ax.set_xlim(-0.5, n_ct - 0.5)
    ax.set_ylim(-0.5, n_cn - 0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.savefig(out_pdf, bbox_inches="tight",dpi=1200)  # PDF: vector
    plt.close(fig)

    return map_size


def save_colorbar_pdf(norm, out_pdf, cmap_name=CMAP_NAME, label="log2FC"):
    """
    Save ONLY the colorbar as a standalone PDF.
    """
    import matplotlib as mpl

    fig = plt.figure(figsize=(1.4, 3.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    cax = fig.add_axes([0.35, 0.12, 0.25, 0.80])  # [left, bottom, width, height]
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(cmap_name))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.ax.set_title(label,  pad=6)
    cbar.ax.tick_params(labelsize=10)

    fig.savefig(out_pdf, bbox_inches="tight",dpi=1200)
    plt.close(fig)


def save_size_legend_pdf(map_size, out_pdf, cap_mlog10=CAP_MLOG10, cmap_name=CMAP_NAME):
    fig = plt.figure(figsize=(2.6, 2.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    # legend axis
    lax = fig.add_axes([0.05, 0.05, 0.90, 0.90])
    lax.axis("off")

    ref_vals = np.array([0, 0.5, 1.3, 2], dtype=float)
    ref_vals = ref_vals[ref_vals <= cap_mlog10]

    warm_color = plt.get_cmap(cmap_name)(0.9)
    handles = [
        plt.scatter([], [], s=map_size(r), color=warm_color, edgecolors="white", linewidths=0.35)
        for r in ref_vals
    ]
    labels = [f"{r:g}" for r in ref_vals]

    lax.legend(
        handles, labels,
        title="-log10(pvalue)",
        frameon=False,
        loc="upper left",
        borderaxespad=0.0,
        labelspacing=1.3,
        handletextpad=1.0
    )

    fig.savefig(out_pdf, bbox_inches="tight",dpi=1200)
    plt.close(fig)


# =========================
# Main
# =========================
def main():
    celltypes = read_celltypes(CELLTYPE_FILE)
    cond_mean = load_condition_means(ES_DIR, celltypes)

    any_cond = list(cond_mean.keys())[0]
    n_cn = cond_mean[any_cond].shape[1]

    # 1) global norm (shared across pairs)
    vmax_global = compute_global_vmax(cond_mean, PAIRS, percentile=98, hard_cap=4.0)
    norm = TwoSlopeNorm(vmin=-vmax_global, vcenter=0.0, vmax=vmax_global)
    print(f"[OK] Global colorbar vmax = {vmax_global:.4g}")

    # 2) generate dotplots (PDF only); capture one map_size
    legend_map_size = None
    any_done = False

    for pair_folder in PAIRS:
        pair_dir = os.path.join(TTEST_BASE, pair_folder)
        if not os.path.isdir(pair_dir):
            continue

        pair = pair_folder.lstrip("_")  # t0-t1
        a_raw, b_raw = pair.split("-")
        a = a_raw.replace("t", "")
        b = b_raw.replace("t", "")

        if a not in cond_mean or b not in cond_mean:
            print(f"[WARN] Missing condition {a} or {b}")
            continue

        log2fc = compute_log2fc(cond_mean[a], cond_mean[b])
        pmat = load_pair_pvals(pair_dir, n_cn, celltypes)

        out_pdf = os.path.join(OUT_DIR, f"DotPlot_{pair}.pdf")
        map_size = plot_dot_only_pdf(log2fc, pmat, celltypes, out_pdf, norm)

        if legend_map_size is None:
            legend_map_size = map_size

        any_done = True
        print(f"[OK] Saved {pair}: {out_pdf}")

    # 3) generate legends ONCE (two separate pdf files)
    if any_done and legend_map_size is not None:
        out_cbar = os.path.join(OUT_DIR, "DotPlot_COLORBAR.pdf")
        out_sleg = os.path.join(OUT_DIR, "DotPlot_SIZELEGEND.pdf")
        save_colorbar_pdf(norm, out_cbar, label="log2FC")
        save_size_legend_pdf(legend_map_size, out_sleg)
        print(f"[OK] Saved legends: {out_cbar} + {out_sleg}")
    else:
        print("[WARN] No dotplot generated; legends not created.")


if __name__ == "__main__":
    main()