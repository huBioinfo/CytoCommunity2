# -*- coding: utf-8 -*-
import os
import glob
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# ===== 配置 =====
SCORE_DIR = "data/EnrichScoreMatrix"
NUM_TCN = 10

# ===== 输出目录 =====
OUT_DIR = os.path.join("data/Communication", "config")
os.makedirs(OUT_DIR, exist_ok=True)

# ===== 1) 读取 CellType 列表 =====
ct_path = os.path.join(SCORE_DIR, "CellTypeVec_List.csv")
if not os.path.exists(ct_path):
    raise FileNotFoundError(f"找不到 {ct_path}")

celltypes = pd.read_csv(ct_path).iloc[:, 0].astype(str).tolist()

# ===== 2) 整理所有样本的 EnrichScore 矩阵为长表 =====
files = sorted(glob.glob(os.path.join(SCORE_DIR, "*_EnrichScoreMatrix.csv")))
if not files:
    raise FileNotFoundError(f"未在 {SCORE_DIR} 找到 *_EnrichScoreMatrix.csv")

rows = []
for fp in files:
    sample = os.path.basename(fp).replace("_EnrichScoreMatrix.csv", "")
    mat = pd.read_csv(fp, header=None).to_numpy(dtype=float)

    if mat.shape[0] != len(celltypes):
        raise ValueError(f"{fp}: 行数 {mat.shape[0]} != CellType 数 {len(celltypes)}")
    if mat.shape[1] != NUM_TCN:
        raise ValueError(f"{fp}: 列数 {mat.shape[1]} != NUM_TCN {NUM_TCN}")

    for i, ct in enumerate(celltypes):
        for k in range(1, NUM_TCN + 1):
            rows.append({
                "Sample": sample,
                "CN": str(k),
                "CellType": ct,
                "Score": float(mat[i, k - 1])
            })

scores_long = pd.DataFrame(rows)
scores_long.to_csv(os.path.join(OUT_DIR, "EnrichScoreMatrix_long.csv"), index=False)

# ===== 3) 计算每个 CN 内所有细胞类型两两 Spearman =====
all_pairs = []

for k in range(1, NUM_TCN + 1):
    sub = scores_long[scores_long["CN"] == str(k)]
    mat_df = sub.pivot(index="Sample", columns="CellType", values="Score")

    valid_cols = [
        ct for ct in mat_df.columns
        if (mat_df[ct].notna().sum() >= 3) and (mat_df[ct].dropna().std() > 0)
    ]
    mat_df = mat_df[valid_cols]

    if mat_df.shape[1] < 2:
        continue

    corr, pval = spearmanr(mat_df.values, axis=0, nan_policy="omit")

    if np.isscalar(corr):
        corr = np.array([[1.0, corr], [corr, 1.0]])
        pval = np.array([[0.0, pval], [pval, 0.0]])

    corr_df = pd.DataFrame(corr, index=mat_df.columns, columns=mat_df.columns)
    pval_df = pd.DataFrame(pval, index=mat_df.columns, columns=mat_df.columns)

    ct_list = list(mat_df.columns)
    for i in range(len(ct_list)):
        for j in range(i + 1, len(ct_list)):   # 直接不含对角线
            ct1, ct2 = ct_list[i], ct_list[j]
            n = len(mat_df[[ct1, ct2]].dropna())
            all_pairs.append({
                "CN": str(k),
                "CellType1": ct1,
                "CellType2": ct2,
                "rho": corr_df.loc[ct1, ct2],
                "p": pval_df.loc[ct1, ct2],
                "n": n
            })

out_df = pd.DataFrame(all_pairs)
out_df.to_csv(os.path.join(OUT_DIR, "SpearmanWithinCNs.csv"), index=False)