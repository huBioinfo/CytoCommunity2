# -*- coding: utf-8 -*-
import os
import pandas as pd
from scipy.stats import spearmanr

CCA_DIR = "data/Communication/config"
OUT_DIR = CCA_DIR

# ===== 1 读取 enrichment long table =====
scores_long = pd.read_csv(os.path.join(CCA_DIR, "EnrichScoreMatrix_long.csv"))
scores_long = scores_long.dropna(subset=["Sample", "CN", "CellType"])
scores_long["CN"] = scores_long["CN"].astype(str).str.strip()

# ===== 2 读取 CCA CN 对（取前10）=====
pairs = pd.read_csv(os.path.join(CCA_DIR, "CCA_config.csv"))
pairs = pairs.sort_values("rho1", ascending=False).head(10)

rows = []

# ===== 3 对每个 CN pair 做 Spearman =====
for _, r in pairs.iterrows():

    cnA = str(int(r["CN_A"]))
    cnB = str(int(r["CN_B"]))

    coef_file = os.path.join(CCA_DIR, f"CCA_coefficients_CN{cnA}_vs_CN{cnB}.csv")
    if not os.path.exists(coef_file):
        print("⚠️ 缺少文件:", coef_file)
        continue

    coef = pd.read_csv(coef_file)

    if not {"CellType","CN_A_coef","CN_B_coef"}.issubset(coef.columns):
        print("⚠️ 文件列缺失:", coef_file)
        continue

    # 取 canonical coefficient 最大的 cell type
    topA = coef.loc[coef["CN_A_coef"].abs().idxmax(),"CellType"]
    topB = coef.loc[coef["CN_B_coef"].abs().idxmax(),"CellType"]

    subA = scores_long.loc[
        (scores_long["CN"]==cnA) & (scores_long["CellType"]==topA),
        ["Sample","Score"]
    ].rename(columns={"Score":"Score_A"})

    subB = scores_long.loc[
        (scores_long["CN"]==cnB) & (scores_long["CellType"]==topB),
        ["Sample","Score"]
    ].rename(columns={"Score":"Score_B"})

    merged = pd.merge(subA, subB, on="Sample", how="inner")

    print(f"🔹 CN{cnA}({topA}) vs CN{cnB}({topB}) overlap={len(merged)}")

    if len(merged) < 3:
        print("⚠️ 样本重叠不足")
        continue

    if merged["Score_A"].nunique()<2 or merged["Score_B"].nunique()<2:
        print("⚠️ 其中一侧几乎常数")
        continue

    rho, p = spearmanr(merged["Score_A"], merged["Score_B"], nan_policy="omit")

    rows.append({
        "CN_A":cnA,
        "CN_B":cnB,
        "rho1_CCA":r["rho1"],
        "pval_CCA":r["pval"],
        "Top_CellType_A":topA,
        "Top_CellType_B":topB,
        "Spearman_rho":rho,
        "Spearman_p":p,
        "N_overlap":len(merged)
    })

# ===== 4 输出 =====
res_df = pd.DataFrame(rows).sort_values("rho1_CCA", ascending=False)

out_csv = os.path.join(OUT_DIR,"Spearman_TopCCA_Pairs.csv")
res_df.to_csv(out_csv, index=False)

print("\n✅ 结果已保存:", out_csv)