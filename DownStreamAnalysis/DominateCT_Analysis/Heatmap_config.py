import os
import glob
import pandas as pd
import numpy as np

base_path    = "data/EnrichScoreMatrix/"
celltypes_fp = os.path.join(base_path, 'CellTypeVec_List.csv')

out_path="data/DominateCT"
os.makedirs(out_path,exist_ok=True)

freq_fp      = os.path.join(out_path, 'CN_CT_Frequency.csv')
dominant_fp  = os.path.join(out_path, 'DominantCT.csv')

SCORE_THRESH = -np.log10(0.05) 
FREQ_THRESH  = 0.5              

celltypes_df = pd.read_csv(celltypes_fp, header=0)
celltypes    = celltypes_df.iloc[:, 0].tolist()
n_types      = len(celltypes)

files = sorted(glob.glob(os.path.join(base_path, '*_EnrichScoreMatrix.csv')))
if not files:
    raise RuntimeError("未找到任何 *_EnrichScoreMatrix.csv 文件")

mats = []
n_tcns = None
for fp in files:
    mat = pd.read_csv(fp, header=None).values
    if mat.shape[0] != n_types:
        print(f"[WARN] 跳过 {os.path.basename(fp)}：行数({mat.shape[0]}) != n_types({n_types})")
        continue
    if n_tcns is None:
        n_tcns = mat.shape[1]
    elif mat.shape[1] != n_tcns:
        print(f"[WARN] 跳过 {os.path.basename(fp)}：TCN 列数不一致({mat.shape[1]} vs {n_tcns})")
        continue
    mats.append(mat)

if not mats:
    raise RuntimeError("没有合法的 EnrichScoreMatrix 文件可用于统计")

stack = np.stack(mats, axis=0)  # [n_samples, n_types, n_tcns]
n_samples = stack.shape[0]
tcns = [f"TCN{i+1}" for i in range(n_tcns)]
print(f"[INFO] 有效样本数={n_samples}, 细胞类型数={n_types}, TCN 数={n_tcns}")

dominant_bool = (stack > SCORE_THRESH).astype(np.int32)   # [S, T, C]

dominant_counts = dominant_bool.sum(axis=0)               # [T, C]
dominant_freq   = dominant_counts / float(n_samples)      # [T, C]

freq_df = pd.DataFrame(dominant_freq, index=celltypes, columns=tcns)
freq_df.to_csv(freq_fp, float_format='%.4f')

dominant_dict = {}
for j, tcn in enumerate(tcns):
    mask = dominant_freq[:, j] > FREQ_THRESH   
    selected = [celltypes[i] for i, ok in enumerate(mask) if ok]
    dominant_dict[tcn] = ", ".join(selected) if selected else ""

dominant_df = pd.DataFrame(list(dominant_dict.items()), columns=["CN", "DominantCellTypes"])
dominant_df.to_csv(dominant_fp, index=False)
print(f"[DONE] 基于样本内主导阈值与跨样本频率的 CN 注释已保存到 {dominant_fp}")
