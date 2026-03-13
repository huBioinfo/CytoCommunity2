import numpy as np
from matplotlib.colors import BoundaryNorm
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


out_path="plot/DominateCT/"
os.makedirs(out_path,exist_ok=True)

mat_fp = "data/DominateCT/CN_CT_Frequency.csv"
pmat = pd.read_csv(mat_fp, header=0, index_col=0)

pmat.columns = [c.replace("TCN", "CN", 1) if str(c).startswith("TCN") else c for c in pmat.columns]

# 转成数值矩阵
pmat = pmat.apply(pd.to_numeric, errors="coerce")

brks = np.concatenate([
    np.linspace(0.00, 0.40, 100),
    np.linspace(0.401, 0.599, 300),
    np.linspace(0.60, 1.00, 100)
])

cmap = plt.get_cmap("RdYlBu_r", len(brks) - 1)
norm = BoundaryNorm(brks, ncolors=cmap.N, clip=True)

fig, ax = plt.subplots(figsize=(9, 11))

# 用 pcolormesh 而不是 imshow，这样更方便画白色网格线
mesh = ax.pcolormesh(
    np.arange(pmat.shape[1] + 1),
    np.arange(pmat.shape[0] + 1),
    pmat.values,
    cmap=cmap,
    norm=norm,
    edgecolors="white",   # 白色分隔线
    linewidth=0.3,        # 分隔线粗细，可调成 0.3 / 0.8
    shading="flat"
)

# 让第一行显示在上方，效果更像 pheatmap
ax.invert_yaxis()
ax.set_aspect('equal')
# 坐标轴标签放在每个格子的中心
ax.set_xticks(np.arange(pmat.shape[1]) + 0.5)
ax.set_xticklabels(pmat.columns, rotation=90,fontsize=10)

ax.set_yticks(np.arange(pmat.shape[0]) + 0.5)
ax.set_yticklabels(pmat.index,fontsize=10)

# 去掉刻度线本体
ax.tick_params(axis="both", length=0)

# 去掉边框
for spine in ax.spines.values():
    spine.set_visible(False)

# 手动指定 colorbar 位置
cax = fig.add_axes([1.0, 0.6, 0.02, 0.15])
cbar = fig.colorbar(mesh, cax=cax)

tick_values = [0.0, 0.5, 1.0]

cbar.set_ticks(tick_values)
cbar.set_ticklabels([f"{t:.1f}" for t in tick_values])
cbar.ax.tick_params(labelsize=10, length=3)
cbar.outline.set_visible(False)
cbar.minorticks_off()

plt.tight_layout()
plt.savefig(out_path+"heatmap.pdf", bbox_inches="tight")
plt.close()