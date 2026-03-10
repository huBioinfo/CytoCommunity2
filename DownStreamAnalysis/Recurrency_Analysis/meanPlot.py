import matplotlib as mpl
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

out_path="plot/OverlapCoefficient"
os.makedirs(out_path,exist_ok=True)

df = pd.read_csv('data/OverlapCoefficient/config.csv')
df['TCN'] = df['TCN'].astype(str).str.strip()
df['OverlapCoefficient'] = pd.to_numeric(df['OverlapCoefficient'], errors='coerce')

pooled = pd.DataFrame({
    'Group': ['All CNs'] * len(df),
    'OverlapCoefficient': df['OverlapCoefficient'].values
})


mean_of_means = (
    df.groupby('TCN')['OverlapCoefficient'].mean().mean()
)

plt.figure(figsize=(5, 6), dpi=1200)
ax = sns.violinplot(
    x='Group', y='OverlapCoefficient', data=pooled,
    inner=None, color='white', linewidth=1, cut=0
)

for pc in ax.collections:
    try:
        pc.set_edgecolor('black')
        pc.set_facecolor('white')
    except Exception:
        pass

ax.set_ylim(0, 1)
ax.set_xlabel('')
ax.set_ylabel('Overlap Coefficient')

x0, x1 = -0.35, 0.35   
ax.plot([x0, x1], [mean_of_means, mean_of_means], color='red', linewidth=2)

ax.text(0.42, mean_of_means, f'{mean_of_means:.3f}', color='red', va='center')

plt.tight_layout()
plt.savefig(out_path+'/OC_mean.pdf', dpi=1200)
plt.close()
