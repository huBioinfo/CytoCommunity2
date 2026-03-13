import os
import glob
import numpy as np
import pandas as pd
from collections import defaultdict

TCN_num = 10

input_folder = "./TNBC_Input"
target_folder = "./data/EnrichScoreMatrix"
score_folder = os.path.join(target_folder, "Score")
ttest_folder = os.path.join(target_folder, "t_test")

os.makedirs(score_folder, exist_ok=True)
os.makedirs(ttest_folder, exist_ok=True)


# 读取样本名
region_list = pd.read_csv(
    os.path.join(input_folder, "ImageNameList.txt"),
    sep="\t",
    header=None,
    names=["Image"]
)["Image"].tolist()

label_to_regions = defaultdict(list)

for region_name in region_list:
    label_path = os.path.join(input_folder, f"{region_name}_GraphLabel.txt")
    graph_label = int(np.loadtxt(label_path, dtype="int64", delimiter="\t"))
    label_to_regions[graph_label].append(region_name)

for label, regions in sorted(label_to_regions.items()):
    print(f"Label {label}: {regions} (n={len(regions)})")


for region_name in region_list:
    label_path = os.path.join(input_folder, f"{region_name}_GraphLabel.txt")
    graph_label = int(np.loadtxt(label_path, dtype="int64", delimiter="\t"))

    old_csv = os.path.join(target_folder, f"{region_name}_EnrichScoreMatrix.csv")
    if not os.path.isfile(old_csv):
        continue

    new_csv = os.path.join(target_folder, f"{graph_label}_{region_name}_EnrichScoreMatrix.csv")
    os.rename(old_csv, new_csv)

celltype_path = os.path.join(target_folder, "CellTypeVec_List.csv")
celltype_list = pd.read_csv(celltype_path).iloc[:, 0].tolist()
print(celltype_list)

matrix_files = glob.glob(os.path.join(target_folder, "*_EnrichScoreMatrix.csv"))


# 为每个 TCN 生成长表 EnrichmentScore.csv
for i in range(TCN_num):
    rows = []

    for file in matrix_files:
        data = pd.read_csv(file, header=None)
        file_name = os.path.basename(file)
        condition = file_name.split("_")[0]   # 文件名第一个字段作为 condition

        for j, cell_type in enumerate(celltype_list):
            enrichment_sco = data.iloc[j, i]
            rows.append({
                "Condition": condition,
                "EnrichSco": enrichment_sco,
                "CellType": cell_type
            })

    tcn_df = pd.DataFrame(rows)

    # 去掉某个 CellType 在所有样本中都为 0 的情况
    tcn_df = tcn_df.groupby("CellType").filter(lambda x: (x["EnrichSco"] != 0).sum() > 0)

    output_file = os.path.join(score_folder, f"TCN{i+1}_EnrichmentScore.csv")
    tcn_df.to_csv(output_file, index=False)
    print(f"[SAVE] {output_file}")

# 转成 t-test 用的向量格式
for tcn in range(1, TCN_num + 1):
    input_file = os.path.join(score_folder, f"TCN{tcn}_EnrichmentScore.csv")
    output_file = os.path.join(ttest_folder, f"TCN{tcn}.csv")

    data = pd.read_csv(input_file)
    result = []

    for cell_type, group in data.groupby("CellType"):
        for condition in group["Condition"].unique():
            enrich_sco_vector = tuple(
                float(v) for v in group[group["Condition"] == condition]["EnrichSco"].values
            )
            enrich_sco_vector = tuple(0 if v == -0.0 else v for v in enrich_sco_vector)

            result.append({
                "CellType": cell_type,
                "Condition": condition,
                "EnrichSco": enrich_sco_vector
            })

    result_df = pd.DataFrame(result)
    result_df.to_csv(output_file, index=False)
    print(f"[SAVE] {output_file}")