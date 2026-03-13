from pathlib import Path
import pandas as pd

# =========================
# Paths
# =========================
source_folder = Path("./data/EnrichScoreMatrix/t_test/bias_result")
final_folder = source_folder / "final"
final_folder.mkdir(parents=True, exist_ok=True)

# =========================
# Step 1: Filter significant rows
# Keep rows where either p_value_left or p_value_right < 0.05
# =========================
for csv_path in source_folder.rglob("*.csv"):
    # Skip files already inside the final folder
    if final_folder in csv_path.parents:
        continue

    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()

    required_cols = {"p_value_left", "p_value_right"}
    if not required_cols.issubset(df.columns):
        continue

    filtered_df = df[
        (df["p_value_left"] < 0.05) | (df["p_value_right"] < 0.05)
    ].copy()

    if filtered_df.empty:
        print(f"[skip empty] {csv_path}")
        continue

    out_path = final_folder / csv_path.name
    filtered_df.to_csv(out_path, index=False)
    print(f"[wrote] {out_path}")

# =========================
# Step 2: Keep only CellType and the minimum p-value
# =========================
for csv_path in final_folder.glob("*.csv"):
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()

    required_cols = {"CellType", "p_value_left", "p_value_right"}
    if not required_cols.issubset(df.columns):
        print(f"[skip missing columns] {csv_path}")
        continue

    df["p_value"] = df[["p_value_left", "p_value_right"]].min(axis=1)
    out_df = df[["CellType", "p_value"]].copy()

    if out_df.empty:
        print(f"[skip empty] {csv_path}")
        continue

    out_df.to_csv(csv_path, index=False)
    print(f"[wrote minimal] {csv_path}")

print(f"Done: filtered results and wrote minimal p_value files into: {final_folder.as_posix()}")