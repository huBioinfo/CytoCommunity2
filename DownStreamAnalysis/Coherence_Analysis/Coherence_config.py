import os, glob
import pandas as pd
import anndata as ad
import numpy as np
from scipy.spatial import distance_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors


result_table_dir = "Step4_Output/ResultTable_File"
output_path = "data/Cohernece"
os.makedirs(output_path,exist_ok=True)
outputfile=os.path.join(output_path,"config.csv")

from sklearn.neighbors import NearestNeighbors

def compute_CHAOS_fast(clusterlabel, location):
    clusterlabel = np.asarray(clusterlabel)
    location = StandardScaler().fit_transform(np.asarray(location))

    chaos_vals = []

    for k in np.unique(clusterlabel):
        pts = location[clusterlabel == k]
        if len(pts) <= 1:
            continue

        nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree').fit(pts)
        dists, _ = nbrs.kneighbors(pts)
        chaos_vals.append(np.sum(dists[:, 1]))

    return np.sum(chaos_vals) / len(clusterlabel)

def compute_PAS_fast(clusterlabel, location, k=10):
    clusterlabel = np.asarray(clusterlabel)
    location = np.asarray(location)

    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='kd_tree').fit(location)
    dists, inds = nbrs.kneighbors(location)
    neigh_inds = inds[:, 1:]
    diff = (clusterlabel[neigh_inds] != clusterlabel[:, None])
    pas_flags = (diff.sum(axis=1) > (k / 2)).astype(int)

    return pas_flags.mean()


def fx_1NN(i, location_in):
    location_in = np.array(location_in)
    dist_array = distance_matrix(location_in[i, :][None, :], location_in)[0, :]
    dist_array[i] = np.inf
    return np.min(dist_array)


def fx_kNN(i, location_in, k, cluster_in):
    location_in = np.array(location_in)
    cluster_in = np.array(cluster_in)
    dist_array = distance_matrix(location_in[i, :][None, :], location_in)[0, :]
    dist_array[i] = np.inf
    ind = np.argsort(dist_array)[:k]
    cluster_use = np.array(cluster_in)
    if np.sum(cluster_use[ind] != cluster_in[i]) > (k / 2):
        return 1
    else:
        return 0


def compute_CHAOS(adata, pred_key, spatial_key='spatial'):
    return compute_CHAOS_fast(adata.obs[pred_key], adata.obsm[spatial_key])

def compute_PAS(adata, pred_key, spatial_key='spatial'):
    return compute_PAS_fast(adata.obs[pred_key], adata.obsm[spatial_key])



rows = []
for fp in glob.glob(os.path.join(result_table_dir, "*.csv")):
    sample = os.path.basename(fp).replace(".csv", "").replace("ResultTable_", "")
    df = pd.read_csv(fp)

    df = df.loc[~df["TCN_Label"].isna()].copy()
    coords = df[["x_coordinate", "y_coordinate"]].to_numpy(float)
    adata = ad.AnnData(X=coords, obs=pd.DataFrame({"pred": df["TCN_Label"].astype(str).values}))
    adata.obsm["spatial"] = coords

    chaos = float(compute_CHAOS(adata, "pred"))
    pas = float(compute_PAS(adata, "pred"))

    rows.append([sample, chaos, pas])
    print(f"[{sample}] CHAOS={chaos:.4f}, PAS={pas:.4f}")


df_out = pd.DataFrame(rows, columns=["Sample", "CHAOS", "PAS"]).round(4)


mean_row = pd.DataFrame([["Average", df_out["CHAOS"].mean(), df_out["PAS"].mean()]],
                        columns=df_out.columns).round(4)
df_out = pd.concat([df_out, mean_row], ignore_index=True)

df_out.to_csv(outputfile, index=False)
print(f"[SAVE] {outputfile}")
