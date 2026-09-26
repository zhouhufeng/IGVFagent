"""GSM3017261 MATLAB DGE -> AnnData with the authors' deposited labels (per their GEO gist)."""
import gzip, shutil, os, json
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, anndata as ad
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "Data", "rosenberg2018"))
if not os.path.exists("GSM3017261_150000_CNS_nuclei.mat"):
    with gzip.open("GSM3017261_150000_CNS_nuclei.mat.gz") as fi, open("GSM3017261_150000_CNS_nuclei.mat", "wb") as fo:
        shutil.copyfileobj(fi, fo)
d = sio.loadmat("GSM3017261_150000_CNS_nuclei.mat")
for k, v in d.items():
    if not k.startswith("__"):
        print(k, type(v).__name__, getattr(v, "shape", None), getattr(v, "dtype", None))
X = sp.csr_matrix(d["DGE"])
st = lambda k: pd.Series(d[k]).str.strip().values
obs = pd.DataFrame({"sample_type": st("sample_type"), "cluster_assignment": st("cluster_assignment"),
                    "spinal_cluster_assignment": st("spinal_cluster_assignment")},
                   index=[f"n{int(b)}" for b in np.ravel(d["barcodes"])])
var = pd.DataFrame(index=pd.Series(d["genes"]).str.strip().values)
var.index = pd.Index(var.index).astype(str)
a = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)
a.obs_names_make_unique(); a.var_names_make_unique()
a.write_h5ad("rosenberg_cns_full.h5ad", compression="gzip")
ca = a.obs.cluster_assignment.value_counts()
sc = a.obs.spinal_cluster_assignment.value_counts()
out = {"n_nuclei": int(a.n_obs), "n_genes": int(a.n_vars),
       "sample_type": a.obs.sample_type.value_counts().to_dict(),
       "cluster_assignment": ca.to_dict(), "spinal_cluster_assignment": sc.to_dict()}
json.dump(out, open("atlas_labels_summary.json", "w"), indent=1)
print(json.dumps({k: (v if not isinstance(v, dict) else len(v)) for k, v in out.items()}))
