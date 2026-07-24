#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Uso: $0 <file.clpy> [outdir]"
  exit 1
fi

CLPY="$1"
OUTDIR="${2:-clpy_extract}"
mkdir -p "$OUTDIR"

python3 - "$CLPY" "$OUTDIR" <<'PY'
import os, sys
import numpy as np
import pandas as pd

fn = sys.argv[1]
outdir = sys.argv[2]

def safe(s: str) -> str:
    s = str(s)
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in s)[:180]

def is_matrix_like(x):
    try:
        a = np.asarray(x)
    except Exception:
        return False
    return a.ndim == 2 and a.size > 0

def write_matrix(mat, name):
    mat = np.asarray(mat, dtype=float)
    tsv = os.path.join(outdir, f"{name}.tsv")
    npy = os.path.join(outdir, f"{name}.npy")
    np.savetxt(tsv, mat, delimiter="\t")
    np.save(npy, mat)
    return tsv, npy, mat

with pd.HDFStore(fn, "r") as store:
    keys = store.keys()
    print("HDF5 keys:", keys)
    if "/annotation" not in keys:
        raise SystemExit("ERRORE: non trovo la key /annotation nel .clpy")

    df = store["/annotation"]

print("== annotation loaded ==")
print("rows, cols:", df.shape)
print("columns:", list(df.columns))

# trova colonne che contengono matrici 2D
matrix_cols = []
for col in df.columns:
    # prova sul primo valore non-null
    v0 = None
    for _, v in df[col].items():
        if v is None:
            continue
        # filtra NaN scalari
        if isinstance(v, float) and np.isnan(v):
            continue
        v0 = v
        break
    if v0 is not None and is_matrix_like(v0):
        matrix_cols.append(col)

if not matrix_cols:
    raise SystemExit("ERRORE: nessuna colonna in /annotation sembra contenere matrici 2D.")

print("== matrix-like columns found ==", matrix_cols)

summary_path = os.path.join(outdir, "summary.tsv")
with open(summary_path, "w") as f:
    f.write("\t".join(["row_id","col","n","shape","min","max","mean","out_tsv","out_npy"]) + "\n")

    for row_id, row in df.iterrows():
        n_val = row["n"] if "n" in row.index else ""
        for col in matrix_cols:
            mat = row[col]
            if mat is None:
                continue
            name = safe(f"{col}_row{row_id}")
            out_tsv, out_npy, m = write_matrix(mat, name)
            f.write("\t".join([
                str(row_id),
                str(col),
                str(n_val),
                f"{m.shape[0]}x{m.shape[1]}",
                f"{np.nanmin(m):.6g}",
                f"{np.nanmax(m):.6g}",
                f"{np.nanmean(m):.6g}",
                os.path.basename(out_tsv),
                os.path.basename(out_npy),
            ]) + "\n")

# se c’è num e n, calcola anche la matrice media num/n (quella più “plot-like”)
if "num" in matrix_cols and "n" in df.columns:
    for row_id, row in df.iterrows():
        num = row["num"]
        n = row["n"]
        if num is None or n in ("", None) or (isinstance(n, float) and np.isnan(n)) or float(n) == 0:
            continue
        avg = np.asarray(num, dtype=float) / float(n)
        name = safe(f"avg_num_over_n_row{row_id}")
        write_matrix(avg, name)

print("== written ==")
print("outdir:", outdir)
print("summary:", summary_path)
PY

echo "Done. Output in: $OUTDIR"
echo " - summary.tsv"
echo " - <col>_row<id>.tsv/.npy (matrici raw)"
echo " - avg_num_over_n_row<id>.tsv/.npy (se presenti num e n)"
