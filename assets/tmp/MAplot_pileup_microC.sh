chmod +x ./git/Lara_MLL2/bin/extract_clpy.sh

python3 - <<'PY'
import pandas as pd
fn="outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/win500_anch3_aLp_WT_day0.Dd.clpy"
with pd.HDFStore(fn, "r") as s:
    print("HDF5 keys:")
    for k in s.keys():
        print(" ", k, "->", type(s.get(k)))
PY

ls outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/*.clpy

ls outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/*.clpy
# outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/win500_anch2_aLp_KO_day0.Dd.clpy
# outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/win500_anch2_aLp_WT_day0.Dd.clpy
# outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/win500_anch3_aLp_KO_day0.Dd.clpy
# outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/win500_anch3_aLp_WT_day0.Dd.clpy
# outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/win500_anch3_aLp_WT_day4.Dd.clpy

ls outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/*.clpy
# outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/win500_anch2_aLp_KO_day4.Dd.clpy
# outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/win500_anch2_aLp_WT_day4.Dd.clpy
# outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/win500_anch3_aLp_KO_day4.Dd.clpy
# outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/win500_anch3_aLp_WT_day4.Dd.clpy

bash ./git/Lara_MLL2/bin/extract_clpy.sh \
  outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/win500_anch3_aLp_WT_day0.Dd.clpy \
  outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/WT_day0.Dd_extracted

bash ./git/Lara_MLL2/bin/extract_clpy.sh \
  outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/win500_anch3_aLp_WT_day0.Dd.clpy \
  outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/WT_day0.Dd_extracted

bash ./git/Lara_MLL2/bin/extract_clpy.sh \
  outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/win500_anch3_aLp_WT_day4.Dd.clpy \
  outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/WT_day4.Dd_extracted

bash ./git/Lara_MLL2/bin/extract_clpy.sh \
  outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/win500_anch3_aLp_KO_day4.Dd.clpy \
  outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/KO_day4.Dd_extracted

head outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/extracted/avg_num_over_n_row0.tsv

python3 - <<'PY'
import pandas as pd

fn="outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/win500_anch3_aLp_KO_day4.Dd.clpy"
df = pd.read_hdf(fn, "/annotation")
cols = [c for c in ["group","n","control_n","num","control_num","features","subset"] if c in df.columns]
print(df[cols])
PY


ls outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/WT_day0.Dd_extracted/avg_num_over_n_row0.tsv

mkdir outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot
mkdir outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp

Rscript git/Lara_MLL2/bin/microC_MAplot_pileup.R
