#!/usr/bin/env bash
set -euo pipefail
usage() {
  cat <<'EOF'
Usage: check_geo_files.sh [TSV] [OUTPUT] [--paths-only] [--tsv]

TSV         Input TSV (default: geo_sub_chipSEQ.tsv)
OUTPUT      Output report (default: missing_files.txt)

Supports two schemas:
  RAW-ONLY:  GEOname1, GEOname2, originalName1, originalName2, fastq_path
  FULL:      + processed_path, GEOprocessed, originalProcessed

Options:
  --paths-only   Print only missing paths (one per line)
  --tsv          Output a TSV with columns: row, geoname1, geoname2, type, path
Note: this version avoids AWK functions for maximum portability (BSD/mawk).
EOF
}
tsv="geo_sub_chipSEQ.tsv"; out="missing_files.txt"; paths_only=0; as_tsv=0
args=()
for a in "${@:-}"; do
  case "$a" in
    -h|--help) usage; exit 0;;
    --paths-only) paths_only=1;;
    --tsv) as_tsv=1;;
    *) args+=("$a");;
  esac
done
set -- "${args[@]:-}"
(( $#>=1 )) && tsv="$1"
(( $#>=2 )) && out="$2"
(( $#>2 )) && { usage; exit 2; }
(( paths_only==0 )) && : > "$out"

awk -v OFS="\t" -v out="$out" -v paths_only="$paths_only" -v as_tsv="$as_tsv" '
BEGIN{ FS="\t" }
NR==1{
  for(i=1;i<=NF;i++){ gsub(/\r/,"",$i); gsub(/^[ \t]+|[ \t]+$/, "", $i); h[$i]=i }
  rawreq="GEOname1 GEOname2 originalName1 originalName2 fastq_path"
  n=split(rawreq,r," "); for(i=1;i<=n;i++){ if(!(r[i] in h)) miss=miss r[i] " " }
  if(length(miss)){ print "ERROR: missing raw columns: " miss > "/dev/stderr"; exit 2 }
  has_proc = ("processed_path" in h) && ("GEOprocessed" in h) && ("originalProcessed" in h)
  if(as_tsv && paths_only==0){ print "row","geoname1","geoname2","type","path" > out }
  next
}
{
  fqdir=$(h["fastq_path"]); on1=$(h["originalName1"]); on2=$(h["originalName2"])
  geon1=$(h["GEOname1"]); geon2=$(h["GEOname2"])

  raw1 = (length(on1)? fqdir "/" on1 : "")
  raw2 = (length(on2)? fqdir "/" on2 : "")

  proc=""
  if (has_proc){
    pp=$(h["processed_path"]); op=$(h["originalProcessed"])
    if (op ~ /^\// || op ~ /\//)      proc = op;
    else if (pp ~ /\// && length(op)) proc = pp "/" op;
  }

  gsub(/^[ \t]+|[ \t]+$/, "", raw1); gsub(/^[ \t]+|[ \t]+$/, "", raw2); gsub(/^[ \t]+|[ \t]+$/, "", proc)

  header_printed=0

  # check RAW1
  if (length(raw1) && system("[ -e \"" raw1 "\" ]")!=0){
    if(paths_only){ print raw1 }
    else if(as_tsv){ print NR, geon1, geon2, "FASTQ1", raw1 >> out }
    else{
      if(!header_printed){ print "ROW", NR, geon1, geon2 >> out; header_printed=1 }
      print "MISSING", "FASTQ1", raw1 >> out
    }
  }

  # check RAW2
  if (length(raw2) && system("[ -e \"" raw2 "\" ]")!=0){
    if(paths_only){ print raw2 }
    else if(as_tsv){ print NR, geon1, geon2, "FASTQ2", raw2 >> out }
    else{
      if(!header_printed){ print "ROW", NR, geon1, geon2 >> out; header_printed=1 }
      print "MISSING", "FASTQ2", raw2 >> out
    }
  }

  # check PROCESSED if present
  if (length(proc) && system("[ -e \"" proc "\" ]")!=0){
    if(paths_only){ print proc }
    else if(as_tsv){ print NR, geon1, geon2, "PROCESSED", proc >> out }
    else{
      if(!header_printed){ print "ROW", NR, geon1, geon2 >> out; header_printed=1 }
      print "MISSING", "PROCESSED", proc >> out
    }
  }

  if (header_printed && !as_tsv && !paths_only){ print "" >> out }
}
' "$tsv"

if (( paths_only==0 )); then echo "Done. Report written to: $out"; fi
