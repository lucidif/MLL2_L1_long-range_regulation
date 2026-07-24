#!/usr/bin/env bash
set -euo pipefail

usage(){ cat <<'EOF'
Usage:
  verify_geo_md5.sh <TSV> <MD5_REPORT> [--tsv] [--quiet]

Compares MD5 of *source* files (computed now) vs MD5s recorded in MD5_REPORT
(which are of the *destination* files written by organize_geo_files.sh).

MD5_REPORT accepted formats (per TSV row):
  (A) 3 cols : md5_1  md5_2  md5_P
  (B) 6 cols : name1  md5_1  name2  md5_2  nameP  md5_P   (with or without header)

Options:
  --tsv      Output TSV with columns:
             row  geoname1  geoname2  kind  src_path  src_md5  dst_md5  status
  --quiet    Print only mismatches/missing (suppresses matches)
EOF
}

(( $# < 2 )) && { usage; exit 2; }

TSV="$1"; MD5="$2"; shift 2
tsv_out=0; quiet=0
for a in "${@:-}"; do
  [[ -z "${a:-}" ]] && continue
  case "$a" in
    --tsv) tsv_out=1;;
    --quiet) quiet=1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $a" >&2; usage; exit 2;;
  esac
done

md5_of(){
  if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'
  elif command -v md5 >/dev/null 2>&1; then md5 -r "$1" | awk '{print $1}'
  else echo "ERROR: neither md5sum nor md5 found" >&2; exit 3
  fi
}

# Read entire MD5 report
mapfile -t MD5_LINES < "$MD5"

# If first line looks like a header in 6-col form, drop it
if (( ${#MD5_LINES[@]} > 0 )); then
  IFS=$'\t' read -r c1 c2 c3 c4 c5 c6 <<< "${MD5_LINES[0]}"
  if [[ "$c1" == "name1" && "$c2" == "md5_1" ]]; then
    MD5_LINES=("${MD5_LINES[@]:1}")
  fi
fi

# Helper: extract destination MD5s for a given row index from either 3-col or 6-col row
get_md5s() {
  local line="$1" m1="" m2="" mp=""
  IFS=$'\t' read -r a b c d e f <<< "$line"
  if [[ -n "${f:-}" ]]; then
    # 6 columns (names + md5s)
    m1="$b"; m2="$d"; mp="$f"
  else
    # 3 columns (just md5s)
    m1="$a"; m2="$b"; mp="$c"
  fi
  printf '%s\t%s\t%s\n' "$m1" "$m2" "$mp"
}

(( tsv_out==1 )) && printf "row\tgeoname1\tgeoname2\tkind\tsrc_path\tsrc_md5\tdst_md5\tstatus\n"

status=0; rowidx=0

awk -v FS="\t" -v OFS="\t" '
NR==1{
  for(i=1;i<=NF;i++){ gsub(/\r/,"",$i); gsub(/^[ \t]+|[ \t]+$/, "", $i); h[$i]=i }
  req="GEOname1 GEOname2 originalName1 originalName2 fastq_path"
  n=split(req,r," "); for(i=1;i<=n;i++){ if(!(r[i] in h)) miss=miss r[i] " " }
  if(length(miss)){ print "ERROR: missing required columns: " miss > "/dev/stderr"; exit 2 }
  has_proc = ("processed_path" in h) && ("GEOprocessed" in h) && ("originalProcessed" in h)
  next
}
{
  fqdir=$(h["fastq_path"]); on1=$(h["originalName1"]); on2=$(h["originalName2"])
  g1=$(h["GEOname1"]); g2=$(h["GEOname2"])
  gsub(/\r/,"",fqdir); gsub(/\r/,"",on1); gsub(/\r/,"",on2); gsub(/\r/,"",g1); gsub(/\r/,"",g2)
  gsub(/^[ \t]+|[ \t]+$/, "", fqdir); gsub(/^[ \t]+|[ \t]+$/, "", on1); gsub(/^[ \t]+|[ \t]+$/, "", on2)
  gsub(/^[ \t]+|[ \t]+$/, "", g1); gsub(/^[ \t]+|[ \t]+$/, "", g2)

  src1=(length(on1)? fqdir "/" on1 : ""); src2=(length(on2)? fqdir "/" on2 : "")

  proc_src=""
  if(has_proc){
    pp=$(h["processed_path"]); op=$(h["originalProcessed"])
    gsub(/\r/,"",pp); gsub(/\r/,"",op); gsub(/^[ \t]+|[ \t]+$/, "", pp); gsub(/^[ \t]+|[ \t]+$/, "", op)
    if (op ~ /^\// || op ~ /\//)      proc_src = op;
    else if (pp ~ /\// && length(op)) proc_src = pp "/" op;
  }

  print NR, g1, g2, src1, src2, proc_src, has_proc
}
' "$TSV" | while IFS=$'\t' read -r ROWNO GEO1 GEO2 SRC1 SRC2 PROCSRC HASP; do
  # pull destination md5s for this row
  md5_1=""; md5_2=""; md5_p=""
  if (( rowidx < ${#MD5_LINES[@]} )); then
    IFS=$'\t' read -r md5_1 md5_2 md5_p <<< "$(get_md5s "${MD5_LINES[$rowidx]}")"
  fi

  for KIND in FASTQ1 FASTQ2 PROCESSED; do
    [[ "$KIND" == "PROCESSED" && "$HASP" != "1" ]] && continue
    case "$KIND" in
      FASTQ1) SRC="$SRC1"; DMD5="$md5_1";; FASTQ2) SRC="$SRC2"; DMD5="$md5_2";; PROCESSED) SRC="$PROCSRC"; DMD5="$md5_p";;
    esac

    # Skip if both sides are empty (no file expected)
    [[ -z "${SRC:-}" && -z "${DMD5:-}" ]] && continue

    # Source MD5 (compute now)
    SMISSING=0
    if [[ -n "${SRC:-}" && -e "$SRC" ]]; then SMD5="$(md5_of "$SRC")"; else SMD5=""; SMISSING=1; fi

    # Status
    if [[ -z "${DMD5:-}" ]]; then ST="MISSING_DST_HASH"
    elif (( SMISSING==1 )); then ST="MISSING_SRC_FILE"
    elif [[ "$SMD5" == "$DMD5" ]]; then ST="MATCH"
    else ST="MISMATCH"; fi

    if (( tsv_out==1 )); then
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$ROWNO" "$GEO1" "$GEO2" "$KIND" "${SRC:-}" "${SMD5:-}" "${DMD5:-}" "$ST"
    else
      (( quiet==0 || ST != "MATCH" )) && {
        printf "[row %s] %-9s | %s\n  src: %s\n  md5 src: %s\n  md5 dst: %s\n\n" "$ROWNO" "$KIND" "$ST" "${SRC:-<empty>}" "${SMD5:-<none>}" "${DMD5:-<none>}"
      }
    fi

    [[ "$ST" == "MISMATCH" || "$ST" == "MISSING_SRC_FILE" ]] && status=5 || true
  done

  rowidx=$((rowidx+1))
done

exit ${status:-0}
