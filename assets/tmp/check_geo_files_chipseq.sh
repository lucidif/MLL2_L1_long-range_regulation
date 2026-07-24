#!/usr/bin/env bash
tsv="${1:-geo_sub_chipSEQ.tsv}"
out="${2:-missing_files.txt}"

# empty the output file
: > "$out"

awk -v OFS="\t" -v out="$out" '
BEGIN{ FS="\t" }
# Map header to column indices
NR==1{
  for(i=1;i<=NF;i++){ h[$i]=i }
  # required columns
  req = "GEOname1 GEOname2 originalName1 originalName2 fastq_path processed_path originalProcessed"
  n=split(req, r, " ")
  for(i=1;i<=n;i++){ if(!(r[i] in h)) { missingcols=missingcols r[i] " " } }
  if(length(missingcols)) {
    printf("ERROR: missing required columns: %s\n", missingcols) > "/dev/stderr"
    exit 2
  }
  next
}
# For each data row, compose full paths and test for existence
{
  fq1 = $(h["fastq_path"]) "/" $(h["originalName1"])
  fq2 = $(h["fastq_path"]) "/" $(h["originalName2"])
  bw  = $(h["processed_path"]) "/" $(h["originalProcessed"])

  # trim spaces (just in case)
  gsub(/^[ \t]+|[ \t]+$/, "", fq1); gsub(/^[ \t]+|[ \t]+$/, "", fq2); gsub(/^[ \t]+|[ \t]+$/, "", bw)

  # Keep a nice row header for context (printed once if anything is missing)
  row_header_printed = 0
  # Check FASTQ1
  if (system("[ -e \"" fq1 "\" ]")==0) {
    # ok
  } else {
    if (!row_header_printed) { print "ROW", NR, $(h["GEOname1"]), $(h["GEOname2"]) >> out; row_header_printed=1 }
    print "MISSING", "FASTQ1", fq1 >> out
  }
  # Check FASTQ2
  if (system("[ -e \"" fq2 "\" ]")==0) {
    # ok
  } else {
    if (!row_header_printed) { print "ROW", NR, $(h["GEOname1"]), $(h["GEOname2"]) >> out; row_header_printed=1 }
    print "MISSING", "FASTQ2", fq2 >> out
  }
  # Check bigWig only if non-empty originalProcessed
  if ($(h["originalProcessed"])!="") {
    if (system("[ -e \"" bw "\" ]")==0) {
      # ok
    } else {
      if (!row_header_printed) { print "ROW", NR, $(h["GEOname1"]), $(h["GEOname2"]) >> out; row_header_printed=1 }
      print "MISSING", "BIGWIG", bw >> out
    }
  }
  if (row_header_printed) { print "" >> out }
}
END{
  # noop
}
' "$tsv"

echo "Done. Missing files (if any) were written to: $out"
