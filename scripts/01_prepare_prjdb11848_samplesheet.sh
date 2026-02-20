#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${1:-validation_run_downstream}"
META_DIR="${BASE_DIR}/metadata"
mkdir -p "${META_DIR}"

echo "[INFO] Fetching PRJDB11848 metadata from ENA ..."
wget -qO "${META_DIR}/prjdb11848.tsv" \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB11848&result=read_run&fields=run_accession,sample_title,scientific_name,library_layout,fastq_ftp"

echo "[INFO] Building samplesheet.csv ..."
awk -F"\t" '
BEGIN{
  OFS=",";
  print "sample_id,run,genotype,condition,time,replicate,layout,fastq_1,fastq_2,batch";
}
NR>1{
  run=$1; title=$2; layout=$4; ftp=$5;
  n=split(title,a,"_");
  genotype=a[1]; condition=a[2]; time=a[3]; rep=a[4];
  split(ftp,f,";");
  gsub(/^ftp\.sra\.ebi\.ac\.uk\//,"",f[1]);
  gsub(/^ftp\.sra\.ebi\.ac\.uk\//,"",f[2]);
  fq1="http://ftp.sra.ebi.ac.uk/" f[1];
  fq2="http://ftp.sra.ebi.ac.uk/" f[2];
  sample_id=run;
  print sample_id,run,genotype,condition,time,rep,layout,fq1,fq2,"B1";
}
' "${META_DIR}/prjdb11848.tsv" > "${META_DIR}/samplesheet.csv"

echo "[INFO] samplesheet rows: $(wc -l < "${META_DIR}/samplesheet.csv")"
echo "[INFO] Preview:"
head -n 8 "${META_DIR}/samplesheet.csv"
