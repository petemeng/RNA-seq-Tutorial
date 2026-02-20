#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${1:-validation_run_downstream}"
SHEET="${BASE_DIR}/metadata/samplesheet.csv"
RAW_DIR="${BASE_DIR}/data/raw_data"
LOG_DIR="${BASE_DIR}/logs"
LOG_FILE="${LOG_DIR}/download_fastq.log"

mkdir -p "${RAW_DIR}" "${LOG_DIR}"

if [[ ! -f "${SHEET}" ]]; then
  echo "[ERROR] Missing samplesheet: ${SHEET}" >&2
  exit 1
fi

: > "${LOG_FILE}"
echo "[INFO] Download start: $(date -Iseconds)" | tee -a "${LOG_FILE}"

awk -F',' 'NR>1 {print $8"\n"$9}' "${SHEET}" | while read -r url; do
  [[ -z "${url}" ]] && continue
  echo "[$(date -Iseconds)] ${url}" | tee -a "${LOG_FILE}"
  ok=0
  for i in $(seq 1 10); do
    if wget -c -nv --tries=1 --timeout=30 "${url}" -P "${RAW_DIR}/" >> "${LOG_FILE}" 2>&1; then
      ok=1
      break
    fi
    echo "retry ${i} failed: ${url}" | tee -a "${LOG_FILE}"
    sleep 2
  done
  if [[ "${ok}" -ne 1 ]]; then
    echo "[ERROR] failed: ${url}" | tee -a "${LOG_FILE}"
    exit 1
  fi
done

echo "[INFO] Download done: $(date -Iseconds)" | tee -a "${LOG_FILE}"
echo "[INFO] FASTQ file count: $(ls "${RAW_DIR}"/*.fastq.gz | wc -l)"
