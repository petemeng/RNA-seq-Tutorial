#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${1:-validation_run_downstream}"
SHEET="${BASE_DIR}/metadata/samplesheet.csv"
RAW_DIR="${BASE_DIR}/data/raw_data"
QUANT_DIR="${BASE_DIR}/data/quant"
REF_DIR="${BASE_DIR}/data/reference"
LOG_DIR="${BASE_DIR}/logs"

SALMON_BIN="${SALMON_BIN:-salmon}"
THREADS="${THREADS:-8}"
SALMON_INDEX="${SALMON_INDEX:-${REF_DIR}/salmon_index}"

mkdir -p "${QUANT_DIR}" "${LOG_DIR}"

if [[ ! -f "${SHEET}" ]]; then
  echo "[ERROR] Missing samplesheet: ${SHEET}" >&2
  exit 1
fi
if [[ ! -d "${SALMON_INDEX}" ]]; then
  echo "[ERROR] Missing salmon index: ${SALMON_INDEX}" >&2
  exit 1
fi

while IFS=',' read -r sample_id run genotype condition time replicate layout fastq1 fastq2 batch; do
  [[ "${sample_id}" == "sample_id" ]] && continue
  f1="${RAW_DIR}/$(basename "${fastq1}")"
  f2="${RAW_DIR}/$(basename "${fastq2}")"

  "${SALMON_BIN}" quant \
    -i "${SALMON_INDEX}" \
    -l A \
    -1 "${f1}" \
    -2 "${f2}" \
    --validateMappings \
    --gcBias \
    --seqBias \
    -p "${THREADS}" \
    -o "${QUANT_DIR}/${sample_id}" \
    > "${LOG_DIR}/salmon_quant_${sample_id}.log" 2>&1

  echo "salmon done: ${sample_id}"
done < "${SHEET}"

echo "[INFO] Quantification complete. Samples: $(ls "${QUANT_DIR}" | wc -l)"
