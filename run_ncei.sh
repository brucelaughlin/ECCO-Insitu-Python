#!/usr/bin/env bash
#
# Controller for the NCEI post-processing chain.
#
# Runs NCEI.py, teeing all stdout+stderr to a timestamped log under z_logs/.
#
# DEFAULTS (run with no arguments):
#   input  = ${DEFAULT_INPUT_DIR}      (see below)
#   output = ${DEFAULT_OUTPUT_BASE}_<YYYYMMDD_HHMMSS>
# The timestamp is appended to the output base automatically, so no run ever
# overwrites a previous one.
#
# OPTIONAL OVERRIDES:
#   -i <input_dir>    use a different input directory
#   -d <dest_dir>     use this EXACT output directory (no timestamp appended;
#                     you are being explicit, so you own the name)
#   -c                run step11_profiles_compact.py (compact is skipped by default).
#
# Examples:
#   ./run_ncei.sh                                   # defaults + timestamp, no compact
#   ./run_ncei.sh -i /path/to/other_inputs          # other input, default+ts output
#   ./run_ncei.sh -d /path/to/exact_output_dir      # default input, exact output
#   ./run_ncei.sh -c                                # also run compact step
#   ./run_ncei.sh -n 8 -c                           # 8 workers + compact
#
# The exit code of NCEI.py is preserved (via PIPESTATUS).

set -euo pipefail

# --- defaults (edit these to change what a bare run does) ---
DEFAULT_INPUT_DIR="/Users/brucel/ecco/yip/profile_data/Interp_Profiles"
DEFAULT_OUTPUT_BASE="/Users/brucel/ecco/yip/processed_by_NCEI_profile_data/profile_files_NCEI_processed"

# --- locate ourselves so paths work regardless of CWD ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/z_logs"
mkdir -p "${LOG_DIR}"

STAMP="$(date +%Y%m%d_%H%M%S)"

# --- parse optional overrides ---
INPUT_DIR=""
DEST_DIR=""
N_WORKERS=""
RUN_COMPACT=0
while getopts "i:d:n:ch" opt; do
  case "${opt}" in
    i) INPUT_DIR="${OPTARG}" ;;
    d) DEST_DIR="${OPTARG}" ;;
    n) N_WORKERS="${OPTARG}" ;;
    c) RUN_COMPACT=1 ;;
    h) grep '^#' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Usage: $0 [-i input_dir] [-d dest_dir] [-n n_workers] [-c]"; exit 2 ;;
  esac
done

# --- apply defaults where not overridden ---
if [[ -z "${INPUT_DIR}" ]]; then
  INPUT_DIR="${DEFAULT_INPUT_DIR}"
fi
if [[ -z "${DEST_DIR}" ]]; then
  # no -d given: default base + timestamp so nothing is ever overwritten
  DEST_DIR="${DEFAULT_OUTPUT_BASE}_${STAMP}"
fi

LOG_FILE="${LOG_DIR}/ncei_log_${STAMP}.log"

# --- safety checks ---
if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "ERROR: input dir does not exist: ${INPUT_DIR}" >&2
  exit 2
fi
if [[ -e "${DEST_DIR}" ]]; then
  echo "ERROR: output dir already exists, refusing to overwrite: ${DEST_DIR}" >&2
  echo "       (pass a different -d, or remove it first)" >&2
  exit 2
fi

# --- banner ---
{
echo "==============================================================="
echo "NCEI chain run started: $(date)"
echo "Input dir:  ${INPUT_DIR}"
echo "Output dir: ${DEST_DIR}"
echo "Log file:   ${LOG_FILE}"
echo "==============================================================="
echo ""
} | tee "${LOG_FILE}"

# python3 -u keeps output unbuffered so the log fills in real time; 2>&1 folds
# stderr into stdout before the pipe so warnings/tracebacks land in the log too.
N_WORKERS_ARG=""
if [[ -n "${N_WORKERS}" ]]; then N_WORKERS_ARG="-n ${N_WORKERS}"; fi
python3 -u "${SCRIPT_DIR}/NCEI.py" -i "${INPUT_DIR}" -d "${DEST_DIR}" ${N_WORKERS_ARG} 2>&1 | tee -a "${LOG_FILE}"
STATUS="${PIPESTATUS[0]}"

{
echo ""
echo "==============================================================="
echo "NCEI chain run finished: $(date) (exit code ${STATUS})"
echo "Output dir: ${DEST_DIR}"
echo "==============================================================="
} | tee -a "${LOG_FILE}"

if [[ "${RUN_COMPACT}" -eq 1 && "${STATUS}" -eq 0 ]]; then
  {
  echo ""
  echo "==============================================================="
  echo "Step 11 (compact): rewriting output files in-place"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"

  COMPACT_FAILED=0
  while IFS= read -r -d '' NC_FILE; do
    STEM="${NC_FILE%.nc}"
    TMP_OUT="${STEM}_compact.nc"
    python3 -u "${SCRIPT_DIR}/step11_profiles_compact.py" \
      "${NC_FILE}" -o "${TMP_OUT}" --complevel 2 --overwrite 2>&1 | tee -a "${LOG_FILE}"
    COMPACT_STATUS="${PIPESTATUS[0]}"
    if [[ "${COMPACT_STATUS}" -eq 0 ]]; then
      mv "${TMP_OUT}" "${NC_FILE}"
    else
      echo "WARNING: step11 failed for ${NC_FILE}; leaving original in place" | tee -a "${LOG_FILE}"
      COMPACT_FAILED=1
    fi
  done < <(find "${DEST_DIR}" -name '*.nc' -print0)

  {
  echo ""
  echo "==============================================================="
  echo "Step 11 (compact) finished: $(date)"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"

  if [[ "${COMPACT_FAILED}" -ne 0 ]]; then
    STATUS=1
  fi
fi

exit "${STATUS}"
