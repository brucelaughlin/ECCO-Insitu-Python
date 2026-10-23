#!/usr/bin/env bash
#
# Controller for the NCEI post-processing chain.
#
# Runs NCEI.py, teeing all stdout+stderr to a timestamped log under z_logs/.
#
# SETTINGS come from config/defaults.yaml, config/local.yaml (this machine's data
# locations; copy config/local.example.yaml) and any -C files, merged in that order.
# Defaults for a bare run:
#   input  = run.input_dir
#   output = run.output_base + _<YYYYMMDD_HHMMSS>
# The timestamp is appended to the output base automatically, so no run ever
# overwrites a previous one. The flags below override the settings for one run.
#
# OPTIONAL OVERRIDES:
#   -C <file.yaml>    extra settings file (repeatable; later files win)
#   -s key=value      override one setting, e.g. -s chain.new_floor.prof_S=0.006 (repeatable)
#   -i <input_dir>    use a different input directory
#   -d <dest_dir>     use this EXACT output directory (no timestamp appended;
#                     you are being explicit, so you own the name)
#   -n <workers>      parallel worker processes
#   -S                skip step11 spatial scaling (run.spatial_scaling).
#   -w [monthly|10day] spatial scaling window mode (spatial_scaling.window).
#   -p                run diagnostics/report_spatial_scaling.py after step11 (run.scaling_plots).
#   -c                run step12_profiles_compact.py (run.compact).
#   -A                skip step13 ECCO-ready file packing (run.pack).
#
# OUTPUT LAYOUT (inside the timestamped output dir):
#   yearly/        one file per source and year (NCEI.py, then steps 11 and 12)
#   ecco_arrays/   step13: the model-ready profile files only (see optimize_profile_files.py)
#   documentation/ step13 summary: manifest, NOBSGLOB / NFILESPROFMAX trade-off,
#                  PROFILES_SIZE.h values and data.profiles file list
#
# Examples:
#   ./run_ncei.sh                                   # defaults + timestamp, scaling on, no compact
#   ./run_ncei.sh -i /path/to/other_inputs          # other input, default+ts output
#   ./run_ncei.sh -d /path/to/exact_output_dir      # default input, exact output
#   ./run_ncei.sh -S                                # skip spatial scaling
#   ./run_ncei.sh -w 10day                          # use 10-day windows instead of monthly
#   ./run_ncei.sh -c                                # also run compact step
#   ./run_ncei.sh -n 8 -c                           # 8 workers + compact
#   ./run_ncei.sh -C my_run.yaml                    # settings for this run from a file
#   ./run_ncei.sh -s chain.new_floor.prof_S=0.006   # change one setting for this run
#
# The settings used are saved in <output>/documentation/run_config.yaml.
#
# The exit code of NCEI.py is preserved (via PIPESTATUS).

set -euo pipefail

# --- locate ourselves so paths work regardless of CWD ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/z_logs"
mkdir -p "${LOG_DIR}"

STAMP="$(date +%Y%m%d_%H%M%S)"

# --- parse optional overrides ---
CONFIG_ARGS=()               # passed to every Python step: --config FILE / --set KEY=VALUE
INPUT_DIR=""
DEST_DIR=""
N_WORKERS=""
OPT_SCALING=""; OPT_WINDOW=""; OPT_PLOTS=""; OPT_COMPACT=""; OPT_PACK=""
while getopts "C:s:i:d:n:w:SpcAh" opt; do
  case "${opt}" in
    C) CONFIG_ARGS+=(--config "${OPTARG}") ;;
    s) CONFIG_ARGS+=(--set "${OPTARG}") ;;
    i) INPUT_DIR="${OPTARG}" ;;
    d) DEST_DIR="${OPTARG}" ;;
    n) N_WORKERS="${OPTARG}" ;;
    w) OPT_WINDOW="${OPTARG}" ;;
    S) OPT_SCALING=0 ;;
    p) OPT_PLOTS=1 ;;
    c) OPT_COMPACT=1 ;;
    A) OPT_PACK=0 ;;
    h) grep '^#' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Usage: $0 [-C settings.yaml] [-s key=value] [-i input_dir] [-d dest_dir] [-n n_workers] [-w monthly|10day] [-S] [-p] [-c] [-A]"; exit 2 ;;
  esac
done

# --- read the settings (one value per line; true/false as 1/0, null as empty) ---
SETTINGS="$(python3 "${SCRIPT_DIR}/config_loader.py" ${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"} --get \
  run.input_dir run.output_base run.n_workers run.spatial_scaling spatial_scaling.window \
  run.scaling_plots run.compact run.pack)" || { echo "ERROR: could not read the settings (see above)" >&2; exit 2; }
{ read -r CFG_INPUT; read -r CFG_OUTPUT_BASE; read -r CFG_WORKERS; read -r CFG_SCALING; read -r CFG_WINDOW
  read -r CFG_PLOTS; read -r CFG_COMPACT; read -r CFG_PACK; } <<< "${SETTINGS}"
RUN_SCALING="${OPT_SCALING:-${CFG_SCALING}}"
SCALING_WINDOW="${OPT_WINDOW:-${CFG_WINDOW}}"
RUN_PLOTS="${OPT_PLOTS:-${CFG_PLOTS}}"
RUN_COMPACT="${OPT_COMPACT:-${CFG_COMPACT}}"
RUN_PACK="${OPT_PACK:-${CFG_PACK}}"
N_WORKERS="${N_WORKERS:-${CFG_WORKERS}}"

# --- apply defaults where not overridden ---
if [[ -z "${INPUT_DIR}" ]]; then
  INPUT_DIR="${CFG_INPUT}"
fi
if [[ -z "${INPUT_DIR}" ]]; then
  echo "ERROR: no input dir: give -i, or set run.input_dir in config/local.yaml" >&2
  exit 2
fi
if [[ -z "${DEST_DIR}" ]]; then
  if [[ -z "${CFG_OUTPUT_BASE}" ]]; then
    echo "ERROR: no output dir: give -d, or set run.output_base in config/local.yaml" >&2
    exit 2
  fi
  # no -d given: default base + timestamp so nothing is ever overwritten
  DEST_DIR="${CFG_OUTPUT_BASE}_${STAMP}"
fi

LOG_FILE="${LOG_DIR}/ncei_log_${STAMP}.log"
YEARLY_DIR="${DEST_DIR}/yearly"
ARRAYS_DIR="${DEST_DIR}/ecco_arrays"

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
python3 -u "${SCRIPT_DIR}/NCEI.py" -i "${INPUT_DIR}" -d "${YEARLY_DIR}" ${N_WORKERS_ARG} \
  ${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"} 2>&1 | tee -a "${LOG_FILE}"
STATUS="${PIPESTATUS[0]}"

{
echo ""
echo "==============================================================="
echo "NCEI chain run finished: $(date) (exit code ${STATUS})"
echo "Output dir: ${DEST_DIR}"
echo "==============================================================="
} | tee -a "${LOG_FILE}"

if [[ "${RUN_SCALING}" -eq 1 && "${STATUS}" -eq 0 ]]; then
  {
  echo ""
  echo "==============================================================="
  echo "Step 11 (spatial scaling): window=${SCALING_WINDOW}"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"

  python3 -u "${SCRIPT_DIR}/step11_spatial_scaling.py" \
    "${YEARLY_DIR}" --window "${SCALING_WINDOW}" ${N_WORKERS_ARG} \
    ${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"} 2>&1 | tee -a "${LOG_FILE}"
  STATUS="${PIPESTATUS[0]}"

  {
  echo ""
  echo "==============================================================="
  echo "Step 11 (spatial scaling) finished: $(date) (exit code ${STATUS})"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"
fi

if [[ "${RUN_PLOTS}" -eq 1 && "${STATUS}" -eq 0 ]]; then
  {
  echo ""
  echo "==============================================================="
  echo "Spatial scaling diagnostic plots"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"

  python3 -u "${SCRIPT_DIR}/diagnostics/report_spatial_scaling.py" \
    "${YEARLY_DIR}" --window "${SCALING_WINDOW}" 2>&1 | tee -a "${LOG_FILE}"

  {
  echo ""
  echo "==============================================================="
  echo "Diagnostic plots finished: $(date)"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"
fi

if [[ "${RUN_COMPACT}" -eq 1 && "${STATUS}" -eq 0 ]]; then
  {
  echo ""
  echo "==============================================================="
  echo "Step 12 (compact): rewriting output files in-place"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"

  COMPACT_FAILED=0
  while IFS= read -r -d '' NC_FILE; do
    STEM="${NC_FILE%.nc}"
    TMP_OUT="${STEM}_compact.nc"
    python3 -u "${SCRIPT_DIR}/step12_profiles_compact.py" \
      "${NC_FILE}" -o "${TMP_OUT}" --complevel 2 --overwrite 2>&1 | tee -a "${LOG_FILE}"
    COMPACT_STATUS="${PIPESTATUS[0]}"
    if [[ "${COMPACT_STATUS}" -eq 0 ]]; then
      # Rename to reflect step_12 and remove the old file
      NEW_NC_FILE="$(echo "${NC_FILE}" | sed 's/__ncei_step_[0-9][0-9]*/__ncei_step_12/')"
      mv "${TMP_OUT}" "${NEW_NC_FILE}"
      if [[ "${NEW_NC_FILE}" != "${NC_FILE}" ]]; then
        rm "${NC_FILE}"
      fi
    else
      echo "WARNING: step12 failed for ${NC_FILE}; leaving original in place" | tee -a "${LOG_FILE}"
      COMPACT_FAILED=1
    fi
  done < <(find "${YEARLY_DIR}" -name '*.nc' -print0)

  {
  echo ""
  echo "==============================================================="
  echo "Step 12 (compact) finished: $(date)"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"

  if [[ "${COMPACT_FAILED}" -ne 0 ]]; then
    STATUS=1
  fi
fi

if [[ "${RUN_PACK}" -eq 1 && "${STATUS}" -eq 0 ]]; then
  {
  echo ""
  echo "==============================================================="
  echo "Step 13 (ECCO-ready files): ${YEARLY_DIR} -> ${ARRAYS_DIR}"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"

  python3 -u "${SCRIPT_DIR}/step13_pack_ecco_arrays.py" \
    "${YEARLY_DIR}" "${ARRAYS_DIR}" ${N_WORKERS_ARG} ${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"} 2>&1 | tee -a "${LOG_FILE}"
  STATUS="${PIPESTATUS[0]}"

  {
  echo ""
  echo "==============================================================="
  echo "Step 13 finished: $(date) (exit code ${STATUS})"
  echo "==============================================================="
  } | tee -a "${LOG_FILE}"
fi

exit "${STATUS}"
