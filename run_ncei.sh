#!/usr/bin/env bash
#
# Controller for the NCEI post-processing chain.
#
# Runs NCEI.py and automatically tees all stdout+stderr to a timestamped
# log under z_logs/, so consecutive runs never overwrite each other and
# you don't have to rename anything by hand.
#
# Usage:
#   ./run_ncei.sh -i <input_dir> -d <dest_dir>
#
# Example:
#   ./run_ncei.sh \
#       -i /Users/brucel/ecco/yip/profile_data/Interp_Profiles \
#       -d /Users/brucel/ecco/yip/profile_files_NCEI_processed
#
# The exit code of NCEI.py is preserved (via PIPESTATUS) so this script
# can be chained or monitored.

set -euo pipefail

# --- locate ourselves so paths work regardless of CWD ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/z_logs"
mkdir -p "${LOG_DIR}"

# --- timestamped log path (seconds precision => no collisions) ---
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/ncei_log_${STAMP}.log"

# --- pass all args straight through to NCEI.py (-i / -d, etc.) ---
echo "===============================================================" | tee    "${LOG_FILE}"
echo "NCEI chain run started: $(date)"                                 | tee -a "${LOG_FILE}"
echo "Arguments: $*"                                                   | tee -a "${LOG_FILE}"
echo "Log file:  ${LOG_FILE}"                                          | tee -a "${LOG_FILE}"
echo "===============================================================" | tee -a "${LOG_FILE}"
echo ""                                                                | tee -a "${LOG_FILE}"

# python3 -u keeps output unbuffered so the log fills in real time (useful
# when tailing a run inside screen) and works on stock macOS without stdbuf.
# 2>&1 folds stderr into stdout BEFORE the pipe so warnings/tracebacks land
# in the log too.
python3 -u "${SCRIPT_DIR}/NCEI.py" "$@" 2>&1 | tee -a "${LOG_FILE}"
STATUS="${PIPESTATUS[0]}"

echo ""                                                                | tee -a "${LOG_FILE}"
echo "===============================================================" | tee -a "${LOG_FILE}"
echo "NCEI chain run finished: $(date) (exit code ${STATUS})"          | tee -a "${LOG_FILE}"
echo "===============================================================" | tee -a "${LOG_FILE}"

exit "${STATUS}"
