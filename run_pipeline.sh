#!/usr/bin/env bash
# Run the strain-tracking pipeline for one or more tests, from sorted (or unsorted) CellScale images to
# tracks, figures and exports. Run from anywhere; paths are relative to this script's folder.
#
#   ./run_pipeline.sh N1739                       all N1739 tests, every step
#   ./run_pipeline.sh N1739/test1 N1739/test2     two tests
#   ./run_pipeline.sh --steps track,jfig N1757    only some steps
#   ./run_pipeline.sh --dry-run N1739             print the commands, run nothing
#   ./run_pipeline.sh --setup                     only build .venv from requirements.txt
#
# Run ./run_pipeline.sh --help for all options.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

ALL_STEPS="sort config scale widths track artifacts jfig exports edge compare"
STEPS="$ALL_STEPS"
DEFS="Def11 Def12 Def13 Def14 Def15 Def16 Def17 Def18 Def19"
EXPORT_DEFS="Def11"
PY="$HERE/.venv/bin/python"
LABEL=""
DRY=0
YES=0
SETUP=0
if command -v nproc >/dev/null 2>&1; then NCPU=$(nproc); else NCPU=$(sysctl -n hw.ncpu 2>/dev/null || echo 4); fi
JOBS=$(( NCPU > 1 ? NCPU - 1 : 1 ))
STATUS_EVERY=30

usage() {
  cat <<EOF
Usage: ./run_pipeline.sh [options] TEST...

TEST is SAMPLE/testN (e.g. N1739/test1) or a whole sample (e.g. N1739 = all its tests). Each part is
any unique prefix of the folder name, as for --test in the scripts.

Steps (in this order; --steps picks a subset, comma separated):
  sort     sort_images.py       sort the images into Def01-Def19 folders (safe to re-run)
  config   make_config.py       find and refine the tine tips (--update keeps every other setting;
                                a test without test_config.json gets a new one, see --label)
  scale    calibrate_scale.py   image scale from Def01 and Def19 (tine pairs in calibration_exclude
                                of test_config.json are left out)
  widths   measure_tines.py     tine spans -> widths_mm in test_config.json
  track    track_protocols.py   Def11-Def19 against the common reference, one job per deformation
  artifacts artifact_check.py   find local tracking errors (sliding particles, glare) at the Def11 peak; when
                                one is confirmed, exclude it in test_config.json and track the test again
  jfig     peak_j_figure.py     peak J/mean figure per test (needs all of Def11-Def19 tracked)
  exports  compare_cellscale.py, export_fields.py, animate_fields.py for the --export-defs
  edge     edge_effects.py      tine-effect check across the tests given (all of Def11-Def19 tracked)
  compare  compare_peak_j.py    cross-test J/mean comparison for the samples given (same)

Options:
  -s, --steps LIST         steps to run (default: all)
  -j, --jobs N             parallel tracking / export jobs (default: CPU count - 1 = $JOBS)
  -d, --defs LIST          deformations to track (default: Def11-Def19)
      --status N           while tracking, print the progress of the running deformations every N s
                           (default: 30; 0 turns it off)
  -e, --export-defs LIST   deformations to export (default: Def11)
  -l, --label TEXT         label for a new config, e.g. "labeled healthy"
  -p, --python PATH        Python to use (default: .venv/bin/python)
      --setup              build .venv from requirements.txt first (needs python3 and python3-venv)
  -n, --dry-run            print the commands only
  -y, --yes                don't ask for confirmation
  -h, --help               this help

A failed sort, config, scale, widths or tracking step skips that test's later steps; a failed figure or
export doesn't. Logs: <test folder>/pipeline_logs/run_<date>_<time>/, one file per step plus pipeline.log
(everything printed here); edge and compare log to comparisons/pipeline_logs/run_<date>_<time>/.
EOF
}

TESTS=()
while [ $# -gt 0 ]; do
  case "$1" in
    -s|--steps) STEPS="${2//,/ }"; shift 2 ;;
    -j|--jobs) JOBS="$2"; shift 2 ;;
    -d|--defs) DEFS="${2//,/ }"; shift 2 ;;
    --status) STATUS_EVERY="$2"; shift 2 ;;
    -e|--export-defs) EXPORT_DEFS="${2//,/ }"; shift 2 ;;
    -l|--label) LABEL="$2"; shift 2 ;;
    -p|--python) PY="$2"; shift 2 ;;
    --setup) SETUP=1; shift ;;
    -n|--dry-run) DRY=1; shift ;;
    -y|--yes) YES=1; shift ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "unknown option $1 (see --help)"; exit 2 ;;
    *) TESTS+=("$1"); shift ;;
  esac
done

for s in $STEPS; do
  case " $ALL_STEPS " in *" $s "*) ;; *) echo "unknown step '$s'; steps are: $ALL_STEPS"; exit 2 ;; esac
done

# ---- environment ---------------------------------------------------------------------------------------
if [ "$SETUP" = 1 ]; then
  if [ -x "$HERE/.venv/bin/python" ]; then
    echo ".venv exists; installing/updating requirements"
  else
    command -v python3 >/dev/null || { echo "python3 not found (Ubuntu: sudo apt install python3 python3-venv)"; exit 1; }
    python3 -m venv "$HERE/.venv" || { echo "could not create .venv (Ubuntu: sudo apt install python3-venv)"; exit 1; }
  fi
  "$HERE/.venv/bin/python" -m pip install --upgrade pip >/dev/null
  "$HERE/.venv/bin/python" -m pip install -r "$HERE/requirements.txt" || exit 1
  [ ${#TESTS[@]} -eq 0 ] && { echo "environment ready"; exit 0; }
fi

[ ${#TESTS[@]} -eq 0 ] && { usage; exit 2; }
if [ "$DRY" = 0 ]; then
  [ -x "$PY" ] || { echo "$PY not found: run ./run_pipeline.sh --setup first"; exit 1; }
  "$PY" - <<'EOF' || exit 1
import sys
missing = []
for m in ("numpy", "scipy", "pandas", "matplotlib", "cv2", "imageio", "imageio_ffmpeg", "skimage"):
    try:
        __import__(m)
    except ImportError:
        missing.append(m)
if missing:
    sys.exit(f"missing Python packages: {', '.join(missing)}; run ./run_pipeline.sh --setup")
if not (3, 9) <= sys.version_info[:2] <= (3, 12):
    print(f"note: tested with Python 3.9-3.12, this is {sys.version.split()[0]}")
EOF
fi

# ---- tests -----------------------------------------------------------------------------------------------
# Expand a bare sample (N1739) to its test folders; keep SAMPLE/testN as given.
EXPANDED=()
for t in "${TESTS[@]}"; do
  if [[ "$t" == */* ]]; then
    ls -d "$HERE"/"${t%%/*}"*/"${t##*/}"_*/ >/dev/null 2>&1 || { echo "no test folder found for '$t'"; exit 1; }
    EXPANDED+=("$t")
  else
    found=0
    for d in "$HERE"/"$t"*/test[0-9]*_*/; do
      [ -d "$d" ] || continue
      s=$(basename "$(dirname "$d")"); n=$(basename "$d"); n=${n%%_*}
      EXPANDED+=("${s%%_*}/$n"); found=1
    done
    [ "$found" = 1 ] || { echo "no test folders found for sample '$t'"; exit 1; }
  fi
done
TESTS=("${EXPANDED[@]}")
SAMPLES=$(for t in "${TESTS[@]}"; do echo "${t%%/*}"; done | awk '!seen[$0]++' | tr '\n' ' ')

has() { case " $STEPS " in *" $1 "*) return 0 ;; esac; return 1; }
key() { echo "${1//\//_}"; }

STAMP=$(date +%Y%m%d_%H%M%S)
# tlogdir TEST: this run's log folder for TEST, or for the cross-test steps when TEST is empty
tlogdir() {
  if [ -z "$1" ]; then echo "$HERE/comparisons/pipeline_logs/run_$STAMP"; return; fi
  local d; d=$(ls -d "$HERE"/"${1%%/*}"*/"${1##*/}"_*/ 2>/dev/null | head -1)
  echo "${d%/}/pipeline_logs/run_$STAMP"
}
ndefs=$(echo $DEFS | wc -w | tr -d ' ')
echo "Tests:   ${TESTS[*]}"
echo "Steps:   $STEPS"
has track && echo "Track:   $DEFS ($(( ${#TESTS[@]} * ndefs )) deformations, $JOBS at a time)"
has exports && echo "Exports: $EXPORT_DEFS"
echo "Python:  $PY"
[ "$DRY" = 1 ] && echo "(dry run: commands are printed, not run)"
if [ "$DRY" = 0 ]; then
  echo "Logs:    <test folder>/pipeline_logs/run_$STAMP/"
  if [ "$YES" = 0 ] && [ -t 0 ]; then
    read -r -p "Proceed? [y/N] " ans
    case "$ans" in y|Y|yes) ;; *) echo "stopped"; exit 0 ;; esac
  fi
  STATE=$(mktemp -d)                    # failures list for the summary; removed at the end
  CONSOLE=()
  for t in "${TESTS[@]}"; do mkdir -p "$(tlogdir "$t")"; CONSOLE+=("$(tlogdir "$t")/pipeline.log"); done
  exec > >(tee "${CONSOLE[@]}") 2>&1   # everything printed from here on also goes to each test's pipeline.log
fi
export DIC_PROGRESS=lines DIC_PROGRESS_EVERY=10 PY DRY HERE STAMP STATE
export -f tlogdir

# run NAME TEST CMD...: run one command with its own log in TEST's log folder (TEST empty: a cross-test step).
# A failure marks the test as failed, so its later steps are skipped, unless SOFT=1 (figures and exports:
# nothing else depends on them).
run() {
  local name="$1" t="$2"; shift 2
  if [ "$DRY" = 1 ]; then echo "  $*"; return 0; fi
  local dir tk="${t//\//_}"; dir=$(tlogdir "$t"); mkdir -p "$dir"
  if [ -n "$t" ] && [ -e "$dir/FAILED" ]; then echo "  skip $name ($tk failed earlier)"; return 1; fi
  local log="$dir/$name.log" t0=$SECONDS
  if "$@" >"$log" 2>&1; then
    echo "  ok    $name ${tk:+$tk }($(( SECONDS - t0 )) s)"
  else
    echo "  FAIL  $name ${tk:+$tk }($(( SECONDS - t0 )) s): last lines of $log"
    tail -n 5 "$log" | sed 's/^/        /'
    [ -n "$t" ] && [ "${SOFT:-0}" != 1 ] && touch "$dir/FAILED"
    echo "$name ${tk:+$tk }($log)" >> "$STATE/failures.txt"
    return 1
  fi
}
export -f run

ST="$HERE/strain_tracking"
step() { echo; echo "== $1  ($(date +%H:%M:%S))"; }

# track_status SECONDS: while tracking runs, print one line every SECONDS with the deformations finished
# and the progress of each running one (read from the progress lines in its track_*.log)
track_status() {
  local every="$1" total=$(( ${#TESTS[@]} * ndefs )) t0=$SECONDS
  [ "$every" -gt 0 ] 2>/dev/null || return 0
  while sleep "$every"; do
    local done=0 running="" f last tk
    for t in "${TESTS[@]}"; do
      tk=""; [ "${#TESTS[@]}" -gt 1 ] && tk="$(key "$t") "
      for f in "$(tlogdir "$t")"/track_*.log; do
        [ -e "$f" ] || continue
        if grep -q "deformations tracked in\|Traceback" "$f"; then done=$(( done + 1 )); continue; fi
        last=$(grep -E "^  Def[0-9]+ +[0-9]+/[0-9]+ +[0-9]+%" "$f" | tail -1 | awk '{print $1, $3}')
        running+="${running:+, }$tk${last:-$(basename "$f" .log | sed 's/track_//') starting}"
      done
    done
    echo "  [$(date +%H:%M:%S), $(( (SECONDS - t0) / 60 )) min] $done/$total done${running:+; running: $running}"
  done
}
# all_tracked TEST: true when Def11-Def19 all have a full-domain track (needed by jfig, edge and compare)
all_tracked() {
  [ "$DRY" = 1 ] && return 0
  local s=${1%%/*} n=${1##*/} d
  for d in Def11 Def12 Def13 Def14 Def15 Def16 Def17 Def18 Def19; do
    ls "$HERE"/"$s"*/"$n"_*/Images/"$d"_*/strain_results/"$d"_*_tine_region.npz >/dev/null 2>&1 || return 1
  done
}

if has sort; then
  step "sort images"
  for t in "${TESTS[@]}"; do run sort "$t" "$PY" "$ST/sort_images.py" --test "$t"; done
fi

if has config; then
  step "test_config.json (tine tips)"
  for t in "${TESTS[@]}"; do
    s=${t%%/*}; n=${t##*/}
    if ls "$HERE"/"$s"*/"$n"_*/test_config.json >/dev/null 2>&1; then
      run config "$t" "$PY" "$ST/make_config.py" --test "$t" --update
    elif [ -n "$LABEL" ]; then
      run config "$t" "$PY" "$ST/make_config.py" --test "$t" --label "$LABEL"
    else
      run config "$t" "$PY" "$ST/make_config.py" --test "$t"
    fi
  done
fi

if has scale; then
  step "scale calibration"
  for t in "${TESTS[@]}"; do run scale "$t" "$PY" "$ST/calibrate_scale.py" --test "$t"; done
fi

if has widths; then
  step "tine spans"
  for t in "${TESTS[@]}"; do run widths "$t" "$PY" "$ST/measure_tines.py" --test "$t" --write; done
fi

if has track; then
  step "tracking ($JOBS jobs at a time; one log per deformation)"
  track_one() {   # TEST DEF
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      run "track_$2" "$1" "$PY" "$HERE/strain_tracking/track_protocols.py" --test "$1" --defs "$2"
  }
  export -f track_one
  if [ "$DRY" = 1 ]; then
    for t in "${TESTS[@]}"; do for d in $DEFS; do track_one "$t" "$d"; done; done
  else
    track_status "$STATUS_EVERY" &
    status_pid=$!
    t_track=$SECONDS
    for t in "${TESTS[@]}"; do for d in $DEFS; do echo "$t $d"; done; done \
      | xargs -P "$JOBS" -L 1 bash -c 'track_one "$0" "$1"'
    kill "$status_pid" 2>/dev/null; wait "$status_pid" 2>/dev/null
    echo "  tracking finished in $(( (SECONDS - t_track) / 60 )) min $(( (SECONDS - t_track) % 60 )) s"
    for t in "${TESTS[@]}"; do      # a test with a failed deformation is not passed on to later steps
      d=$(tlogdir "$t")
      ls "$d"/track_*.log >/dev/null 2>&1 || continue
      grep -l "Traceback\|Error" "$d"/track_*.log >/dev/null 2>&1 && touch "$d/FAILED"
    done
  fi
fi

if has artifacts; then
  step "artifact check"
  RETRACK=()
  for t in "${TESTS[@]}"; do
    tk=$(key "$t")
    if [ "$DRY" = 1 ]; then echo "  $PY $ST/artifact_check.py --test $t   (exit 3: exclusions added, track again)"; continue; fi
    [ -e "$(tlogdir "$t")/FAILED" ] && { echo "  skip artifacts ($tk failed earlier)"; continue; }
    ls "$HERE"/"${t%%/*}"*/"${t##*/}"_*/Images/Def11_*/strain_results/Def11_*_tine_region.npz >/dev/null 2>&1 \
      || { echo "  skip  artifacts $tk: no Def11 track"; continue; }
    log="$(tlogdir "$t")/artifacts.log"
    "$PY" "$ST/artifact_check.py" --test "$t" > "$log" 2>&1; rc=$?
    case $rc in
      0) echo "  ok    artifacts $tk: nothing to change" ;;
      3) echo "  found artifacts $tk: $(grep -h 'added to dic.exclude_px' "$log" | sed 's/^ *//')"; RETRACK+=("$t") ;;
      *) echo "  FAIL  artifacts $tk (see $log)"; echo "artifacts $tk ($log)" >> "$STATE/failures.txt" ;;
    esac
  done
  if [ ${#RETRACK[@]} -gt 0 ]; then
    step "tracking again with the new exclusions: ${RETRACK[*]}"
    for t in "${RETRACK[@]}"; do for d in $DEFS; do echo "$t $d"; done; done \
      | xargs -P "$JOBS" -L 1 bash -c 'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 run "retrack_$1" "$0" "$PY" "$HERE/strain_tracking/track_protocols.py" --test "$0" --defs "$1"'
    for t in "${RETRACK[@]}"; do
      SOFT=1 run artifacts_recheck "$t" "$PY" "$ST/artifact_check.py" --test "$t" --report-only
      grep -h "CONFIRMED" "$(tlogdir "$t")/artifacts_recheck.log" >/dev/null 2>&1 \
        && echo "  note: $(key "$t") still has a confirmed artifact after the exclusion; see its artifact_check.png"
    done
  fi
fi

if has jfig; then
  step "peak J/mean figures"
  for t in "${TESTS[@]}"; do
    if all_tracked "$t"; then SOFT=1 run jfig "$t" "$PY" "$ST/peak_j_figure.py" --test "$t"
    else echo "  skip  jfig $(key "$t"): needs the tracks of all of Def11-Def19"; fi
  done
fi

if has exports; then
  step "exports ($EXPORT_DEFS)"
  export_one() {   # TEST DEF: CellScale comparison first, so the animations include the force traces
    local t="$1" d="$2" s=${1%%/*} n=${1##*/} npz
    npz=$(ls "$HERE"/"$s"*/"$n"_*/Images/"$d"_*/strain_results/"$d"_*_tine_region.npz 2>/dev/null | head -1)
    if [ -z "$npz" ] && [ "$DRY" = 0 ]; then echo "  FAIL  exports ${t//\//_} $d: no track"; return 1; fi
    npz=${npz:-<track of $t $d>}
    SOFT=1 run "cellscale_$d" "$t" "$PY" "$HERE/strain_tracking/compare_cellscale.py" "$npz" &&
    SOFT=1 run "fields_$d" "$t" "$PY" "$HERE/strain_tracking/export_fields.py" "$npz" &&
    SOFT=1 run "animate_$d" "$t" "$PY" "$HERE/strain_tracking/animate_fields.py" "$npz"
  }
  export -f export_one
  if [ "$DRY" = 1 ]; then
    for t in "${TESTS[@]}"; do for d in $EXPORT_DEFS; do export_one "$t" "$d"; done; done
  else
    for t in "${TESTS[@]}"; do for d in $EXPORT_DEFS; do echo "$t $d"; done; done \
      | xargs -P "$JOBS" -L 1 bash -c 'export_one "$0" "$1"'
  fi
fi

OK_TESTS=()      # tests that passed every prerequisite step and have all of Def11-Def19 tracked
for t in "${TESTS[@]}"; do
  [ "$DRY" = 0 ] && [ -e "$(tlogdir "$t")/FAILED" ] && continue
  all_tracked "$t" && OK_TESTS+=("$t")
done
if { has edge || has compare; } && [ ${#OK_TESTS[@]} -eq 0 ]; then
  echo; echo "== edge-effect check and comparison skipped: they need tests with all of Def11-Def19 tracked"
fi

if has edge && [ ${#OK_TESTS[@]} -gt 0 ]; then
  step "edge-effect check"
  SOFT=1 run edge "" "$PY" "$ST/edge_effects.py" --tests "${OK_TESTS[@]}"
fi

if has compare && [ ${#OK_TESTS[@]} -gt 0 ]; then
  step "cross-test comparison"
  OK_SAMPLES=$(for t in "${OK_TESTS[@]}"; do echo "${t%%/*}"; done | awk '!seen[$0]++' | tr '\n' ' ')
  SOFT=1 run compare "" "$PY" "$ST/compare_peak_j.py" --samples $OK_SAMPLES
fi

[ "$DRY" = 1 ] && exit 0
echo
echo "== summary  ($(date +%H:%M:%S))"
if has track || has artifacts; then
  for t in "${TESTS[@]}"; do
    for f in "$(tlogdir "$t")"/track_*.log "$(tlogdir "$t")"/retrack_*.log; do
      [ -e "$f" ] || continue
      line=$(grep -h "valid per frame" "$f" | tail -1 | sed 's/; saved.*//;s/^ *//')
      echo "  $(key "$t") $(basename "$f" .log): ${line:-no result}"
    done
  done
fi
echo "logs: <test folder>/pipeline_logs/run_$STAMP/"
if [ -s "$STATE/failures.txt" ]; then
  echo "FAILED steps:"; sed 's/^/  /' "$STATE/failures.txt"; rm -rf "$STATE"; exit 1
fi
rm -rf "$STATE"
echo "all steps finished"
