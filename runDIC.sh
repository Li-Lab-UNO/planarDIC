#!/usr/bin/env bash
# Run the strain-tracking pipeline for one or more tests, from sorted (or unsorted) CellScale images to
# tracks, figures and exports. Run from anywhere; paths are relative to this script's folder.
#
#   ./runDIC.sh N1739                       all N1739 tests, every step
#   ./runDIC.sh N1739/test1 N1739/test2     two tests
#   ./runDIC.sh --steps track,jfig N1757    only some steps
#   ./runDIC.sh --dry-run N1739             print the commands, run nothing
#   ./runDIC.sh --setup                     only build .venv from requirements.txt
#
# Run ./runDIC.sh --help for all options.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

ALL_STEPS="sort config scale widths thickness track artifacts jfig exports edge compare"
STEPS="$ALL_STEPS"
DEFS="Def11 Def12 Def13 Def14 Def15 Def16 Def17 Def18 Def19"
EXPORT_DEFS="Def11"
PY="$HERE/.venv/bin/python"
LABEL=""
THICKNESS=auto
DRY=0
YES=0
SETUP=0
if command -v nproc >/dev/null 2>&1; then NCPU=$(nproc); else NCPU=$(sysctl -n hw.ncpu 2>/dev/null || echo 4); fi
JOBS=$(( NCPU > 1 ? NCPU - 1 : 1 ))
STATUS_EVERY=1
LOADING_ONLY=""

usage() {
  cat <<EOF
Usage: ./runDIC.sh [options] [TEST...]

Defaults come from this computer's settings file, runDIC.mac.ini on macOS and runDIC.ubuntu.ini on
Linux (both tracked by git); options given here override it for one run. Without TEST, the tests listed
in the settings file are run.

TEST is SAMPLE/testN (e.g. N1739/test1) or a whole sample (e.g. N1739 = all its tests). Each part is
any unique prefix of the folder name, as for --test in the scripts.

Steps (in this order; --steps picks a subset, comma separated):
  sort     sort_images.py       sort the images into Def01-Def19 folders (safe to re-run)
  config   make_config.py       find and refine the tine tips (--update keeps every other setting;
                                a test without test_config.json gets a new one, see --label)
  scale    calibrate_scale.py   image scale from Def01 and Def19 (tine pairs in calibration_exclude
                                of test_config.json are left out)
  widths   measure_tines.py     tine spans -> widths_mm in test_config.json
  thickness measure_thickness.py thickness from the photos in <SAMPLE>/Thickness*/test N -> thickness_mm,
                                in place of the workbook value (thickness = auto), or only for tests
                                without one (thickness = workbook / --thickness workbook)
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
  -d, --defs LIST          deformations to track, e.g. Def11,Def12 or Def11-Def15 (default: Def11-Def19)
      --loading-only       track each deformation to the end of loading only (about 55% of the frames;
                           the loading frames are the same as in a full-cycle track)
      --full-cycle         track the full load-unload cycle (overrides loading_only = yes)
      --status N           while tracking, print the progress of the running deformations every N s
                           (default: 1; 0 turns it off)
  -e, --export-defs LIST   deformations to export (default: Def11)
  -l, --label TEXT         label for a new config, e.g. "labeled healthy"
      --thickness SOURCE   auto (default): thickness measured from the photos for every test (a test whose
                           photos give none keeps the workbook value); workbook: the sample's workbook,
                           measured from the photos only where there is none
  -p, --python PATH        Python to use (default: .venv/bin/python)
      --setup              build .venv from requirements.txt first (needs python3 and python3-venv)
  -n, --dry-run            print the commands only
  -y, --yes                don't ask for confirmation
      --ask                ask for confirmation (overrides yes = yes)
      --config FILE        read the settings from FILE instead of runDIC.mac.ini / runDIC.ubuntu.ini
  -h, --help               this help

A failed sort, config, scale, widths or tracking step skips that test's later steps; a failed figure or
export doesn't. Logs: <test folder>/DIC_logs/run_<date>_<time>/, one file per step plus runDIC.log
(everything printed here); edge and compare log to comparisons/DIC_logs/run_<date>_<time>/.
EOF
}

# ---- settings file: built-in defaults < runDIC.mac.ini / runDIC.ubuntu.ini < command line ----------------------------------------
# expand_defs LIST: "Def11-Def19" -> "Def11 Def12 ... Def19"; commas -> spaces
expand_defs() {
  local out="" d
  for d in ${1//,/ }; do
    if [[ "$d" =~ ^[Dd]ef([0-9]+)-[Dd]ef([0-9]+)$ ]]; then
      for ((i = 10#${BASH_REMATCH[1]}; i <= 10#${BASH_REMATCH[2]}; i++)); do out+=" $(printf 'Def%02d' "$i")"; done
    else out+=" $d"; fi
  done
  echo "${out# }"
}
yesno() {   # yesno FILE LINE KEY VALUE -> 1 / 0
  case "$(echo "$4" | tr '[:upper:]' '[:lower:]')" in
    yes|y|true|on|1) echo 1 ;; no|n|false|off|0) echo 0 ;;
    *) echo "$1 line $2: $3 must be yes or no, not '$4'" >&2; exit 2 ;;
  esac
}
CONF_TESTS=()
read_config() {
  local f="$1" line key val n=0
  while IFS= read -r line || [ -n "$line" ]; do
    n=$((n + 1))
    [[ "$line" =~ ^[[:space:]]*\; ]] && continue                # ; comment line (INI style)
    [[ "$line" =~ ^[[:space:]]*\[.*\][[:space:]]*$ ]] && continue  # [section] headers are ignored
    line="${line%%#*}"
    [[ "$line" =~ ^[[:space:]]*$ ]] && continue
    [[ "$line" == *=* ]] || { echo "$f line $n: expected 'key = value'"; exit 2; }
    key=$(echo "${line%%=*}" | tr -d '[:space:]')
    val=$(echo "${line#*=}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    [ -z "$val" ] && continue                                  # empty: keep the built-in default
    case "$key" in
      tests) read -r -a CONF_TESTS <<< "${val//,/ }" ;;
      steps) STEPS="${val//,/ }" ;;
      jobs) if [ "$val" != auto ]; then
              [[ "$val" =~ ^[0-9]+$ ]] && [ "$val" -gt 0 ] || { echo "$f line $n: jobs must be auto or a number"; exit 2; }
              JOBS="$val"
            fi ;;
      defs) DEFS="$(expand_defs "$val")" ;;
      export_defs) EXPORT_DEFS="$(expand_defs "$val")" ;;
      status) [[ "$val" =~ ^[0-9]+$ ]] || { echo "$f line $n: status must be a number of seconds"; exit 2; }
              STATUS_EVERY="$val" ;;
      loading_only) v=$(yesno "$f" "$n" "$key" "$val") || exit 2
                    if [ "$v" = 1 ]; then LOADING_ONLY="--loading-only"; else LOADING_ONLY=""; fi ;;
      yes) YES="$(yesno "$f" "$n" "$key" "$val")" || exit 2 ;;
      label) LABEL="$val" ;;
      thickness) case "$val" in workbook|auto) THICKNESS="$val" ;;
                   *) echo "$f line $n: thickness must be workbook or auto"; exit 2 ;; esac ;;
      python) case "$val" in /*) PY="$val" ;; *) PY="$HERE/$val" ;; esac ;;
      *) echo "$f line $n: unknown setting '$key' (see runDIC.mac.ini or runDIC.ubuntu.ini)"; exit 2 ;;
    esac
  done < "$f"
}
case "$(uname -s)" in                         # this computer's settings file
  Darwin) CONF="$HERE/runDIC.mac.ini" ;;
  Linux) CONF="$HERE/runDIC.ubuntu.ini" ;;
  *) CONF="" ;;
esac
for ((i = 1; i <= $#; i++)); do              # --config FILE is read before the other options
  if [ "${!i}" = "--config" ]; then j=$((i + 1)); CONF="${!j}"; fi
done
[ -n "$CONF" ] || { echo "no settings file for $(uname -s); give one with --config FILE"; exit 2; }
if [ -e "$CONF" ]; then read_config "$CONF"; else echo "settings file $CONF not found"; exit 2; fi

TESTS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --config) shift 2 ;;
    -s|--steps) STEPS="${2//,/ }"; shift 2 ;;
    -j|--jobs) JOBS="$2"; shift 2 ;;
    -d|--defs) DEFS="$(expand_defs "$2")"; shift 2 ;;
    --status) STATUS_EVERY="$2"; shift 2 ;;
    --loading-only) LOADING_ONLY="--loading-only"; shift ;;
    --full-cycle) LOADING_ONLY=""; shift ;;
    -e|--export-defs) EXPORT_DEFS="$(expand_defs "$2")"; shift 2 ;;
    -l|--label) LABEL="$2"; shift 2 ;;
    --thickness) case "$2" in workbook|auto) THICKNESS="$2" ;; *) echo "--thickness must be workbook or auto"; exit 2 ;; esac; shift 2 ;;
    -p|--python) PY="$2"; shift 2 ;;
    --setup) SETUP=1; shift ;;
    -n|--dry-run) DRY=1; shift ;;
    -y|--yes) YES=1; shift ;;
    --ask) YES=0; shift ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "unknown option $1 (see --help)"; exit 2 ;;
    *) TESTS+=("$1"); shift ;;
  esac
done
[ ${#TESTS[@]} -eq 0 ] && TESTS=("${CONF_TESTS[@]+"${CONF_TESTS[@]}"}")   # no tests given: the settings file's

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

[ ${#TESTS[@]} -eq 0 ] && { echo "no tests: give them on the command line or set 'tests' in ${CONF#$HERE/} (--help for all options)"; exit 2; }
if [ "$DRY" = 0 ]; then
  [ -x "$PY" ] || { echo "$PY not found: run ./runDIC.sh --setup first"; exit 1; }
  "$PY" - <<'EOF' || exit 1
import sys
missing = []
for m in ("numpy", "scipy", "pandas", "matplotlib", "cv2", "imageio", "imageio_ffmpeg", "skimage"):
    try:
        __import__(m)
    except ImportError:
        missing.append(m)
if missing:
    sys.exit(f"missing Python packages: {', '.join(missing)}; run ./runDIC.sh --setup")
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
  if [ -z "$1" ]; then echo "$HERE/comparisons/DIC_logs/run_$STAMP"; return; fi
  local d; d=$(ls -d "$HERE"/"${1%%/*}"*/"${1##*/}"_*/ 2>/dev/null | head -1)
  echo "${d%/}/DIC_logs/run_$STAMP"
}
ndefs=$(echo $DEFS | wc -w | tr -d ' ')
echo "Tests:   ${TESTS[*]}"
echo "Config:  ${CONF#$HERE/}"
echo "Steps:   $STEPS"
has thickness && echo "Thickness: $([ "$THICKNESS" = auto ] && echo "measured from the photos (in place of workbook values)" || echo "workbook; measured from the photos where there is none")"
has track && echo "Track:   $DEFS ($(( ${#TESTS[@]} * ndefs )) deformations, $JOBS at a time, $(( NCPU / JOBS > 1 ? NCPU / JOBS : 1 )) threads each${LOADING_ONLY:+, loading branch only})"
has exports && echo "Exports: $EXPORT_DEFS"
echo "Python:  $PY"
[ "$DRY" = 1 ] && echo "(dry run: commands are printed, not run)"
if [ "$DRY" = 0 ]; then
  echo "Logs:    <test folder>/DIC_logs/run_$STAMP/"
  if [ "$YES" = 0 ] && [ -t 0 ]; then
    read -r -p "Proceed? [y/N] " ans
    case "$ans" in y|Y|yes) ;; *) echo "stopped"; exit 0 ;; esac
  fi
  STATE=$(mktemp -d)                    # failures list for the summary; removed at the end
  CONSOLE=()
  for t in "${TESTS[@]}"; do mkdir -p "$(tlogdir "$t")"; CONSOLE+=("$(tlogdir "$t")/runDIC.log"); done
  exec > >(tee "${CONSOLE[@]}") 2>&1   # everything printed from here on also goes to each test's runDIC.log
  T_RUN=$SECONDS                         # runtime is counted from here (after the confirmation question)
fi
# each tracking job writes a progress line to its log this often, so the status line can be as current
DIC_PROGRESS_EVERY=$(( STATUS_EVERY > 0 && STATUS_EVERY < 10 ? STATUS_EVERY : 10 ))
# threads per job: the CPUs shared among the parallel jobs (dic2d.THREADS; each tracking job matches chunks
# of its subsets in parallel, and OpenCV is limited to match)
DIC_THREADS=$(( NCPU / JOBS > 1 ? NCPU / JOBS : 1 ))
export DIC_PROGRESS=lines DIC_PROGRESS_EVERY DIC_THREADS PY DRY HERE STAMP STATE LOADING_ONLY
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
# hms SECONDS: "1 h 02 min 05 s", "12 min 05 s" or "5 s"
hms() {
  local s=$1
  if [ "$s" -ge 3600 ]; then printf "%d h %02d min %02d s" $((s / 3600)) $((s % 3600 / 60)) $((s % 60))
  elif [ "$s" -ge 60 ]; then printf "%d min %02d s" $((s / 60)) $((s % 60))
  else printf "%d s" "$s"; fi
}

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
# track_of TEST DEF: path of the deformation's full-domain track, <test>/DIC_results/<deformation folder>.npz
track_of() {
  local s=${1%%/*} n=${1##*/} f
  for f in "$HERE"/"$s"*/"$n"_*/Images/"$2"_*/; do
    f=${f%/}; f="${f%/Images/*}/DIC_results/${f##*/}.npz"
    [ -e "$f" ] && { echo "$f"; return 0; }
  done
  return 1
}
export -f track_of
# all_tracked TEST: true when Def11-Def19 all have a full-domain track (needed by jfig, edge and compare)
all_tracked() {
  [ "$DRY" = 1 ] && return 0
  local d
  for d in Def11 Def12 Def13 Def14 Def15 Def16 Def17 Def18 Def19; do
    track_of "$1" "$d" >/dev/null || return 1
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

if has thickness; then
  step "thickness from the photos ($([ "$THICKNESS" = auto ] && echo "in place of workbook values" || echo "where there is no workbook value"))"
  for t in "${TESTS[@]}"; do
    n=${t##*/test}; photos=""
    for d in "$HERE"/"${t%%/*}"*/Thickness*/*/; do      # <SAMPLE>/Thickness, <ID>/test <N> or test<N>
      b=$(basename "$d"); [[ "$b" =~ ^[Tt]est[[:space:]]*0*$n$ ]] && photos="${d%/}"
    done
    if [ -z "$photos" ]; then echo "  skip  thickness $(key "$t"): no thickness photo folder"; continue; fi
    SOFT=1 run thickness "$t" "$PY" "$HERE/thickness/measure_thickness.py" "$photos" --write \
      $([ "$THICKNESS" = workbook ] && echo --keep-workbook)
  done
fi

if has track; then
  step "tracking ($JOBS jobs at a time; one log per deformation)"
  track_one() {   # TEST DEF
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      run "track_$2" "$1" "$PY" "$HERE/strain_tracking/track_protocols.py" --test "$1" --defs "$2" $LOADING_ONLY
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
    echo "  tracking finished in $(hms $(( SECONDS - t_track )))"
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
    track_of "$t" Def11 >/dev/null || { echo "  skip  artifacts $tk: no Def11 track"; continue; }
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
      | xargs -P "$JOBS" -L 1 bash -c 'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 run "retrack_$1" "$0" "$PY" "$HERE/strain_tracking/track_protocols.py" --test "$0" --defs "$1" $LOADING_ONLY'
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
    local t="$1" d="$2" npz
    npz=$(track_of "$t" "$d")
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
echo "logs: <test folder>/DIC_logs/run_$STAMP/"
if [ -s "$STATE/failures.txt" ]; then
  echo "FAILED steps:"; sed 's/^/  /' "$STATE/failures.txt"; rm -rf "$STATE"
  echo "total runtime: $(hms $(( SECONDS - T_RUN )))  (finished $(date '+%Y-%m-%d %H:%M:%S'))"; exit 1
fi
rm -rf "$STATE"
echo "all steps finished"
echo "total runtime: $(hms $(( SECONDS - T_RUN )))  (finished $(date '+%Y-%m-%d %H:%M:%S'))"
