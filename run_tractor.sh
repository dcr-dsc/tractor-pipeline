#!/bin/bash
set -euo pipefail

## RUNNING TRACTOR

start_time=$(date +%s)
exec >> run_tractor_serie1.log 2>&1
echo "Running Tractor started at $(date)"

CHR_RAW=""
MAX_JOBS=4
THREADS=10

PHENOFILE=""
COVARCOLLIST=""
METHOD=""
SAMPLEIDCOL=""
PHENOCOL=""

usage() {
  echo "Usage:"
  echo "  bash run_tractor.sh -chr \"1,6,7,12,13,18\"|\"1-5,21\" [-jobs 4] [-threads 10] \\"
  echo "    -phenofile covarlist.txt \\"
  echo "    -covarcollist \"sex,age,rnd_source,pc_1,pc_2\" \\"
  echo "    -method logistic|linear \\"
  echo "    -sampleidcol id \\"
  echo "    -phenocol ad"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -chr)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -chr requires an argument"; usage; }
      CHR_RAW="$2"
      shift 2
      ;;
    -jobs)
      MAX_JOBS="${2:-}"
      shift 2
      ;;
    -threads)
      THREADS="${2:-}"
      shift 2
      ;;
    -phenofile)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -phenofile requires an argument"; usage; }
      PHENOFILE="$2"
      shift 2
      ;;
    -covarcollist)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -covarcollist requires an argument"; usage; }
      COVARCOLLIST="$2"
      shift 2
      ;;
    -method)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -method requires an argument"; usage; }
      METHOD="$2"
      shift 2
      ;;
    -sampleidcol)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -sampleidcol requires an argument"; usage; }
      SAMPLEIDCOL="$2"
      shift 2
      ;;
    -phenocol)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -phenocol requires an argument"; usage; }
      PHENOCOL="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      ;;
  esac
done

[[ -z "$CHR_RAW" ]] && { echo "ERROR: -chr is required"; usage; }
[[ -z "$PHENOFILE" ]] && { echo "ERROR: -phenofile is required"; usage; }
[[ -z "$COVARCOLLIST" ]] && { echo "ERROR: -covarcollist is required"; usage; }
[[ -z "$METHOD" ]] && { echo "ERROR: -method is required"; usage; }
[[ -z "$SAMPLEIDCOL" ]] && { echo "ERROR: -sampleidcol is required"; usage; }
[[ -z "$PHENOCOL" ]] && { echo "ERROR: -phenocol is required"; usage; }

if ! [[ "$MAX_JOBS" =~ ^[0-9]+$ ]] || (( MAX_JOBS < 1 )); then
  echo "ERROR: -jobs must be an integer >= 1"
  exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS < 1 )); then
  echo "ERROR: -threads must be an integer >= 1"
  exit 1
fi

expand_chr_list() {
  local raw="$1"
  raw="$(echo "$raw" | tr -d '[:space:]')"
  local IFS=',' parts token
  read -r -a parts <<< "$raw"

  local out=()
  for token in "${parts[@]}"; do
    [[ -z "$token" ]] && continue
    if [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      local a="${BASH_REMATCH[1]}"
      local b="${BASH_REMATCH[2]}"
      if (( a <= b )); then
        for ((c=a; c<=b; c++)); do out+=("$c"); done
      else
        for ((c=a; c>=b; c--)); do out+=("$c"); done
      fi
    elif [[ "$token" =~ ^[0-9]+$ ]]; then
      out+=("$token")
    else
      echo "ERROR: Invalid token in -chr: '$token'"
      exit 1
    fi
  done

  printf "%s\n" "${out[@]}" | awk '!seen[$0]++' | sort -n
}

mapfile -t CHRS < <(expand_chr_list "$CHR_RAW")

echo "Chromosomes: ${CHRS[*]}"
echo "Parallel jobs (-jobs): $MAX_JOBS"
echo "Threads per chromosome (-threads): $THREADS"
echo "Phenofile: $PHENOFILE"
echo "Covariates: $COVARCOLLIST"
echo "Method: $METHOD"
echo "Sample ID column: $SAMPLEIDCOL"
echo "Pheno column: $PHENOCOL"

mkdir -p geno_out

running=0
for chr in "${CHRS[@]}"; do
  (
    Rscript run_tractor.R \
      --hapdose "geno_${chr}_ph" \
      --phenofile "$PHENOFILE" \
      --covarcollist "$COVARCOLLIST" \
      --method "$METHOD" \
      --output "geno_out/geno_${chr}_out" \
      --sampleidcol "$SAMPLEIDCOL" \
      --phenocol "$PHENOCOL" \
      --nthreads "$THREADS"
  ) &

  ((running+=1))
  if (( running >= MAX_JOBS )); then
    wait
    running=0
  fi
done

wait

end_time=$(date +%s)
elapsed_time=$((end_time-start_time))
elapsed_h=$((elapsed_time/3600))
elapsed_m=$(((elapsed_time % 3600)/60))
elapsed_s=$((elapsed_time % 60))

echo "Running Tractor finished"
printf "Total execution time: %02d:%02d:%02d\n"  $elapsed_h $elapsed_m $elapsed_s
echo "=^.^="









