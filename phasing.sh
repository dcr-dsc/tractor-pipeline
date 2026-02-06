#!/bin/bash
set -euo pipefail

## PHASING STEP

start_time=$(date +%s)

exec >> phasing.log 2>&1
echo "Phasing started at $(date)"

CHR_RAW=""
THREADS=8
MAX_JOBS=4

usage() {
  echo "Usage: bash phasing.sh -chr \"1-22\"|\"1-3,5,22\" [-threads 8] [-jobs 4]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -chr)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -chr requires an argument"; usage; }
      CHR_RAW="$2"
      shift 2
      ;;
    -threads)
      THREADS="${2:-}"
      shift 2
      ;;
    -jobs)
      MAX_JOBS="${2:-}"
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

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS < 1 )); then
  echo "ERROR: -threads must be an integer >= 1"
  exit 1
fi

if ! [[ "$MAX_JOBS" =~ ^[0-9]+$ ]] || (( MAX_JOBS < 1 )); then
  echo "ERROR: -jobs must be an integer >= 1"
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

mkdir -p geno_ph

running=0
for chr in "${CHRS[@]}"; do
  (
    in_bcf="geno_chr/geno_${chr}_com.bcf"
    gmap="gmap_1kGP/chr${chr}.b38.gmap"
    ref_bcf="ref_1kGP/ref_${chr}_com.bcf"
    out_bcf="geno_ph/geno_${chr}_ph.bcf"

    shapeit4 \
      --input "$in_bcf" \
      --map "$gmap" \
      --region "$chr" \
      --reference "$ref_bcf" \
      --output "$out_bcf" \
      --thread "$THREADS"

    bcftools index -f "$out_bcf"
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

echo "Phasing finished"
echo "$(ls geno_ph/*ph.bcf.csi 2>/dev/null | wc -l) genotype bcf files were phased and saved in the folder geno_ph/:"
echo "$(ls geno_ph/*ph.bcf 2>/dev/null || true)"

printf "Total execution time: %02d:%02d:%02d\n"  $elapsed_h $elapsed_m $elapsed_s
echo "=^.^="


