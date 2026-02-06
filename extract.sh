#!/bin/bash
set -euo pipefail

## EXTRACTING TRACTS

start_time=$(date +%s)

exec >> extract.log 2>&1
echo "Extracting tracts started at $(date)"

RAW_ANC=""
CHR_RAW=""
MAX_JOBS=22

usage() {
  echo "Usage: bash extract.sh -ancs anc1,anc2,... -chr \"1,4,22\"|\"1-3,5,22\" [-jobs 22]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -ancs)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -ancs requires an argument"; usage; }
      RAW_ANC="$2"
      shift 2
      ;;
    -chr)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -chr requires an argument"; usage; }
      CHR_RAW="$2"
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

[[ -z "$RAW_ANC" ]] && { echo "ERROR: -ancs is required"; usage; }
[[ -z "$CHR_RAW" ]] && { echo "ERROR: -chr is required"; usage; }

if ! [[ "$MAX_JOBS" =~ ^[0-9]+$ ]] || (( MAX_JOBS < 1 )); then
  echo "ERROR: -jobs must be an integer >= 1"
  exit 1
fi

RAW_ANC="$(echo "$RAW_ANC" | tr -d '[:space:]')"
IFS=',' read -r -a ANC_ITEMS <<< "$RAW_ANC"
NUM_ANCS="${#ANC_ITEMS[@]}"

if (( NUM_ANCS < 1 )); then
  echo "ERROR: -ancs must contain at least one ancestry (e.g., eur,afr)"
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
echo "Ancestries (-ancs): $RAW_ANC"
echo "Derived --num-ancs: $NUM_ANCS"

running=0
for chr in "${CHRS[@]}"; do
  (
    in_bcf="geno_ph/geno_${chr}_ph.bcf"
    final_vcf="geno_ph/geno_${chr}_ph.vcf.gz"

    # 1) BCF -> VCF.GZ
    if [[ ! -s "$final_vcf" ]]; then
      echo "[chr${chr}] Converting BCF -> VCF.GZ"
      bcftools view "$in_bcf" -Oz -o "$final_vcf"
    fi

    # 2) Index VCF.GZ
    if [[ ! -s "${final_vcf}.csi" && ! -s "${final_vcf}.tbi" ]]; then
      echo "[chr${chr}] Indexing $final_vcf"
      bcftools index -f "$final_vcf"
    fi

    # 3) Strip INFO and FORMAT (overwrite same file safely)
    echo "[chr${chr}] Stripping INFO,FORMAT"
    bcftools annotate -x INFO,FORMAT "$final_vcf" -Oz -o "${final_vcf}.tmp"
    mv "${final_vcf}.tmp" "$final_vcf"

    # 4) Re-index after overwrite
    bcftools index -f "$final_vcf"

    # 5) Extract tracts
    python3 extract_tracts.py \
      --vcf "geno_ph/geno_${chr}_ph.vcf.gz" \
      --msp "dec_chr/dec_${chr}.msp.tsv" \
      --num-ancs "${NUM_ANCS}"
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

echo "Extracting tracts finished"
printf "Total execution time: %02d:%02d:%02d\n"  $elapsed_h $elapsed_m $elapsed_s
echo "=^.^="











