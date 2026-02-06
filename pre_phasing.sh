#!/bin/bash
set -euo pipefail

## PREPHASING (build sample/reference per-chr BCFs + intersect variants)

start_time=$(date +%s)
exec >> pre_phasing.log 2>&1
echo "Prephasing started at $(date)"

CHR_RAW=""
THREADS=8
MAX_JOBS=4
SAMPLE_VCF=""
CONTIGS_MAP="contigs.txt"
REF_DIR="ref_1kGP"

usage() {
  echo "Usage: bash pre_phasing.sh -sample geno.vcf.gz -chr \"1-22\"|\"1-3,5,22\" [-threads 8] [-jobs 4]"
  exit 1
}

# Parse CLI arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -sample)
      [[ -z "${2:-}" || "${2:-}" == -* ]] && { echo "ERROR: -sample requires an argument"; usage; }
      SAMPLE_VCF="$2"
      shift 2
      ;;
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

[[ -z "$SAMPLE_VCF" ]] && { echo "ERROR: -sample is required"; usage; }
[[ -z "$CHR_RAW" ]] && { echo "ERROR: -chr is required"; usage; }

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS < 1 )); then
  echo "ERROR: -threads must be an integer >= 1"
  exit 1
fi

if ! [[ "$MAX_JOBS" =~ ^[0-9]+$ ]] || (( MAX_JOBS < 1 )); then
  echo "ERROR: -jobs must be an integer >= 1"
  exit 1
fi

if [[ ! -f "$CONTIGS_MAP" ]]; then
  echo "ERROR: contigs mapping file not found: $CONTIGS_MAP"
  exit 1
fi

if [[ ! -d "$REF_DIR" ]]; then
  echo "ERROR: reference directory not found: $REF_DIR"
  exit 1
fi

if [[ ! -f "$SAMPLE_VCF" ]]; then
  echo "ERROR: sample file not found: $SAMPLE_VCF"
  exit 1
fi

# Expand chromosome list
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
echo "Threads (-threads): $THREADS"
echo "Sample VCF: $SAMPLE_VCF"
echo "Reference dir: $REF_DIR"
echo "Contigs map: $CONTIGS_MAP"

# Create output dirs
mkdir -p geno_chr

# Index sample VCF (force to be safe across machines)
bcftools index -f "$SAMPLE_VCF" --threads "$THREADS"

# Run per chromosome
running=0
for chr in "${CHRS[@]}"; do
  (
    sample_bcf="geno_chr/geno_${chr}.bcf"
    sample_com_bcf="geno_chr/geno_${chr}_com.bcf"

    ref_vcfgz="${REF_DIR}/ref_${chr}.vcf.gz"
    ref_bcf="${REF_DIR}/ref_${chr}.bcf"
    ref_com_bcf="${REF_DIR}/ref_${chr}_com.bcf"

    # Check reference VCF exists
    if [[ ! -f "$ref_vcfgz" ]]; then
      echo "ERROR: reference file not found: $ref_vcfgz"
      exit 1
    fi

    # 1) Extract chromosome from sample -> BCF + index
    bcftools view --threads "$THREADS" -r "$chr" "$SAMPLE_VCF" -O b -o "$sample_bcf"
    bcftools index -f "$sample_bcf" --threads "$THREADS"

    # 2) Rename contigs in reference -> BCF + index
    # (Reads ref_${chr}.vcf.gz and produces ref_${chr}.bcf)
    bcftools annotate --rename-chrs "$CONTIGS_MAP" --threads "$THREADS" \
      -O b -o "$ref_bcf" "$ref_vcfgz"
    bcftools index -f "$ref_bcf" --threads "$THREADS"

    # 3) Intersect variants (create *_com files) + index
    bcftools isec -n=2 -w1 --threads "$THREADS" -O b -o "$sample_com_bcf" \
      "$sample_bcf" "$ref_bcf"
    bcftools index -f "$sample_com_bcf" --threads "$THREADS"

    bcftools isec -n=2 -w2 --threads "$THREADS" -O b -o "$ref_com_bcf" \
      "$sample_bcf" "$ref_bcf"
    bcftools index -f "$ref_com_bcf" --threads "$THREADS"

  ) &

  ((running+=1))
  if (( running >= MAX_JOBS )); then
    wait
    running=0
  fi
done

wait

end_time=$(date +%s)
elapsed=$((end_time-start_time))

printf "Prephasing finished in %02d:%02d:%02d\n" \
  $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

echo "=^.^="
