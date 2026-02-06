#!/bin/bash

set -euo pipefail

## LOCAL ANCESTRY INFERENCE

# Startime of execution for computing elapsed time
start_time=$(date +%s)

# Creating the la_inf.log file, which will contain all stdout and stderr lines
exec >> la_inf.log 2>&1

# Starting date and time
echo "Local Ancestry Inference started at $(date)"

# Default values
RAW_ANC=""
CHR_RAW=""          # REQUIRED
MAX_JOBS=4
THREADS=8

usage() {
  echo "Usage: bash la_inf.sh -ancs anc1,anc2,... -chr \"1,4,22\"|\"1-3,5,22\" [-jobs 4] [-threads 8]"
  exit 1
}

# Parse command line arguments (single pass)
while [[ $# -gt 0 ]]; do
  case "$1" in
    -ancs)
      if [[ -z "${2:-}" || "${2:-}" == -* ]]; then
        echo "Error: '-ancs' requires a non-empty argument."
        usage
      fi
      RAW_ANC="$2"
      shift 2
      ;;
    -chr)
      if [[ -z "${2:-}" || "${2:-}" == -* ]]; then
        echo "Error: '-chr' requires a non-empty argument."
        usage
      fi
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
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown parameter passed: $1"
      usage
      ;;
  esac
done

# Required args checks
if [[ -z "$RAW_ANC" ]]; then
  echo "Error: -ancs argument is required."
  usage
fi

if [[ -z "$CHR_RAW" ]]; then
  echo "Error: -chr argument is required."
  usage
fi

# Numeric validation
if ! [[ "$MAX_JOBS" =~ ^[0-9]+$ ]] || (( MAX_JOBS < 1 )); then
  echo "ERROR: -jobs must be an integer >= 1"
  exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS < 1 )); then
  echo "ERROR: -threads must be an integer >= 1"
  exit 1
fi

# Convert to uppercase and wrap each element in single quotes
IFS=',' read -ra ANC_ITEMS <<< "$RAW_ANC"

NORMALIZED_ANC_LIST=""
for item in "${ANC_ITEMS[@]}"; do
    UPPER=$(echo "$item" | tr '[:lower:]' '[:upper:]')
    if [[ -n "$NORMALIZED_ANC_LIST" ]]; then
        NORMALIZED_ANC_LIST="$NORMALIZED_ANC_LIST,'$UPPER'"
    else
        NORMALIZED_ANC_LIST="'$UPPER'"
    fi
done

# Call the Python script ($NORMALIZED_ANC_LIST is the list of variables in a standarized way to be red by sel_anc.py 
python3 sel_anc.py -ancs "$NORMALIZED_ANC_LIST"

# Samplemap file
cleaned=$(echo "$NORMALIZED_ANC_LIST" | tr -d "'")
ancs=$(echo "$cleaned" | tr '[:upper:]' '[:lower:]' | tr ',' '.')
filename="b38_${ancs}_map.txt"

# Correcting gmaps, LA inference

# Expand chromosome list (supports ranges) Example: "1-3,5,22" -> 1 2 3 5 22
expand_chr_list() {
  local raw="$1"
  raw="$(echo "$raw" | tr -d '[:space:]')"  # Remove whitespace
  local IFS=','
  read -r -a parts <<< "$raw"

  local out=()
  for token in "${parts[@]}"; do
    [[ -z "$token" ]] && continue

    if [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      local start="${BASH_REMATCH[1]}"
      local end="${BASH_REMATCH[2]}"

      if (( start <= end )); then
        for ((c=start; c<=end; c++)); do out+=("$c"); done
      else
        for ((c=start; c>=end; c--)); do out+=("$c"); done
      fi

    elif [[ "$token" =~ ^[0-9]+$ ]]; then
      out+=("$token")
    else
      echo "ERROR: Invalid token in -chr: '$token'"
      exit 1
    fi
  done

  # Remove duplicates and sort numerically
  printf "%s\n" "${out[@]}" | awk '!seen[$0]++' | sort -n
}

# Convert expanded list into array
mapfile -t CHRS < <(expand_chr_list "$CHR_RAW")

echo "Chromosomes: ${CHRS[*]}"
echo "Parallel jobs (-jobs): $MAX_JOBS"
echo "Threads per rfmix run (-threads): $THREADS"

# Creating the dec_chr folder
mkdir dec_chr/

# Run chromosomes in parallel
running=0
for chr in "${CHRS[@]}"; do
  (
    # Prepare genetic map
    awk '{print $2, $1, $3}' "gmap_1kGP/chr${chr}.b38.gmap" > "chr${chr}.b38.tmp.txt"
    tail -n +2 "chr${chr}.b38.tmp.txt" > "gmap_1kGP/chr${chr}.b38.txt"
    rm -f "chr${chr}.b38.tmp.txt"

    # Run RFMix
    rfmix -f "geno_ph/geno_${chr}_ph.bcf" \
          -r "ref_1kGP/ref_${chr}_com.bcf" \
          -m "${filename}" \
          -g "gmap_1kGP/chr${chr}.b38.txt" \
          -o "dec_chr/dec_${chr}" \
          --n-threads="${THREADS}" \
          --chromosome="${chr}"
  ) &

  ((running+=1))

  # Wait when reaching max parallel jobs
  if (( running >= MAX_JOBS )); then
    wait
    running=0
  fi
done

wait



# Defining end time of execution and computing elapsed time in hh:mm:ss format 

end_time=$(date +%s)

elapsed_time=$((end_time-start_time))
elapsed_h=$((elapsed_time/3600))
elapsed_m=$(((elapsed_time % 3600)/60))
elapsed_s=$((elapsed_time % 60))

# Final messages
echo "Local Ancestry Inference finished"

echo "$(ls dec_chr/*dec_chr* | wc -l) output files were created and saved in the folder dec_chr/:"
echo "$(ls dec_chr/*dec_chr*)"

printf "Total execution time: %02d:%02d:%02d\n"  $elapsed_h $elapsed_m $elapsed_s
echo "=^.^="








