#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  align_bowtie2.sh -x <bt2_index_prefix> -1 <R1.fastq.gz> -2 <R2.fastq.gz> -o <out_prefix> [-t <threads>] [--rg <readgroup>] [--preset <name>]

Required:
  -x  Bowtie2 index prefix (e.g., /path/genome or genome if genome.1.bt2 etc exist)
  -1  R1 FASTQ (optionally gzipped)
  -2  R2 FASTQ (optionally gzipped)
  -o  Output prefix (no extension). Creates:
        <out_prefix>.sorted.bam
        <out_prefix>.sorted.bam.bai
        <out_prefix>.bowtie2.log
        <out_prefix>.flagstat.txt

Optional:
  -t        Threads (default: 4)
  --rg      Read group string for bowtie2, e.g.:
              "ID:lane1\tSM:sample1\tPL:ILLUMINA"
  --preset  Bowtie2 preset (passed through directly):
              very-fast, fast, sensitive (default), very-sensitive

Notes:
  - No MAPQ filtering is applied
EOF
}

# Defaults
threads=4
rg=""
preset="sensitive"

# Parse args
genome=""
R1=""
R2=""
out_prefix=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -x) genome="$2"; shift 2 ;;
    -1) R1="$2"; shift 2 ;;
    -2) R2="$2"; shift 2 ;;
    -o) out_prefix="$2"; shift 2 ;;
    -t) threads="$2"; shift 2 ;;
    --rg) rg="$2"; shift 2 ;;
    --preset) preset="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

# Validate
if [[ -z "${genome}" || -z "${R1}" || -z "${R2}" || -z "${out_prefix}" ]]; then
  echo "ERROR: Missing required arguments." >&2
  usage
  exit 2
fi

for f in "${R1}" "${R2}"; do
  [[ -r "$f" ]] || { echo "ERROR: Cannot read file: $f" >&2; exit 2; }
done

# Check bowtie2 index existence (bt2 or bt2l)
if ! ls "${genome}".*.bt2* >/dev/null 2>&1; then
  echo "ERROR: Bowtie2 index not found for prefix: ${genome}" >&2
  echo "Expected files like: ${genome}.1.bt2 (or .bt2l)" >&2
  exit 2
fi

# Tool checks
command -v bowtie2 >/dev/null || { echo "ERROR: bowtie2 not in PATH" >&2; exit 2; }
command -v samtools >/dev/null || { echo "ERROR: samtools not in PATH" >&2; exit 2; }

# Outputs
log="${out_prefix}.bowtie2.log"
sorted_bam="${out_prefix}.sorted.bam"
flagstat="${out_prefix}.flagstat.txt"

# Temp dir (auto-clean)
tmpdir="$(mktemp -d)"
cleanup() { rm -rf "${tmpdir}"; }
trap cleanup EXIT

# Map our preset name to Bowtie2's preset flag
preset_flag=()
case "${preset}" in
  very-fast)      preset_flag=(--very-fast) ;;
  fast)           preset_flag=(--fast) ;;
  sensitive)      preset_flag=(--sensitive) ;;
  very-sensitive) preset_flag=(--very-sensitive) ;;
  *)
    echo "ERROR: Unknown --preset '${preset}'" >&2
    echo "Valid: very-fast, fast, sensitive, very-sensitive" >&2
    exit 2
    ;;
esac

echo "[$(date)] Starting alignment"
echo "  Index:    ${genome}"
echo "  R1:       ${R1}"
echo "  R2:       ${R2}"
echo "  Threads:  ${threads}"
echo "  Out:      ${sorted_bam}"
echo "  Preset:   ${preset} (${preset_flag[*]})"

# Optional read group
rg_args=()
if [[ -n "${rg}" ]]; then
  rg_args=(--rg "${rg}")
fi

# Align -> coordinate sort (streaming) -> index -> stats
bowtie2 -x "${genome}" \
  -1 "${R1}" \
  -2 "${R2}" \
  --threads "${threads}" \
  "${preset_flag[@]}" \
  "${rg_args[@]}" \
  2> "${log}" \
| samtools sort \
    -@ "${threads}" \
    -T "${tmpdir}/sorttmp" \
    -o "${sorted_bam}" -

samtools index -@ "${threads}" "${sorted_bam}"
samtools flagstat -@ "${threads}" "${sorted_bam}" > "${flagstat}"

echo "[$(date)] Done"
echo "  BAM:      ${sorted_bam}"
echo "  BAI:      ${sorted_bam}.bai"
echo "  Bowtie2:  ${log}"
echo "  Flagstat: ${flagstat}"
