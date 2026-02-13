#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME="$(basename "$0")"

# ------------------------------------------------------------
# Build ONE FASTQ set containing the UNION of:
#   A) "unmapped fragments": read unmapped OR mate unmapped
#   B) "decoy-associated fragments": either end maps to any decoy contig
#
# Decoy contigs are provided as a newline-separated list of contig names.
#
# Key decisions (explicit):
#   - Exclude from all selection steps: secondary, QC-fail, duplicates, supplementary (0xF00)
#   - Unmapped criterion uses: include-flags 0x4 OR 0x8 (12)
#   - Decoy criterion: include a pair if either end maps to a decoy contig
#     (samtools view --fetch-pairs on decoy regions)
#   - To avoid duplicates when concatenating FASTQs, the decoy extraction excludes
#     read-unmapped (0x4) AND mate-unmapped (0x8). Decoy pairs with an unmapped
#     mate are still included via criterion (A).
# ------------------------------------------------------------

usage() {
  cat <<EOF
Usage:
  ${SCRIPT_NAME} --bam FILE.bam --decoys contigs.txt --prefix NAME [options]

Required:
  --bam FILE.bam          Input BAM (indexed recommended)
  --decoys contigs.txt    Newline-separated list of decoy contig names
  --prefix NAME           Output file prefix (written into --outdir)

Options:
  -t, --threads N         Threads for samtools (-@). Default: 1
  -o, --outdir DIR        Output directory. Default: .
  --keep-temp             Keep intermediate FASTQs (unmapped_only.* and decoy.*). Default: remove them

Behavior:
  Produces a single FASTQ set containing all read pairs/fragments that satisfy either:
    (A) read unmapped OR mate unmapped (flags 0x4 or 0x8)
    (B) either end maps to a decoy contig (pair-complete via --fetch-pairs)

  Selection excludes 0xF00 (secondary/QC-fail/duplicate/supplementary).

Requirements:
  samtools and bgzip must be available in PATH.
EOF
}

die() {
  echo "ERROR: $*" >&2
  usage >&2
  exit 1
}

# Defaults
threads=1
outdir="."
bam=""
decoys=""
prefix=""
keep_temp=0

# -------------------------
# Parse arguments
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
  --bam)
    bam="${2:-}"
    shift 2
    ;;
  --decoys)
    decoys="${2:-}"
    shift 2
    ;;
  --prefix)
    prefix="${2:-}"
    shift 2
    ;;
  -t | --threads)
    threads="${2:-}"
    shift 2
    ;;
  -o | --outdir)
    outdir="${2:-}"
    shift 2
    ;;
  --keep-temp)
    keep_temp=1
    shift
    ;;
  -h | --help)
    usage
    exit 0
    ;;
  *) die "Unknown argument: $1 (use --help)" ;;
  esac
done

# -------------------------
# Validate inputs
# -------------------------
command -v samtools >/dev/null 2>&1 || die "samtools not found in PATH"
command -v bgzip >/dev/null 2>&1 || die "bgzip not found in PATH"

[[ -n "$bam" ]] || die "--bam is required"
[[ -n "$decoys" ]] || die "--decoys is required"
[[ -n "$prefix" ]] || die "--prefix is required"

[[ -f "$bam" ]] || die "BAM not found: $bam"
[[ -f "$decoys" ]] || die "Decoy contig file not found: $decoys"

mkdir -p "$outdir"

[[ "$threads" =~ ^[0-9]+$ ]] || die "--threads must be an integer"
[[ "$threads" -ge 0 ]] || die "--threads must be >= 0"

# Load decoy contigs: ignore blanks and comments
decoy_contigs=()
while IFS= read -r line; do
  decoy_contigs+=("$line")
done < <(grep -vE '^\s*($|#)' "$decoys" || true)

[[ "${#decoy_contigs[@]}" -gt 0 ]] || die "No contig names found in $decoys"

# Warn if BAM index missing
if [[ ! -f "${bam}.bai" && ! -f "${bam%.*}.bai" ]]; then
  echo "WARNING: BAM index (.bai) not found. Please ensure bam is indexed and try again" >&2
  exit 1
fi

# Output files
out_r1_unmapped="${outdir}/${prefix}.unmapped_only.R1.fq"
out_r2_unmapped="${outdir}/${prefix}.unmapped_only.R2.fq"
out_s_unmapped="${outdir}/${prefix}.unmapped_only.singleton.fq"
out_r1_decoy="${outdir}/${prefix}.decoy.R1.fq"
out_r2_decoy="${outdir}/${prefix}.decoy.R2.fq"
out_s_decoy="${outdir}/${prefix}.decoy.singleton.fq"

out_r1_final="${outdir}/${prefix}.R1.fq.gz"
out_r2_final="${outdir}/${prefix}.R2.fq.gz"
out_s_final="${outdir}/${prefix}.singleton.fq.gz"

# -------------------------
# Unmapped fragments
# -------------------------
echo "[1/3] Writing unmapped fragments (read unmapped OR mate unmapped)" >&2
# Note we don't need keep-pairs because the --include-flags includes read unmapped or mate unmapped

samtools view \
  -u \
  --exclude-flags 3840 \
  --include-flags 12 \
  "$bam" |
  samtools collate -O - |
  samtools fastq \
    -1 "${out_r1_unmapped}" \
    -2 "${out_r2_unmapped}" \
    -s "${out_s_unmapped}" \
    -@ "${threads}"

# -------------------------
# Decoy-associated fragments
# -------------------------
echo "[2/3] Writing decoy-associated fragments (either end on decoy contigs)" >&2
samtools view \
  -u \
  --fetch-pairs \
  --exclude-flags 3852 \
  "$bam" \
  "${decoy_contigs[@]}" |
  samtools collate -O - |
  samtools fastq \
    -1 "${out_r1_decoy}" \
    -2 "${out_r2_decoy}" \
    -s "${out_s_decoy}" \
    -@ "${threads}"

# -------------------------
# Combine FASTQs
# -------------------------
echo "[3/3] Combining FASTQs into a single bgzipped record..." >&2
cat "${out_r1_unmapped}" "${out_r1_decoy}" | bgzip -c >"${out_r1_final}"
cat "${out_r2_unmapped}" "${out_r2_decoy}" | bgzip -c >"${out_r2_final}"
cat "${out_s_unmapped}" "${out_s_decoy}" | bgzip -c >"${out_s_final}"

# Cleanup
if [[ "${keep_temp}" -eq 0 ]]; then
  rm -f "${out_r1_unmapped}" "${out_r2_unmapped}" "${out_s_unmapped}" \
    "${out_r1_decoy}" "${out_r2_decoy}" "${out_s_decoy}"
else
  echo "Keeping intermediate FASTQs (--keep-temp set)." >&2
fi

echo "Done." >&2
