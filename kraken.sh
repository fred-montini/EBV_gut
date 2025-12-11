
set -euo pipefail
SAMPLE="${1:-GSM7434453_HC_14}"
R1="${2}"
R2="${3}"
OUT="${4:-out}"
DB="${5:-kraken2_db}"
HUMAN_FA="${6:-hg38.fa}"
THREADS="${THREADS:-16}"
mkdir -p "$OUT" "$DB" "$(dirname "$HUMAN_FA")" "${OUT}/tmp_${SAMPLE}"
if [ ! -f "${HUMAN_FA}.bitmask" ]; then bmtool -d "$HUMAN_FA" -o "${HUMAN_FA}.bitmask" -A 0 -w 18; fi
if [ ! -f "${HUMAN_FA}.srprism.idx" ] && [ ! -f "${HUMAN_FA}.srprism" ]; then srprism mkindex -i "$HUMAN_FA" -o "${HUMAN_FA}.srprism" -M 7168; fi
if [ ! -f "${HUMAN_FA}.nsq" ]; then makeblastdb -in "$HUMAN_FA" -dbtype nucl; fi
bmtagger.sh -q 1 -1 "$R1" -2 "$R2" -b "${HUMAN_FA}.bitmask" -x "${HUMAN_FA}.srprism" -d "$HUMAN_FA" -T "${OUT}/tmp_${SAMPLE}" -o "${OUT}/${SAMPLE}.human_blacklist.txt"
FILTER_FASTQ() { BL="$1"; IN="$2"; OUTF="$3"; if printf "%s" "$IN" | grep -qE "\.gz$"; then DECOMP="zcat"; else DECOMP="cat"; fi; awk 'NR==FNR{a[$1]=1;next} /^@/{h=$1;sub(/^@/,"",h);split(h,s," ");id=s[1];getline seq;getline plus;getline qual;if(!(id in a)){print $0;print seq;print plus;print qual}}' "$BL" <($DECOMP "$IN") | gzip -c > "$OUTF"; }
FILTER_FASTQ "${OUT}/${SAMPLE}.human_blacklist.txt" "$R1" "${OUT}/${SAMPLE}.nonhuman_R1.fastq.gz"
FILTER_FASTQ "${OUT}/${SAMPLE}.human_blacklist.txt" "$R2" "${OUT}/${SAMPLE}.nonhuman_R2.fastq.gz"
if [ ! -f "${DB}/taxo.k2d" ]; then kraken2-build --download-taxonomy --db "$DB" --threads "$THREADS"; for L in archaea bacteria plasmid viral human UniVec_Core; do kraken2-build --download-library "$L" --db "$DB" --threads "$THREADS"; done; kraken2-build --build --db "$DB" --threads "$THREADS"; fi
kraken2 --db "$DB" --threads "$THREADS" --use-names --paired --report "${OUT}/${SAMPLE}.kraken2.report" "${OUT}/${SAMPLE}.nonhuman_R1.fastq.gz" "${OUT}/${SAMPLE}.nonhuman_R2.fastq.gz" > "${OUT}/${SAMPLE}_results.txt"
gzip -c "${OUT}/${SAMPLE}_results.txt" > "${OUT}/${SAMPLE}_results.txt.gz"
awk 'BEGIN{FS="\t";OFS="\t"}{n=$6;gsub(/^ +/,"",n);if($4=="D"){dom=n}if(dom=="Viruses"&&$4=="O"){print n,$2}}' "${OUT}/${SAMPLE}.kraken2.report" > "${OUT}/${SAMPLE}.virome_order.tsv"
awk 'BEGIN{FS="\t";OFS="\t"}{n=$6;gsub(/^ +/,"",n);if($4=="D"){dom=n}if(dom=="Viruses"&&$4=="F"){print n,$2}}' "${OUT}/${SAMPLE}.kraken2.report" > "${OUT}/${SAMPLE}.virome_family.tsv"
