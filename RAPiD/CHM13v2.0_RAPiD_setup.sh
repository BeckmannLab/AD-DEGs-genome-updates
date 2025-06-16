#!/bin/bash

##########
# Build STAR index (CDS version)
##########
cd "$t2t_root"
ml star/2.7.3a

STAR --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir "$star_dir" \
  --genomeFastaFiles "$cdna_file" \
  --limitGenomeGenerateRAM 73287799609 \
  --sjdbGTFfile "$gtf_file" \
  --sjdbOverhang 100 \
  --outFileNamePrefix "${t2t_root}/GCA_009914755.4"

##########
# Build STAR index (unmasked genome version)
##########
STAR --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir "$star_dir" \
  --genomeFastaFiles "$fasta_file" \
  --limitGenomeGenerateRAM 73287799609 \
  --sjdbGTFfile "$gtf_file" \
  --sjdbOverhang 100 \
  --outFileNamePrefix "${t2t_root}/GCA_009914755.4_cdna"

##########
# Generate GenePred and reflat
##########
cd "$t2t_root"
"$gtftogenepred" -genePredExt "$gtf_file" "${output_name}.GenePred"

awk '($2 !~ /^H/){print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' \
  "${output_name}.GenePred" > "${output_name}.reflat"

##########
# Create rRNA BED and interval list
##########
awk 'OFS="\t" { if ($3 == "CDS") print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" | tr -d '";' > "${output_name}.bed"
awk 'OFS="\t" { if ($2 == "rRNA" && $1 !~ /^H/) print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" | tr -d '";' > "${output_name}.rRNA.bed"

ml gatk
gatk CreateSequenceDictionary -R "$fasta_file"

module load picard/2.18.4
java -jar "$PICARD" BedToIntervalList \
  I="${output_name}.rRNA.bed" \
  O="${output_name}.rRNA.interval_list" \
  SD="$dict_file"

ln -s "${output_name}.rRNA.interval_list" "${output_name}.rRNA.interval.list"

##########
# Create MtRNA gene annotation (no globin in T2T)
##########
grep "^MT" "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "CDS") print $16 }' | tr -d '";' > MtRNA_gene

grep "^MT" "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "CDS") print $1, $4-1, $5, $10, 0, $7 }' | tr -d '";' > MtRNA_gene.bed

ml bedtools
bedtools sort -i MtRNA_gene.bed | \
  bedtools merge -s -c 4,5,6 -delim "|" -o distinct > MtRNA_gene.merge.bed

# Additional processing (e.g., convert to GTF if needed)
