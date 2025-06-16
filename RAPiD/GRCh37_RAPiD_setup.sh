#!/bin/bash

##########
# Download Ensembl v70 FASTA and GTF
##########
mkdir -p "$fasta_dir" "$gtf_dir"

wget -P "$gtf_dir" \
  "http://ftp.ensembl.org/pub/release-70/gtf/homo_sapiens/Homo_sapiens.GRCh37.70.gtf.gz"

wget -P "$fasta_dir" \
  "http://ftp.ensembl.org/pub/release-70/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.70.dna_sm.primary_assembly.fa.gz"

##########
# Generate STAR genome index
##########
cd "$base_dir"
ml star/2.7.3a
mkdir -p "$star_dir"

STAR --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir "$star_dir" \
  --genomeFastaFiles "$fasta_file" \
  --sjdbGTFfile "$gtf_file" \
  --sjdbOverhang 100 \
  --outFileNamePrefix "${base_dir}/${output_name}"

##########
# Generate GenePred and reflat
##########
cd "$base_dir"
"$gtftogenepred" -genePredExt "$gtf_file" "${output_name}.GenePred"

awk '($2 !~ /^H/){print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' \
  "${output_name}.GenePred" > "${output_name}.reflat"

##########
# Create rRNA BED and interval list
##########
awk 'OFS="\t" { if ($3 == "CDS") print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" | tr -d '";' > "${output_name}.bed"
awk 'OFS="\t" { print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" | tr -d '";' > "check${output_name}.bed"
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
# Create MtRNA gene list and annotation
##########
grep "^MT" "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "CDS") print $16 }' | tr -d '";' > MtRNA_gene

grep "^MT" "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "CDS") print $1, $4-1, $5, $10, 0, $7 }' | tr -d '";' > MtRNA_gene.bed

ml bedtools
bedtools merge -s -c 4,5,6 -delim "|" -o distinct -i MtRNA_gene.bed > MtRNA_gene.merge.bed

ml ucsc-utils/2020-03-17
bedToGenePred MtRNA_gene.merge.bed MtRNA_gene.merge.genePred
genePredToGtf file MtRNA_gene.merge.genePred MtRNA_gene.merge.gtf
awk '$3 == "exon"' MtRNA_gene.merge.gtf > MtRNA_gene.gtf

##########
# Create globin gene annotation
##########
cp "$globin_gene_file" ./

cat globin_gene | xargs -n1 -I% grep % "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "CDS") { if ($1 !~ /^H/) print $1, $4-1, $5, $10, 0, $7; else print $1, $4-1, $5, $10, 0, $7 } }' \
  | tr -d '";' | sort -u > globin_gene.bed

bedtools sort -i globin_gene.bed > globin_gene.sort.bed
bedtools merge -s -c 4,5,6 -delim "|" -o distinct -i globin_gene.sort.bed > globin_gene.merge.bed
bedToGenePred globin_gene.merge.bed globin_gene.merge.genePred
genePredToGtf file globin_gene.merge.genePred globin_gene.merge.gtf
awk '$3 == "exon"' globin_gene.merge.gtf > globin_gene.gtf

##########
# Clean up intermediate files
##########
rm -f MtRNA_gene.bed MtRNA_gene.merge.* \
      globin_gene.bed globin_gene.sort.bed globin_gene.merge.*
