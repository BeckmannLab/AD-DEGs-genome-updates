#!/bin/bash

########################################
# Setup GRCh38 Gencode v43 Reference
########################################

##########
# Define top-level directories and input paths
##########
base_dir="/sc/arion/projects/mscic1/results/anina/references/GRCh38/v43"
fasta_dir="${base_dir}/FASTA"
gtf_dir="${base_dir}/gencode.v43"
star_dir="${base_dir}/chr_primary"

fasta_file="${fasta_dir}/GRCh38.primary_assembly.genome.fa"
gtf_file="${gtf_dir}/gencode.v43.primary_assembly.annotation.gtf"
genepred_file="${gtf_dir}/gencode.v43.primary_assembly.annotation.GenePred"
reflat_file="${gtf_dir}/gencode.v43.primary_assembly.annotation.reflat"
bed_file="${gtf_dir}/gencode.v43.primary_assembly.annotation.bed"
rrna_bed="${gtf_dir}/gencode.v43.rRNA.bed"
dict_file="/sc/arion/projects/mscic1/results/anina/Noam_testing/FASTA/GRCh38.primary_assembly.genome.dict"
globin_gene_list="/sc/arion/projects/H_PBG/REFERENCES/GRCh38/Gencode/release_30/globin_gene"

##########
# Download Gencode v43 FASTA and GTF
##########
mkdir -p "$fasta_dir" "$gtf_dir"
wget -P "$fasta_dir" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz"
wget -P "$gtf_dir" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz"

##########
# Generate STAR index
##########
ml star/2.7.3a
mkdir -p "$star_dir"

STAR --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir "$star_dir" \
  --genomeFastaFiles "$fasta_file" \
  --sjdbGTFfile "$gtf_file" \
  --sjdbOverhang 100 \
  --outFileNamePrefix "${star_dir}/chr_primary"

##########
# Create GenePred and reflat
##########
cd "$base_dir"
/sc/arion/projects/mscic1/results/anina/gtfToGenePred -genePredExt "$gtf_file" "$genepred_file"

awk '{ print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10 }' OFS='\t' \
  "$genepred_file" > "$reflat_file"

##########
# Create rRNA BED and interval list
##########
awk 'OFS="\t" { if ($3 == "gene") print $1, $4-1, $5, $10, $12, $7 }' "$gtf_file" | tr -d '";' > "$bed_file"
grep rRNA "$bed_file" > "$rrna_bed"

ml gatk
gatk CreateSequenceDictionary -R "$fasta_file"

module load picard/2.18.4
java -jar "$PICARD" BedToIntervalList \
  I="$rrna_bed" \
  O="${gtf_dir}/gencode.v43.rRNA.interval_list" \
  SD="$dict_file"

ln -s "${gtf_dir}/gencode.v43.rRNA.interval_list" "${gtf_dir}/gencode.v43.rRNA.interval.list"

##########
# Create MtRNA gene list and annotation
##########
grep chrM "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "gene") print $14 }' | tr -d '";' > MtRNA_gene

grep chrM "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "gene") print $1, $4-1, $5, $10, 0, $7 }' | tr -d '";' > MtRNA_gene.bed

ml bedtools
bedtools merge -s -c 4,5,6 -delim "|" -o distinct -i MtRNA_gene.bed > MtRNA_gene.merge.bed

ml ucsc-utils/2020-03-17
bedToGenePred MtRNA_gene.merge.bed MtRNA_gene.merge.genePred
genePredToGtf file MtRNA_gene.merge.genePred MtRNA_gene.merge.gtf
awk '$3 == "exon"' MtRNA_gene.merge.gtf > MtRNA_gene.gtf

##########
# Create globin gene annotation
##########
cp "$globin_gene_list" ./

cat globin_gene | xargs -n1 -I% grep % "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "gene") print $1, $4-1, $5, $10, 0, $7 }' \
  | tr -d '";' | sort -u > globin_gene.bed

bedtools sort -i globin_gene.bed > globin_gene.sort.bed
bedtools merge -s -c 4,5,6 -delim "|" -o distinct -i globin_gene.sort.bed > globin_gene.merge.bed
bedToGenePred globin_gene.merge.bed globin_gene.merge.genePred
genePredToGtf file globin_gene.merge.genePred globin_gene.merge.gtf
awk '$3 == "exon"' globin_gene.merge.gtf > globin_gene.gtf

##########
# Cleanup intermediate files
##########
rm -f MtRNA_gene.bed MtRNA_gene.merge.* \
       globin_gene.bed globin_gene.sort.bed globin_gene.merge.*
