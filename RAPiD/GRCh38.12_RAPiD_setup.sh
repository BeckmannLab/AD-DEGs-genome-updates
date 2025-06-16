#!/bin/bash

##########
# Copy FASTA and GTF files
##########
mkdir -p "$fasta_dir" "$gtf_dir"
cp /sc/arion/projects/H_PBG/REFERENCES/GRCh38/FASTA/GRCh38.primary_assembly.genome.fa "$fasta_file"
cp /sc/arion/projects/H_PBG/REFERENCES/GRCh38/Gencode/release_30/gencode.v30.primary_assembly.annotation.gtf "$gtf_file"
wget -P "$fasta_dir" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh38.primary_assembly.genome.fa.gz"

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
  --outFileNamePrefix "${base_dir}/GRCh38.gencode.v30"

##########
# Compare STAR output with reference
##########
ref_star_dir="/sc/arion/projects/mscic1/results/anina/Noam_testing/chr_primary"
cd "$star_dir"
for file in *; do
  sum1=$(md5sum "$file" | awk '{print $1}')
  sum2=$(md5sum "${ref_star_dir}/$file" | awk '{print $1}')
  if [ "$sum1" = "$sum2" ]; then
    echo "files matching $file"
  else
    echo "files not matching $file"
  fi
done

##########
# Convert GTF to GenePred and reflat
##########
cd "$gtf_dir"
/sc/arion/projects/mscic1/results/anina/gtfToGenePred -genePredExt "$gtf_file" "$genepred_file"

awk '($2 !~ /^H/){print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' \
  "$genepred_file" > "$reflat_file"

##########
# Create rRNA BED and interval list
##########
awk 'OFS="\t" { if ($3 == "CDS") print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" | tr -d '";' > "$bed_file"
awk 'OFS="\t" { if ($2 == "rRNA" && $1 !~ /^H/) print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" | tr -d '";' > "$rrna_bed"

ml gatk
gatk CreateSequenceDictionary -R "$fasta_file"

module load picard/2.18.4
java -jar "$PICARD" BedToIntervalList \
  I="$rrna_bed" \
  O="${gtf_dir}/gencode.v30.rRNA.interval_list" \
  SD="$dict_file"

ln -s "${gtf_dir}/gencode.v30.rRNA.interval_list" "${gtf_dir}/gencode.v30.rRNA.interval.list"

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
cp "$globin_gene_list" ./

cat globin_gene | xargs -n1 -I% grep % "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "CDS") { if ($1 !~ /^H/) print $1, $4-1, $5, $10, 0, $7; else print $1, $4-1, $5, $10, 0, $7 } }' \
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
