#!/bin/bash

######################################################
# Define all paths and files
######################################################
base_dir="/sc/arion/projects/mscic1/results/anina/references/GRCh38/v30"
ref_dir="/sc/arion/projects/H_PBG/REFERENCES/GRCh38"
gtf_src="${ref_dir}/Gencode/release_30/gencode.v30.primary_assembly.annotation.gtf"
fasta_src="${ref_dir}/FASTA/GRCh38.primary_assembly.genome.fa"
download_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh38.primary_assembly.genome.fa.gz"

gtf_dest="${base_dir}/gencode.v30.primary_assembly.annotation.gtf"
fasta_dest="${base_dir}/FASTA/GRCh38.primary_assembly.genome.fa"
star_index_dir="${base_dir}/chr_primary"
md5_compare_dir="/sc/arion/projects/mscic1/results/anina/Noam_testing/chr_primary"

gtf_file="/sc/arion/projects/mscic1/results/anina/references/GRCh37/v70/ensembl.v30/Homo_sapiens.GRCh37.70.gtf"
gtf_dir="/sc/arion/projects/mscic1/results/anina/references/GRCh37/v70"
gtf_base=$(basename "$gtf_file" .gtf)
gtftogenepred="/sc/arion/projects/mscic1/results/anina/gtfToGenePred"
dict_file="/sc/arion/projects/mscic1/results/anina/references/GRCh37/v70/FASTA/Homo_sapiens.GRCh37.70.dna_sm.primary_assembly.dict"
globin_gene_list="${ref_dir}/Gencode/release_30/globin_gene"

######################################################
# Fetch FASTA and GTF
######################################################
cd "$base_dir"
cp "$fasta_src" "$fasta_dest"
cp "$gtf_src" "$gtf_dest"
wget -P "${base_dir}/FASTA" "$download_url"

######################################################
# Generate STAR index
######################################################
ml star/2.7.3a
mkdir -p "$star_index_dir"

STAR --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir "$star_index_dir" \
  --genomeFastaFiles "$fasta_dest" \
  --sjdbGTFfile "$gtf_dest" \
  --sjdbOverhang 100 \
  --outFileNamePrefix GRCh38.gencode.v30

######################################################
# Compare STAR-generated files with trusted versions
######################################################
cd "$star_index_dir"
for file in *; do
  sum1=$(md5sum "$file" | awk '{print $1}')
  sum2=$(md5sum "${md5_compare_dir}/$file" | awk '{print $1}')
  if [ "$sum1" = "$sum2" ]; then
    echo "files matching $file"
  else
    echo "files not matching $file"
  fi
done

######################################################
# Convert GTF to GenePred and Reflat
######################################################
cd "$gtf_dir"
"$gtftogenepred" -genePredExt "$gtf_file" "${gtf_base}.GenePred"

awk '($2 !~ /^H/){print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' \
  "${gtf_base}.GenePred" > "${gtf_base}.reflat"

######################################################
# Generate rRNA BED and interval list
######################################################
awk 'OFS="\t" { if ($3 == "CDS") print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" \
  | tr -d '";' > "${gtf_base}.bed"

awk 'OFS="\t" { if ($2 == "rRNA" && $1 !~ /^H/) print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" \
  | tr -d '";' > "${gtf_base}.rRNA.bed"

ml gatk
gatk CreateSequenceDictionary -R "${gtf_dir}/Homo_sapiens.GRCh37.70.dna_sm.primary_assembly.fa"

module load picard/2.18.4
cd "$gtf_dir"
java -jar "$PICARD" BedToIntervalList \
  I="${gtf_base}.rRNA.bed" \
  O="${gtf_base}.rRNA.interval_list" \
  SD="$dict_file"

ln -s "${gtf_base}.rRNA.interval_list" "${gtf_base}.rRNA.interval.list"

######################################################
# Create MtRNA gene list and annotation
######################################################
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

######################################################
# Create globin gene annotation
######################################################
cp "$globin_gene_list" ./

cat globin_gene | xargs -n1 -I% grep % "$gtf_file" \
  | awk 'OFS="\t" { if ($3 == "CDS") { if ($1 !~ /^H/) print $1, $4-1, $5, $10, 0, $7; else print $1, $4-1, $5, $10, 0, $7 } }' \
  | tr -d '";' | sort | uniq > globin_gene.bed

bedtools sort -i globin_gene.bed > globin_gene.sort.bed
bedtools merge -s -c 4,5,6 -delim "|" -o distinct -i globin_gene.sort.bed > globin_gene.merge.bed

bedToGenePred globin_gene.merge.bed globin_gene.merge.genePred
genePredToGtf file globin_gene.merge.genePred globin_gene.merge.gtf
awk '$3 == "exon"' globin_gene.merge.gtf > globin_gene.gtf

######################################################
# Cleanup intermediate files
######################################################
rm MtRNA_gene.bed MtRNA_gene.merge.* globin_gene.bed globin_gene.sort.bed globin_gene.merge.*
