
#######################################
# Download Ensembl release 53 GTF and FASTA
#######################################
cd "$base_dir"
wget http://ftp.ensembl.org/pub/release-53/gtf/homo_sapiens/Homo_sapiens.NCBI36.53.gtf.gz
wget http://ftp.ensembl.org/pub/release-53/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.53.dna.toplevel.fa.gz


#######################################
# Create STAR genome reference
#######################################
ml star/2.7.3a
mkdir -p "$base_dir/chr_primary"

STAR --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir "$base_dir/chr_primary" \
  --genomeFastaFiles "$fasta_file" \
  --sjdbGTFfile "$gtf_file" \
  --sjdbOverhang 100 \
  --outFileNamePrefix GRCh37.ensembl.v53

#######################################
# Generate GenePred and Reflat files
#######################################
cd "$base_dir"
"$gtftogenepred" -genePredExt "$gtf_file" "$output_name.GenePred"

awk '($2 !~ /^H/){print $12 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' \
  "$output_name.GenePred" > "$output_name.reflat"

#######################################
# Create rRNA BED and interval list
#######################################
awk -F'\t' 'OFS="\t" { if ($3 == "CDS") print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" | tr -d '";' > "$output_name.bed"
awk -F'\t' 'OFS="\t" { print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" | tr -d '";' > "check_$output_name.bed"
awk -F'\t' 'OFS="\t" { if ($2 == "rRNA" && $1 !~ /^H/) print $1, $4-1, $5, $10, $2, $7 }' "$gtf_file" | tr -d '";' > "$output_name.rRNA.bed"

ml gatk
gatk CreateSequenceDictionary -R "$fasta_file"

module load picard/2.18.4
java -jar "$PICARD" BedToIntervalList \
  I="$output_name.rRNA.bed" \
  O="$output_name.rRNA.interval_list" \
  SD="$dict_file"

ln -s "$output_name.rRNA.interval_list" "$output_name.rRNA.interval.list"

#######################################
# Extract MtRNA and globin genes
#######################################
grep "^MT" "$gtf_file" | awk 'OFS="\t" { if ($3 == "CDS") print $16 }' | tr -d '";' > MtRNA_gene
grep "^MT" "$gtf_file" | awk 'OFS="\t" { if ($3 == "CDS") print $1, $4-1, $5, $10, 0, $7 }' | tr -d '";' > MtRNA_gene.bed

ml bedtools
bedtools sort -i MtRNA_gene.bed | \
  bedtools merge -s -c 4,5,6 -delim "|" -o distinct > MtRNA_gene.merge.bed

ml ucsc-utils/2020-03-17
bedToGenePred MtRNA_gene.merge.bed MtRNA_gene.merge.genePred
genePredToGtf file MtRNA_gene.merge.genePred MtRNA_gene.merge.gtf
awk '$3=="exon"' MtRNA_gene.merge.gtf > MtRNA_gene.gtf

cp /sc/arion/projects/H_PBG/REFERENCES/GRCh38/Gencode/release_30/globin_gene ./

cat globin_gene | xargs -n1 -I% grep % "$gtf_file" | \
  awk 'OFS="\t" { if ($3 == "CDS" && $1 !~ /^H/) print $1, $4-1, $5, $10, 0, $7 }' | \
  tr -d '";' | sort -u > globin_gene.bed

bedtools sort -i globin_gene.bed > globin_gene.sort.bed
bedtools merge -s -c 4,5,6 -delim "|" -o distinct -i globin_gene.sort.bed > globin_gene.merge.bed
bedToGenePred globin_gene.merge.bed globin_gene.merge.genePred
genePredToGtf file globin_gene.merge.genePred globin_gene.merge.gtf
awk '$3=="exon"' globin_gene.merge.gtf > globin_gene.gtf

#######################################
# Clean up intermediate BED/GTF files
#######################################
rm -f MtRNA_gene.bed MtRNA_gene.merge.* globin_gene.bed globin_gene.sort.bed globin_gene.merge.*
