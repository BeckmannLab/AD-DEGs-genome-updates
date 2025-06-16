#!/bin/bash

######################################################
# Define all path and file variables
######################################################
genome=hg19
assembly=GRCh37.Ensembl.v70

msbb_dir="${base_dir}/msbb"
fastq_dir="${msbb_dir}/fastq"
rapid_dir="${msbb_dir}/rapid_run/${genome}"
run_dir="${rapid_dir}/run"
ids_file="${msbb_dir}/all_ids.txt"

batch_size=100

######################################################
# Link FASTQ files by sample
######################################################
while read -r NAME; do
  sample_dir="${run_dir}/${NAME}"
  mkdir -p "$sample_dir"
  ln -sf "${fastq_dir}/${NAME}.sorted.fastq.gz" "${sample_dir}/${NAME}.fastq.gz"
done < "$ids_file"

######################################################
# Split into batches and create RAPiD run scripts
######################################################
cd "$run_dir"
sample_dirs=($(ls -d */ | grep -v "batch"))
numbatches=$(( (${#sample_dirs[@]} + batch_size - 1) / batch_size ))

for batchID in $(seq 0 $((numbatches - 1))); do
  batch_name="batch${batchID}"
  batch_path="${run_dir}/${batch_name}"
  mkdir -p "$batch_path"

  start=$((batchID * batch_size))
  batch_samples=("${sample_dirs[@]:$start:$batch_size}")
  cp -r "${batch_samples[@]}" "$batch_path"
  rm -rf "$batch_path"/*/RAPiD/

  script_name="${run_dir}/rapid_run_${genome}_batch${batchID}.sh"
  cat <<EOF > "$script_name"
#BSUB -J rapid_run_${genome}
#BSUB -P 
#BSUB -q 
#BSUB -n 
#BSUB -W 
#BSUB -R 
#BSUB -o ${scratch_dir}/%J.stdout
#BSUB -eo ${scratch_dir}/%J.stderr
#BSUB -L /bin/bash

scratch_rapid1="${scratch_dir}"
mkdir -p "\$scratch_rapid1"
run_folder1="${batch_path}"
script_rapid="${pipeline_script}"

cd "\$scratch_rapid1"
module purge
module load java R

${nextflow_bin} run "\$script_rapid" \\
  --run "\$run_folder1" \\
  --singleEnd -profile RiboZero \\
  --genome "${assembly}" \\
  --stranded none \\
  --twopass1readsN 18446744073709551615 \\
  --twopassMode Basic \\
  --rawPath . \\
  --outPath RAPiD \\
  --fastqc --featureCounts --qc -resume
EOF
done

######################################################
# Submit selected batch jobs manually
######################################################
# Example:
# bsub < ${run_dir}/rapid_run_${genome}_batch0.sh
# bsub < ${run_dir}/rapid_run_${genome}_batch1.sh
# ...

######################################################
# Copy logs and output from scratch to project directory
######################################################
batchID=3  # Adjust manually
dest_dir="${run_dir}/batch${batchID}"
for f in "${trace_files[@]}"; do
  cp "${scratch_dir}/$f" "$dest_dir"
done

######################################################
# Final cleanup after job completion
######################################################
rm -rf "${scratch_dir:?}"/*
rm -rf "${run_dir}/${genome}"

cd "$scratch_dir" && rm -rf *
cd "${run_dir}"
