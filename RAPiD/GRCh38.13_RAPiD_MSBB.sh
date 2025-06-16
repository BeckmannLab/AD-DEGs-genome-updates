#!/bin/bash

genome="hg38"
transcriptome="v43"
assembly="GRCh38.Gencode.v43"

batch_size=100

########################################
# Link FASTQ files
########################################
while read -r NAME; do
  sample_dir="${run_dir}/${NAME}"
  mkdir -p "$sample_dir"
  ln -sf "${fastq_dir}/${NAME}.sorted.fastq.gz" "${sample_dir}/${NAME}.fastq.gz"
done < "$ids_file"

########################################
# Split into batches and generate scripts
########################################
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

  script_name="rapid_run_${genome}${transcriptome}_batch${batchID}.sh"
  cat <<EOF > "$script_name"
#!/bin/bash

scratch_rapid1="${scratch_dir}"
mkdir -p "\$scratch_rapid1"
run_folder1="${batch_path}"
script_rapid="${script_path}"

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

########################################
# Copy logs and outputs from scratch
########################################
batchID=3  # adjust as needed
batch_dest="${run_dir}/batch${batchID}"
cp "${scratch_dir}/pipeline_trace.txt"* "$batch_dest"
cp "${scratch_dir}"/*std* "$batch_dest"
cp "${scratch_dir}/report.html" "$batch_dest"

########################################
# Clean up scratch space
########################################
rm -rf "${scratch_dir:?}"/*
cd "$run_dir"
