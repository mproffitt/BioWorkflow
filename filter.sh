output_dir=${HOME}/Sequences/run
start_time=$(date +"%Y-%m-%d_%H_%M_%S")
merged_bam=${output_dir}/bam_merged
merged_bam_list=($(ls ${merged_bam}/*.bam))
log_dir=${output_dir}/log
bam_filtered_dir="${merged_bam}/filtered"
command_log="${log_dir}/filtered-${start_time}.log"
touch ${command_log}

for file in ${merged_bam_list[*]}; do
    name=$(basename $file)
    echo "samtools view -b ${file} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ${bam_filtered_dir}/${name}" > ${command_log}
    samtools view -b ${file} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ${bam_filtered_dir}/${name}
    
    samtools index "${bam_filtered_dir}/${name}"
done