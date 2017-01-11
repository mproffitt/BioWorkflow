#!/bin/bash

output_dir=${HOME}/Sequences/run
sam_dir=${output_dir}/sam
files=($(ls ${sam_dir}/*.sam))
bam_dir=${output_dir}/bam_lanewise
tmp_dir=${output_dir}/tmp
log_dir=${output_dir}/log
start_time=$(date +"%Y-%m-%d_%H_%M_%S")

command_log=${log_dir}/sam_to_bam_log_${start_time}.txt
touch ${command_log}

for file in ${files[*]}; do
    sample_name=$(basename ${file} _.sam)
    output=${bam_dir}/${sample_name}.bam
    echo "Triggering 'samtools view -bS -q 10 -F 260 ${file} -o ${output}'" | tee -a ${command_log}
    time samtools view -bS -q 10 -F 260 ${file} -o ${output} &>>${command_log}

# Mapping statistics
#    output=${tmp_dir}/${sample_name}.bam
#    echo "Triggering 'samtools view -bS -F 260 ${file} -o ${output}'" | tee -a ${command_log}
#    time samtools view -bS -F 260 ${file} -o ${output} &>>${command_log}
done

