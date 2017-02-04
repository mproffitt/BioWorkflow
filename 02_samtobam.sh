#!/bin/bash

##
# Convert sam file format to BAM format
function sam_to_bam()
{
    local sam_dir=${RUN_DIR}/sam
    local files=($(ls ${sam_dir}/*.sam))
    local bam_dir=${RUN_DIR}/bam_lanewise
    local tmp_dir=${RUN_DIR}/tmp
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/sam_to_bam_log_${start_time}.txt"

    touch ${command_log}
    for file in ${files[*]}; do
        local sample_name=$(basename ${file} _.sam)
        local output="${bam_dir}/${sample_name}.bam"
        local command="samtools view -bS -q 10 -F 260 ${file} -o ${output}"
        echo "Triggering '${command}'" | tee -a ${command_log}
        time $command | tee -a ${command_log}
    done
}
