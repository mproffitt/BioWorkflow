#!/bin/bash

function filter_chromasones()
{
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local merged_bam=${RUN_DIR}/bam_merged
    local merged_bam_list=($(ls ${merged_bam}/*.bam))
    local bam_filtered_dir="${merged_bam}/filtered"
    local command_log="${LOG_DIR}/filtered-${start_time}.log"

    touch ${command_log}
    for file in ${merged_bam_list[*]}; do
        local name=$(basename $file)
        local command="samtools view -b ${file} ${CHROM_LIST[@]} > ${bam_filtered_dir}/${name}"
        echo "Triggering '${command}'" | tee -a ${command_log}
        time $command | tee ${command_log}

        command="samtools index ${bam_filtered_dir}/${name}"
        echo "Triggering '${command}'" | tee -a ${command_log}
        time $command
    done
}

