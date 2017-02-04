#!/bin/bash

##
# Remove duplications from all bam files in the input directory
#
# Arguments
# $1 Input directory
function remove_duplicates()
{
    input_dir=''
    if [ -z "$1" ] ; then
        input_dir="/$1"
    fi
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local merged_bam="${RUN_DIR}/bam_merged${input_dir}"
    local temp_bam_folder="${merged_bam}/de_duped_bam_tmp"
    local picard_metrics_dir="${temp_bam_folder}/picard_metrics"
    local command_log="${LOG_DIR}/de-duped_${start_time}.txt"
    local files=($(ls ${merged_bam}/*.bam))

    touch ${command_log}
    for file in ${files[*]}; do
        local sample_name=$(basename $file .bam)
        local mark_output=${sample_name}_remove_duplicates.bam
        local output="OUTPUT=${temp_bam_folder}/${mark_output}"
        local metrics="METRICS_FILE=${picard_metrics_dir}/$(basename ${file}).txt"

        local com="java -jar ${PICARD_LOCATION} MarkDuplicates INPUT=${file} ${output} $metrics REMOVE_DUPLICATES=TRUE"
        echo "Triggering $com" | tee -a ${command_log}
        time $com &>>${command_log}
    done
}
