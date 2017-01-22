#!/bin/bash

output_dir=${HOME}/Sequences/run
start_time=$(date +"%Y-%m-%d_%H_%M_%S")
merged_bam=${output_dir}/bam_merged/filtered
merged_bam_list=($(ls ${merged_bam}/*.bam))
log_dir=${output_dir}/log


temp_bam_folder=${merged_bam}/de_duped_bam_tmp
picard_metrics_dir=${temp_bam_folder}/picard_metrics

command_log=${log_dir}/de-duped_${start_time}.txt
touch ${command_log}
for file in ${merged_bam_list[*]}; do
    sample_name=$(basename $file .bam)
    mark_output=${sample_name}_remove_duplicates.bam

    command="java -jar ${output_dir}/../tools/picard-tools-1.129/picard.jar MarkDuplicates INPUT=${file} OUTPUT=${temp_bam_folder}/${mark_output} METRICS_FILE=${picard_metrics_dir}/$(basename ${file}).txt REMOVE_DUPLICATES=TRUE"
    echo "Triggering $command" | tee -a ${command_log}

    time $command &>>${command_log}
done
