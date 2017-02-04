#!/bin/bash

function read_counts()
{
    local input_dir=''
    if [ -z "$1" ] ; then
        input_dir="/$1"
    fi
    local bam_dir=${RUN_DIR}/bam_merged${input_dir}
    local merge_dir=${RUN_DIR}/merge_pvalue

    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local read_count=${merge_dir}/read_count
    local command_log=${LOG_DIR}/read_counts_${start_time}.txt
    local files=($(ls ${merge_dir}/*.bed))

    touch ${command_log}
    for file in ${bed_list[*]}; do
        if ! echo $file | grep -q peaks ; then
            continue
        fi
        local input_file=$(basename $file _remove_duplicates_${PRECISION}_peaks_extended_${EXTEND_RANGE}.bed)
        if [ -f ${bam_dir}/${input_file}.bam ] && [ ! -f "${read_count}/${input_file}.txt" ] ; then
            echo "Creating read counts for ${input_file}" | tee -a ${command_log}
            for peak_value in $(sed 's/\t/:/;s/\t/-/;' ${merge_dir}/$(basename $file)); do
                echo -n ${peak_value}$'\t'
                samtools view ${bam_dir}/${input_file}.bam ${peak_value} | wc -l
            done > ${read_count}/${input_file}.txt
        fi
    done
}
