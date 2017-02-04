#!/bin/bash

function merge_pvalue()
{
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log=${LOG_DIR}/merge_broad_${start_time}.txt
    local extended_bed_dir=${RUN_DIR}/macs_pvalue/extended
    local known_sizes=${LOCATION}/${GENOME}/${STRAIN}/${STRAIN}.chrom.sizes
    local merge_dir=${RUN_DIR}/merge_pvalue
    local files=($(ls ${extended_bed_dir}/*.bed))

    touch ${command_log}
    for file in ${files[*]}; do
        local command="bedtools merge -d ${MERGE_RANGE} -i ${file} > ${merge_dir}/$(basename $file)"
        echo "Triggering '${command}'" | tee -a ${command_log}
        time $command | tee -a ${command_log}
    done
}

