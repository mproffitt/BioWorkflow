#!/bin/bash

function _extend()
{
    local peak_type=$1
    local bed_dir=${RUN_DIR}/macs_pvalue/clean
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log=${LOG_DIR}/extend_${peak_type}_${start_time}.txt
    local known_sizes=${RUN_DIR}/../${GENOME}/${STRAIN}/${STRAIN}.chrom.sizes
    local extended=""

    touch ${command_log}
    local files=($(ls ${bed_dir}/*_${peak_type}_clean.bed))
    for file in ${files[*]}; do
        local filename=$(basename $file _${peak_type}_clean.bed)
        local edited=${filename}_extended_${EXTEND_RANGE}.bed

        local command="bedtools slop -i ${file} -b ${EXTEND_RANGE} -g ${known_sizes} |"
        command="${command} bedtools merge -i - > ${extended}/${edited}"

        echo "Triggering '${command}'" | tee -a ${command_log}
        time $command | tee -a ${command_log}
    done
}

function extend_pvalue()
{
    _extend "narrowPeak"
    _extend "summits"
}

