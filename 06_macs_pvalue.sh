#!/bin/bash

function macs_pvalue()
{
    local precision=$1
    local input_dir=''
    if [ -z "$2" ] ; then
        input_dir="/$2"
    fi
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log=${LOG_DIR}/macs_pvalue_${start_time}.txt
    local bam_merged_dir=${RUN_DIR}/bam_merged${input_dir}
    local de_duped_dir=${bam_merged_dir}/de_duped_bam_tmp
    local macs_pvalue_dir=${RUN_DIR}/macs_pvalue
    local files=($(ls ${de_duped_dir}/*.bam))

    touch ${command_log}
    for file in ${files[*]}; do
        if echo $file | grep -qi input; then
            continue
        fi
        local input="${de_duped_dir}/$(ls $de_duped_dir | grep -i input | grep -i $(basename $file | cut -d_ -f1))"
        local name=$(basename $file .bam)
        local command="macs2 callpeak -t ${file} -c ${input} --format BAM -g mm -n"
        command="$command ${name}_${precision} --outdir ${macs_pvalue_dir} -p ${precision}'"

        echo "Triggering '$command'" | tee -a ${command_log}
        time $command | tee -a ${command_log}
    done
}
