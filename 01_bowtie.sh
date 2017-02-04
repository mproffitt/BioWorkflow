#!/bin/bash

# Script to execute Bowtie and align to the requested genome using bowtie2
#
function run_bowtie()
{
    if [ -n "$1" ] ; then
        echo "No input pattern provided"
        exit 1
    fi
    local pattern=$1

    local files=($(ls ${FASTQ_DIR}/${pattern}))
    local sam_dir="${RUN_DIR}/sam"
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/bowtie2Commands_${start_time}.txt"

    touch ${command_log}
    for file in ${files[*]}; do
        local noext=${file%%_R1_001.fastq.gz}
        local sample_name=$(basename ${noext})
        local output=${sam_dir}/${sample_name}.sam

        command="bowtie2 -p ${PROCESSORS} -q -x ${BOWTIE_INDEX} -U ${file}  -S ${output}"
        echo "Triggering ${command}" | tee -a ${command_log}
        time $command | tee -a ${command_log}
        sleep 5
    done
}
