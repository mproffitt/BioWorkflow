#!/bin/bash

# Script to execute Bowtie and align to mouse genome mm10 using bowtie2
#
bowtie_index=${HOME}/Sequences/mouse/mm9/mm9
fastq_dir=${HOME}'/Sequences/2016-07-18/Unaligned'
files=($(ls ${fastq_dir}/GL30*_RT*.gz))

output_dir=${HOME}/Sequences/run
sam_dir=${output_dir}/sam
log_dir=${output_dir}/log

start_time=$(date +"%Y-%m-%d_%H_%M_%S")

command_log=${log_dir}/bowtie2Commands_${start_time}.txt
touch ${command_log}

for file in ${files[*]}; do
    noext=${file%%_R1_001.fastq.gz}
    sample_name=$(basename ${noext})
    output=${sam_dir}/${sample_name}.sam

    command="bowtie2 -p 4 -q -x ${bowtie_index} -U ${file}  -S ${output}"
    echo "Triggering ${command}" | tee -a ${command_log}
    time bowtie2 -p 16 -q -x ${bowtie_index} -U ${file} -S ${output} &>>${command_log}
    sleep 5
done
