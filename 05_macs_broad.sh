#!/bin/bash

output_dir=${HOME}/Sequences/run
start_time=$(date +"%Y-%m-%d_%H_%M_%S")
log_dir=${output_dir}/log

bam_merged_dir=${output_dir}/bam_merged
de_duped_dir=${bam_merged_dir}/de_duped_bam_tmp

macs_broad_dir=${output_dir}/macs_broad

input=$(ls $de_duped_dir | grep -qi input)
files=($(ls ${de_duped_dir}/*.bam))

command_log=${log_dir}/macs_broad_${start_time}.txt
touch ${command_log}
for file in ${files[*]}; do
    if echo $file | grep -qi input; then
        continue
    fi
    name=$(basename $file .bam)
    echo "Triggering 'macs2 callpeak -t ${file} -c ${input} --broad --format BAM -g mm -n ${name}_0.001 --outdir ${macs_broad_dir} --broad-cutoff 0.001'" | tee -a ${command_log}
    time macs2 callpeak -t ${file} -c ${input} --broad --format BAM -g mm -n ${name}_0.01 --outdir ${macs_broad_dir} --broad-cutoff 0.001 &>>${command_log}
done
