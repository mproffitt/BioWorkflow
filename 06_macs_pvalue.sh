#!/bin/bash

output_dir=${HOME}/Sequences/run
log_dir=${output_dir}/log
start_time=$(date +"%Y-%m-%d_%H_%M_%S")
command_log=${log_dir}/macs_pvalue_${start_time}.txt
touch ${command_log}

bam_merged_dir=${output_dir}/bam_merged
de_duped_dir=${bam_merged_dir}/de_duped_bam_tmp

macs_pvalue_dir=${output_dir}/macs_pvalue

input=$(ls $de_duped_dir | grep -qi input)
files=($(ls ${de_duped_dir}/*.bam))
for file in ${files[*]}; do
    if echo $file | grep -qi input; then
        continue
    fi
    name=$(basename $file .bam)
    echo "Triggering 'macs2 callpeak -t ${file} -c ${input} --format BAM -g mm -n ${name}_0.01 --outdir ${macs_pvalue_dir} -p 0.01'" | tee -a ${command_log}
    time macs2 callpeak -t ${file} -c ${input} --format BAM -g mm -n ${name}_0.01 --outdir ${macs_pvalue_dir} -p 0.01 &>>${command_log}
done
