#!/bin/bash

output_dir=${HOME}/Sequences/run
log_dir=${output_dir}/log

start_time=$(date +"%Y-%m-%d_%H_%M_%S")
command_log=${log_dir}/merge_broad_${start_time}.txt
touch ${command_log}

extended_bed_dir=${output_dir}/macs_broad/extended
bed_list=($(ls ${extended_bed_dir}/*.bed))
range=100
known_sizes=${output_dir}/../mouse/mm10.chrom.sizes
merge_dir=${output_dir}/merge_broad

for file in ${bed_list[*]}; do
    echo "Triggering 'bedtools merge -d ${range} -i ${file} > ${merge_dir}/$(basename $file)'" | tee -a ${command_log}
    time bedtools merge -d ${range} -i ${file} > ${merge_dir}/$(basename $file)
done
