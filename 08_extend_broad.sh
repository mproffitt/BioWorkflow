#!/bin/bash

output_dir=${HOME}/Sequences/run
log_dir=${output_dir}/log

bed_dir=${output_dir}/macs_broad/clean

start_time=$(date +"%Y-%m-%d_%H_%M_%S")
command_log=${log_dir}/extend_broadPeak_${start_time}.txt
touch ${command_log}

bed_list=($(ls ${bed_dir}/*_broadPeak_clean.bed))

range=250

known_sizes=${output_dir}/../mouse/mm9/mm9.chrom.sizes

extended=${output_dir}/macs_broad/extended
mkdir $extended

for file in ${bed_list[*]}; do
    filename=$(basename $file _broadPeak_clean.bed)
    edited=${filename}_extended_${range}.bed
    echo "Triggering 'bedtools slop -i ${file} -b ${range} -g ${known_sizes} | bedtools merge -i - > ${extended}/${edited}'" | tee -a ${command_log}
    time bedtools slop -i ${file} -b ${range} -g ${known_sizes} | bedtools merge -i - > ${extended}/${edited}
done
