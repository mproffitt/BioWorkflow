#!/bin/bash

output_dir=${HOME}/Sequences/run
log_dir=${output_dir}/log

bed_dir=${output_dir}/macs_pvalue/clean

start_time=$(date +"%Y-%m-%d_%H_%M_%S")
command_log=${log_dir}/extend_narrowPeak_${start_time}.txt
touch ${command_log}

range=250

known_sizes=${output_dir}/../mouse/mm10/mm10.chrom.sizes

extended=${output_dir}/macs_pvalue/extended
mkdir $extended

bed_list=($(ls ${bed_dir}/*_narrowPeak_clean.bed))
for file in ${bed_list[*]}; do
    filename=$(basename $file _narrowPeak_clean.bed)
    edited=${filename}_extended_${range}.bed
    echo "Triggering 'bedtools slop -i ${file} -b ${range} -g ${known_sizes} | bedtools merge -i - > ${extended}/${edited}'" | tee -a ${command_log}
    time bedtools slop -i ${file} -b ${range} -g ${known_sizes} | bedtools merge -i - > ${extended}/${edited}
done

command_log=${log_dir}/extend_summits_${start_time}.txt
bed_list=($(ls ${bed_dir}/*_summits_clean.bed))
for file in ${bed_list[*]}; do
    filename=$(basename $file _summits_clean.bed)
    edited=${filename}_extended_${range}.bed
    echo "Triggering 'bedtools slop -i ${file} -b ${range} -g ${known_sizes} | bedtools merge -i - > ${extended}/${edited}'" | tee -a ${command_log}
    time bedtools slop -i ${file} -b ${range} -g ${known_sizes} | bedtools merge -i - > ${extended}/${edited}
done
