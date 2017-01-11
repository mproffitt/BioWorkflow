#!/bin/bash

output_dir=${HOME}/Sequences/run
log_dir=${output_dir}/log

start_time=$(date +"%Y-%m-%d_%H_%M_%S")
command_log=${log_dir}/remove_blacklist_broad_${start_time}.txt
touch ${command_log}

bed_folder=${output_dir}/macs_broad
bed_list=($(ls ${bed_folder}/*.broadPeak))

broad_clean_folder=${bed_folder}/clean
mkdir $broad_clean_folder

blacklist=${output_dir}/../mouse/blacklist/mm9-blacklist.bed

# 1) remove blacklisted positions from broadPeak files
for file in ${bed_list[*]}; do
    filename=$(basename ${file} .broadPeak)
    edited=${filename}_broadPeak_clean.bed
    echo "Triggering 'bedtools intersect -a ${file} -b ${blacklist} -v > ${broad_clean_folder}/${edited}'" | tee -a ${command_log}
    time bedtools intersect -a ${file} -b ${blacklist} -v > ${broad_clean_folder}/${edited}
done

# 2) Remove blacklisted position in xls output files
#    Also removes redundant header lines
xls_list=($(ls ${bed_folder}/*.xls))
command_log=${log_dir}/remove_blacklist_broad_xls_${start_time}.txt
touch ${command_log}
for file in ${xls_list[*]}; do
    filename=$(basename ${file} .xls)
    # sed breakdown:
    # The following sed command deletes everything that starts with # or is empty line
    #     x -> exchange pattern and hold buffer
    #     /^#\|^$/ -> search for lines beginning with # symbol
    #              -> or empty lines
    #     !g -> Copy hold space into pattern buffer
    #     !p -> don't print anything in current pattern space (e.g, delete everything that is in hold space)
    sed -n 'x;/^#\|^$/!g;//!p' ${file} > ${broad_clean_folder}/${filename}_edited.xls
    echo "Triggering 'bedtools intersect -a ${broad_clean_folder}/${filename}_edited.xls -b ${blacklist} -v > ${broad_clean_folder}/${filename}.xls'" | tee -a ${command_log}
    time bedtools intersect -a ${broad_clean_folder}/${filename}_edited.xls -b ${blacklist} -v > ${broad_clean_folder}/${filename}.xls
    rm ${broad_clean_folder}/${filename}_edited.xls
done

