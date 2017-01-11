#!/bin/bash

output_dir=${HOME}/Sequences/run
log_dir=${output_dir}/log

bed_folder=${output_dir}/macs_pvalue
bed_list=($(ls ${bed_folder}/*.narrowPeak))

pvalue_clean_folder=${bed_folder}/clean
mkdir ${pvalue_clean_folder}

blacklist=${output_dir}/../mouse/blacklist/mm9-blacklist.bed

# 1) remove blacklisted positions from narrowPeak files
start_time=$(date +"%Y-%m-%d_%H_%M_%S")
command_log=${log_dir}/remove_blacklist_narrowPeak_${start_time}.txt
touch ${command_log}
for file in ${bed_list[*]}; do
    filename=$(basename ${file} .narrowPeak)
    edited=${filename}_narrowPeak_clean.bed
    echo "Triggering 'bedtools intersect -a ${file} -b ${blacklist} -v > ${pvalue_clean_folder}/${edited}'" | tee -a ${command_log}
    time bedtools intersect -a ${file} -b ${blacklist} -v > ${pvalue_clean_folder}/${edited}
done

bed_list=($(ls ${bed_folder}/*_summits.bed))
# 2) remove blacklisted positions from summits files
command_log=${log_dir}/remove_blacklist_summits_${start_time}.txt
touch ${command_log}
for file in ${bed_list[*]}; do
    filename=$(basename ${file} _summits.bed)
    edited=${filename}_summits_clean.bed
    echo "Triggering 'bedtools intersect -a ${file} -b ${blacklist} -v > ${pvalue_clean_folder}/${edited}'" | tee ${command_log}
    time bedtools intersect -a ${file} -b ${blacklist} -v > ${pvalue_clean_folder}/${edited}
done
# 2) Remove blacklisted position in xls output files
#    Also removes redundant header lines
xls_list=($(ls ${bed_folder}/*.xls))

command_log=${log_dir}/remove_blacklist_peaks_and_summits_xls_${start_time}.txt
touch ${command_log}s
for file in ${xls_list[*]}; do
    filename=$(basename ${file} .xls)
    # sed breakdown:
    # The following sed command deletes everything that starts with # or is empty line
    #     x -> exchange pattern and hold buffer
    #     /^#\|^$/ -> search for lines beginning with # symbol
    #              -> or empty lines
    #     !g -> Copy hold space into pattern buffer
    #     !p -> don't print anything in current pattern space (e.g, delete everything that is in hold space)
    sed -n 'x;/^#\|^$/!g;//!p' ${file} > ${pvalue_clean_folder}/${filename}_edited.xls
    
    echo "Triggering 'bedtools intersect -a ${pvalue_clean_folder}/${filename}_edited.xls -b ${blacklist} -v > ${pvalue_clean_folder}/${filename}.xls'" | tee -a ${command_log}
    time bedtools intersect -a ${pvalue_clean_folder}/${filename}_edited.xls -b ${blacklist} -v > ${pvalue_clean_folder}/${filename}.xls
    rm ${pvalue_clean_folder}/${filename}_edited.xls
done

