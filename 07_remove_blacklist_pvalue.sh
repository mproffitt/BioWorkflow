#!/bin/bash

function _blacklist()
{
    local peak_type=$1
    local bed_folder=${RUN_DIR}/macs_pvalue
    local pvalue_clean_folder=${bed_folder}/clean

    local blacklist=${RUN_DIR}/../${GENOME}/blacklist/${STRAIN}-blacklist.bed

    local bed_list=($(ls ${bed_folder}/*.${peak_type}))
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")

    local command_log=${LOG_DIR}/remove_blacklist_${peak_type}_${start_time}.txt
    touch ${command_log}
    for file in ${bed_list[*]}; do
        local filename=$(basename ${file} .${peak_type})
        local edited=${filename}_${peak_type}_clean.bed

        command="bedtools intersect -a ${file} -b ${blacklist} -v > ${pvalue_clean_folder}/${edited}"
        echo "Triggering '$command'" | tee -a ${command_log}
        time $command
    done
}

function remove_blacklist()
{
    _blacklist narrowPeak
    _blacklist summits

    local bed_folder=${RUN_DIR}/macs_pvalue
    local xls_list=($(ls ${bed_folder}/*.xls))

    command_log=${LOG_DIR}/remove_blacklist_peaks_and_summits_xls_${start_time}.txt
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

        local command="bedtools intersect -a ${pvalue_clean_folder}/${filename}_edited.xls"
        command="$comamnd -b ${blacklist} -v > ${pvalue_clean_folder}/${filename}.xls"

        echo "Triggering '${command}'" | tee -a ${command_log}
        time $command && rm ${pvalue_clean_folder}/${filename}_edited.xls
    done
}

