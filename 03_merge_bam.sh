#!/bin/bash

##
# Merges multiple bam files into a single bam, one per track.
#
# Arguments:
# This function takes a set of patterns/replacements as its arguments.
#
# For each group of lanes, a pattern and replacement must be provided
# in the form of "pattern^replacement"
#
# For example, if the output file is "AB123_Merged.bam" and the list of files constructing this merged is:
#     AB123_ENTR_1_RT_12.bam
#     AB123_ENTR_2_RT_12.bam
#     AB123_ENTR_3_RT_12.bam
#     AB123_ENTR_4_RT_12.bam
#
# then the argument matching this would be: "AB123_*_RT*^AB123_Merged"
#
function merge_bam_files()
{
    local arguments=($@)
    local patterns=(
        $(
            for pattern in ${arguments[@]}; do
                echo ${pattern} | cut -d^ -f1
            done
        )
    )

    local replacements=(
        $(
            for pattern in ${arguments[@]}; do
                echo ${pattern} | cut -d^ -f2
            done
        )
    )
    local bam_merged=${RUN_DIR}/bam_merged
    local input_dir=${RUN_DIR}/bam_lanewise
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/merge_bam_${start_time}.txt"

    touch command_log
    for ((i=0; i < ${#patterns[@]}; i++)); do
        local pattern="${patterns[$i]}"
        local replacement="${bam_merged}/${replacements[$i]}.bam"

        local files=($(ls $input_dir/${pattern}.bam))

        command="samtools merge - ${files[@]} | samtools sort -o ${bam_merged}/${replacement} -"

        echo "Triggering '$command'" | tee -a ${command_log}
        time $command

        echo "Triggering 'samtools index ${replacement}'" | tee -a ${command_log}
        time samtools index ${replacement}
    done
}
