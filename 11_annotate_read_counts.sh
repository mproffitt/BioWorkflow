#!/bin/bash


function annotate_read_counts()
{
    local merge_dir=${RUN_DIR}/merge_pvalue
    local read_count=${merge_dir}/read_count
    local annotated=${merge_dir}/annotated

    local sorted=${annotated}/sorted
    local partial=${annotated}/part
    local combined=${annotated}/combined
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/annotate_read_counts_${start_time}.log"

    touch ${command_log}
    for file in $(ls ${merge_dir}/*_${PRECISION}_peaks_extended_${EXTEND_RANGE}.bed); do
        local command="annotatePeaks.pl $file ${STRAIN} > ${annotated}/$(basename $file)";
        echo "Triggering '${command}'" | tee -a ${command_log}
        time $command | tee -a ${command_log}
    done

    for file in $(ls ${read_count}/*.txt); do
        local sname=$(basename $file .txt)
        local name=${sname}_remove_duplicates_${PRECISION}_peaks_extended_${EXTEND_RANGE}.bed

        sed 1d ${annotated}/${name} | sort -k2,2 -k3,3n -k4,4n $file > ${sorted}/${sname}
        awk -v filename=$file '{ getline v < filename ; split( v, a ); print a[2],"\t"$0 }' ${sorted}/${sname}.bed > ${combined}/${sname}.bed
    done
}
