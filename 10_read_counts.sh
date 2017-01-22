#!/bin/bash

output_dir=${HOME}/Sequences/run
bam_dir=${output_dir}/bam_merged/filtered

merge_dir=${output_dir}/merge_pvalue
bed_list=($(ls ${merge_dir}/*.bed))

read_count=${merge_dir}/read_count
mkdir -p ${read_count}

for file in ${bed_list[*]}; do
    if ! echo $file | grep -q peaks ; then
        continue
    fi
    input_file=$(basename $file _remove_duplicates_0.01_peaks_extended_250.bed)
    if [ -f ${bam_dir}/${input_file}.bam ] && [ ! -f "${read_count}/${input_file}.txt" ] ; then
        echo "Creating read counts for ${input_file}"
        for peak_value in $(sed 's/\t/:/;s/\t/-/;' ${merge_dir}/$(basename $file)); do
	    echo -n ${peak_value}$'\t'
            samtools view ${bam_dir}/${input_file}.bam ${peak_value} | wc -l
        done > ${read_count}/${input_file}.txt
    fi
done

