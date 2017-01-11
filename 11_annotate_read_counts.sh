#!/bin/bash

output_dir=${HOME}/Sequences/run

merge_dir=${output_dir}/merge_pvalue
read_count=${merge_dir}/read_count
annotated=${merge_dir}/annotated

sorted=${annotated}/sorted
partial=${annotated}/part
combined=${annotated}/combined

mkdir ${read_count}
mkdir ${annotated}
mkdir ${sorted}
mkdir ${partial}
mkdir ${combined}

for file in $(ls ${merge_dir}/*_0.01_peaks_extended_250.bed); do
    annotatePeaks.pl $file mm9 > ${annotated}/$(basename $file);
done

for file in $(ls ${read_count}/*.txt); do
    name=$(basename $file .txt);
    sed 1d ${annotated}/${name}_0.01_peaks_extended_250.bed > ${partial}/${name}.bed;
done

for file in $(ls ${partial}/*.bed); do
    sort -k2,2 -k3,3n -k4,4n $file > ${sorted}/$file;
done

for file in $(ls ${read_count}/*.txt); do 
    name=$(basename $file .txt);
    awk -v filename=$file '{ getline v < filename ; split( v, a ); print a[2],"\t"$0 }' ${sorted}/${name}.bed > ${combined}/${name}.bed;
done