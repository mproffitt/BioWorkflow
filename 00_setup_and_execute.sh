#!/bin/bash
export XDG_COFNIG_HOME=~/.config

run_dir=${HOME}/Sequences/run

[ -d ~/pyccata/src/images ] && mv ~/pyccata/src/images ${run_dir}/overlap_output
[ ! -d ~/pyccata/src/images ] && mkdir ~/pyccata/src/images

# Sets the backup name to the earliest modified date of files under the run directory
[ -d "${run_dir}" ] && run_backup='run-'$(date -d "$(ls -al ${run_dir} | awk '{print $7,$6,2017}' | tail -n +2 | sort | uniq | head -1)" +"%Y-%m-%d") 

# Set up the run directory structure
[ -d ${run_dir} ] && mv ${run_dir} ${HOME}/Sequences/${run_backup}
mkdir ${run_dir}

mkdir ${HOME}/Sequences/run/alignment_stats
mkdir ${HOME}/Sequences/run/bam_lanewise
mkdir ${HOME}/Sequences/run/bam_merged
mkdir ${HOME}/Sequences/run/bam_merged/de_duped_bam_tmp
mkdir ${HOME}/Sequences/run/bam_merged/de_duped_bam_tmp/picard_metrics

mkdir ${HOME}/Sequences/run/bam_merged/filtered
mkdir ${HOME}/Sequences/run/bam_merged/filtered/de_duped_bam_tmp
mkdir ${HOME}/Sequences/run/bam_merged/filtered/de_duped_bam_tmp/picard_metrics

mkdir ${HOME}/Sequences/run/blacklist
mkdir ${HOME}/Sequences/run/log
mkdir ${HOME}/Sequences/run/macs_broad
mkdir ${HOME}/Sequences/run/macs_pvalue
mkdir ${HOME}/Sequences/run/macs_pvalue/clean
mkdir ${HOME}/Sequences/run/merge_broad
mkdir ${HOME}/Sequences/run/merge_broad/clean
mkdir ${HOME}/Sequences/run/merge_broad/extended

mkdir ${HOME}/Sequences/run/merge_pvalue
mkdir ${HOME}/Sequences/run/merge_pvalue/clean
mkdir ${HOME}/Sequences/run/merge_pvalue/extended
mkdir ${HOME}/Sequences/run/merge_pvalue/read_count
mkdir ${HOME}/Sequences/run/merge_pvalue/annotated
mkdir ${HOME}/Sequences/run/merge_pvalue/annotated/sorted
mkdir ${HOME}/Sequences/run/merge_pvalue/annotated/part
mkdir ${HOME}/Sequences/run/merge_pvalue/annotated/combined

mkdir ${HOME}/Sequences/run/sam
mkdir ${HOME}/Sequences/run/scripts
mkdir ${HOME}/Sequences/run/tmp

cp ${HOME}/scripts/* ${HOME}/Sequences/run/scripts/

rm -rf ${run_dir}/bam_merged
cp -r ${HOME}/Sequences/${run_backup}/bam_merged ${run_dir}/

rm -f ${run_dir}/bam_merged/de_duped_bam_tmp/*.bam
rm -f ${run_dir}/bam_merged/de_duped_bam_tmp/picard_metrics/*.txt
rm -f ${run_dir}/bam_merged/de_duped_bam_tmp/filtered/*.bam
rm -f ${run_dir}/bam_merged/de_duped_bam_tmp/filtered/picard_metrics/*.txt

cd ${HOME}

./scripts/04_remove_duplicates.sh &&
    ./scripts/06_macs_pvalue.sh && 
    ./scripts/07_remove_blacklist_pvalue.sh && 
    ./scripts/08_extend_pvalue.sh && 
    ./scripts/09_merge_pvalue.sh && 
    ./scripts/10_read_counts.sh && 
    ./scripts/11_annotate_read_counts.sh

if [ $? -eq 0 ] ; then
    cd ~/pyccata/src
    [ ! -d images ] && mkdir images
    time python3 annotation.py
    mv images ${run_dir}/overlap_output
fi
