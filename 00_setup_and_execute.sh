#!/bin/bash

run_dir=${HOME}/Sequences/run

# Sets the backup name to the earliest modified date of files under the run directory
[ -d "${run_dir}" ] && run_backup='run-'$(date -d "$(ls -al ${run_dir} | awk '{print $7,$6,2016}' | tail -n +2 | sort | uniq | head -1)" +"%Y-%m-%d") 

# Set up the run directory structure
[ -d ${run_dir} ] && mv ${run_dir} ${HOME}/Sequences/${run_backup}
mkdir ${run_dir}

mkdir ${HOME}/Sequences/run/alignment_stats
mkdir ${HOME}/Sequences/run/bam_lanewise
mkdir ${HOME}/Sequences/run/bam_merged
mkdir ${HOME}/Sequences/run/bam_merged/de_duped_bam_tmp
mkdir ${HOME}/Sequences/run/bam_merged/de_duped_bam_tmp/picard_metrics

mkdir ${HOME}/Sequences/run/blacklist
mkdir ${HOME}/Sequences/run/log
mkdir ${HOME}/Sequences/run/macs_broad
mkdir ${HOME}/Sequences/run/macs_pvalue
mkdir ${HOME}/Sequences/run/merge_broad
mkdir ${HOME}/Sequences/run/merge_broad/clean
mkdir ${HOME}/Sequences/run/merge_broad/extended

mkdir ${HOME}/Sequences/run/merge_pvalue
mkdir ${HOME}/Sequences/run/merge_pvalue/clean
mkdir ${HOME}/Sequences/run/merge_pvalue/extended

mkdir ${HOME}/Sequences/run/sam
mkdir ${HOME}/Sequences/run/scripts
mkdir ${HOME}/Sequences/run/tmp

