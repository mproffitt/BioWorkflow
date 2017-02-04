#!/bin/bash

export GENOME='mouse'
export STRAIN='mm9'
export RUN_DATE='2017-02-04'
export PICARD_VERSION='2.5.0'
export XDG_CONFIG_HOME=~/.config
export EXTEND_RANGE=250
export MERGE_RANGE=100

export LOCATION="${HOME}/Sequences"

export BOWTIE_INDEX="${LOCATION}/${GENOME}/${STRAIN}"
export FASTQ_DIR="${LOCATION}/${RUN_DATE}/Unaligned"

export RUN_DIR=${LOCATION}/run

export LOG_DIR=${RUN_DIR}/log
export PICARD_LOCATION="${RUN_DIR}/../tools/picard-tools-${PICARD_VERSION}/picard.jar"
export PROCESSORS=$(grep processor /proc/cpuinfo | wc -l)
export CHROM_LIST=(
    'chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9'
    'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16'
    'chr17' 'chr18' 'chr19' 'chrX' 'chrY'
)

if [ "${GENOME}" = 'human' ] ; then
    CHROM_LIST+=( 'chr20' 'chr21')
fi

pattern="GL30*_RT*.gz"

[ -d ~/pyccata/src/images ] && mv ~/pyccata/src/images ${RUN_DIR}/overlap_output
[ ! -d ~/pyccata/src/images ] && mkdir ~/pyccata/src/images

# Sets the backup name to the earliest modified date of files under the run directory
if [ -d "${RUN_DIR}" ]; then
    run_backup='run-'$(
        date -d "$(ls -al ${RUN_DIR} | awk '{print $7,$6,2017}' | tail -n +2 | sort | uniq | head -1)" +"%Y-%m-%d"
    )
    mv ${RUN_DIR} ${LOCATION}/${run_backup}
fi

# Set up the run directory structure
mkdir ${RUN_DIR}
mkdir ${RUN_DIR}/alignment_stats
mkdir ${RUN_DIR}/bam_lanewise
mkdir ${RUN_DIR}/bam_merged
mkdir ${RUN_DIR}/bam_merged/de_duped_bam_tmp
mkdir ${RUN_DIR}/bam_merged/de_duped_bam_tmp/picard_metrics

mkdir ${RUN_DIR}/bam_merged/filtered
mkdir ${RUN_DIR}/bam_merged/filtered/de_duped_bam_tmp
mkdir ${RUN_DIR}/bam_merged/filtered/de_duped_bam_tmp/picard_metrics

mkdir ${RUN_DIR}/blacklist
mkdir ${RUN_DIR}/log
mkdir ${RUN_DIR}/macs_broad
mkdir ${RUN_DIR}/macs_pvalue
mkdir ${RUN_DIR}/macs_pvalue/extended
mkdir ${RUN_DIR}/macs_pvalue/clean
mkdir ${RUN_DIR}/merge_broad
mkdir ${RUN_DIR}/merge_broad/clean
mkdir ${RUN_DIR}/merge_broad/extended

mkdir ${RUN_DIR}/merge_pvalue
mkdir ${RUN_DIR}/merge_pvalue/clean
mkdir ${RUN_DIR}/merge_pvalue/extended
mkdir ${RUN_DIR}/merge_pvalue/read_count
mkdir ${RUN_DIR}/merge_pvalue/annotated
mkdir ${RUN_DIR}/merge_pvalue/annotated/sorted
mkdir ${RUN_DIR}/merge_pvalue/annotated/part
mkdir ${RUN_DIR}/merge_pvalue/annotated/combined

mkdir ${RUN_DIR}/sam
mkdir ${RUN_DIR}/scripts
mkdir ${RUN_DIR}/tmp

cp ${HOME}/scripts/* ${RUN_DIR}/scripts/

rm -rf ${RUN_DIR}/bam_merged
cp -r ${LOCATION}/${run_backup}/bam_merged ${RUN_DIR}/

rm -f ${RUN_DIR}/bam_merged/de_duped_bam_tmp/*.bam
rm -f ${RUN_DIR}/bam_merged/de_duped_bam_tmp/picard_metrics/*.txt
rm -f ${RUN_DIR}/bam_merged/de_duped_bam_tmp/filtered/*.bam
rm -f ${RUN_DIR}/bam_merged/de_duped_bam_tmp/filtered/picard_metrics/*.txt

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
    mv images ${RUN_DIR}/overlap_output
fi

