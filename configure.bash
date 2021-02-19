#!/bin/bash
#
# Configuration for the BioInformatics workflow
export DISABLE_PYCCATA=1

export GENOME='mouse'
export GENOME_VERSION='mm10'
export RUN_DATE='2021-02-01'
export PICARD_VERSION='2.24.1'
export XDG_CONFIG_HOME=~/.config
export EXTEND_RANGE=1
export MERGE_RANGE=1
export PRECISION='0.0085'
export INPUTFILTER='IgG'

export LOCATION="/home/mproffitt/2021-01-26_ACD5RFANXX"

export BOWTIE_INDEX="${LOCATION}/${GENOME}/${GENOME_VERSION}/$(echo ${GENOME_VERSION} | cut -d- -f1)"
export FASTQ_DIR="${LOCATION}/${RUN_DATE}/Unaligned"

export RUN_DIR=${LOCATION}/run-${RUN_DATE}
export LOG_DIR=${RUN_DIR}/log
export PICARD_LOCATION="${LOCATION}/tools/picard-tools-${PICARD_VERSION}/picard.jar"
export PROCESSORS=$(grep processor /proc/cpuinfo | wc -l)
export THREAD_SLEEP=1
export CHROM_LIST=(
    'chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9'
    'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16'
    'chr17' 'chr18' 'chr19' 'chrX' 'chrY'
)
if [ "${GENOME}" = 'human' ] ; then
    CHROM_LIST+=( 'chr20' 'chr21')
fi

export QUEUE=()
export PROCESSES=()
export STATUSES=()
export ACTIVE=0
export MANAGER_RUNNING=1
export CONFIG_DEFINED=1
