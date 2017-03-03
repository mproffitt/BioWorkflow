#!/bin/bash
#
# Configuration for the BioInformatics workflow

export GENOME='mouse'
export GENOME_VERSION='mm10'
export RUN_DATE='2017-02-04'
export PICARD_VERSION='2.5.0'
export XDG_CONFIG_HOME=~/.config
export EXTEND_RANGE=250
export MERGE_RANGE=100
export PRECISION='0.01'

export LOCATION="${HOME}/Sequences"

export BOWTIE_INDEX="${LOCATION}/${GENOME}/${GENOME_VERSION}/$(echo ${GENOME_VERSION} | cut -d- -f1)"
export FASTQ_DIR="${LOCATION}/${RUN_DATE}/Unaligned"

export RUN_DIR=${LOCATION}/run
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
