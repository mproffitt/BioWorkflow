#!/bin/bash
source ~/scripts/functions.bash

# Temporary until pyccata is complete
function setup_pyccata()
{
    cd ${RUN_DIR}
    ln -sf $(python3 -c 'import site; print(site.getsitepackages()[0])')/pyccata pyccata
    [ ! -d overlap_output ] && mkdir overlap_output
    INPUT_FILES=($1)
    LIMITS=($2)
    source $__SCRIPTS_DIR/configure_pyccata.bash
    echo "$CONFIG" > bio.json
    cd -
}

# Sets the backup name to the earliest modified date of files under the run directory
# if $1 is provided, will rename the folder as "run-<name>" instead
function backup_and_recreate()
{
    local run_backup=''
    if [ -d "${RUN_DIR}" ]; then
        if [ ! -z $1 ] ; then
            run_backup="run-$1"
        else
            local year=$(date +%Y)
            run_backup='run-'$(
                date -d "$(ls -al ${RUN_DIR} | awk -v year=$year '{print $7,$6,year}' | tail -n +2 | sort | uniq | head -1)" +"%Y-%m-%d"
            )
        fi

        mv ${RUN_DIR} ${LOCATION}/${run_backup}
    fi
    create_structure
    cp ${HOME}/scripts/* ${RUN_DIR}/scripts/
}

function check_bowtie_pattern()
{
    if ! find -P ${FASTQ_DIR} -maxdepth 2 -type f | grep -q $1; then
        echo "Input pattern for bowtie does not recognise any files" >&2
        echo "Please check the pattern and try again"
        exit 1
    fi
}

function check_bam_filter_patterns()
{
    local sam_patterns=$@
    local cwd=$(pwd)
    local valid=false
    cd ${FASTQ_DIR}
    for ((i=0; i < ${#sam_patterns[@]}; i++)); do
        local pattern=${sam_patterns[@]}
        if ! $(ls ${pattern} &>/dev/null); then
            valid=true
        fi
    done
    cd ${cwd}
    if ! $valid; then
        echo 'One or more input patterns for the bam_merge function does not match any files.' >&2
        echo 'Please check your patterns and try again.' >&2
        echo 'Patterns should not contain any file extensions and can be tested with `ls '${FASTQ_DIR}'/<pattern>`' >&2
        echo 'where <pattern> is the input pattern you wish to test' >&2
        exit 1
    fi
}

function execute()
{
    unset INPUT_DIR
    local filter=''
    if [ ! -z $1 ] && [ "$1" = '-f' ]; then
        echo "In filter config"
        # reset filter location
        export FILTERED='filtered/'
        create_structure
        shift
    fi

    local bowtie_pattern=$1
    shift
    if [ ! -z "$bowtie_pattern" ]; then
        check_bowtie_pattern $bowtie_pattern
    fi

    local bam_filter=$@
    if [ ${#bam_filter[@]} -gt 0 ]; then
        check_bam_filter_patterns $@
    fi

    if [ -z $FILTERED ] || ([ ! -z $FILTERED ] && [ $(ls ${RUN_DIR}/sam | wc -l) -eq 0 ]) ; then
        run_bowtie $bowtie_pattern
        sam_to_bam
        merge_bam_files $bam_filter
    fi

    if [ ! -z $FILTERED ] ; then
        filter
    fi
    remove_duplicates
    macs_pvalue ${PRECISION}
    remove_blacklist
    extend_pvalue
    merge_pvalue
    read_counts
    annotate_read_counts
    process_queue
    setup_pyccata

    cd ${RUN_DIR}
    ${HOME}/$__SCRIPTS_DIR/annotation.py
    cd ${HOME}
}

gl21_filters=(
    "GL21_*_Input*#GL21_Input"
    "GL21_*_H3K4me3*#GL21_2_H3K4me3_CTTGTA"
    "GL21_*_Hdac1*#GL21_3_Hdac1_GGCTAC"
    "GL21_*_Hdac2*#GL21_4_Hdac2_TAGCTT"
)

function gl21()
{
    #backup_and_recreate
    export RUN_DIR=${LOCATION}/run-gl21
    export LOG_DIR=${RUN_DIR}/log
    execute 'GL21_.*.fastq*' ${gl21_filters[@]}
}

function gl21_filter()
{
    export RUN_DATE='2017-02-04'
    export FASTQ_DIR="${LOCATION}/${RUN_DATE}/Unaligned"
    export RUN_DIR=${LOCATION}/run-gl21
    export LOG_DIR=${RUN_DIR}/log
    execute -f 'GL21_.*.fastq*' ${gl21_filters[@]}
}

gl25_filters=(
    "GL25*input*#GL25_Input"
    "GL25*loxlox-H3K4me3*#GL25_loxlox_H3k4me3"
    "GL25*loxlox-H3K27me3*#GL25_loxlox_H3k27me3"
    "GL25*loxlox-Hd1*#GL25_loxlox_Hd1"
    "GL25*loxlox-Hd2*#GL25_loxlox_Hd2"
    "GL25*4oht-H3K4me3*#GL25_4oht_H3k4me3"
    "GL25*4oht-H3K27me3*#GL25_4oht_H3k27me3"
    "GL25*4oht-Hd1*#GL25_4oht_Hd1"
    "GL25*4oht-Hd2*#GL25_4oht_Hd2"
)

function gl25()
{
    #backup_and_recreate
    export RUN_DIR=${LOCATION}/run-gl25
    export LOG_DIR=${RUN_DIR}/log
    execute 'GL25-.*.fastq*' ${gl25_filters[@]}
}

function gl25_filter()
{

    export RUN_DATE='2017-02-04'
    export FASTQ_DIR="${LOCATION}/${RUN_DATE}/Unaligned"
    export RUN_DIR=${LOCATION}/run-gl25
    export LOG_DIR=${RUN_DIR}/log
    execute -f 'GL25-.*.fastq*' ${gl25_filters[@]}
}

gl30_filters=(
    "GL30*Input*#GL30_Input"
    "GL30*Hd2lox_Hd1*#GL30_Hd2lox_Hd1"
    "GL30*Hd2lox-Hd2*#GL30_Hd2lox_Hd2"
    "GL30*Hd2delta_Hd1*#GL30_Hd2delta_Hd1"
    "GL30*Hd2delta_Hd2*#GL30_Hd2delta_Hd2"
    "GL30*_RT7_S8*#GL30_CoREST"
)

function gl30()
{
    export RUN_DATE='2016-07-18'
    export FASTQ_DIR="${LOCATION}/${RUN_DATE}/Unaligned"

    export RUN_DIR=${LOCATION}/run-gl30
    export LOG_DIR=${RUN_DIR}/log
    #backup_and_recreate 'gl25'
    execute 'GL30.*RT.*fastq*' ${gl30_filters[@]}
}

function gl30_filter()
{
    export RUN_DATE='2016-07-18'
    export FASTQ_DIR="${LOCATION}/${RUN_DATE}/Unaligned"
    export RUN_DIR=${LOCATION}/run-gl30
    export LOG_DIR=${RUN_DIR}/log
    execute -f 'GL30.*RT.*fastq*' ${gl30_filters[@]}
}

if [ "$1" = '-f' ] ; then
    gl30_filter
    gl25_filter
    gl21_filter
else
    echo ''
    gl30
    gl25
    gl21
fi

