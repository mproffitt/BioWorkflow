#!/bin/bash
#
# High performance pipeline for analysing ChIPSeq data
#
# This pipeline sets up and executes a pipeline to convert ChIPSeq data
# taken from the sequencer in fastq.gz format all the way to final analysis.
#
# To execute this pipeline, the following dependencies need to exist in $PATH:
#
#   * bowtie2  (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
#   * samtools (http://samtools.sourceforge.net/)
#   * bedtools (http://bedtools.readthedocs.io/en/latest/)
#   * homer    (http://homer.ucsd.edu/homer/)
#   * java
#   * python3
#
# Python3 is required for executing the pyccata overlap flow. This software is still
# in an alpha development stage and is due to be published in the summer of 2017.
#   A Preliminary copy of pyccata may be obtained by emailing <mproffitt@jitsc.co.uk>
#   with the subject line 'pyccata prelim request'
#
# pyccata can be disabled by defining an environment variable of `DISABLE_PYCCATA=1`
#
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
                date -d "$(
                    ls -al ${RUN_DIR} | awk -v year=$year '{print $7,$6,year}' |
                        tail -n +2 | sort | uniq | head -1)" +"%Y-%m-%d"
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
    if [ ! -z $DISABLE_PYCCATA ]; then
        setup_pyccata
        cd ${RUN_DIR}
        ${HOME}/$__SCRIPTS_DIR/annotation.py
        cd ${HOME}
    fi
}

