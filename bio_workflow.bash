#!/bin/bash
source ~/scripts/configure.bash
source ~/scripts/process_manager.bash
source ~/scripts/functions.bash

pattern="GL30*_RT*.gz"

# Temporary until pyccata is complete
function setup_pyccata()
{
    cd ${RUN_DIR}
    ln -s $(python3 -c 'import site; print(site.getsitepackages()[0])')/pyccata pyccata
    mkdir overlap_output
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
            run_backup=$1
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

function execute()
{
    local filter=''
    if [ ! -z $1 ] && [ "$1" = '-f' ]; then
        filter='filtered'
        shift
    fi

    local bowtie_pattern=$1
    shift
    local bam_filter=$@

    run_bowtie $bowtie_pattern
    wait_for 'last'
    sam_to_bam
    wait_for 'last'
    merge_bam_files $bam_filter
    wait_for 'last'
    remove_duplicates $filter
    wait_for 'last'
    macs_pvalue ${PRECISION} $filter
    wait_for 'last'
    remove_blacklist $filter
    wait_for 'last'
    extend_pvalue
    wait_for 'last'
    merge_pvalue
    wait_for 'last'
    read_counts $filter
    annotate_read_counts
}

function gl21()
{
    #backup_and_recreate
    local bam_filters=(
        "GL21_*_Input*.sam#GL21_Input"
        "GL21_*_H3K4me3*.sam#GL21_2_H3K4me3_CTTGTA"
        "GL21_*_Hdac1*.sam#GL21_3_Hdac1_GGCTAC"
        "GL21_*_Hdac2*.sam#GL21_4_Hdac2_TAGCTT"
    )
    execute 'GL21_.*.fastq*' ${bam_filters[@]}
    print_queue
    process_queue
    #setup_pyccata
}

function gl25()
{
    backup_and_recreate
    local bam_filters=(
        "GL25*input*.sam#GL25_Input"
        "GL25*loxlox-H3K4me3*.sam#GL25_loxlox_H3k4me3"
        "GL25*loxlox-H3K27me3*.sam#GL25_loxlox_H3k27me3"
        "GL25*loxlox-Hd1*.sam#GL25_loxlox_H3k4me3"
        "GL25*loxlox-Hd2*.sam#GL25_loxlox_H3k27me3"
        "GL25*4oht-H3K4me3*.sam#GL25_4oht_H3k4me3"
        "GL25*4oht-H3K27me3*.sam#GL25_4oht_H3k27me3"
        "GL25*4oht-Hd1*.sam#GL25_4oht_H3k4me3"
        "GL25*4oht-Hd2*.sam#GL25_4oht_H3k4me3"
    )
    execute 'GL25-.*.fastq*' ${bam_filters[@]}
    print_queue
    process_queue
}

gl25
