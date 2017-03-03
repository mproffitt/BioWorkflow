#
# Wrapper script for GL21, GL25 and GL30
#
source bio_workflow.bash

gl21_filters=(
    "GL21_*_Input*#GL21_Input"
    "GL21_*_H3K4me3*#GL21_2_H3K4me3_CTTGTA"
    "GL21_*_Hdac1*#GL21_3_Hdac1_GGCTAC"
    "GL21_*_Hdac2*#GL21_4_Hdac2_TAGCTT"
)

function gl21()
{
    backup_and_recreate
    execute 'GL21_.*.fastq*' ${gl21_filters[@]}
    mv $RUN_DIR ${LOCATION}/run-gl21
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
    backup_and_recreate
    execute 'GL25-.*.fastq*' ${gl25_filters[@]}
    mv $RUN_DIR ${LOCATION}/run-gl25
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

    backup_and_recreate
    execute 'GL30.*RT.*fastq*' ${gl30_filters[@]}
    mv $RUN_DIR ${LOCATION}/run-gl30
}

function gl30_filter()
{
    export RUN_DATE='2016-07-18'
    export FASTQ_DIR="${LOCATION}/${RUN_DATE}/Unaligned"
    export RUN_DIR=${LOCATION}/run-gl30
    export LOG_DIR=${RUN_DIR}/log
    execute -f 'GL30.*RT.*fastq*' ${gl30_filters[@]}
}

gl21
gl21_filter
gl25
gl25_filter
gl30
gl30_filter

