#
# Wrapper script for GL21, GL25 and GL30
#
source configure.bash
source bio_workflow.bash

set1_filters=(
    XXA_10IK*#XXA_10IK
    XXA_115NIC*#XXA_115NIC
    XXA_115Su*#XXA_115Su
    XXA_11CD31*#XXA_11CD31
    XXA_15IK*#XXA_15IK
    XXA_15NI*#XXA_15NI
    XXA_15SU*#XXA_15SU
    XXA_IgG*#XXA_IgG
    XXA_KOFull*#XXA_KOFull
    XXA_KOH3*#XXA_KOH3
    XXA_Wtfull*#XXA_Wtfull
    XXA_WTH3*#XXA_WTH3
)

function set1_filter()
{
    #create_structure
    export FASTQ_DIR="${LOCATION}/${RUN_DATE}/Unaligned"
    export RUN_DATE='2021-02-15'
    export RUN_DIR=${LOCATION}/run-${RUN_DATE}
    export LOG_DIR=${RUN_DIR}/log
    execute 'XXA_.*.fastq*' ${set1_filters[@]}
}

set1_filter
