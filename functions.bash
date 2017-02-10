#!/bin/bash

function _find_parent_dir()
{
    if [ "$0" = '-bash' ]; then
        __SCRIPTS_DIR=$(dirname ${BASH_SOURCE[$(expr ${#BASH_SOURCE[@]} - 1)]})
    fi
    if [ -z $__SCRIPTS_DIR ] ; then
        [ "$(dirname $0)" = '.' ] && __SCRIPTS_DIR=$(pwd) || __SCRIPTS_DIR=$(dirname $0);
        ls -l $__SCRIPTS_DIR | grep -q ^l && __SCRIPTS_DIR=$(dirname `ls -l $__SCRIPTS_DIR | awk '{print $NF}'`);
    fi
    export __SCRIPTS_DIR=$(echo "$__SCRIPTS_DIR" | sed -e 's,\\,/,g');
}

if [ -z ${PROCESSORS} ] ; then
    [ -z $__SCRIPTS_DIR ] && _find_parent_dir
    source $__SCRIPTS_DIR/configure.bash
fi

if ! typeset -f queue ; then
    [ -z $__SCRIPTS_DIR ] && _find_parent_dir
    source $__SCRIPTS_DIR/process_manager.bash
fi

##
# Execute Bowtie2 to align input files to the requested genome
#
function run_bowtie()
{
    if [ -z "$1" ] ; then
        echo "No input pattern provided"
        exit 1
    fi

    local pattern=$1
    local files=($(find -P ${FASTQ_DIR} -maxdepth 2 -type f | grep ${pattern}))
    local sam_dir="${RUN_DIR}/sam"
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/bowtie2Commands_${start_time}.txt"

    for file in ${files[*]}; do
        local noext=${file%%_R1_001.fastq.gz}
        local sample_name=$(basename ${noext})
        local output=${sam_dir}/${sample_name}.sam

        command="bowtie2 -p ${PROCESSORS} -q -x ${BOWTIE_INDEX} -U ${file}  -S ${output}"
        queue "${command}" "${command_log}" "block=true"
    done
    process_queue
}

##
# Convert sam file format to BAM format
function sam_to_bam()
{
    local sam_dir=${RUN_DIR}/sam
    local files=($(ls ${sam_dir}/*.sam))
    local bam_dir=${RUN_DIR}/bam_lanewise
    local tmp_dir=${RUN_DIR}/tmp
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/sam_to_bam_log_${start_time}.txt"

    for file in ${files[*]}; do
        local sample_name=$(basename ${file} _.sam)
        local output="${bam_dir}/${sample_name}.bam"
        local command="samtools view -bS -q 10 -F 260 ${file} -o ${output}"
        queue "${command}" "${command_log}"
    done
    process_queue
}

##
# Merges multiple bam files into a single bam, one per track.
#
# Arguments:
# This function takes a set of patterns/replacements as its arguments.
#
# For each group of lanes, a pattern and replacement must be provided
# in the form of "pattern#replacement"
#
# For example, if the output file is "AB123_Merged.bam" and the list of files constructing this merged is:
#     AB123_ENTR_1_RT_12.bam
#     AB123_ENTR_2_RT_12.bam
#     AB123_ENTR_3_RT_12.bam
#     AB123_ENTR_4_RT_12.bam
#
# then the argument matching this would be: "AB123_*_RT*#AB123_Merged"
#
function merge_bam_files()
{
    local arguments=($@)
    local patterns=(
        $(
            for pattern in ${arguments[@]}; do
                echo ${pattern} | cut -d\# -f1
            done
        )
    )

    local replacements=(
        $(
            for pattern in ${arguments[@]}; do
                echo ${pattern} | cut -d\# -f2
            done
        )
    )
    local bam_merged=${RUN_DIR}/bam_merged
    local input_dir=${RUN_DIR}/bam_lanewise
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/merge_bam_${start_time}.txt"

    for (( i=0; i < ${#patterns[@]}; i++ )); do
        local pattern="${patterns[$i]}"
        local replacement="${bam_merged}/${replacements[$i]}.bam"
        local files=($(ls $input_dir/${pattern}.bam))

        command="samtools merge - ${files[@]} | samtools sort -o ${replacement} -"
        queue "${command}" "${command_log}"
        queue "samtools index ${replacement}" "${command_log}" "waitfor=$?"
    done
    process_queue
}

function filter()
{
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local merged_bam=${RUN_DIR}/bam_merged
    local merged_bam_list=($(ls ${merged_bam}/*.bam))
    local bam_filtered_dir="${merged_bam}/${FILTERED}"
    local command_log="${LOG_DIR}/${FILTERED}${start_time}.log"
    local chromasones="${CHROM_LIST[@]}"

    for file in ${merged_bam_list[*]}; do
        name=$(basename $file)
        queue "samtools view -b ${file} $chromasones | samtools sort -o ${bam_filtered_dir}/${name} -" "${command_log}"
        queue "samtools index ${bam_filtered_dir}/${name}" "${command_log}" "waitfor=$?"
    done
    process_queue
}

##
# Remove duplications from all bam files in the input directory
#
# Arguments
# $1 Input directory
function remove_duplicates()
{
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local merged_bam="${RUN_DIR}/bam_merged/${FILTERED}"
    local temp_bam_folder="${merged_bam}/de_duped_bam_tmp"
    local picard_metrics_dir="${temp_bam_folder}/picard_metrics"
    local command_log="${LOG_DIR}/${FILTERED}de-duped_${start_time}.txt"
    local files=($(ls ${merged_bam}/*.bam))

    for file in ${files[*]}; do
        local sample_name=$(basename $file .bam)
        local mark_output=${sample_name}_remove_duplicates.bam
        local output="OUTPUT=${temp_bam_folder}/${mark_output}"
        local metrics="METRICS_FILE=${picard_metrics_dir}/$(basename ${file}).txt"

        local command="java -jar ${PICARD_LOCATION} MarkDuplicates INPUT=${file} ${output} $metrics REMOVE_DUPLICATES=TRUE"
        queue "$command" "${command_log}"
    done
    process_queue
}

function macs_pvalue()
{
    local precision=$1
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log=${LOG_DIR}/${FILTERED}/macs_pvalue_${start_time}.txt
    local bam_merged_dir=${RUN_DIR}/bam_merged/${FILTERED}
    local de_duped_dir=${bam_merged_dir}/de_duped_bam_tmp
    local macs_pvalue_dir=${RUN_DIR}/macs_pvalue/${FILTERED}
    local files=($(ls ${de_duped_dir}/*.bam))

    for file in ${files[*]}; do
        if echo $file | grep -qi input; then
            continue
        fi
        local input="${de_duped_dir}/$(ls $de_duped_dir | grep -i input | grep -i $(basename $file | cut -d_ -f1))"
        local name=$(basename $file .bam)
        local command="macs2 callpeak -t ${file} -c ${input} --format BAM -g mm -n"
        command="$command ${name}_${precision} --outdir ${macs_pvalue_dir} -p ${precision}"
        queue "${command}" "${command_log}"
    done
    process_queue
}

function remove_blacklist()
{
    _blacklist narrowPeak
    _blacklist summits

    local bed_folder=${RUN_DIR}/macs_pvalue/${FILTERED}
    local pvalue_clean_folder=${bed_folder}/clean
    local xls_list=($(ls ${bed_folder}/*.xls))
    local blacklist=${RUN_DIR}/../${GENOME}/blacklist/${GENOME_VERSION}-blacklist.bed

    command_log=${LOG_DIR}/${FILTERED}remove_blacklist_peaks_and_summits_xls_${start_time}.txt
    for file in ${xls_list[*]}; do
        filename=$(basename ${file} .xls)
        # sed breakdown:
        # The following sed command deletes everything that starts with # or is empty line
        #     x -> exchange pattern and hold buffer
        #     /^#\|^$/ -> search for lines beginning with # symbol
        #              -> or empty lines
        #     !g -> Copy hold space into pattern buffer
        #     !p -> don't print anything in current pattern space (e.g, delete everything that is in hold space)
        sed -n 'x;/^#\\|^$/!g;//!p' ${file} > ${pvalue_clean_folder}/${filename}_edited.xls

        local command="bedtools intersect -a ${pvalue_clean_folder}/${filename}_edited.xls"
        command="$command -b ${blacklist} -v > ${pvalue_clean_folder}/${filename}.xls"

        queue "$command" "${command_log}"
        index=$?
        echo $index
        queue "rm ${pvalue_clean_folder}/${filename}_edited.xls" "${command_log}" "waitfor=$index"
    done
    process_queue
}

function extend_pvalue()
{
    _extend "narrowPeak"
    _extend "summits"
    process_queue
}

function merge_pvalue()
{
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log=${LOG_DIR}/${FILTERED}merge_pvalue_${start_time}.txt
    local extended_bed_dir=${RUN_DIR}/macs_pvalue/${FILTERED}extended
    local known_sizes=${LOCATION}/${GENOME}/${GENOME_VERSION}/${GENOME_VERSION}.chrom.sizes
    local merge_dir=${RUN_DIR}/merge_pvalue/${FILTERED}
    local files=($(ls ${extended_bed_dir}/*.bed))

    for file in ${files[*]}; do
        local command="bedtools merge -d ${MERGE_RANGE} -i ${file} > ${merge_dir}/$(basename $file)"
        queue "${command}" "${command_log}"
    done
    process_queue
}

function read_counts()
{
    echo "Entering read_counts function"
    local bam_dir=${RUN_DIR}/bam_merged/${FILTERED}
    local merge_dir=${RUN_DIR}/merge_pvalue/${FILTERED}

    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local read_count=${merge_dir}/read_count
    local command_log=${LOG_DIR}/${FILTERED}read_counts_${start_time}.txt
    local files=($(ls ${merge_dir}/*.bed))

    for file in ${files[@]}; do
        if ! echo $file | grep -q peaks ; then
            continue
        fi
        local input_file=$(basename $file _remove_duplicates_${PRECISION}_peaks_extended_${EXTEND_RANGE}.bed)
        if [ -f ${bam_dir}/${input_file}.bam ] && [ ! -f "${read_count}/${input_file}.txt" ] ; then
            echo "Creating read counts for ${input_file}" | tee -a ${command_log}
            for peak_value in $(sed 's/\t/:/;s/\t/-/;' ${merge_dir}/$(basename $file)); do
                echo -n ${peak_value}$'\t'
                samtools view ${bam_dir}/${input_file}.bam ${peak_value} | wc -l
            done > ${read_count}/${input_file}.txt
        fi
    done
}

function annotate_read_counts()
{
    local merge_dir=${RUN_DIR}/merge_pvalue/${FILTERED}
    local read_count=${merge_dir}/read_count
    local annotated=${merge_dir}/annotated

    local sorted=${annotated}/sorted
    local partial=${annotated}/part
    local combined=${annotated}/combined
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/annotate_read_counts_${start_time}.log"

    for file in $(ls ${merge_dir}/*_${PRECISION}_peaks_extended_${EXTEND_RANGE}.bed); do
        local command="annotatePeaks.pl $file ${GENOME_VERSION} > ${annotated}/$(basename $file)";
        echo "Triggering '${command}'"
        queue "${command}" "${command_log}"
    done
    process_queue

    for file in $(ls ${read_count}/*.txt); do
        name=$(basename $file .txt);
        sed 1d ${annotated}/${name}_remove_duplicates_${PRECISION}_peaks_extended_${EXTEND_RANGE}.bed > ${partial}/${name}.bed;
    done

    for file in $(ls ${partial}/*.bed); do
        sort -k2,2 -k3,3n -k4,4n $file > ${sorted}/$(basename $file);
    done

    for file in $(ls ${read_count}/*.txt); do
        name=$(basename $file .txt);
        awk -v filename=$file '{ getline v < filename ; split( v, a ); print a[2],"\t"$0 }' ${sorted}/${name}.bed > ${combined}/${name}.bed;
    done
}

# Set up the run directory structure
function create_structure()
{
    [ ! -z $FILTERED ] && FILTERED='filtered/'
    mkdir -p ${RUN_DIR}

    [ ! -d ${RUN_DIR}/bam_lanewise ] &&
        mkdir -p ${RUN_DIR}/bam_lanewise

    [ ! -d ${RUN_DIR}/bam_merged/${FILTERED} ] &&
        mkdir -p ${RUN_DIR}/bam_merged${FILTERED}

    [ ! -d ${RUN_DIR}/bam_merged/${FILTERED}de_duped_bam_tmp ] &&
        mkdir -p ${RUN_DIR}/bam_merged/${FILTERED}de_duped_bam_tmp

    [ ! -d ${RUN_DIR}/bam_merged/${FILTERED}de_duped_bam_tmp/picard_metrics ] &&
        mkdir -p ${RUN_DIR}/bam_merged/${FILTERED}de_duped_bam_tmp/picard_metrics

    [ ! -d ${RUN_DIR}/blacklist/${FILTERED} ] &&
        mkdir -p ${RUN_DIR}/blacklist/${FILTERED}

    [ ! -d ${RUN_DIR}/log/${FILTERED} ] &&
        mkdir -p ${RUN_DIR}/log/${FILTERED}

    [ ! -d ${RUN_DIR}/macs_pvalue/${FILTERED} ] &&
        mkdir -p ${RUN_DIR}/macs_pvalue/${FILTERED}

    [ ! -d ${RUN_DIR}/macs_pvalue/${FILTERED}extended ] &&
        mkdir -p ${RUN_DIR}/macs_pvalue/${FILTERED}extended

    [ ! -d ${RUN_DIR}/macs_pvalue/${FILTERED}clean ] &&
        mkdir -p ${RUN_DIR}/macs_pvalue/${FILTERED}clean

    [ ! -d ${RUN_DIR}/merge_pvalue/${FILTERED} ] &&
        mkdir -p ${RUN_DIR}/merge_pvalue/${FILTERED}

    [ ! -d ${RUN_DIR}/merge_pvalue/${FILTERED}clean ] &&
        mkdir -p ${RUN_DIR}/merge_pvalue/${FILTERED}clean

    [ ! -d ${RUN_DIR}/merge_pvalue/${FILTERED}extended ] &&
        mkdir -p ${RUN_DIR}/merge_pvalue/${FILTERED}extended

    [ ! -d ${RUN_DIR}/merge_pvalue/${FILTERED}read_count ] &&
        mkdir -p ${RUN_DIR}/merge_pvalue/${FILTERED}read_count

    [ ! -d ${RUN_DIR}/merge_pvalue/${FILTERED}annotated ] &&
        mkdir -p ${RUN_DIR}/merge_pvalue/${FILTERED}annotated

    [ ! -d ${RUN_DIR}/merge_pvalue/${FILTERED}annotated/sorted ] &&
        mkdir -p ${RUN_DIR}/merge_pvalue/${FILTERED}annotated/sorted

    [ ! -d ${RUN_DIR}/merge_pvalue/${FILTERED}annotated/part ] &&
        mkdir -p ${RUN_DIR}/merge_pvalue/${FILTERED}annotated/part

    [ ! -d ${RUN_DIR}/merge_pvalue/${FILTERED}annotated/combined ] &&
        mkdir -p ${RUN_DIR}/merge_pvalue/${FILTERED}annotated/combined

    [ ! -d ${RUN_DIR}/sam ] &&
        mkdir -p ${RUN_DIR}/sam

    [ ! -d ${RUN_DIR}/scripts ] &&
        mkdir -p ${RUN_DIR}/scripts

    [ ! -d ${RUN_DIR}/tmp ] &&
        mkdir -p ${RUN_DIR}/tmp
}

function _blacklist()
{
    local peak_type=$1
    local bed_folder=${RUN_DIR}/macs_pvalue/${FILTERED}
    local pvalue_clean_folder=${bed_folder}/clean

    local blacklist=${RUN_DIR}/../${GENOME}/blacklist/${GENOME_VERSION}-blacklist.bed

    local bed_list=($(ls ${bed_folder}/*${peak_type}*))
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")

    local command_log=${LOG_DIR}/${FILTERED}remove_blacklist_${peak_type}_${start_time}.txt
    for file in ${bed_list[*]}; do
        local filename=$(basename ${file} .${peak_type})
        local edited=${filename}_${peak_type}_clean.bed

        command="bedtools intersect -a ${file} -b ${blacklist} -v > ${pvalue_clean_folder}/${edited}"
        queue "$command" "${command_log}"
    done
}

function _extend()
{
    local peak_type=$1
    local bed_dir=${RUN_DIR}/macs_pvalue/${FILTERED}clean
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log=${LOG_DIR}/${FILTERED}extend_${peak_type}_${start_time}.txt
    local known_sizes=${RUN_DIR}/../${GENOME}/${GENOME_VERSION}/${GENOME_VERSION}.chrom.sizes
    local extended="${bed_dir}/../extended"
    local files=($(ls ${bed_dir}/*_${peak_type}_clean.bed))

    for file in ${files[*]}; do
        local filename=$(basename $file _${peak_type}_clean.bed)
        local edited=${filename}_extended_${EXTEND_RANGE}.bed

        local command="bedtools slop -i ${file} -b ${EXTEND_RANGE} -g ${known_sizes} |"
        command="${command} bedtools merge -i - > ${extended}/${edited}"
        queue "$command" "${command_log}"
    done
}

function _clean()
{
    rm -rf ${RUN_DIR}/bam_merged
    cp -r ${LOCATION}/${run_backup}/bam_merged ${RUN_DIR}/

    rm -f ${RUN_DIR}/bam_merged/de_duped_bam_tmp/*.bam
    rm -f ${RUN_DIR}/bam_merged/de_duped_bam_tmp/picard_metrics/*.txt
    rm -f ${RUN_DIR}/bam_merged/de_duped_bam_tmp/filtered/*.bam
    rm -f ${RUN_DIR}/bam_merged/de_duped_bam_tmp/filtered/picard_metrics/*.txt
}

