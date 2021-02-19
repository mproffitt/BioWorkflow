#!/bin/bash

function _find_parent_dir()
{
    if [ "$0" = '-bash' ]; then
        __SCRIPTS_DIR=$(dirname ${BASH_SOURCE[$(expr ${#BASH_SOURCE[@]} - 1)]})
    fi
    echo "1 = __SCRIPTS_DIR='$__SCRIPTS_DIR'"
    if [ -z $__SCRIPTS_DIR ] ; then
        [ "$(dirname $0)" = '.' ] && __SCRIPTS_DIR=$(pwd) || __SCRIPTS_DIR=$(dirname $0);
        echo "2 = __SCRIPTS_DIR='$__SCRIPTS_DIR'"
        #ls -l $__SCRIPTS_DIR | grep -q ^l && __SCRIPTS_DIR=$(dirname `ls -l $__SCRIPTS_DIR | awk '{print $NF}'`);
    fi
    echo "3 = __SCRIPTS_DIR='$__SCRIPTS_DIR'"
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

function run_bowtie()
{
    echo "PATTERN IS $1"; echo; echo;
    if [ -z "$1" ] ; then
        echo "No input pattern provided"
        return 1
    fi

    local pattern=$1
    local files=($(find -P ${FASTQ_DIR} -maxdepth 2 -type f | grep ${pattern}))
    local sam_dir="${RUN_DIR}/00sam"
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/bowtie2Commands_${start_time}"

    mkdir -p ${sam_dir}

    echo "Building list ${files[@]}"
    for sample in $(for file in ${files[@]}; do echo $(basename ${file%%_read*.fastq.gz}); done | sort | uniq); do
        local command="bowtie2 --local --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${PROCESSORS} -q -x ${BOWTIE_INDEX}"

        let i=1
        for match in $(ls ${FASTQ_DIR}/$sample*); do
            command="${command} -${i} ${match}"
            let i=$i+1
        done;
        command="${command} -S ${sam_dir}/${sample}.sam"
        echo "Triggering ${command}"; echo
        queue "${command}" "${command_log}_${sample}.txt" "block=true"
    done
    process_queue
}

##
# Convert sam file format to BAM format
function sam_to_bam()
{
    local sam_dir=${RUN_DIR}/00sam
    local files=($(ls ${sam_dir}/*.sam))
    local bam_dir=${RUN_DIR}/01bam
    local tmp_dir=${RUN_DIR}/tmp
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/sam_to_bam_log_${start_time}.txt"

    mkdir -p ${bam_dir}
    mkdir -p ${tmp_dir}

    for file in ${files[*]}; do
        local sample_name=$(basename ${file} .sam)
        local output="${bam_dir}/${sample_name}.bam"
        local command="samtools view -bS -q 10 -F 4 ${file} | samtools sort -o ${output} -"
        echo "Triggering ${command}"; echo
        queue "${command}" "${command_log}"
        local ret=$?

        command="samtools index ${bam_dir}/${sample_name}.bam"
        echo "Triggering ${command}"; echo
        queue "$command" "${command_log}" "waitfor=${ret}"
    done
    process_queue
}

function filter()
{
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local bam_dir=${RUN_DIR}/01bam
    local merged_bam_list=($(ls ${bam_dir}/*.bam))

    local filtered_dir=${RUN_DIR}/02filtered
    local command_log="${LOG_DIR}/${FILTERED}${start_time}.log"
    local chromasones="${CHROM_LIST[@]}"

    mkdir -p ${filtered_dir}
    for file in ${merged_bam_list[*]}; do
        name=$(basename $file)

        command="samtools view -b ${file} $chromasones | samtools sort -o ${filtered_dir}/${name} -"
        echo "Triggering ${command}"; echo
        queue "$command" "${command_log}"
        local ret=$?
        command="samtools index ${filtered_dir}/${name}"
        echo "Triggering ${command}"; echo
        queue "$command" "${command_log}" "waitfor=${ret}"
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
    local filtered_dir=${RUN_DIR}/02filtered

    local unique_dir=${RUN_DIR}/03unique
    local tmp_dir=${RUN_DIR}/tmp

    local picard_metrics_dir="${unique_dir}/picard_metrics"
    local command_log="${LOG_DIR}/${FILTERED}de-duped_${start_time}.txt"

    local files=($(ls ${filtered_dir}/*.bam))
    mkdir -p ${unique_dir}
    mkdir -p ${picard_metrics_dir}

    for file in ${files[*]}; do
        local sample_name=$(basename $file .bam)
        local mark_output=${sample_name}_remove_duplicates.bam
        local output="OUTPUT=${unique_dir}/${mark_output}"
        local metrics="METRICS_FILE=${picard_metrics_dir}/$(basename ${file}).txt"

        local command="java -jar ${PICARD_LOCATION} MarkDuplicates INPUT=${file} ${output} $metrics REMOVE_DUPLICATES=TRUE"
        echo "Triggering ${command}"; echo
        queue "$command" "${command_log}"
    done
    process_queue
}

function macs_pvalue()
{
    local precision=$1
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log=${LOG_DIR}/${FILTERED}/macs_pvalue_${start_time}.txt
    local unique_dir=${RUN_DIR}/03unique

    local macs_pvalue_dir=${RUN_DIR}/04macs_pvalue
    local files=($(ls ${unique_dir}/*.bam))

    mkdir -p ${macs_pvalue_dir}

    for file in ${files[*]}; do
        if echo $file | grep -qi ${INPUTFILTER}; then
            continue
        fi
        local input="${unique_dir}/$(ls ${unique_dir} | grep -i ${INPUTFILTER} | grep -i $(basename $file | cut -d_ -f1))"
        local name=$(basename $file .bam)
        #local command="macs2 callpeak -t ${file} -c ${input} --format BAM --broad --min-length 2000 --broad-cutoff 0.1 -g mm -n"
        local command="macs2 callpeak -t ${file} -c ${input} --format BAM -g mm -n"
        command="$command ${name}_${precision} --outdir ${macs_pvalue_dir} -p ${precision}"
        queue "${command}" "${command_log}"
    done
    process_queue
}

function remove_blacklist()
{
    local bed_folder=${RUN_DIR}/04macs_pvalue
    local pvalue_clean_folder=${RUN_DIR}/05clean
    local xls_list=($(ls ${bed_folder}/*.xls))
    local blacklist=${RUN_DIR}/../${GENOME}/blacklist/${GENOME_VERSION}-blacklist.bed

    mkdir -p ${pvalue_clean_folder}

    _blacklist narrowPeak
    _blacklist summits
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
    local extended_bed_dir=${RUN_DIR}/06extended
    local known_sizes=${LOCATION}/${GENOME}/${GENOME_VERSION}/${GENOME_VERSION}.chrom.sizes
    local merge_dir=${RUN_DIR}/07merge_pvalue
    local files=($(ls ${extended_bed_dir}/*.bed))

    mkdir -p ${merge_dir}

    for file in ${files[*]}; do
        local command="bedtools merge -d ${MERGE_RANGE} -i ${file} > ${merge_dir}/$(basename $file)"
        queue "${command}" "${command_log}"
    done
    process_queue
}

function read_counts()
{
    echo "Entering read_counts function"
    local bam_dir=${RUN_DIR}/02filtered
    local merge_dir=${RUN_DIR}/06extended

    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local read_count_dir=${RUN_DIR}/08read_count
    local command_log=${LOG_DIR}/${FILTERED}read_counts_${start_time}.txt
    local files=($(ls ${merge_dir}/*.bed))

    mkdir -p ${read_count_dir}

    for file in ${files[@]}; do
        if ! echo $file | grep -q peaks ; then
            continue
        fi
        local input_file=$(basename $file _remove_duplicates_${PRECISION}_peaks_extended_${EXTEND_RANGE}.bed)
        if [ -f ${bam_dir}/${input_file}.bam ] && [ ! -f "${read_count_dir}/${input_file}.txt" ] ; then
            echo "Creating read counts for ${input_file}" | tee -a ${command_log}
            for peak_value in $(sed 's/\t/:/;s/\t/-/;' ${merge_dir}/$(basename $file)); do
                echo -n ${peak_value}$'\t'
                samtools view ${bam_dir}/${input_file}.bam ${peak_value} | wc -l
            done > ${read_count_dir}/${input_file}.txt
        fi
    done
}

function annotate_read_counts()
{
    local merge_dir=${RUN_DIR}/07merge_pvalue
    local read_count=${RUN_DIR}/08read_count
    local annotated=${RUN_DIR}/09annotated

    local sorted=${annotated}/sorted
    local partial=${annotated}/part
    local combined=${annotated}/combined
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log="${LOG_DIR}/annotate_read_counts_${start_time}.log"

    mkdir -p ${annotated} ${sorted} ${partial} ${combined}

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
    mkdir -p ${RUN_DIR}
}

function _blacklist()
{
    local peak_type=$1
    local bed_folder=${RUN_DIR}/04macs_pvalue
    local pvalue_clean_folder=${RUN_DIR}/05clean

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
    local bed_dir=${RUN_DIR}/05clean
    local start_time=$(date +"%Y-%m-%d_%H_%M_%S")
    local command_log=${LOG_DIR}/${FILTERED}extend_${peak_type}_${start_time}.txt
    local known_sizes=${RUN_DIR}/../${GENOME}/${GENOME_VERSION}/${GENOME_VERSION}.chrom.sizes
    local extended="${RUN_DIR}/06extended"
    local files=($(ls ${bed_dir}/*_${peak_type}_clean.bed))

    mkdir -p ${RUN_DIR}/06extended

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

