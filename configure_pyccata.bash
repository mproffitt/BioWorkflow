
[ -z "${INPUT_DIR}" ] && INPUT_DIR=${RUN_DIR}/merge_pvalue/${FILTERED}annotated/combined
INPUT_FILES=(
    $(find -P ${INPUT_DIR} -type f -exec basename {} \; )
)

LIMITS=(
    $(find -P ${INPUT_DIR} -type f -exec basename {} \; |
        tr '[:upper:]' '[:lower:]' | cut -d. -f1 | awk -F_ '{print $NF}' |
            xargs | awk -v extended="250" '{ for(i = 1; i <= NF; i++) { print $i":" extended; } }'
    )
)

[ ! -z "$FILTERED" ] && FILTERED="/$FILTERED"
echo "Creating config file"
read -r -d '' CONFIG <<EOF
{
    "manager": "csv",
    "reporting": "null",
    "modules":[
        "csv"
    ],
    "csv": {
        "datapath": "${INPUT_DIR}",
        "input_files": [
            $(
                for (( i=0; i < ${#INPUT_FILES[@]}; i++ )); do
                    file="${INPUT_FILES[$i]}"
                    if [ $i -lt $(expr ${#INPUT_FILES[@]} - 1) ]; then
                        echo "\"$file\","
                    else
                        echo "\"$file\""
                    fi
                done
            )
        ],
        "output_directory": "overlap_output${FILTERED}",
        "namespace": "bioinformatics"
    },
    "replacements":[],
    "report": {
        "datapath": "data",
        "title": "",
        "subtitle": "",
        "abstract": "",
        "sections":[
            {
                "title": "Hd2Lox with Hd1",
                "abstract":"",
                "level": 1,
                "structure": [
                    {
                        "type": "graph",
                        "title": "Up-Set of Chromasone summary for 2 or more files",
                        "content": {
                            "width": 100,
                            "graphtype": "upset",
                            "structure": null,
                            "query": {
                                "query": "",
                                "fields": ["gene_name", "start", "end", "read_count", "peak_id"],
                                "max_results": false,
                                "group_by": "gene_name",
                                "collate": {
                                    "method": "combinatorics",
                                    "fields": ["gene_name", "start", "end", "read_count", "peak_id"],
                                    "query": {
                                        "inclusive": "(start_x >= (start_y - {left_limit}) and (start_x <= (end_y + {right_limit})))",
                                        "exclusive": "(start_y >= (start_x - {left_limit}) and (start_y <= (end_x + {right_limit})))"
                                    },
                                    "join": {
                                        "method": "outer",
                                        "column": "gene_name"
                                    },
                                    "limits": {
                                        $(
                                            for (( i=0; i < ${#LIMITS[@]}; i++ )); do
                                                limit=${LIMITS[$i]}
                                                if [ $i -lt $(expr ${#LIMITS[@]} - 1) ]; then
                                                    echo "\"$(echo $limit | cut -d: -f1)\": $(echo $limit | cut -d: -f2),"
                                                else
                                                    echo "\"$(echo $limit | cut -d: -f1)\": $(echo $limit | cut -d: -f2)"
                                                fi
                                            done
                                        )
                                    }
                                }
                            }
                        }
                    }
                ]
            }
        ]
    }
}
EOF
