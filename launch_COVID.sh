#!/bin/bash

#./scriptTestDatasets ../dataset ../Tools ./out_table 8

dataset_directory=$1;
tools_directory=$2;
output_directory=$3;
alpha=$4

###### ADD TEXT LENGTH N, NUMBER OF RUNS R,  R' 

pattern_length=("8" "16" "32" "64")
##########
#needs to change the input data types (len seqs, BWT len, ...) in Parameters.h
#dataset_array=("pdb_seqres.fa"); 
##########

echo -e "M_LF_build\ndataset,%CPU,WALL_CLOCK,RAM,WALL_CLOCK_BS,n,r,r',DEG_SYMBOLS,STRINGS" &> ${output_directory}/table_build_MLF.txt
echo -e "MLF_search\ndataset,%CPU,WALL_CLOCK,RAM,WALL_CLOCK_BS,FOUND,NOT_FOUND" &> ${output_directory}/table_MLF_search.txt
echo -e "MLF_search\ndataset,%CPU,WALL_CLOCK,RAM,WALL_CLOCK_BS,FOUND,NOT_FOUND" &> ${output_directory}/table_EDSBWT_search.txt

for eds_file in ${dataset_directory}/datasets/*; do

    dataset=$(basename "$eds_file")
    dataset="${dataset%.eds}"

    echo $dataset
  
    echo "Building EDS ..."

    ${tools_directory}/EDS-BWTransform.sh "$dataset_directory/datasets/$dataset.eds" example_${dataset} 0 > /dev/null 2>&1 &
    wait

    #CHECK EXECUTION
    if [ $? -eq 1 ]; then
    echo "    Command terminated by signal"
    else
    echo "    Done"
    fi

    #BUILF MLF for each dataset
    echo "Building M_LF"

        #${tools_directory}/build/build_MLF example_${dataset}_10 ${alpha}

    (/usr/bin/time -v ${tools_directory}/build/build_MLF example_${dataset} ${alpha}) &> ${output_directory}/data_${dataset}.txt


    #CHECK EXECUTION
    if [ $? -eq 0 ]; then
    echo "  Command terminated by signal"
    CPU=0
    CLOCK=0
    RAM=0
    WALL_CLOCK_BS=0
    LENGTH=0
    RUNS=0
    R_=0 
    else
    echo "    Done"
    CPU=$(grep "CPU" ${output_directory}/data_${dataset}.txt | cut -f 7 -d " ")
    CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}.txt | cut -f 8 -d " ")
    RAM=$(grep "Maximum" ${output_directory}/data_${dataset}.txt | cut -f 6 -d " ") 
    WALL_CLOCK_BS=$(grep "M_LF build took:" ${output_directory}/data_${dataset}.txt | cut -f 4 -d " ")
    LENGTH=$(grep "text length" ${output_directory}/data_${dataset}.txt | cut -f 3 -d " ")
    RUNS=$(grep "new size" ${output_directory}/data_${dataset}.txt | cut -f 4 -d " ")
    R_=$(grep "r'" ${output_directory}/data_${dataset}.txt | cut -f 6 -d " ")

    fi

    python3 ../junctions/scripts/msatoeds/get_stats.py ${eds_file} eds &> ${output_directory}/stats_${dataset}.txt
    echo "Getting stats"

    echo "    Done"
    DEG_SYMBOLS=$(grep "Number of segments:" ${output_directory}/stats_${dataset}.txt | cut -f 4 -d " ")
    STRINGS=$(grep "Number of strings:" ${output_directory}/stats_${dataset}.txt | cut -f 4 -d " ")


    #STORE INFORMATION
    echo "${dataset},${CPU},${CLOCK},${RAM},${WALL_CLOCK_BS},${LENGTH},${RUNS},${R_},${DEG_SYMBOLS},${STRINGS}" >> ${output_directory}/table_build_MLF.txt
    wait

    #rm ${output_directory}/example_${dataset}*
    rm ${output_directory}/data_${dataset}.txt

    for length in ${pattern_length[@]}; do

            echo " Searching patterns in" ${dataset}"_"${length}.txt
            (/usr/bin/time -v  ${tools_directory}/build/MOVE_EDSBWTSearch example_${dataset} "${dataset_directory}/randomPatterns/patterns_${dataset}"_"${length}.txt")  &> "${output_directory}/data_${dataset}_${length}.txt"
            wait

            #CHECK EXECUTION
            if [ $? -eq 1 ]; then
            echo "  Command terminated by signal"
            CPU=0
            CLOCK=0
            RAM=0
            WALL_CLOCK_BS=0
            FOUND=0
            NOT_FOUND=0 
            else
            echo "    Done"
            CPU=$(grep "CPU" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 7 -d " ")
            CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 8 -d " ")
            RAM=$(grep "Maximum" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 6 -d " ")
            WALL_CLOCK_BS=$(grep "bs took:" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 2 -d " ")
            FOUND=$(grep "count_found" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 3 -d " ")
            NOT_FOUND=$(grep "count_not_found" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 3 -d " ")
            fi

            #STORE INFORMATION
            echo "${dataset}"_"${length},${CPU},${CLOCK},${RAM},${WALL_CLOCK_BS},${FOUND},${NOT_FOUND}" >> ${output_directory}/table_MLF_search.txt

            #cp ${output_directory}/data_${dataset}"_"${length}.txt ${output_dirStdERROUT}/data_RLO_${dataset}.txt
            rm "${output_directory}/data_${dataset}_${length}.txt"
            #rm ${output_directory}/${dataset}"_"${length}.txtoutput_M_LF.csv

    done

    for length in ${pattern_length[@]}; do

        echo " Searching patterns in" ${dataset}"_"${length}.txt
        echo "$patterns_${dataset}"_"${length}.txt"
        (/usr/bin/time -v  ${tools_directory}/EDSBWTsearch example_${dataset} "${dataset_directory}/randomPatterns/patterns_${dataset}"_"${length}.txt")  &> "${output_directory}/data_EDSBWT_${dataset}_${length}.txt"
        #${tools_directory}/build/MOVE_EDSBWTSearch example_${dataset}_10 ${dataset_directory}/randomPatterns/${dataset}"_"${length}.txt &
        wait

        #CHECK EXECUTION
        if [ $? -eq 1 ]; then
        echo "  Command terminated by signal"
        CPU=0
        CLOCK=0
        RAM=0
        WALL_CLOCK_BS=0
        FOUND=0
        NOT_FOUND=0 
        else
        echo "    Done"
        CPU=$(grep "CPU" ${output_directory}/data_EDSBWT_${dataset}"_"${length}.txt | cut -f 7 -d " ")
        CLOCK=$(grep "Elapsed" ${output_directory}/data_EDSBWT_${dataset}"_"${length}.txt | cut -f 8 -d " ")
        RAM=$(grep "Maximum" ${output_directory}/data_EDSBWT_${dataset}"_"${length}.txt | cut -f 6 -d " ")
        WALL_CLOCK_BS=$(grep "bs took:" ${output_directory}/data_EDSBWT_${dataset}"_"${length}.txt | cut -f 2 -d " ")
        FOUND=$(grep "count_found" ${output_directory}/data_EDSBWT_${dataset}"_"${length}.txt | cut -f 3 -d " ")
        NOT_FOUND=$(grep "count_not_found" ${output_directory}/data_EDSBWT_${dataset}"_"${length}.txt | cut -f 3 -d " ")
        fi

        #STORE INFORMATION
        echo "${dataset}"_"${length},${CPU},${CLOCK},${RAM},${WALL_CLOCK_BS},${FOUND},${NOT_FOUND}" >> ${output_directory}/table_EDSBWT_search.txt

        #cp ${output_directory}/data_${dataset}"_"${length}.txt ${output_dirStdERROUT}/data_RLO_${dataset}.txt
        rm "${output_directory}/data_EDSBWT_${dataset}_${length}.txt"
        #rm ${output_directory}/${dataset}"_"${length}.txtoutput_M_LF.csv

    done

    rm example_*

done