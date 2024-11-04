#!/bin/bash

#./scriptTestDatasets ../dataset ../Tools ./out_table ../Tools ./out_BCRdata

dataset_directory=$1;
tools_directory=$2;
output_directory=$3;
GSUFPATH=$4
alpha=$5

###### ADD TEXT LENGTH N, NUMBER OF RUNS R,  R' 

dataset_array=("100000" "200000" "400000" "800000" "1600000")
pattern_length=("8" "16" "32" "64")
##########
#needs to change the input data types (len seqs, BWT len, ...) in Parameters.h
#dataset_array=("pdb_seqres.fa"); 
##########

echo -e "M_LF_build\ndataset,%CPU,WALL_CLOCK,RAM,n,r,r'" &> ${output_directory}/table_build_MLF.txt
echo -e "MLF_search\ndataset,%CPU,WALL_CLOCK,RAM,WALL_CLOCK_BS" &> ${output_directory}/table_MLF_search.txt
echo -e "MLF_search\ndataset,%CPU,WALL_CLOCK,RAM,WALL_CLOCK_BS" &> ${output_directory}/table_EDSBWT_search.txt

for dataset in ${dataset_array[@]}; do

  echo $dataset
  
  echo "Building EDS ..."

    ${tools_directory}/EDS-BWTransform.sh ${dataset_directory}/syntheticDatasets/${dataset}_10 example_${dataset}_10 > /dev/null 2>&1 &
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

    (/usr/bin/time -v ${tools_directory}/build/build_MLF example_${dataset}_10 ${alpha}) &> ${output_directory}/data_${dataset}.txt

    #CHECK EXECUTION
    if [ $? -eq 0 ]; then
    echo "  Command terminated by signal"
    CPU=0
    CLOCK=0
    RAM=0
    LENGTH=0
    RUNS=0
    R_=0 
    else
    echo "    Done"
    CPU=$(grep "CPU" ${output_directory}/data_${dataset}.txt | cut -f 7 -d " ")
    CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}.txt | cut -f 8 -d " ")
    RAM=$(grep "Maximum" ${output_directory}/data_${dataset}.txt | cut -f 6 -d " ") 
    LENGTH=$(grep "text length" ${output_directory}/data_${dataset}.txt | cut -f 3 -d " ")
    RUNS=$(grep "new size" ${output_directory}/data_${dataset}.txt | cut -f 3 -d " ")
    R_=$(grep "r'" ${output_directory}/data_${dataset}.txt | cut -f 2 -d " ")

    fi

    #STORE INFORMATION
    echo "${dataset},${CPU},${CLOCK},${RAM},${LENGTH},${RUNS},${R_}" >> ${output_directory}/table_build_MLF.txt

    wait

    rm ${output_directory}/data_${dataset}.txt

    for length in ${pattern_length[@]}; do

        echo " Searching patterns in" ${dataset}"_"${length}.txt
        (/usr/bin/time -v  ${tools_directory}/build/MOVE_EDSBWTSearch example_${dataset}_10 ${dataset_directory}/randomPatterns/${dataset}"_"${length}.txt)  &> ${output_directory}/data_${dataset}"_"${length}.txt
        #${tools_directory}/build/MOVE_EDSBWTSearch example_${dataset}_10 ${dataset_directory}/randomPatterns/${dataset}"_"${length}.txt &
        wait

        #CHECK EXECUTION
        if [ $? -eq 1 ]; then
        echo "  Command terminated by signal"
        CPU=0
        CLOCK=0
        RAM=0
        WALL_CLOCK_BS=0 
        else
        echo "    Done"
        CPU=$(grep "CPU" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 7 -d " ")
        CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 8 -d " ")
        RAM=$(grep "Maximum" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 6 -d " ")
        WALL_CLOCK_BS=$(grep "bs took:" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 2 -d " ")
        fi

        #STORE INFORMATION
        echo "${dataset}"_"${length},${CPU},${CLOCK},${RAM},${WALL_CLOCK_BS}" >> ${output_directory}/table_MLF_search.txt

        #cp ${output_directory}/data_${dataset}"_"${length}.txt ${output_dirStdERROUT}/data_RLO_${dataset}.txt
        rm ${output_directory}/data_${dataset}"_"${length}.txt
        #rm ${output_directory}/${dataset}"_"${length}.txtoutput_M_LF.csv

    done


    for length in ${pattern_length[@]}; do

        echo " Searching patterns in" ${dataset}"_"${length}.txt
        (/usr/bin/time -v  ${tools_directory}/EDSBWTsearch example_${dataset}_10 ${dataset_directory}/randomPatterns/${dataset}"_"${length}.txt)  &> ${output_directory}/data_${dataset}"_"${length}.txt
        #${tools_directory}/build/MOVE_EDSBWTSearch example_${dataset}_10 ${dataset_directory}/randomPatterns/${dataset}"_"${length}.txt &
        wait

        #CHECK EXECUTION
        if [ $? -eq 1 ]; then
        echo "  Command terminated by signal"
        CPU=0
        CLOCK=0
        RAM=0
        WALL_CLOCK_BS=0
        else
        echo "    Done"
        CPU=$(grep "CPU" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 7 -d " ")
        CLOCK=$(grep "Elapsed" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 8 -d " ")
        RAM=$(grep "Maximum" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 6 -d " ")
        WALL_CLOCK_BS=$(grep "bs took:" ${output_directory}/data_${dataset}"_"${length}.txt | cut -f 2 -d " ")
        fi

        #STORE INFORMATION
        echo "${dataset}"_"${length},${CPU},${CLOCK},${RAM},${WALL_CLOCK_BS}" >> ${output_directory}/table_EDSBWT_search.txt

        #cp ${output_directory}/data_${dataset}"_"${length}.txt ${output_dirStdERROUT}/data_RLO_${dataset}.txt
        rm ${output_directory}/data_${dataset}"_"${length}.txt
        #rm ${output_directory}/${dataset}"_"${length}.txtoutput_M_LF.csv

    done

    rm ${output_directory}/example_*

done