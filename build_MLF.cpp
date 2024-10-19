#include <iostream>
#include <assert.h>
#include <string.h>   
#include <sstream>
#include <stdio.h>
#include <math.h>

using namespace std;
using std::cout;
using std::endl;

#include <move_r/move_r.hpp>
#include "Parameters.h"


int build_MLF(std::string inputFileName);


int main(int argc, char *argv[]){

    if( argc != 2 ) {
		std::cerr << "usage: " << argv[0] << " inputEBWTfileName" << std::endl;
		std::cerr << "where:" << std::endl;
		std::cerr << "  inputEBWTfile is the BWT filename without the extension .ebwt (and .ebwt.qs for the QS string)" << std::endl;
        std::cerr << "  The M_LF data structure will be saved in inputEBWTfile_MLF.aux" << std::endl;
		exit(1);
    }

	std::cout << "The input ebwt file is " << argv[1] << std::endl;

	std::string InputFileName=argv[1];

    std::cout<<"BUILDING M_LF ... "<<std::endl;

    time_t startI,endI;
		double difI;
			time (&startI);
    build_MLF(InputFileName);   
    	
        time (&endI);
        difI = difftime (endI,startI);  
        

    std::cout<<"M_LF BUILT in "<<difI<<" seconds;"<<std::endl;

    return 1;
}


int build_MLF(std::string inputFileName){

	string ebwtFileName = string(inputFileName) + ".ebwt";
	FILE *ebwtFile = fopen(ebwtFileName.c_str(), "r");
	if (ebwtFile == NULL) {
		std::cerr << "Error opening \"" << ebwtFile << "\" file"<< std::endl;
		exit (1);
	}
	fseek(ebwtFile, 0L, SEEK_END);
	uint32_t n=ftell(ebwtFile); //retrieve text length n
	cout<<"text length "<<n<<endl;
	fclose(ebwtFile);


	string runsFileName = string(inputFileName) + "_runs.txt";
	FILE *runsFile = fopen(runsFileName.c_str(), "rb");
	if (runsFile == NULL) {
		std::cerr << "Error opening \"" << runsFile << "\" file"<< std::endl;
		exit (1);
	}

	uint32_t LPos;
	uint32_t LFPos;
	uint32_t I_LF_index=0;
    std::vector<std::pair<uint32_t,uint32_t>> I_LF;//disjoint interval sequence to build M_LF
	no_init_resize(I_LF,n);

	int z=0;

	while (fscanf(runsFile, "%u,%u\n", &LPos, &LFPos) != EOF) {
	 	/* Write the pair (i',LF(i')) to the next position i in I_LF*/
		I_LF[I_LF_index++]=std::make_pair(LPos,LFPos);
		//cout<<I_LF[I_LF_index-1].first<<endl;
		z++;
		/* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i'-1),
            for each c in [0..255] 
            C[i_p][run_sym(i_p,i)] += run_len(i_p,i);*/
    }

	fclose(runsFile);

    // Iterate through the vector
	cout<<"size= "<<I_LF.size()<<endl;
	I_LF.resize(z);
	cout<<"new size= "<<I_LF.size()<<endl;
	cout<<"actual size= "<<z<<endl;

    move_data_structure_l_<> M_LF;//M_LF data structure, initialized

	M_LF=move_data_structure_l_<>(I_LF,n,{
		.a = 6 //gestire a poi
    }, 8); 

    I_LF.clear();
    I_LF.shrink_to_fit();

    uint32_t r_=M_LF.num_intervals();

    for (uint32_t i=0;i<r_;i++){
        //std::cout<<M_LF.p(i)<<" "<<M_LF.q(i)<<std::endl;
    }

	cout<<"r'="<<r_<<endl;
	cout<<"n= "<<n<<endl;

	string runsAuxFileName = string(inputFileName) + "_runs.aux";
	FILE *runsAuxFile = fopen(runsAuxFileName.c_str(), "rb");
	if (runsAuxFile == NULL) {
		std::cerr << "Error opening \"" << runsAuxFile << "\" file"<< std::endl;
		exit (1);
	}

	uint64_t i=0;
	uint32_t initPos;
    dataTypedimAlpha let;

	dataTypedimAlpha prev;

	while (fscanf(runsAuxFile, "%u,%c\n", &initPos, &let) != EOF) {

		while(M_LF.p(i)<initPos){
			M_LF.set_L_(i,prev);
			//cout<<M_LF.p(i)<<" minore di "<<initPos<<endl;
			//cout<<"old letter "<< prev << " = " << M_LF.L_(i-1)<<" new letter "<< M_LF.L_(i);
			i++;
		}

		if (let==EMPTY_CHAR){
			M_LF.set_L_(i,TERMINATE_CHAR);
		}
		else{
			M_LF.set_L_(i,let);
		}

		prev=let;
		i++;
    }

	fclose(runsAuxFile);

    string outputM_LF_Name = string(inputFileName) + "_MLF.aux";
    std::ofstream ofs_MLF(outputM_LF_Name);
    if (!ofs_MLF) {
        std::cerr << "Error opening file: " << outputM_LF_Name << std::endl;
        return 1;
    }
    M_LF.serialize(ofs_MLF);

    ofs_MLF.close();

	return 1;
}