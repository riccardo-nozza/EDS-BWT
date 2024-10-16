#include <iostream>
#include <assert.h>

#include "Parameters.h"
#include "Sorting.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <move_r/move_r.hpp>
#include "MOVE_EDSBWTSearch.hpp"



using namespace std;
using namespace sdsl;


MOVE_EDSBWT::MOVE_EDSBWT(string inputFileName, string filepatterns)
{
    
	std::cerr << "Backward Search\n";

	cout << "DEBUG: " << DEBUG << endl;

	cout << "BUILDING M_LF..."<< endl;
	build_MLF(inputFileName);
	cout << "M_LF built"<< endl;

	std::function<uint32_t(uint32_t)> read = [this](uint32_t i){return M_LF.L_(i);};//function that given an index i, returns M_LF.L_(i)
	//usefuk for building rank-select data structure
	_RS_L_=rsl_t(read,0,r_);


	string searchOutput_s = filepatterns + "output_M_LF.csv";
	std::ofstream searchOutput;
	searchOutput.open(searchOutput_s,ios::out);
	if(searchOutput.is_open()){
//			searchOutput << "Pattern_#" << "\t" << "word_index" << "\t" << "segment_index" << "\t" << "WordInSeg_index" << "\t" << "position_in_string\n";
		searchOutput << "#Pat" << "\t" << "$_i" << "\t" << "D[i]" << "\t" << "S_j" << "\t" << "S_j[r]\n";
		searchOutput.close();
	}
	else{
		cerr << "ERROR opening file " << searchOutput_s << " to write output\n";
		exit(1);
	}

	dataTypeNSeq count_found=0, count_not_found=0;

	std::ifstream InFileKmer(filepatterns);
	std::string kmer; 
	dataTypeNChar i=0;
	dataTypeNChar lenKmer=0;

	while (std::getline(InFileKmer, kmer)) {		
		lenKmer = kmer.length();
		cout << "Pattern: " << kmer << " of length " << lenKmer << endl;

		#if DEBUG==1
		time_t startI,endI;
		double difI;
			time (&startI);
		#endif

		init_backward_search();
		std::cout<<b<<","<<e<<","<<","<<b<<","<<e_<<std::endl;
		
		if(backwardSearch(inputFileName.c_str(), filepatterns.c_str(), i+1, kmer, lenKmer) > 0){
			//std::cerr << "1" << endl;
			count_found++;
		}
		else{
			//std::cerr << "0" << endl;
			count_not_found++;
		}
		
		#if DEBUG==1
			time (&endI);
			difI = difftime (endI,startI);
			std::cerr << "End backwardSearch " << endI << " seconds\n";
			std::cerr << "backwardSearch tooks " << difI << " seconds\n";
		#endif
		i++;
	}

	InFileKmer.close();
	

	std::cerr << endl;
	std::cerr << "count_found = " << count_found << endl;
	std::cerr << "count_not_found = " << count_not_found << endl;

	/*
    //Free the memory
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        delete [] tableOcc[j];
        tableOcc[j] = NULL;
    }
    delete [] tableOcc;
	
	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        delete [] EOF_ID[j];
        EOF_ID[j] = NULL;
    }
    delete [] EOF_ID;

    delete[] alphaInverse;
	*/
	
	
	//va in bad_optional_access se gli passi una lettera che non è nel testo	
	/*try {
	std::cout<<"Numero run"<<_RS_L_.select('C',1)<<endl;


	//std::cout<<"Numero input interval"<<M_LF.p(_RS_L_.select('C',1))<<endl;
	//std::cout << to_string<>(M_LF.move({M_LF.p(_RS_L_.select('C',1)),_RS_L_.select('C',1)})) << std::endl;
	std::cout<<"Numero run"<<_RS_L_.select(_RS_L_.rank('C',0)+1)<<endl;
	std::cout<<"Posizione"<<M_LF.p(select(_RS_L_.rank('C',0)+1))<<endl;

	int numC=_RS_L_.rank('C',r_);
	int run=_RS_L_.select('C',numC);
	std::cout<<"run number"<<run<<std::endl;
	int runPosition=M_LF.p(run+1)-1;
	std::cout << to_string<>(M_LF.move({runPosition, run})) << std::endl;

    } catch (const std::bad_optional_access& e) {
        // Handle the exception
        std::cerr << "Caught exception: " << e.what() << std::endl;
		std::cerr << "Pattern Non presente"<< std::endl;
    }*/


}


int MOVE_EDSBWT::build_MLF(std::string inputFileName){

	string ebwtFileName = string(inputFileName) + ".ebwt";
	FILE *ebwtFile = fopen(ebwtFileName.c_str(), "r");
	if (ebwtFile == NULL) {
		std::cerr << "Error opening \"" << ebwtFile << "\" file"<< std::endl;
		exit (1);
	}
	fseek(ebwtFile, 0L, SEEK_END);
	n=ftell(ebwtFile); //retrieve text length n

	fclose(ebwtFile);


	string runsFileName = string(inputFileName) + "_runs.txt";
	FILE *runsFile = fopen(runsFileName.c_str(), "r");
	if (runsFile == NULL) {
		std::cerr << "Error opening \"" << runsFile << "\" file"<< std::endl;
		exit (1);
	}

	int LPos;
	int LFPos;
	int I_LF_index=0;
	no_init_resize(I_LF,n);

	while (fscanf(runsFile, "%d,%d\n", &LPos, &LFPos) != EOF) {
		std::cout<<LPos<<" "<<LFPos<<std::endl; //questo si può poi cancellare, ora è per controllo

	 	/* Write the pair (i',LF(i')) to the next position i in I_LF*/
		I_LF[I_LF_index++]=std::make_pair(LPos,LFPos);
		/* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i'-1),
            for each c in [0..255] 
            C[i_p][run_sym(i_p,i)] += run_len(i_p,i);*/
    }
	fclose(runsFile);


	M_LF=move_data_structure_l_<>(I_LF,n,{
        .num_threads = 4, .a = 2
    }, 8); 

    r_=M_LF.num_intervals();

    for (uint32_t i=0; i<r_; i++) {
        std::cout << to_string<>({M_LF.p(i),M_LF.q(i)})<<std::endl;
    }

	string runsAuxFileName = string(inputFileName) + "_runs.aux";
	FILE *runsAuxFile = fopen(runsAuxFileName.c_str(), "r");
	if (runsAuxFile == NULL) {
		std::cerr << "Error opening \"" << runsAuxFile << "\" file"<< std::endl;
		exit (1);
	}

	int i=0;
	int initPos;
    dataTypedimAlpha let;

	while (fscanf(runsAuxFile, "%d,%c\n", &initPos, &let) != EOF) {
		if (let==EMPTY_CHAR){
			M_LF.set_L_(i,TERMINATE_CHAR);
		}
		else{
			M_LF.set_L_(i,let);
		}
		i++;
    }

	fclose(runsAuxFile);

	return 1;
}


void MOVE_EDSBWT::init_backward_search(){//check dei parametri qui
	b = 0;
	e = n-1;
	b_ = 0;
	e_ = r_-1;
}


int MOVE_EDSBWT::backwardSearch(std::string,std::string, dataTypeNSeq n_kmer, std::string kmers, dataTypelenSeq lenKmer){

	

	return 1;
}



MOVE_EDSBWT::~MOVE_EDSBWT()
{

}