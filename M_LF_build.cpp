#include <iostream>
#include <assert.h>

#include "Parameters.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <move_r/move_r.hpp>


using namespace std;



//CONTROLLARE TIPATURE VARIE
#define sym_t char;
using i_sym_t = constexpr_switch_t<
        constexpr_case<sizeof(char) == 1,    uint8_t>,
        constexpr_case<sizeof(char) == 2,    uint16_t>,
        constexpr_case<sizeof(char) == 4,    uint32_t>,
     /* constexpr_case<sizeof(sym_t) == 8, */ uint64_t
    >;
std::vector<std::pair<uint32_t,uint32_t>> I_LF;
using rsl_t = rank_select_support<i_sym_t,uint32_t,true,true>; // type of RS_L'
rsl_t _RS_L_;


int main(int argc, char *argv[]){
    
	if(argc != 2 ) {
		std::cerr << "usage: " << argv[0] << " input" << std::endl;
		exit(1);
	}

	
	if(strcmp(argv[1], "--help")==0 || strcmp(argv[1], "-help")==0 || strcmp(argv[1], "-h")==0) {
		std::cerr << "usage: " << argv[0] << " input" << std::endl;
		exit(1);
	}
	
	string inputName = argv[1];	

	string runsFileName = string(inputName) + "_runs.txt";
	FILE *runsFile = fopen(runsFileName.c_str(), "r");
	if (runsFile == NULL) {
		std::cerr << "Error opening \"" << runsFile << "\" file"<< std::endl;
		exit (1);
	}

	int initPos;
	int LFPos;
	int I_LF_index=0;
	no_init_resize(I_LF,33);

	while (fscanf(runsFile, "%d,%d\n", &initPos, &LFPos) != EOF) {
		std::cout<<initPos<<" "<<LFPos<<std::endl;

	 	/* Write the pair (i',LF(i')) to the next position i in I_LF, where
		LF(i') = C[L[i']] + rank(L,L[i'],i'-1) = C[p'][L[i']] + C[i_p][L[i']]. */
		I_LF[I_LF_index++]=std::make_pair(initPos,LFPos);
		/* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i'-1),
            for each c in [0..255] 
            C[i_p][run_sym(i_p,i)] += run_len(i_p,i);*/
    }
	 move_data_structure_l_<> M_LF(I_LF,33,{//33 DA CAMBIARE E METTERE IN BASE ALLA LUNGHEZZA DI T
        .num_threads = 4, .a = 2
    }, 8);
	//: (uint8_t)(std::ceil(std::log2(idx.sigma+1)/(double)8)*8) se non usiamo un byte alphabet

    uint32_t r_=M_LF.num_intervals();

    for (uint32_t i=0; i<r_; i++) {
        std::cout << to_string<>({M_LF.p(i),M_LF.q(i)})<<std::endl;
    }

	string runsAuxFileName = string(inputName) + "_runs.aux";
	FILE *runsAuxFile = fopen(runsAuxFileName.c_str(), "r");
	if (runsAuxFile == NULL) {
		std::cerr << "Error opening \"" << runsAuxFile << "\" file"<< std::endl;
		exit (1);
	}

	int i=0;
	int initPos2;
    dataTypedimAlpha let2;

	while (fscanf(runsAuxFile, "%d,%c\n", &initPos2, &let2) != EOF) {
		if (let2==EMPTY_CHAR){
			M_LF.set_L_(i,TERMINATE_CHAR);
		}
		else{
			M_LF.set_L_(i,let2);
		}
		i++;
    }

    // print the pairs of the resulting disjoint interval sequence
    /*for (uint32_t i=0; i<r_; i++) {
        std::cout << M_LF.L_(i)<<std::endl;
    }*/
	//std::function<char(uint32_t)> read =[&M_LF](uint32_t i){return M_LF.L_(i);};

	/*_RS_L_ = rsl_t(read,SIZE_ALPHA,(uint32_t)0,r_-1);
	std::cout<<_RS_L_.frequency('T')<<std::cout;*/
	//std::string ciao="ciao";
	/*_RS_L_ = rsl_t(ciao);
	_RS_L_.contains('c');*/
	/*
    std::pair<uint32_t,uint32_t> ix1{5,0};
    std::cout<<""<<std::endl;
    std::cout << to_string<>(ix1 = mds_str.move(ix1)) << std::endl;
    std::cout << to_string<>(ix1 = mds_str.move(ix1)) << std::endl;
	*/
	return 1;
}