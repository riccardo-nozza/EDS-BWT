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
using rsl_t = rank_select_support<char>; // type of RS_L'
rsl_t _RS_L_;
std::string L_;


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

	string ebwtFileName = string(inputName) + ".ebwt";
	FILE *ebwtFile = fopen(ebwtFileName.c_str(), "r");
	if (ebwtFile == NULL) {
		std::cerr << "Error opening \"" << ebwtFile << "\" file"<< std::endl;
		exit (1);
	}
	fseek(ebwtFile, 0L, SEEK_END);
	int n=ftell(ebwtFile); //retrieve text length n

	fclose(ebwtFile);


	string runsFileName = string(inputName) + "_runs.txt";
	FILE *runsFile = fopen(runsFileName.c_str(), "r");
	if (runsFile == NULL) {
		std::cerr << "Error opening \"" << runsFile << "\" file"<< std::endl;
		exit (1);
	}

	int initPos;
	int LFPos;
	int I_LF_index=0;
	no_init_resize(I_LF,n);

	while (fscanf(runsFile, "%d,%d\n", &initPos, &LFPos) != EOF) {
		std::cout<<initPos<<" "<<LFPos<<std::endl;

	 	/* Write the pair (i',LF(i')) to the next position i in I_LF, where
		LF(i') = C[L[i']] + rank(L,L[i'],i'-1) = C[p'][L[i']] + C[i_p][L[i']]. */
		I_LF[I_LF_index++]=std::make_pair(initPos,LFPos);
		/* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i'-1),
            for each c in [0..255] 
            C[i_p][run_sym(i_p,i)] += run_len(i_p,i);*/
    }
	fclose(runsFile);


	move_data_structure_l_<> M_LF(I_LF,n,{
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

	fclose(runsAuxFile);

    // print the L'
    for (uint32_t i=0; i<r_; i++) {
		L_ += M_LF.L_(i);
	}

	std::cout<<L_<<std::endl;

	_RS_L_ = rsl_t("da");
	//std::cout<<_RS_L_.frequency('T')<<std::cout;
	//va in error se gli passi una lettera sbagliata
	try {
		std::cout<<"Numero run"<<_RS_L_.select('C',1)<<endl;
    } catch (const std::bad_optional_access& e) {
        // Handle the exception
        std::cerr << "Caught exception: " << e.what() << std::endl;
		std::cerr << "Pattern Non presente"<< std::endl;
    }
	//std::cout<<"Numero input interval"<<M_LF.p(_RS_L_.select('C',1))<<endl;
	//std::cout << to_string<>(M_LF.move({M_LF.p(_RS_L_.select('C',1)),_RS_L_.select('C',1)})) << std::endl;
	/*std::cout<<"Numero run"<<_RS_L_.select(_RS_L_.rank('C',0)+1)<<endl;
	std::cout<<"Posizione"<<M_LF.p(select(_RS_L_.rank('C',0)+1))<<endl;*/

	/*int numC=_RS_L_.rank('C',r_);
	int run=_RS_L_.select('C',numC);
	std::cout<<"run number"<<run<<std::endl;
	int runPosition=M_LF.p(run+1)-1;
	std::cout << to_string<>(M_LF.move({runPosition, run})) << std::endl;
	/*std::cout<<"Contiene"<<_RS_L_.rank('C',r_)<<endl;
	std::cout<<"Numero run"<<_RS_L_.select(_RS_L_.rank('C',r_))+1<<endl;
	std::cout<<"Posizione"<<M_LF.p(select(_RS_L_.rank('C',0)+1)+1)-1<<endl;*/




	return 1;
}