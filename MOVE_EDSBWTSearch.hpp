
#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
//#include <stdlib.h>

#include "Parameters.h" // Defines ulong and uchar.
#include "Sorting.h" // Defines ulong and uchar.
#include <string.h>

#include <move_r/move_r.hpp>


#include <deque>
/*
#include <map>
#include <sdsl/bit_vectors.hpp>


#define BUFFERSIZE 1024

#define SIZE_ALPHA 256

#define DIMBLOCK  1024     // 1048576

using std::string;
using std::vector;
using std::map;*/

//CONTROLLARE TIPATURE VARIE
typedef char sym_t;
typedef uint32_t pos_t;


using i_sym_t = constexpr_switch_t<
        constexpr_case<sizeof(char) == 1,    uint8_t>,
        constexpr_case<sizeof(char) == 2,    uint16_t>,
        constexpr_case<sizeof(char) == 4,    uint32_t>,
     /* constexpr_case<sizeof(sym_t) == 8, */ uint64_t
    >;
using rsl_t = rank_select_support<char>; // type of RS_L'

using namespace sdsl;


class MOVE_EDSBWT
{
public:

    MOVE_EDSBWT(std::string, std::string);
    ~MOVE_EDSBWT();

    int build_MLF(std::string);
    pos_t findInputInterval(uint32_t i);

    void init_backward_search(rangeElement& firstInterval);

    bool backward_search_step(sym_t sym, std::vector<rangeElement>& vectRangeOtherPile);
    int updateSingleInterval(sym_t sym, rangeElement& interval);

    //int backwardSearch(string fileInput, string fileOutDecode, dataTypeNSeq n_kmer, string kmers, dataTypelenSeq lenKmer, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);    
    int backwardSearch(std::string,std::string, dataTypeNSeq n_kmer, std::string kmers, dataTypelenSeq lenKmer,rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);

    int link(rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);

    void dollars_in_interval(std::deque<dataTypeNSeq> &d_out,dataTypeNChar i,dataTypeNChar j);
    
    int recoverInfo(std::string filename);

    rangeElement preceding_dollars_finder(dataTypeNSeq i, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);

    void MergeAndRemove (std::vector< rangeElement > &vectRange, dataTypeNSeq &k,dataTypeNSeq &k_tmp);



    void print(std::vector<rangeElement> &vectRange);

    void print_interval_number();

    /*
    string  fileOutBwt;
    string  fileOutCyc;*/
    std::string  ext;
    
	std::vector<rangeElement> vectRangeDollarPile;
	std::vector<rangeElement> vectRangeOtherPile;
	dataTypedimAlpha symbPile;
    
    //dataTypeNChar SIZEBUFFERcycFiles;
    
    //dataTypelenSeq lengthRead;    //Length of each text
    //dataTypeNChar lengthTot;   //Total length of all texts without $-symbols
    dataTypeNChar lengthTot_plus_eof;   //Total length of all texts with $-symbols

    dataTypeNSeq nText;   //number total of texts in filename1
    //dataTypeNChar freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters
    
    int numthreads;
    
    dataTypeNChar** tableOcc; //contains the number of occurrences of each symbol
    dataTypedimAlpha alpha[SIZE_ALPHA]; //Corresponding between the alphabet, the piles and tableOcc
    dataTypedimAlpha sizeAlpha;  //number of the different symbols in the input texts
    dataTypedimAlpha *alphaInverse;  //Corresponding between alpha[i] and the symbol as char
    
	//Added/modified by Lapo 06-06-2024
	bit_vector rrrb;
	dataTypeNSeq first_symbol_index;

	//Added by Vero 15_05_2024
	std::vector< rrr_vector<> > EOFpos; 
	std::vector< rrr_vector<>::rank_1_type > EOFpos_rank; 
	dataTypeNChar** EOF_ID;


    dataTypeNChar* EOF_ID_Copy;
    uint32_t* M_LF_Dollar_Input_interval;
    /*
    std::vector <bool> vectInsTexts;
    dataTypeNSeq textToInsert;
    vector <dataTypeNChar> vectSizeCurrentPile;

    vector< vector< vector<dataTypeNChar> > > vectorOcc;
    vector <dataTypeNChar> numBlocksInPartialBWT;
    
    int recoverInfo(string);

    int splitIntoPartial(string, int);
    
    int buildFreq(string);
  
    int decodeBCRmultipleReverse(string);
    
    //dataTypeNChar rankManySymbolsByVector(FILE & , dataTypeNChar *, dataTypeNChar, uchar *, uchar *, FILE &);
    
    int computeVectorOcc(string filename);
    
    int convertFromCycFileToFastaOrFastq( string );
	
	#if RECOVERBW==1
		int findMultipleDollarsBackward(std::vector<rangeElementBW> &vectRange, string fileInput,string fileOutDecode,dataTypeNSeq n_kmer,rank_support_v<1> &rb_1,bit_vector::select_1_type &bsel_1);
		int updateSingleIntervalBW(std::vector<rangeElementBW> &vectRange, FILE *InFileBWT, dataTypeNSeq k, dataTypedimAlpha currentPile, dataTypeNChar * counters, dataTypeNChar *numBlockCounter, dataTypeNChar * contInCurrentBlock, dataTypeNChar toRead, uchar *bufferBlock);
		int findBlockToReadBWT(dataTypedimAlpha currentPile, dataTypeNChar *toRead, dataTypeNChar *numBlock);
		void printElementBW(std::vector<rangeElementBW> &vectRange);
	#endif
	void indices_of_dollars_in_interval(std::ostream& searchOutput,dataTypelenSeq symbPile, dataTypeNChar i, dataTypeNChar j, dataTypeNSeq n_kmer, dataTypeNChar occPos, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1);
	
	int rankInverseManyByVector (string fileInput, dataTypeNSeq numKmersInput, uchar *toFindSymbols, uchar *bufferBlock);
	int update_Pos_Pile_Blocks(dataTypeNChar *toRead, dataTypeNChar *numBlock, dataTypedimAlpha currentPile, uchar toFindSymbol);
	int update_Pos_Pile(rangeElement *tripla);*/
	////////////////////
    
private:

     // ############ VARIABLES ############
    int n; //text length
    int r_; //r', number of runs in M_LF after balancing
    move_data_structure_l_<> M_LF;//M_LF data structure, initialized
    rsl_t _RS_L_;//Rank-select data structure for L'
    std::vector<std::pair<uint32_t,uint32_t>> I_LF;//disjoint interval sequence to build M_LF
    int num_of_eof;    
};