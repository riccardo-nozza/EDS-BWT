#include <iostream>
#include <assert.h>
#include <string.h>   
#include <sstream>
#include <stdio.h>
#include <math.h>

#include <map>

using namespace std;
using std::cout;
using std::endl;

#include "MOVE_EDSBWTSearch.hpp"
#include "Parameters.h"

int main(int argc, char *argv[]){

    if( argc != 3 ) {
		std::cerr << "usage: " << argv[0] << " inputEBWTfileName inputPATTERNfile" << std::endl;
		std::cerr << "where:" << std::endl;
		std::cerr << "  inputEBWTfile is the BWT filename without the extension .ebwt (and .ebwt.qs for the QS string)" << std::endl;
		std::cerr << "  inputPATTERNfile is the pattern file" << std::endl;
		exit(1);
    }

	std::cout << "BCR_eds: " << argv[0] << std::endl;
	std::cout << "BCR_eds: The input ebwt file is " << argv[1] << std::endl;
	std::cout << "BCR_eds: The pattern file is " << argv[2] << std::endl;

	std::string InputFileName=argv[1];
	std::string filePattern=argv[2];

    MOVE_EDSBWT *MOVE_Search;

	MOVE_Search = new MOVE_EDSBWT(InputFileName, filePattern); //mode e num threads altri param

    /*
	#if DEBUG == 1
	if (MODE == 1)
		std::cout << "unBCR_QS: MODE is 1 --> unBCR " << std::endl;
	else if (MODE == 2)
		std::cout << "unBCR_QS: MODE is 2 --> unBCR by using existing partial files" << std::endl;
	else
		std::cout << "unBCR_QS: MODE is 3 --> unBCR by using .table file" << std::endl;
	#endif*/


	//int num_threads = 1; //decidere poi
    
	#if RECOVERBW==1
		std::cerr << "\nThe csv file is ready! \n";
	#else
		std::cerr << "\nThe search is finished! \n";
	#endif

	delete MOVE_Search;

	std::cerr << "The End!\n";

    return 1;
}


