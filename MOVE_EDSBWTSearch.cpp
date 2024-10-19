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

bool compare(rangeElement a, rangeElement b);

MOVE_EDSBWT::MOVE_EDSBWT(string inputFileName, string filepatterns){
	
	std::cerr << "Backward Search\n";
	ext = ".aux";
	
	cout << "DEBUG: " << DEBUG << endl;

	num_of_eof=0;

	recoverInfo(inputFileName);

	M_LF_Dollar_Input_interval = new uint32_t[num_of_eof];
	
	n=lengthTot_plus_eof;//text length

	#if RECOVER_INFO
		#if DEBUG==1
			fprintf(stderr, "##BEFORE recoverInfo\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
		#endif
	#endif

	cout << "RETRIEVING M_LF..."<< endl;
	if (retrieve_MLF(inputFileName)==-1){
		cout<<"Error retrievinf MLF, exit"<<endl;
		exit(0);
	}
	cout << "M_LF RETRIEVED"<< endl;

	std::function<uint32_t(uint32_t)> read = [this](uint32_t i){return M_LF.L_(i);};//function that given an index i, returns M_LF.L_(i)
	//useful for building rank-select data structure
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


	char *fileBitVector = new char[256];
	sprintf (fileBitVector,"%s%s",inputFileName.c_str(),".bitvector");
	
	#if DEBUG == 1
	cout<<"The bitvector is in: "<<fileBitVector<<endl;
	fprintf(stderr, "##BEFORE rank_1_type\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	#endif
	
	if (load_from_file(rrrb,fileBitVector) != 1) {
		std::cerr << "Error opening \"" << fileBitVector << "\" file"<< std::endl;
		exit (1);
	}
	delete[] fileBitVector;

	std::function<int(uint32_t)> read_bitvector = [this](uint32_t i){return (int)(rrrb[i]);};//function that given an index i, rrrb[i]
	//IMPORTANTE SI POTREBE LIBERARE RRRB GIÀ QUI. (o forse no)

	rrrb_size = rrrb.size();
	_RS_bitvector_=rsl_t(read_bitvector,0,rrrb_size); //mi pare non trovi gli 0
	cout<<"BitVector size: "<<_RS_bitvector_.size()<<endl;

	//rank_support_v<1> rb_1(&rrrb);
	
	#if DEBUG == 1
		fprintf(stderr, "##BEFORE select_1_type\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	#endif
	bit_vector::select_1_type bsel_1(&rrrb);
	first_symbol_index = bsel_1(2)-1; //IMPORTANTE CAMBIARE POI


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
		
		if(backwardSearch(inputFileName.c_str(), filepatterns.c_str(), i+1, kmer, lenKmer) > 0){
			//std::cerr << "1" << endl;
			count_found++;
			std::cerr<<"OCCORRENZA DI: "<<kmer<<" TROVATA"<<endl;
		}
		else{
			//std::cerr << "0" << endl;
			count_not_found++;
			std::cerr<<"OCCORRENZA DI: "<<kmer<<" NON TROVATA"<<endl;
		}
		
		#if DEBUG==1
			time (&endI);
			difI = difftime (endI,startI);
			std::cerr << "End backwardSearch " << endI << " seconds\n";
			std::cerr << "backwardSearch tooks " << difI << " seconds\n";
		#endif
		i++;
	}

	fprintf(stderr, "##\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());

	InFileKmer.close();
	
	std::cerr << endl;
	std::cerr << "count_found = " << count_found << endl;
	std::cerr << "count_not_found = " << count_not_found << endl;


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
	delete[] EOF_ID_Copy;
	delete[] M_LF_Dollar_Input_interval;

	//IMPORTANTE LIBERARE TUTTO
}

int MOVE_EDSBWT::retrieve_MLF(std::string inputFileName){

	no_init_resize(I_LF,n);

	string input_MLF_Name = string(inputFileName) + "_MLF.aux";
    std::ifstream ifs_MLF(input_MLF_Name);
    if (!ifs_MLF) {
        std::cerr << "Error opening file: " << input_MLF_Name << std::endl;
        return -1;
    }

	M_LF.load(ifs_MLF);

	cout<<"size= "<<I_LF.size()<<endl;
	r_=M_LF.num_intervals();


	I_LF.resize(r_);
	cout<<"new size= "<<I_LF.size()<<endl;
	cout<<"actual size= "<<r_<<endl;

	ifs_MLF.close();

	//build M_LF__Dollar_Input_Interval
	M_LF_Dollar_Input_interval[0]=0;
	for (uint32_t i=1;i<num_of_eof;i++){
		pos_t x=findInputInterval(i);
		M_LF_Dollar_Input_interval[i]=x;
	}

	return 1;
}


pos_t MOVE_EDSBWT::findInputInterval(uint32_t i){
	pos_t x = 0;
	while (i >= M_LF.p(x)) {
		x++;
	}
	return x-1;
}

void MOVE_EDSBWT::init_backward_search(rangeElement& firstInterval){//check dei parametri qui
	firstInterval.startPosN = 0; //b
	firstInterval.endPosN = n-1; //e
	firstInterval.startPosNII = 0; //b', input interval of b
	firstInterval.endPosNII = r_-1; //e', input interval of e
}


int MOVE_EDSBWT::backwardSearch(std::string fileInput, string fileOutDecode, dataTypeNSeq n_kmer, std::string kmer, dataTypelenSeq lenKmer){ 

	bool res;
	bool found_one_occ;

#if DEBUG==1
	std::cerr << "Pattern: " << kmer << "\n";
#endif	

	//Initialization
	uchar symbol = kmer[lenKmer-1];
	cout <<symbol<<endl;
	vectRangeOtherPile.resize(1);

	init_backward_search(vectRangeOtherPile[0]);

	/*
	uint32_t b_=_RS_L_.rank('A',0);
	b_=_RS_L_.select('A',b_+1);
	cout<<"positionII of first A "<<b_<<endl;
	uint32_t b = M_LF.p(b_);
	uint32_t e_=_RS_L_.rank('A',r_);
	cout<<"number of A end"<<e_<<endl;
	e_=_RS_L_.select('A',e_);
	uint32_t e=M_LF.p(e_+1)-1;
	cout<<"positionII of last A "<<e_<<endl;
	M_LF.move(b,b_);
		cout<<"new position of first A "<<b<<endl;
	M_LF.move(e,e_);
		cout<<"new position of last A "<<e<<endl;

	b_=_RS_L_.rank('#',0);
	b_=_RS_L_.select('#',b_+1);
	cout<<"positionII of first # "<<b_<<endl;
	b = M_LF.p(b_);
	e_=_RS_L_.rank('#',r_);
	cout<<"number of # end"<<e_<<endl;
	e_=_RS_L_.select('#',e_);
	cout<<"positionII of last # "<<e_<<endl;
	e=M_LF.p(e_+1)-1;
	M_LF.move(b,b_);
		cout<<"new position of first # "<<b<<endl;
	M_LF.move(e,e_);
		cout<<"new position of last # "<<e<<endl;


	b_=_RS_L_.rank('C',0);
	b_=_RS_L_.select('C',b_+1);
	cout<<"positionII of first C "<<b_<<endl;
	b = M_LF.p(b_);
	e_=_RS_L_.rank('C',r_);
	cout<<"number of C end"<<e_<<endl;
	e_=_RS_L_.select('C',e_);
	cout<<"positionII of last C "<<e_<<endl;
		e=M_LF.p(e_+1)-1;

	M_LF.move(b,b_);
		cout<<"new position of first C "<<b<<endl;
	M_LF.move(e,e_);
		cout<<"new position of last C "<<e<<endl;

	b_=_RS_L_.rank('G',0);
	b_=_RS_L_.select('G',b_+1);
	cout<<"positionII of first G "<<b_<<endl;
	b = M_LF.p(b_);
	e_=_RS_L_.rank('G',r_);
	cout<<"number of G end"<<e_<<endl;
	e_=_RS_L_.select('G',e_);
	cout<<"positionII of last G "<<e_<<endl;
		e=M_LF.p(e_+1)-1;

	M_LF.move(b,b_);
		cout<<"new position of first G "<<b<<endl;
	M_LF.move(e,e_);
		cout<<"new position of last G "<<e<<endl;

	b_=_RS_L_.rank('T',0);
	b_=_RS_L_.select('T',b_+1);
	cout<<"positionII of first T "<<b_<<endl;
	b = M_LF.p(b_);
	e_=_RS_L_.rank('T',r_);
	cout<<"number of T end"<<e_<<endl;
	e_=_RS_L_.select('T',e_);
	cout<<"positionII of last T "<<e_<<endl;
		e=M_LF.p(e_+1)-1;

	M_LF.move(b,b_);
		cout<<"new position of first T "<<b<<endl;
	M_LF.move(e,e_);
		cout<<"new position of last T "<<e<<endl;
	
	*/


	res=backward_search_step(symbol, vectRangeOtherPile);
	int link_res;

	if (!res){
		std::cout<<"pattern non presente"<<std::endl;
		return 0;
	}
	for (dataTypelenSeq posSymb=lenKmer-1; posSymb>0; posSymb--) {   //For each symbol of the kmer

		found_one_occ=false;

		link_res = link();
		if (link_res!=1){
			std::cerr << "Error link"<< std::endl;
        	exit (EXIT_FAILURE);
		}

		/*cerr<<"backwardSearch - dopo link (including merge)"<<"\n";
		print_interval_number();
		cerr << "vectRangeDollarPile (computed): ";
		print(vectRangeDollarPile);		
		cerr << "-symbPile: " << (int)symbPile << "\n";
		cerr << "vectRangeOtherPile (original): ";
		print(vectRangeOtherPile);	*/

		symbol= kmer[posSymb-1];

		if (!vectRangeDollarPile.empty()){
			res = backward_search_step (symbol,vectRangeDollarPile);
			if (!res){
			std::cout<<"pattern non presente in vectRangeDollarPile"<<std::endl;
			}
			else found_one_occ=true;
		}

		res=backward_search_step(symbol,vectRangeOtherPile);
		if (!res){
			std::cout<<"pattern non presente in vectRangeOtherPile"<<std::endl;
		}
		else {
			found_one_occ=true;
		}

		if (!found_one_occ){
			return 0;
		}

		//Append elements of vectRangeOtherPile to vectRangeDollarPile
		vectRangeDollarPile.insert(std::end(vectRangeDollarPile), std::begin(vectRangeOtherPile), std::end(vectRangeOtherPile));
		
		vectRangeDollarPile.shrink_to_fit();  //Requests the container to reduce its capacity to fit its size.
		vectRangeOtherPile.swap(vectRangeDollarPile);
		vector<rangeElement>().swap(vectRangeDollarPile);   // clear output reallocating 

		std::sort(vectRangeOtherPile.begin(),vectRangeOtherPile.end(),compare);

		//print(vectRangeOtherPile);

		//merge part 
		if (vectRangeOtherPile.size()>1){		
		for (uint32_t i=1;i<vectRangeOtherPile.size();i++){
			if ( (vectRangeOtherPile[i].startPosN == vectRangeOtherPile[i-1].endPosN +1)){
			vectRangeOtherPile[i-1].endPosN = vectRangeOtherPile[i].endPosN;
			vectRangeOtherPile[i-1].endPosNII = vectRangeOtherPile[i].endPosNII;
			vectRangeOtherPile.erase(vectRangeOtherPile.begin()+i);
			if (vectRangeOtherPile.size()>0){
				i--;
					}
				}
			}
		}

		//print(vectRangeOtherPile);

		#if DEBUG == 1
		std::cerr << "\n Iteration: posSymb in Pattern = "  << (int) posSymb << " symbol "<< symbol  << "\t";
		#endif
	}
		//assert(link(rb_1,bsel_1)==1);

		//For each symbol in the kmer we have to update both vectRangeDollarPile (if not empty) and vectRangeOtherPile 

//		std::cerr << "--------------------------------------backward search: the cycle for " << "symbol in position " << (int)posSymb << " took " << difftime(end,start) << " seconds\n\n";

	return 1;
}

bool MOVE_EDSBWT::backward_search_step(sym_t sym, std::vector<rangeElement>& vectRange){

	dataTypeNSeq sizeVectRange = vectRange.size();

	bool res1=false; //da cambiare poi

	cout<<"size:" <<sizeVectRange<<endl;

	for (uint32_t k=0; k < sizeVectRange; k++) {	
		cout<<vectRange[k].startPosN<<endl;
		bool res = updateSingleInterval(sym,vectRange[k]);
		//cout<<"Lettera "<<sym<<"presente nell'intervallo "<	<vectRange[k].startPosN<<" "<<vectRange[k].endPosN<<endl;
		if (!res) {
			cout<<"Lettera "<<sym<<"non presente nell'intervallo "<<vectRange[k].startPosN<<" "<<vectRange[k].endPosN<<endl; //direi di fare remove qui
			//cout<<"rimozione intervallo"<<endl;
			vectRange.erase(vectRange.begin()+k);
			if (vectRange.size()==0 || k>=vectRange.size()){
				break;
			}
				k--;
		}
		else res1=true;
	}

	return res1;

}

int MOVE_EDSBWT::updateSingleInterval(sym_t sym, rangeElement& interval){

	dataTypeNChar b,e,b_,e_;
	b = interval.startPosN;	
	e = interval.endPosN;
	b_ = interval.startPosNII;	
	e_ = interval.endPosNII;
	cout<<" start b= "<<b<<" e= "<<e<<" b'= "<<b_<<" e'="<<e_<<endl;	

   try{
	if (!_RS_L_.contains(sym)) return false;

	if (sym != M_LF.L_(b_)){
		b_ = _RS_L_.rank(sym,b_);
		if (b_ == _RS_L_.frequency(sym)) return false;
		b_ = _RS_L_.select(sym,b_+1);
		if (b_ > e_) return false;
		b = M_LF.p(b_);
	}

    // Find the lexicographically largest suffix in the current suffix array interval that is prefixed by P[i]
    if (sym != M_LF.L_(e_)) {
		e_ = _RS_L_.select(sym,_RS_L_.rank(sym,e_));
		e = M_LF.p(e_+1)-1;
	}   
	//cout<<"b"<<b<<"e"<<e<<"b'"<<b_<<"e'"<<e_<<endl;	
	

    // If the suffix array interval is empty, P does not occur in T, so return false.
    if (b > e) return false;
    
    /* Else, set b <- LF(b) and e <- LF(e). The following two optimizations increase query throughput slightly
        if there are only few occurrences */
    if (b_ == e_) {
        if (b == e) {
            /* If \hat{b'}_i == \hat{e'}_i and b'_i = e'_i, then computing
            (e_i,\hat{e}_i) <- M_LF.move(e'_i,\hat{e'}_i) is redundant */
            M_LF.move(b,b_);
            e = b;
            e_ = b_;

			interval.startPosN=b;	
			interval.endPosN=e;
			interval.startPosNII=b_;	
			interval.endPosNII=e_;	
			cout<<" end b= "<<b<<" e= "<<e<<" b'= "<<b_<<" e'="<<e_<<endl;	

        } else {
            /* If \hat{b'}_i == \hat{e'}_i, but b'_i != e'_i, then e_i = b_i + e'_i - b'_i and therefore
            \hat{b'}_i < \hat{e'}_i, hence we can compute \hat{e'}_i by setting e_ <- \hat{b'}_i = b_ and
            incrementing e_ until e < M_LF.p[e_+1] holds; This takes O(a) time because of the a-balancedness property */
            pos_t diff_eb = e - b;
            M_LF.move(b,b_);
            e = b + diff_eb;
            e_ = b_;
            
            while (e >= M_LF.p(e_+1)) {
                e_++;
            }
			interval.startPosN=b;	
			interval.endPosN=e;
			interval.startPosNII=b_;	
			interval.endPosNII=e_;	
			cout<<" end b= "<<b<<" e= "<<e<<" b'= "<<b_<<" e'="<<e_<<endl;	
        }

    } else {
        M_LF.move(b,b_);
        M_LF.move(e,e_);
		interval.startPosN=b;	
		interval.endPosN=e;
		interval.startPosNII=b_;	
		interval.endPosNII=e_;	
		cout<<" end b= "<<b<<" e= "<<e<<" b'= "<<b_<<" e'="<<e_<<endl;	
    }
   } catch (const std::bad_optional_access& e) {
        // Handle the exception
        std::cerr << "Caught exception: " << e.what() << std::endl;
		std::cerr << "Pattern Non presente"<< std::endl;
		return false;
    }

    return true;
}

int MOVE_EDSBWT::link(){

	// find the indexes of each dollar contained in each interval [startPosN,endPosN] of vectRangeOtherPile and store it in d.
	deque<dataTypeNSeq> d;
	dataTypeNSeq sizeVectRange = vectRangeOtherPile.size();

	dataTypeNChar b_,e_;

	for (dataTypeNSeq itVectRange=0; itVectRange < sizeVectRange; itVectRange++) {

		b_ = vectRangeOtherPile[itVectRange].startPosNII;	
		e_ = vectRangeOtherPile[itVectRange].endPosNII;
		//cout<<"start"<<vectRangeOtherPile[itVectRange].startPosN<<" end"<<vectRangeOtherPile[itVectRange].endPosN<<endl;


		if(b_ <= e_){
			dollars_in_interval(d,b_,e_);
		}
	}

	//after obtaining the vector of the indexes (note that there are no repetitions), sort it.
	std::sort(d.begin(),d.end());
	//cout<<d.size()<<endl;
	//cout<<"num of eof"<<num_of_eof<<endl;


	//TO DO: limitation
	vectRangeDollarPile.reserve(d.size()*2/3);

	rangeElement current;


	while(!(d.empty())){

		//cout<<d.back()<<endl;
		current = preceding_dollars_finder(d.back());
		d.pop_back();

		if (!(vectRangeDollarPile.empty()) and (current.endPosN + 1 >= (vectRangeDollarPile.back()).startPosN)){
			if(current.startPosN < (vectRangeDollarPile.back()).startPosN){	
			
				(vectRangeDollarPile.back()).startPosN = current.startPosN;
				(vectRangeDollarPile.back()).startPosNII = current.startPosNII;
			}
		}
		else {
			vectRangeDollarPile.insert(vectRangeDollarPile.end(),current);
		}
		cout<<"start "<<current.startPosN<<" end "<<current.endPosN<<" startII "<<current.startPosNII<<" endII "<<current.endPosNII<<endl;
		dollars_in_interval(d,current.startPosNII,current.endPosNII);
		//d.shrink_to_fit();
	}

	reverse(vectRangeDollarPile.begin(),vectRangeDollarPile.end());
	return 1;
}

bool compare(rangeElement a, rangeElement b) {
    return a.startPosN < b.startPosN; // Sorts in ascending order
}


rangeElement MOVE_EDSBWT::preceding_dollars_finder(dataTypeNSeq i){
	
	dataTypeNChar a;

	//a = rb_1(i+1); //a is the index of the segment containing the i-th dollar
	if (i+1>rrrb_size){
		a=_RS_bitvector_.rank(1,rrrb_size);
	}
	else{
		a=_RS_bitvector_.rank(1,i+1);
	}

	cout<<"# dollars before"<<i<<" = "<<a<<endl;
	dataTypeNChar start=_RS_bitvector_.select(1,a-1);

	//dataTypeNChar start = bsel_1(a-1); //bsel_1(k) gives the index of the k-th 1. That is, the index of the first word of the k-th segment.
	
	
	//dataTypeNChar end = bsel_1(a)-1; //same as before. We subtract 1 because we want to find the index of the last word of the (a-1)-th segment.
	dataTypeNChar end =_RS_bitvector_.select(1,a)-1;
		//cout<<"start"<<start<<" end"<<end<<endl;

	#if DEBUG == 1
	std::cerr << "preceding_dollars_finder - start " << start << " end " << end << "\n";
	#endif

	std::cerr << "preceding_dollars_finder - start " << start << " end " << end << "\n";

	//IMPORTANTE SORT NECESSARIA?
	rangeElement output;
	output.startPosN = start;//provato a togliere il +1 a start e end, va fatto così
	output.startPosNII = M_LF_Dollar_Input_interval[start];//?
	output.endPosN = end;
	output.endPosNII = M_LF_Dollar_Input_interval[end];//?
										
	return output;
}

void MOVE_EDSBWT::dollars_in_interval(deque<dataTypeNSeq> &d_out,dataTypeNChar i,dataTypeNChar j){

	dataTypeNSeq l = _RS_L_.rank('#',i);	
	dataTypeNSeq u = _RS_L_.rank('#',j+1);
	dataTypeNSeq index;
	//cout<<"ci sono "<<u-l<<" dollari nell'intervallo l="<<i<<", u="<<j<<endl;

	for (dataTypeNSeq k = l; k < u; k++){

		index = EOF_ID_Copy[k];
		//cout<<index<<endl;

		if(index > first_symbol_index){	
			if (std::find(d_out.begin(), d_out.end(), (dataTypeNSeq) index) != d_out.end()){//occhio che d contiene ushort (index sono uint invece)
			cout<<"dollar already in pile, skip"<<endl;
			return;
			}
			d_out.push_back(index);	
		}
	}
	
}


int MOVE_EDSBWT::recoverInfo(string filename) {         
    
	//Read FileInfo
    string fnInfoFile = filename + "_info" + ext;
    FILE* InfoFile = fopen(fnInfoFile.c_str(), "rb");
    if (InfoFile==NULL) {
        std::cerr << "Error opening " << fnInfoFile << "." << std::endl;
        exit (EXIT_FAILURE);
    }
	//Read BWT-length, num. sequences and sizeAlpha
    int res = fread(&lengthTot_plus_eof,sizeof(dataTypeNChar),1,InfoFile); 
	if (res!=1){
		std::cerr << "Error reading lengthTot_plus_eof" << InfoFile << "." << std::endl;
        exit (EXIT_FAILURE);
	}
    res = fread(&nText,sizeof(dataTypeNSeq),1,InfoFile);
	if (res!=1){
		std::cerr << "Error reading nText" << InfoFile << "." << std::endl;
        exit (EXIT_FAILURE);
	}
    res = fread(&sizeAlpha,sizeof(dataTypedimAlpha),1,InfoFile);
	if (res!=1){
		std::cerr << "Error reading sizeAlpha" << InfoFile << "." << std::endl;
        exit (EXIT_FAILURE);
	}
	//set alpha and alphaInverse
	alphaInverse = new dataTypedimAlpha[sizeAlpha];
    for (dataTypedimAlpha i = 0; i < sizeAlpha; ++i){
		res = fread(&(alphaInverse[i]), sizeof(dataTypedimAlpha), 1, InfoFile);
		if (res!=1){
		std::cerr << "Error reading alphaInverse" << InfoFile << "." << std::endl;
        exit (EXIT_FAILURE);
	}
		alpha[alphaInverse[i]] = i;
	}
	
	std::cout << "\nFrom " << fnInfoFile << " file:\n";
    std::cout << "\tNumber of sequences: " << nText << "\n";
    std::cout << "\tTotal length (with $): " << lengthTot_plus_eof << "\n";
    std::cout << "\tSize alpha: " << (int)sizeAlpha << "\n";
	std::cout << "\tAlphabet: ";
	for (dataTypedimAlpha i = 0; i < sizeAlpha; ++i)
		std::cout << (char)alphaInverse[i] << "\t";
	std::cout << "\n";
	
	//Upload partial bitvectors with EOF positions
	EOFpos.resize(sizeAlpha);
	EOFpos_rank.resize(sizeAlpha);
	
	for(dataTypedimAlpha j=0; j<sizeAlpha; j++){
		char *fileBitVector = new char[256];
		sprintf (fileBitVector,"%s%s%u%s",filename.c_str(),"_bv_",j,ext.c_str());
		
		#if DEBUG==1
			cerr<<"Uploaded bitvector "<< fileBitVector << endl;
		#endif

		if (load_from_file(EOFpos[j],fileBitVector) != 1) {
			std::cerr << "Error loading bitvector from " << fileBitVector << std::endl;
			exit (1);
		}
		delete [] fileBitVector;
		rrr_vector<>::rank_1_type tmp(&EOFpos[j]);
		EOFpos_rank.insert(EOFpos_rank.begin()+j,tmp);
	}
	#if DEBUG==1
		fprintf(stderr, "##AFTER uploading bitvectors with EOF pos\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	#endif
	//for(dataTypedimAlpha j=0; j<sizeAlpha; j++){
	//	cerr << EOFpos_rank[j].size() << endl;
	//}
	
	//Read dollars ID from FileInfo
	//set EOF_ID
    EOF_ID = new dataTypeNChar*[sizeAlpha];
	dataTypeNChar *numEOF = new dataTypeNChar[sizeAlpha];
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
		numEOF[j] = EOFpos_rank[j](EOFpos_rank[j].size());
		//cout << numEOF[j] << endl;
        EOF_ID[j] = new dataTypeNChar[numEOF[j]];
    }

	for(dataTypedimAlpha j=0; j<sizeAlpha; j++){
		for (dataTypeNChar h = 0 ; h < numEOF[j]; h++) {
        }
	}
	
	for(dataTypedimAlpha j=0; j<sizeAlpha; j++){
		for (dataTypeNChar h = 0 ; h < numEOF[j]; h++) {
            res = fread(&EOF_ID[j][h],sizeof(dataTypeNSeq),1,InfoFile);
			if (res!=1){
			std::cerr << "Error reading EOF_ID" << InfoFile << "." << std::endl;
        	exit (EXIT_FAILURE);
			}
			num_of_eof++;
        }
	}
	cout<<"NUM OF EOF"<<num_of_eof<<endl;
	#if DEBUG==1
		fprintf(stderr, "##AFTER uploading EOF_ID\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());
	#endif

	int z=0;
	EOF_ID_Copy = new dataTypeNChar[num_of_eof];
	for(dataTypedimAlpha j=0; j<sizeAlpha; j++){//copy into eof_ID_copy
		for (dataTypeNChar h = 0 ; h < numEOF[j]; h++) {
			memcpy(&EOF_ID_Copy[z++],&EOF_ID[j][h], sizeof(dataTypeNChar));
        }
	}

	//free(numEOF);
	delete[] numEOF;
	
    //Read tableOcc from FileInfo
	//set tableOcc
    tableOcc = new dataTypeNChar*[sizeAlpha];
    //Counting for each pile, es. $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        tableOcc[j] = new dataTypeNChar[sizeAlpha];
    }
	for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        //Update tableOcc with the symbols of the previous partial BWT files
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++) {
            res = fread(&tableOcc[j][h],sizeof(dataTypeNChar),1,InfoFile);
			if (res!=1){
			std::cerr << "Error reading tableOcc" << InfoFile << "." << std::endl;
        	exit (EXIT_FAILURE);
	}
        }
    }
    fclose(InfoFile);
	
    //#if DEBUG == 1
    std::cout << "\nFrom " << fnInfoFile << " file (TableOcc):\n";
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
            std::cout << tableOcc[j][h] << "\t";
		std::cout << std::endl;
    }
    //#endif

    return 1;
}


void MOVE_EDSBWT::print(std::vector<rangeElement> &vectRange){
			std::cerr << "size = " << vectRange.size() << std::endl;
			std::cerr << "startPosN, endPosN, startII, endII: ";
			for (dataTypeNSeq g = 0 ; g < vectRange.size(); g++) {
					std::cerr << vectRange[g].startPosN  << " " << vectRange[g].endPosN<<" " << vectRange[g].startPosNII<<" " << vectRange[g].endPosNII  << "\t";
			}
			std::cerr << std::endl;
}

void MOVE_EDSBWT::print_interval_number (){
			std::cerr << "Number of intervals in vectRangeOtherPile: " << vectRangeOtherPile.size() << std::endl;			
}


MOVE_EDSBWT::~MOVE_EDSBWT()
{

}