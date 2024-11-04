#include <iostream>
#include <assert.h>
#include <sys/time.h>

#include "Parameters.h"
#include "Sorting.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <move_r/move_r.hpp>
#include "MOVE_EDSBWTSearch.hpp"



using namespace std;
using namespace sdsl;

//std::vector<uint32_t> EOF_RANK;
uint32_t* EOF_RANK; 

uint32_t number_updates;
uint32_t endpos;

std::vector<unsigned char> symbols = {
        static_cast<unsigned char>('#'),
        static_cast<unsigned char>('A'),
        static_cast<unsigned char>('C'),
        static_cast<unsigned char>('G'),
        static_cast<unsigned char>('T')
    };

std::unordered_map<unsigned char, std::vector<uint32_t>> rank_results;
std::unordered_map<unsigned char, std::vector<uint32_t>> select_results;




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

	//cout << "RETRIEVING M_LF..."<< endl;
	if (retrieve_MLF(inputFileName)==-1){
		cout<<"Error retrievinf MLF, exit"<<endl;
		exit(0);
	}
	//cout << "M_LF RETRIEVED"<< endl;

	std::function<uint32_t(uint32_t)> read = [this](uint32_t i){return M_LF.L_(i);};//function that given an index i, returns M_LF.L_(i)
	//useful for building rank-select data structure
	_RS_L_=rsl_t(read,0,r_);


	for (uchar sym : symbols) {
        std::vector<uint32_t> ranks;
        for (int b_ = 0; b_ <= r_; ++b_) {
            int rank_value = _RS_L_.rank(sym, b_);
            ranks.push_back(rank_value);
        }
        rank_results[sym] = ranks; // Store the ranks for the current symbol
    }

		for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
			std::vector<uint32_t> selects;
			uchar sym = symbols[j];
			cout<<sym<<endl;
			uint32_t num_occ = 0;
			num_occ = rank_results[sym][r_];
			for (int b_ = 1; b_ <= num_occ; b_++) {
            int select_value = _RS_L_.select(sym, b_);
            selects.push_back(select_value);
        	}
        	select_results[sym] = selects; // Store the ranks for the current symbol
    	}
	
	//EOF_RANK.reserve(n);
	EOF_RANK= new uint32_t[n];

	//migliorare questa roba
	/*for(int i=0;i<r_;i++){
		uint32_t num_dollars=_RS_L_.rank('#',M_LF.p(i));
		for (uint32_t j=M_LF.p(i);j<M_LF.p(i+1);j++){
		if (num_dollars>num_of_eof){
			EOF_RANK[j]=num_of_eof;
		}
		else{
		EOF_RANK[j]=num_dollars;
		}
		}
		//cout<<i<<" "<<M_LF.p(i)<<" "<<num_dollars<<endl;

		
	}*/

	//cout<<_RS_L_.rank('#',475017)<<endl;
	//cout<<EOF_RANK[475017]<<endl;


	
	string searchOutput_s = filepatterns + "output_M_LF.csv";
	searchOutput.open(searchOutput_s,ios::out);
	if(searchOutput.is_open()){
//			searchOutput << "Pattern_#" << "\t" << "word_index" << "\t" << "segment_index" << "\t" << "WordInSeg_index" << "\t" << "position_in_string\n";
		searchOutput << "#Pat" << "\t" << "$_i" << "\t" << "D[i]" << "\t" << "S_j" << "\t" << "S_j[r] \n";
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
	first_symbol_index = bsel_1(2)-1; 


	dataTypeNSeq count_found=0, count_not_found=0;

	std::ifstream InFileKmer(filepatterns);
	std::string kmer; 
	dataTypeNChar i=0;
	dataTypeNChar lenKmer=0;

	/*#if DEBUG==1
		clock_t startI,endI;
		double difI;
			time (&startI);
		#endif*/
	const clock_t begin_time = clock();

	while (std::getline(InFileKmer, kmer)) {		
		lenKmer = kmer.length();
		//cout << "Pattern: " << kmer << " of length " << lenKmer << endl;

		/*#if DEBUG==1
		time_t startI,endI;
		double difI;
			time (&startI);
		#endif*/
		number_updates=0;
		if(backwardSearch(inputFileName.c_str(), filepatterns.c_str(), i+1, kmer, lenKmer) > 0){
			//std::cerr << "1" << endl;
			count_found++;
			//std::cerr<<"OCCORRENZA DI: "<<kmer<<" TROVATA"<<endl;
		}
		else{
			//std::cerr << "0" << endl;
			count_not_found++;
			//std::cerr<<"OCCORRENZA DI: "<<kmer<<" NON TROVATA"<<endl;
		}
		
			//std::cerr << "End backwardSearch " << endI << " seconds\n";
			//std::cerr << "backwardSearch tooks " << difI << " seconds\n";
		i++;
		//cout<<"number of updates: "<<number_updates<<endl;
	}

	/*#if DEBUG==1
			time (&endI);
			difI = difftime (endI,startI);
			long long milliseconds = (long long)(difI * 1000);
	#endif
	std::cerr << "Total time " << milliseconds<< " milliSeconds\n";		
	std::cerr << "Total time " << difI<< " Seconds\n";*/	
	std::cout << "bs took:"<<float( clock () - begin_time ) /  CLOCKS_PER_SEC;


	//fprintf(stderr, "##\nmalloc_count ### current peak: %zu\n##\n", malloc_count_peak());

	InFileKmer.close();
	searchOutput.close();
	
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
	//cout<<"new size= "<<I_LF.size()<<endl;
	//cout<<"actual size= "<<r_<<endl;

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
	pos_t x = M_LF_Dollar_Input_interval[i-1]; //great
	while (i >= M_LF.p(x)) {
		x++;
	}
	return x-1;
}

void MOVE_EDSBWT::init_backward_search(rangeElement& firstInterval){
	firstInterval.startPosN = 0; //b
	firstInterval.endPosN = n-1; //e
	firstInterval.startPosNII = 0; //b', input interval of b
	firstInterval.endPosNII = r_-1; //e', input interval of e
}


int MOVE_EDSBWT::backwardSearch(std::string fileInput, string fileOutDecode, dataTypeNSeq n_kmer, std::string kmer, dataTypelenSeq lenKmer){ 

	bool res;
	bool found_one_occ;

	uint32_t num_occ=0;
#if DEBUG==1
	//std::cerr << "Pattern: " << kmer << "\n";
#endif	

	//Initialization
	uchar symbol = kmer[lenKmer-1];
	//cout <<symbol<<endl;
	vectRangeOtherPile.resize(1);

	init_backward_search(vectRangeOtherPile[0]);

	//const clock_t begin_time = clock();
	res=backward_search_step(symbol, vectRangeOtherPile);
	int link_res;
	//std::cout << "bs took:"<<float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl;

	if (!res){
		//std::cout<<"pattern non presente"<<std::endl;
		return 0;
	}
	for (dataTypelenSeq posSymb=lenKmer-1; posSymb>0; posSymb--) {   //For each symbol of the kmer

		found_one_occ=false;
		//const clock_t ex = clock();
		link_res = link();
		//std::cout << "link took:"<<float( clock () - ex ) /  CLOCKS_PER_SEC<<endl;
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
		print(vectRangeOtherPile);*/

		symbol= kmer[posSymb-1];

		if (!vectRangeDollarPile.empty()){
			const clock_t begin_time0 = clock();
			//cout<<"dollar pile size before"<<vectRangeDollarPile.size()<<endl;
			res = backward_search_step (symbol,vectRangeDollarPile);
			//cout<<"dollar pile size after"<<vectRangeDollarPile.size()<<endl;
			//std::cout << "bs took:"<<float( clock () - begin_time0 ) /  CLOCKS_PER_SEC<<endl;
			if (!res){
			//std::cout<<"pattern non presente in vectRangeDollarPile"<<std::endl;
			}
			else found_one_occ=true;
		}

		const clock_t begin_time1 = clock();
		//cout<<"other pile size before"<<vectRangeOtherPile.size()<<endl;
		res=backward_search_step(symbol,vectRangeOtherPile);
		//cout<<"other pile size after"<<vectRangeOtherPile.size()<<endl;
		//std::cout << "bs took:"<<float( clock () - begin_time1 ) /  CLOCKS_PER_SEC<<endl;

		if (res) found_one_occ=true;

		if (!found_one_occ){
			return 0;
		}

		//Append elements of vectRangeOtherPile to vectRangeDollarPile
		vectRangeDollarPile.insert(std::end(vectRangeDollarPile), std::begin(vectRangeOtherPile), std::end(vectRangeOtherPile));
		
		vectRangeDollarPile.shrink_to_fit();  //Requests the container to reduce its capacity to fit its size.
		vectRangeOtherPile.swap(vectRangeDollarPile);
		vector<rangeElement>().swap(vectRangeDollarPile);   // clear output reallocating 

		//std::sort(vectRangeOtherPile.begin(),vectRangeOtherPile.end(),compare);

		//merge 
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
	}
	// LOCATE
	
	for (uint32_t i=0;i<vectRangeOtherPile.size();i++){
		rangeElement interval = vectRangeOtherPile[i];
		uint32_t prevII_copy = interval.startPosNII;
		uint32_t prevII = interval.startPosNII;

		uint32_t start = interval.startPosN;
		uint32_t end = interval.endPosN;

		num_occ=num_occ+end-start+1;

		for (uint32_t j=start;j<=end;j++){
			uint32_t position = j;
			uint32_t pos_in_string=0;

			if (!(M_LF.p(prevII_copy+1)-position>0)){//next input interval
				prevII_copy=prevII_copy+1;
				prevII = prevII_copy;	
			}
			while (M_LF.L_(prevII)!=TERMINATE_CHAR){
				M_LF.move(position,prevII);
				pos_in_string++;
				//cout<<position<<endl;
				//cout<<prevII<<endl;
			}	
			uint32_t num_of_dollar = _RS_L_.rank('#',prevII);
			uint32_t index = EOF_ID_Copy[num_of_dollar];	

			//uint32_t dollar_II = M_LF_Dollar_Input_interval[index];
			dataTypeNChar degenerate_symbol=_RS_bitvector_.rank(1,index+1);			
			//cout<<degenerate_symbol<<endl;
			dataTypeNChar start =_RS_bitvector_.select(1,degenerate_symbol);		
			uint32_t offset_in_symbol = index - start;
			
			//dataTypeNChar end =_RS_bitvector_.select(1,a+1)-1;

			searchOutput << n_kmer << "\t" << index << "\t" << degenerate_symbol << "\t" << offset_in_symbol << "\t" << pos_in_string << endl;
			prevII = prevII_copy;
			
				}
			}

			//cout<<"num occ "<<num_occ<<endl;

	return 1;
}

bool MOVE_EDSBWT::backward_search_step(sym_t sym, std::vector<rangeElement>& vectRange){

	dataTypeNSeq sizeVectRange = vectRange.size();

	std::vector<rangeElement> validRanges;
	//cout<<"number of intervals before update "<<vectRange.size()<<endl;

	//bool res1=false;

	//cout<<"size:" <<sizeVectRange<<endl;

	for (uint32_t k=0; k < sizeVectRange; k++) {
		//cout<<vectRange[k].startPosN<<endl;
		//const clock_t begin_time1 = clock();
		bool res = updateSingleInterval(sym,vectRange[k]);
		//cout<<vectRange[k].endPosN;
		number_updates+=1;
		//std::cout << "single interval step took:"<<float( clock () - begin_time1 ) /  CLOCKS_PER_SEC<<endl;
		//cout<<"Lettera "<<sym<<"presente nell'intervallo "<	<vectRange[k].startPosN<<" "<<vectRange[k].endPosN<<endl;
		/*if (!res){ //|| (k>0 && vectRange[k].startPosN<vectRange[k-1].startPosN)) {
			//cout<<"Lettera "<<sym<<"non presente nell'intervallo "<<vectRange[k].startPosN<<" "<<vectRange[k].endPosN<<endl; //remove qui
			//cout<<"rimozione intervallo"<<endl;
			vectRange.erase(vectRange.begin()+k);
			if (vectRange.size()==0) return false;
			if (k>=vectRange.size()){
				break;
			}
				k--;	
		}*/
		if (res){
			//res1=true;
			validRanges.push_back(vectRange[k]);	
		}
		//else res1=true;
	}
	if (vectRange.empty()) {
        return false;
    }

	vectRange.swap(validRanges);
	//vectRange.shrink_to_fit();
	//cout<<"number of intervals after update"<<vectRange.size()<<endl;

	//return res1;
	return true;

}

int MOVE_EDSBWT::updateSingleInterval(sym_t sym, rangeElement& interval){

	dataTypeNChar b,e,b_,e_;
	b = interval.startPosN;	
	e = interval.endPosN;
	b_ = interval.startPosNII;	
	e_ = interval.endPosNII;
	//cout<<" start b= "<<b<<" e= "<<e<<" b'= "<<b_<<" e'="<<e_<<endl;	

   try{
	if (!_RS_L_.contains(sym)) return false;

	if (sym != M_LF.L_(b_)){
		//b_ = _RS_L_.rank(sym,b_);
		b_=rank_results[sym][b_];
		//if (b_ == _RS_L_.frequency(sym)) return false;
		if (b_ ==rank_results[sym][r_]) return false;
		//b_ = _RS_L_.select(sym,b_+1);
		b_=select_results[sym][b_];
		if (b_ > e_) return false;
		b = M_LF.p(b_);
	}

    // Find the lexicographically largest suffix in the current suffix array interval that is prefixed by P[i]
    if (sym != M_LF.L_(e_)) {
		e_=rank_results[sym][e_];
		//e_ = _RS_L_.select(sym,e_);
		e_=select_results[sym][e_-1];
		e = M_LF.p(e_+1)-1;
	}   
	

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
			//cout<<" end b= "<<b<<" e= "<<e<<" b'= "<<b_<<" e'="<<e_<<endl;	

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
			//cout<<" end b= "<<b<<" e= "<<e<<" b'= "<<b_<<" e'="<<e_<<endl;	
        }

    } else {
        M_LF.move(b,b_);
		
        M_LF.move(e,e_);
		
		interval.startPosN=b;	
		interval.endPosN=e;
		interval.startPosNII=b_;	
		interval.endPosNII=e_;	
		//cout<<" end b= "<<b<<" e= "<<e<<" b'= "<<b_<<" e'="<<e_<<endl;	
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

	//cout<<"# dollars before"<<i<<" = "<<a<<endl;
	dataTypeNChar start=_RS_bitvector_.select(1,a-1);

	//dataTypeNChar start = bsel_1(a-1); //bsel_1(k) gives the index of the k-th 1. That is, the index of the first word of the k-th segment.
	
	
	//dataTypeNChar end = bsel_1(a)-1; //same as before. We subtract 1 because we want to find the index of the last word of the (a-1)-th segment.
	dataTypeNChar end =_RS_bitvector_.select(1,a)-1;
		//cout<<"start"<<start<<" end"<<end<<endl;

	#if DEBUG == 1
	//std::cerr << "preceding_dollars_finder - start " << start << " end " << end << "\n";
	#endif

	//std::cerr << "preceding_dollars_finder - start " << start << " end " << end << "\n";

	rangeElement output;
	output.startPosN = start;//provato a togliere il +1 a start e end, va fatto così
	output.startPosNII = M_LF_Dollar_Input_interval[start];
	output.endPosN = end;
	output.endPosNII = M_LF_Dollar_Input_interval[end];
										
	return output;
}

void MOVE_EDSBWT::dollars_in_interval(deque<dataTypeNSeq> &d_out,dataTypeNChar i,dataTypeNChar j){

	dataTypeNSeq l = _RS_L_.rank('#',i);
	//cout<<i<<" "<<l<<" "<<EOF_RANK[i]<<endl;
	dataTypeNSeq u = _RS_L_.rank('#',j+1);
	//cout<<u<<" "<<EOF_RANK[j+1]<<endl;
	//dataTypeNSeq l =EOF_RANK[i];
	//dataTypeNSeq u =EOF_RANK[j+1];
	dataTypeNSeq index;

	//cout<<"ci sono "<<u-l<<" dollari nell'intervallo l="<<i<<", u="<<j<<endl;

	for (dataTypeNSeq k = l; k < u; k++){

		index = EOF_ID_Copy[k];
		//cout<<index<<endl;
		
		if(index > first_symbol_index){	
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