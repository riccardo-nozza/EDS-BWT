#include <iostream>
#include <fstream>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;


int main (int argc, char *argv[]){
	if(argc != 3 ) {
		std::cerr << "usage: " << argv[0] << " input output" << std::endl;
		std::cerr << "input is the full filename, while output is the filename without any extension.\n";
	}

	char current;
	char next;
	char null_char_input='E';
	char null_char='Z';	//must be lexicographically greater than every other in the input file.
	char end_pos_char='#';
	bit_vector b;
	int empty_counter=0;


	string inputName = argv[1];
	string outputName = argv[2];

	string edsName = inputName;
	ifstream edsFileForLength (edsName);
	cerr<<"Checking number of characters in "<<edsName<<" to initialize the bitvector\n";
	if (edsFileForLength.is_open()){
		edsFileForLength.seekg(0,std::ios_base::end);
		std::ios_base::streampos end_pos = edsFileForLength.tellg();
		b = bit_vector(end_pos,0);
		cerr<<edsName<<" contains "<<end_pos<<" characters\n";
	}
	else{
		cerr << "Errore: impossibile aprire il file "<<edsName<<endl;
		exit(1);
	}
	edsFileForLength.close();


	ifstream edsFile (edsName);

	string newName = outputName + ".fasta";
	ofstream lista (newName);

//	int k=1;	
	long int i=0;	
	long int chars=0;	

	cerr<<"Start reading "<<edsName<<endl;
	
	if (edsFile.is_open() && lista.is_open()){

		edsFile.get(current);
		next=edsFile.peek();
		if(current=='{'){
			b[i]=1;
			i=i+1;
//			lista<<">"<<k<<endl;
			lista<<">1"<<endl;
			if (next == ','){
				//cout<<"Ho trovato la parola vuota, scrivo:"<<null_char<<endl;
				lista<<null_char;
				empty_counter++;
				chars++;
			}
		}
		else {
			cout<<"Error: the string does not start with {."<<endl;
			return 0;
		}


		while (edsFile.get(current)){
			if(current=='{'){
//				lista<<endl<<">"<<k<<endl;
				lista<<endl<<">1"<<endl;
				b[i]=1;
				i=i+1;
				next=edsFile.peek();
				if (next == ','){
					//cerr<<"Ho trovato la parola vuota, scrivo:"<<null_char<<endl;
					lista<<null_char;
					empty_counter++;
					chars++;
				}
				if(i%1000000 == 0){
					cerr<<i/1000000<<" million words read"<<"\r";
				}
			}
			else if(current == '}'){
//				k=k+1;
//				if(max_l < l_tmp){
//					max_l = l_tmp;
//				}
//				l_tmp = 0;
			}
			else if(current == ','){
				b[i]=0;
				i=i+1;
				lista<<endl<<">2"<<endl;
				next=edsFile.peek();
				if (next == ',' or next == '}'){
					lista<<null_char;
					empty_counter++;
					chars++;
				}
				//if(i%1000000 == 0){
				//	cerr<<i/1000000<<" million words read"<<"\r";
				//}
			}
			else if(current == null_char_input){
				//cerr<<"Ho trovato la parola vuota, scrivo:"<<null_char<<endl;
				lista<<null_char;
				empty_counter++;
				chars++;
			}
			else {
//				lista<<toupper(current);
				lista<<current;
//				l_tmp++;
				chars++;
			}
		}
		lista<<endl;
		edsFile.close();
		lista.close();
	}
	else {
		cerr << "Error: open file: " << edsName << " or "<< newName << endl;
		return 0;
	}

cerr<<"\nfine lettura"<<endl;


	b.resize(i);


	string bvName = outputName + ".bitvector";
	store_to_file(b,bvName);

	string empty_name = outputName+".empty.sh";
	ofstream fix_bwt_command_txt (empty_name);
	if (fix_bwt_command_txt.is_open()){
		fix_bwt_command_txt << "#!/bin/bash\n\n";
		fix_bwt_command_txt << "head -c " << (chars + i - empty_counter) << " " << outputName << ".fasta.bwt ";
		fix_bwt_command_txt << "| awk 'BEGIN {ORS=\"\"} {gsub(/" << null_char <<"/,\"" << end_pos_char << "\")}1' ";	//ORS="" make the output on a single line, not generating new unexpected characters in the bwt.
		fix_bwt_command_txt	<< "> " << outputName << ".ebwt";
		fix_bwt_command_txt.close();
	}
	else{
		cerr << "Errore: impossibile aprire il file "<<empty_name<<endl;
	}

	
	return 1;
}