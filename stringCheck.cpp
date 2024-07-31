#include <iostream>
#include <fstream>
#include "Parameters.h"

using namespace std;


// StringCheck takes in input a text file which represent an elastic degenerate string, but which do not have curly brackets around each degenerate symbol, and adds them.
// It also checks if the eds contains the character, specified in the parameters file, that EDS-BWT uses to internally represent empty words, and warns user if that happen.

int stringCheck (string originalName, string newName){
	char current;
	char next;

	ifstream originalFile (originalName);

	// create file which will contain the output
	string nuovoNome = newName + ".eds";
	cout << "stringCheck on " << originalName <<", output " << nuovoNome << endl;
	ofstream lista (nuovoNome);



	// read original file, adding parenthesys whenever needed
	if (originalFile.is_open() && lista.is_open()){

		//check first character
		originalFile.get(current);
		if(current!=EMPTY_CHAR){
			if(current=='{'){
				lista<<current;
			}
			else {
				lista<<"{"<<current;
			}
		}
		else{
			cerr << "Error: the input contains character " << EMPTY_CHAR << ", which is used by this tool to represent an empty string. Please change EMPTY" << "_CHAR in Parameters.h or change the character in your input." << endl;
			exit(1);
		}

		//check subsequent characters
		while (originalFile.get(current)){
			if(current!=EMPTY_CHAR){
				next=originalFile.peek();
				if(current=='}' && next!='{' && next!=EOF){
					lista<<current<<'{';
				}
				else if(current!='}' && next=='{'){
					lista<<current<<'}';
				}
				else if(current=='<'){
					originalFile.get(current);
					while(current!='>'){
						originalFile.get(current);
					}
					next=originalFile.peek();
					if(next=='{'){
						lista<<'}';
					}
				}
				else if(next==EOF && current!='}'){
					lista<<current<<'}';
				}
				else {
					lista<<current;
				}
			}
			else{
				cerr << "Error: the input contains character " << EMPTY_CHAR << ", which is used by this tool to represent an empty string. Please change EMPTY" << "_CHAR in Parameters.h or change the character in your input." << endl;
				exit(1);
			}
		}

		originalFile.close();
		lista.close();
		cout<<"Done.\n";
	}

	else {
		cerr << "Error: cannot open file " << originalName << " or file " << nuovoNome << endl;
		exit(1);
	}

	return 1;
}









int main (int argc, char *argv[]){
	if(argc != 3 ) {
		std::cerr << "usage: " << argv[0] << " input" << "output" << std::endl;
		std::cerr << "input is the full filename, output will have \".eds\" appended.\n";
		exit(1);
	}

	string inputName = argv[1];
	string outputName = argv[2];
	stringCheck (inputName, outputName);
}

