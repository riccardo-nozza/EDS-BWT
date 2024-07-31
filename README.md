# EDS-BWT
The Burrows-Wheeler transform of an elastic-degenerate string

## Install

You need to install both EDS-BWT and gsufsort for the tool to function
```sh
git clone --recursive https://github.com/giovannarosone/EDS-BWT.git 
cd EDS-BWT
make
git clone https://github.com/felipelouza/gsufsort.git
cd gsufsort
make
```

If you want the tool to just search for the number of occurrences of a pattern, without finding the positions, then 
```sh
git clone --recursive https://github.com/giovannarosone/EDS-BWT.git 
cd EDS-BWT
make RECOVERBW=0
git clone https://github.com/felipelouza/gsufsort.git
cd gsufsort
make
```


### Test install
To check if the tool is properly installed, execute
```sh
./eds_to_ebwt.sh test/test/test
./EDS-BWT test/test/test test/test/kmers.txt
```


## Run
EDS-BWT can be used to compute the EDS-BWT of an elastic degenerate string and use it to search patterns in an eds.

To compute the EDS-BWT of an elastic degenerate string written in file input.eds:

```sh
./eds_to_ebwt.sh input
```

Afterwards, to search one or more patterns, contained in file kmers:
```sh
EDS-BWT input kmers
```


If you want to choose the name of the output files that encode the EDS-BWT, execute instead the following istructions:
```sh
./eds_to_fasta input.eds output
gsufsort/gsufsort output.fasta --da --bwt --output output.fasta
./da_to_everything output
```


## Input

### Elastic degenerate string
The input eds file must be a text file with extension .eds, where each degenerate symbol is contained in curly bracket and each string in a degenerate symbol is separated by comma.
Empty strings can be represented both by empty strings themselves, or by using the character for empty string EMPTY_CHAR_EDS (default is E, can be changed in Parameters.h).
Example:

```sh
{ATTGCT}{CTA,TA,A}{CTACGGACT}{A,}{CTGT}
```
and
```sh
{ATTGCT}{CTA,TA,A}{CTACGGACT}{A,E}{CTGT}
```
represent the same eds.

If your input file input.txt does not have curly brackets around all degenerate symbols, such as the following example
```sh
ATTGCT{CTA,TA,A}CTACGGACT{A,}CTGT
```
then please execute
```sh
stringCheck input.txt output
```

This will produce the file output.eds, which will have the brackets as required.


Finally, your input must not contain the character Z, which is internally used to represent empty strings.
If your input does contain this character, please change the parameter EMPTY_CHAR in Parameters.h to a character that is not used in your input and is lexicographically greater than all the characters that you are using. 

### Patterns
The patterns file is a text file where each pattern is in a different line, with no empty lines, for example

```sh
TATT
ACT
TTAT
```

Note that UNIX end-of-line markers (LF) should be used, as opposed for example to Windows ones (CR+LF).


