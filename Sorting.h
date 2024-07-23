#ifndef SORTED_INCLUDED
#define SORTED_INCLUDED
		
#include <stdlib.h>
#include <vector>
#include "Parameters.h"

struct IntervalElement{		
	dataTypeNSeq start;
	dataTypeNSeq end;	
};

#if BUILD_LCP == 1
	struct __attribute__((__packed__)) rangeElement {
//	struct rangeElement {	
	rangeElement() {};
	rangeElement( dataTypeNChar x,dataTypeNSeq y, dataTypelenSeq l1, dataTypelenSeq l2) { posN = x; seqN = y; lcpCurN = l1; lcpSucN = l2; };
	~rangeElement() {};
	dataTypelenSeq lcpCurN;
	dataTypelenSeq lcpSucN;
	dataTypeNSeq seqN;	
	dataTypeNChar posN;
	};
	#if RECOVERBW == 1
	struct __attribute__((__packed__)) rangeElementBW {
//	struct rangeElement {	
	rangeElementBW() {};
	rangeElementBW( dataTypedimAlpha z, dataTypeNChar x,dataTypeNSeq y, dataTypelenSeq l1, dataTypelenSeq l2) { pileN = z; posN = x; seqN = y; lcpCurN = l1; lcpSucN = l2; };
	~rangeElementBW() {};
	dataTypelenSeq lcpCurN;
	dataTypelenSeq lcpSucN;
	dataTypedimAlpha pileN;
	dataTypeNSeq seqN;	
	dataTypeNChar posN;
	};
	#endif
#else
	struct __attribute__((__packed__)) rangeElement {
//	struct rangeElement {	
	rangeElement() {};	
	rangeElement( dataTypeNChar x, dataTypeNChar w) { startPosN = x; endPosN = w; };
	~rangeElement() {};
	dataTypeNChar startPosN;
	dataTypeNChar endPosN;
	};
	#if RECOVERBW == 1
	struct __attribute__((__packed__)) rangeElementBW {
//	struct rangeElement {	
	rangeElementBW() {};
	rangeElementBW( dataTypedimAlpha z, dataTypeNChar x, dataTypeNChar w) { pileN = z; startPosN = x; endPosN = w; };
	~rangeElementBW() {};
	dataTypedimAlpha pileN;
	dataTypeNChar startPosN;
	dataTypeNChar endPosN;
	};
	#endif
#endif

#if RECOVERBW == 1
void quickSort(std::vector< rangeElementBW > &v);
#endif

#endif
