#include "Sorting.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

#if RECOVERBW == 1
bool cmpSortEl (rangeElementBW a,rangeElementBW b) { 
	if (a.pileN == b.pileN)
		return (a.startPosN < b.startPosN);
	else
		return (a.pileN<b.pileN); 
}

void quickSort(std::vector< rangeElementBW > &v)
{
  sort( v.begin(),v.end(),cmpSortEl);  
}	
#endif