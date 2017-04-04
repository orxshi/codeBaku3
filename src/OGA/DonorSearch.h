#ifndef DONORSEARCH_H
#define	DONORSEARCH_H

#include "../Grid.h"

struct DonorSearch
{
    enum cellSelection_t {WALL=0, SIZE=1};    
    cellSelectinon_t cellSelectionCriterion; 

    void searchWithADT();
};

#endif	/* DONORSEARCH_H */