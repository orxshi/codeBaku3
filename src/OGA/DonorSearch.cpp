#include "OGA.h"

using namespace OGA;

void searchDonorUsingADT(MeshBlock& mbA, queryPoints)
{
    // find donor cells of mesh-block, 'mbA' which contain query points, 'queryPoints'.
    
    // make sure that input mesh-blocks are valid.

    // lambda function.
    function<int(Cell&)> getIndex = [&] (Cell& cll)
    {
        ADT::ADTPoint vec;
        CVector cnt;
        
        cnt = cll.cnt;
        
        for (int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = cnt[i];
            vec.dim[i*2+1] = cnt[i];
        }             
        
        return grPas.cellADT.search (vec);
    };
    
    // loop through all interior elements. would be nice to work on a reduced set of elements.
    for (int c=grAct.n_bou_elm; c<grAct.cell.size(); ++c)
    {
        int index;        
        Cell& cll = grAct.cell[c];

        index = getIndex (cll);
                
        // if an overlap found.
        if (index != -1)
        {
            double var1;
            double var2;
            
            if (cellSelectionCriterion == cellSelection_t::WALL)
            {
                var1 = grPas.cell[index].wallDistance;
                var2 = cll.wallDistance;
            }
            else if (cellSelectionCriterion == cellSelection_t::SIZE)
            {
                var1 = grPas.cell[index].vol;
                var2 = cll.vol;
            }
            
            if (var1 < var2)
            {
                cll.iBlank = iBlank_t::FRINGE;
                cll.donor = &grPas.cell[index];
            }
            else
            {
                cll.iBlank = iBlank_t::FIELD;
            }
        }
        else
        {
            // oprhan
        }
    }
}