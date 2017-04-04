#include "OGA.h"

HoleCutter::HoleCutter()
{
    holeCuttingMethod = holeCutter_t::EXPLICIT_AABB;
}

void HoleCutter::make_AABB_of_hole(MeshBlock& mb)
{
    // make AABB of hole (if exists).
    
    /*
    Steps:
    - Loop through ghost elements.
    - Pick the ones which has wall boundary condition.
    - Create AABB of hole with respect to minimum and maximum coordinates of elements.
    
    Modifies:
    - HoleCutter::aabb
    */
}

void HoleCutter::explicit_AABB (Grid& g, CVector& min, CVector& max)
{
    /*
    Identify hole cells in grid 'g' by using AABB with user-defined parameters.
    
    min: minimum coordinates (x, y, z) of AABB.
    max: maximum coordinates (x, y, z) of AABB.
    */
    
    // loop through all interior elements. would be nice to work on a reduced set of elements.
    for (int c=g.n_bou_elm; c<g.cell.size(); ++c)
    {
        Cell& cll = g.cell[c];    
    
        cll.iBlank = iBlank_t::HOLE;
        
        for (int j=0; j<N_DIM; ++j)
        {
            if (cll.cnt[j] < min[j] || cll.cnt[j] > max[j])
            {
                cll.iBlank = iBlank_t::UNDEFINED;
                break;
            }            
        }
    }
}