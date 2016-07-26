#ifndef OGA_H
#define	OGA_H

#include "../Grid.h"

// For now only Grid is used instead of MeshBlock.
struct MeshBlock
{
    /* Mesh-block is different than grid.
       Grid refers to complete individual grid while mesh-block could be itself or a part of grid. */
    
    // constructor.
    MeshBlock();       
    // If mesh-block has a hole then 'hm' will be allocated.
    HoleMap* hm;
};

struct HoleMap
{        
    /* Hole map is a map which represents hole region with Cartesian grid. */
    
    // aabb is axis-aligned bounding box (rectangle) which has 4 corners.
    aabb[4];
    // constructor.
    HoleMap();    
    // make AABB of hole (if exists).
    void create_AABB(vector <Cell>& cell);
    // sub-divide hole map so that none of the sub-cells does not intersect both hole and outer boundary.
    void subdivide();
    // fill hole map starting from exterior.
    void fill();
};

// assemble overset mesh-blocks.
void assemble();
// do donor search using ADT.
void searchDonorUsingADT();

namespace OGA
{
    // criterion to evaluate better cell of two overlapping cells.
    // options: wall/volume.
    cellSelectionCriterion;
};

#endif	/* OGA_H */