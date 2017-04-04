#ifndef HOLECUTTER_H
#define	HOLECUTTER_H

#include "../Grid.h"

struct HoleCutter
{
    enum holeCutter_t {EXPLICIT_AABB=0, CARTESIAN_HOLE_MAP=1, DIRECT=2}; // define hole cutters.
    holeCutter_t holeCuttingMethod; // create a hole cutter.
    
    HoleCutter();
    
    void explicit_AABB (Grid& g, CVector& min, CVector& max);
};

#endif	/* HOLECUTTER_H */