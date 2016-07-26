#include "Face.h"

Face::Face(int tag, vector<int> vtx, string gshape, pt, const vector<Point>& pt);
{
    vb[0] = 0.;
    vb[1] = 0.;
    vb[2] = 0.;
    
    this.tag = tag;
    this.vtx = vtx;
    
    if (gshape.compare("tri"))
    {
        face_.shape = new GeometricShape::Tri();
    }
    else if (gshape.compare("quad"))
    {
        face_.shape = new GeometricShape::Quad();
    }
    else
    {
        cout << "undefined geometric shape in constructor of f" << endl;
    }
    
    set_area(pt);
    set_centroid(pt);
}

void Face::set_area (const vector<Point>& pt)
{
    // For area of triangle and quad

    unsigned int tri  = static_cast<int>(nVerticesFace_t::TRI);
    unsigned int quad = static_cast<int>(nVerticesFace_t::QUAD);

    if (vtx.size() == tri || vtx.size() == quad)
    {
        area = crossP (pt[vtx[1]].dim - pt[vtx[0]].dim, pt[vtx.back()].dim - pt[vtx[0]].dim);

        if (vtx.size() == tri)
        {
            area *= 0.5;
        }
    }
}

void Face::set_centroid (const vector<Point>& pt)
{
    for (unsigned int i=0; i<cnt.size(); ++i)
    {
        cnt[i] = 0.;
    }

    for (const int p: vtx)
    {        
        cnt += pt[p].dim;
    }

    cnt /= vtx.size();
}