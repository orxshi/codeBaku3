/* 
 * File:   Face.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 12:30 AM
 */

#ifndef FACE_H
#define	FACE_H

#include <string>
#include <vector>
#include <array>
#include "../Cell/Cell.h"

#define N_NEI_FACE 2

using std::vector;
using std::array;
using std::string;
using std::reference_wrapper;

enum class BC {UNDEFINED=-2, NA=-1, EMPTY=0, SLIP_WALL=1, DIRICHLET=2};

struct Face
{
    GeometricShape::Shape* shape;
    int tag;
    vector <int> vtx;
    vector <int> nei;
    CVector area;
    CVector cnt;
    CVector vb;    
    geometric_shape_t geometric_shape;
    
    // constructor.
    Face(int tag, vector<int> vtx, string gshape, const vector<Point>& pt);
    
    void set_area (const vector<Point>& pt);
    void set_centroid (const vector<Point>& pt);
};

struct InteriorFace:Face
{
};

struct BoundaryFace:Face
{ 
    BC bc;
};

#endif	/* FACE_H */

