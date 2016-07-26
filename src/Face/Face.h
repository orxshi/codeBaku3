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

enum class nVerticesFace_t {UNDEFINED=-1, TRI=3, QUAD=4};

struct Face
{
    int tag;
    int phys; // transform this to something else otherwise not clear.
    int nVertices;
    int nFaces;
    vector <int> vtx;
    vector <int> nei;
    CVector area;
    CVector cnt;
    CVector vb;    
    geometric_shape_t geometric_shape;
    
    // constructor.
    Face();
    
    void set_area (const vector<Point>& pt);
    void set_centroid (const vector<Point>& pt);
};

struct InteriorFace:Face
{
};

struct BoundaryFace:Face
{
};

#endif	/* FACE_H */

