/* 
 * File:   Cell.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 3:40 AM
 */

#ifndef CELL_H
#define	CELL_H

#include <vector>
#include <memory>
#include <array>
#include "../Vector/Vector.h"
#include "../Constants.h"
#include "../Matrix5/Matrix5.h"
#include "../Point/Point.h"

using std::vector;
using std::array;
using std::reference_wrapper;
using std::move;

class Face; // forward declaration

enum class iBlank_t {UNDEFINED, NA, HOLE, FRINGE, FIELD};

enum class vtkCellType_t {UNDEFINED=-1, TET=10, HEX=12, WEDGE=13, TRI=5};
enum class fringeBou_t {UNDEFINED=-1, NO=0, YES=1};





struct Cell
{
    GeometricShape::Shape shape*;
    int phys; // transform this to something else otherwise not clear.
    Cell* donor;
    vector<Cell*> receiver;
    iBlank_t iBlank;
    int belonging;
    int nTrims;    
    
    vector <int> nei; // only cell neighbors.
    
    double Mach;
    double sigma;
    double wallDistance;
    double vol;
    Vector<N_VAR> R;
    Vector<N_VAR> dQ, old_dQ;
    Vector<N_VAR> prim, cons, old_cons, oldold_cons;
    Vector<N_VAR> resInner, resOuter;
    Matrixd<N_VAR, N_VAR> D;
    bool trim; 
    bool trimmedNow;
    bool newlyCreated;
    bool ghost;
    Vector<3> cnt;
    //Vector2D <N_DIM,N_VAR> grad;
    Vector2D <N_DIM,N_VAR> emin;
    Vector2D <N_DIM,N_VAR> emax;
    vector <int> face; // stores indices of faces.
    geometric_shape_t geometric_shape;
    vector <int> vtx;
    vector <int> vtxBelo;
    BC bc;    
    fringeBou_t fringeBou;

    // Constructors
    Cell();
    Cell (const Cell& other);
    Cell& operator= (const Cell& other);
    Cell (Cell&& other);

    // Methods
    void set_centroid (const vector<Point>& pt);
    void prim_to_cons();
    void cons_to_prim();
    //void interpolate();
    
};

#endif	/* ELEMENT_H */

