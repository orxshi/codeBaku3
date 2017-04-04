/* 
 * File:   Grid.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 5:41 AM
 */

#ifndef GRID_H
#define	GRID_H

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <functional>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <mpi.h>
#include "../Vector/Vector.h"
#include "../Face/Face.h"
#include "../ADT/ADT.h"
#include "../LinearAlgebra/LinearAlgebra.h"
#include "../Point/Point.h"
#include "../Solid/Solid.h"

using std::ifstream;
using std::ofstream;
using std::vector;
using std::string;
using std::ref;
using std::cout;
using std::endl;
using std::exit;
using std::move;
using std::to_string;
using std::min;
using std::max;
using std::setprecision;
using std::fixed;
using std::function;
using std::flush;
using std::left;
using std::setw;
using std::stringstream;
using std::scientific;

// Decide whether this struct or HoleMap will be used.
struct Hole
{
    array <double, N_DIM> min;
    array <double, N_DIM> max;
};

struct Grid
{
    int tag; // tag of grid.    
    int nHoles; // number of holes in grid.
    double wallDistance; // not clear.    
    string outputDir; // not clear.
    string mainDir; // not clear.
    string logDir; // path of log file.     
    vector <Point> pt; // points.    
    vector <BoundaryFace> bface; // boundary faces.
    vector <InteriorFace> iface; // interior faces.
    vector <Cell> cell; // cells.    
    vector <Hole> holes; // vector of holes.    
    
    // ADT specific for cells. do not know if I can make it generic.
    struct CellADT: ADT
    {
        bool compareFunction (const Node *node, const ADTPoint& targetPoint) override;
        bool doCubesOverlap (const Node *node, const ADTPoint& targetPoint) override;
        void build (const Grid& gr);
    } cellADT;

    // constructor
    Grid (string mainDir, int id);
    //Grid (const Grid&);
    Grid (Grid&& other);
    Grid& operator=(Grid&& other);
    
    // Destructor
    ~Grid ();

    // Private Methods
    void read_ptSize ();
    void read_pt ();
    void read_elmSize ();
    void read_elm ();
    void set_faceVertices (Face& face, const Cell& elm, const int index);
    void cell_to_cell_connectivity ();
    void updateVars();
    double setExpRes();
    void printInput();
    void printMeshInfo();
    void read_grid_GMSH(string meshFile);

    // Public Methods
    void readInput();    
    void read_grid();
    void set_grid();
    void apply_BCs();
    void set_BCs();
    //void createOutputDir(const string mainDir);
    void identifyIBlank(Grid& gr);
    //void leastSquaresCoeffs();
    //void leastSquaresGrad();
    //void interpolate (Gradient& gradient);
    void trimWhoHasFringeNeighbor();
    void trimWhoHasTrimNeighbor (int threshold);
    void trimToUntrim (int crt);
    void nontrimToTrim (int crt);
    void fieldToFringe (int crt);
    void fringeToField (int crt);
    //void blankWithPhys (int phys);
    void getMovingFaceVelocity();
    void rotateMesh (double delAlpha, CVector& centerAirfoil);    
    void outAllTecplot();
    void outAllUnsteady(string title);
    void log();
    void setWallDistance (int phys);
    void outAllVTK (int time);    
};

#endif	/* GRID_H */

