#include <petscsys.h>
#include "../../../../src/Grid/Grid.h"
#include "../../../../src/IBlank/IBlank.h"
#include "../../../../src/Output/Output.h"
#include "../../../../src/Solid/Solid.h"
#include "../../../../src/Coeffs/Coeffs.h"
#include "../../../../src/Solver/Solver.h"
#include "../../../../src/AFT/AFT.h"
#include "Init/Init.h"
#include "MovingGrid/OscAirfoil.h"
#include "MovingGrid/StraightMovingAirfoil.h"
#include "Output/Output.h"

int main(int argc, char** argv)
{
    PetscInitialize (&argc, &argv, NULL, NULL);
    
    Watch watchSteady;
    Watch watchOscAirfoil;
    
    PetscMPIInt rank, n_procs;
    MPI_Comm world = PETSC_COMM_WORLD;
    MPI_Comm_rank (world, &rank);
    MPI_Comm_size (world, &n_procs);
    
    string mainDir = createOutputDir();
    
    // background grid
    Grid bg (mainDir, 1);
    bg.read_grid();
    bg.set_grid();    
    
    // airfoil grid
    Grid ag (mainDir, 0);
    ag.read_grid();
    ag.set_grid();    
    
    // initialize grids
    OscInit oscInit;
    oscInit.read();
    oscInit.init (ag);
    oscInit.init (bg);
    
    // push grids to vector
    vector<Grid> grs;    
    grs.push_back(move(ag));
    grs.push_back(move(bg));
    
    // set wall distances
    grs[0].setWallDistance(3);
    grs[1].setWallDistance(2);
    
    // build cell trees
    grs[0].cellADT.build (grs[0]);
    grs[1].cellADT.build (grs[1]);        
    
    // hole cutting
    Iblank iblank;
    iblank.identify (grs[0], grs[1]);
    iblank.identify (grs[1], grs[0]);
    iblank.treatFieldIslands (grs[0]);
    iblank.treatFieldIslands (grs[1]);
    iblank.treatFringeIslands (grs[0]);
    iblank.treatFringeIslands (grs[1]);
    iblank.treatVoidAreas (grs[0]);
    iblank.treatVoidAreas (grs[1]);
    
    /*for (int g=0; g<grs.size(); ++g)
    {
        for (int c=grs[g].n_bou_elm; c<grs[g].cell.size(); ++c)
        {
            grs[g].cell[c].iBlank = iBlank_t::FIELD;
        }
    }*/
    
    grs[0].outAllVTK (0);
    grs[1].outAllVTK (0);
    
    exit (-2);
    
    int nActiveElms[2];
    nActiveElms[0] = 0;
    nActiveElms[1] = 0;
    for (int g=0; g<grs.size(); ++g)
    {
        for (const Cell& cll: grs[g].cell)
        {
            if (cll.iBlank == iBlank_t::FIELD)
            {
                ++nActiveElms[g];
            }
        }
    }
    
    // solvers
    array<Solver,2> solver = { Solver(grs[0], "SOLVER-STEADY-AG", grs[0].n_in_elm), Solver(grs[1], "SOLVER-STEADY-BG", grs[1].n_in_elm) };
    solver[0].read ("Solver/solSteady.dat");
    solver[1].read ("Solver/solSteady.dat");
    
    // moving airfoil object
    SMAirfoil sma (solver[0].dt);
    sma.read ("MovingGrid/smAirfoil.dat");
    
    log (mainDir, nActiveElms[0] + nActiveElms[1], "nActiveElms", "");
    cout << "nActiveElms[0] = " << nActiveElms[0] << endl;
    cout << "nActiveElms[1] = " << nActiveElms[1] << endl;
    cout << "nActiveElms = " << nActiveElms[0] + nActiveElms[1] << endl;
    
    watchSteady.start();
    double err = BIG_POS_NUM;
    //solver[0].maxTimeStep = 1;
    //solver[1].maxTimeStep = 1;
    double oversetTol = 1e-12;
    while (err > oversetTol)    
    {
        for (int g=0; g<grs.size(); ++g)
        {
            sma.getAllFaceVelocities (grs[g]);
            
            (solver[g].implicit) ? solver[g].impl(grs[g]) : solver[g].expl(grs[g]);
            
            iblank.interpolate (grs[1-g], solver[g].gradient, grs[g]);
            //grs[1-g].interpolate();
            sma.getAllFaceVelocities (grs[1-g]);
            
            //solver[1-g].getRes (grs[1-g]);
            solver[1-g].getResiduals (grs[1-g]);
            err = solver[1-g].rmsRes[0];
            //err = solver[1-g].setExpRes (grs[1-g]);
        }
        
        cout << "err = " << err << endl;
    }
    watchSteady.stop();
    log (mainDir, watchSteady.elapsedTime, "elapsedTimeSteady", watchSteady.unit);
    
    solver[0].petsc.finalize();
    solver[1].petsc.finalize();
    
    grs[0].outAllVTK (0);
    grs[1].outAllVTK (0);
    
    PetscFinalize();

    return 0;
}