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
    PetscMPIInt rank, n_procs;
    MPI_Comm world = PETSC_COMM_WORLD;
    MPI_Comm_rank (world, &rank);
    MPI_Comm_size (world, &n_procs);
    
    Watch watchSteady;
    Watch watchOscAirfoil;
    Watch watchAFT;
    Watch watchPre;
    Watch watchIblank;    
    
    //watchPre.start();
    double prestart = MPI_Wtime();
    
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
    grs.push_back(move(bg));
    grs.push_back(move(ag));
    
    // set wall distances
    grs[0].setWallDistance(3);
    grs[1].setWallDistance(2);
    
    grs[0].cellADT.build (grs[0]);
    grs[1].cellADT.build (grs[1]);
    
    //watchPre.stop();
    double preend = MPI_Wtime();
    cout << "pre = " << preend - prestart << endl;
    //log (mainDir, watchPre.elapsedTime, "elapsedTimePre", watchPre.unit);
    
    //watchIblank.start();
    double iblankstart = MPI_Wtime();
    
    /*Iblank iBlank;
    iBlank.identify (grs[0], grs[1]);
    iBlank.identify (grs[1], grs[0]);*/
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
    
    //watchIblank.stop();
    //log (mainDir, watchIblank.elapsedTime, "elapsedTimeIblank", watchIblank.unit);
    double iblankend = MPI_Wtime();
    cout << "iblank = " << iblankend - iblankstart << endl;

    for (int g=0; g<grs.size(); ++g)
    {
        for (int c=0; c<grs[g].cell.size(); ++c)
        {
            if (grs[g].cell[c].iBlank == iBlank_t::UNDEFINED)
            {
                cout << "undefined iblank" << endl;
                cout << "g = " << g << endl;
                cout << "c = " << c << endl;
                cout << "d = " << grs[g].n_bou_elm << endl;
                cout << "d = " << grs[g].cell.size() << endl;
                exit(-2);
            }
        }
    }
    
    grs[0].outAllVTK (0);
    grs[1].outAllVTK (0);
    
    
    
    Grid finalGrid (mainDir, 3);
    
    
    //watchAFT.start();
    double aftstart = MPI_Wtime();
    AFT::aft (grs, finalGrid);    
    //cout << "out of AFT" << endl;
    //watchAFT.stop();    
    double aftend = MPI_Wtime();
    cout << "aft = " << aftend - aftstart << endl;
    //cout << "stopped aft watch" << endl;
    //log (mainDir, watchAFT.elapsedTime, "elapsedTimeAFT", watchAFT.unit);
    //cout << "logged AFT" << endl;
    
    finalGrid.outAllVTK (0);
    //cout << "output final grid" << endl;
    
    exit(-2);
    
    finalGrid.readInput();
    //cout << "read final grid" << endl;
    //finalGrid.leastSquaresCoeffs();    
    finalGrid.cellADT.build (finalGrid);
    //cout << "built final grid" << endl;
    oscInit.init (finalGrid);
    //cout << "osc init" << endl;
    
    Solver solSteady (finalGrid, "SOLVER-STEADY", finalGrid.n_in_elm);
    //cout << "made solSteady" << endl;
    solSteady.read ("Solver/solSteady.dat");
    //cout << "read solSteady" << endl;

    // solve steady state
    SMAirfoil sma (solSteady.dt);
    OscAirfoil oa (1.); // 1 is time step
    sma.read ("MovingGrid/smAirfoil.dat");
    oa.read ("MovingGrid/oscAirfoil.dat");
    
    //cout << "ma read" << endl;

    //Coeffs coeffs (finalGrid, oscInit.rhoInf, oscInit.pInf, oscInit.Mach, oa.MachAirfoil);
    
    sma.getAllFaceVelocities (finalGrid);
    //cout << "sma read" << endl;
    watchSteady.start();
    (solSteady.implicit) ? solSteady.impl(finalGrid) : solSteady.expl(finalGrid);
    solSteady.petsc.finalize();
    watchSteady.stop();
    
    //finalGrid.outAllVTK (0);
    exit(-2);
        
    // solve osc airfoil
    //Grid oldGrid (mainDir, 4);
    //Grid oldGrid = move(finalGrid);
    int countr = 0;
    watchOscAirfoil.start();
    for (double time=0.; time<50.; time+=1.) // 1 is dt
    {
        cout << "time = " << time << endl;
        
        grs[0].cellADT.build (grs[0]);
        grs[1].cellADT.build (grs[1]);   
        grs[0].identifyIBlank (grs[1]);
        grs[1].identifyIBlank (grs[0]);
        
        grs[0].outAllVTK (countr);
        grs[1].outAllVTK (countr);
        
        Grid finalGrid (mainDir, 3);
        AFT::aft (grs, finalGrid);        
        finalGrid.cellADT.build (finalGrid);        
        finalGrid.readInput();
        //finalGrid.leastSquaresCoeffs();
        
        if (time == 0.)
        {            
            /*oa.delAlpha = 0.;
            oscInit.init (finalGrid);
            oa.interFromOldTS (finalGrid, oldGrid);*/
            //finalGrid = move(oldGrid);
        }
        else
        {   
            /*oscInit.init (finalGrid);
            oa.interFromOldTS (finalGrid, oldGrid);*/
            //finalGrid.set_BCs();
            //finalGrid.apply_BCs();
        }
        
        Solver solOscAirfoil (finalGrid, "SOLVER-OSC-AIRFOIL", finalGrid.n_in_elm);
        solOscAirfoil.read ("Solver/solOscAirfoil.dat");
        solOscAirfoil.time = time;
        
        oa.setAngles (time);
        /*oa.getAllFaceVelocities (finalGrid);
        (solOscAirfoil.implicit) ? solOscAirfoil.impl(finalGrid) : solOscAirfoil.expl(finalGrid);
        coeffs.getCoeffs (finalGrid);
        outLiftCoef (coeffs, oa.alpha, solOscAirfoil.time);
        coeffs.outPresCoef (countr);*/
        finalGrid.outAllVTK (countr);
        oa.moveGrid (grs[1]);
        
        //oldGrid = move(finalGrid);
        
        ++countr;
    }
    watchOscAirfoil.stop();
    
    /*if (rank == MASTER_RANK)
    {
        //gr.outAllTecplot();
        finalGrid.outAllVTK (0);
        //coeffs.out.close();
        log (mainDir, watchSteady.elapsedTime, "elapsedTimeSteady", watchSteady.unit);
        //log (mainDir, watchOscAirfoil.elapsedTime, "elapsedTimeOscAirfoil", watchOscAirfoil.unit);
        solSteady.log (finalGrid.logDir);
        //solOscAirfoil.log (gr.logDir);
        sma.log (finalGrid.logDir);
        //oa.log (gr.logDir);
    }*/
    
    PetscFinalize();

    return 0;
}