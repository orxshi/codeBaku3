#include "Solver.h"

void Solver::impl (Grid& gr)
{       
    wImpl.start();

    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    
    preSolverCheck (gr);
    
    string dir = gr.outputDir;
    string temps = "res.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);
    int printThres = 1;
    
    if (sOrder == 2) { gradient.leastSquaresGrad (gr); } // parallel
    roe.roeflx (gr, limiter, M0, M1, gradient); // parallel
    
    for (nTimeStep=0; nTimeStep<maxTimeStep; ++nTimeStep)
    {
        wInne.start();        
        interflux(gr); // serial
        
        switch (linearSolverType)
        {
            case 1:
                //gauss_seidel (gr); // only fields
                cout << "this solver type is supported any more" << endl;
                exit(-2);
                break;
            case 2:                
                petsc.solveAxb (gr, M0, M1);                
                break;
            default:
                cout << "undefined linear solver in Solver::impl(...)" << endl;
                exit(-2);
                break;
        }
        
        diff_to_cons_prim (gr); // only fields
        gr.apply_BCs();
        if (sOrder == 2) { gradient.leastSquaresGrad (gr); } // parallel
        roe.roeflx (gr, limiter, M0, M1, gradient); // parallel
        getResiduals (gr); // only fields
        
        if (rank == MASTER_RANK) { outRes(gr.outputDir); }
        
        if (verbose && rank == MASTER_RANK && nTimeStep%printThres==0)
        {
            cout << left << setw(10) << fixed << setprecision(5) << time;
            cout << setw(10) << nTimeStep;
            cout << setw(10) << scientific << setprecision(6) << rmsRes[0] << endl;
        }        
        
        if (fabs(rmsRes[0]) < tol) { break; }
        
        ++glo_nTimeStep;
        
        if (!isSampledInne)
        {
            wInne.stop();
            if (nTimeStep < nSampleInne)
            {
                eTimeInne += wInne.elapsedTime;
            }
            else
            {
                isSampledInne = true;
                printTimeInne();
            }
        }
    }
    
    for (Cell& cll: gr.cell)
    {
        cll.oldold_cons = cll.old_cons;
        cll.old_cons = cll.cons;
    }    
    
    ++nImplicitCalls;
    
    if (!isSampledImpl)
    {
        wImpl.stop();
        if (nImplicitCalls <= nSampleImpl)
        {
            eTimeImpl += wImpl.elapsedTime;
        }
        else
        {
            isSampledImpl = true;
            printTimeImpl();
        }
    }
}
