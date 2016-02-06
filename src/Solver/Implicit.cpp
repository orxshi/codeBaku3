#include "Solver.h"

void Solver::impl (Grid& gr)
{       
    wImpl.start();

    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    
    //preSolverCheck (gr);
    
    string dir = gr.outputDir;
    string temps = "res.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);
    int printThres = 1;
    
    /*if (sOrder == 2)
    {
        gradient.leastSquaresGrad (gr); // parallel
    }
    
    roe.roeflx (gr, limiter, M0, M1, gradient); // parallel*/
    
    
    for (nTimeStep=0; nTimeStep<maxTimeStep; ++nTimeStep)
    {
        wInne.start();
        
        if (sOrder == 2)
        {
            gradient.leastSquaresGrad (gr); // parallel
        }
        
        roe.roeflx (gr, limiter, M0, M1, gradient); // parallel
        //cout << "R[0] = " << gr.cell[gr.n_bou_elm].R[0] << endl;
        /*cout << "R[0] = " << gradient.grad[0][0][0] << endl;
        cout << "R[0] = " << gradient.grad[0][1][0] << endl;
        cout << "R[0] = " << gradient.grad[0][2][0] << endl;
        cout << "R[0] = " << gradient.grad[0][3][0] << endl;
        cout << "R[0] = " << gradient.grad[0][4][0] << endl;
        cin.ignore();*/
        
        interflux(gr); // serial
        
        /*double sumr = 0.;
        for (int c=0; c<gr.cell.size(); ++c)
        {
            sumr += gr.cell[c].R[0];
        }
        cout << sumr << endl;
        exit(-2);*/
        
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
        //getMaxRes (gr); // only fields // includes roe
        getRmsRes (gr); // only fields // includes roe
        
        if (rank == MASTER_RANK) { outRes(gr.outputDir); }
        
        
        if (verbose && rank == MASTER_RANK && nTimeStep%printThres==0)
        {
            cout << left << setw(10) << fixed << setprecision(5) << time;
            cout << setw(10) << nTimeStep;
            cout << setw(10) << scientific << setprecision(6) << (rmsRes[0] + rmsRes[1] + rmsRes[2] + rmsRes[3] + rmsRes[4])/5. << endl;
            //cout << setw(10) << scientific << setprecision(3) << maxRes[0] << endl;
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
