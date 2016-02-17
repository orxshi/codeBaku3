#include "Solver.h"

void Solver::Petsc::solveAxb (Grid& gr, vector <Matrixd<N_VAR,N_VAR>>& M0, vector <Matrixd<N_VAR,N_VAR>>& M1)
{
    //cout << "inside petsc" << endl;

    int rank, nProcs;

    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);

    // set values of b
    int ind[vecLocalSize];
    double val[vecLocalSize];
    
    for (PetscInt gp=vecLocBeg; gp<vecLocEnd; ++gp)
    {
        PetscInt brow = static_cast <int> (floor(gp/bs));
        PetscInt c = brow + gr.n_bou_elm;
        //PetscInt c = ic[brow];
        PetscInt i = gp % bs;
        
        Cell& cll = gr.cell[c];        
        
        ind[gp-vecLocBeg] = gp;
        
        /*if (cll.iBlank == iBlank_t::FIELD)
        {
            ind[gp-vecLocBeg] = gp;
        }
        else if (cll.iBlank == iBlank_t::FRINGE)
        {
            ind[gp-vecLocBeg] = -gp;
        }*/
        
        val[gp-vecLocBeg] = cll.R[i];
    }
    
    //cout << "loaded RHS vector in petsc" << endl;
        
    VecSetValues (b, vecLocalSize, ind, val, INSERT_VALUES);
    
    VecAssemblyBegin (b);
    VecAssemblyEnd (b);
    
    // set values of A
    for (PetscInt brow=matLocBeg/bs; brow<matLocEnd/bs; ++brow)
    {
        PetscInt c = brow + gr.n_bou_elm;
        //PetscInt c = ic[brow];
        
        Cell& cll = gr.cell[c];
        
        
        
        /*PetscInt idxm;
        PetscInt idxn;
        
        if (cll.iBlank == iBlank_t::FIELD)
        {
            idxm = brow;
            idxn = brow;
        }
        else if (cll.iBlank == iBlank_t::FRINGE)
        {
            idxm = -brow;
            idxn = -brow;
        }*/
        
        
        
        PetscInt idxm = brow;
        PetscInt idxn = brow;
        double v[bs*bs];
        
        for (int q=0; q<bs; ++q)
        {
            for (int w=0; w<bs; ++w)
            {
                v[q*bs+w] = cll.D(q,w);
            }
        }
        
        MatSetValuesBlocked (A, 1, &idxm, 1, &idxn, v, INSERT_VALUES);                
        
        for (int nn=0; nn<cll.nei.size(); ++nn)
        {
            if (cll.nei[nn] >= gr.n_bou_elm)
            {
                //if (gr.cell[cll.nei[nn]].iBlank == iBlank_t::FIELD)
                if (true)
                {
                    /*if (gr.cell[cll.nei[nn]].iBlank == iBlank_t::FIELD)
                    {
                        idxn = brow;
                    }
                    else if (gr.cell[cll.nei[nn]].iBlank == iBlank_t::FRINGE)
                    {
                        idxn = -brow;
                    }*/
                    
                    
                
                    idxn = cll.nei[nn] - gr.n_bou_elm;
                    //idxn = icn[brow][nn] - gr.n_bou_elm;
                    //idxn = nn - gr.n_bou_elm;
                    
                    Face& f = gr.face[cll.face[nn]];
                    if (c == f.nei[0])
                    {
                        for (int q=0; q<bs; ++q)
                        {
                            for (int w=0; w<bs; ++w)
                            {
                                v[q*bs+w] = M1[cll.face[nn]](q,w);
                            }
                        }
                    }
                    else
                    {
                        for (int q=0; q<bs; ++q)
                        {
                            for (int w=0; w<bs; ++w)
                            {
                                v[q*bs+w] = -M0[cll.face[nn]](q,w);
                            }
                        }
                    }
                    
                    MatSetValuesBlocked (A, 1, &idxm, 1, &idxn, v, INSERT_VALUES);
                }
            }
        }
    }
    
    MatAssemblyBegin (A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY);    
    
    //cout << "loaded matrix in petsc" << endl;
    
    KSPSolve (ksp, b, x);
    
    /*double* bb = NULL;
    bb = new double [5];
    VecGetArray (b, &bb);
    cout << bb[0] << endl;
    cout << bb[1] << endl;
    cout << bb[2] << endl;
    cout << bb[3] << endl;
    cout << bb[4] << endl;
    VecRestoreArray (b, &bb);*/
    
    int recvcounts[nProcs];
    int displs[nProcs];
    displs[0] = 0;
    recvcounts[0] = localSizes[0];
    
    for (int i=1; i<nProcs; ++i)
    {
        recvcounts[i] = localSizes[i];
        displs[i] = displs[i-1] + recvcounts[i-1];
    }
    
    double *dx = NULL;
    VecGetArray (x, &dx);
    MPI_Allgatherv (dx, localSizes[rank], MPI_DOUBLE, DX, localSizes, displs, MPI_DOUBLE, world);
    VecRestoreArray (x, &dx);
    //
    
    for (PetscInt gp=0; gp<vecGlobalSize; ++gp)
    {
        PetscInt brow = static_cast <int> (floor(gp/bs));
        PetscInt c = brow + gr.n_bou_elm;
        //PetscInt c = ic[brow];
        PetscInt i = gp % bs;        
        Cell& cll = gr.cell[c];
        
        cll.dQ[i] = DX[gp];
        
        
        /*if (c == gr.n_bou_elm && i == 3)
        {
            cout << "dx[gp] = " << dx[gp] << endl;
        }*/
    }
    
    
}

void Solver::Petsc::finalize()
{
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    //PetscFree (DX);
    delete DX;
    //free (DX);
    DX = NULL;
    delete [] localSizes;
    ic.clear();
}

void Solver::Petsc::IC (Grid& gr)
{
    for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
    {
        if (gr.cell[c].iBlank == iBlank_t::FIELD)
        {
            ic.push_back (c);
            
            for (int nn=0; nn<gr.cell[c].nei.size(); ++nn)
            {
                if (gr.cell[gr.cell[c].nei[nn]].iBlank == iBlank_t::FIELD && gr.cell[c].nei[nn]>=gr.n_bou_elm)
                {
                    icn[ic.size()-1].push_back (gr.cell[c].nei[nn]);
                }
            }
        }
    }
}