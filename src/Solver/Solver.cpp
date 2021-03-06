#include "Solver.h"

Solver::Solver (Grid& gr, string instanceName, int nActiveElms) : petsc(nActiveElms, gr), roe(gr), gradient (gr), limiter (gr)
{
    //default    
    tOrder = 2;
    sOrder = 1;
    linearSolverType = 1; // MYGS
    nGaussIter = 5;
    maxTimeStep = 10000;
    cfl = 5;
    dt = 1.;
    finalTime = 10.;
    tol = 1e-12;    
    steady = false;
    implicit = true;
    verbose = true;    
    nSampleImpl = 5;
    nSampleInne = 5;
    
    // initialize
    time = 0.;
    glo_nTimeStep = 0;
    nImplicitCalls = 0;    
    eTimeImpl = 0.;
    eTimeInne = 0.;
    isSampledImpl = false;
    isSampledInne = false;
    this->instanceName = instanceName;
    M0.resize (gr.face.size());
    M1.resize (gr.face.size());
}

Solver::Petsc::Petsc (int nActiveElms, Grid& gr)
{    
    //cout << "inside petsc constructor" << endl;

    int nProcs;

    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    /*n = 0;
    for (Cell& cll: gr.cell)
    {
        if (cll.iBlank == iBlank_t::FIELD)
        {
            ++n;
        }
    }*/
    
    //n = gr.n_in_elm;
    n = nActiveElms;
    bs = N_VAR;
    vecGlobalSize = n*bs;
    //DX = (double*) malloc (vecGlobalSize*sizeof(double));
    //PetscMalloc1 (xGlobalSize, &DX);
    DX = new double [vecGlobalSize];
    
    
    
    // set x
    VecCreate (world, &x);
    VecSetType (x, VECSTANDARD);
    VecSetBlockSize(x, bs);
    VecSetSizes (x, PETSC_DECIDE, vecGlobalSize);
    VecGetLocalSize (x, &vecLocalSize);
    VecGetOwnershipRange (x, &vecLocBeg, &vecLocEnd);
    
    //cout << "set x" << endl;

    // set b
    VecDuplicate (x, &b);
    
    //cout << "set b" << endl;
    
    // set A
    MatCreate (world, &A);    
    MatSetType (A, MATMPIBAIJ);
    //MatSetType (A, MATMPIAIJ);
    MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE, vecGlobalSize, vecGlobalSize);
    //MatSeqBAIJSetPreallocation (A, bs, 4, NULL);
    //MatSeqAIJSetPreallocation (A, 4*bs, NULL);    
    //MatMPIAIJSetPreallocation (A, 4*bs, NULL, 4*bs, NULL);
    MatMPIBAIJSetPreallocation (A, bs, 5, NULL, 5, NULL);
    //MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    //MatCreateBAIJ (world, bs, PETSC_DECIDE, PETSC_DECIDE, n*bs, n*bs, 1, NULL, 3, NULL, &A);
    // specific to pentagonal mesh . change later
    MatGetOwnershipRange (A, &matLocBeg, &matLocEnd);
    MatGetLocalSize (A, &matLocalSize, NULL);
    
    //cout << "set matrix" << endl;
    
    KSPCreate (world, &ksp);
    KSPSetOperators (ksp, A, A);
    KSPGetPC (ksp, &pc);
    PCSetType (pc, PCSOR);
    KSPSetType (ksp, KSPGMRES);
    KSPSetFromOptions (ksp);
        
    localSizes = new int [nProcs];
    int recvcounts[nProcs];
    int displs[nProcs];
    displs[0] = 0;    
    for (int i=0; i<nProcs; ++i)
    {
        recvcounts[i] = 1;
    }
    for (int i=1; i<nProcs; ++i)
    {
        displs[i] = displs[i-1] + recvcounts[i-1];
    }
    
    MPI_Allgatherv (&vecLocalSize, 1, MPI_INT, localSizes, recvcounts, displs, MPI_INT, world);
    
    //ic.reserve(n);
    //icn.reserve(n);
    
    //IC (gr);
}

void Solver::preSolverCheck (const Grid& gr)
{
    for (int c=0; c<gr.cell.size(); ++c)
    //for (const Cell& cll: gr.cell)
    {
        const Cell& cll = gr.cell[c];
    
        if (cll.iBlank == iBlank_t::UNDEFINED)
        {
            cout << "iBlank is undefined in Solver::preSolverCheck(...)" << endl;
            cout << "grid id = " << gr.id << endl;
            cout << "c = " << c << endl;
            exit(-2);
        }
    }
}

bool Solver::cm (string s, ifstream& in)
{    
    bool found = false;
    
    cmh (s, STRINGTIFY(nGaussIter), nGaussIter, in, found);
    cmh (s, STRINGTIFY(maxTimeStep), maxTimeStep, in, found);
    cmh (s, STRINGTIFY(cfl), cfl, in, found);
    cmh (s, STRINGTIFY(dt), dt, in, found);
    cmh (s, STRINGTIFY(finalTime), finalTime, in, found);
    cmh (s, STRINGTIFY(tol), tol, in, found);    
    cmh (s, STRINGTIFY(steady), steady, in, found);
    cmh (s, STRINGTIFY(implicit), implicit, in, found);
    cmh (s, STRINGTIFY(verbose), verbose, in, found);
    cmh (s, STRINGTIFY(tOrder), tOrder, in, found);
    cmh (s, STRINGTIFY(sOrder), sOrder, in, found);
    cmh (s, STRINGTIFY(linearSolverType), linearSolverType, in, found);
    cmh (s, STRINGTIFY(nSampleImpl), nSampleImpl, in, found);
    cmh (s, STRINGTIFY(nSampleInne), nSampleInne, in, found);
    
    return found;
}

void Solver::read (string fileName)
{
    string tmps;
    ifstream in;
    in.open (fileName);
    
    if (in.is_open())
    {
        in >> tmps;
        while ( !in.eof() )
        {
            if ( cm (tmps, in) == false )
            {
                cout << "undefined input in Solver::read(...)" << endl;
                exit(-2);
            }

            in >> tmps;
        }
    }
    else
    {
        cout << "could not open file in Solver::read(...)" << endl;
        exit(-2);
    }    

    in.close();
}

void Solver::diff_to_cons_prim(Grid& g)
{
    //double maxDq = BIG_NEG_NUM;

    for (int ic=g.n_bou_elm; ic<g.cell.size(); ++ic)
    {
        Cell& e = g.cell[ic];
        
        if (e.iBlank == iBlank_t::FIELD)
        {
            e.cons += e.dQ;
            e.cons_to_prim();
            
            //maxDq = max (maxDq,e.dQ[1]);
            /*if (ic == g.n_bou_elm)
            {
                cout << "e.dQ[3] = " << e.dQ[3] << endl;
            }*/
        }
    }
    
    //cout << "maxDq = " << maxDq << endl;
}