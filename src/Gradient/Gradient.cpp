#include "Gradient.h"

Gradient::Gradient (Grid& gr)
{
    grad.resize (gr.n_in_elm);
    
    for (int ic=0; ic<gr.n_in_elm; ++ic)
    {
        for (int i=0; i<N_VAR; ++i)
        {
            grad[ic][i].fill(0.);
        }
    }
    
    Wx.resize (gr.n_in_elm);
    Wy.resize (gr.n_in_elm);
    Wz.resize (gr.n_in_elm);

    leastSquaresCoeffs (gr);
    initParallelVars (gr);
}

void Gradient::initParallelVars (Grid& gr)
{
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    localSizes  = new int [nProcs];
    localSizesNVARNDIM  = new int [nProcs];
    displs = new int [nProcs];
    displsNVARNDIM = new int [nProcs];

    int rem;

    // get localSizeFace
    rem = gr.n_in_elm % nProcs;
    
    localSize = (gr.n_in_elm - rem) / nProcs;
    
    if (rank == nProcs-1)
    {
        localSize += rem;
    }
    
    localSizeNVARNDIM = localSize * N_VAR * N_DIM;
    
    // gather localSizesFace and localSizesFaceM, F, V
    MPI_Allgather (&localSize, 1, MPI_INT, localSizes, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather (&localSizeNVARNDIM, 1, MPI_INT, localSizesNVARNDIM, 1, MPI_INT, MPI_COMM_WORLD);
    
    displs[0] = gr.n_bou_elm;
    for (int i=1; i<nProcs; ++i)
    {
        displs[i] = displs[i-1] + localSizes[i-1];
    }
    
    displsNVARNDIM[0] = 0;
    for (int i=1; i<nProcs; ++i)
    {
        displsNVARNDIM[i] = displsNVARNDIM[i-1] + localSizesNVARNDIM[i-1];
    }
}

void Gradient::leastSquaresCoeffs (Grid& gr)
{    
    double dx, dy, dz, a1, a2, a3, psi;
    double r_11;
    double r_12;
    double r_13;
    double r_22;
    double r_23;
    double r_33;
    CVector d;

    for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
    {
        Cell& cll = gr.cell[c];
        
        r_11 = 0.;
        r_12 = 0.;
        r_13 = 0.;
        r_22 = 0.;
        r_23 = 0.;
        r_33 = 0.;
        
        for (const int f: gr.cell[c].nei)        
        {            
            d = gr.cell[f].cnt - cll.cnt;

            dx = d[0];
            dy = d[1];
            dz = d[2];

            r_11 += dx * dx;
            r_12 += dx * dy;
            r_13 += dx * dz;
            r_22 += dy * dy;
            r_23 += dy * dz;
            r_33 += dz * dz;
        }
        
        
        
        //cout << "dx = " << dx << endl;
        //cout << "dy = " << dy << endl;
        //cout << "dz = " << dz << endl;
        
        //cout << "r_12 = " << r_12 << endl;
        //cout << "r_22 = " << r_22 << endl;

        r_11 = sqrt (r_11);
        r_12 /= r_11;
        r_13 /= r_11;
        r_22 = sqrt (r_22 - pow(r_12,2));
        r_23 = (r_23 - r_12 * r_13) / r_22;
        r_33 = sqrt (r_33 - pow(r_13,2.) - pow(r_23,2.));
        
        
            
        
        for (const int f: cll.nei)
        {
            d = gr.cell[f].cnt - cll.cnt;

            dx = d[0];
            dy = d[1];
            dz = d[2];

            a1 = dx / pow (r_11,2.);
            a2 = (dy - dx * r_12 / r_11) / pow (r_22,2.);
            psi = (r_12 * r_23 - r_13 * r_22) / (r_11 * r_22);
            a3 = (dz - dy * r_23 / r_22 + psi * dx) / pow (r_33,2.);            
            
            Wx[c-gr.n_bou_elm].push_back (a1 - a2 * r_12 / r_11 + psi * a3);
            Wy[c-gr.n_bou_elm].push_back (a2 - a3 * r_23 / r_22);
            Wz[c-gr.n_bou_elm].push_back (a3);
            
            //cout << "a3 = " << a3 << endl;
             
        }
    }    
}

void Gradient::leastSquaresGrad (Grid& gr)
{    
    double tempf;
    //CVector d;

    
    double dx, dy, dz, a1, a2, a3, psi;
    double r_11;
    double r_12;
    double r_13;
    double r_22;
    double r_23;
    double r_33;
    CVector d;
    double Wxx, Wyy, Wzz;
    
    
    
    
    
    for (int ic=displs[rank]; ic<displs[rank]+localSize; ++ic)
    {
        Cell& cll = gr.cell[ic];
        
        for (int i=0; i<N_VAR; ++i)
        {
            grad[ic-gr.n_bou_elm][i].fill(0.);
        }
        
        /*r_11 = 0.;
            r_12 = 0.;
            r_13 = 0.;
            r_22 = 0.;
            r_23 = 0.;
            r_33 = 0.;
            
        
        for (const int f: cll.nei)        
        {            
            d = gr.cell[f].cnt - cll.cnt;

            dx = d[0];
            dy = d[1];
            dz = d[2];

            r_11 += dx * dx;
            r_12 += dx * dy;
            r_13 += dx * dz;
            r_22 += dy * dy;
            r_23 += dy * dz;
            r_33 += dz * dz;
        }
        
      
        
        r_11 = sqrt (r_11);
        r_12 /= r_11;
        r_13 /= r_11;
        r_22 = sqrt (r_22 - pow(r_12,2));
        r_23 = (r_23 - r_12 * r_13) / r_22;
        r_33 = sqrt (r_33 - pow(r_13,2.) - pow(r_23,2.));
        
        /*if (isnan(r_33))
        {
            cout << "r11 = " << r_11 << endl;
            cout << "r12 = " << r_12 << endl;
            cout << "r13 = " << r_13 << endl;
            cout << "r22 = " << r_22 << endl;
            cout << "r23 = " << r_23 << endl;
            cout << "r33 = " << r_33 << endl;
            exit(-2);
        }*/
        
        
        int ff=0;
        for (const int f: cll.nei)
        {
            for (int n=0; n<N_VAR; ++n)
            {
                tempf = gr.cell[f].prim[n] - cll.prim[n];
                
                
                
                
                /*d = gr.cell[f].cnt - cll.cnt;

            dx = d[0];
            dy = d[1];
            dz = d[2];

            a1 = dx / pow (r_11,2.);
            a2 = (dy - dx * r_12 / r_11) / pow (r_22,2.);
            psi = (r_12 * r_23 - r_13 * r_22) / (r_11 * r_22);
            a3 = (dz - dy * r_23 / r_22 + psi * dx) / pow (r_33,2.);            
            
            Wxx = (a1 - a2 * r_12 / r_11 + psi * a3);
            Wyy =  (a2 - a3 * r_23 / r_22);
            Wzz = (a3);
            
            grad[ic-gr.n_bou_elm][n][0] += Wxx * tempf;
                grad[ic-gr.n_bou_elm][n][1] += Wyy * tempf;
                grad[ic-gr.n_bou_elm][n][2] += Wzz * tempf;*/
            
            
            
                
                
                
                
                

                grad[ic-gr.n_bou_elm][n][0] += Wx[ic-gr.n_bou_elm][ff] * tempf;
                grad[ic-gr.n_bou_elm][n][1] += Wy[ic-gr.n_bou_elm][ff] * tempf;
                grad[ic-gr.n_bou_elm][n][2] += Wz[ic-gr.n_bou_elm][ff] * tempf;
                
                /*if (ic-gr.n_bou_elm == 0 && n == 3)
                {
                    cout << "nei = " << f << endl;
                    cout << "Wx = " << Wx[ic-gr.n_bou_elm][ff] << endl;
                    cout << "tempf = " << tempf << endl;
                    cout << "grad = " << grad[ic-gr.n_bou_elm][n][0] << endl;
                    cout << "gr.cell[f].prim[n] = " << gr.cell[f].prim[n] << endl;
                    cout << "cll.prim[n] = " << cll.prim[n] << endl;
                    
                }*/
                
                
                //cout << "Wz[ic-gr.n_bou_elm][ff] * tempf = " << Wz[ic-gr.n_bou_elm][ff] * tempf << endl;
                //cin.ignore();
                
                
                /*if (tempf != 0.)
                {
                
                
                cout << "tempf = " << tempf << endl;
                cout << "Wx[ic-gr.n_bou_elm][ff] = " << Wx[ic-gr.n_bou_elm][ff] << endl;
                exit(-2);
                }*/
                
            }
            
            ++ff;
        }
    }
    
    MPI_Allgatherv (MPI_IN_PLACE, localSizesNVARNDIM[rank], MPI_DOUBLE, &grad[0][0][0], localSizesNVARNDIM, displsNVARNDIM, MPI_DOUBLE, MPI_COMM_WORLD);
    
    
}