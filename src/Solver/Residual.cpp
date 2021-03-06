#include "Solver.h"

void Solver::getMaxRes (Grid& gr)
{
    Vector<N_VAR> res;
    //res.fill(BIG_NEG_NUM); // for max criterion
    maxRes.fill(BIG_NEG_NUM); // for max criterion
    
    // calculate cll.R with updated values. updateVars should be called before this
    /*if (sOrder == 2)
    {
        gradient.leastSquaresGrad (gr);
    }
    
    roe.roeflx (gr, limiter, M0, M1, gradient);
    cout << "called roe in maxres" << endl;*/
    //cin.ignore();
    
    if (tOrder == 1)
    {
        if (steady)
        {
            for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
            {
                Cell& cll = gr.cell[ic];
            
                if (cll.iBlank == iBlank_t::FIELD)
                //if (cll.iBlank == iBlank_t::FRINGE)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        //double RHS = cll.R[i];
                        double RHS = cll.dQ[i];
                        maxRes[i] = max( fabs(RHS), maxRes[i] );
                        
                        
                    }
                }
            }
        }
        else
        {
            for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
            {
                Cell& cll = gr.cell[ic];
            
                if (cll.iBlank == iBlank_t::FIELD)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        double LHS = (cll.cons[i] - cll.old_cons[i]) * cll.vol / dt;
                        double RHS = cll.R[i];

                        //res[i] += pow(LHS - RHS, 2.);
                        maxRes[i] = max( fabs(LHS - RHS), maxRes[i] );
                    }
                }
            }
        }
    }
    else
    {
        cout << "time derivative is not 1 or 2 in Solver::setExpRes(...)" << endl;
        exit(-2);
    }
    
    /*maxRes = BIG_NEG_NUM;
    for (int i=0; i<N_VAR; ++i)
    {
        maxRes = max(res[i],maxRes);
    }*/
    
    for (int i=0; i<N_VAR; ++i)
    {
        if ( isnan(maxRes[i]) ) { cout << "maxRes " << i << " is NAN in Solver::getRes()" << endl; exit(-2); }
        if ( isinf(maxRes[i]) ) { cout << "maxRes " << i << " is INF in Solver::getRes()" << endl; exit(-2); }
    }
}

void Solver::getRmsRes (Grid& gr)
{
    Vector<N_VAR> res;
    //res.fill(0.);
    rmsRes.fill(0.);
    int nField = 0;
    
    /*if (sOrder == 2)
    {
        gradient.leastSquaresGrad (gr);
    }
    
    roe.roeflx (gr, limiter, M0, M1, gradient);
    cout << "called roe in rmsres" << endl;*/
    //cin.ignore();
    
    if (tOrder == 1)
    {
        if (steady)
        {
            for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
            {
                Cell& cll = gr.cell[ic];
            
                if (cll.iBlank == iBlank_t::FIELD)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        //double RHS = cll.R[i];
                        double RHS = cll.dQ[i];
                        rmsRes[i] += pow(RHS,2.);
                        
                        
                    }
                    
                    ++nField;
                }
            }
        }
        else
        {
            for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
            {
                Cell& cll = gr.cell[ic];
            
                if (cll.iBlank == iBlank_t::FIELD)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        double LHS = (cll.cons[i] - cll.old_cons[i]) * cll.vol / dt;
                        double RHS = cll.R[i];

                        rmsRes[i] += pow(LHS - RHS, 2.);
                    }
                    
                    ++nField;
                }
            }
        }
    }
    else
    {
        cout << "time derivative is not 1 or 2 in Solver::setExpRes(...)" << endl;
        exit(-2);
    }
    
    for (int i=0; i<N_VAR; ++i)
    {
        rmsRes[i] /= nField;
        rmsRes[i] = sqrt (rmsRes[i]);
    }
    
    /*rmsRes = BIG_NEG_NUM;
    for (int i=0; i<N_VAR; ++i)
    {
        rmsRes = max(res[i],rmsRes);
    }*/
    
    for (int i=0; i<N_VAR; ++i)
    {
        if ( isnan(rmsRes[i]) ) { cout << "rmsRes " << i << " is NAN in Solver::getRes()" << endl; exit(-2); }
        if ( isinf(rmsRes[i]) ) { cout << "rmsRes " << i << " is INF in Solver::getRes()" << endl; exit(-2); }
    }
}

void Solver::getResiduals (Grid& gr)
{
    rmsRes.fill(0.);
    maxRes.fill(BIG_NEG_NUM);
    int nField = 0;
    
    if (tOrder == 1)
    {
        if (steady)
        {
            for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
            {
                Cell& cll = gr.cell[ic];
            
                if (cll.iBlank == iBlank_t::FIELD)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        double RHS = cll.R[i];
                        rmsRes[i] += pow(RHS,2.);
                        maxRes[i] = max( fabs(RHS), maxRes[i] );
                    }
                    
                    ++nField;
                }
            }
        }
        else
        {
            for (int ic=gr.n_bou_elm; ic<gr.cell.size(); ++ic)
            {
                Cell& cll = gr.cell[ic];
            
                if (cll.iBlank == iBlank_t::FIELD)
                {
                    for (int i=0; i<N_VAR; ++i)
                    {
                        double LHS = (cll.cons[i] - cll.old_cons[i]) * cll.vol / dt;
                        double RHS = cll.R[i];

                        rmsRes[i] += pow(LHS - RHS, 2.);
                        maxRes[i] = max( fabs(LHS - RHS), maxRes[i] );
                    }
                    
                    ++nField;
                }
            }
        }
    }
    else
    {
        cout << "only tOrder=1 is supported in Solver::getResiduals (Grid& gr)" << endl;
        exit(-2);
    }
    
    for (int i=0; i<N_VAR; ++i)
    {
        rmsRes[i] /= nField;
        rmsRes[i] = sqrt (rmsRes[i]);
    }
    
    for (int i=0; i<N_VAR; ++i)
    {
        if ( isnan(rmsRes[i]) ) { cout << "rmsRes " << i << " is NAN in Solver::getResiduals()" << endl; exit(-2); }
        if ( isnan(maxRes[i]) ) { cout << "maxRes " << i << " is NAN in Solver::getResiduals()" << endl; exit(-2); }
        if ( isinf(rmsRes[i]) ) { cout << "rmsRes " << i << " is INF in Solver::getResiduals()" << endl; exit(-2); }
        if ( isinf(maxRes[i]) ) { cout << "rmsRes " << i << " is INF in Solver::getResiduals()" << endl; exit(-2); }
    }
}