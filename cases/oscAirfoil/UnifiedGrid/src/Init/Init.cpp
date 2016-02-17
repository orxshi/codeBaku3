#include "Init.h"

void OscInit::init (Grid& gr)
{    
    //cout << "in oscinit" << endl;

    double aoaRad = aoa * DEG_TO_RAD;
    
    double p[5];
    
    p[0] = rhoInf;
    p[1] = Mach * cos(aoaRad);
    p[2] = Mach * sin(aoaRad);
    p[3] = 0.;
    p[4] = pInf;
    
    //cout << "cp1 oscinit" << endl;

    for (Cell& cll: gr.cell)
    {
        for (int i=0; i<N_VAR; ++i)
        {
            cll.prim[i] = p[i];
        }
        
        cll.prim_to_cons();
        
        cll.old_cons = cll.cons;
        cll.oldold_cons = cll.cons;
        
        //cll.iBlank = iBlank_t::FIELD;
    }
    
    //cout << "cp2 oscinit" << endl;
    
    gr.set_BCs();
    
    //cout << "cp3 oscinit" << endl;
    
    gr.apply_BCs();
    
    //cout << "cp4 oscinit" << endl;
    
    for (int c=0; c<gr.n_bou_elm; ++c)
    {
        Cell& cll = gr.cell[c];
        
        if (cll.bc == BC::DIRICHLET)
        {
            cll.fringeBou = fringeBou_t::YES;
        }
        else
        {
            cll.fringeBou = fringeBou_t::NO;
        }
    }
    
    //cout << "cp5 oscinit" << endl;
}