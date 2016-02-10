#include "../Grid/Grid.h"
#include "IBlank.h"

void Grid::CellADT::build (const Grid& gr)
{
    points.resize (gr.n_in_elm);
    
    for (unsigned int c=0; c<points.size(); ++c)
    {
        for (int i=0; i<ADT_DIM; ++i)
        {
            points[c].dim[i*2]   = BIG_POS_NUM;
            points[c].dim[i*2+1] = BIG_NEG_NUM;
        }

        for (const int ip: gr.cell[c+gr.n_bou_elm].vtx)
        {
            const Point& p = gr.pt[ip];
            
            for (int i=0; i<ADT_DIM; ++i)
            {
                points[c].dim[i*2]   = min (p.dim[i], points[c].dim[i*2]);
                points[c].dim[i*2+1] = max (p.dim[i], points[c].dim[i*2+1]);
            }

            points[c].vertices.push_back (p.dim);
        }

        points[c].idx = c+gr.n_bou_elm;
    }
    
    ADT::build();
}

void Grid::setWallDistance (int phys)
{
    for (int c=0; c<cell.size(); ++c)
    {
        double d = BIG_POS_NUM;
        
        for (int g=0; g<n_bou_elm; ++g)
        {
            if (cell[g].phys == phys)
            {
                CVector tmpv = cell[g].cnt - cell[c].cnt;

                d = min ( d, mag (tmpv) );
            }
        }

        cell[c].wallDistance = d;
    }
}

Iblank::Iblank ()
{
    // default
    cellCriter = cellCriter_t::WALL;

    ifstream in;
    in.open("iblank.dat");
    
    if (in.is_open())
    {
        string tmps;
        int tmpi;
        in >> tmps; in >> tmpi;
        
        cellCriter = static_cast<cellCriter_t> (tmpi);
    }
    
    in.close();
}

void Iblank::identify (Grid& grAct, Grid& grPas)
{
    function<int(Cell&)> getIndex = [&] (Cell& cll)
    {
        ADT::ADTPoint vec;
        CVector cnt;
        
        cnt = cll.cnt;
        
        for (int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = cnt[i];
            vec.dim[i*2+1] = cnt[i];
        }

        /*vec.dim[0] = cnt[0];
        vec.dim[2] = cnt[1];
        vec.dim[4] = cnt[2];

        vec.dim[1] = cnt[0];
        vec.dim[3] = cnt[1];
        vec.dim[5] = cnt[2];*/        
        
        return grPas.cellADT.search (vec);
    };

    for (int c=grAct.n_bou_elm; c<grAct.cell.size(); ++c)
    {
        
    
        int index;        
        Cell& cll = grAct.cell[c];
        
        
        
        if (cll.iBlank == iBlank_t::UNDEFINED)
        {            
            
        
            index = getIndex (cll);
            
            if (index == -1)
            {
                bool insideHole;
                
                if (grPas.nHoles > 0)
                {
                    for (int i=0; i<grPas.nHoles; ++i)
                    {
                        for (int j=0; j<N_DIM; ++j)
                        {
                            if (cll.cnt[j] >= grPas.holes[i].min[j] && cll.cnt[j] <= grPas.holes[i].max[j])
                            {
                                insideHole = true;
                            }
                            else
                            {
                                insideHole = false;
                                break;
                            }
                        }
                    }                    
                }
                
                if (grPas.nHoles > 0 && insideHole)
                {
                    cll.iBlank = iBlank_t::HOLE;                    
                }
                else
                {
                    cll.iBlank = iBlank_t::FIELD;
                }
            }
            else
            {
            
            
                if (grPas.cell[index].iBlank == iBlank_t::FIELD)
                {
                    cll.iBlank = iBlank_t::FRINGE;
                    
                    /*if (grPas.cell[index].iBlank == iBlank_t::FIELD && grPas.id == 0)
                    {
                        cout << "case B" << endl;
                        exit(-2);
                    }*/
                    
                    cll.donor = &grPas.cell[index];
                    grPas.cell[index].receiver.push_back (&cll);
                }
                else if (grPas.cell[index].iBlank == iBlank_t::FRINGE || grPas.cell[index].iBlank == iBlank_t::HOLE)
                {
                    
                
                    cll.iBlank = iBlank_t::FIELD;
                    
                    /*if (cll.iBlank == iBlank_t::FIELD && grAct.id == 0)
                    {
                        cout << "case D" << endl;
                        cout << "state = " << static_cast<int> (grPas.cell[index].iBlank) << endl;
                        exit(-2);
                    }*/
                }
                else if (grPas.cell[index].iBlank == iBlank_t::UNDEFINED)
                {
                    double var1;
                    double var2;
                    
                    
                    
                    if (cellCriter == cellCriter_t::WALL)
                    {
                        var1 = grPas.cell[index].wallDistance;
                        var2 = cll.wallDistance;
                        
                            
                    }
                    else if (cellCriter == cellCriter_t::SIZE)
                    {
                        var1 = grPas.cell[index].vol;
                        var2 = cll.vol;
                    }
                    
                    if (var1 < var2)
                    {
                        cll.iBlank = iBlank_t::FRINGE;
                        grPas.cell[index].iBlank = iBlank_t::FIELD;
                        cll.donor = &grPas.cell[index];
                        grPas.cell[index].receiver.push_back (&cll);
                    }
                    else
                    {
                        cll.iBlank = iBlank_t::FIELD;
                        grPas.cell[index].iBlank = iBlank_t::FRINGE;
                        grPas.cell[index].donor = &cll;
                        cll.receiver.push_back (&grPas.cell[index]);
                    }
                    
                    /*if (grPas.cell[index].iBlank == iBlank_t::FRINGE && grPas.id == 1)
                    {
                        cout << "case A" << endl;
                        exit(-2);
                    }
                    
                    if (cll.iBlank == iBlank_t::FIELD && grPas.id == 0)
                    {
                        cout << "case C" << endl;
                        exit(-2);
                    }*/
                    
                }
                else
                {
                    cout << "undefined behavior in Grid::identifyIBlank(...)" << endl;
                    cout << "state = " << static_cast<int> (grPas.cell[index].iBlank) << endl;
                    exit(-2);
                }
            }
        }
        
    }
    
    for (int c=0; c<grAct.n_bou_elm; ++c)
    {
        int index;
        Cell& cll = grAct.cell[c];
        
        if (cll.iBlank == iBlank_t::UNDEFINED)
        {
            if (cll.fringeBou == fringeBou_t::YES)
            {
                index = getIndex (cll);

                if (index == -1)
                {
                    cll.iBlank == iBlank_t::NA;
                }
                else
                {
                    if (grPas.cell[index].iBlank == iBlank_t::FIELD)
                    {
                        cll.iBlank = iBlank_t::FRINGE;
                        cll.donor = &grPas.cell[index];
                        grPas.cell[index].receiver.push_back (&cll);
                    }
                    else if (grPas.cell[index].iBlank == iBlank_t::FRINGE)
                    {                        
                        cll.iBlank = iBlank_t::FRINGE;
                        grPas.cell[index].iBlank = iBlank_t::FIELD;
                        
                        for (int r=0; r<grPas.cell[index].donor->receiver.size(); ++r)
                        {
                            if (grPas.cell[index].donor->receiver[r] == &grPas.cell[index])
                            {
                                grPas.cell[index].donor->receiver[r] = NULL;
                                grPas.cell[index].donor->receiver.erase (grPas.cell[index].donor->receiver.begin() + r);
                                break;
                            }
                        }
                        
                        grPas.cell[index].donor = NULL;
                        grPas.cell[index].receiver.push_back (&cll);
                        cll.donor = &grPas.cell[index];
                        //cout << "unset situation in Grid::identifyIBlank(...)" << endl;
                        //exit(-2);
                    }
                    else if (grPas.cell[index].iBlank == iBlank_t::UNDEFINED)
                    {
                        cll.iBlank = iBlank_t::FRINGE;
                        cll.donor = &grPas.cell[index];
                        grPas.cell[index].receiver.push_back (&cll);
                        grPas.cell[index].iBlank = iBlank_t::FIELD;
                    }
                    else
                    {
                        cout << "second undefined behavior in Grid::identifyIBlank(...)" << endl;
                        cout << "iBlank = " << static_cast<int>(grPas.cell[index].iBlank) << endl;
                        cout << "bc = " << static_cast<int>(cll.bc) << endl;
                        cout << "cell[c].cnt[0] = " << cll.cnt[0] << endl;
                        cout << "cell[c].cnt[1] = " << cll.cnt[1] << endl;
                        cout << "cell[c].cnt[2] = " << cll.cnt[2] << endl;
                        cout << "cell[c].belonging = " << cll.belonging << endl;
                        cout << "c = " << c << endl;
                        exit(-2);
                    }
                }
            }
            else if (cll.fringeBou == fringeBou_t::NO)
            {
                cll.iBlank == iBlank_t::NA;
            }
            else if (cll.fringeBou == fringeBou_t::UNDEFINED)
            {
                cout << "undefined fringeBou_t in Grid::identifyIBlank(...)" << endl;
                exit(-2);
            }
        }
    }
}

bool Grid::CellADT::compareFunction (const Node* node, const ADTPoint& targetPoint)
{
    bool inside;
    CVector tempVector;
    unsigned int size = node->p->vertices.size();
    double frac[size];
    vector<CVector> xc(size);

    for (unsigned int v=0; v<size; ++v)
    {
        xc[v] = node->p->vertices[v];
    }
    
    for (int i=0; i<ADT_DIM; ++i)
    {
        tempVector[i] = targetPoint.dim[i*2];
    }
    
    /*tempVector[0] = targetPoint.dim[0];
    tempVector[1] = targetPoint.dim[2];
    tempVector[2] = targetPoint.dim[4];*/
    
    osInterpolants (xc, tempVector, size, frac);

    inside = true;
    for (unsigned int v=0; v<size; ++v)
    {
        if (frac[v]<=0 || frac[v]>=1)
        {
            inside = false;
            break;
        }
    }

    return inside;
}

bool Grid::CellADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
{
    bool insideCube = true;

    for (int d=0; d<ADT_DIM; ++d)
    {
        if (!(node->p->dim[d*2] <= targetPoint.dim[d*2]) || !(node->p->dim[d*2+1] >= targetPoint.dim[d*2+1]))
        {
            insideCube = false;
            break;
        }
    }

    return insideCube;
}

void Iblank::interpolate (Grid& gr, Gradient& gradient, Grid& ogr)
{
    for (int ic=0; ic<gr.cell.size(); ++ic)
    //for (Cell& cll: cell)
    {
        Cell& cll = gr.cell[ic];
    
        //cll.interpolate();        
        
        CVector dis;        
    
        if ( cll.iBlank == iBlank_t::FRINGE )
        {
            for (int i=0; i<N_VAR; ++i)
            {
                dis = cll.cnt - cll.donor->cnt;
                
                if (ic<gr.n_bou_elm)
                {
                    cll.prim[i] = cll.donor->prim[i];
                    if ( isnan(cll.prim[i]) ) { cout << "cll.prim[i] is NAN in Solver::getRes()" << endl; exit(-2); }
                }
                else
                {
                    //cll.prim[i] = cll.donor->prim[i] + dotP(gradient.grad[cll.donor-&ogr.cell[0]][i], dis);
                    cll.prim[i] = cll.donor->prim[i];                    
                    if ( isnan(cll.donor->prim[i]) ) { cout << "cll.prim[i] is NAN in Solver::getRes() for cells" << endl; exit(-2); }
                    //if ( isnan(dotP(gradient.grad[cll.donor-&ogr.cell[0]][i], dis)) ) { cout << "grad is NAN in Solver::getRes() for cells" << endl; exit(-2); }
                }
            }

            cll.prim_to_cons();
        }
        
        
    }
    gr.apply_BCs();
}

void Iblank::treatFieldIslands (Grid& grAct)
{
    for (int c=grAct.n_bou_elm; c<grAct.cell.size(); ++c)
    {        
        int nFieldNeis = 0;
    
        Cell& cll = grAct.cell[c];
        
        if (cll.iBlank == iBlank_t::FIELD)
        {
            for (int n:cll.nei)
            {
                if (grAct.cell[n].iBlank == iBlank_t::FIELD)
                {
                    ++nFieldNeis;
                }
            }
            
            if (nFieldNeis == 0) // island found
            {
                int iR;
                double dis = BIG_POS_NUM;
            
                cll.iBlank = iBlank_t::FRINGE;
                
                for (int r=0; r<cll.receiver.size(); ++r)
                {                
                    cll.receiver[r]->iBlank = iBlank_t::FIELD;
                    cll.receiver[r]->donor = NULL;
                    double tdis = mag(cll.cnt - cll.receiver[r]->cnt);
                    if (tdis < dis)
                    {
                        iR = r;
                        dis = tdis;
                    }
                }
                
                if (cll.receiver.size() > 0)
                {
                    cll.receiver[iR]->receiver.push_back (&cll);
                    cll.donor = cll.receiver[iR];
                    
                    for (int r=0; r<cll.receiver.size(); ++r)
                    {
                        cll.receiver[r] = NULL;
                    }            
                    cll.receiver.clear();                
                }
            }
        }
    }
}

void Iblank::treatFringeIslands (Grid& grAct)
{
    for (int c=grAct.n_bou_elm; c<grAct.cell.size(); ++c)
    {        
        int nFringeNeis = 0;
    
        Cell& cll = grAct.cell[c];
        
        if (cll.iBlank == iBlank_t::FRINGE)
        {
            for (int n:cll.nei)
            {
                if (grAct.cell[n].iBlank == iBlank_t::FRINGE)
                {
                    ++nFringeNeis;
                }
            }
            
            if (nFringeNeis == 0) // island found
            {
                int iR;
                double dis = BIG_POS_NUM;
            
                cll.iBlank = iBlank_t::FIELD;
                
                if (cll.donor != NULL)
                {   
                    for (int r=0; r<cll.donor->receiver.size(); ++r)
                    {
                        if (cll.donor->receiver[r] == &cll)
                        {
                            cll.donor->receiver.erase(cll.donor->receiver.begin() + r);
                            break;
                        }
                    }
                    
                    cll.donor = NULL;
                }
                else
                {
                    cout << "has no donor" << endl;
                    exit(-2);
                }
            }
        }
    }
}

void Iblank::treatVoidAreas (Grid& grAct)
{
    for (int c=grAct.n_bou_elm; c<grAct.cell.size(); ++c)
    {        
        int nFringeNeis = 0;
        int nHoleNeis = 0;
    
        Cell& cll = grAct.cell[c];
        
        if (cll.iBlank == iBlank_t::FRINGE)
        {
            for (int n:cll.nei)
            {
                if (grAct.cell[n].iBlank == iBlank_t::FRINGE)
                {
                    ++nFringeNeis;
                }
                else if (grAct.cell[n].iBlank == iBlank_t::HOLE)
                {
                    ++nHoleNeis;
                }
            }
            
            if (nFringeNeis == 1 && nHoleNeis == 0)
            {
                int iR;                
            
                cll.iBlank = iBlank_t::FIELD;
                
                if (cll.donor != NULL)
                {   
                    for (int r=0; r<cll.donor->receiver.size(); ++r)
                    {
                        if (cll.donor->receiver[r] == &cll)
                        {
                            cll.donor->receiver.erase(cll.donor->receiver.begin() + r);
                            break;
                        }
                    }
                    
                    cll.donor = NULL;
                }
                else
                {
                    cout << "has no donor" << endl;
                    exit(-2);
                }
            }
        }
    }
}