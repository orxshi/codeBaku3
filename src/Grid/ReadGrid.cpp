#include "Grid.h"

void Grid::read_grid_GMSH()
{
    /* Read grid from a file which is created by GMSH. */
    
    string temps; // for storing redundant a string.
    int tempi; // for storing redundant a number.    
    
    // read number of points.
    // open mesh file.
    in.open (meshFile);

    if (in.is_open())
    {
        // read line until "$Nodes" is found.
        do
        {
            getline (in, temps);
        }
        while (temps != "$Nodes");

        // number under "$Nodes" is the number of points.
        in >> tempi;
        pt.resize (tempi);        
    }
    else
    {
        cout << "Could not open mesh file" << endl;
        exit(-2);
    }
    
    // read points.
    for (unsigned int i=0; i<pt.size(); ++i)
    {
        in >> temps;
        in >> pt[i].dim[0];
        in >> pt[i].dim[1];
        in >> pt[i].dim[2];
        
        pt[i].belonging = id;
    }
    
    // read number of elements.
    // read line until "$Elements" is found.
    do
    {
        getline (in, temps);
    }
    while (temps != "$Elements");

    // the number under "$Elements" is total number of elements which includes ghosts and interiors.
    in >> totalNElms;
    
    // read elements.
    bool stillBoundary = true;
    int n_tags; // related to GMSH.
    int tag_count; // related to GMSH.
    int n_part_belongs;

    vector <btree> bt; // for fast cell-to-cell connectivity.
    bt.resize (pt.size() + 1);

    for (int e=0; e<totalNElms; ++e)
    {
        in >> tag_; // read tag of boundary face or cell.
        in >> gs_; // read geometric shape of boundary face or cell.        
        in >> n_tags_; // read number of GMSH tags.

        tag_count_ = 0;

        // read GMSH tags.
        while (n_tags_ > 0)
        {
            // read physical number.
            in >> phys_;
            ++ tag_count_;
            if (tag_count_ == n_tags_) break;
            in >> geo_; // read geometrical number.
            ++ tag_count_;
            if (tag_count_ == n_tags_) break;
            in >> n_part_; // read number of partitions to which element belongs.
            for (int i=0; i<n_part_; ++i)
            {
                in >> temps;
            }            
            
            break;
        }
        
        if (gs_ <= geometric_shape_t::QUAD)
        {            
            BoundaryFace bface_; // create a boundary face.
            bface_.tag = tag_; // set tag of boundary face to that of grid.
            bface_.phys = phys_;            
            
            // set geometric shape of boundary face.
            switch (tempi)
            {
                case 2:
                    bface_.geometric_shape = geometric_shape_t::TRI;
                    break;
                case 3:
                    bface_.geometric_shape = geometric_shape_t::QUAD;
                    break;
                case 4:
                    bface_.geometric_shape = geometric_shape_t::TET;
                    break;
                case 5:
                    bface_.geometric_shape = geometric_shape_t::HEX;
                    break;
                case 6:
                    bface_.geometric_shape = geometric_shape_t::PEN;
                    break;
                default:
                    cout << "unknown boundary face type in read grid" << endl;
                    exit(-2);
                    break;
            }            
            
            switch (geometric_shape)
            {
                case geometric_shape_t::TRI:            
                    cell_.nVertices = 3;
                    cell_.nFaces = 1;
                    break;
                case geometric_shape_t::HEX:            
                    cell_.nVertices = 8;
                    cell_.nFaces = 6;
                    break;
                case geometric_shape_t::TET:           
                    cell_.nVertices = 4;
                    cell_.nFaces = 4;
                    break;
                case geometric_shape_t::QUAD:            
                    cell_.nVertices = 1;
                    cell_.nFaces = 4;
                    break;
                case geometric_shape_t::PEN:            
                    cell_.nVertices = 6;
                    cell_.nFaces = 5;
                    break;                
            }
            
            bface_.vtx.reserve (bface_.nVertices);
            bface_.vtxBelo.reserve (bface_.nVertices);
        }
        else
        {            
            Cell cell_; // create a cell.
            cell_.tag = tag; // set tag of cell to that of grid.
            cell_.phys = phys_;
            
            switch (tempi)
            {                
                case 2:
                    cell_.geometric_shape = geometric_shape_t::TRI;
                    break;
                case 3:
                    cell_.geometric_shape = geometric_shape_t::QUAD;
                    break; 
                case 4:
                    cell_.geometric_shape = geometric_shape_t::TET;
                    break;
                case 5:
                    cell_.geometric_shape = geometric_shape_t::HEX;
                    break;
                case 6:
                    cell_.geometric_shape = geometric_shape_t::PEN;
                    break;
                default:
                    cout << "unknown boundary face in read grid" << endl;
                    exit(-2);
                    break;
            }
            
            switch (geometric_shape)
            {
                case geometric_shape_t::TRI:            
                    cell_.nVertices = 3;
                    cell_.nFaces = 1;
                    break;
                case geometric_shape_t::HEX:            
                    cell_.nVertices = 8;
                    cell_.nFaces = 6;
                    break;
                case geometric_shape_t::TET:           
                    cell_.nVertices = 4;
                    cell_.nFaces = 4;
                    break;
                case geometric_shape_t::QUAD:            
                    cell_.nVertices = 1;
                    cell_.nFaces = 4;
                    break;
                case geometric_shape_t::PEN:            
                    cell_.nVertices = 6;
                    cell_.nFaces = 5;
                    break;                
            }
            
            cell_.vtx.reserve (cell_.nVertices);
            cell_.vtxBelo.reserve (cell_.nVertices);
            
            // vertices.
            for (int i=0; i<cell_.nVertices; ++i)
            {
                in >> tempi;
                cell_.vtx.push_back(tempi-1);
                cell_.vtxBelo.push_back (id);                
                bt[tempi-1].insert(e);
            }
        }        
        
        // was using bt for connectivity but this time boundary faces and cells are separated. so redesign bt or connectivity function.
        
        // vertices.
        for (unsigned int i=0; i<tmpElm.nVertices; ++i)
        {
            in >> tempi;
            tmpElm.vtx.push_back(tempi-1);
            tmpElm.vtxBelo.push_back (id);
            //tmpElm.vtx.push_back ( ref(pt[tempi-1]) );
            //tmpElm.vtxIndex.push_back (tempi-1);
            bt[tempi-1].insert(e);
        }
        
        cell.push_back( move(tmpElm) );        
    }
    
    for (int c=0; c<n_bou_elm; ++c)
    {
        cell[c].ghost = true;
    }
    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        cell[c].ghost = false;
    }

    n_in_elm = totalNElms - n_bou_elm;
}

void Grid::read_ptSize ()
{
    string temps;
    int tempi;
        
    in.open (meshFile);

    if (in.is_open())
    {
        do
        {
            getline (in, temps);
        }
        while (temps != "$Nodes");

        in >> tempi;
        pt.resize (tempi);
        //ptIndex.resize (tempi);
    }
    else
    {
        cout << "Could not open mesh file" << endl;
    }
}

void Grid::read_pt()
{
    string temps;
    
    for (unsigned int i=0; i<pt.size(); ++i)
    {
        in >> temps;
        in >> pt[i].dim[0];
        in >> pt[i].dim[1];
        in >> pt[i].dim[2];
        
        pt[i].belonging = id;
    }
}

void Grid::read_elmSize ()
{
    std::string temps;
    int tempi;

    do
    {
        getline (in, temps);
    }
    while (temps != "$Elements");

    in >> totalNElms;
}

void Grid::read_elm ()
{
    bool stillBoundary = true;
    int n_tags;
    string temps;
    int tempi;
    int tag_count;
    int n_part_belongs;

    bt.resize (pt.size() + 1);

    for (int e=0; e<totalNElms; ++e)
    {
        in >> temps; // ID
        in >> tempi; // type
        
        Cell tmpElm;
        tmpElm.belonging = id;

        switch (tempi)
        {
            case 2:
                tmpElm.type = elmType_t::TRI;
                break;
            case 3:
                tmpElm.type = elmType_t::QUAD;
                break;
            case 4:
                tmpElm.type = elmType_t::TET;
                break;
            case 5:
                tmpElm.type = elmType_t::HEX;
                break;
            case 6:
                tmpElm.type = elmType_t::PEN;
                break;
            default:
                cout << "unknown boundary element type in read_elm()" << endl;
                exit(-2);
                break;
        }

        if (stillBoundary)
        {
            if (tmpElm.type > elmType_t::QUAD)
            {
                n_bou_elm = e;
                stillBoundary = false;
            }
        }

        in >> n_tags; // num of tags

        tag_count = 0;

        while (n_tags > 0)
        {
            in >> tmpElm.phys;
            ++ tag_count;
            if (tag_count == n_tags) break;
            in >> temps; // geometrical num
            ++ tag_count;
            if (tag_count == n_tags) break;
            in >> n_part_belongs; // num of partitions to which element belongs

            for (int i=0; i<n_part_belongs; ++i)
            {
                in >> temps;
            }

            break;
        }

        switch ( tmpElm.type )
        {
            case elmType_t::TRI:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::TRI) );
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::TRI));
                tmpElm.nVertices = static_cast<int>(nVertices_t::TRI);
                tmpElm.nFaces = static_cast<int>(nFaces_t::TRI);
                break;

            case elmType_t::HEX:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::HEX));
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::HEX));
                tmpElm.nVertices = static_cast<int>(nVertices_t::HEX);
                tmpElm.nFaces = static_cast<int>(nFaces_t::HEX);
                break;

            case elmType_t::TET:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::TET));
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::TET));
                tmpElm.nVertices = static_cast<int>(nVertices_t::TET);
                tmpElm.nFaces = static_cast<int>(nFaces_t::TET);
                break;

            case elmType_t::QUAD:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::QUAD));
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::QUAD));
                tmpElm.nVertices = static_cast<int>(nVertices_t::QUAD);
                tmpElm.nFaces = static_cast<int>(nFaces_t::QUAD);
                break;

            case elmType_t::PEN:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::PEN));
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::PEN));
                tmpElm.nVertices = static_cast<int>(nVertices_t::PEN);
                tmpElm.nFaces = static_cast<int>(nFaces_t::PEN);
                break;
            default:
                cout << "elm[" << e << "].type = " << static_cast<int>(tmpElm.type) << endl;
                exit(-2);
                break;
        }

        tmpElm.vtx.reserve (tmpElm.nVertices);
        tmpElm.vtxBelo.reserve (tmpElm.nVertices);
        for (unsigned int i=0; i<tmpElm.nVertices; ++i)
        {
            in >> tempi;
            tmpElm.vtx.push_back(tempi-1);
            tmpElm.vtxBelo.push_back (id);
            //tmpElm.vtx.push_back ( ref(pt[tempi-1]) );
            //tmpElm.vtxIndex.push_back (tempi-1);
            bt[tempi-1].insert(e);
        }
        
        cell.push_back( move(tmpElm) );
        //cell.push_back(tmpElm);
    }
    
    for (int c=0; c<n_bou_elm; ++c)
    {
        cell[c].ghost = true;
    }
    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        cell[c].ghost = false;
    }

    n_in_elm = totalNElms - n_bou_elm;
}

/*void Grid::read_input()
{
    int tmp;
    string tmps;
    
    ifstream in;
    string inputDir;
    string ids = to_string (id);
    inputDir = "../input_";
    inputDir.append (ids);
    inputDir.append (".dat");
    in.open(inputDir);
    
    in >> tmps; in >> meshFile;
    in >> tmps; in >> steady;
    in >> tmps; in >> cfl;
    in >> tmps; in >> useCFL;
    in >> tmps; in >> dt;
    in >> tmps; in >> tmp;
    
    if (tmp == 1)
    {
        tOrder = tsOrder_t::FIRST;
    }
    else if (tmp == 2)
    {
        tOrder = tsOrder_t::SECOND;
    }
    else
    {
        cout << "temporal order is given as: " << tmp << endl;
        exit(-2);
    }
    
    in >> tmps; in >> tmp;
    
    if (tmp == 1)
    {
        sOrder = tsOrder_t::FIRST;
    }
    else if (tmp == 2)
    {
        sOrder = tsOrder_t::SECOND;
    }
    else
    {
        cout << "spatial order is given as: " << tmp << endl;
        exit(-2);
    }
    
    in >> tmps; in >> implicit;
    in >> tmps; in >> nGaussIter;
    in >> tmps; in >> tol;
    in >> tmps; in >> maxTimeStep;
    in >> tmps; in >> nSolidBoundaries;
    
    for (int i=0; i<nSolidBoundaries; ++i)
    {
        in >> tmpd;
        solidBoundaryPhys.push_back (tmpd);
    }
    
    in >> tmps; in >> solidBoundaryPhys;
    in >> tmps; in >> phys_count;
    
    phys.resize(phys_count);
    bc.resize(phys_count);
    
    for (int i=0; i<phys_count; ++i)
    {
        in >> tmps; in >> phys[i];
    }
    
    for (int i=0; i<phys_count; ++i)
    {
        in >> tmps; in >> bc[i];
    }
    
    in >> tmps; in >> rhoInf;
    in >> tmps; in >> pInf;
    in >> tmps; in >> MachAirfoil;
    in >> tmps; in >> MachAir;
    in >> tmps; in >> aoa;
    in >> tmps; in >> alphaMean;
    in >> tmps; in >> alphaMax;
    in >> tmps; in >> kc;
    
    for (int d=0; d<N_DIM; ++d)
    {
        in >> tmps; in >> centerAirfoil[d];
    }
    
    in.close();
}*/

/*void Grid::printInput()
{
    cout << ">>> Grid: " << id << endl;
    cout << ">>> Mesh file: " << meshFile << endl;
    cout << ">>> CFL: " << cfl << endl;

    if (steady)
    {
        cout << ">>> Steady: true" << endl;
    }
    else
    {
        cout << ">>> Steady: false" << endl;
    }

    if (useCFL)
    {
        cout << ">>> Use CFL: true" << endl;
    }
    else
    {
        cout << ">>> Use CFL: false" << endl;
        cout << ">>> time step: " << dt << endl;
    }

    cout << ">>> Temporal order: " << static_cast<int>(tOrder) << endl;
    cout << ">>> Spatial order: " << static_cast<int>(sOrder) << endl;

    if (implicit)
    {
        cout << ">>> Implicit: true" << endl;
        cout << ">>> nGaussIter: " << nGaussIter << endl;
    }
    else
    {
        cout << ">>> Implicit: false" << endl;
    }

    cout << ">>> Tolerance: " << tol << endl;
    cout << ">>> maxTimeStep: " << maxTimeStep << endl;
    cout << ">>> solidBoundaryPhys: " << solidBoundaryPhys << endl;
    cout << ">>> Number of phys: " << phys_count << endl;
    
    for (int i=0; i<phys_count; ++i)
    {
        cout << ">>> phys_" << phys[i] << " : " << bcVerbose[bc[i]] << endl;
    }
    
    cout << ">>> Ref density: " << rhoInf << endl;
    cout << ">>> Ref pressure: " << pInf << endl;
    cout << ">>> MachAirfoil: " << MachAirfoil << endl;
    cout << ">>> MachAir: " << MachAir << endl;
    cout << ">>> aoa: " << aoa << endl;
    cout << ">>> alphaMean: " << alphaMean << endl;
    cout << ">>> alphaMax: " << alphaMax << endl;
    cout << ">>> kc: " << kc << endl;
    
    for (int d=0; d<N_DIM; ++d)
    {
        cout << ">>> centerAirfoil[" << d << "] : " << centerAirfoil[d] << endl;
    }
}*/

void Grid::printMeshInfo()
{
    cout << ">>> Number of cells: " << n_in_elm << endl;
}

void Grid::readInput()
{
    string tmps;
    
    ifstream in;
    string inputDir;
    string ids = to_string (id);
    inputDir = "./Grid_";
    inputDir.append (ids);
    inputDir.append ("/input.dat");
    in.open(inputDir);
    
    in >> tmps; in >> meshFile;
    in >> tmps; in >> phys_count;
    
    phys.resize(phys_count);
    bc.resize(phys_count);

    for (int i=0; i<phys_count; ++i)
    {
        in >> tmps; in >> phys[i];
    }
    
    for (int i=0; i<phys_count; ++i)
    {
        in >> tmps; in >> bc[i];
    }
    
    in >> tmps; in >> nHoles;
    holes.resize(nHoles);
    
    for (int i=0; i<nHoles; ++i)
    {
        for (int j=0; j<N_DIM; ++j)
        {
            in >> tmps; in >> holes[i].min[j];
            in >> tmps; in >> holes[i].max[j];
        }
    }
    
    in.close();
}

void Grid::printInput()
{    
    cout << meshFile << endl;
}

void Grid::read_grid ()
{
    readInput();
    printInput();
    read_ptSize ();
    read_pt();
    read_elmSize ();
    read_elm ();    
    printMeshInfo();
    in.close();
}



