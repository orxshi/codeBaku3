#include "Grid.h"

void Face::set_faceVertices (vector <int>& facevtx, const Cell& cell, int index)
{
    switch (cell.shape)
    {
        case GeometricShape::Hex:
            
            switch (index)
            {
                case 0:                    
                    facevtx.push_back ( cell.vtx[0] );
                    facevtx.push_back ( cell.vtx[3] );
                    facevtx.push_back ( cell.vtx[2] );
                    facevtx.push_back ( cell.vtx[1] );
                    break;

                case 1:
                    facevtx.push_back ( cell.vtx[0] );
                    facevtx.push_back ( cell.vtx[1] );
                    facevtx.push_back ( cell.vtx[5] );
                    facevtx.push_back ( cell.vtx[4] );
                    break;

                case 2:
                    facevtx.push_back ( cell.vtx[1] );
                    facevtx.push_back ( cell.vtx[2] );
                    facevtx.push_back ( cell.vtx[6] );
                    facevtx.push_back ( cell.vtx[5] );
                    break;

                case 3:
                    facevtx.push_back ( cell.vtx[2] );
                    facevtx.push_back ( cell.vtx[3] );
                    facevtx.push_back ( cell.vtx[7] );
                    facevtx.push_back ( cell.vtx[6] );
                    break;

                case 4:
                    facevtx.push_back ( cell.vtx[0] ); 
                    facevtx.push_back ( cell.vtx[4] ); 
                    facevtx.push_back ( cell.vtx[7] ); 
                    facevtx.push_back ( cell.vtx[3] ); 
                    break;

                case 5:
                    facevtx.push_back ( cell.vtx[4] ); 
                    facevtx.push_back ( cell.vtx[5] ); 
                    facevtx.push_back ( cell.vtx[6] ); 
                    facevtx.push_back ( cell.vtx[7] ); //face.vtxIndex.push_back (elm.vtxIndex[7]);
                    break;
            }
            break;

        case GeometricShape::Tet:

            switch (index)
            {
                
                
                case 0:
                    facevtx.push_back ( cell.vtx[0] ); 
                    facevtx.push_back ( cell.vtx[2] ); 
                    facevtx.push_back ( cell.vtx[1] ); 
                    break;

                case 1:
                    facevtx.push_back ( cell.vtx[0] ); 
                    facevtx.push_back ( cell.vtx[1] ); 
                    facevtx.push_back ( cell.vtx[3] ); 
                break;

                case 2:
                    facevtx.push_back ( cell.vtx[1] ); 
                    facevtx.push_back ( cell.vtx[2] ); 
                    facevtx.push_back ( cell.vtx[3] ); 
                break;

                case 3:
                    facevtx.push_back ( cell.vtx[2] ); 
                    facevtx.push_back ( cell.vtx[0] ); 
                    facevtx.push_back ( cell.vtx[3] ); 
                break;
            }
            break;

        case GeometricShape::Pen:

        switch(index)
        {
            
            
            case 0:
                face.vtxpush_back ( cell.vtx[0] ); 
                face.vtxpush_back ( cell.vtx[1] ); 
                face.vtxpush_back ( cell.vtx[4] ); 
                face.vtxpush_back ( cell.vtx[3] ); 
                break;

            case 1:
                face.vtxpush_back ( cell.vtx[1] ); 
                face.vtxpush_back ( cell.vtx[2] ); 
                face.vtxpush_back ( cell.vtx[5] ); 
                face.vtxpush_back ( cell.vtx[4] ); 
                break;

            case 2:
                face.vtxpush_back ( cell.vtx[2] ); 
                face.vtxpush_back ( cell.vtx[0] ); 
                face.vtxpush_back ( cell.vtx[3] ); 
                face.vtxpush_back ( cell.vtx[5] ); 
                break;

            case 3:
                face.vtxpush_back ( cell.vtx[0] ); 
                face.vtxpush_back ( cell.vtx[2] );
                face.vtxpush_back ( cell.vtx[1] );
                break;

            case 4:
                face.vtxpush_back ( cell.vtx[3] ); 
                face.vtxpush_back ( cell.vtx[4] ); 
                face.vtxpush_back ( cell.vtx[5] ); 
                break;
        }
        break;
        
        default:
            cout << "undefined cell shape" << endl;
            exit(-2);
            break;
    }
}

void Grid::cell_to_cell_connectivity()
{    
    vector <faceNumCheck> fnc (totalNElms);    
    double tmpd;
    
    // loop through cells.
    for (int c=0; c<cell.size(); ++c)   
    {        
        cell[c].face.reserve (cell[c].shape.nFaces);
        
        // loop through cell face indices.
        for (int j=0; j<cell[c].shape.nFaces; ++j)
        {
            if (!pt[face_.vtx[0]].bfaceTag.empty())
            {
                vector<int> facevtx;                
                set_faceVertices (facevtx, cell[c], j);
                BoundaryFace face_(move(facevtx)); // bface is already created in readgrid.
        
                // loop through boundary face tags of first vertex of face.
                for (int nei: pt[face_.vtx[0]].bfaceTag)
                {
                    if (nei == e)
                    {
                        continue;
                    }
                    
                    int sig = 0;
                    for (int n=face_.vtx[1]; n<face_.vtx.back(); ++n)
                    {
                        if (std::find(pt[n].bfaceTag.begin(), pt[n].bfaceTag.end(), nei) != pt[n].bfaceTag.end())
                        {
                            ++sig;
                        }
                    }
                    
                    if (sig == face_.vtx.size()-1)
                    {                            
                        cell[c].nei.push_back(-1);
                        face_.nei.push_back(c);
            
                        // set volumes
                        tmpd = face_.cnt[0] * face_.area[0];
                        cell[c].vol += tmpd;                                                    
                        
                        bface.push_back(face_);
                        
                        cell[c].face.push_back( bface.size()-1 );                            

                        break;                        
                    }
                }
                
                continue; // continue with the next face of cell.            
            }
            
            InteriorFace face_;
            set_faceVertices (face_, cell[e], j);
            
            // loop through cell tags of first vertex of face.
            for (int nei: pt[face_.vtx[0]].cellTag)
            {            
                if (nei == e || fnc[c].nmap[nei] == 1)
                {
                    continue;
                }
                
                int sig = 0;
                for (int n=face_.vtx[1]; n<face_.vtx.back(); ++n)
                {
                    if (std::find(pt[n].cellTag.begin(), pt[n].cellTag.end(), nei) != pt[n].cellTag.end())
                    {
                        ++sig;
                    }
                }
                
                if (sig == face_.vtx.size()-1)
                {
                    fnc[c].nmap[nei] = 1;
                    fnc[nei].nmap[e] = 1;
                    
                    cell[c]  .nei.push_back (nei);
                    cell[nei].nei.push_back (c);
                    
                    face_.nei.push_back (c);
                    face_.nei.push_back (nei);
                    
                    face_.set_area (pt);
                    face_.set_centroid(pt);
        
                    // set volumes
                    tmpd = face_.cnt[0] * face_.area[0];
                    cell[c]  .vol += tmpd;
                    cell[nei].vol -= tmpd;
                    
                    iface.push_back(face_);
                    
                    cell[c]  .face.push_back( face.size()-1 );
                    cell[nei].face.push_back( face.size()-1 );

                    break;                    
                }
            }
        }
    }    
}

void Cell::calcVolume ()
{
    double tmpd;
    
    for (const Face& f: face)
    {
        tmpd = f.cnt[1] * f.area[1];
        
        if (f.nei.size() == 1)
        {
            cell[f.nei[0]].vol += tmpd;
        }
        else if (f.nei.size() == 2)
        {
            cell[f.nei[0]].vol += tmpd;
            cell[f.nei[1]].vol -= tmpd;
        }
        else
        {
            cout << "!!! Error: f.nei.size() != 1 or 2" << endl;
            exit(-2);
        }
    }
}
