#ifndef FACE_H
#define FACE_H

#include "Types.h"

class Face {
public:
    // one of the halfedges associated with this face
    HalfEdgeIter he;
    
    // id between 0 and |F|-1
    int index;
    
    // parent for tree-cotree decomposition
    FaceCIter parent;
    
    // field angle
    double beta;
    
    // field angle
    Eigen::Vector3d field;
    
    // flag for setting beta
    bool visited;
    
    // returns centroid
    Eigen::Vector3d centroid() const;
        
    // returns shared halfedge
    HalfEdgeIter sharedEdge(FaceCIter& f) const;
    
    // returns local axis associated with the triangle
    void axis(Eigen::Vector3d& x, Eigen::Vector3d& y) const;
    
    // returns face area
    double area() const;
    
    // returns normal to face
    Eigen::Vector3d normal() const;
};

#endif
