#ifndef VERTEX_H
#define VERTEX_H

#include "Types.h"

class Vertex {
public:
    // outgoing halfedge
    HalfEdgeIter he;
    
    // location in 3d
    Eigen::Vector3d position;
    
    // id between 0 and |V|-1
    int index;
    
    // target holonomy / 2pi
    double k;
    
    // parent for tree-cotree decomposition
    VertexCIter parent;
    
    // checks if vertex is contained in any edge or face
    bool isIsolated() const;
    
    // 2pi - ∑ø
    double angleDefect() const;
};

#endif
