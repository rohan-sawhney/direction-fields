#include "HalfEdge.h"
#include "Vertex.h"
#include "Face.h"

double HalfEdge::cotan() const
{    
    Eigen::Vector3d p0 = vertex->position;
    Eigen::Vector3d p1 = next->vertex->position;
    Eigen::Vector3d p2 = next->next->vertex->position;
    
    Eigen::Vector3d v1 = p2 - p1;
    Eigen::Vector3d v2 = p2 - p0;
    
    return v1.dot(v2) / v1.cross(v2).norm();
}

bool HalfEdge::inPrimalSpanningTree() const
{
    VertexCIter v1 = vertex;
    VertexCIter v2 = flip->vertex;
    
    return v1->parent == v2 || v2->parent == v1;
}

bool HalfEdge::inDualSpanningTree() const
{
    FaceCIter f1 = face;
    FaceCIter f2 = flip->face;
    
    return f1->parent == f2 || f2->parent == f1;
}