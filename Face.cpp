#include "Face.h"
#include "HalfEdge.h"
#include "Vertex.h"

Eigen::Vector3d Face::centroid() const
{
    return (he->vertex->position +
            he->next->vertex->position +
            he->next->next->vertex->position) / 3.0;
}

HalfEdgeIter Face::sharedEdge(FaceCIter& f) const
{
    HalfEdgeIter h = he;
    do {
        if (h->flip->face == f) return h;
        h = h->next;
        
    } while (h != he);
    
    std::cout << "should not happen" << std::endl;
    // this should not happen
    return h;
}

void Face::axis(Eigen::Vector3d& x, Eigen::Vector3d& y) const
{
    Eigen::Vector3d a = he->vertex->position;
    Eigen::Vector3d b = he->next->vertex->position;
    Eigen::Vector3d c = he->next->next->vertex->position;
    
    x = (b - a).normalized();
    
    // Gram-Schmidt
    y = c - a;
    y = (y - y.dot(x)*x).normalized();
}

double Face::area() const
{
    return 0.5 * normal().norm();
}

Eigen::Vector3d Face::normal() const
{
    Eigen::Vector3d a = he->vertex->position;
    Eigen::Vector3d b = he->next->vertex->position;
    Eigen::Vector3d c = he->next->next->vertex->position;
    
    Eigen::Vector3d v1 = a - b;
    Eigen::Vector3d v2 = c - b;
    
    return v1.cross(v2);
}
