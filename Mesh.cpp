#include "Mesh.h"
#include "MeshIO.h"
#include <queue>

#define EPSILON 10-6

Mesh::Mesh()
{
    
}

Mesh::Mesh(const Mesh& mesh)
{
    *this = mesh;
}

bool Mesh::read(const std::string& fileName)
{
    std::ifstream in(fileName.c_str());

    if (!in.is_open()) {
        std::cerr << "Error: Could not open file for reading" << std::endl;
        return false;
    }
    
    bool readSuccessful = false;
    if ((readSuccessful = MeshIO::read(in, *this))) {
        normalize();
        setup();
    }
    
    return readSuccessful;
}

bool Mesh::write(const std::string& fileName) const
{
    std::ofstream out(fileName.c_str());
    
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return false;
    }
    
    MeshIO::write(out, *this);
    
    return false;
}

int Mesh::eulerCharacteristic() const
{
    return (int)(vertices.size() - edges.size() + faces.size());
}

void Mesh::buildContractibleCycles(std::vector<Cycle>& basisCycles)
{
    for (VertexCIter v = vertices.begin(); v!= vertices.end(); v++) {
        // add halfedges to cycle
        Cycle cycle;
        HalfEdgeIter he = v->he;
        do {
            cycle.push_back(he);
            he = he->flip->next;
                
        } while (he != v->he);
        
        basisCycles.push_back(cycle);
    }
}

void Mesh::buildPrimalSpanningTree()
{
    // initialize parents
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->parent = v;
    }
    
    // build tree
    VertexIter root = vertices.begin();
    std::queue<VertexIter> queue;
    queue.push(root);
    
    while (!queue.empty()) {
        VertexCIter v1 = queue.front();
        queue.pop();
        
        HalfEdgeCIter he = v1->he;
        do {
            VertexIter v2 = he->flip->vertex;
            
            if (v2 != root && v2->parent == v2) {
                v2->parent = v1;
                queue.push(v2);
            }
            
            he = he->flip->next;
            
        } while (he != v1->he);
    }
}

void Mesh::buildDualSpanningTree()
{
    // initialize parents
    for (FaceIter f = faces.begin(); f != faces.end(); f++) {
        f->parent = f;
    }
    
    // build tree
    FaceIter root = faces.begin();
    std::queue<FaceIter> queue;
    queue.push(root);
    
    while (!queue.empty()) {
        FaceCIter f1 = queue.front();
        queue.pop();
        
        HalfEdgeCIter he = f1->he;
        do {
            if (!he->inPrimalSpanningTree()) {
                FaceIter f2 = he->flip->face;
                
                if (f2 != root && f2->parent == f2) {
                    f2->parent = f1;
                    queue.push(f2);
                }
            }
            
            he = he->next;
            
        } while (he != f1->he);
    }
}

void Mesh::buildNonContractibleCycles(std::vector<Cycle>& basisCycles)
{
    generators.clear();
    buildPrimalSpanningTree();
    buildDualSpanningTree();
    
    for (EdgeCIter e = edges.begin(); e != edges.end(); e++) {
        
        HalfEdgeIter he = e->he;
        if (!he->inPrimalSpanningTree() && !he->inDualSpanningTree()) {
            Cycle cycle;
            cycle.push_back(he);
                
            // trace loop in both directions
            Cycle temp1;
            FaceCIter f = he->flip->face;
            while (f != f->parent) {
                FaceCIter fp = f->parent;
                temp1.push_back(f->sharedEdge(fp));
                f = fp;
            }
            
            Cycle temp2;
            f = he->face;
            while (f != f->parent) {
                FaceCIter fp = f->parent;
                temp2.push_back(f->sharedEdge(fp));
                f = fp;
            }

            // remove common edges
            int t1 = (int)temp1.size()-1;
            int t2 = (int)temp2.size()-1;
            while (temp1[t1] == temp2[t2]) {
                t1--;
                t2--;
            }
            
            for (int i = 0; i <= t1; i++) {
                cycle.push_back(temp1[i]);
            }
            
            for (int i = t2; i >= 0; i--) {
                cycle.push_back(temp2[i]->flip);
            }
            
            generators.push_back(cycle);
        }
    }
    
    basisCycles.insert(basisCycles.end(), generators.begin(), generators.end());
}

void Mesh::buildBasisCycles()
{
    std::vector<Cycle> basisCycles;
    buildContractibleCycles(basisCycles);
    buildNonContractibleCycles(basisCycles);
    
    // build A
    A.resize((int)basisCycles.size(), (int)edges.size());
    A.setZero();
    
    std::vector<Eigen::Triplet<double>> ATriplet;
    for (int i = 0; i < (int)basisCycles.size(); i++) {
        for (int j = 0; j < (int)basisCycles[i].size(); j++) {
        
            HalfEdgeIter he = basisCycles[i][j];
            int e = he->edge->index;
            double coefficient = sqrt((he->cotan() + he->flip->cotan()) * 0.5);
            if (coefficient != coefficient) coefficient = 0; // check for negative cotan weights
            
            if (he == he->edge->he) {
                ATriplet.push_back(Eigen::Triplet<double>(i, e, coefficient));
                
            } else {
                ATriplet.push_back(Eigen::Triplet<double>(i, e, -coefficient));
            }
        }
    }

    A.setFromTriplets(ATriplet.begin(), ATriplet.end());
    solver.compute(A);
}

double parallelTransport(double angle, HalfEdgeIter& he)
{
    Eigen::Vector3d x, y;
    Eigen::Vector3d v = he->flip->vertex->position - he->vertex->position;
    if (he != he->edge->he) v = -v;
    
    // subtract delta ij
    he->face->axis(x, y);
    angle -= atan2(v.dot(y), v.dot(x));
    
    // subtract delta ji
    he->flip->face->axis(x, y);
    angle += atan2(v.dot(y), v.dot(x));
    
    return angle;
}

void Mesh::setup()
{
    // build basis cycles
    buildBasisCycles();

    // compute angle defects
    K.resize(A.rows());
    size_t v = vertices.size();
    for (size_t i = 0; i < A.rows(); i++) {
        if (i < v) {
            K(i) = -vertices[i].angleDefect();
            
        } else {
            // compute holonomy of the discrete Levi-Civita connection
            double alpha = 0;
            for (size_t j = 0; j < generators[i-v].size(); j++) {
                alpha = parallelTransport(alpha, generators[i-v][j]);
            }
            
            while (alpha >=  M_PI) alpha -= 2 * M_PI;
            while (alpha <  -M_PI) alpha += 2 * M_PI;
            
            K(i) = -alpha;
        }
    }
}

void Mesh::solve(const std::vector<double>& generatorKs)
{
    // set singularities
    double sum = 0;
    size_t v = vertices.size();
    Eigen::VectorXd b(K.rows()); 
    for (size_t i = 0; i < A.rows(); i++) {
        if (i < v) {
            b(i) = K(i) + 2 * M_PI * vertices[i].k;
            sum += vertices[i].k;
        }
        else b(i) = K(i) + 2 * M_PI * generatorKs[i-v];
    }
    
    if (std::abs(eulerCharacteristic() - sum) > EPSILON) {
        std::cerr << "Error: Indices must add to Euler Characteristic" << std::endl;
        return;
    }
    
    // solve for adjustment angles
    Eigen::VectorXd x = solver.solve(b);
    
    // set adjustment angles
    for (EdgeIter e = edges.begin(); e != edges.end(); e++) {
        double coefficient = sqrt((e->he->cotan() + e->he->flip->cotan()) * 0.5);
        if (coefficient != coefficient) coefficient = 0; // check for negative cotan weights
        e->phi = x(e->index) * coefficient;
    }
    
    // construct direction fields
    FaceIter root = faces.begin() + (rand() % faces.size());
    root->visited = true;
    root->beta = M_PI;
    std::queue<FaceIter> queue;
    queue.push(root);
    
    while (!queue.empty()) {
        FaceIter f = queue.front();
        queue.pop();
        
        HalfEdgeIter he = f->he;
        do {
            FaceIter neighbor = he->flip->face;
            if (!neighbor->visited) {
                neighbor->beta = parallelTransport(f->beta, he->flip) - he->edge->phi;
                neighbor->visited = true;
                queue.push(neighbor);
            }
            
            he = he->next;
            
        } while (he != f->he);
    }
    
    for (FaceIter f = faces.begin(); f != faces.end(); f++) {
        while (f->beta >=  M_PI) f->beta -= 2 * M_PI;
        while (f->beta <  -M_PI) f->beta += 2 * M_PI;
        
        Eigen::Vector3d n = f->normal().normalized();
        Eigen::Vector3d x = (f->he->next->vertex->position - f->he->vertex->position).normalized();
        if (f->he != f->he->edge->he) n = -n;
        f->field = x*cos(f->beta) + n.cross(x)*sin(f->beta) + n*(n.dot(x))*(1-cos(f->beta));
    }
}

void Mesh::normalize()
{
    // compute center of mass
    Eigen::Vector3d cm = Eigen::Vector3d::Zero();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        cm += v->position;
    }
    cm /= (double)vertices.size();
    
    // translate to origin
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position -= cm;
    }
    
    // determine radius
    double rMax = 0;
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        rMax = std::max(rMax, v->position.norm());
    }
    
    // rescale to unit sphere
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position /= rMax;
    }
}
