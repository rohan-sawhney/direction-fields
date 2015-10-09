#ifndef MESH_H
#define MESH_H

#include "Types.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"
#include <Eigen/SparseCore>

typedef std::vector<HalfEdgeIter> Cycle;

class Mesh {
public:
    // default constructor
    Mesh();
    
    // copy constructor
    Mesh(const Mesh& mesh);
        
    // read mesh from file
    bool read(const std::string& fileName);
    
    // write mesh to file
    bool write(const std::string& fileName) const;
    
    // returns euler characteristic
    int eulerCharacteristic() const;
    
    // solves for adjustment angles
    void solve(const std::vector<double>& generatorKs);
    
    // member variables
    std::vector<HalfEdge> halfEdges;
    std::vector<Vertex> vertices;
    std::vector<Eigen::Vector3d> uvs;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<HalfEdgeIter> boundaries;
    std::vector<Cycle> generators;
    
private:
    // center mesh about origin and rescale to unit radius
    void normalize();
    
    // builds contractible cycles
    void buildContractibleCycles(std::vector<Cycle>& basisCycles);
    
    // builds primal spanning tree for tree-cotree decomoposition
    void buildPrimalSpanningTree();
    
    // builds dual spanning tree for tree-cotree decomoposition
    void buildDualSpanningTree();
    
    // builds non contractible cycles
    void buildNonContractibleCycles(std::vector<Cycle>& basisCycles);
    
    // builds basis cycles
    void buildBasisCycles();
    
    // sets up A and K 
    void setup();

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd K;
};

#endif