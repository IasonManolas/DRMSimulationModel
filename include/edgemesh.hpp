#ifndef EDGEMESH_HPP
#define EDGEMESH_HPP
#include <vcg/complex/complex.h>
#include <array>
#include <optional>
#include <unordered_set>
#include <vector>
#include "beam.hpp"
#include "mesh.hpp"
#include "utilities.hpp"

#ifdef POLYSCOPE_DEFINED
#include <polyscope/curve_network.h>
#endif

using EdgeIndex = size_t;
using VertexIndex = size_t;

class VCGEdgeMeshVertexType;
class VCGEdgeMeshEdgeType;
class VCGEdgeMeshFaceType;

struct VCGEdgeMeshUsedTypes
    : public vcg::UsedTypes<vcg::Use<VCGEdgeMeshVertexType>::AsVertexType,
                            vcg::Use<VCGEdgeMeshEdgeType>::AsEdgeType,
                            vcg::Use<VCGEdgeMeshFaceType>::AsFaceType> {};

class VCGEdgeMeshVertexType : public vcg::Vertex<VCGEdgeMeshUsedTypes,
                                                 vcg::vertex::Coord3d,
                                                 vcg::vertex::Normal3d,
                                                 vcg::vertex::BitFlags,
                                                 vcg::vertex::Color4b,
                                                 vcg::vertex::VEAdj> {};
class VCGEdgeMeshEdgeType : public vcg::Edge<VCGEdgeMeshUsedTypes,
                                             vcg::edge::VertexRef,
                                             vcg::edge::BitFlags,
                                             vcg::edge::EEAdj,
                                             vcg::edge::VEAdj> {};

class VCGEdgeMeshFaceType : public vcg::Face<VCGEdgeMeshUsedTypes,
                                             vcg::face::FFAdj,
                                             vcg::face::VFAdj,
                                             vcg::face::VertexRef,
                                             vcg::face::BitFlags,
                                             vcg::face::Normal3d> {};

class VCGEdgeMesh : public vcg::tri::TriMesh<std::vector<VCGEdgeMeshVertexType>,
                                             std::vector<VCGEdgeMeshEdgeType>,
                                             std::vector<VCGEdgeMeshFaceType>>,
                    public Mesh {
 protected:
  Eigen::MatrixX2i eigenEdges;
  Eigen::MatrixX3d eigenVertices;
  Eigen::MatrixX3d eigenEdgeNormals;

  void computeEdges(Eigen::MatrixX2i& edges);
  void computeVertices(Eigen::MatrixX3d& vertices);

 public:
  VCGEdgeMesh();
  template <typename MeshElement>
  size_t getIndex(const MeshElement& meshElement) {
    return vcg::tri::Index<VCGEdgeMesh>(*this, meshElement);
  }
  void updateEigenEdgeAndVertices();
  /*
   * The copy function shold be a virtual function of the base interface Mesh
   * class.
   * https://stackoverflow.com/questions/2354210/can-a-class-member-function-template-be-virtual
   * use type erasure (?)
   * */
  bool copy(const VCGEdgeMesh& mesh);

  void set(const std::vector<double>& vertexPositions,
           const std::vector<int>& edges);

  void removeDuplicateVertices();
  virtual void deleteDanglingVertices();
  virtual void deleteDanglingVertices(
      vcg::tri::Allocator<VCGEdgeMesh>::PointerUpdater<VertexPointer>& pu);

  void computeEdges(Eigen::MatrixX3d& edgeStartingPoints,
                    Eigen::MatrixX3d& edgeEndingPoints) const;

  Eigen::MatrixX3d getNormals() const;

  bool plyFileHasAllRequiredFields(const std::string& plyFilename);

  //    bool loadUsingNanoply(const std::string &plyFilename);

  bool load(const std::filesystem::path& meshFilePath) override;
  bool save(const std::filesystem::path& meshFilePath =
                std::filesystem::path()) override;

  bool createSpanGrid(const size_t squareGridDimension);
  bool createSpanGrid(const size_t desiredWidth, const size_t desiredHeight);
  void createSpiral(const float& degreesOfArm, const size_t& numberOfSamples);

  Eigen::MatrixX2i getEigenEdges() const;
  std::vector<vcg::Point2i> computeEdges();

  Eigen::MatrixX3d getEigenVertices() const;
  Eigen::MatrixX3d getEigenEdgeNormals() const;
  void printVertexCoordinates(const size_t& vi) const;
#ifdef POLYSCOPE_DEFINED
  polyscope::CurveNetwork* registerForDrawing(
      const std::optional<std::array<float, 3>>& desiredColor = std::nullopt,
      const float& desiredRadius = 0.001,
      const bool& shouldEnable = true);
  void unregister() const;
  void drawInitialFrames(
      polyscope::CurveNetwork* polyscopeHandle_initialMesh) const;
#endif
  void removeDuplicateVertices(
      vcg::tri::Allocator<VCGEdgeMesh>::PointerUpdater<VertexPointer>&
          pu_vertices,
      vcg::tri::Allocator<VCGEdgeMesh>::PointerUpdater<EdgePointer>& pu_edges);

  void centerMesh() {
    CoordType centerOfMass(0, 0, 0);

    for (const auto& v : vert) {
      centerOfMass = centerOfMass + v.cP();
    }
    centerOfMass = centerOfMass / VN();
    for (auto& v : vert) {
      v.P() = v.cP() - centerOfMass;
    }
  }

  static std::string ToSvgText(const VCGEdgeMesh& edgeMesh);
  void writeToSvg(
      const std::filesystem::path& writeToPath = std::filesystem::path()) const;

  void markVertices(const std::vector<size_t>& vertsToMark);
  /*
   * Returns the first maxNumberOfNeighbors vertices of a breadth-first search
   * starting at vertex startVi.
   * NOTE: Should be const but due to the way vcg is written and since i modify
   * flags of the mesh during the bfs it isnt const.
   * */
  std::vector<VertexIndex> getBFSNeighbours(const VertexIndex& startVi,
                                            const int& maxNumberOfNeighbors);
  /*
   * Returns the unordered set of vertices belonging to the same border as vi.
   * Assumes that any two borders are allways more than 2-edges apart
   * Clears selection flags
   * */
  std::unordered_set<VertexIndex> getVerticesOfSameBorder(
      const VertexIndex& vi);

  void sampleMesh(VCGEdgeMesh& sampledEdgeMesh) const;

 private:
  void GeneratedRegularSquaredPattern(
      const double angleDeg,
      std::vector<std::vector<vcg::Point2d>>& pattern,
      const size_t& desiredNumberOfSamples);
  bool loadUsingDefaultLoader(const std::filesystem::path& meshFilePath);
};

using VectorType = VCGEdgeMesh::CoordType;
using CoordType = VCGEdgeMesh::CoordType;
using VertexPointer = VCGEdgeMesh::VertexPointer;
using EdgePointer = VCGEdgeMesh::EdgePointer;
using ConstVCGEdgeMesh = VCGEdgeMesh;

#endif  // EDGEMESH_HPP
