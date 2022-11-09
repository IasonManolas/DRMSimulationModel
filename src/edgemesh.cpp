#include "edgemesh.hpp"
#include <wrap/io_trimesh/import.h>
#include "vcg/simplex/face/topology.h"
//#include <wrap/nanoply/include/nanoplyWrapper.hpp>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/import.h>
#include "polymesh.hpp"

Eigen::MatrixX2i VCGEdgeMesh::getEigenEdges() const {
  return eigenEdges;
}

std::vector<vcg::Point2i> VCGEdgeMesh::computeEdges() {
  computeEdges(eigenEdges);
  std::vector<vcg::Point2i> edges(eigenEdges.rows());
  for (int ei = 0; ei < eigenEdges.rows(); ei++) {
    edges[ei] = vcg::Point2i(eigenEdges(ei, 0), eigenEdges(ei, 1));
  }
  return edges;
}

Eigen::MatrixX3d VCGEdgeMesh::getEigenVertices() const {
  //    getVertices(eigenVertices);
  return eigenVertices;
}

Eigen::MatrixX3d VCGEdgeMesh::getEigenEdgeNormals() const {
  return eigenEdgeNormals;
}

bool VCGEdgeMesh::save(const std::filesystem::path& meshFilePath) {
  const std::string outputFileExtension = ".ply";

  std::string filename = meshFilePath.string();
  if (filename.empty()) {
    filename = std::filesystem::current_path()
                   .append(getLabel() + outputFileExtension)
                   .string();
  } else if (std::filesystem::is_directory(
                 std::filesystem::path(meshFilePath))) {
    filename = std::filesystem::path(meshFilePath)
                   .append(getLabel() + outputFileExtension)
                   .string();
  }
  //  assert(std::filesystem::path(filename).extension().string() == ".ply");
  unsigned int mask = 0;
  mask |= vcg::tri::io::Mask::IOM_VERTCOORD;
  mask |= vcg::tri::io::Mask::IOM_EDGEINDEX;
  mask |= vcg::tri::io::Mask::IOM_VERTNORMAL;
  mask |= vcg::tri::io::Mask::IOM_VERTCOLOR;
  //    if (nanoply::NanoPlyWrapper<VCGEdgeMesh>::SaveModel(filename.c_str(),
  //    *this, mask, false) != 0) {
  if (std::filesystem::is_directory(meshFilePath.parent_path())) {
    std::filesystem::create_directories(meshFilePath.parent_path());
  }
  if (vcg::tri::io::Exporter<VCGEdgeMesh>::Save(*this, filename.c_str(),
                                                mask) != 0) {
    return false;
  }
  return true;
}

void VCGEdgeMesh::GeneratedRegularSquaredPattern(
    const double angleDeg,
    std::vector<std::vector<vcg::Point2d>>& pattern,
    const size_t& desiredNumberOfSamples) {
  static const size_t piSamples = 10;

  // generate a pattern in a 1x1 quad
  const vcg::Point2d offset(0, 0);

  const size_t samplesNo = desiredNumberOfSamples;
  //      std::max(desiredNumberOfSamples, size_t(piSamples * (angleDeg /
  //      180)));
  const double angle = vcg::math::ToRad(angleDeg);

  pattern.clear();

  // first arm
  std::vector<vcg::Point2d> arm;
  {
    for (int k = 0; k <= samplesNo; k++) {
      const double t = double(k) / samplesNo;
      const double a = (1 - t) * angle;
      // const double r = vcg::math::Sin(t*M_PI_2) /*(1-((1-t)*(1-t)))*/ *
      // 0.5;
      const double r = t * 0.5;  // linear

      vcg::Point2d p(vcg::math::Cos(a), vcg::math::Sin(a));

      arm.push_back((p * r));
    }
    pattern.push_back(arm);
  }

  // other arms
  for (int i = 0; i < 3; i++) {
    for (vcg::Point2d& p : arm) {
      p = vcg::Point2d(-p.Y(), p.X());
    }
    pattern.push_back(arm);
  }

  assert(pattern.size() == 4);

  // offset all
  for (auto& arm : pattern) {
    for (vcg::Point2d& p : arm) {
      p += offset;
    }
  }
}

void VCGEdgeMesh::createSpiral(const float& degreesOfArm,
                               const size_t& numberOfSamples) {
  std::vector<std::vector<vcg::Point2d>> spiralPoints;
  GeneratedRegularSquaredPattern(degreesOfArm, spiralPoints, numberOfSamples);

  for (size_t armIndex = 0; armIndex < spiralPoints.size(); armIndex++) {
    for (size_t pointIndex = 0; pointIndex < spiralPoints[armIndex].size() - 1;
         pointIndex++) {
      const vcg::Point2d p0 = spiralPoints[armIndex][pointIndex];
      const vcg::Point2d p1 = spiralPoints[armIndex][pointIndex + 1];
      CoordType n(0, 0, 1);
      auto ei = vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(
          *this, VCGEdgeMesh::CoordType(p0.X(), p0.Y(), 0),
          VCGEdgeMesh::CoordType(p1.X(), p1.Y(), 0));
      ei->cV(0)->N() = n;
      ei->cV(1)->N() = n;
    }
  }

  //  setDefaultAttributes();
}

bool VCGEdgeMesh::createSpanGrid(const size_t squareGridDimension) {
  return createSpanGrid(squareGridDimension, squareGridDimension);
}

bool VCGEdgeMesh::createSpanGrid(const size_t desiredWidth,
                                 const size_t desiredHeight) {
  std::cout << "Grid of dimensions:" << desiredWidth << "," << desiredHeight
            << std::endl;
  const VCGEdgeMesh::CoordType n(0, 0, 1);
  int x = 0;
  int y = 0;
  //  for (size_t vi = 0; vi < numberOfVertices; vi++) {
  while (y <= desiredHeight) {
    //    std::cout << x << " " << y << std::endl;
    auto p = VCGEdgeMesh::CoordType(x, y, 0);
    vcg::tri::Allocator<VCGEdgeMesh>::AddVertex(*this, p, n);
    x++;
    if (x > desiredWidth) {
      x = 0;
      y++;
    }
  }

  for (size_t vi = 0; vi < VN(); vi++) {
    int x = vi % (desiredWidth + 1);
    int y = vi / (desiredWidth + 1);
    const bool isCornerNode = (y == 0 && x == 0) ||
                              (y == 0 && x == desiredWidth) ||
                              (y == desiredHeight && x == 0) ||
                              (y == desiredHeight && x == desiredWidth);
    if (isCornerNode) {
      continue;
    }
    if (y == 0) {  // row 0.Connect with node above
      vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(*this, vi,
                                                vi + desiredWidth + 1);
      continue;
    } else if (x == 0) {  // col 0.Connect with node to the right
      vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(*this, vi, vi + 1);
      continue;
    } else if (y == desiredHeight) {  // row 0.Connect with node below
      //      vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(*this, vi,
      //                                                vi - (desiredWidth +
      //                                                1));
      continue;
    } else if (x == desiredWidth) {  // row 0.Connect with node to the left
      //      vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(*this, vi, vi - 1);
      continue;
    }

    vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(*this, vi, vi + desiredWidth + 1);
    vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(*this, vi, vi + 1);
    //    vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(*this, vi,
    //                                              vi - (desiredWidth + 1));
    //    vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(*this, vi, vi - 1);
  }

  vcg::tri::Allocator<VCGEdgeMesh>::DeleteVertex(*this, vert[0]);
  vcg::tri::Allocator<VCGEdgeMesh>::DeleteVertex(*this, vert[desiredWidth]);
  vcg::tri::Allocator<VCGEdgeMesh>::DeleteVertex(
      *this, vert[desiredHeight * (desiredWidth + 1)]);
  vcg::tri::Allocator<VCGEdgeMesh>::DeleteVertex(
      *this, vert[(desiredHeight + 1) * (desiredWidth + 1) - 1]);
  vcg::tri::Allocator<VCGEdgeMesh>::CompactVertexVector(*this);
  computeEdges(eigenEdges);
  computeVertices(eigenVertices);
  //  vcg::tri::Allocator<VCGEdgeMesh>::CompactEdgeVector(*this);

  //  const size_t numberOfEdges =
  //      desiredHeight * (desiredWidth - 1) + desiredWidth * (desiredHeight -
  //      1);
  //  handleBeamDimensions._handle->data.resize(
  //      numberOfEdges, CylindricalElementDimensions(0.03, 0.026));
  //  handleBeamMaterial._handle->data.resize(numberOfEdges,
  //                                          ElementMaterial(0.3, 200));

  return true;
}

bool VCGEdgeMesh::load(const std::filesystem::path& meshFilePath) {
  std::string usedPath = meshFilePath.string();
  if (std::filesystem::path(meshFilePath).is_relative()) {
    usedPath = std::filesystem::absolute(meshFilePath).string();
  }
  assert(std::filesystem::exists(usedPath));
  Clear();
  //  if (!loadUsingNanoply(usedPath)) {
  //    std::cerr << "Error: Unable to open " + usedPath << std::endl;
  //    return false;
  //  }
  if (!loadUsingDefaultLoader(usedPath)) {
    std::cerr << "Error: Unable to open " + usedPath << std::endl;
    return false;
  }

  computeEdges(eigenEdges);
  computeVertices(eigenVertices);
  vcg::tri::UpdateTopology<VCGEdgeMesh>::VertexEdge(*this);
  label = std::filesystem::path(meshFilePath).stem().string();

  const bool printDebugInfo = false;
  if (printDebugInfo) {
    std::cout << meshFilePath << " was loaded successfuly." << std::endl;
    std::cout << "Mesh has " << EN() << " edges." << std::endl;
  }

  label = std::filesystem::path(meshFilePath).stem().string();
  return true;
}

// bool VCGEdgeMesh::loadUsingNanoply(const std::string &plyFilename) {

//  this->Clear();
//  //  assert(plyFileHasAllRequiredFields(plyFilename));
//  // Load the ply file
//  unsigned int mask = 0;
//  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_VERTCOORD;
//  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_VERTNORMAL;
//  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_VERTCOLOR;
//  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_EDGEINDEX;
//  if (nanoply::NanoPlyWrapper<VCGEdgeMesh>::LoadModel(plyFilename.c_str(),
//                                                      *this, mask) != 0) {
//    return false;
//  }
//  return true;
//}

// bool VCGEdgeMesh::loadUsingDefaultLoader(const std::string& filePath) {
//  Clear();
//  //  assert(plyFileHasAllRequiredFields(plyFilename));
//  // Load the ply file
//  int mask = 0;
//  //  if (nanoply::NanoPlyWrapper<VCGEdgeMesh>::LoadModel(plyFilename.c_str(),
//  //                                                      *this, mask) != 0) {
//  const int loadErrorCode =
//      vcg::tri::io::Importer<VCGEdgeMesh>::Open(*this, filePath.c_str(),
//      mask);
//  if (mask & Mask::IOM_BITPOLYGONAL) {
//    VCGPolyMesh poly;
//    poly.load(filePath);
//  }
//  if (loadErrorCode != 0) {
//    std::cerr << vcg::tri::io::Importer<VCGEdgeMesh>::ErrorMsg(loadErrorCode)
//              << std::endl;
//    return false;
//  }
//  return true;
//}

bool VCGEdgeMesh::loadUsingDefaultLoader(
    const std::filesystem::path& meshFilePath) {
  Clear();
  // Load the ply file
  int mask = 0;
  mask |= vcg::tri::io::Mask::IOM_VERTCOORD;
  mask |= vcg::tri::io::Mask::IOM_VERTNORMAL;
  mask |= vcg::tri::io::Mask::IOM_VERTCOLOR;
  mask |= vcg::tri::io::Mask::IOM_EDGEINDEX;
  mask |= vcg::tri::io::Mask::IOM_VERTQUALITY;
  //  if (nanoply::NanoPlyWrapper<VCGEdgeMesh>::L
  //  mask |= /*vcg::tri::io::Mask::IOM_BITPOLYGONAL;*/
  //      vcg::tri::io::Importer<VCGPolyMesh>::LoadMask(meshFilePath.c_str(),
  //      mask);
  constexpr bool shouldComputeBorderVertices = true;
  const int loadErrorCode = vcg::tri::io::Importer<VCGEdgeMesh>::Open(
      *this, meshFilePath.string().c_str(), mask);
  if (mask & vcg::tri::io::Mask::IOM_BITPOLYGONAL) {
    VCGPolyMesh polyMesh;
    polyMesh.load(meshFilePath);
    vcg::tri::UpdateNormal<VCGPolyMesh>::PerVertexNormalized(polyMesh);
    //    std::cout << "vertices before:" << polyMesh.VN() << std::endl;
    int mergedVertices =
        vcg::tri::Clean<VCGPolyMesh>::MergeCloseVertex(polyMesh, 0.000061524);
    //    std::cout << "Merged vertices:" << mergedVertices << std::endl;
    vcg::tri::Allocator<VCGPolyMesh>::CompactEveryVector(polyMesh);
    vcg::tri::UpdateTopology<VCGPolyMesh>::VertexEdge(polyMesh);
    vcg::tri::UpdateTopology<VCGPolyMesh>::VertexFace(polyMesh);
    vcg::tri::UpdateTopology<VCGPolyMesh>::EdgeEdge(polyMesh);
    vcg::tri::UpdateTopology<VCGPolyMesh>::FaceFace(polyMesh);
    //    std::cout << "vertices after:" << polyMesh.VN() << std::endl;
    //    std::vector<VCGPolyMesh::FaceType*> selfIntersectingFaces;
    //    vcg::tri::Clean<VCGPolyMesh>::SelfIntersections(polyMesh,
    //                                                    selfIntersectingFaces);
    //    if (!selfIntersectingFaces.empty()) {
    //      std::cerr << "polygonal mesh contains self intersecting
    //      faces"
    //                << std::endl;
    //      std::for_each(selfIntersectingFaces.begin(),
    //      selfIntersectingFaces.end(),
    //                    [](VCGPolyMesh::FaceType* pf) { pf->SetD(); });
    //    }
    //    vcg::tri::Allocator<VCGPolyMesh>::CompactVertexVector(polyMesh);
    //    vcg::tri::Allocator<VCGPolyMesh>::CompactEdgeVector(polyMesh);
    //    vcg::tri::Allocator<VCGPolyMesh>::CompactFaceVector(polyMesh);
    //    vcg::tri::UpdateTopology<VCGPolyMesh>::VertexEdge(polyMesh);
    //    vcg::tri::UpdateTopology<VCGPolyMesh>::VertexFace(polyMesh);
    //    vcg::tri::UpdateTopology<VCGPolyMesh>::EdgeEdge(polyMesh);
    //    vcg::tri::UpdateTopology<VCGPolyMesh>::FaceFace(polyMesh);
    vcg::tri::UpdateFlags<VCGPolyMesh>::VertexBorderFromFaceAdj(polyMesh);
    //    polyMesh.registerForDrawing();
    // TODO: create a getPositions return std::vector<double>
    const std::vector<double> flattenedPositions = [&]() {
      Eigen::MatrixX3d positions = polyMesh.getPositions();
      //      std::cout << "positions size:" << positions.size() / 3 <<
      //      std::endl;
      Eigen::MatrixXd At = positions.transpose();
      Eigen::VectorXd v =
          Eigen::Map<const Eigen::VectorXd>(At.data(), At.size());
      std::vector<double> vec(v.data(), v.data() + v.size());
      return vec;
    }();
    const auto edges = polyMesh.getEdges();
    //    set(flattenedPositions, edges);
    //    std::cout << "edge mesh verts after set:" << VN() << std::endl;
    vcg::tri::Append<VCGEdgeMesh, VCGPolyMesh>::MeshCopyConst(*this, polyMesh);
    //    /*int merged =*/removeDuplicateVertices(/*pu_vertices,
    //        pu_edges*/);
    //    std::cout << "edge mesh verts after remove dups:" << VN() <<
    //    std::endl;
    assert(polyMesh.VN() == VN());
    //    for (int vi = 0; vi < polyMesh.VN(); vi++) {
    //      if (polyMesh.vert[vi].IsB())
    //        vert[vi].SetB();
    //    }
    //    registerForDrawing();
    //    polyscope::show();
    //    vcg::tri::Allocator<VCGEdgeMesh>::PointerUpdater<VertexPointer>
    //    pu_vertices;
    //    vcg::tri::Allocator<VCGEdgeMesh>::PointerUpdater<EdgePointer>
    //    pu_edges;
    //    for (int vi = 0; vi < polyMesh.VN(); vi++) {
    //      if (pu_vertices[vi].IsB())
    //      std::cout<<
    //    }
    //    save("./this.ply");
    //  } else if (FN() != 0) {
    //    vcg::tri::UpdateNormal<VCGEdgeMesh>::PerVertexNormalized(*this);
    //    //
    //    vcg::tri::UpdateNormal<VCGEdgeMesh>::NormalizePerVertex(*this);

  } else if (shouldComputeBorderVertices) {
    removeDuplicateVertices();
    vcg::tri::UpdateTopology<VCGEdgeMesh>::VertexEdge(*this);
    //    vcg::tri::UpdateTopology<VCGEdgeMesh>::VertexFace(*this);
    vcg::tri::UpdateTopology<VCGEdgeMesh>::EdgeEdge(*this);
    vcg::tri::UpdateTopology<VCGEdgeMesh>::FaceFace(*this);
    if (VN() != 0) {
      vcg::tri::UpdateFlags<VCGEdgeMesh>::VertexBorderFromFaceAdj(*this);

    } else {
      vcg::tri::UpdateFlags<VCGEdgeMesh>::VertexBorderFromEdgeAdj(*this);
    }
  }
  vcg::tri::UpdateTopology<VCGEdgeMesh>::VertexEdge(*this);
  //    vcg::tri::UpdateTopology<VCGEdgeMesh>::VertexFace(*this);
  vcg::tri::UpdateTopology<VCGEdgeMesh>::EdgeEdge(*this);
  if (EN() == 0) {
    vcg::tri::UpdateTopology<VCGEdgeMesh>::AllocateEdge(*this);
  }

  std::for_each(vert.begin(), vert.end(),
                [](VertexType& v) { v.N() = VectorType(0, 0, 1); });
  if (loadErrorCode != 0) {
    std::cerr << vcg::tri::io::Importer<VCGEdgeMesh>::ErrorMsg(loadErrorCode)
              << std::endl;
    return false;
  }
  // clear faces such that they are not exported when the mesh is saved
  //  this->face.clear();
  //  for (FaceIterator fi = face.begin(); fi != face.end(); ++fi)
  //    (*fi).Dealloc();
  //  fn = 0;

  assert(EN() != 0);
  return true;
}
// bool VCGEdgeMesh::plyFileHasAllRequiredFields(const std::string
// &plyFilename)
// {
//  const nanoply::Info info(plyFilename);
//  const std::vector<nanoply::PlyElement>::const_iterator edgeElemIt =
//      std::find_if(info.elemVec.begin(), info.elemVec.end(),
//                   [&](const nanoply::PlyElement &plyElem) {
//                     return plyElem.plyElem == nanoply::NNP_EDGE_ELEM;
//                   });
//  if (edgeElemIt == info.elemVec.end()) {
//    std::cerr << "Ply file is missing edge elements." << std::endl;
//    return false;
//  }

//  const std::vector<nanoply::PlyProperty> &edgePropertyVector =
//      edgeElemIt->propVec;
//  return hasPlyEdgeProperty(plyFilename, edgePropertyVector,
//                            plyPropertyBeamDimensionsID) &&
//         hasPlyEdgeProperty(plyFilename, edgePropertyVector,
//                            plyPropertyBeamMaterialID);
//}

Eigen::MatrixX3d VCGEdgeMesh::getNormals() const {
  Eigen::MatrixX3d vertexNormals;
  vertexNormals.resize(VN(), 3);

  for (int vertexIndex = 0; vertexIndex < VN(); vertexIndex++) {
    VCGEdgeMesh::CoordType vertexNormal =
        vert[static_cast<size_t>(vertexIndex)].cN();
    vertexNormals.row(vertexIndex) =
        vertexNormal.ToEigenVector<Eigen::Vector3d>();
  }

  return vertexNormals;
}
void VCGEdgeMesh::computeEdges(Eigen::MatrixX3d& edgeStartingPoints,
                               Eigen::MatrixX3d& edgeEndingPoints) const {
  edgeStartingPoints.resize(EN(), 3);
  edgeEndingPoints.resize(EN(), 3);
  for (int edgeIndex = 0; edgeIndex < EN(); edgeIndex++) {
    const VCGEdgeMesh::EdgeType& edge = this->edge[edgeIndex];
    edgeStartingPoints.row(edgeIndex) =
        edge.cP(0).ToEigenVector<Eigen::Vector3d>();
    edgeEndingPoints.row(edgeIndex) =
        edge.cP(1).ToEigenVector<Eigen::Vector3d>();
  }
}

VCGEdgeMesh::VCGEdgeMesh() {}

void VCGEdgeMesh::updateEigenEdgeAndVertices() {
#ifdef POLYSCOPE_DEFINED
  computeEdges(eigenEdges);
  computeVertices(eigenVertices);
#endif
}

bool VCGEdgeMesh::copy(const VCGEdgeMesh& mesh) {
  vcg::tri::Append<VCGEdgeMesh, VCGEdgeMesh>::MeshCopyConst(*this, mesh);
  label = mesh.getLabel();
  eigenEdges = mesh.getEigenEdges();
  //  assert(eigenEdges.rows() != 0);
  //  if (eigenEdges.rows() == 0) {
  //    getEdges(eigenEdges);
  //  }
  eigenVertices = mesh.getEigenVertices();
  //  assert(eigenVertices.rows() != 0);
  //  if (eigenVertices.rows() == 0) {
  //    getVertices(eigenVertices);
  //  }
  vcg::tri::UpdateTopology<VCGEdgeMesh>::VertexEdge(*this);

  return true;
}

void VCGEdgeMesh::set(const std::vector<double>& vertexPositions,
                      const std::vector<int>& edges) {
  Clear();
  for (int ei = 0; ei < edges.size(); ei += 2) {
    const int vi0 = edges[ei];
    const int vi1 = edges[ei + 1];
    const int vi0_offset = 3 * vi0;
    const int vi1_offset = 3 * vi1;
    const CoordType p0(vertexPositions[vi0_offset],
                       vertexPositions[vi0_offset + 1],
                       vertexPositions[vi0_offset + 2]);
    const CoordType p1(vertexPositions[vi1_offset],
                       vertexPositions[vi1_offset + 1],
                       vertexPositions[vi1_offset + 2]);
    auto eIt = vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(*this, p0, p1);
    CoordType n(0, 0, 1);
    eIt->cV(0)->N() = n;
    eIt->cV(1)->N() = n;
  }
  //    removeDuplicateVertices();

  updateEigenEdgeAndVertices();
}

void VCGEdgeMesh::removeDuplicateVertices() {
  vcg::tri::Clean<VCGEdgeMesh>::MergeCloseVertex(*this, 0.000061524);
  vcg::tri::Allocator<VCGEdgeMesh>::CompactVertexVector(*this);
  vcg::tri::Allocator<VCGEdgeMesh>::CompactEdgeVector(*this);
  vcg::tri::UpdateTopology<VCGEdgeMesh>::VertexEdge(*this);
  vcg::tri::UpdateTopology<VCGEdgeMesh>::EdgeEdge(*this);
}

void VCGEdgeMesh::removeDuplicateVertices(
    vcg::tri::Allocator<VCGEdgeMesh>::PointerUpdater<VertexPointer>&
        pu_vertices,
    vcg::tri::Allocator<VCGEdgeMesh>::PointerUpdater<EdgePointer>& pu_edges) {
  vcg::tri::Clean<VCGEdgeMesh>::MergeCloseVertex(*this, 0.000061524);
  vcg::tri::Allocator<VCGEdgeMesh>::CompactVertexVector(*this, pu_vertices);
  vcg::tri::Allocator<VCGEdgeMesh>::CompactEdgeVector(*this, pu_edges);
  vcg::tri::UpdateTopology<VCGEdgeMesh>::VertexEdge(*this);
}

std::string VCGEdgeMesh::ToSvgText(const VCGEdgeMesh& edgeMesh) {
  // Measures in mm (?)
  constexpr double ImgPadding = 5;
  // constexpr double PatternSize = 200 * 19.230769231;
  // constexpr double PatternSize = 200 * 17.913461538;
  constexpr double PatternSize = 1000;
  constexpr double ImgSize = 2 * ImgPadding + PatternSize;

  constexpr double StrokeWidth = 5;
  VCGEdgeMesh edgeMeshCopy;
  edgeMeshCopy.copy(edgeMesh);

  vcg::tri::UpdateBounding<VCGEdgeMesh>::Box(edgeMeshCopy);
  const double maxDim = edgeMeshCopy.bbox.Dim()[edgeMeshCopy.bbox.MaxDim()];
  vcg::tri::UpdatePosition<VCGEdgeMesh>::Translate(
      edgeMeshCopy,
      VCGEdgeMesh::CoordType(maxDim / 2.0 - edgeMeshCopy.bbox.Center().X(),
                             maxDim / 2.0 - edgeMeshCopy.bbox.Center().Y(), 0));
  vcg::tri::UpdatePosition<VCGEdgeMesh>::Scale(edgeMeshCopy,
                                               PatternSize / maxDim);

  // debug
  //	vcg::tri::UpdateBounding<EdgeMesh>::Box(em);
  //	std::cout << "pattern size "
  //	          << em.bbox.min.X() << " "
  //	          << em.bbox.max.X() << " "
  //	          << em.bbox.min.Y() << " "
  //	          << em.bbox.max.Y() << " " << std::endl;

  std::stringstream ss;

  // svg header
  ss << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl
     << "<!DOCTYPE svg  PUBLIC '-//W3C//DTD SVG 1.1//EN'  "
        "'http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd'>"
     << std::endl;

  // size & background
  ss << "<svg height=\"" << ImgSize << "\" width=\"" << ImgSize
     << "\" xmlns=\"http://www.w3.org/2000/svg\" "
        "xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
     << std::endl;
  ss << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>" << std::endl;

  for (const auto& e : edgeMeshCopy.edge) {
    const auto& p0 = e.cP(0);
    const auto& p1 = e.cP(1);
    ss << "<line "
       << "x1=\"" << ImgPadding + p0.X() << "\" y1=\""
       << -ImgPadding - p0.Y() + ImgSize << "\" "
       << "x2=\"" << ImgPadding + p1.X() << "\" y2=\""
       << -ImgPadding - p1.Y() + ImgSize << "\" "
       << "style=\"stroke:rgb(67,160,232);stroke-width:" << StrokeWidth
       << "\" stroke-linecap=\"round\"/>" << std::endl;
  }

  ss << "</svg>" << std::endl;
  return ss.str();
}

void VCGEdgeMesh::writeToSvg(const std::filesystem::path& writeToPath) const {
  // retrieve filepath for saving svg
  const std::filesystem::path svgPath = [=]() {
    if (writeToPath.empty()) {
      return std::filesystem::current_path().append(getLabel()).concat(".svg");
    }
    return std::filesystem::absolute(writeToPath);
  }();
  // save to svg file
  std::cout << "saving to " << svgPath << std::endl;

  std::ofstream ofs(svgPath);
  if (!ofs.is_open()) {
    std::cout << "unable to save to " << svgPath << std::endl;
    assert(false);
    return;
  }

  ofs << ToSvgText(*this);
  ofs.close();
}

void VCGEdgeMesh::deleteDanglingVertices() {
  vcg::tri::Allocator<VCGEdgeMesh>::PointerUpdater<VertexPointer> pu;
  deleteDanglingVertices(pu);
}

void VCGEdgeMesh::deleteDanglingVertices(
    vcg::tri::Allocator<VCGEdgeMesh>::PointerUpdater<TriMesh::VertexPointer>&
        pu) {
  for (VertexType& v : vert) {
    std::vector<VCGEdgeMesh::EdgePointer> incidentElements;
    vcg::edge::VEStarVE(&v, incidentElements);
    if (incidentElements.size() == 0 && !v.IsD()) {
      vcg::tri::Allocator<VCGEdgeMesh>::DeleteVertex(*this, v);
    }
  }

  vcg::tri::Allocator<VCGEdgeMesh>::CompactVertexVector(*this, pu);
  vcg::tri::Allocator<VCGEdgeMesh>::CompactEdgeVector(*this);

  updateEigenEdgeAndVertices();
}

void VCGEdgeMesh::computeVertices(Eigen::MatrixX3d& vertices) {
  vertices = Eigen::MatrixX3d();
  vertices.resize(VN(), 3);
  for (int vi = 0; vi < VN(); vi++) {
    if (vert[vi].IsD()) {
      continue;
    }
    VCGEdgeMesh::CoordType vertexCoordinates =
        vert[static_cast<size_t>(vi)].cP();
    vertices.row(vi) = vertexCoordinates.ToEigenVector<Eigen::Vector3d>();
  }
}

void VCGEdgeMesh::computeEdges(Eigen::MatrixX2i& edges) {
  edges = Eigen::MatrixX2i();
  edges.resize(EN(), 2);
  for (int edgeIndex = 0; edgeIndex < EN(); edgeIndex++) {
    const VCGEdgeMesh::EdgeType& edge = this->edge[edgeIndex];
    assert(!edge.IsD());
    auto vp0 = edge.cV(0);
    auto vp1 = edge.cV(1);
    assert(vcg::tri::IsValidPointer(*this, vp0) &&
           vcg::tri::IsValidPointer(*this, vp1));
    const size_t vi0 = vcg::tri::Index<VCGEdgeMesh>(*this, vp0);
    const size_t vi1 = vcg::tri::Index<VCGEdgeMesh>(*this, vp1);
    assert(vi0 != -1 && vi1 != -1);
    edges.row(edgeIndex) = Eigen::Vector2i(vi0, vi1);
  }
}

void VCGEdgeMesh::printVertexCoordinates(const size_t& vi) const {
  std::cout << "vi:" << vi << " " << vert[vi].cP()[0] << " " << vert[vi].cP()[1]
            << " " << vert[vi].cP()[2] << std::endl;
}

#ifdef POLYSCOPE_DEFINED
void VCGEdgeMesh::markVertices(const std::vector<size_t>& vertsToMark) {
  if (vertsToMark.empty()) {
    return;
  }

  std::vector<std::array<double, 3>> nodeColors(VN(), {0, 0, 0});
  for (const size_t vi : vertsToMark) {
    nodeColors[vi] = {1, 0, 0};
  }

  polyscope::getCurveNetwork(getLabel())
      ->addNodeColorQuantity("Marked vertices" + getLabel(), nodeColors)
      ->setEnabled(true);
}

// TODO: make const getEigenVertices is not
polyscope::CurveNetwork* VCGEdgeMesh::registerForDrawing(
    const std::optional<std::array<float, 3>>& desiredColor,
    const float& desiredRadius,
    const bool& shouldEnable) {
  PolyscopeInterface::init();
  const double drawingRadius = desiredRadius;
  updateEigenEdgeAndVertices();
  polyscope::CurveNetwork* polyscopeHandle_edgeMesh =
      polyscope::registerCurveNetwork(label, getEigenVertices(),
                                      getEigenEdges());
  //    std::cout << "EDGES:" << polyscopeHandle_edgeMesh->nEdges() <<
  //    std::endl;
  assert(polyscopeHandle_edgeMesh->nEdges() == getEigenEdges().rows() &&
         polyscopeHandle_edgeMesh->nNodes() == getEigenVertices().rows());

  polyscopeHandle_edgeMesh->setEnabled(shouldEnable);
  polyscopeHandle_edgeMesh->setRadius(drawingRadius, false);
  if (desiredColor.has_value()) {
    const glm::vec3 desiredColor_glm(desiredColor.value()[0],
                                     desiredColor.value()[1],
                                     desiredColor.value()[2]);
    polyscopeHandle_edgeMesh->setColor(
        /*glm::normalize(*/ desiredColor_glm /*)*/);
  }

  return polyscopeHandle_edgeMesh;
}

void VCGEdgeMesh::unregister() const {
  if (!polyscope::hasCurveNetwork(label)) {
    std::cerr << "No curve network registered with a name: " << getLabel()
              << std::endl;
    std::cerr << "Nothing to remove." << std::endl;
    return;
  }
  polyscope::removeCurveNetwork(label);
}

void VCGEdgeMesh::drawInitialFrames(
    polyscope::CurveNetwork* polyscopeHandle_initialMesh) const {
  Eigen::MatrixX3d frameInitialX(VN(), 3);
  Eigen::MatrixX3d frameInitialY(VN(), 3);
  Eigen::MatrixX3d frameInitialZ(VN(), 3);
  for (int vi = 0; vi < VN(); vi++) {
    frameInitialX.row(vi) = Eigen::Vector3d(1, 0, 0);
    frameInitialY.row(vi) = Eigen::Vector3d(0, 1, 0);
    frameInitialZ.row(vi) = Eigen::Vector3d(0, 0, 1);
  }
  polyscopeHandle_initialMesh->addNodeVectorQuantity("FrameX", frameInitialX)
      ->setVectorColor(glm::vec3(1, 0, 0))
      ->setEnabled(true);
  polyscopeHandle_initialMesh->addNodeVectorQuantity("FrameY", frameInitialY)
      ->setVectorColor(glm::vec3(0, 1, 0))
      ->setEnabled(true);
  polyscopeHandle_initialMesh->addNodeVectorQuantity("FrameZ", frameInitialZ)
      ->setVectorColor(glm::vec3(0, 0, 1))
      ->setEnabled(true);
}
#endif

std::vector<VertexIndex> VCGEdgeMesh::getBFSNeighbours(
    const VertexIndex& startVi,
    const int& maxNumberOfNeighbors) {
  std::vector<VertexIndex> bfsNeighbors;
  bfsNeighbors.reserve(maxNumberOfNeighbors);
  vcg::tri::UpdateFlags<VCGEdgeMesh>::VertexClearV(*this);
  // Create a queue for BFS
  std::deque<VertexPointer> queue;

  // Mark the current node as visited and enqueue it
  vert[startVi].SetV();
  queue.push_back(&vert[startVi]);

  while (!queue.empty() && bfsNeighbors.size() < maxNumberOfNeighbors) {
    // Dequeue a vertex from queue and print it
    const VertexPointer& vpOnExpandedBoarder = queue.front();
    queue.pop_front();

    // Get all adjacent vertices of the dequeued
    // vertex s. If a adjacent has not been visited,
    // then mark it visited and enqueue it
    std::vector<VertexPointer> adjacentVp;
    vcg::edge::VVStarVE(vpOnExpandedBoarder, adjacentVp);
    for (const VertexPointer& neighbour : adjacentVp) {
      if (!neighbour->IsV()) {
        neighbour->SetV();
        queue.push_back(neighbour);
        const int neighbourVi = getIndex(neighbour);
        bfsNeighbors.push_back(neighbourVi);
        if (bfsNeighbors.size() == maxNumberOfNeighbors) {
          break;
        }
      }
    }
  }

  if (bfsNeighbors.size() < maxNumberOfNeighbors) {
    std::cout << "Number of neighbours found is not the desired one:"
              << bfsNeighbors.size() << std::endl;
  }
  return bfsNeighbors;
}

std::unordered_set<VertexIndex> VCGEdgeMesh::getVerticesOfSameBorder(
    const VertexIndex& vi) {
  std::unordered_set<VertexIndex> borderVertices;

  vcg::tri::RequireVEAdjacency(*this);
  vcg::tri::UpdateTopology<MeshType>::VertexEdge(*this);
  vcg::tri::UpdateSelection<VCGEdgeMesh>::VertexClear(*this);
  VertexIndex currentVi = vi;
  do {
    vert[currentVi].SetS();
    borderVertices.insert(currentVi);
    std::vector<VertexPointer> VVStarVec;
    vcg::edge::VVStarVE(&(vert[currentVi]), VVStarVec);

    for (const VertexPointer& vp : VVStarVec) {
      const bool nonVisitedSameBorderVertex = !vp->IsS() && vp->IsB();
      if (nonVisitedSameBorderVertex) {
        currentVi = getIndex(vp);
        break;  // found a vertex on the same border
      }
    }
    if (vert[currentVi].IsS()) {
      break;  // traversed the whole border
    }
  } while (true);
  return borderVertices;
}

void VCGEdgeMesh::sampleMesh(VCGEdgeMesh& sampledEdgeMesh) const {
  //  sampledEdgeMesh.copy(*this);
  constexpr int numSamples = 10;
  for (const EdgeType& e : edge) {
    const auto& v0 = e.cV(0);
    const auto& p0 = v0->cP();
    const auto& v1 = e.cV(1);
    const auto& p1 = v1->cP();

    //    std::cout << "p0:" << p0[0] << " " << p0[1] << " " << p0[2] <<
    //    std::endl; std::cout << "p1:" << p1[0] << " " << p1[1] << " " << p1[2]
    //    << std::endl;

    VectorType sampleDisplacementVector = (p1 - p0) / (numSamples + 1);
    for (int sampleIndex = 0; sampleIndex < numSamples + 2; sampleIndex++) {
      CoordType samplePos =
          static_cast<double>(sampleIndex) * sampleDisplacementVector + p0;
      CoordType previousSamplePos =
          sampleIndex != 0 ? samplePos - sampleDisplacementVector : p0;
      vcg::tri::Allocator<VCGEdgeMesh>::AddEdge(sampledEdgeMesh, samplePos,
                                                previousSamplePos);
      //      std::cout << "samplePos:" << samplePos[0] << " " << samplePos[1]
      //      << " "
      //                << samplePos[2] << std::endl;
    }
  }

  sampledEdgeMesh.removeDuplicateVertices();
  //  sampledEdgeMesh.registerForDrawing();
  //  polyscope::show();
  //  return sampledEdgeMesh;
}
