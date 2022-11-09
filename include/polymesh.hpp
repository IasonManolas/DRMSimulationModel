#ifndef POLYMESH_HPP
#define POLYMESH_HPP
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/import.h>
#include <filesystem>
#include "mesh.hpp"
#include "utilities.hpp"
#include "vcg/complex/complex.h"

#ifdef POLYSCOPE_DEFINED
#include <polyscope/surface_mesh.h>
#endif
class PFace;
class PVertex;
class PEdge;

struct PUsedTypes : public vcg::UsedTypes<vcg::Use<PVertex>::AsVertexType,
                                          vcg::Use<PFace>::AsFaceType,
                                          vcg::Use<PEdge>::AsEdgeType> {};

class PVertex : public vcg::Vertex<PUsedTypes,
                                   vcg::vertex::Coord3d,
                                   vcg::vertex::Normal3d,
                                   vcg::vertex::Mark,
                                   vcg::vertex::Qualityf,
                                   vcg::vertex::BitFlags,
                                   vcg::vertex::VFAdj,
                                   vcg::vertex::VEAdj> {};
class PEdge : public vcg::Edge<PUsedTypes,
                               vcg::edge::VertexRef,
                               vcg::edge::BitFlags,
                               vcg::edge::EEAdj,
                               vcg::edge::EFAdj,
                               vcg::edge::VEAdj,
                               vcg::edge::EVAdj> {};

class PFace
    : public vcg::Face<
          PUsedTypes,
          vcg::face::PolyInfo  // this is necessary  if you use component in
                               // vcg/simplex/face/component_polygon.h
          // It says "this class is a polygon and the memory for its components
          // (e.g. pointer to its vertices will be allocated dynamically")
          ,
          //                               vcg::face::FHAdj,
          vcg::face::PVFAdj,
          vcg::face::PFEAdj,
          vcg::face::PFVAdj,
          vcg::face::PVFAdj,
          //                               vcg::face::PVFAdj,
          vcg::face::PFFAdj  // Pointer to edge-adjacent face (just like FFAdj )
          ,
          vcg::face::BitFlags  // bit flags
          ,
          vcg::face::Mark,
          vcg::face::Qualityf  // quality
          ,
          vcg::face::Normal3f  // normal
          ,
          vcg::face::Color4b> {};

class VCGPolyMesh : public vcg::tri::TriMesh<std::vector<PVertex>,
                                             std::vector<PFace>,
                                             std::vector<PEdge>>,
                    public Mesh {
 public:
  virtual bool load(const std::filesystem::path& meshFilePath) override {
    int mask;
    vcg::tri::io::Importer<VCGPolyMesh>::LoadMask(meshFilePath.string().c_str(),
                                                  mask);
    int error = vcg::tri::io::Importer<VCGPolyMesh>::Open(
        *this, meshFilePath.string().c_str(), mask);
    if (error != 0) {
      std::cerr << "Could not load polygonal mesh:" << meshFilePath
                << std::endl;
      return false;
    }

    vcg::tri::UpdateTopology<VCGPolyMesh>::VertexEdge(*this);
    vcg::tri::UpdateTopology<VCGPolyMesh>::VertexFace(*this);
    vcg::tri::UpdateTopology<VCGPolyMesh>::FaceFace(*this);
    vcg::tri::UpdateTopology<VCGPolyMesh>::AllocateEdge(*this);
    vcg::tri::UpdateTopology<VCGPolyMesh>::EdgeEdge(*this);
    vcg::tri::UpdateNormal<VCGPolyMesh>::PerVertexNormalized(*this);
    vcg::tri::Clean<VCGPolyMesh>::RemoveUnreferencedVertex(*this);

    // finally remove valence 1 vertices on the border
    //        vcg::PolygonalAlgorithm<PolyMeshType>::RemoveValence2Vertices(dual);
    return true;
  }
  //    //    vcg::tri::io::ImporterOBJ<VCGPolyMesh>::Open();
  //    //    unsigned int mask = 0;
  //    //    mask |= nanoply::NanoPlyWrapper<VCGPolyMesh>::IO_VERTCOORD;
  //    //    mask |= nanoply::NanoPlyWrapper<VCGPolyMesh>::IO_VERTNORMAL;
  //    //    mask |= nanoply::NanoPlyWrapper<VCGPolyMesh>::IO_VERTCOLOR;
  //    //    mask |= nanoply::NanoPlyWrapper<VCGPolyMesh>::IO_EDGEINDEX;
  //    //    mask |= nanoply::NanoPlyWrapper<VCGPolyMesh>::IO_FACEINDEX;
  //    //    if (nanoply::NanoPlyWrapper<VCGPolyMesh>::LoadModel(
  //    //            std::filesystem::absolute(filename).c_str(), *this, mask)
  //    !=
  //    //            0) {
  //    //      std::cout << "Could not load tri mesh" << std::endl;
  //    //    }
  //  }
  Eigen::MatrixX3d getPositions() const {
    Eigen::MatrixX3d vertices(VN(), 3);
    for (size_t vi = 0; vi < VN(); vi++) {
      VCGPolyMesh::CoordType vertexCoordinates = vert[vi].cP();
      vertices.row(vi) = vertexCoordinates.ToEigenVector<Eigen::Vector3d>();
    }
    return vertices;
  }

  std::vector<std::vector<int>> getFaces() const {
    std::vector<std::vector<int>> faces(FN());
    for (const VCGPolyMesh::FaceType& f : this->face) {
      const int fi = vcg::tri::Index<VCGPolyMesh>(*this, f);
      for (size_t vi = 0; vi < f.VN(); vi++) {
        const size_t viGlobal = vcg::tri::Index<VCGPolyMesh>(*this, f.cV(vi));
        faces[fi].push_back(viGlobal);
      }
    }

    return faces;
  }

  std::vector<int> getEdges() const {
    std::vector<int> edges(2 * EN(), -1);
    for (const VCGPolyMesh::EdgeType& e : this->edge) {
      const int ei = vcg::tri::Index<VCGPolyMesh>(*this, e);
      for (size_t vi = 0; vi < 2; vi++) {
        const size_t viGlobal = vcg::tri::Index<VCGPolyMesh>(*this, e.cV(vi));
        edges[2 * ei + vi] = viGlobal;
      }
    }

    return edges;
  }

  //  bool load(const std::filesystem::__cxx11::path &meshFilePath)
  //  {
  //      const std::string extension = ".ply";
  //      std::filesystem::path filePath = meshFilePath;
  //      assert(std::filesystem::path(filePath).extension().string() ==
  //      extension); unsigned int mask = 0; mask |=
  //      vcg::tri::io::Mask::IOM_VERTCOORD; mask |=
  //      vcg::tri::io::Mask::IOM_VERTNORMAL; mask |=
  //      vcg::tri::io::Mask::IOM_FACEINDEX; mask |=
  //      vcg::tri::io::Mask::IOM_FACECOLOR; if
  //      (vcg::tri::io::Importer<VCGPolyMesh>::Open(*this, filePath.c_str()) !=
  //      0) {
  //          return false;
  //      }
  //      label = meshFilePath.filename();
  //      return true;
  //  }

  bool save(const std::filesystem::path& meshFilePath =
                std::filesystem::path()) override {
    if (meshFilePath.extension() == ".obj") {
      return saveOBJ(meshFilePath);
    } else if (meshFilePath.extension() == ".ply") {
      return savePLY(meshFilePath);
    }

    return false;
  }

  bool saveOBJ(
      const std::filesystem::path& objFilePath = std::filesystem::path()) {
    const std::string extension = ".obj";
    std::filesystem::path filePath = objFilePath;
    if (filePath.empty()) {
      filePath = std::filesystem::current_path()
                     .append(getLabel() + extension)
                     .string();
    } else if (std::filesystem::is_directory(
                   std::filesystem::path(objFilePath))) {
      filePath = std::filesystem::path(objFilePath)
                     .append(getLabel() + extension)
                     .string();
    }
    assert(std::filesystem::path(filePath).extension().string() == extension);
    unsigned int mask = 0;
    mask |= vcg::tri::io::Mask::IOM_VERTCOORD;
    mask |= vcg::tri::io::Mask::IOM_VERTNORMAL;
    mask |= vcg::tri::io::Mask::IOM_FACEINDEX;
    mask |= vcg::tri::io::Mask::IOM_FACECOLOR;
    if (vcg::tri::io::ExporterOBJ<VCGPolyMesh>::Save(
            *this, filePath.string().c_str(), mask) != 0) {
      return false;
    }
    return true;
  }

  bool savePLY(
      const std::filesystem::path& objFilePath = std::filesystem::path()) {
    const std::string extension = ".ply";
    std::filesystem::path filePath = objFilePath;
    if (filePath.empty()) {
      filePath = std::filesystem::current_path()
                     .append(getLabel() + extension)
                     .string();
    } else if (std::filesystem::is_directory(
                   std::filesystem::path(objFilePath))) {
      filePath = std::filesystem::path(objFilePath)
                     .append(getLabel() + extension)
                     .string();
    }
    assert(std::filesystem::path(filePath).extension().string() == extension);
    unsigned int mask = 0;
    mask |= vcg::tri::io::Mask::IOM_VERTCOORD;
    mask |= vcg::tri::io::Mask::IOM_VERTNORMAL;
    mask |= vcg::tri::io::Mask::IOM_FACEINDEX;
    mask |= vcg::tri::io::Mask::IOM_FACECOLOR;
    if (vcg::tri::io::ExporterPLY<VCGPolyMesh>::Save(
            *this, filePath.string().c_str(), mask, false) != 0) {
      return false;
    }
    return true;
  }

#ifdef POLYSCOPE_DEFINED
  using RGBColor = std::array<float, 3>;
  polyscope::SurfaceMesh* registerForDrawing(
      const std::optional<RGBColor>& desiredColor = std::nullopt,
      const bool& shouldEnable = true) {
    auto vertices = getPositions();
    auto faces = getFaces();
    PolyscopeInterface::init();

    polyscope::SurfaceMesh* polyscopeHandle_mesh =
        polyscope::registerSurfaceMesh(label, vertices, faces);

    const double drawingRadius = 0.002;
    polyscopeHandle_mesh->setEnabled(shouldEnable);
    polyscopeHandle_mesh->setEdgeWidth(drawingRadius);

    if (desiredColor.has_value()) {
      const glm::vec3 desiredColor_glm(desiredColor.value()[0],
                                       desiredColor.value()[1],
                                       desiredColor.value()[2]);
      polyscopeHandle_mesh->setSurfaceColor(desiredColor_glm);
    }
    std::vector<std::array<double, 3>> colors(VN(),
                                              std::array<double, 3>({0, 0, 0}));
    for (int vi = 0; vi < VN(); vi++) {
      if (vert[vi].IsB()) {
        colors[vi] = std::array<double, 3>({1, 0, 0});
      }
    }
    polyscopeHandle_mesh->addVertexColorQuantity("Selected Vertices", colors)
        ->setEnabled(true);

    return polyscopeHandle_mesh;
  }
#endif
  void moveToCenter() {
    CoordType centerOfMass(0, 0, 0);

    for (int vi = 0; vi < VN(); vi++) {
      centerOfMass += vert[vi].cP();
    }
    centerOfMass /= VN();
    vcg::tri::UpdatePosition<VCGPolyMesh>::Translate(*this, -centerOfMass);
  }

  /*
   * Returns the average distance from the center of each edge to the center of
   * its face over the whole mesh
   * */
  double getAverageFaceRadius() const {
    double averageFaceRadius = 0;
    for (int fi = 0; fi < FN(); fi++) {
      const VCGPolyMesh::FaceType& f = face[fi];
      CoordType centerOfFace(0, 0, 0);
      for (int vi = 0; vi < f.VN(); vi++) {
        centerOfFace = centerOfFace + f.cP(vi);
      }
      centerOfFace /= f.VN();

      double faceRadius = 0;
      //                    for (int face_ei = 0; face_ei < f.EN(); face_ei++) {
      //          std::cout << "fi:" << getIndex(f) << std::endl;
      //          auto vps = f.FVp(0);
      //          auto vpe = vps;
      for (int i = 0; i < f.VN(); i++) {
        faceRadius += vcg::Distance(centerOfFace, (f.cP0(i) + f.cP1(i)) / 2);
      }

      //          }
      const int faceEdges = f.VN();  // NOTE: When does this not hold?
      faceRadius /= faceEdges;
      averageFaceRadius += faceRadius;
    }
    averageFaceRadius /= FN();
    return averageFaceRadius;
  }

  bool copy(VCGPolyMesh& copyFrom) {
    vcg::tri::Append<VCGPolyMesh, VCGPolyMesh>::MeshCopy(*this, copyFrom);
    label = copyFrom.getLabel();
    //  eigenEdges = mesh.getEigenEdges();
    //  if (eigenEdges.rows() == 0) {
    //    getEdges(eigenEdges);
    //  }
    //  eigenVertices = mesh.getEigenVertices();
    //  if (eigenVertices.rows() == 0) {
    //    getVertices();
    //  }
    vcg::tri::UpdateTopology<VCGPolyMesh>::VertexEdge(*this);
    vcg::tri::UpdateTopology<VCGPolyMesh>::VertexFace(*this);
    vcg::tri::UpdateTopology<VCGPolyMesh>::FaceFace(*this);
    vcg::tri::UpdateTopology<VCGPolyMesh>::AllocateEdge(*this);
    vcg::tri::UpdateTopology<VCGPolyMesh>::EdgeEdge(*this);
    //      vcg::tri::UpdateTopology<VCGPolyMesh>::VertexFace(*this);

    return true;
  }

  VCGPolyMesh(VCGPolyMesh& copyFrom) { copy(copyFrom); }
  VCGPolyMesh() {}
  template <typename MeshElement>
  size_t getIndex(const MeshElement& meshElement) const {
    return vcg::tri::Index<VCGPolyMesh>(*this, meshElement);
  }
};

using ConstVCGPolyMesh = VCGPolyMesh;

#endif  // POLYMESH_HPP
