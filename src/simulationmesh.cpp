#include "simulationmesh.hpp"
#include <Eigen/Dense>

SimulationEdgeMesh::SimulationEdgeMesh() {
  elements =
      vcg::tri::Allocator<SimulationEdgeMesh>::GetPerEdgeAttribute<Element>(
          *this, std::string("Elements"));
  nodes = vcg::tri::Allocator<SimulationEdgeMesh>::GetPerVertexAttribute<Node>(
      *this, std::string("Nodes"));
}

SimulationEdgeMesh::SimulationEdgeMesh(VCGEdgeMesh& mesh) {
  vcg::tri::MeshAssert<VCGEdgeMesh>::VertexNormalNormalized(mesh);

  VCGEdgeMesh::copy(mesh);

  elements =
      vcg::tri::Allocator<SimulationEdgeMesh>::GetPerEdgeAttribute<Element>(
          *this, std::string("Elements"));
  nodes = vcg::tri::Allocator<SimulationEdgeMesh>::GetPerVertexAttribute<Node>(
      *this, std::string("Nodes"));
  initializeNodes();
  initializeElements();
}

SimulationEdgeMesh::~SimulationEdgeMesh() {
  vcg::tri::Allocator<SimulationEdgeMesh>::DeletePerEdgeAttribute<Element>(
      *this, elements);
  vcg::tri::Allocator<SimulationEdgeMesh>::DeletePerVertexAttribute<Node>(
      *this, nodes);
}

SimulationEdgeMesh::SimulationEdgeMesh(SimulationEdgeMesh& m) {
  vcg::tri::MeshAssert<SimulationEdgeMesh>::VertexNormalNormalized(m);
  VCGEdgeMesh::copy(m);

  elements =
      vcg::tri::Allocator<SimulationEdgeMesh>::GetPerEdgeAttribute<Element>(
          *this, std::string("Elements"));
  nodes = vcg::tri::Allocator<SimulationEdgeMesh>::GetPerVertexAttribute<Node>(
      *this, std::string("Nodes"));
  initializeNodes();

  for (size_t ei = 0; ei < EN(); ei++) {
    elements[ei] = m.elements[ei];
  }
  reset();
}

void SimulationEdgeMesh::computeElementalProperties() {
  const std::vector<CrossSectionType> elementalDimensions = getBeamDimensions();
  const std::vector<ElementMaterial> elementalMaterials = getBeamMaterial();
  assert(EN() == elementalDimensions.size() &&
         elementalDimensions.size() == elementalMaterials.size());

  for (const EdgeType& e : edge) {
    const EdgeIndex ei = getIndex(e);
    elements[e].setDimensions(elementalDimensions[ei]);
    elements[e].setMaterial(elementalMaterials[ei]);
  }
}

void SimulationEdgeMesh::initializeNodes() {
  // set initial and previous locations,
  for (const VertexType& v : vert) {
    const VertexIndex vi = getIndex(v);
    Node& node = nodes[v];
    node.vi = vi;
    node.initialLocation = v.cP();
    node.initialNormal = v.cN();
    node.derivativeOfNormal.resize(6, VectorType(0, 0, 0));

    node.displacements[3] =
        v.cN()[0];  // initialize nx diplacement with vertex normal x
                    //        component
    node.displacements[4] =
        v.cN()[1];  // initialize ny(in the paper) diplacement with vertex
    //        normal
    // y component.

    // Initialize incident elements
    std::vector<VCGEdgeMesh::EdgePointer> incidentElements;
    vcg::edge::VEStarVE(&v, incidentElements);
    //    assert(
    //        vcg::tri::IsValidPointer<SimulationEdgeMesh>(*this,
    //        incidentElements[0]) && incidentElements.size() > 0);
    if (incidentElements.size() != 0) {
      nodes[v].incidentElements = incidentElements;
      node.referenceElement = getReferenceElement(v);
      //        std::vector<int>
      //        incidentElementsIndices(node.incidentElements.size()); if
      //        (drawGlobal && vi == 5) {
      //            std::vector<glm::vec3> edgeColors(EN(), glm::vec3(0, 1, 0));
      //            std::vector<glm::vec3> vertexColors(VN(), glm::vec3(0, 1,
      //            0)); vertexColors[vi] = glm::vec3(0, 0, 1); for (int iei =
      //            0; iei < incidentElementsIndices.size(); iei++) {
      //                incidentElementsIndices[iei] =
      //                this->getIndex(node.incidentElements[iei]);
      //                edgeColors[incidentElementsIndices[iei]] = glm::vec3(1,
      //                0, 0);
      //            }
      //            polyHandle->addEdgeColorQuantity("chosenE",
      //            edgeColors)->setEnabled(true);
      //            polyHandle->addNodeColorQuantity("chosenV",
      //            vertexColors)->setEnabled(true); draw();
      //        }
      //        const int referenceElementIndex =
      //        getIndex(node.referenceElement);
      // Initialze alpha angles

      const EdgeType& referenceElement = *node.referenceElement;
      const VectorType t01 =
          computeT1Vector(referenceElement.cP(0), referenceElement.cP(1));
      const VectorType f01 = (t01 - (v.cN() * (t01.dot(v.cN())))).Normalize();
      node.alphaAngles.reserve(incidentElements.size());

      for (const VCGEdgeMesh::EdgePointer& ep : nodes[v].incidentElements) {
        assert(referenceElement.cV(0) == ep->cV(0) ||
               referenceElement.cV(0) == ep->cV(1) ||
               referenceElement.cV(1) == ep->cV(0) ||
               referenceElement.cV(1) == ep->cV(1));
        const VectorType t1 = computeT1Vector(*ep);
        const VectorType f1 = t1 - (v.cN() * (t1.dot(v.cN()))).Normalize();
        const EdgeIndex ei = getIndex(ep);
        const double alphaAngle = computeAngle(f01, f1, v.cN());
        node.alphaAngles.emplace_back(std::make_pair(ei, alphaAngle));
      }
    }
  }
}

void SimulationEdgeMesh::reset() {
  for (const EdgeType& e : edge) {
    Element& element = elements[e];
    element.ei = getIndex(e);
    const VCGEdgeMesh::CoordType p0 = e.cP(0);
    const VCGEdgeMesh::CoordType p1 = e.cP(1);
    const vcg::Segment3<double> s(p0, p1);
    element.initialLength = s.Length();
    element.length = element.initialLength;
    element.updateRigidity();
  }

  for (const VertexType& v : vert) {
    Node& node = nodes[v];
    node.vi = getIndex(v);
    node.initialLocation = v.cP();
    node.initialNormal = v.cN();
    node.derivativeOfNormal.resize(6, VectorType(0, 0, 0));

    node.displacements[3] =
        v.cN()[0];  // initialize nx diplacement with vertex normal x
                    //        component
    node.displacements[4] =
        v.cN()[1];  // initialize ny(in the paper) diplacement with vertex
    //        normal
    // y component.

    const EdgeType& referenceElement = *getReferenceElement(v);
    const VectorType t01 =
        computeT1Vector(referenceElement.cP(0), referenceElement.cP(1));
    const VectorType f01 = (t01 - (v.cN() * (t01.dot(v.cN())))).Normalize();
    node.alphaAngles.clear();
    node.alphaAngles.reserve(node.incidentElements.size());
    for (const VCGEdgeMesh::EdgePointer& ep : nodes[v].incidentElements) {
      assert(referenceElement.cV(0) == ep->cV(0) ||
             referenceElement.cV(0) == ep->cV(1) ||
             referenceElement.cV(1) == ep->cV(0) ||
             referenceElement.cV(1) == ep->cV(1));
      const VectorType t1 = computeT1Vector(*ep);
      const VectorType f1 = t1 - (v.cN() * (t1.dot(v.cN()))).Normalize();
      const EdgeIndex ei = getIndex(ep);
      const double alphaAngle = computeAngle(f01, f1, v.cN());
      node.alphaAngles.emplace_back(std::make_pair(ei, alphaAngle));
    }
  }
}

void SimulationEdgeMesh::initializeElements() {
  for (const EdgeType& e : edge) {
    Element& element = elements[e];
    element.ei = getIndex(e);
    // Initialize dimensions
    element.dimensions = CrossSectionType();
    // Initialize material
    element.material = ElementMaterial();
    // Initialize lengths
    const VCGEdgeMesh::CoordType p0 = e.cP(0);
    const VCGEdgeMesh::CoordType p1 = e.cP(1);
    const vcg::Segment3<double> s(p0, p1);
    element.initialLength = s.Length();
    element.length = element.initialLength;
    // Initialize const factors
    element.updateRigidity();
    element.derivativeT1.resize(
        2, std::vector<VectorType>(6, VectorType(0, 0, 0)));
    element.derivativeT2.resize(
        2, std::vector<VectorType>(6, VectorType(0, 0, 0)));
    element.derivativeT3.resize(
        2, std::vector<VectorType>(6, VectorType(0, 0, 0)));
    element.derivativeR.resize(2,
                               std::vector<VectorType>(6, VectorType(0, 0, 0)));
  }
  updateElementalFrames();
}

void SimulationEdgeMesh::updateElementalLengths() {
  for (const EdgeType& e : edge) {
    const EdgeIndex ei = getIndex(e);
    const VertexIndex vi0 = getIndex(e.cV(0));
    const VCGEdgeMesh::CoordType p0 = e.cP(0);
    const VertexIndex vi1 = getIndex(e.cV(1));
    const VCGEdgeMesh::CoordType p1 = e.cP(1);
    const vcg::Segment3<double> s(p0, p1);
    const double elementLength = s.Length();
    elements[e].length = elementLength;
  }
}

void SimulationEdgeMesh::updateElementalFrames() {
  for (const EdgeType& e : edge) {
    const VectorType elementNormal =
        (e.cV(0)->cN() + e.cV(1)->cN()).Normalize();
    elements[e].frame = computeElementFrame(e.cP(0), e.cP(1), elementNormal);
  }
}

#ifdef POLYSCOPE_DEFINED
polyscope::CurveNetwork* SimulationEdgeMesh::registerForDrawing(
    const std::optional<std::array<float, 3>>& desiredColor,
    const bool& shouldEnable,
    const std::optional<float>& inputDrawingRadius) {
  //    std::cout << __FUNCTION__ << " revert this" << std::endl;
  const float drawingRadius = [&]() {
    if (inputDrawingRadius.has_value()) {
      return inputDrawingRadius.value();
    } else {
      return getBeamDimensions()[0].getDrawingRadius();
    }
  }();
  return VCGEdgeMesh::registerForDrawing(desiredColor, drawingRadius,
                                         shouldEnable);
}

void SimulationEdgeMesh::unregister() const {
  VCGEdgeMesh::unregister();
}
#endif

void SimulationEdgeMesh::setBeamCrossSection(
    const CrossSectionType& beamDimensions) {
  for (size_t ei = 0; ei < EN(); ei++) {
    elements[ei].dimensions = beamDimensions;
    elements[ei].updateRigidity();
  }
}

void SimulationEdgeMesh::setBeamMaterial(const double& pr, const double& ym) {
  for (size_t ei = 0; ei < EN(); ei++) {
    elements[ei].setMaterial(ElementMaterial(pr, ym));
    //    elements[ei].updateRigidity();
  }
}

void SimulationEdgeMesh::setBeamMaterial(const ElementMaterial& material) {
  for (size_t ei = 0; ei < EN(); ei++) {
    elements[ei].setMaterial(material);
    //    elements[ei].updateRigidity();
  }
}

std::vector<CrossSectionType> SimulationEdgeMesh::getBeamDimensions() {
  std::vector<CrossSectionType> beamDimensions(EN());
  for (size_t ei = 0; ei < EN(); ei++) {
    beamDimensions[ei] = elements[ei].dimensions;
  }
  return beamDimensions;
}

std::vector<ElementMaterial> SimulationEdgeMesh::getBeamMaterial() {
  std::vector<ElementMaterial> beamMaterial(EN());
  for (size_t ei = 0; ei < EN(); ei++) {
    beamMaterial[ei] = elements[ei].material;
  }
  return beamMaterial;
}

bool SimulationEdgeMesh::load(const std::filesystem::path& plyFilePath) {
  this->Clear();
  //  assert(plyFileHasAllRequiredFields(plyFilename));
  // Load the ply file
  //  VCGEdgeMesh::PerEdgeAttributeHandle<CrossSectionType> handleBeamDimensions
  //  =
  //      vcg::tri::Allocator<SimulationEdgeMesh>::AddPerEdgeAttribute<
  //          CrossSectionType>(*this, plyPropertyBeamDimensionsID);
  //  VCGEdgeMesh::PerEdgeAttributeHandle<ElementMaterial> handleBeamMaterial =
  //      vcg::tri::Allocator<SimulationEdgeMesh>::AddPerEdgeAttribute<ElementMaterial>(
  //          *this, plyPropertyBeamMaterialID);
  //  nanoply::NanoPlyWrapper<SimulationEdgeMesh>::CustomAttributeDescriptor
  //      customAttrib;
  //  customAttrib.GetMeshAttrib(plyFilename);
  //  customAttrib.AddEdgeAttribDescriptor<CrossSectionType, double, 2>(
  //      plyPropertyBeamDimensionsID, nanoply::NNP_LIST_INT8_FLOAT64, nullptr);
  //  /*FIXME: Since I allow CrossSectionType to take two types I should export
  //  the
  //   * type as well such that that when loaded the correct type of cross
  //   section
  //   * is used.
  //   */
  //  customAttrib.AddEdgeAttribDescriptor<vcg::Point2d, double, 2>(
  //      plyPropertyBeamMaterialID, nanoply::NNP_LIST_INT8_FLOAT64, nullptr);

  //  VCGEdgeMesh::PerEdgeAttributeHandle<std::array<double, 6>>
  //      handleBeamProperties =
  //          vcg::tri::Allocator<SimulationEdgeMesh>::AddPerEdgeAttribute<
  //              std::array<double, 6>>(*this, plyPropertyBeamProperties);
  //  customAttrib.AddEdgeAttribDescriptor<std::array<double, 6>, double, 6>(
  //      plyPropertyBeamProperties, nanoply::NNP_LIST_INT8_FLOAT64, nullptr);

  //  VCGEdgeMesh::PerEdgeAttributeHandle<ElementMaterial>
  //  handleBeamRigidityContants;
  //  customAttrib.AddEdgeAttribDescriptor<vcg::Point4f, float, 4>(
  //      plyPropertyBeamRigidityConstantsID, nanoply::NNP_LIST_INT8_FLOAT32,
  //      nullptr);
  //  unsigned int mask = 0;
  //  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_VERTCOORD;
  //  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_VERTNORMAL;
  //  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_EDGEINDEX;
  //  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_EDGEATTRIB;
  //  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_MESHATTRIB;
  if (!VCGEdgeMesh::load(plyFilePath)) {
    return false;
  }
  removeDuplicateVertices();

  //  elements =
  //  vcg::tri::Allocator<SimulationEdgeMesh>::GetPerEdgeAttribute<Element>(
  //      *this, std::string("Elements"));
  //  nodes =
  //  vcg::tri::Allocator<SimulationEdgeMesh>::GetPerVertexAttribute<Node>(
  //      *this, std::string("Nodes"));
  vcg::tri::UpdateTopology<SimulationEdgeMesh>::VertexEdge(*this);
  initializeNodes();
  initializeElements();
  //  setBeamMaterial(0.3, 1 * 1e9);
  updateEigenEdgeAndVertices();

  //  if (!handleBeamProperties._handle->data.empty()) {
  //    for (size_t ei = 0; ei < EN(); ei++) {
  //      elements[ei] =
  //      Element::Properties(handleBeamProperties[ei]);
  //      elements[ei].updateRigidity();
  //    }
  //  }
  //  for (size_t ei = 0; ei < EN(); ei++) {
  //    elements[ei].setDimensions(handleBeamDimensions[ei]);
  //    elements[ei].setMaterial(handleBeamMaterial[ei]);
  //    elements[ei].updateRigidity();
  //  }

  bool normalsAreAbsent = vert[0].cN().Norm() < 0.000001;
  if (normalsAreAbsent) {
    CoordType normalVector(0, 0, 1);
    std::cout << "Warning: Normals are missing from " << plyFilePath
              << ". Added normal vector:" << toString(normalVector)
              << std::endl;
    for (auto& v : vert) {
      v.N() = normalVector;
    }
  }

  return true;
}

bool SimulationEdgeMesh::save(const std::filesystem::path& plyFilePath) {
  //  std::filesystem::path filePath = plyFilePath;
  //  if (filePath.empty()) {
  //    filePath =
  //        std::filesystem::current_path().append(getLabel() +
  //        ".ply").string();
  //  }
  //  nanoply::NanoPlyWrapper<VCGEdgeMesh>::CustomAttributeDescriptor
  //  customAttrib; customAttrib.GetMeshAttrib(filePath);

  //  std::vector<CrossSectionType> dimensions = getBeamDimensions();
  //  customAttrib.AddEdgeAttribDescriptor<CrossSectionType, double, 2>(
  //      plyPropertyBeamDimensionsID, nanoply::NNP_LIST_INT8_FLOAT64,
  //      dimensions.data());
  //  std::vector<ElementMaterial> material = getBeamMaterial();
  //  customAttrib.AddEdgeAttribDescriptor<vcg::Point2d, double, 2>(
  //      plyPropertyBeamMaterialID, nanoply::NNP_LIST_INT8_FLOAT64,
  //      material.data());
  //  unsigned int mask = 0;
  //  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_VERTCOORD;
  //  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_EDGEINDEX;
  //  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_EDGEATTRIB;
  //  mask |= nanoply::NanoPlyWrapper<VCGEdgeMesh>::IO_VERTNORMAL;
  //  if (nanoply::NanoPlyWrapper<VCGEdgeMesh>::SaveModel(
  //          filePath.c_str(), *this, mask, customAttrib, false) != 1) {
  //    return false;
  //  }
  if (!VCGEdgeMesh::save(plyFilePath)) {
    return false;
  }
  return true;
}

SimulationEdgeMesh::EdgePointer SimulationEdgeMesh::getReferenceElement(
    const VCGEdgeMesh::VertexType& v) {
  const VertexIndex vi = getIndex(v);
  return nodes[v].incidentElements[0];
}

VectorType computeT1Vector(const SimulationEdgeMesh::EdgeType& e) {
  return computeT1Vector(e.cP(0), e.cP(1));
}

VectorType computeT1Vector(const CoordType& p0, const CoordType& p1) {
  const VectorType t1 = (p1 - p0).Normalize();
  return t1;
}

Element::LocalFrame computeElementFrame(const CoordType& p0,
                                        const CoordType& p1,
                                        const VectorType& elementNormal) {
  const VectorType t1 = computeT1Vector(p0, p1);
  const VectorType t2 = (elementNormal ^ t1).Normalize();
  const VectorType t3 = (t1 ^ t2).Normalize();

  return Element::LocalFrame{t1, t2, t3};
}

double computeAngle(const VectorType& v0,
                    const VectorType& v1,
                    const VectorType& normalVector) {
  //  double cosAngle = vector0.dot(vector1);
  //  const double epsilon = std::pow(10, -6);
  //  if (abs(cosAngle) > 1 && abs(cosAngle) < 1 + epsilon) {
  //    if (cosAngle > 0) {
  //      cosAngle = 1;

  //    } else {
  //      cosAngle = -1;
  //    }
  //  }
  //  assert(abs(cosAngle) <= 1);
  //  const double angle =
  //      acos(cosAngle);  // NOTE: I compute the alpha angle not between
  //                       // two consecutive elements but rather between
  //                       // the first and the ith. Is this correct?
  //  assert(!std::isnan(angle));

  //  const VectorType cp = vector0 ^ vector1;
  //  if (cp.dot(normalVector) < 0) {
  //    return -angle;
  //  }
  //  return angle;

  double angle = std::atan2(v1[1], v1[0]) - std::atan2(v0[1], v0[0]);
  if (angle < 0) {
    angle += 2 * M_PI;
  }
  return angle;
}

// void Element::computeMaterialProperties(const ElementMaterial &material) {
//    G = material.youngsModulus / (2 * (1 + material.poissonsRatio));
//}

// void Element::computeCrossSectionArea(const RectangularBeamDimensions
// &dimensions, double &A)
//{
//    A = dimensions.b * dimensions.h;
//}

// void Element::computeDimensionsProperties(
//    const RectangularBeamDimensions &dimensions) {
//  assert(typeid(CrossSectionType) == typeid(RectangularBeamDimensions));
//  computeCrossSectionArea(dimensions, A);
//  computeMomentsOfInertia(dimensions, dimensions.inertia);
//}

// void Element::computeDimensionsProperties(
//    const CylindricalBeamDimensions &dimensions) {
//  assert(typeid(CrossSectionType) == typeid(CylindricalBeamDimensions));
//  A = M_PI * (std::pow(dimensions.od / 2, 2) - std::pow(dimensions.id / 2,
//  2)); dimensions.inertia.I2 = M_PI * (std::pow(dimensions.od, 4) -
//  std::pow(dimensions.id, 4)) / 64; dimensions.inertia.I3 =
//  dimensions.inertia.I2; dimensions.inertia.J = dimensions.inertia.I2 +
//  dimensions.inertia.I3;
//}

void Element::setDimensions(const CrossSectionType& dimensions) {
  this->dimensions = dimensions;
  assert(this->dimensions.A == dimensions.A);
  updateRigidity();
}

void Element::setMaterial(const ElementMaterial& material) {
  this->material = material;
  updateRigidity();
}

double Element::getMass(const double& materialDensity) {
  const double beamVolume = dimensions.A * length;
  return beamVolume * materialDensity;
}

void Element::updateRigidity() {
  //    assert(initialLength != 0);
  rigidity.axial = material.youngsModulus * dimensions.A / initialLength;
  //    assert(rigidity.axial != 0);
  rigidity.torsional = material.G * dimensions.inertia.J / initialLength;
  //    assert(rigidity.torsional != 0);
  rigidity.firstBending =
      2 * material.youngsModulus * dimensions.inertia.I2 / initialLength;
  //    assert(rigidity.firstBending != 0);
  rigidity.secondBending =
      2 * material.youngsModulus * dimensions.inertia.I3 / initialLength;
  //    assert(rigidity.secondBending != 0);
}
