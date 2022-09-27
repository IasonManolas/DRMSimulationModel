#ifndef SIMULATIONMESH_HPP
#define SIMULATIONMESH_HPP

#include "edgemesh.hpp"

struct Element;
struct Node;
// TODO: use templates for different cross section types
using CrossSectionType = RectangularBeamDimensions;
// using CrossSectionType = CylindricalBeamDimensions;

class SimulationEdgeMesh : public VCGEdgeMesh {
 private:
  void computeElementalProperties();
  void initializeElements();
  void initializeNodes();
  EdgePointer getReferenceElement(const VertexType& v);

  const std::string plyPropertyBeamDimensionsID{"beam_dimensions"};
  const std::string plyPropertyBeamMaterialID{"beam_material"};

 public:
  ~SimulationEdgeMesh();
  SimulationEdgeMesh();
  SimulationEdgeMesh(ConstVCGEdgeMesh& edgeMesh);
  SimulationEdgeMesh(SimulationEdgeMesh& elementalMesh);
  void updateElementalLengths();
  void updateIncidentElements();
  void updateElementalFrames();
  bool save(const std::filesystem::path& plyFilePath = std::string());
  virtual bool load(const std::filesystem::path& plyFilePath);
  void reset();
  std::vector<VCGEdgeMesh::EdgePointer> getIncidentElements(
      const VCGEdgeMesh::VertexType& v);
  std::vector<CrossSectionType> getBeamDimensions();
  std::vector<ElementMaterial> getBeamMaterial();
  void setBeamCrossSection(const CrossSectionType& beamDimensions);
  void setBeamMaterial(const double& pr, const double& ym);
  void setBeamMaterial(const ElementMaterial& material);

  PerEdgeAttributeHandle<Element> elements;
  PerVertexAttributeHandle<Node> nodes;
  double pre_previousTotalKineticEnergy{0};
  double pre_previousTotalTranslationalKineticEnergy{0};
  double pre_previousTotalRotationalKineticEnergy{0};
  double previousTotalKineticEnergy{0};
  double previousTotalTranslationalKineticEnergy{0};
  double previousTotalRotationalKineticEnergy{0};
  double previousTotalResidualForcesNorm{0};
  double currentTotalKineticEnergy{0};
  double currentTotalTranslationalKineticEnergy{0};
  double currentTotalRotationalKineticEnergy{0};
  double totalResidualForcesNorm{0};
  double totalExternalForcesNorm{0};
  double averageResidualForcesNorm{0};
  double currentTotalPotentialEnergykN{0};
  double previousTotalPotentialEnergykN{0};
  double residualForcesMovingAverageDerivativeNorm{0};
  double residualForcesMovingAverage{0};
  double perVertexAverageNormalizedDisplacementNorm{0};

#ifdef POLYSCOPE_DEFINED
  polyscope::CurveNetwork* registerForDrawing(
      const std::optional<std::array<float, 3>>& desiredColor = std::nullopt,
      const bool& shouldEnable = true,
      const std::optional<float>& inputDrawingRadius = std::nullopt);
  void unregister() const;
#endif
};

struct Element {
  CrossSectionType dimensions;
  ElementMaterial material;

  void computeMaterialProperties(const ElementMaterial& material);
  //  void computeDimensionsProperties(const RectangularBeamDimensions
  //  &dimensions); void computeDimensionsProperties(const
  //  CylindricalBeamDimensions &dimensions);
  void setDimensions(const CrossSectionType& dimensions);
  void setMaterial(const ElementMaterial& material);
  double getMass(const double& matrialDensity);

  struct LocalFrame {
    VectorType t1;
    VectorType t2;
    VectorType t3;
  };

  EdgeIndex ei;
  double length{0};
  double initialLength;
  LocalFrame frame;

  struct Rigidity {
    double axial;
    double torsional;
    double firstBending;
    double secondBending;
    std::string toString() const {
      return std::string("Rigidity:") + std::string("\nAxial=") +
             std::to_string(axial) + std::string("\nTorsional=") +
             std::to_string(torsional) + std::string("\nFirstBending=") +
             std::to_string(firstBending) + std::string("\nSecondBending=") +
             std::to_string(secondBending);
    }
  };
  Rigidity rigidity;
  void updateRigidity();

  VectorType f1_j;
  VectorType f1_jplus1;
  VectorType f2_j;
  VectorType f2_jplus1;
  VectorType f3_j;
  VectorType f3_jplus1;
  double cosRotationAngle_j;
  double cosRotationAngle_jplus1;
  double sinRotationAngle_j;
  double sinRotationAngle_jplus1;
  std::vector<std::vector<VectorType>> derivativeT1;
  std::vector<std::vector<VectorType>> derivativeT2;
  std::vector<std::vector<VectorType>> derivativeT3;
  std::vector<std::vector<VectorType>> derivativeR;
  struct RotationalDisplacements {
    double theta1{0}, theta2{0}, theta3{0};
  };
  RotationalDisplacements rotationalDisplacements_j;
  RotationalDisplacements rotationalDisplacements_jplus1;

  static void computeCrossSectionArea(
      const RectangularBeamDimensions& dimensions,
      double& A);
};

struct Node {
  struct Forces {
    Vector6d external{0};
    Vector6d internal{0};
    Vector6d residual{0};
    bool hasExternalForce() const { return external.isZero(); }
  };

  struct Mass {
    double translational;
    double rotational;
  };

  Mass mass;
  //  Vector6d mass_6d;
  Vector6d damping_6d;
  VertexIndex vi;
  CoordType initialLocation;
  CoordType initialNormal;
  Vector6d acceleration{0};
  Forces force;
  Vector6d previousVelocity{0};
  Vector6d velocity{0};
  double kineticEnergy{0};
  Vector6d displacements{0};
  double nR{0};
  //  std::unordered_map<EdgeIndex, double>
  //      alphaAngles; // contains the initial angles between the first star
  //      element
  //                   // incident to this node and the other elements of the
  //                   star
  //  // has size equal to the valence of the vertex
  std::vector<std::pair<EdgeIndex, double>> alphaAngles;

  std::vector<VCGEdgeMesh::EdgePointer> incidentElements;
  std::vector<VectorType> derivativeOfNormal;
  SimulationEdgeMesh::EdgePointer referenceElement;
};

Element::LocalFrame computeElementFrame(const CoordType& p0,
                                        const CoordType& p1,
                                        const VectorType& elementNormal);
VectorType computeT1Vector(const SimulationEdgeMesh::EdgeType& e);

VectorType computeT1Vector(const CoordType& p0, const CoordType& p1);
double computeAngle(const VectorType& vector0,
                    const VectorType& vector1,
                    const VectorType& normalVector);

#endif  // ELEMENTALMESH_HPP
