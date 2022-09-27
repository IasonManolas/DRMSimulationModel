#ifndef BEAMFORMFINDER_HPP
#define BEAMFORMFINDER_HPP

#include <chrono>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <unordered_set>
#include "edgemesh.hpp"
#include "simulation_structs.hpp"

class SimulationJob;

using EdgeType = VCGEdgeMesh::EdgeType;
using VertexType = VCGEdgeMesh::VertexType;

struct DifferentiateWithRespectTo {
  const VertexType& v;
  const int& dofi;
};

class DRMSimulationModel {
 public:
  inline const static std::string label{"DRM"};
  enum DoF { Ux = 0, Uy, Uz, Nx, Ny, Nr, NumDoF };
  using DoFType = int;
  struct Settings {
    inline static std::string jsonDefaultFileName{"DRMSettings.json"};
    double Dtini{0.1};
    double xi{0.9969};
    double gamma{1};
    std::optional<double> threshold_residualForcesNorm;
    std::optional<double> threshold_averageResidualToExternalForcesNorm;
    std::optional<double> threshold_totalTranslationalKineticEnergy;
    std::optional<int> maxDRMIterations;
    std::optional<TemporaryLoad> initialDistortion;
    std::optional<StepsDuration> gradualForcedDisplacementSteps;
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(
        Settings,
        Dtini,
        xi,
        gamma,
        threshold_residualForcesNorm,
        threshold_averageResidualToExternalForcesNorm,
        threshold_totalTranslationalKineticEnergy,
        maxDRMIterations,
        initialDistortion,
        gradualForcedDisplacementSteps)

    Settings() {}
    bool hasConvergenceCriterion() const;
    void printConvergenceCriteria() const;
    void save(const std::filesystem::path& jsonPath) const;
    bool load(const std::filesystem::path& jsonFilePath);
    std::optional<int> debugModeStep;
    bool shouldDraw{false};
    bool beVerbose{false};
    bool shouldCreatePlots{false};
  };

 private:
  std::string postConvergenceMessage;
#ifdef POLYSCOPE_DEFINED
  RGBColor color_deformed{67.0 / 255, 160.00 / 255, 232.0 / 255};
  RGBColor color_initial{24.0 / 255, 23.0 / 255, 23.0 / 255};
#endif
  std::chrono::steady_clock::time_point beginTime;

  SimulationHistory history;

  std::vector<double> plotYValues;

  std::vector<bool> isVertexConstrained;
  std::vector<bool> isRigidSupport;
  bool shouldApplyInitialDistortion{false};
  Settings mSettings;
  double Dt{mSettings.Dtini};
  bool kineticDampingWasAppliedInThisIteration{false};
  bool checkedForMaximumMoment{false};
  double externalMomentsNorm{0};
  size_t mCurrentSimulationStep{0};
  size_t numOfDampings{0};
  double minTotalResidualForcesNorm{std::numeric_limits<double>::max()};
  double minLocalMaximaOfTranslationaKineticEnergy{
      std::numeric_limits<double>::max()};
  double firstTranslationalKineticEnergyPeak{-1};
  bool firstTranslationalKineticEnergyPeakOccured{false};
  const std::string polyscopeLabel_deformed{"Simulation mesh"};
  std::unique_ptr<SimulationEdgeMesh> pMesh;
  std::unordered_map<VertexIndex,
                     std::unordered_set<DRMSimulationModel::DoFType>>
      constrainedVertices;
  void reset(const std::shared_ptr<SimulationJob>& pJob,
             const Settings& settings);
  void updateNodalExternalForces(
      const std::unordered_map<VertexIndex, Vector6d>& nodalForces);
  void updateRotationalDisplacements();
  void updateElementalLengths();

  void updateNodalMasses();

  void updateNodalAccelerations();

  void updateNodalVelocities();

  void updateNodalDisplacements();

  void applyForcedDisplacements(
      const std::unordered_map<VertexIndex, Eigen::Vector3d>&
          nodalForcedDisplacements);

  void updateKineticEnergy();

  void resetVelocities();

  SimulationResults computeResults(const std::shared_ptr<SimulationJob>& pJob);

  void updateNodePosition(VertexType& v);

  void applyDisplacements();

#ifdef POLYSCOPE_DEFINED
  void draw(const std::shared_ptr<SimulationJob>& pJob,
            const std::string& screenshotsFolder = {});
#endif
  void updateNodalInternalForce(
      Vector6d& nodalInternalForce,
      const Vector6d& elementInternalForce,
      const std::unordered_set<DRMSimulationModel::DoFType>& nodalFixedDof);

  Vector6d computeElementInternalForce(
      const Element& elem,
      const Node& n0,
      const Node& n1,
      const std::unordered_set<DRMSimulationModel::DoFType>& n0ConstrainedDof,
      const std::unordered_set<DRMSimulationModel::DoFType>& n1ConstrainedDof);

  Vector6d computeElementAxialForce(const ::EdgeType& e) const;
  VectorType computeDisplacementDifferenceDerivative(
      const EdgeType& e,
      const DifferentiateWithRespectTo& dui) const;
  double computeDerivativeElementLength(
      const EdgeType& e,
      const DifferentiateWithRespectTo& dui) const;

  VectorType computeDerivativeT1(const EdgeType& e,
                                 const DifferentiateWithRespectTo& dui) const;

  VectorType computeDerivativeOfNormal(
      const VertexType& v,
      const DifferentiateWithRespectTo& dui) const;

  VectorType computeDerivativeT3(const EdgeType& e,
                                 const DifferentiateWithRespectTo& dui) const;

  VectorType computeDerivativeT2(const EdgeType& e,
                                 const DifferentiateWithRespectTo& dui) const;

  double computeDerivativeTheta2(
      const EdgeType& e,
      const VertexIndex& evi,
      const VertexIndex& dwrt_evi,
      const DRMSimulationModel::DoFType& dwrt_dofi) const;

  void updateElementalFrames();

  VectorType computeDerivativeOfR(const EdgeType& e,
                                  const DifferentiateWithRespectTo& dui) const;

  static double computeDerivativeOfNorm(const VectorType& x,
                                        const VectorType& derivativeOfX);
  static VectorType computeDerivativeOfCrossProduct(
      const VectorType& a,
      const VectorType& derivativeOfA,
      const VectorType& b,
      const VectorType& derivativeOfB);

  double computeTheta3(const EdgeType& e, const VertexType& v);
  double computeDerivativeTheta3(const EdgeType& e,
                                 const VertexType& v,
                                 const DifferentiateWithRespectTo& dui) const;
  double computeTotalPotentialEnergy();
  void computeRigidSupports();
  void updateNormalDerivatives();
  void updateT1Derivatives();
  void updateT2Derivatives();
  void updateT3Derivatives();
  void updateRDerivatives();

  double computeDerivativeTheta1(
      const EdgeType& e,
      const VertexIndex& evi,
      const VertexIndex& dwrt_evi,
      const DRMSimulationModel::DoFType& dwrt_dofi) const;

  void updateNodeNormal(VertexType& v);

  // TODO: move to experimental branch
  void applyForcedNormals(
      const std::unordered_map<VertexIndex, VectorType> nodalForcedRotations);

  void printCurrentState() const;

  void printDebugInfo() const;

  void applySolutionGuess(const SimulationResults& solutionGuess,
                          const std::shared_ptr<SimulationJob>& pJob);

  void updateNodeNr(VertexType& v);

  void applyKineticDamping(const std::shared_ptr<SimulationJob>& pJob);

  void reset(const std::shared_ptr<SimulationJob>& pJob);

  void updateDerivatives();

  bool convergedUsingResidualForcesCriteria();
  void updateResidualForces();

 public:
  DRMSimulationModel();
  SimulationResults executeSimulation(
      const std::shared_ptr<SimulationJob>& pJob,
      const Settings& settings);
  //#ifdef USE_ENSMALLEN
  //  std::shared_ptr<SimulationJob> pJob;
  //  double EvaluateWithGradient(const arma::mat &x, arma::mat &g);
  //  void setJob(const std::shared_ptr<SimulationJob> &pJob);
  //  SimulationEdgeMesh *getDeformedMesh(const arma::mat &x, const
  //  std::shared_ptr<SimulationJob> &pJob); double Evaluate(const arma::mat
  //  &x);
  //#endif

  static void runUnitTests();
  void reset(const Settings& settings);
  void updateNodalExternalForces(
      const std::unordered_map<VertexIndex, Vector6d>& nodalForces,
      const std::unordered_map<VertexIndex, std::unordered_set<DoFType>>&
          fixedVertices);
  void updateNodePosition(
      VertexType& v,
      const std::unordered_map<VertexIndex, std::unordered_set<DoFType>>&
          fixedVertices);
  void updateNodeNormal(
      VertexType& v,
      const std::unordered_map<VertexIndex, std::unordered_set<DoFType>>&
          fixedVertices);
  void applyDisplacements(
      const std::unordered_map<VertexIndex, std::unordered_set<DoFType>>&
          fixedVertices);
};

namespace nlohmann {}  // namespace nlohmann

template <typename PointType>
PointType Cross(PointType p1, PointType p2) {
  return p1 ^ p2;
}

inline size_t currentStep{0};
#endif  // BEAMFORMFINDER_HPP
