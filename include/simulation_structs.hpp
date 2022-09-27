#ifndef SIMULATIONSTRUCTS_HPP
#define SIMULATIONSTRUCTS_HPP
//#include "csvfile.hpp"
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>
#include "nlohmann/json.hpp"
#include "simulationmesh.hpp"
#include "utilities.hpp"

#ifdef POLYSCOPE_DEFINED
#include <polyscope/point_cloud.h>
#endif
namespace Eigen {
template <class Matrix>
void writeBinary(const std::filesystem::path& filePath, const Matrix& matrix) {
  std::ofstream out(filePath,
                    std::ios::out | std::ios::binary | std::ios::trunc);
  typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
  out.write((char*)(&rows), sizeof(typename Matrix::Index));
  out.write((char*)(&cols), sizeof(typename Matrix::Index));
  out.write((char*)matrix.data(),
            rows * cols * sizeof(typename Matrix::Scalar));
  out.close();
}
template <class Matrix>
void readBinary(const std::filesystem::path& filePath, Matrix& matrix) {
  std::ifstream in(filePath, std::ios::in | std::ios::binary);
  typename Matrix::Index rows = 0, cols = 0;
  in.read((char*)(&rows), sizeof(typename Matrix::Index));
  in.read((char*)(&cols), sizeof(typename Matrix::Index));
  matrix.resize(rows, cols);
  in.read((char*)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
  in.close();
}
// const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
// template<class Matrix>
// void writeToCSV(const std::string &filename, Matrix &matrix)
//{
//    ofstream file(filename.c_str());
//    file << matrix.format(CSVFormat);
//}
}  // namespace Eigen

struct SimulationHistory {
  SimulationHistory() {}

  size_t numberOfSteps{0};
  std::string label;
  std::vector<double> logResidualForces;
  std::vector<double> kineticEnergy;
  std::vector<double> potentialEnergies;
  std::vector<size_t> redMarks;
  std::vector<double> greenMarks;

  void markRed(const size_t& stepNumber) { redMarks.push_back(stepNumber); }

  void markGreen(const size_t& stepNumber) { greenMarks.push_back(stepNumber); }

  void stepPulse(const SimulationEdgeMesh& mesh) {
    kineticEnergy.push_back(std::log10(mesh.currentTotalKineticEnergy));
    logResidualForces.push_back(std::log10(mesh.totalResidualForcesNorm));
  }

  void clear() {
    logResidualForces.clear();
    kineticEnergy.clear();
  }
};

struct StepsDuration {
  int startStep{0};
  int endStep{0};
  bool isActive(const int& currentStep) const {
    return currentStep >= startStep && currentStep < endStep;
  }
  int getSize() const { return endStep - startStep; }
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(StepsDuration, startStep, endStep)
};

struct TemporaryLoad {
  StepsDuration duration;
  Vector6d load;

  bool isActive(const int& currentStep) const {
    return duration.isActive(currentStep);
  }
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(TemporaryLoad, duration, load)
};

namespace nlohmann {

template <>
struct adl_serializer<std::unordered_map<VertexIndex, Vector6d>> {
  static void to_json(json& j,
                      const std::unordered_map<VertexIndex, Vector6d>& value) {
    // calls the "to_json" method in T's namespace
  }

  static void from_json(const nlohmann::json& j,
                        std::unordered_map<VertexIndex, Vector6d>& m) {
    std::cout << "Entered." << std::endl;
    for (const auto& p : j) {
      m.emplace(p.at(0).template get<VertexIndex>(),
                p.at(1).template get<std::array<double, 6>>());
    }
  }
};
}  // namespace nlohmann

class SimulationJob {
  //  const std::unordered_map<VertexIndex, VectorType> nodalForcedNormals;
  // json labels
  struct JSONLabels {
    inline static std::string meshFilename{"mesh filename"};
    inline static std::string forcedDisplacements{"forced displacements"};
    inline static std::string constrainedVertices{"fixed vertices"};
    inline static std::string nodalForces{"forces"};
    inline static std::string label{"label"};
    inline static std::string meshLabel{
        "meshLabel"};  // TODO: should be in the savePly function of the
    // simulation mesh class
    inline static std::string crossSectionDimensionLabel{"cross section dim"};
    inline static std::string material{"beam material"};

  } jsonLabels;

#ifdef POLYSCOPE_DEFINED
  const std::string polyscopeLabel_bcAsPointCloud{"BC_spheres"};
#endif

 public:
  inline static std::string jsonDefaultFileName = "SimulationJob.json";
  std::shared_ptr<SimulationEdgeMesh> pMesh;
  std::string label{"empty_job"};
  std::unordered_map<VertexIndex, std::unordered_set<int>> constrainedVertices;
  std::unordered_map<VertexIndex, Vector6d> nodalExternalForces;
  std::unordered_map<VertexIndex, Eigen::Vector3d> nodalForcedDisplacements;
  SimulationJob(
      const std::shared_ptr<SimulationEdgeMesh>& m,
      const std::string& label,
      const std::unordered_map<VertexIndex, std::unordered_set<int>>& cv,
      const std::unordered_map<VertexIndex, Vector6d>& ef = {},
      const std::unordered_map<VertexIndex, Eigen::Vector3d>& fd = {})
      : pMesh(m),
        label(label),
        constrainedVertices(cv),
        nodalExternalForces(ef),
        nodalForcedDisplacements(fd) {}
  SimulationJob() { pMesh = std::make_shared<SimulationEdgeMesh>(); }
  SimulationJob(const std::string& jsonFilename) { load(jsonFilename); }

  bool isEmpty() {
    return constrainedVertices.empty() && nodalExternalForces.empty() &&
           nodalForcedDisplacements.empty() && pMesh == nullptr;
  }
  void remap(const std::unordered_map<size_t, size_t>& sourceToDestinationViMap,
             SimulationJob& destination_simulationJob) const {
    std::unordered_map<VertexIndex, std::unordered_set<int>>
        destination_fixedVertices;
    for (const auto& source_fixedVertex : this->constrainedVertices) {
      destination_fixedVertices[sourceToDestinationViMap.at(
          source_fixedVertex.first)] = source_fixedVertex.second;
    }

    std::unordered_map<VertexIndex, Vector6d> destination_nodalForces;
    for (const auto& source_nodalForces : this->nodalExternalForces) {
      destination_nodalForces[sourceToDestinationViMap.at(
          source_nodalForces.first)] = source_nodalForces.second;
    }

    std::unordered_map<VertexIndex, Eigen::Vector3d>
        destination_forcedDisplacements;
    for (const auto& source_forcedDisplacements :
         this->nodalForcedDisplacements) {
      destination_forcedDisplacements[sourceToDestinationViMap.at(
          source_forcedDisplacements.first)] =
          source_forcedDisplacements.second;
    }

    destination_simulationJob.constrainedVertices = destination_fixedVertices;
    destination_simulationJob.nodalExternalForces = destination_nodalForces;
    destination_simulationJob.label = this->getLabel();
    destination_simulationJob.nodalForcedDisplacements =
        destination_forcedDisplacements;
  }

  SimulationJob getCopy() const {
    SimulationJob jobCopy;
    jobCopy.pMesh = std::make_shared<SimulationEdgeMesh>();
    jobCopy.pMesh->copy(*pMesh);
    jobCopy.label = label;
    jobCopy.constrainedVertices = constrainedVertices;
    jobCopy.nodalExternalForces = nodalExternalForces;
    jobCopy.nodalForcedDisplacements = nodalForcedDisplacements;

    return jobCopy;
  }

  std::string getLabel() const { return label; }

  std::string toString() const {
    nlohmann::json json;
    if (!constrainedVertices.empty()) {
      json[jsonLabels.constrainedVertices] = constrainedVertices;
    }
    if (!nodalExternalForces.empty()) {
      std::unordered_map<VertexIndex, std::array<double, 6>> arrForces;
      for (const auto& f : nodalExternalForces) {
        arrForces[f.first] = f.second;
      }
      json[jsonLabels.nodalForces] = arrForces;
    }

    return json.dump();
  }
  bool operator==(const SimulationJob& otherSimulationJob) {
    return this->toString() == otherSimulationJob.toString();
  }

  void clear() {
    label = "empty_job";
    constrainedVertices.clear();
    nodalExternalForces.clear();
    nodalForcedDisplacements.clear();
    if (pMesh.use_count() == 1) {
      std::cout << "Job mesh is deleted" << std::endl;
    }
    pMesh.reset();
  }

  bool load(const std::string& jsonFilename,
            const bool& shouldLoadMesh = true) {
    label = "empty_job";
    constrainedVertices.clear();
    nodalExternalForces.clear();
    nodalForcedDisplacements.clear();
    const bool beVerbose = false;
    if (std::filesystem::path(jsonFilename).extension() != ".json") {
      std::cerr << "A json file is expected as input. The given file has the "
                   "following extension:"
                << std::filesystem::path(jsonFilename).extension() << std::endl;
      assert(false);
      return false;
    }

    if (!std::filesystem::exists(std::filesystem::path(jsonFilename))) {
      std::cerr << "The json file does not exist. Json file provided:"
                << jsonFilename << std::endl;
      assert(false);
      return false;
    }

    if (beVerbose) {
      std::cout << "Loading json file:" << jsonFilename << std::endl;
    }
    nlohmann::json json;
    std::ifstream ifs(jsonFilename);
    ifs >> json;

    if (shouldLoadMesh) {
      pMesh.reset();
      pMesh = std::make_shared<SimulationEdgeMesh>();
      if (json.contains(jsonLabels.meshFilename)) {
        const std::string relativeFilepath = json[jsonLabels.meshFilename];
        const auto meshFilepath =
            std::filesystem::path(
                std::filesystem::path(jsonFilename).parent_path())
                .append(relativeFilepath);
        pMesh->load(meshFilepath.string());
        pMesh->setLabel(
            json[jsonLabels.meshLabel]);  // FIXME: This should be exported
                                          // using nanoply but nanoply might not
                                          // be able to write a string(??)
      }

      if (json.contains(jsonLabels.meshLabel)) {
        pMesh->setLabel(json[jsonLabels.meshLabel]);
      }
    }

    if (json.contains(jsonLabels.constrainedVertices)) {
      constrainedVertices =
          //      auto conV =
          json[jsonLabels.constrainedVertices]
              .get<std::unordered_map<VertexIndex, std::unordered_set<int>>>();
      if (beVerbose) {
        std::cout << "Loaded constrained vertices. Number of constrained "
                     "vertices found:"
                  << constrainedVertices.size() << std::endl;
      }
    }

    if (json.contains(jsonLabels.nodalForces)) {
      auto f =
          json[jsonLabels.nodalForces]
              .get<std::unordered_map<VertexIndex, std::array<double, 6>>>();
      for (const auto& force : f) {
        nodalExternalForces[force.first] = Vector6d(force.second);
      }
      if (beVerbose) {
        std::cout << "Loaded forces. Number of forces found:"
                  << nodalExternalForces.size() << std::endl;
      }
    }

    if (json.contains(jsonLabels.forcedDisplacements)) {
      //      auto conV =
      std::unordered_map<VertexIndex, std::array<double, 3>> forcedDisp =
          json[jsonLabels.forcedDisplacements]
              .get<std::unordered_map<VertexIndex, std::array<double, 3>>>();

      for (const auto& fd : forcedDisp) {
        nodalForcedDisplacements[fd.first] =
            Eigen::Vector3d(fd.second[0], fd.second[1], fd.second[2]);
      }
      if (beVerbose) {
        std::cout << "Loaded forced displacements. Number of forced displaced"
                     "vertices found:"
                  << nodalForcedDisplacements.size() << std::endl;
      }
    }

    if (json.contains(jsonLabels.label)) {
      label = json[jsonLabels.label];
    }

    if (json.contains(jsonLabels.crossSectionDimensionLabel)) {
      std::pair<float, float> crossSectionDimensions =
          json[jsonLabels.crossSectionDimensionLabel]
              .get<std::pair<float, float>>();
      pMesh->setBeamCrossSection(RectangularBeamDimensions(
          {crossSectionDimensions.first, crossSectionDimensions.second}));
    }

    if (json.contains(jsonLabels.material)) {
      std::pair<float, float> poissonsRatioYoungsModulusPair =
          json[jsonLabels.material].get<std::pair<float, float>>();
      pMesh->setBeamMaterial(
          ElementMaterial({poissonsRatioYoungsModulusPair.first,
                           poissonsRatioYoungsModulusPair.second}));
    }

    return true;
  }

  bool save(const std::string& folderDirectory) const {
    const std::filesystem::path pathFolderDirectory(folderDirectory);
    if (!std::filesystem::is_directory(pathFolderDirectory)) {
      std::cerr << "A folder directory is expected for saving the simulation "
                   "job. Exiting.."
                << std::endl;
      return false;
    }

    bool returnValue = true;
    std::string jsonFilename(
        std::filesystem::path(pathFolderDirectory)
            //            .append(label + "_" + pMesh->getLabel() +
            //            "_simulationJob.json")
            .append("SimulationJob.json")
            .string());

    const std::string meshFilename =
        std::filesystem::absolute(
            std::filesystem::canonical(
                std::filesystem::path(pathFolderDirectory)))
            .append(pMesh->getLabel() + ".ply")
            .string();
    returnValue = pMesh->save(meshFilename);
    nlohmann::json json;
    json[jsonLabels.meshFilename] =
        std::filesystem::relative(
            std::filesystem::path(meshFilename),
            std::filesystem::path(jsonFilename).parent_path())
            .string();
    json[jsonLabels.meshLabel] =
        pMesh->getLabel();  // FIXME: This should be exported using nanoply but
                            // nanoply might not be able to write a string(??)
    if (!constrainedVertices.empty()) {
      json[jsonLabels.constrainedVertices] = constrainedVertices;
    }
    if (!nodalExternalForces.empty()) {
      std::unordered_map<VertexIndex, std::array<double, 6>> arrForces;
      for (const auto& f : nodalExternalForces) {
        arrForces[f.first] = f.second;
      }
      json[jsonLabels.nodalForces] = arrForces;
    }
    if (!nodalForcedDisplacements.empty()) {
      std::unordered_map<VertexIndex, std::array<double, 3>> forcedDisp;
      for (const auto& fd : nodalForcedDisplacements) {
        forcedDisp[fd.first] = {fd.second[0], fd.second[1], fd.second[2]};
      }

      json[jsonLabels.forcedDisplacements] = forcedDisp;
    }

    if (!label.empty()) {
      json[jsonLabels.label] = label;
    }
    if (!pMesh->getLabel().empty()) {
      json[jsonLabels.meshLabel] = pMesh->getLabel();
    }
    //    if (!pMesh->getLabel().empty()) {
    json[jsonLabels.crossSectionDimensionLabel] =
        std::pair<float, float>({pMesh->getBeamDimensions()[0].dim1,
                                 pMesh->getBeamDimensions()[0].dim2});
    //    }

    json[jsonLabels.material] =
        std::pair<float, float>({pMesh->getBeamMaterial()[0].poissonsRatio,
                                 pMesh->getBeamMaterial()[0].youngsModulus});

    std::ofstream jsonFile(jsonFilename);
    jsonFile << json;
    jsonFile.close();
    //    std::cout << "Saved simulation job as:" << jsonFilename << std::endl;

    return returnValue;
  }
#ifdef POLYSCOPE_DEFINED
  void registerForDrawing(const std::string& meshLabel,
                          const bool& shouldEnable = false) const {
    PolyscopeInterface::init();
    if (meshLabel.empty()) {
      assert(false);
      std::cerr << "Expects a mesh label on which to draw the simulation job."
                << std::endl;
      return;
    }

    if (!polyscope::hasCurveNetwork(meshLabel)) {
      assert(false);
      std::cerr << "Expects mesh already being registered to draw the "
                   "simulation job. No struct named " +
                       meshLabel
                << std::endl;
      return;
    }
    std::vector<std::array<double, 3>> nodeColors(pMesh->VN());
    for (const auto& fixedVertex : constrainedVertices) {
      const bool hasRotationalDoFConstrained = fixedVertex.second.contains(3) ||
                                               fixedVertex.second.contains(4) ||
                                               fixedVertex.second.contains(5);
      const bool hasTranslationalDoFConstrained =
          fixedVertex.second.contains(0) || fixedVertex.second.contains(1) ||
          fixedVertex.second.contains(2);
      if (hasTranslationalDoFConstrained && !hasRotationalDoFConstrained) {
        nodeColors[fixedVertex.first] = {0, 0, 1};
      } else if (!hasTranslationalDoFConstrained &&
                 hasRotationalDoFConstrained) {
        nodeColors[fixedVertex.first] = {0, 1, 0};
      } else {
        nodeColors[fixedVertex.first] = {0, 1, 1};
      }
    }
    if (!nodalForcedDisplacements.empty()) {
      for (const std::pair<VertexIndex, Eigen::Vector3d>& viDisplPair :
           nodalForcedDisplacements) {
        const VertexIndex vi = viDisplPair.first;
        nodeColors[vi][0] += 1;
        nodeColors[vi][0] /= 2;
        nodeColors[vi][1] += 0;
        nodeColors[vi][1] /= 2;
        nodeColors[vi][2] += 0;
        nodeColors[vi][2] /= 2;
      }
    }
    std::for_each(nodeColors.begin(), nodeColors.end(),
                  [](std::array<double, 3>& color) {
                    const double norm =
                        sqrt(std::pow(color[0], 2) + std::pow(color[1], 2) +
                             std::pow(color[2], 2));
                    if (norm > std::pow(10, -7)) {
                      color[0] /= norm;
                      color[1] /= norm;
                      color[2] /= norm;
                    }
                  });

    if (!nodeColors.empty()) {
      polyscope::getCurveNetwork(meshLabel)
          ->addNodeColorQuantity("Boundary conditions", nodeColors)
          ->setEnabled(shouldEnable);
    }
    drawBcAsSpheres(polyscope::getCurveNetwork(meshLabel));

    // per node external forces
    std::vector<std::array<double, 3>> externalForces(pMesh->VN());
    for (const auto& forcePair : nodalExternalForces) {
      auto index = forcePair.first;
      auto force = forcePair.second;
      externalForces[index] = {force[0], force[1], force[2]};
    }

    if (!externalForces.empty()) {
      const std::string polyscopeLabel_externalForces = "External force";
      polyscope::getCurveNetwork(meshLabel)->removeQuantity(
          polyscopeLabel_externalForces);
      polyscope::CurveNetworkNodeVectorQuantity* externalForcesVectors =
          polyscope::getCurveNetwork(meshLabel)->addNodeVectorQuantity(
              polyscopeLabel_externalForces, externalForces);

      const std::array<double, 3> color_loads{1.0, 0, 0};
      externalForcesVectors->setVectorColor(
          glm::vec3(color_loads[0], color_loads[1], color_loads[2]));
      externalForcesVectors->setEnabled(shouldEnable);
    }
    // per node external moments
    bool hasExternalMoments = false;
    std::vector<std::array<double, 3>> externalMoments(pMesh->VN());
    for (const auto& forcePair : nodalExternalForces) {
      auto index = forcePair.first;
      const Vector6d& load = forcePair.second;
      if (load.getRotation().norm() != 0) {
        hasExternalMoments = true;
      }
      externalMoments[index] = {load[3], load[4], load[5]};
    }

    if (hasExternalMoments) {
      polyscope::getCurveNetwork(meshLabel)
          ->addNodeVectorQuantity("External moment", externalMoments)
          ->setEnabled(shouldEnable);
    }
  }
  void unregister(const std::string& meshLabel) const {
    if (polyscope::getCurveNetwork(meshLabel) == nullptr) {
      return;
    }
    if (!nodalExternalForces.empty()) {
      polyscope::getCurveNetwork(meshLabel)->removeQuantity("External force");
    }
    if (!constrainedVertices.empty()) {
      polyscope::getCurveNetwork(meshLabel)->removeQuantity(
          "Boundary conditions");
      polyscope::removeStructure(polyscopeLabel_bcAsPointCloud, false);
    }

    // per node external moments
    bool hasExternalMoments = false;
    for (const auto& forcePair : nodalExternalForces) {
      const Vector6d& load = forcePair.second;
      if (load.getRotation().norm() != 0) {
        hasExternalMoments = true;
        break;
      }
    }
    if (hasExternalMoments) {
      polyscope::getCurveNetwork(meshLabel)->removeQuantity("External moment");
    }
  }

  void drawBcAsSpheres(polyscope::CurveNetwork* polyscopeHandle_drawOn) const {
    polyscope::removeStructure(polyscopeLabel_bcAsPointCloud, false);
    std::vector<glm::vec3> bcPos;
    std::vector<glm::vec3> bcColors;
    const std::array<double, 3> color_rigid{63.0 / 255, 191.0 / 255,
                                            127.0 / 255};
    const std::array<double, 3> color_fixedTranslation({93.0 / 255, 0, 0});
    for (std::pair<VertexIndex, std::unordered_set<int>> bc :
         constrainedVertices) {
      bcPos.push_back(polyscopeHandle_drawOn->nodes[bc.first]);
      const bool hasRotationalDoFConstrained = bc.second.contains(3) ||
                                               bc.second.contains(4) ||
                                               bc.second.contains(5);
      const bool hasTranslationalDoFConstrained = bc.second.contains(0) ||
                                                  bc.second.contains(1) ||
                                                  bc.second.contains(2);

      if (hasTranslationalDoFConstrained && !hasRotationalDoFConstrained) {
        bcColors.push_back(glm::vec3(color_fixedTranslation[0],
                                     color_fixedTranslation[1],
                                     color_fixedTranslation[2]));
      } else if (!hasTranslationalDoFConstrained &&
                 hasRotationalDoFConstrained) {
        bcColors.push_back(glm::vec3(0, 0, 1));
      } else {
        bcColors.push_back(
            glm::vec3(color_rigid[0], color_rigid[1], color_rigid[2]));
      }
    }

    for (std::pair<VertexIndex, Eigen::Vector3d> bc :
         nodalForcedDisplacements) {
      bcPos.push_back(polyscopeHandle_drawOn->nodes[bc.first]);
      const bool hasTranslationalDoFConstrained = true;
      if (hasTranslationalDoFConstrained) {
        bcColors.push_back(glm::vec3(color_fixedTranslation[0],
                                     color_fixedTranslation[1],
                                     color_fixedTranslation[2]));
      }
    }

    auto bcPolyscopePointCloud =
        polyscope::registerPointCloud(polyscopeLabel_bcAsPointCloud, bcPos);
    bcPolyscopePointCloud->setPointRadius(
        polyscopeHandle_drawOn->getRadius() * 0.05, true);
    bcPolyscopePointCloud->addColorQuantity("bc_colors", bcColors)
        ->setEnabled(true);
  }
#endif  // POLYSCOPE_DEFINED
};
struct SimulationResults {
  /*TODO: remove rotationalDisplacementQuaternion since the last three
   * components of the displacments vector contains the same info using euler
   * angles
   */
  inline const static std::string defaultJsonFilename{"SimulationResults.json"};
  bool converged{false};
  std::shared_ptr<SimulationJob> pJob;
  SimulationHistory history;
  std::vector<Vector6d> debug_drmDisplacements;
  std::vector<Eigen::Quaternion<double>> debug_q_f1;      // per vertex
  std::vector<Eigen::Quaternion<double>> debug_q_normal;  // per vertex
  std::vector<Eigen::Quaternion<double>> debug_q_nr;      // per vertex
  std::vector<Vector6d> displacements;
  std::vector<Eigen::Quaternion<double>>
      rotationalDisplacementQuaternion;  // per vertex
  double internalPotentialEnergy{0};
  double executionTime{0};
  std::vector<std::array<Vector6d, 4>>
      perVertexInternalForces;  // axial,torsion,bending1,bending2
  std::string simulationModelUsed{""};
  std::string labelPrefix{"deformed"};
  inline static char deliminator{' '};
  SimulationResults() { pJob = std::make_shared<SimulationJob>(); }

  std::vector<VectorType> getTranslationalDisplacements() const {
    std::vector<VectorType> translationalDisplacements(displacements.size());
    std::transform(displacements.begin(), displacements.end(),
                   translationalDisplacements.begin(), [&](const Vector6d& d) {
                     return VectorType(d[0], d[1], d[2]);
                   });

    return translationalDisplacements;
  }

  void setLabelPrefix(const std::string& lp) {
    labelPrefix += deliminator + lp;
  }
  std::string getLabel() const {
    return labelPrefix + deliminator + simulationModelUsed + deliminator +
           pJob->pMesh->getLabel() + deliminator + pJob->getLabel();
  }

  bool saveDeformedModel(const std::string& outputFolder = std::string()) {
    VCGEdgeMesh m;
    vcg::tri::Append<VCGEdgeMesh, SimulationEdgeMesh>::MeshCopy(m,
                                                                *pJob->pMesh);
    for (int vi = 0; vi < m.VN(); vi++) {
      m.vert[vi].P() =
          m.vert[vi].P() + CoordType(displacements[vi][0], displacements[vi][1],
                                     displacements[vi][2]);
    }

    // clear faces such that they are not exported when the mesh is saved
    m.face.clear();
    for (VCGEdgeMesh::FaceIterator fi = m.face.begin(); fi != m.face.end();
         ++fi)
      (*fi).Dealloc();
    m.fn = 0;

    return m.save(std::filesystem::path(outputFolder)
                      .append(getLabel() + ".ply")
                      .string());
  }

  //    void saveInternalForces(const std::filesystem::path &outputDirPath)
  //    {
  //        std::cout << "out to:" << outputDirPath << std::endl;
  //        const std::filesystem::path internalForcesDirPath =
  //        std::filesystem::path(outputDirPath);
  //        std::filesystem::create_directories(internalForcesDirPath);
  //        csvFile csv_axial6d(std::filesystem::path(internalForcesDirPath)
  //                                .append("forces_axial_6d.csv"),
  //                            true);
  //        csvFile csv_axialMagn(std::filesystem::path(internalForcesDirPath)
  //                                  .append("forces_axial_magn.csv"),
  //                              true);
  //        csvFile csv_torsion6d(std::filesystem::path(internalForcesDirPath)
  //                                  .append("forces_torsion_6d.csv"),
  //                              true);
  //        csvFile csv_torsionMagn(std::filesystem::path(internalForcesDirPath)
  //                                    .append("forces_torsion_magn.csv"),
  //                                true);
  //        csvFile
  //        csv_firstBending6d(std::filesystem::path(internalForcesDirPath)
  //                                       .append("forces_firstBending_6d.csv"),
  //                                   true);
  //        csvFile
  //        csv_firstBendingMagn(std::filesystem::path(internalForcesDirPath)
  //                                         .append("forces_firstBending_magn.csv"),
  //                                     true);
  //        csvFile
  //        csv_secondBending6d(std::filesystem::path(internalForcesDirPath)
  //                                        .append("forces_secondBending_6d.csv"),
  //                                    true);
  //        csvFile
  //        csv_secondBendingMagn(std::filesystem::path(internalForcesDirPath)
  //                                          .append("forces_secondBending_magn.csv"),
  //                                      true);
  //        for (const std::array<Vector6d, 4> &internalForce :
  //        perVertexInternalForces) {
  //            for (int dofi = 0; dofi < 6; dofi++) {
  //                csv_axial6d << internalForce[0][dofi];
  //                csv_torsion6d << internalForce[1][dofi];
  //                csv_firstBending6d << internalForce[2][dofi];
  //                csv_secondBending6d << internalForce[3][dofi];
  //            }
  //            csv_axial6d << endrow;
  //            csv_torsion6d << endrow;
  //            csv_firstBending6d << endrow;
  //            csv_secondBending6d << endrow;
  //            csv_axialMagn << internalForce[0].norm() << endrow;
  //            csv_torsionMagn << internalForce[1].norm() << endrow;
  //            csv_firstBendingMagn << internalForce[2].norm() << endrow;
  //            csv_secondBendingMagn << internalForce[3].norm() << endrow;
  //        }
  //    }

  void save(const std::string& outputFolder = std::string()) {
    const std::filesystem::path outputFolderPath =
        outputFolder.empty() ? std::filesystem::current_path()
                             : std::filesystem::path(outputFolder);
    //        std::cout << "Saving results to:" << outputFolderPath <<
    //        std::endl;
    std::filesystem::path simulationJobOutputFolderPath =
        std::filesystem::path(outputFolderPath).append("SimulationJob");
    std::filesystem::create_directories(simulationJobOutputFolderPath);
    pJob->save(simulationJobOutputFolderPath.string());
    //    const std::string filename(getLabel() + "_displacements.eigenBin");

    //    Eigen::MatrixXd m = Utilities::toEigenMatrix(displacements);
    const std::filesystem::path resultsFolderPath(
        std::filesystem::path(outputFolder).append("Results"));
    std::filesystem::create_directories(resultsFolderPath);
    //    Eigen::writeBinary(
    //        std::filesystem::path(resultsFolderPath).append(filename).string(),
    //        m);

    //    nlohmann::json json;

    //    json[GET_VARIABLE_NAME(internalPotentialEnergy)] =
    //    internalPotentialEnergy;
    //    // Write internal forces
    //    if (!perVertexInternalForces.empty()) {
    //      std::vector<Vector6d> internalForces_axial(
    //          perVertexInternalForces.size());
    //      std::vector<Vector6d> internalForces_torsion(
    //          perVertexInternalForces.size());
    //      std::vector<Vector6d> internalForces_firstBending(
    //          perVertexInternalForces.size());
    //      std::vector<Vector6d> internalForces_secondBending(
    //          perVertexInternalForces.size());
    //      for (int vi = 0; vi < pJob->pMesh->VN(); vi++) {
    //        internalForces_axial[vi] = perVertexInternalForces[vi][0];
    //        internalForces_torsion[vi] = perVertexInternalForces[vi][1];
    //        internalForces_firstBending[vi] = perVertexInternalForces[vi][2];
    //        internalForces_secondBending[vi] = perVertexInternalForces[vi][3];
    //      }
    //      json[std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
    //      "_axial"] =
    //          internalForces_axial;
    //      json[std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
    //           "_torsion"] = internalForces_torsion;
    //      json[std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
    //           "_firstBending"] = internalForces_firstBending;
    //      json[std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
    //           "_secondBending"] = internalForces_secondBending;
    //    }
    //    std::filesystem::path jsonFilePath(
    //        std::filesystem::path(resultsFolderPath).append(defaultJsonFilename));
    //    std::ofstream jsonFile(jsonFilePath.string());
    //    jsonFile << json;
    //    jsonFile.close();

    saveDeformedModel(resultsFolderPath.string());
  }

  bool isEqual(const SimulationResults& otherSimulationResults, double& error) {
    assert(displacements.size() == otherSimulationResults.displacements.size());
    const double translationalError = [&]() {
      Eigen::MatrixXd thisTranslations =
          Utilities::toEigenMatrix(this->displacements);
      thisTranslations.conservativeResize(pJob->pMesh->VN(), 3);
      Eigen::MatrixXd otherTranslations =
          Utilities::toEigenMatrix(otherSimulationResults.displacements);
      otherTranslations.conservativeResize(pJob->pMesh->VN(), 3);
      return (thisTranslations - otherTranslations).norm();
    }();

    const double rotationalError = [&]() {
      double rotationalError = 0;
      for (int vi = 0; vi < pJob->pMesh->VN(); vi++) {
        rotationalError +=
            std::abs(rotationalDisplacementQuaternion[vi].angularDistance(
                otherSimulationResults.rotationalDisplacementQuaternion[vi]));
        //                const double rotationalError_x =
        //                std::abs(sin(displacements[vi][3] -), cos()); const
        //                double rotationalError_y = displacements[vi][4]; const
        //                double rotationalError_z = displacements[vi][5];
      }
      return rotationalError;
    }();

    error = translationalError + rotationalError;
    return error < 1e-1 ? true : false;
  }

  // The comparison of the results happens comparing the 6-dof nodal
  // displacements
  bool isEqual(const Eigen::MatrixXd& otherMatrix, double& error) {
    assert(otherMatrix.cols() == 6);
    Eigen::MatrixXd thisMatrix = Utilities::toEigenMatrix(this->displacements);
    //        eigenDisplacements.conservativeResize(eigenDisplacements.rows(),
    //        3); Eigen::MatrixXd debug_nodalDisplacements = nodalDisplacements;
    //        debug_nodalDisplacements.conservativeResize(eigenDisplacements.rows(),
    //        3); const double errorNorm = (eigenDisplacements -
    //        debug_nodalDisplacements).norm();
    const Eigen::MatrixXd errorMatrix = (thisMatrix - otherMatrix);
    const double errorNorm = errorMatrix.norm();
    error = errorNorm;
    return errorNorm < 1e-10;
    //    return eigenDisplacements.isApprox(nodalDisplacements);
  }

  void load(const std::filesystem::path& loadFromPath,
            const std::filesystem::path& loadJobFrom) {
    pJob->load(std::filesystem::path(loadJobFrom)
                   .append("SimulationJob.json")
                   .string());
    load(loadFromPath);
  }
  void load(const std::filesystem::path& loadFromPath,
            const std::shared_ptr<SimulationJob>& pJob) {
    this->pJob = pJob;
    load(loadFromPath);
  }

  template <
      typename Container,
      typename T = std::decay_t<decltype(*begin(std::declval<Container>()))>,
      typename = std::enable_if_t<
          std::is_convertible_v<T, std::pair<VertexIndex, VertexIndex>>>>
  static double computeDistance(
      const SimulationResults& resultsA,
      const SimulationResults& resultsB,
      const Container& resultsAToResultsBViMap
      /*,const std::unordered_map<VertexIndex, VertexIndex> */) {
    double distance = 0;
    for (std::pair<int, int> resultsAToResultsBViPair :
         resultsAToResultsBViMap) {
      const double vertexToVertexDistance =
          (resultsA.displacements[resultsAToResultsBViPair.first]
               .getTranslation() -
           resultsB.displacements[resultsAToResultsBViPair.second]
               .getTranslation())
              .norm();
      distance += vertexToVertexDistance;
    }
    return distance;
  }

  double computeDistance(const SimulationResults& other,
                         const std::unordered_map<VertexIndex, VertexIndex>&
                             thisToOtherViMap) const {
    return computeDistance(*this, other, thisToOtherViMap);
  }

#ifdef POLYSCOPE_DEFINED
  void unregister() const {
    if (!polyscope::hasCurveNetwork(getLabel())) {
      std::cerr << "No curve network registered with a name: " << getLabel()
                << std::endl;
      std::cerr << "Nothing to remove." << std::endl;
      return;
    }
    pJob->unregister(getLabel());
    polyscope::removeCurveNetwork(getLabel());
  }

  polyscope::CurveNetwork* registerForDrawing(
      const std::optional<std::array<float, 3>>& desiredColor = std::nullopt,
      const bool& shouldEnable = true,
      const std::optional<float>& inputDrawingRadius = std::nullopt,
      const bool& shouldShowFrames = false) const {
    if (!converged) {
      std::cerr << "Simulation did not converge. Drawing nothing." << std::endl;
      return nullptr;
    }
    PolyscopeInterface::init();
    const std::shared_ptr<SimulationEdgeMesh>& simulationEdgeMesh = pJob->pMesh;
    polyscope::CurveNetwork* polyscopeHandle_deformedEdmeMesh;
    const std::string& meshLabel = getLabel();
    if (!polyscope::hasStructure(meshLabel)) {
      const auto& verts = simulationEdgeMesh->getEigenVertices();
      const auto& edges = simulationEdgeMesh->getEigenEdges();
      polyscopeHandle_deformedEdmeMesh =
          polyscope::registerCurveNetwork(meshLabel, verts, edges);
      //         simulationEdgeMesh->registerForDrawing()
    } else {
      polyscopeHandle_deformedEdmeMesh = polyscope::getCurveNetwork(meshLabel);
    }
    polyscopeHandle_deformedEdmeMesh->setEnabled(shouldEnable);
    const float drawingRadius = [&]() {
      if (inputDrawingRadius.has_value()) {
        return inputDrawingRadius.value();
      } else {
        return simulationEdgeMesh->getBeamDimensions()[0].getDrawingRadius();
      }
    }();
    polyscopeHandle_deformedEdmeMesh->setRadius(drawingRadius /*0.00025*/,
                                                false);
    if (desiredColor.has_value()) {
      const glm::vec3 desiredColor_glm(desiredColor.value()[0],
                                       desiredColor.value()[1],
                                       desiredColor.value()[2]);
      polyscopeHandle_deformedEdmeMesh->setColor(desiredColor_glm);
    }
    Eigen::MatrixX3d nodalDisplacements(simulationEdgeMesh->VN(), 3);
    Eigen::MatrixX3d framesX(simulationEdgeMesh->VN(), 3);
    Eigen::MatrixX3d framesY(simulationEdgeMesh->VN(), 3);
    Eigen::MatrixX3d framesZ(simulationEdgeMesh->VN(), 3);

    Eigen::MatrixX3d framesX_initial(simulationEdgeMesh->VN(), 3);
    Eigen::MatrixX3d framesY_initial(simulationEdgeMesh->VN(), 3);
    Eigen::MatrixX3d framesZ_initial(simulationEdgeMesh->VN(), 3);

    //        std::unordered_set<int> interfaceNodes{1, 3, 5, 7, 9, 11};
    std::unordered_set<int> interfaceNodes{3, 7, 11, 15, 19, 23};
    //      std::unordered_set<int> interfaceNodes{};
    for (VertexIndex vi = 0; vi < simulationEdgeMesh->VN(); vi++) {
      const Vector6d& nodalDisplacement = displacements[vi];
      nodalDisplacements.row(vi) = Eigen::Vector3d(
          nodalDisplacement[0], nodalDisplacement[1], nodalDisplacement[2]);
      //           Eigen::Quaternion<double>
      //           Rx(Eigen::AngleAxis(nodalDisplacement[2],Eigen::Vector3d(1,
      //           0, 0))); Eigen::Quaternion<double>
      //           Ry(Eigen::AngleAxis(nodalDisplacement[4],Eigen::Vector3d(0,
      //           1, 0))); Eigen::Quaternion<double>
      //           Rz(Eigen::AngleAxis(nodalDisplacement[5],Eigen::Vector3d(0,
      //           0, 1))); Eigen::Quaternion<double> R=Rx*Ry*Rz;
      //            if (interfaceNodes.contains(vi)) {
      assert(!rotationalDisplacementQuaternion.empty());
      auto deformedNormal =
          rotationalDisplacementQuaternion[vi] * Eigen::Vector3d(0, 0, 1);
      auto deformedFrameY =
          rotationalDisplacementQuaternion[vi] * Eigen::Vector3d(0, 1, 0);
      auto deformedFrameX =
          rotationalDisplacementQuaternion[vi] * Eigen::Vector3d(1, 0, 0);
      framesX.row(vi) = Eigen::Vector3d(deformedFrameX[0], deformedFrameX[1],
                                        deformedFrameX[2]);
      framesY.row(vi) = Eigen::Vector3d(deformedFrameY[0], deformedFrameY[1],
                                        deformedFrameY[2]);
      framesZ.row(vi) = Eigen::Vector3d(deformedNormal[0], deformedNormal[1],
                                        deformedNormal[2]);
      framesX_initial.row(vi) = Eigen::Vector3d(1, 0, 0);
      framesY_initial.row(vi) = Eigen::Vector3d(0, 1, 0);
      framesZ_initial.row(vi) = Eigen::Vector3d(0, 0, 1);
      //            } else {
      //                framesX.row(vi) = Eigen::Vector3d(0, 0, 0);
      //                framesY.row(vi) = Eigen::Vector3d(0, 0, 0);
      //                framesZ.row(vi) = Eigen::Vector3d(0, 0, 0);
      //                framesX_initial.row(vi) = Eigen::Vector3d(0, 0, 0);
      //                framesY_initial.row(vi) = Eigen::Vector3d(0, 0, 0);
      //                framesZ_initial.row(vi) = Eigen::Vector3d(0, 0, 0);
      //            }
    }
    polyscopeHandle_deformedEdmeMesh->updateNodePositions(
        simulationEdgeMesh->getEigenVertices() + nodalDisplacements);

    const double frameRadius_default = 0.035;
    const double frameLength_default = 0.035;
    const bool shouldEnable_default = true;
    // if(showFramesOn.contains(vi)){
    auto polyscopeHandle_frameX =
        polyscopeHandle_deformedEdmeMesh->addNodeVectorQuantity("FrameX",
                                                                framesX);
    polyscopeHandle_frameX->setVectorLengthScale(frameLength_default);
    polyscopeHandle_frameX->setVectorRadius(frameRadius_default);
    polyscopeHandle_frameX->setVectorColor(
        /*polyscope::getNextUniqueColor()*/ glm::vec3(1, 0, 0));
    auto polyscopeHandle_frameY =
        polyscopeHandle_deformedEdmeMesh->addNodeVectorQuantity("FrameY",
                                                                framesY);
    polyscopeHandle_frameY->setVectorLengthScale(frameLength_default);
    polyscopeHandle_frameY->setVectorRadius(frameRadius_default);
    polyscopeHandle_frameY->setVectorColor(
        /*polyscope::getNextUniqueColor()*/ glm::vec3(0, 1, 0));
    auto polyscopeHandle_frameZ =
        polyscopeHandle_deformedEdmeMesh->addNodeVectorQuantity("FrameZ",
                                                                framesZ);
    polyscopeHandle_frameZ->setVectorLengthScale(frameLength_default);
    polyscopeHandle_frameZ->setVectorRadius(frameRadius_default);
    polyscopeHandle_frameZ->setVectorColor(
        /*polyscope::getNextUniqueColor()*/ glm::vec3(0, 0, 1));

    //        if (!polyscope::hasCurveNetwork(mesh->getLabel())) {
    //            const std::array<double, 3> initialColor({0, 0, 0});
    //            /*auto polyscopeHandle_initialMesh
    //            =*/mesh->registerForDrawing(initialColor);
    //        }
    pJob->registerForDrawing(polyscopeHandle_deformedEdmeMesh->name);

    //      auto polyscopeHandle_frameX_initial = polyscopeHandle_initialMesh
    //                                                ->addNodeVectorQuantity("FrameX",
    //                                                framesX_initial);
    //      polyscopeHandle_frameX_initial->setVectorLengthScale(frameLength_default);
    //      polyscopeHandle_frameX_initial->setVectorRadius(frameRadius_default);
    //      polyscopeHandle_frameX_initial->setVectorColor(
    //          /*polyscope::getNextUniqueColor()*/ glm::vec3(1, 0, 0));
    //      auto polyscopeHandle_frameY_initial = polyscopeHandle_initialMesh
    //                                                ->addNodeVectorQuantity("FrameY",
    //                                                framesY_initial);
    //      polyscopeHandle_frameY_initial->setVectorLengthScale(frameLength_default);
    //      polyscopeHandle_frameY_initial->setVectorRadius(frameRadius_default);
    //      polyscopeHandle_frameY_initial->setVectorColor(
    //          /*polyscope::getNextUniqueColor()*/ glm::vec3(0, 1, 0));
    //      auto polyscopeHandle_frameZ_initial = polyscopeHandle_initialMesh
    //                                                ->addNodeVectorQuantity("FrameZ",
    //                                                framesZ_initial);
    //      polyscopeHandle_frameZ_initial->setVectorLengthScale(frameLength_default);
    //      polyscopeHandle_frameZ_initial->setVectorRadius(frameRadius_default);
    //      polyscopeHandle_frameZ_initial->setVectorColor(
    //          /*polyscope::getNextUniqueColor()*/ glm::vec3(0, 0, 1));
    //      //}
    //    pJob->registerForDrawing(, false);
    //        static bool wasExecuted =false;
    //        if (!wasExecuted) {
    //            std::function<void()> callback = [&]() {
    //                static bool showFrames = shouldShowFrames;

    //                if (ImGui::Checkbox("Show Frames", &showFrames) &&
    //                showFrames) {
    polyscopeHandle_frameX->setEnabled(shouldShowFrames);
    polyscopeHandle_frameY->setEnabled(shouldShowFrames);
    polyscopeHandle_frameZ->setEnabled(shouldShowFrames);
    //            };

    //            PolyscopeInterface::addUserCallback(callback);
    //            wasExecuted = true;
    //        }
    return polyscopeHandle_deformedEdmeMesh;
  }
#endif
 private:
  bool load(const std::filesystem::path& loadFromPath) {
    converged = true;  // assuming it has converged
    assert(pJob != nullptr);
    // load job
    // Use the first .eigenBin file for loading the displacements
    for (auto const& entry :
         std::filesystem::recursive_directory_iterator(loadFromPath)) {
      if (std::filesystem::is_regular_file(entry) &&
          entry.path().extension() == ".eigenBin") {
        Eigen::MatrixXd displacements_eigen;
        Eigen::readBinary(entry.path().string(), displacements_eigen);
        displacements = Utilities::fromEigenMatrix(displacements_eigen);
        break;
      }
    }

    rotationalDisplacementQuaternion.resize(pJob->pMesh->VN());
    for (int vi = 0; vi < pJob->pMesh->VN(); vi++) {
      rotationalDisplacementQuaternion[vi] =
          Eigen::AngleAxisd(displacements[vi][3], Eigen::Vector3d::UnitX()) *
          Eigen::AngleAxisd(displacements[vi][4], Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(displacements[vi][5], Eigen::Vector3d::UnitZ());
    }

    const std::filesystem::path jsonFilepath =
        std::filesystem::path(loadFromPath).append(defaultJsonFilename);
    if (std::filesystem::exists(jsonFilepath)) {
      std::ifstream ifs(jsonFilepath);
      nlohmann::json json;
      ifs >> json;
      //        if (json.contains(GET_VARIABLE_NAME(internalPotentialEnergy))) {
      //            internalPotentialEnergy =
      //            json.at(GET_VARIABLE_NAME(internalPotentialEnergy));
      //        }
      if (json.contains(
              std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
              "_axial")) {
        perVertexInternalForces.resize(pJob->pMesh->VN());
        std::vector<Vector6d> perVertexInternalForces_axial =
            static_cast<std::vector<Vector6d>>(json.at(
                std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
                "_axial"));
        for (int vi = 0; vi < pJob->pMesh->VN(); vi++) {
          perVertexInternalForces[vi][0] = perVertexInternalForces_axial[vi];
        }
      }
      if (json.contains(
              std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
              "_torsion")) {
        perVertexInternalForces.resize(pJob->pMesh->VN());
        std::vector<Vector6d> perVertexInternalForces_axial =
            static_cast<std::vector<Vector6d>>(json.at(
                std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
                "_torsion"));
        for (int vi = 0; vi < pJob->pMesh->VN(); vi++) {
          perVertexInternalForces[vi][0] = perVertexInternalForces_axial[vi];
        }
      }
      if (json.contains(
              std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
              "_firstBending")) {
        perVertexInternalForces.resize(pJob->pMesh->VN());
        std::vector<Vector6d> perVertexInternalForces_axial =
            static_cast<std::vector<Vector6d>>(json.at(
                std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
                "_firstBending"));
        for (int vi = 0; vi < pJob->pMesh->VN(); vi++) {
          perVertexInternalForces[vi][0] = perVertexInternalForces_axial[vi];
        }
      }
      if (json.contains(
              std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
              "_secondBending")) {
        perVertexInternalForces.resize(pJob->pMesh->VN());
        std::vector<Vector6d> perVertexInternalForces_axial =
            static_cast<std::vector<Vector6d>>(json.at(
                std::string(GET_VARIABLE_NAME(perVertexInternalForces)) +
                "_secondBending"));
        for (int vi = 0; vi < pJob->pMesh->VN(); vi++) {
          perVertexInternalForces[vi][0] = perVertexInternalForces_axial[vi];
        }
      }
    }
    return true;
  }
};

#endif  // SIMULATIONHISTORY_HPP
