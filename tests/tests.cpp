#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include "drmsimulationmodel.hpp"
#include "simulation_structs.hpp"

static SimulationResults simulate(
    const std::filesystem::path& folderPath_simulationJob) {
  if (!std::filesystem::exists(folderPath_simulationJob)) {
    std::cerr << "The input simulation job folder path does not exist:"
              << folderPath_simulationJob << std::endl;
    return SimulationResults();
  }
  std::shared_ptr<SimulationJob> pJob = std::make_shared<SimulationJob>(
      std::filesystem::path(folderPath_simulationJob)
          .append(SimulationJob::jsonDefaultFileName));
  DRMSimulationModel::Settings settings;
  settings.load(std::filesystem::path(folderPath_simulationJob)
                    .append(DRMSimulationModel::Settings::jsonDefaultFileName));
  settings.beVerbose = false;
  DRMSimulationModel simulationModel_drm;
  SimulationResults simulationResults_drm =
      simulationModel_drm.executeSimulation(pJob, settings);

  if (!simulationResults_drm.converged) {
    std::cerr << "Simulation did not converge." << std::endl;
    return SimulationResults();
  }
  return simulationResults_drm;
}

constexpr double accuracyThreshold = 1e-5;

TEST_CASE("Simulation results of cantilever", "[single-file]") {
  const std::filesystem::path folderPath_scenario(
      "/home/iason/Coding/Libraries/DRMSimualtionModel_github/"
      "demo/PaperResults/cantilever");
  const std::filesystem::path folderPath_scenarioSimulationJob(
      std::filesystem::path(folderPath_scenario).append("SimulationJob"));
  const std::filesystem::path folderPath_scenarioGroundTruthResults(
      std::filesystem::path(folderPath_scenario).append("Results"));
  SimulationResults simulationResults_groundTruth;
  simulationResults_groundTruth.load(folderPath_scenarioGroundTruthResults,
                                     folderPath_scenarioSimulationJob);
  SimulationResults simulationResults_test =
      simulate(folderPath_scenarioSimulationJob);

  REQUIRE(simulationResults_test.converged);
  REQUIRE(simulationResults_groundTruth.absDifference(simulationResults_test) <
          accuracyThreshold);
}

TEST_CASE("Simulation results of torus", "[single-file]") {
  const std::filesystem::path folderPath_scenario(
      "/home/iason/Coding/Libraries/DRMSimualtionModel_github/"
      "demo/PaperResults/torus");
  const std::filesystem::path folderPath_scenarioSimulationJob(
      std::filesystem::path(folderPath_scenario).append("SimulationJob"));
  const std::filesystem::path folderPath_scenarioGroundTruthResults(
      std::filesystem::path(folderPath_scenario).append("Results"));
  SimulationResults simulationResults_groundTruth;
  simulationResults_groundTruth.load(folderPath_scenarioGroundTruthResults,
                                     folderPath_scenarioSimulationJob);
  SimulationResults simulationResults_test =
      simulate(folderPath_scenarioSimulationJob);

  REQUIRE(simulationResults_test.converged);
  REQUIRE(simulationResults_groundTruth.absDifference(simulationResults_test) <
          accuracyThreshold);
}
