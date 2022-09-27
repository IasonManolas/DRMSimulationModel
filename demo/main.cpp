#include <memory>
#include "drmsimulationmodel.hpp"
#include "simulation_structs.hpp"

int main(int argc, char** argv) {
  constexpr int expectedNumberOfArgs = 1;
  if (argc != expectedNumberOfArgs + 1) {
    std::cerr << "Wrong number of input arguments." << expectedNumberOfArgs
              << " are expected." << std::endl;
    const std::string usageString =
        "Input arguments:\n"
        "1)Simulation job folder path\n"
        "The tool assumes that a .json file containing the drm settings is "
        "found in the same folder as the simulation job. If such a file is not "
        "found the default drm settings are used."
        "Exiting..";
    std::cout << usageString << std::endl;
    std::cerr << "Input arguments are:" << std::endl;
    std::copy(argv + 1, argv + argc,
              std::ostream_iterator<const char*>(std::cout, "\n"));
    return 1;
  }

  const std::filesystem::path folderPath_simulationJob(argv[1]);
  if (!std::filesystem::exists(folderPath_simulationJob)) {
    std::cerr << "The input simulation job folder path does not exist:"
              << folderPath_simulationJob << std::endl;
    return 1;
  }
  std::shared_ptr<SimulationJob> pJob = std::make_shared<SimulationJob>(
      std::filesystem::path(folderPath_simulationJob)
          .append(SimulationJob::jsonDefaultFileName));
  DRMSimulationModel::Settings settings;
  settings.load(std::filesystem::path(folderPath_simulationJob)
                    .append(DRMSimulationModel::Settings::jsonDefaultFileName));
  settings.beVerbose = true;
  DRMSimulationModel simulationModel_drm;
  SimulationResults simulationResults_drm =
      simulationModel_drm.executeSimulation(pJob, settings);

  if (!simulationResults_drm.converged) {
    return 1;
  }

  const std::string folderName_results = "ComputedSimulationResults";
  const std::filesystem::path resultsFolderPath_these =
      std::filesystem::current_path().append(folderName_results);
  std::cout << "Saving results to:" << resultsFolderPath_these << std::endl;
  std::filesystem::remove_all(resultsFolderPath_these);
  std::filesystem::create_directories(resultsFolderPath_these);
  simulationResults_drm.save(resultsFolderPath_these);

#ifdef POLYSCOPE_DEFINED
  //  settings.debugModeStep = 1e4;
  RGBColor color_deformed{67.0 / 255, 160.00 / 255, 232.0 / 255};
  simulationResults_drm.registerForDrawing(color_deformed);
  RGBColor color_initial{24.0 / 255, 23.0 / 255, 23.0 / 255};
  pJob->pMesh->registerForDrawing(color_initial);
  polyscope::show();
#endif

  return 0;
}
