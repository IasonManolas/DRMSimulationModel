#ifndef SIMULATIONHISTORYPLOTTER_HPP
#define SIMULATIONHISTORYPLOTTER_HPP

#include <matplot/matplot.h>
#include <algorithm>
#include "simulation_structs.hpp"
#include "simulationmesh.hpp"
#include "utilities.hpp"

struct SimulationResultsReporter {
  using VertexType = VCGEdgeMesh::VertexType;
  using CoordType = VCGEdgeMesh::CoordType;

  SimulationResultsReporter() {}

  void writeStatistics(const SimulationResults& results,
                       const std::string& reportFolderPath) {
    std::ofstream file;
    file.open(
        std::filesystem::path(reportFolderPath).append("results.txt").string());

    const size_t numberOfSteps = results.history.numberOfSteps;
    file << "Number of steps " << numberOfSteps << "\n";

    //    file << "Force threshold used " << 1000 << "\n";

    //    assert(numberOfSteps == results.history.potentialEnergy.size() &&
    //           numberOfSteps == results.history.residualForces.size());
    // Write kinetic energies
    const SimulationHistory& history = results.history;
    if (!history.kineticEnergy.empty()) {
      file << "Kinetic energies"
           << "\n";
      for (size_t step = 0; step < numberOfSteps; step++) {
        file << history.kineticEnergy[step] << "\n";
      }
      file << "\n";
    }

    if (!history.logResidualForces.empty()) {
      file << "Residual forces"
           << "\n";
      for (size_t step = 0; step < numberOfSteps; step++) {
        file << history.logResidualForces[step] << "\n";
      }
      file << "\n";
    }

    if (!history.potentialEnergies.empty()) {
      file << "Potential energies"
           << "\n";
      for (size_t step = 0; step < numberOfSteps; step++) {
        file << history.potentialEnergies[step] << "\n";
      }
      file << "\n";
    }
    file.close();
  }

  void reportHistory(const SimulationHistory& history,
                     const std::string& reportFolderPath,
                     const std::string& graphSuffix = std::string()) {
    const auto simulationResultPath =
        std::filesystem::path(reportFolderPath).append(history.label);
    std::filesystem::create_directories(simulationResultPath);
    createPlots(history, simulationResultPath.string(), graphSuffix);
  }

  void reportResults(const std::vector<SimulationResults>& results,
                     const std::string& reportFolderPath,
                     const std::string& graphSuffix = std::string()) {
    if (results.empty()) {
      return;
    }

    //    std::filesystem::remove_all(debuggerFolder);
    std::filesystem::create_directory(reportFolderPath);
    for (const SimulationResults& simulationResult : results) {
      const auto simulationResultPath =
          std::filesystem::path(reportFolderPath)
              .append(simulationResult.getLabel());
      std::filesystem::create_directory(simulationResultPath.string());

      createPlots(simulationResult.history, simulationResultPath.string(),
                  graphSuffix);
      writeStatistics(simulationResult, simulationResultPath.string());
    }
  }

  static void createPlot(const std::string& xLabel,
                         const std::string& yLabel,
                         const std::vector<double>& YvaluesToPlot,
                         const std::string& saveTo = {},
                         const std::vector<size_t>& markPoints = {}) {
    if (YvaluesToPlot.size() < 2) {
      return;
    }
    std::vector<double> colors(YvaluesToPlot.size(), 0.5);
    std::vector<double> markerSizes(YvaluesToPlot.size(), 5);
    if (!markPoints.empty()) {
      for (const auto pointIndex : markPoints) {
        colors[pointIndex] = 0.9;
        markerSizes[pointIndex] = 10;
      }
    }
    std::vector<double> x =
        matplot::linspace(0, YvaluesToPlot.size() - 1, YvaluesToPlot.size());
    Utilities::createPlot(xLabel, yLabel, x, YvaluesToPlot, markerSizes, colors,
                          saveTo);
  }

  void createPlots(const SimulationHistory& history,
                   const std::string& reportFolderPath,
                   const std::string& graphSuffix) {
    const auto graphsFolder =
        std::filesystem::path(reportFolderPath).append("Graphs");
    std::filesystem::remove_all(graphsFolder);
    std::filesystem::create_directory(graphsFolder.string());

    if (!history.kineticEnergy.empty()) {
      createPlot("Number of Iterations", "Log of Kinetic Energy log",
                 history.kineticEnergy,
                 std::filesystem::path(graphsFolder)
                     .append("KineticEnergyLog_" + graphSuffix + ".png")
                     .string(),
                 history.redMarks);
    }

    if (!history.logResidualForces.empty()) {
      createPlot("Number of Iterations", "Residual Forces norm log",
                 history.logResidualForces,
                 std::filesystem::path(graphsFolder)
                     .append("ResidualForcesLog_" + graphSuffix + ".png")
                     .string(),
                 history.redMarks);
    }

    if (!history.potentialEnergies.empty()) {
      createPlot("Number of Iterations", "Potential energy",
                 history.potentialEnergies,
                 std::filesystem::path(graphsFolder)
                     .append("PotentialEnergy_" + graphSuffix + ".png")
                     .string(),
                 history.redMarks);
    }
  }
};

#endif  // SIMULATIONHISTORYPLOTTER_HPP
