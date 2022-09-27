#include "utilities.hpp"
#include "matplot/core/figure_registry.h"
#include "matplot/matplot.h"

void Utilities::createPlot(const std::string& xLabel,
                           const std::string& yLabel,
                           const std::vector<double>& x,
                           const std::vector<double>& y,
                           const std::vector<double>& markerSizes,
                           const std::vector<double>& c,
                           const std::string& saveTo) {
  matplot::xlabel(xLabel);
  matplot::ylabel(yLabel);
  matplot::grid(matplot::on);
  auto matplotHandle_scatter = matplot::scatter(x, y, markerSizes, c);
  auto fig = matplot::gcf();
  matplotHandle_scatter->marker_face(true);
  fig->font_size(20.0);
  if (!saveTo.empty()) {
    matplot::save(saveTo);
  }
}
#ifdef POLYSCOPE_DEFINED
#include <functional>
#include "polyscope/curve_network.h"
#include "polyscope/pick.h"
#include "polyscope/polyscope.h"

void PolyscopeInterface::mainCallback() {
  ImGui::PushItemWidth(100);
  for (std::function<void()>& userCallback :
       globalPolyscopeData.userCallbacks) {
    userCallback();
  }

  ImGui::PopItemWidth();
}

void PolyscopeInterface::addUserCallback(
    const std::function<void()>& userCallback) {
  globalPolyscopeData.userCallbacks.push_back(userCallback);
}

void PolyscopeInterface::deinitPolyscope() {
  if (!polyscope::state::initialized) {
    return;
  }

  polyscope::render::engine->shutdownImGui();
}

void PolyscopeInterface::init() {
  if (polyscope::state::initialized) {
    return;
  }
  polyscope::init();
  polyscope::options::groundPlaneEnabled = false;
  polyscope::view::upDir = polyscope::view::UpDir::ZUp;

  polyscope::state::userCallback = &mainCallback;
  polyscope::options::autocenterStructures = false;
  polyscope::options::autoscaleStructures = false;
}

std::pair<PolyscopeInterface::PolyscopeLabel, size_t>
PolyscopeInterface::getSelection() {
  std::pair<polyscope::Structure*, size_t> selection =
      polyscope::pick::getSelection();
  if (selection.first == nullptr) {
    return std::make_pair(std::string(), 0);
  }
  return std::make_pair(selection.first->name, selection.second);
}

void PolyscopeInterface::registerWorldAxes() {
  PolyscopeInterface::init();

  Eigen::MatrixX3d axesPositions(4, 3);
  axesPositions.row(0) = Eigen::Vector3d(0, 0, 0);
  axesPositions.row(1) = Eigen::Vector3d(1, 0, 0);
  axesPositions.row(2) = Eigen::Vector3d(0, 1, 0);
  axesPositions.row(3) = Eigen::Vector3d(0, 0, 1);

  Eigen::MatrixX2i axesEdges(3, 2);
  axesEdges.row(0) = Eigen::Vector2i(0, 1);
  axesEdges.row(1) = Eigen::Vector2i(0, 2);
  axesEdges.row(2) = Eigen::Vector2i(0, 3);
  Eigen::MatrixX3d axesColors(3, 3);
  axesColors.row(0) = Eigen::Vector3d(1, 0, 0);
  axesColors.row(1) = Eigen::Vector3d(0, 1, 0);
  axesColors.row(2) = Eigen::Vector3d(0, 0, 1);

  const std::string worldAxesName = "World Axes";
  polyscope::registerCurveNetwork(worldAxesName, axesPositions, axesEdges);
  polyscope::getCurveNetwork(worldAxesName)->setRadius(0.0001, false);
  const std::string worldAxesColorName = worldAxesName + " Color";
  polyscope::getCurveNetwork(worldAxesName)
      ->addEdgeColorQuantity(worldAxesColorName, axesColors)
      ->setEnabled(true);
}

#endif
