#include "drmsimulationmodel.hpp"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <algorithm>
#include <chrono>
#include <execution>
#include <stdexcept>
#include "simulationhistoryplotter.hpp"
#include "simulationmesh.hpp"

#define RUN_OPT

void DRMSimulationModel::reset(const std::shared_ptr<SimulationJob>& pJob) {
  pMesh.reset();
  pMesh = std::make_unique<SimulationEdgeMesh>(*pJob->pMesh);
  vcg::tri::UpdateBounding<SimulationEdgeMesh>::Box(*pMesh);

  constrainedVertices.clear();
  constrainedVertices = pJob->constrainedVertices;
  for (const std::pair<VertexIndex, Eigen::Vector3d>
           vertexIndexDisplacementPair : pJob->nodalForcedDisplacements) {
    constrainedVertices[vertexIndexDisplacementPair.first].insert(
        {DoF::Ux, DoF::Uy, DoF::Uz});
  }

  computeRigidSupports();
  isVertexConstrained.resize(pMesh->VN(), false);
  for (auto fixedVertex : constrainedVertices) {
    isVertexConstrained[fixedVertex.first] = true;
  }
}

void DRMSimulationModel::reset(const std::shared_ptr<SimulationJob>& pJob,
                               const Settings& settings) {
  mSettings = settings;
  mCurrentSimulationStep = 0;
  if (settings.shouldCreatePlots) {
    history.clear();
    history.label = pJob->pMesh->getLabel() + "_" + pJob->getLabel();
    plotYValues.clear();
  }
  checkedForMaximumMoment = false;
  externalMomentsNorm = 0;
  Dt = settings.Dtini;
  numOfDampings = 0;
  isVertexConstrained.clear();
  postConvergenceMessage.clear();

  reset(pJob);

#ifdef POLYSCOPE_DEFINED
  if (mSettings.shouldDraw) {
    pJob->pMesh->registerForDrawing(color_initial);
    pMesh->prependToLabel(polyscopeLabel_deformed + "_");
    pMesh->registerForDrawing(color_deformed);
  }
#endif

  if (!pJob->nodalForcedDisplacements.empty() &&
      pJob->nodalExternalForces.empty()) {
    if (!mSettings.initialDistortion.has_value()) {
      std::cerr << "An initial out of plane distortion should be applied due "
                   "to bifurcation buckling"
                << std::endl;
    }
  }
  updateElementalFrames();
}

VectorType DRMSimulationModel::computeDisplacementDifferenceDerivative(
    const EdgeType& e,
    const DifferentiateWithRespectTo& dui) const {
  VectorType displacementDiffDeriv(0, 0, 0);
  const DoFType& dofi = dui.dofi;
  const bool differentiateWithRespectToANonEdgeNode =
      e.cV(0) != &dui.v && e.cV(1) != &dui.v;
  if (differentiateWithRespectToANonEdgeNode || dofi > 2) {
    return displacementDiffDeriv;
  }

  if (e.cV(0) == &dui.v) {
    displacementDiffDeriv[dofi] = -1;
  } else if (e.cV(1) == &dui.v) {
    displacementDiffDeriv[dofi] = 1;
  }

  return displacementDiffDeriv;
}

VectorType DRMSimulationModel::computeDerivativeOfNormal(
    const VertexType& v,
    const DifferentiateWithRespectTo& dui) const {
  const size_t vi = pMesh->getIndex(v);
  VectorType normalDerivative(0, 0, 0);
  if (&dui.v != &v ||
      (dui.dofi == 0 || dui.dofi == 1 || dui.dofi == 2 || dui.dofi == 5)) {
    return normalDerivative;
  }
  const VectorType& n = v.cN();
  const double &nx = n[0], ny = n[1];
  const double nxnyMagnitude = std::pow(nx, 2) + std::pow(ny, 2);

  if (dui.dofi == 3) {
    if (nxnyMagnitude + 1e-5 >= 1) {
      const double normalDerivativeX =
          1 / sqrt(nxnyMagnitude) -
          std::pow(nx, 2) / std::pow(nxnyMagnitude, 1.5);
      const double normalDerivativeY = -nx * ny / std::pow(nxnyMagnitude, 1.5);
      const double normalDerivativeZ = 0;
      normalDerivative =
          VectorType(normalDerivativeX, normalDerivativeY, normalDerivativeZ);
    } else {
      const double normalDerivativeX = 1;
      const double normalDerivativeY = 0;
      const double normalDerivativeZ = -nx / std::sqrt(1 - nxnyMagnitude);
      normalDerivative =
          VectorType(normalDerivativeX, normalDerivativeY, normalDerivativeZ);
    }
  } else if (dui.dofi == 4) {
    if (nxnyMagnitude + 1e-5 >= 1) {
      const double normalDerivativeX = -nx * ny / std::pow(nxnyMagnitude, 1.5);
      const double normalDerivativeY =
          1 / sqrt(nxnyMagnitude) -
          std::pow(ny, 2) / std::pow(nxnyMagnitude, 1.5);
      const double normalDerivativeZ = 0;
      normalDerivative =
          VectorType(normalDerivativeX, normalDerivativeY, normalDerivativeZ);
    } else {
      const double normalDerivativeX = 0;
      const double normalDerivativeY = 1;
      const double normalDerivativeZ = -ny / std::sqrt(1 - nxnyMagnitude);
      normalDerivative =
          VectorType(normalDerivativeX, normalDerivativeY, normalDerivativeZ);
    }
  }

  return normalDerivative;
}

double DRMSimulationModel::computeDerivativeElementLength(
    const EdgeType& e,
    const DifferentiateWithRespectTo& dui) const {
  if (e.cV(0) != &dui.v && e.cV(1) != &dui.v) {
    return 0;
  }

  const VectorType& X_j = e.cP(0);
  const VectorType& X_jplus1 = e.cP(1);
  const VectorType positionVectorDiff = X_jplus1 - X_j;
  const VectorType displacementDiffDeriv =
      computeDisplacementDifferenceDerivative(e, dui);
  const double edgeLength = pMesh->elements[e].length;
  const double L_kDeriv =
      positionVectorDiff * displacementDiffDeriv / edgeLength;
  return L_kDeriv;
}

double DRMSimulationModel::computeDerivativeOfNorm(
    const VectorType& x,
    const VectorType& derivativeOfX) {
  return x.dot(derivativeOfX) / x.Norm();
}

VectorType DRMSimulationModel::computeDerivativeOfCrossProduct(
    const VectorType& a,
    const VectorType& derivativeOfA,
    const VectorType& b,
    const VectorType& derivativeOfB) {
  const auto firstTerm = Cross(derivativeOfA, b);
  const auto secondTerm = Cross(a, derivativeOfB);
  return firstTerm + secondTerm;
}

VectorType DRMSimulationModel::computeDerivativeOfR(
    const EdgeType& e,
    const DifferentiateWithRespectTo& dui) const {
  const VertexType& v_j = *e.cV(0);
  const VertexType& v_jplus1 = *e.cV(1);
  const VectorType normal_j = v_j.cN();
  const VectorType normal_jplus1 = v_jplus1.cN();
  const VectorType derivativeOfNormal_j =
      &v_j == &dui.v && dui.dofi > 2
          ? pMesh->nodes[v_j].derivativeOfNormal[dui.dofi]
          : VectorType(0, 0, 0);
  const VectorType derivativeOfNormal_jplus1 =
      &v_jplus1 == &dui.v && dui.dofi > 2
          ? pMesh->nodes[v_jplus1].derivativeOfNormal[dui.dofi]
          : VectorType(0, 0, 0);

  const VectorType derivativeOfSumOfNormals =
      derivativeOfNormal_j + derivativeOfNormal_jplus1;
  const VectorType sumOfNormals = normal_j + normal_jplus1;
  const double normOfSumOfNormals = sumOfNormals.Norm();
  assert(normOfSumOfNormals != 0);
  const double derivativeOfNormOfSumOfNormals =
      computeDerivativeOfNorm(sumOfNormals, derivativeOfSumOfNormals);

  const VectorType derivativeOfR_firstTerm = -sumOfNormals *
                                             derivativeOfNormOfSumOfNormals /
                                             std::pow(normOfSumOfNormals, 2);
  const VectorType derivativeOfR_secondTerm =
      derivativeOfSumOfNormals / normOfSumOfNormals;
  const VectorType derivativeOfR =
      derivativeOfR_firstTerm + derivativeOfR_secondTerm;

  return derivativeOfR;
}

VectorType DRMSimulationModel::computeDerivativeT1(
    const EdgeType& e,
    const DifferentiateWithRespectTo& dui) const {
  const VectorType& X_j = e.cP(0);
  const VectorType& X_jplus1 = e.cP(1);
  const VectorType edgeVector = X_jplus1 - X_j;
  const double L_kDerivative = computeDerivativeElementLength(e, dui);
  const double edgeLength = pMesh->elements[e].length;
  const VectorType firstTerm =
      -edgeVector * L_kDerivative / std::pow(edgeLength, 2);
  const VectorType secondTerm =
      computeDisplacementDifferenceDerivative(e, dui) / edgeLength;
  const VectorType t1Derivative = firstTerm + secondTerm;

  return t1Derivative;
}

VectorType DRMSimulationModel::computeDerivativeT2(
    const EdgeType& e,
    const DifferentiateWithRespectTo& dui) const {
  const DoFType dofi = dui.dofi;

  const VertexType& v_j = *e.cV(0);
  const VertexType& v_jplus1 = *e.cV(1);

  const VectorType r = (v_j.cN() + v_jplus1.cN()).Normalize();
  const VectorType derivativeR_j =
      dofi > 2 && &v_j == &dui.v ? pMesh->elements[e].derivativeR[0][dui.dofi]
                                 : VectorType(0, 0, 0);
  const VectorType derivativeR_jplus1 =
      dofi > 2 && &v_jplus1 == &dui.v
          ? pMesh->elements[e].derivativeR[1][dui.dofi]
          : VectorType(0, 0, 0);
  const VectorType derivativeR = derivativeR_j + derivativeR_jplus1;

  const VectorType& t1 = pMesh->elements[e].frame.t1;
  const VectorType derivativeT1_j =
      dofi < 3 && &v_j == &dui.v ? pMesh->elements[e].derivativeT1[0][dui.dofi]
                                 : VectorType(0, 0, 0);
  const VectorType derivativeT1_jplus1 =
      dofi < 3 && &v_jplus1 == &dui.v
          ? pMesh->elements[e].derivativeT1[1][dui.dofi]
          : VectorType(0, 0, 0);
  const VectorType derivativeT1 = derivativeT1_j + derivativeT1_jplus1;

  const VectorType derivativeOfRCrossT1 =
      computeDerivativeOfCrossProduct(r, derivativeR, t1, derivativeT1);
  const VectorType rCrossT1 = Cross(r, t1);
  const double normOfRCrossT1 = rCrossT1.Norm();
  const double derivativeNormRCrossT1 =
      computeDerivativeOfNorm(rCrossT1, derivativeOfRCrossT1);

  const VectorType t2Deriv_firstTerm =
      -(rCrossT1 * derivativeNormRCrossT1) / std::pow(normOfRCrossT1, 2);
  const VectorType t2Deriv_secondTerm = derivativeOfRCrossT1 / normOfRCrossT1;
  const VectorType t2Deriv = t2Deriv_firstTerm + t2Deriv_secondTerm;

  return t2Deriv;
}

VectorType DRMSimulationModel::computeDerivativeT3(
    const EdgeType& e,
    const DifferentiateWithRespectTo& dui) const {
  const Element& element = pMesh->elements[e];
  const VectorType& t1 = element.frame.t1;
  const VectorType& t2 = element.frame.t2;
  const VectorType t1CrossT2 = Cross(t1, t2);
  const VertexType& v_j = *e.cV(0);
  const size_t vi_j = pMesh->getIndex(v_j);
  const VertexType& v_jplus1 = *e.cV(1);
  const size_t vi_jplus1 = pMesh->getIndex(v_jplus1);
  const VectorType derivativeT1_j =
      dui.dofi < 3 && &v_j == &dui.v
          ? pMesh->elements[e].derivativeT1[0][dui.dofi]
          : VectorType(0, 0, 0);
  const VectorType derivativeT1_jplus1 =
      dui.dofi < 3 && &v_jplus1 == &dui.v
          ? pMesh->elements[e].derivativeT1[1][dui.dofi]
          : VectorType(0, 0, 0);
  const VectorType derivativeT1 = derivativeT1_j + derivativeT1_jplus1;

  VectorType derivativeT2(0, 0, 0);
  if (&v_j == &dui.v) {
    derivativeT2 = pMesh->elements[e].derivativeT2[0][dui.dofi];
  } else if (&v_jplus1 == &dui.v) {
    derivativeT2 = pMesh->elements[e].derivativeT2[1][dui.dofi];
  }

  const VectorType derivativeT1CrossT2 =
      computeDerivativeOfCrossProduct(t1, derivativeT1, t2, derivativeT2);
  const double derivativeOfNormT1CrossT2 =
      computeDerivativeOfNorm(t1CrossT2, derivativeT1CrossT2);
  const double normT1CrossT2 = t1CrossT2.Norm();

  const VectorType t3Deriv_firstTerm =
      -(t1CrossT2 * derivativeOfNormT1CrossT2) / std::pow(normT1CrossT2, 2);
  const VectorType t3Deriv_secondTerm = derivativeT1CrossT2 / normT1CrossT2;
  const VectorType t3Deriv = t3Deriv_firstTerm + t3Deriv_secondTerm;

  return t3Deriv;
}

double DRMSimulationModel::computeDerivativeTheta1(
    const EdgeType& e,
    const VertexIndex& evi,
    const VertexIndex& dwrt_evi,
    const DoFType& dwrt_dofi) const {
  const VertexType& v = *e.cV(evi);
  const size_t vi = pMesh->getIndex(v);
  const Element& element = pMesh->elements[e];
  const VectorType derivativeT1 = element.derivativeT1[dwrt_evi][dwrt_dofi];
  const VectorType derivativeT3 = element.derivativeT3[dwrt_evi][dwrt_dofi];
  const VectorType nDerivative =
      evi != dwrt_evi ? VectorType(0, 0, 0)
                      : pMesh->nodes[v].derivativeOfNormal[dwrt_dofi];
  const VectorType n = v.cN();
  const VectorType& t1 = element.frame.t1;
  const VectorType& t3 = element.frame.t3;
  const double theta1Derivative =
      derivativeT1 * Cross(t3, n) +
      t1 * (Cross(derivativeT3, n) + Cross(t3, nDerivative));

  return theta1Derivative;
}

double DRMSimulationModel::computeDerivativeTheta2(
    const EdgeType& e,
    const VertexIndex& evi,
    const VertexIndex& dwrt_evi,
    const DoFType& dwrt_dofi) const {
  const VertexType& v = *e.cV(evi);

  const Element& element = pMesh->elements[e];
  const VectorType derivativeT2 = element.derivativeT2[dwrt_evi][dwrt_dofi];
  const VectorType derivativeT3 = element.derivativeT3[dwrt_evi][dwrt_dofi];

  const VectorType n = v.cN();
  const VectorType& t2 = element.frame.t2;
  const VectorType& t3 = element.frame.t3;
  const VectorType nDerivative =
      dwrt_evi == evi ? pMesh->nodes[v].derivativeOfNormal[dwrt_dofi]
                      : VectorType(0, 0, 0);
  const double theta2Derivative =
      derivativeT2 * Cross(t3, n) +
      t2 * (Cross(derivativeT3, n) + Cross(t3, nDerivative));

  return theta2Derivative;
}

double DRMSimulationModel::computeTheta3(const EdgeType& e,
                                         const VertexType& v) {
  const VertexIndex& vi = pMesh->nodes[v].vi;

  const Element& elem = pMesh->elements[e];
  const EdgeIndex& ei = elem.ei;
  const Element::LocalFrame& ef = elem.frame;
  const VectorType& t1 = ef.t1;
  const VectorType& n = v.cN();
  const Node& node = pMesh->nodes[v];
  const double& nR = node.nR;

  // Use nR as theta3 only for the first star edge
  if (&e == node.referenceElement) {
    return nR;
  }
  const double alphaAngle =
      std::find_if(node.alphaAngles.begin(), node.alphaAngles.end(),
                   [&](const std::pair<EdgeIndex, double>& p) {
                     return elem.ei == p.first;
                   })
          ->second;
  const EdgeType& refElem = *node.referenceElement;
  const double rotationAngle = nR + alphaAngle;

  const VectorType& t1_k = pMesh->elements[refElem].frame.t1;
  const VectorType f1 = (t1_k - (n * (t1_k * n))).Normalize();
  vcg::Matrix44<double> rotationMatrix;
  rotationMatrix.SetRotateRad(rotationAngle, n);
  const double cosRotationAngle = cos(rotationAngle);
  const double sinRotationAngle = sin(rotationAngle);
  const VectorType f2 =
      (f1 * cosRotationAngle + Cross(n, f1) * sinRotationAngle).Normalize();
  const VectorType& t1Current = t1;
  const VectorType f3 = Cross(t1Current, f2);

  Element& element = pMesh->elements[e];
  // Save for computing theta3Derivative
  if (&v == e.cV(0)) {
    element.f1_j = f1;
    element.f2_j = f2;
    element.f3_j = f3;
    element.cosRotationAngle_j = cosRotationAngle;
    element.sinRotationAngle_j = sinRotationAngle;
  } else {
    element.f1_jplus1 = f1;
    element.f2_jplus1 = f2;
    element.f3_jplus1 = f3;
    element.cosRotationAngle_jplus1 = cosRotationAngle;
    element.sinRotationAngle_jplus1 = sinRotationAngle;
  }
  const double theta3 = f3.dot(n);

  return theta3;
}

double DRMSimulationModel::computeDerivativeTheta3(
    const EdgeType& e,
    const VertexType& v,
    const DifferentiateWithRespectTo& dui) const {
  const Node& node = pMesh->nodes[v];
  const VertexIndex& vi = pMesh->nodes[v].vi;
  if (&e == node.referenceElement && !isRigidSupport[vi]) {
    if (dui.dofi == DoF::Nr && &dui.v == &v) {
      return 1;
    } else {
      return 0;
    }
  }

  const Element& element = pMesh->elements[e];
  const Element::LocalFrame& ef = element.frame;
  const VectorType& t1 = ef.t1;
  const VectorType& n = v.cN();
  const DoFType& dofi = dui.dofi;
  const VertexPointer& vp_j = e.cV(0);
  const VertexPointer& vp_jplus1 = e.cV(1);

  double derivativeTheta3_dofi = 0;
  if (isRigidSupport[vi]) {
    const VectorType& t1Initial =
        computeT1Vector(pMesh->nodes[vp_j].initialLocation,
                        pMesh->nodes[vp_jplus1].initialLocation);
    VectorType g1 = Cross(t1, t1Initial);

    const VectorType derivativeInitialT1_dofi(0, 0, 0);
    const VectorType derivativeT1_j = dui.dofi < 3 && vp_j == &dui.v
                                          ? element.derivativeT1[0][dui.dofi]
                                          : VectorType(0, 0, 0);
    const VectorType derivativeT1_jplus1 =
        dui.dofi < 3 && vp_jplus1 == &dui.v ? element.derivativeT1[1][dui.dofi]
                                            : VectorType(0, 0, 0);
    const VectorType derivativeT1_dofi = derivativeT1_j + derivativeT1_jplus1;

    const VectorType derivativeG1_firstTerm =
        Cross(derivativeT1_dofi, t1Initial);
    const VectorType derivativeG1_secondTerm =
        Cross(t1, derivativeInitialT1_dofi);
    const VectorType derivativeG1 =
        derivativeG1_firstTerm + derivativeG1_secondTerm;
    const VectorType derivativeNormal = &v == &dui.v && dui.dofi > 2
                                            ? node.derivativeOfNormal[dui.dofi]
                                            : VectorType(0, 0, 0);
    derivativeTheta3_dofi = derivativeG1 * n + g1 * derivativeNormal;
    return derivativeTheta3_dofi;
  }
  EdgeType& refElem = *node.referenceElement;
  const VectorType& t1_k = pMesh->elements[refElem].frame.t1;
  VectorType f1, f2, f3;
  double cosRotationAngle, sinRotationAngle;
  if (e.cV(0) == &v) {
    f1 = element.f1_j;
    cosRotationAngle = element.cosRotationAngle_j;
    sinRotationAngle = element.sinRotationAngle_j;
    f2 = element.f2_j;
    f3 = element.f3_j;
  } else {
    f1 = element.f1_jplus1;
    cosRotationAngle = element.cosRotationAngle_jplus1;
    sinRotationAngle = element.sinRotationAngle_jplus1;
    f2 = element.f2_jplus1;
    f3 = element.f3_jplus1;
  }
  const VectorType& t1_kplus1 = t1;
  const VectorType f1Normalized = f1 / f1.Norm();

  VectorType derivativeF1(0, 0, 0);
  VectorType derivativeF2(0, 0, 0);
  VectorType derivativeF3(0, 0, 0);
  if (dui.dofi < 3) {
    const VectorType derivativeT1_kplus1_j =
        vp_j == &dui.v ? element.derivativeT1[0][dui.dofi]
                       : VectorType(0, 0, 0);
    const VectorType derivativeT1_kplus1_jplus1 =
        vp_jplus1 == &dui.v ? element.derivativeT1[1][dui.dofi]
                            : VectorType(0, 0, 0);
    const VectorType derivativeT1_kplus1 =
        derivativeT1_kplus1_j + derivativeT1_kplus1_jplus1;

    const VectorType derivativeT1_k_j =
        refElem.cV(0) == &dui.v
            ? pMesh->elements[refElem].derivativeT1[0][dui.dofi]
            : VectorType(0, 0, 0);
    const VectorType derivativeT1_k_jplus1 =
        refElem.cV(1) == &dui.v
            ? pMesh->elements[refElem].derivativeT1[1][dui.dofi]
            : VectorType(0, 0, 0);
    const VectorType derivativeT1_k = derivativeT1_k_j + derivativeT1_k_jplus1;

    derivativeF1 = derivativeT1_k - (n * (derivativeT1_k * n));

    const double f1Norm = f1.Norm();
    const double derivativeF1Norm = f1 * derivativeF1 / f1Norm;
    const VectorType derivativeF1Normalized =
        -f1 * derivativeF1Norm / (f1Norm * f1Norm) + derivativeF1 / f1Norm;

    derivativeF2 = derivativeF1Normalized * cosRotationAngle +
                   Cross(n, derivativeF1Normalized) * sinRotationAngle;
    const VectorType derivativeF3_firstTerm = Cross(derivativeT1_kplus1, f2);
    const VectorType derivativeF3_secondTerm = Cross(t1_kplus1, derivativeF2);
    derivativeF3 = derivativeF3_firstTerm + derivativeF3_secondTerm;
    derivativeTheta3_dofi = derivativeF3 * n;

  } else if (dui.dofi == DoF::Nr && &dui.v == &v) {
    derivativeF2 = f1Normalized * (-sinRotationAngle) +
                   Cross(n, f1Normalized) * cosRotationAngle;
    derivativeF3 = Cross(t1_kplus1, derivativeF2);
    derivativeTheta3_dofi = derivativeF3 * n;
  } else {  // 2<dofi<5
    if (&v == &dui.v) {
      const VectorType& derivativeOfNormal = node.derivativeOfNormal[dofi];
      derivativeF1 =
          -(n * (t1_k * derivativeOfNormal) + derivativeOfNormal * (t1_k * n));
      const double f1Norm = f1.Norm();
      const double derivativeF1Norm = f1 * derivativeF1 / f1Norm;
      const VectorType derivativeF1Normalized =
          -f1 * derivativeF1Norm / (f1Norm * f1Norm) + derivativeF1 / f1Norm;

      derivativeF2 = derivativeF1Normalized * cosRotationAngle +

                     (Cross(derivativeOfNormal, f1Normalized) +
                      Cross(n, derivativeF1Normalized)) *
                         sinRotationAngle;
      derivativeF3 = Cross(t1_kplus1, derivativeF2);
      derivativeTheta3_dofi = derivativeF3 * n + f3 * derivativeOfNormal;
    }
  }
  return derivativeTheta3_dofi;
}

void DRMSimulationModel::updateResidualForces() {
  std::vector<std::vector<std::pair<int, Vector6d>>>
      internalForcesContributionsFromEachEdge(
          pMesh->EN(),
          std::vector<std::pair<int, Vector6d>>(4, {-1, Vector6d()}));

  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->edge.begin(), pMesh->edge.end(), [&](const EdgeType& e) {
        const int ei = pMesh->getIndex(e);
        const SimulationEdgeMesh::VertexType& ev_j = *e.cV(0);
        const SimulationEdgeMesh::VertexType& ev_jplus1 = *e.cV(1);
        const Element& element = pMesh->elements[e];
        const Element::LocalFrame& ef = element.frame;
        const auto& n_a = ev_j.cN();
        const auto& n_b = ev_jplus1.cN();
        const double theta2_j = n_a * ef.t1;
        const double theta2_jplus1 = n_b * ef.t1;
        const double theta3_j = computeTheta3(e, ev_j);
        const double theta3_jplus1 = computeTheta3(e, ev_jplus1);
        std::vector<std::pair<int, Vector6d>>
            internalForcesContributionFromThisEdge(4, {-1, Vector6d()});
        for (VertexIndex evi = 0; evi < 2; evi++) {
          const SimulationEdgeMesh::VertexType& ev = *e.cV(evi);
          const Node& edgeNode = pMesh->nodes[ev];
          internalForcesContributionFromThisEdge[evi].first = edgeNode.vi;

          const VertexPointer& rev_j = edgeNode.referenceElement->cV(0);
          const VertexPointer& rev_jplus1 = edgeNode.referenceElement->cV(1);
          const VertexPointer& refElemOtherVertex =
              rev_j == &ev ? rev_jplus1 : rev_j;
          const Node& refElemOtherVertexNode = pMesh->nodes[refElemOtherVertex];
          if (edgeNode.referenceElement != &e) {
            internalForcesContributionFromThisEdge[evi + 2].first =
                refElemOtherVertexNode.vi;
          }
          const size_t vi = edgeNode.vi;
          for (DoFType dofi = DoF::Ux; dofi < DoF::NumDoF; dofi++) {
            const bool isDofConstrainedFor_ev =
                (isVertexConstrained[edgeNode.vi] &&
                 constrainedVertices.at(vi).contains(dofi));
            if (!isDofConstrainedFor_ev) {
              DifferentiateWithRespectTo dui{ev, dofi};
              // Axial force computation
              const double axialForce_dofi = [&]() {
                if (dofi > Uz) {
                  return 0.0;
                }
                double e_kDeriv_opt = (evi == 0 ? -element.frame.t1[dofi]
                                                : element.frame.t1[dofi]);
                const double e_k = element.length - element.initialLength;
                const double axialForce_dofi =
                    e_kDeriv_opt * e_k * element.rigidity.axial;
                return axialForce_dofi;
              }();

              const auto& d_n_a =
                  &dui.v != &ev_j ? VectorType(0, 0, 0)
                                  : pMesh->nodes[ev_j].derivativeOfNormal[dofi];
              const auto& d_n_b =
                  &dui.v != &ev_jplus1
                      ? VectorType(0, 0, 0)
                      : pMesh->nodes[ev_jplus1].derivativeOfNormal[dofi];
              const auto& t2Deriv = element.derivativeT2[evi][dofi];
              const double theta1DiffDerivative_opt =
                  (t2Deriv * (n_a - n_b)) +
                  (element.frame.t2 * (d_n_a - d_n_b));
              const double torsionalForce_dofi =
                  element.rigidity.torsional *
                  (element.frame.t2 * (n_a - n_b)) * theta1DiffDerivative_opt;

              // First bending force computation
              ////theta2_j derivative
              const double theta2_j_deriv =
                  d_n_a * element.frame.t1 +
                  n_a * element.derivativeT1[evi][dofi];
              ////theta2_jplus1 derivative
              const double theta2_jplus1_deriv =
                  d_n_b * element.frame.t1 +
                  n_b * element.derivativeT1[evi][dofi];
              ////1st in bracket term
              const double firstBendingForce_inBracketsTerm_0 =
                  theta2_j_deriv * 2 * theta2_j;
              ////2nd in bracket term
              const double firstBendingForce_inBracketsTerm_1 =
                  theta2_jplus1_deriv * theta2_j;
              ////3rd in bracket term
              const double firstBendingForce_inBracketsTerm_2 =
                  theta2_j_deriv * theta2_jplus1;
              ////4th in bracket term
              const double firstBendingForce_inBracketsTerm_3 =
                  2 * theta2_jplus1_deriv * theta2_jplus1;
              // 3rd term computation
              const double firstBendingForce_inBracketsTerm =
                  firstBendingForce_inBracketsTerm_0 +
                  firstBendingForce_inBracketsTerm_1 +
                  firstBendingForce_inBracketsTerm_2 +
                  firstBendingForce_inBracketsTerm_3;
              const double firstBendingForce_dofi =
                  firstBendingForce_inBracketsTerm *
                  element.rigidity.firstBending;

              // Second bending force computation
              ////theta2_j derivative
              const double theta3_j_deriv =
                  computeDerivativeTheta3(e, ev_j, dui);
              ////theta2_jplus1 derivative
              const double theta3_jplus1_deriv =
                  computeDerivativeTheta3(e, ev_jplus1, dui);
              ////1st in bracket term
              const double secondBendingForce_inBracketsTerm_0 =
                  theta3_j_deriv * 2 * theta3_j;
              ////2nd in bracket term
              const double secondBendingForce_inBracketsTerm_1 =
                  theta3_jplus1_deriv * theta3_j;
              ////3rd in bracket term
              const double secondBendingForce_inBracketsTerm_2 =
                  theta3_j_deriv * theta3_jplus1;
              ////4th in bracket term
              const double secondBendingForce_inBracketsTerm_3 =
                  theta3_jplus1_deriv * 2 * theta3_jplus1;
              // 3rd term computation
              const double secondBendingForce_inBracketsTerm =
                  secondBendingForce_inBracketsTerm_0 +
                  secondBendingForce_inBracketsTerm_1 +
                  secondBendingForce_inBracketsTerm_2 +
                  secondBendingForce_inBracketsTerm_3;
              double secondBendingForce_dofi =
                  secondBendingForce_inBracketsTerm *
                  element.rigidity.secondBending;
              internalForcesContributionFromThisEdge[evi].second[dofi] =
                  axialForce_dofi + firstBendingForce_dofi +
                  secondBendingForce_dofi + torsionalForce_dofi;
            }
            if (edgeNode.referenceElement != &e) {
              const bool isDofConstrainedFor_refElemOtherVertex =
                  (isVertexConstrained[refElemOtherVertexNode.vi] &&
                   constrainedVertices.at(refElemOtherVertexNode.vi)
                       .contains(dofi));
              if (!isDofConstrainedFor_refElemOtherVertex) {
                DifferentiateWithRespectTo dui{*refElemOtherVertex, dofi};
                ////theta3_j derivative
                const double theta3_j_deriv =
                    computeDerivativeTheta3(e, ev_j, dui);
                ////theta3_jplus1 derivative
                const double theta3_jplus1_deriv =
                    computeDerivativeTheta3(e, ev_jplus1, dui);
                ////1st in bracket term
                const double secondBendingForce_inBracketsTerm_0 =
                    theta3_j_deriv * 2 * theta3_j;
                ////2nd in bracket term
                const double secondBendingForce_inBracketsTerm_1 =
                    theta3_jplus1_deriv * theta3_j;
                ////3rd in bracket term
                const double secondBendingForce_inBracketsTerm_2 =
                    theta3_j_deriv * theta3_jplus1;
                ////4th in bracket term
                const double secondBendingForce_inBracketsTerm_3 =
                    theta3_jplus1_deriv * 2 * theta3_jplus1;

                // 4th term computation
                const double secondBendingForce_inBracketsTerm =
                    secondBendingForce_inBracketsTerm_0 +
                    secondBendingForce_inBracketsTerm_1 +
                    secondBendingForce_inBracketsTerm_2 +
                    secondBendingForce_inBracketsTerm_3;
                const double secondBendingForce_dofi =
                    secondBendingForce_inBracketsTerm *
                    element.rigidity.secondBending;
                internalForcesContributionFromThisEdge[evi + 2].second[dofi] =
                    secondBendingForce_dofi;
              }
            }
          }
        }
        internalForcesContributionsFromEachEdge[ei] =
            internalForcesContributionFromThisEdge;
      });

  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->nodes._handle->data.begin(), pMesh->nodes._handle->data.end(),
      [](Node& node) {
        Node::Forces& force = node.force;
        force.residual = force.external;
        force.internal = 0;
      });

  for (size_t ei = 0; ei < pMesh->EN(); ei++) {
    for (int i = 0; i < 4; i++) {
      std::pair<int, Vector6d> internalForcePair =
          internalForcesContributionsFromEachEdge[ei][i];
      int vi = internalForcePair.first;
      if (i > 1 && vi == -1) {
        continue;
      }
      Node::Forces& force = pMesh->nodes[vi].force;
      force.internal = force.internal + internalForcePair.second;
      force.residual = force.residual + (internalForcePair.second * -1);
    }
  }

  if (mSettings.initialDistortion.has_value() &&
      mCurrentSimulationStep >=
          mSettings.initialDistortion->duration.startStep &&
      mCurrentSimulationStep < mSettings.initialDistortion->duration.endStep) {
    std::for_each(
#ifdef ENABLE_PARALLEL_DRM
        std::execution::par_unseq,
#endif
        pMesh->nodes._handle->data.begin(), pMesh->nodes._handle->data.end(),
        [&](Node& node) {
          Node::Forces& force = node.force;
          if (!constrainedVertices.contains(node.vi) ||
              !constrainedVertices.at(node.vi).contains(static_cast<int>(Uz))) {
            force.residual = force.residual + mSettings.initialDistortion->load;
          }
        });
  }

  double residualForceNorm_structure = 0;
  Vector6d sumOfResidualForces(0);
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->vert.begin(), pMesh->vert.end(), [&](const VertexType& v) {
        Node& node = pMesh->nodes[v];
        Node::Forces& force = node.force;
        Vector6d& nodeResidualForce = force.residual;
        sumOfResidualForces = sumOfResidualForces + nodeResidualForce;
        const double residualForceNorm_vertex = nodeResidualForce.norm();
        residualForceNorm_structure += residualForceNorm_vertex;
      });
  pMesh->totalResidualForcesNorm = residualForceNorm_structure;
  pMesh->averageResidualForcesNorm = residualForceNorm_structure / pMesh->VN();
}

void DRMSimulationModel::updateNodalExternalForces(
    const std::unordered_map<VertexIndex, Vector6d>& nodalForces,
    const std::unordered_map<VertexIndex, std::unordered_set<DoFType>>&
        fixedVertices) {
  externalMomentsNorm = 0;
  double totalExternalForcesNorm = 0;
  for (const std::pair<VertexIndex, Vector6d>& nodalForce : nodalForces) {
    const VertexIndex nodeIndex = nodalForce.first;
    const bool isNodeConstrained = fixedVertices.contains(nodeIndex);
    Node& node = pMesh->nodes[nodeIndex];
    Vector6d nodalExternalForce(0);
    for (DoFType dofi = DoF::Ux; dofi < DoF::NumDoF; dofi++) {
      const bool isDofConstrained =
          isNodeConstrained && fixedVertices.at(nodeIndex).contains(dofi);
      if (isDofConstrained) {
        continue;
      }
      nodalExternalForce[dofi] = nodalForce.second[dofi];
    }
    externalMomentsNorm += std::sqrt(pow(nodalExternalForce[3], 2) +
                                     pow(nodalExternalForce[4], 2) +
                                     pow(nodalExternalForce[5], 2));

    /*
     * The external moments are given as a rotation around an axis.
     * In this implementation we model moments as rotation of the normal vector
     * and because of that we need to transform the moments.
     */

    if (externalMomentsNorm != 0) {
      VectorType momentAxis(
          nodalExternalForce[3], nodalExternalForce[4],
          nodalExternalForce[5]);  // rotation around this vector
      VectorType transformedVector =
          vcg::RotationMatrix(VectorType(0, 0, 1), vcg::math::ToRad(-90.0)) *
          momentAxis;
      nodalExternalForce[3] = transformedVector[0];
      nodalExternalForce[4] = transformedVector[1];
      nodalExternalForce[5] = transformedVector[2];
    }

    node.force.external = nodalExternalForce;
    totalExternalForcesNorm += node.force.external.norm();
  }

  pMesh->totalExternalForcesNorm = totalExternalForcesNorm;
}

void DRMSimulationModel::computeRigidSupports() {
  isRigidSupport.clear();
  isRigidSupport.resize(pMesh->VN(), false);
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->vert.begin(), pMesh->vert.end(), [&](const VertexType& v) {
        const VertexIndex vi = pMesh->nodes[v].vi;
        const bool isVertexConstrained = constrainedVertices.contains(vi);
        if (isVertexConstrained) {
          auto constrainedDoFType = constrainedVertices.at(vi);
          const bool hasAllDoFTypeConstrained =
              constrainedDoFType.contains(DoF::Ux) &&
              constrainedDoFType.contains(DoF::Uy) &&
              constrainedDoFType.contains(DoF::Uz) &&
              constrainedDoFType.contains(DoF::Nx) &&
              constrainedDoFType.contains(DoF::Ny) &&
              constrainedDoFType.contains(DoF::Nr);
          if (hasAllDoFTypeConstrained) {
            isRigidSupport[vi] = true;
          }
        }
      });
}

void DRMSimulationModel::updateNormalDerivatives() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->vert.begin(), pMesh->vert.end(), [&](const VertexType& v) {
        for (DoFType dofi = DoF::Nx; dofi < DoF::NumDoF; dofi++) {
          DifferentiateWithRespectTo dui{v, dofi};
          pMesh->nodes[v].derivativeOfNormal[dofi] =
              computeDerivativeOfNormal(v, dui);
        }
      });
}

void DRMSimulationModel::updateT1Derivatives() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->edge.begin(), pMesh->edge.end(), [&](EdgeType& e) {
        for (DoFType dofi = DoF::Ux; dofi < DoF::Nx; dofi++) {
          Element& element = pMesh->elements[e];
          DifferentiateWithRespectTo dui_v0{*e.cV(0), dofi};
          element.derivativeT1[0][dofi] = computeDerivativeT1(e, dui_v0);
          DifferentiateWithRespectTo dui_v1{*e.cV(1), dofi};
          element.derivativeT1[1][dofi] = computeDerivativeT1(e, dui_v1);
        }
      });
}

void DRMSimulationModel::updateT2Derivatives() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->edge.begin(), pMesh->edge.end(), [&](EdgeType& e) {
        for (DoFType dofi = DoF::Ux; dofi < DoF::NumDoF; dofi++) {
          DifferentiateWithRespectTo dui_v0{*e.cV(0), dofi};
          DifferentiateWithRespectTo dui_v1{*e.cV(1), dofi};
          Element& element = pMesh->elements[e];
          element.derivativeT2[0][dofi] = computeDerivativeT2(e, dui_v0);
          element.derivativeT2[1][dofi] = computeDerivativeT2(e, dui_v1);
        }
      });
}

void DRMSimulationModel::updateT3Derivatives() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->edge.begin(), pMesh->edge.end(), [&](EdgeType& e) {
        for (DoFType dofi = DoF::Ux; dofi < DoF::NumDoF; dofi++) {
          Element& element = pMesh->elements[e];
          DifferentiateWithRespectTo dui_v0{*e.cV(0), dofi};
          element.derivativeT3[0][dofi] = computeDerivativeT3(e, dui_v0);
          DifferentiateWithRespectTo dui_v1{*e.cV(1), dofi};
          element.derivativeT3[1][dofi] = computeDerivativeT3(e, dui_v1);
        }
      });
}

void DRMSimulationModel::updateRDerivatives() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->edge.begin(), pMesh->edge.end(), [&](EdgeType& e) {
        for (DoFType dofi = DoF::Nx; dofi < DoF::NumDoF; dofi++) {
          DifferentiateWithRespectTo dui_v0{*e.cV(0), dofi};
          pMesh->elements[e].derivativeR[0][dofi] =
              computeDerivativeOfR(e, dui_v0);
          DifferentiateWithRespectTo dui_v1{*e.cV(1), dofi};
          pMesh->elements[e].derivativeR[1][dofi] =
              computeDerivativeOfR(e, dui_v1);
        }
      });
}

void DRMSimulationModel::updateElementalLengths() {
  pMesh->updateElementalLengths();
}

DRMSimulationModel::DRMSimulationModel() {}

void DRMSimulationModel::updateNodalMasses() {
  double gamma = mSettings.gamma;
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->vert.begin(), pMesh->vert.end(), [&](const VertexType& v) {
        double translationalSumSk = 0;
        double rotationalSumSk = 0;
        for (const EdgePointer& ep : pMesh->nodes[v].incidentElements) {
          const Element& elem = pMesh->elements[ep];
          const double SkTranslational =
              elem.material.youngsModulus * elem.dimensions.A / elem.length;
          translationalSumSk += SkTranslational;
          const double lengthToThe3 = std::pow(elem.length, 3);
          const long double SkRotational =
              elem.material.youngsModulus *
              (elem.dimensions.inertia.I2 + elem.dimensions.inertia.I3) /
              lengthToThe3;
          rotationalSumSk += SkRotational;
        }
        pMesh->nodes[v].mass.translational =
            gamma * pow(mSettings.Dtini, 2) * 2 * translationalSumSk;
        pMesh->nodes[v].mass.rotational =
            gamma * pow(mSettings.Dtini, 2) * 8 * rotationalSumSk;

        // Check that Courant–Friedrichs–Lewy's condition holds for the masses
        if (std::pow(mSettings.Dtini, 2.0) * translationalSumSk /
                    pMesh->nodes[v].mass.translational >=
                2 ||
            std::pow(mSettings.Dtini, 2.0) * rotationalSumSk /
                    pMesh->nodes[v].mass.rotational >=
                2) {
          throw std::runtime_error(
              "Courant–Friedrichs–Lewy's condition does not hold");
        }
      });
}

void DRMSimulationModel::updateNodalAccelerations() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->vert.begin(), pMesh->vert.end(), [&](const VertexType& v) {
        Node& node = pMesh->nodes[v];
        node.acceleration.setTranslation(node.force.residual.getTranslation() /
                                         node.mass.translational);
        node.acceleration.setRotation(node.force.residual.getRotation() /
                                      node.mass.rotational);
      });
}

void DRMSimulationModel::updateNodalVelocities() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->vert.begin(), pMesh->vert.end(), [&](const VertexType& v) {
        Node& node = pMesh->nodes[v];
        node.velocity = node.velocity + node.acceleration * Dt;
      });
  updateKineticEnergy();
}

void DRMSimulationModel::updateNodalDisplacements() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->vert.begin(), pMesh->vert.end(), [&](const VertexType& v) {
        Node& node = pMesh->nodes[v];
        Vector6d stepDisplacement = node.velocity * Dt;
        node.displacements = node.displacements + stepDisplacement;
      });
}

void DRMSimulationModel::updateNodePosition(
    VertexType& v,
    const std::unordered_map<VertexIndex, std::unordered_set<DoFType>>&
        fixedVertices) {
  Node& node = pMesh->nodes[v];
  const VertexIndex& vi = pMesh->nodes[v].vi;

  VectorType displacementVector(0, 0, 0);
  displacementVector += VectorType(node.displacements[0], 0, 0);
  displacementVector += VectorType(0, node.displacements[1], 0);
  displacementVector += VectorType(0, 0, node.displacements[2]);

  v.P() = node.initialLocation + displacementVector;
}

void DRMSimulationModel::updateNodeNr(VertexType& v) {
  const VertexIndex& vi = pMesh->nodes[v].vi;
  Node& node = pMesh->nodes[v];
  if (!isRigidSupport[vi]) {
    node.nR = node.displacements[5];
  } else {
    const EdgePointer& refElem = node.referenceElement;
    const VectorType& refT1 = pMesh->elements[refElem].frame.t1;

    const VectorType& t1Initial =
        computeT1Vector(pMesh->nodes[refElem->cV(0)].initialLocation,
                        pMesh->nodes[refElem->cV(1)].initialLocation);
    VectorType g1 = Cross(refT1, t1Initial);
    node.nR = g1.dot(v.cN());
    node.displacements[Nr] = node.nR;
  }
}

void DRMSimulationModel::updateNodeNormal(
    VertexType& v,
    const std::unordered_map<VertexIndex, std::unordered_set<DoFType>>&
        fixedVertices) {
  Node& node = pMesh->nodes[v];
  const VertexIndex& vi = node.vi;
  VectorType normalDisplacementVector(0, 0, 0);
  normalDisplacementVector += VectorType(node.displacements[3], 0, 0);
  normalDisplacementVector += VectorType(0, node.displacements[4], 0);
  v.N() = node.initialNormal + normalDisplacementVector;
  const double &nx = v.N()[0], ny = v.N()[1];
  const double nxnyMagnitude = std::pow(nx, 2) + std::pow(ny, 2);
  if (nxnyMagnitude > 1) {
    VectorType newNormal(nx / std::sqrt(nxnyMagnitude),
                         ny / std::sqrt(nxnyMagnitude), 0);
    v.N() = newNormal;

    /*If an external moment caused the normal to lay on the xy plane this
     * means that in order to disable its effect a greater internal force is
     * needed than what is possible (the constraint on the z of the normals
     * imposes a constraint on the maximum internal force). Because of that
     * the totalResidualForcesNorm can't drop below the magnitude of external
     * moment applied on vertex vi. In order to allow termination of the
     * simulation when the described phenomenon happens we allow the
     * termination of the algorithm if the kinetic energy of the system drops
     * below the set threshold.
     * */
    const bool viHasMoments =
        node.force.external[3] != 0 || node.force.external[4] != 0;
    if (!checkedForMaximumMoment && viHasMoments) {
      checkedForMaximumMoment = true;
    }

  } else {
    const double nzSquared = 1.0 - nxnyMagnitude;
    const double nz = std::sqrt(nzSquared);
    VectorType newNormal(nx, ny, nz);
    v.N() = newNormal;
  }
  node.displacements[DoF::Nx] = v.cN()[0];
  node.displacements[DoF::Ny] = v.cN()[1];
}

void DRMSimulationModel::applyDisplacements(
    const std::unordered_map<VertexIndex, std::unordered_set<DoFType>>&
        fixedVertices) {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->vert.begin(), pMesh->vert.end(), [&](auto& v) {
        updateNodePosition(v, fixedVertices);
        updateNodeNormal(v, fixedVertices);
        updateNodeNr(v);
      });

  updateElementalFrames();
  if (mSettings.shouldDraw) {
    pMesh->updateEigenEdgeAndVertices();
  }
}

void DRMSimulationModel::updateElementalFrames() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->edge.begin(), pMesh->edge.end(), [&](const EdgeType& e) {
        const VectorType elementNormal =
            (e.cV(0)->cN() + e.cV(1)->cN()).Normalize();
        pMesh->elements[e].frame =
            computeElementFrame(e.cP(0), e.cP(1), elementNormal);
      });
}

void DRMSimulationModel::applyForcedDisplacements(
    const std::unordered_map<VertexIndex, Eigen::Vector3d>&
        nodalForcedDisplacements) {
  if (nodalForcedDisplacements.empty() ||
      (mSettings.gradualForcedDisplacementSteps.has_value() &&
       mSettings.gradualForcedDisplacementSteps->startStep >
           mCurrentSimulationStep)) {
    return;
  }
  for (const std::pair<VertexIndex, Eigen::Vector3d>
           vertexIndexDisplacementPair : nodalForcedDisplacements) {
    const VertexIndex vi = vertexIndexDisplacementPair.first;
    const Eigen::Vector3d vertexDisplacement =
        vertexIndexDisplacementPair.second;
    Node& node = pMesh->nodes[vi];
    VectorType displacementVector(vertexDisplacement(0), vertexDisplacement(1),
                                  vertexDisplacement(2));

    if (mSettings.gradualForcedDisplacementSteps.has_value() &&
        mSettings.gradualForcedDisplacementSteps->isActive(
            mCurrentSimulationStep)) {
      const double forcedDisplacementScalingFactor =
          (mCurrentSimulationStep -
           mSettings.gradualForcedDisplacementSteps->startStep) /
          static_cast<double>(
              mSettings.gradualForcedDisplacementSteps->getSize());
      displacementVector *= forcedDisplacementScalingFactor;
    }
    node.displacements = Vector6d(
        {displacementVector[0], displacementVector[1], displacementVector[2],
         node.displacements[3], node.displacements[4], node.displacements[5]});
  }

  if (mSettings.shouldDraw) {
    pMesh->updateEigenEdgeAndVertices();
  }
}

void DRMSimulationModel::updateKineticEnergy() {
  pMesh->previousTotalKineticEnergy = pMesh->currentTotalKineticEnergy;
  pMesh->previousTotalTranslationalKineticEnergy =
      pMesh->currentTotalTranslationalKineticEnergy;
  pMesh->previousTotalRotationalKineticEnergy =
      pMesh->currentTotalRotationalKineticEnergy;
  pMesh->currentTotalKineticEnergy = 0;
  pMesh->currentTotalTranslationalKineticEnergy = 0;
  pMesh->currentTotalRotationalKineticEnergy = 0;

  for (const VertexType& v : pMesh->vert) {
    Node& node = pMesh->nodes[v];
    node.kineticEnergy = 0;

    const double nodeTranslationalKineticEnergy =
        0.5 * node.mass.translational *
        std::pow(node.velocity.getTranslation().norm(), 2);

    const double nodeRotationalKineticEnergy =
        0.5 * node.mass.rotational *
        std::pow(node.velocity.getRotation().norm(), 2);

    node.kineticEnergy =
        nodeTranslationalKineticEnergy + nodeRotationalKineticEnergy;
    if (std::isinf(nodeTranslationalKineticEnergy) /*> 1e15*/) {
      throw std::runtime_error(
          "Infinite kinetic energy detected. Aborting simulation..");
    }

    pMesh->currentTotalKineticEnergy += node.kineticEnergy;
    pMesh->currentTotalTranslationalKineticEnergy +=
        nodeTranslationalKineticEnergy;
    pMesh->currentTotalRotationalKineticEnergy += nodeRotationalKineticEnergy;
  }
}

void DRMSimulationModel::resetVelocities() {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->nodes._handle->data.begin(), pMesh->nodes._handle->data.end(),
      [&](Node& node) { node.velocity = node.acceleration * Dt * 0.5; });
  updateKineticEnergy();

  updateNodalDisplacements();
  applyDisplacements(constrainedVertices);
  updateElementalLengths();
}

SimulationResults DRMSimulationModel::computeResults(
    const std::shared_ptr<SimulationJob>& pJob) {
  std::vector<Vector6d> displacements(pMesh->VN());
  for (size_t vi = 0; vi < pMesh->VN(); vi++) {
    displacements[vi] = pMesh->nodes[vi].displacements;
  }
  history.numberOfSteps = mCurrentSimulationStep;
  SimulationResults results;
  results.converged = true;
  results.pJob = pJob;
  results.history = history;
  results.displacements = displacements;

  if (mSettings.maxDRMIterations.has_value() &&
      mCurrentSimulationStep == mSettings.maxDRMIterations &&
      mCurrentSimulationStep != 0) {
    if (mSettings.beVerbose) {
      std::cout
          << "Did not reach equilibrium before reaching the maximum number "
             "of DRM steps ("
          << mSettings.maxDRMIterations.value() << "). Breaking simulation"
          << std::endl;
    }
    results.converged = false;
  }
#ifdef POLYSCOPE_DEFINED
  if (mSettings.shouldDraw && !mSettings.debugModeStep.has_value()) {
    draw(pJob);
  }
#endif
  if (!std::isnan(pMesh->currentTotalKineticEnergy)) {
    results.debug_drmDisplacements = results.displacements;
    results.rotationalDisplacementQuaternion.resize(pMesh->VN());
    results.debug_q_f1.resize(pMesh->VN());
    results.debug_q_normal.resize(pMesh->VN());
    results.debug_q_nr.resize(pMesh->VN());
    for (int vi = 0; vi < pMesh->VN(); vi++) {
      const Node& node = pMesh->nodes[vi];
      const Eigen::Vector3d nInitial_eigen =
          node.initialNormal.ToEigenVector<Eigen::Vector3d>();
      const Eigen::Vector3d nDeformed_eigen =
          pMesh->vert[vi].cN().ToEigenVector<Eigen::Vector3d>();

      Eigen::Quaternion<double> q_normal;
      q_normal.setFromTwoVectors(nInitial_eigen, nDeformed_eigen);
      Eigen::Quaternion<double> q_nr_nDeformed;
      q_nr_nDeformed =
          Eigen::AngleAxis<double>(pMesh->nodes[vi].nR, nDeformed_eigen);
      Eigen::Quaternion<double> q_nr_nInit;
      q_nr_nInit =
          Eigen::AngleAxis<double>(pMesh->nodes[vi].nR, nInitial_eigen);
      const auto theta = 2 * acos(q_nr_nDeformed.w());
      Eigen::Vector3d deformedNormal_debug(q_nr_nDeformed.x() * sin(theta / 2),
                                           q_nr_nDeformed.y() * sin(theta / 2),
                                           q_nr_nDeformed.z() * sin(theta / 2));
      deformedNormal_debug.normalize();
      const double nr_debug =
          deformedNormal_debug.dot(nDeformed_eigen) > 0 ? theta : -theta;
      assert(pMesh->nodes[vi].nR - nr_debug < 1e-6);
      VectorType referenceT1_deformed =
          pMesh->elements[node.referenceElement].frame.t1;
      const VectorType& nDeformed = pMesh->vert[vi].cN();
      const VectorType referenceF1_deformed =
          (referenceT1_deformed -
           (node.initialNormal * (referenceT1_deformed * node.initialNormal)))
              .Normalize();

      const VectorType referenceT1_initial = computeT1Vector(
          pMesh->nodes[node.referenceElement->cV(0)].initialLocation,
          pMesh->nodes[node.referenceElement->cV(1)].initialLocation);
      const VectorType referenceF1_initial =
          (referenceT1_initial -
           (node.initialNormal * (referenceT1_initial * node.initialNormal)))
              .Normalize();
      Eigen::Quaternion<double> q_f1_nInit;  // nr is with respect to f1
      q_f1_nInit.setFromTwoVectors(
          referenceF1_initial.ToEigenVector<Eigen::Vector3d>(),
          referenceF1_deformed.ToEigenVector<Eigen::Vector3d>());

      Eigen::Quaternion<double> q_f1_nDeformed;  // nr is with respect to f1
      const VectorType referenceF1_initial_def =
          (referenceT1_initial -
           (nDeformed * (referenceT1_initial * nDeformed)))
              .Normalize();
      const VectorType referenceF1_deformed_def =
          (referenceT1_deformed -
           (nDeformed * (referenceT1_deformed * nDeformed)))
              .Normalize();
      q_f1_nDeformed.setFromTwoVectors(
          referenceF1_initial_def.ToEigenVector<Eigen::Vector3d>(),
          referenceF1_deformed_def.ToEigenVector<Eigen::Vector3d>());
      results.debug_q_f1[vi] = q_f1_nInit;
      results.debug_q_normal[vi] = q_normal;
      results.debug_q_nr[vi] = q_nr_nInit;
      results.rotationalDisplacementQuaternion[vi] =
          (q_normal * (q_f1_nInit * q_nr_nInit));
      // Update the displacement vector to contain the euler angles
      const Eigen::Vector3d eulerAngles =
          results
              .rotationalDisplacementQuaternion[vi]
              //            R
              .toRotationMatrix()
              .eulerAngles(0, 1, 2);
      results.displacements[vi][3] = eulerAngles[0];
      results.displacements[vi][4] = eulerAngles[1];
      results.displacements[vi][5] = eulerAngles[2];
    }
  }

  results.simulationModelUsed = label;

  return results;
}

void DRMSimulationModel::printCurrentState() const {
  std::cout << "Simulation steps executed:" << mCurrentSimulationStep
            << std::endl;
  std::cout << "Residual forces norm: " << pMesh->totalResidualForcesNorm
            << std::endl;
  if (pMesh->totalExternalForcesNorm > 1e-15) {
    std::cout << "Average Residual forces norm/extForcesNorm: "
              << pMesh->totalResidualForcesNorm / pMesh->VN() /
                     pMesh->totalExternalForcesNorm
              << std::endl;
  }

  std::cout << "Kinetic energy:" << pMesh->currentTotalKineticEnergy
            << std::endl;

  const double timePerNodePerIteration =
      std::chrono::duration_cast<std::chrono::microseconds>(
          std::chrono::steady_clock::now() - beginTime)
          .count() *
      1e-6 / (static_cast<double>(mCurrentSimulationStep) * pMesh->VN());
  std::cout << "time(s)/(iterations*node) = " << timePerNodePerIteration
            << std::endl;
  std::cout << "Dt:" << Dt << std::endl;
}

#ifdef POLYSCOPE_DEFINED
void DRMSimulationModel::draw(const std::shared_ptr<SimulationJob>& pJob,
                              const std::string& screenshotsFolder) {
  // update positions
  auto deformedMeshPolyscopeLabel = pMesh->getLabel();
  polyscope::CurveNetwork* meshPolyscopeHandle =
      polyscope::getCurveNetwork(deformedMeshPolyscopeLabel);
  meshPolyscopeHandle->updateNodePositions(pMesh->getEigenVertices());

  // Vertex quantities
  std::vector<double> kineticEnergies(pMesh->VN());
  std::vector<std::array<double, 3>> nodeNormals(pMesh->VN());
  std::vector<std::array<double, 3>> internalForces(pMesh->VN());
  std::vector<std::array<double, 3>> externalForces(pMesh->VN());
  std::vector<std::array<double, 3>> externalMoments(pMesh->VN());
  std::vector<double> internalForcesNorm(pMesh->VN());
  std::vector<double> nRs(pMesh->VN());
  std::vector<std::array<double, 3>> residualForces(pMesh->VN());
  std::vector<double> residualForcesNorm(pMesh->VN());
  std::vector<double> accelerationX(pMesh->VN());
  for (const VertexType& v : pMesh->vert) {
    kineticEnergies[pMesh->getIndex(v)] = pMesh->nodes[v].kineticEnergy;
    const VectorType n = v.cN();
    nodeNormals[pMesh->getIndex(v)] = {n[0], n[1], n[2]};
    // per node internal forces
    const Vector6d nodeforce = pMesh->nodes[v].force.internal * (-1);
    internalForces[pMesh->getIndex(v)] = {nodeforce[0], nodeforce[1],
                                          nodeforce[2]};
    internalForcesNorm[pMesh->getIndex(v)] = nodeforce.norm();
    // External force
    const Vector6d nodeExternalForce = pMesh->nodes[v].force.external;
    externalForces[pMesh->getIndex(v)] = {
        nodeExternalForce[0], nodeExternalForce[1], nodeExternalForce[2]};
    externalMoments[pMesh->getIndex(v)] = {nodeExternalForce[3],
                                           nodeExternalForce[4], 0};
    const Node& node = pMesh->nodes[v];
    nRs[pMesh->getIndex(v)] = node.nR;
    const Vector6d nodeResidualForce = pMesh->nodes[v].force.residual;
    residualForces[pMesh->getIndex(v)] = {
        nodeResidualForce[0], nodeResidualForce[1], nodeResidualForce[2]};
    residualForcesNorm[pMesh->getIndex(v)] = nodeResidualForce.norm();
    accelerationX[pMesh->getIndex(v)] = pMesh->nodes[v].acceleration[0];
  }
  meshPolyscopeHandle->addNodeScalarQuantity("Kinetic Energy", kineticEnergies)
      ->setEnabled(false);
  meshPolyscopeHandle->addNodeVectorQuantity("Node normals", nodeNormals)
      ->setEnabled(false);
  meshPolyscopeHandle->addNodeVectorQuantity("Internal force", internalForces)
      ->setEnabled(false);
  meshPolyscopeHandle
      ->addNodeScalarQuantity("Internal force norm", internalForcesNorm)
      ->setEnabled(true);
  pJob->registerForDrawing(pMesh->getLabel(), true);
  polyscope::getCurveNetwork(deformedMeshPolyscopeLabel)
      ->addNodeScalarQuantity("nR", nRs)
      ->setEnabled(false);
  polyscope::getCurveNetwork(deformedMeshPolyscopeLabel)
      ->addNodeVectorQuantity("Residual force", residualForces)
      ->setEnabled(true);
  polyscope::getCurveNetwork(deformedMeshPolyscopeLabel)
      ->addNodeScalarQuantity("Residual force norm", residualForcesNorm)
      ->setEnabled(false);

  // Edge quantities
  std::vector<double> A(pMesh->EN());
  std::vector<double> J(pMesh->EN());
  std::vector<double> I2(pMesh->EN());
  std::vector<double> I3(pMesh->EN());
  for (const EdgeType& e : pMesh->edge) {
    const size_t ei = pMesh->getIndex(e);
    A[ei] = pMesh->elements[e].dimensions.A;
    J[ei] = pMesh->elements[e].dimensions.inertia.J;
    I2[ei] = pMesh->elements[e].dimensions.inertia.I2;
    I3[ei] = pMesh->elements[e].dimensions.inertia.I3;
  }

  polyscope::getCurveNetwork(deformedMeshPolyscopeLabel)
      ->addEdgeScalarQuantity("A", A);
  polyscope::getCurveNetwork(deformedMeshPolyscopeLabel)
      ->addEdgeScalarQuantity("J", J);
  polyscope::getCurveNetwork(deformedMeshPolyscopeLabel)
      ->addEdgeScalarQuantity("I2", I2);
  polyscope::getCurveNetwork(deformedMeshPolyscopeLabel)
      ->addEdgeScalarQuantity("I3", I3);

  // Specify the callback
  static bool calledOnce = false;
  if (!calledOnce) {
    PolyscopeInterface::addUserCallback([&]() {
      ImGui::PushItemWidth(100);
      static int debugModeStep = mSettings.debugModeStep.has_value()
                                     ? mSettings.debugModeStep.value()
                                     : 0;
      if (ImGui::InputInt("Simulation debug step",
                          &debugModeStep)) {  // set a int variable
        if (debugModeStep != 0) {
          *mSettings.debugModeStep = debugModeStep;
        }
      }
      ImGui::Checkbox("Enable drawing",
                      &mSettings.shouldDraw);  // set a int variable
      ImGui::Text("Number of simulation steps: %zu", mCurrentSimulationStep);

      ImGui::PopItemWidth();
    });
    calledOnce = true;
  }

  if (!screenshotsFolder.empty()) {
    static bool firstDraw = true;
    if (firstDraw) {
      for (const auto& entry :
           std::filesystem::directory_iterator(screenshotsFolder))
        std::filesystem::remove_all(entry.path());
      firstDraw = false;
    }
    polyscope::screenshot(
        std::filesystem::path(screenshotsFolder)
            .append(std::to_string(mCurrentSimulationStep) + ".png")
            .string(),
        false);
  }
  polyscope::show();
}
#endif

void DRMSimulationModel::applyKineticDamping(
    const std::shared_ptr<SimulationJob>& pJob) {
  std::for_each(
#ifdef ENABLE_PARALLEL_DRM
      std::execution::par_unseq,
#endif
      pMesh->nodes._handle->data.begin(), pMesh->nodes._handle->data.end(),
      [&](Node& node) {
        const Vector6d stepDisplacement = node.velocity * Dt;
        node.displacements = node.displacements - stepDisplacement;
      });
  if (!pJob->nodalForcedDisplacements.empty()) {
    applyForcedDisplacements(pJob->nodalForcedDisplacements);
  }
  applyDisplacements(constrainedVertices);
  updateElementalLengths();
  resetVelocities();
  kineticDampingWasAppliedInThisIteration = true;
  ++numOfDampings;
}

void DRMSimulationModel::updateDerivatives() {
  updateNormalDerivatives();
  updateT1Derivatives();
  updateRDerivatives();
  updateT2Derivatives();
  updateT3Derivatives();
}

bool DRMSimulationModel::convergedUsingResidualForcesCriteria() {
  const bool fullfillsAverageResidualForcesNormTerminationCriterion =
      mSettings.threshold_averageResidualToExternalForcesNorm.has_value() &&
      (pMesh->totalResidualForcesNorm / pMesh->VN()) /
              pMesh->totalExternalForcesNorm <
          mSettings.threshold_averageResidualToExternalForcesNorm.value();

  const bool fullfillsResidualForcesNormThreshold =
      mSettings.threshold_residualForcesNorm.has_value() &&
      pMesh->totalResidualForcesNorm < mSettings.threshold_residualForcesNorm;
  const bool converged =
      (fullfillsAverageResidualForcesNormTerminationCriterion ||
       fullfillsResidualForcesNormThreshold);

  if (converged && mSettings.beVerbose) {
    postConvergenceMessage = "Simulation converged.\n";
    if (fullfillsAverageResidualForcesNormTerminationCriterion) {
      postConvergenceMessage +=
          "Converged using average residual forces norm threshold "
          "criterion";
    } else {
      assert(fullfillsResidualForcesNormThreshold);
      postConvergenceMessage +=
          "Converged using residual forces norm threshold "
          "criterion";
    }
  }

  return converged;
}

SimulationResults DRMSimulationModel::executeSimulation(
    const std::shared_ptr<SimulationJob>& pJob,
    const DRMSimulationModel::Settings& settings) {
  beginTime = std::chrono::steady_clock::now();
  reset(pJob, settings);

  if (mSettings.beVerbose) {
    mSettings.printConvergenceCriteria();
  }

  if (mSettings.beVerbose) {
    std::cout << "Executing simulation for mesh " << pMesh->getLabel()
              << " which has " << pMesh->VN() << " nodes and " << pMesh->EN()
              << " elements." << std::endl;
  }

  updateNodalMasses();
  updateNodalExternalForces(pJob->nodalExternalForces, constrainedVertices);
  while (!mSettings.maxDRMIterations.has_value() ||
         mCurrentSimulationStep < mSettings.maxDRMIterations.value()) {
    updateDerivatives();
    updateResidualForces();
    updateNodalMasses();
    updateNodalAccelerations();
    updateNodalVelocities();
    updateNodalDisplacements();
    if (!pJob->nodalForcedDisplacements.empty()) {
      applyForcedDisplacements(pJob->nodalForcedDisplacements);
    }
    applyDisplacements(constrainedVertices);
    updateElementalLengths();
    mCurrentSimulationStep++;
    if (mSettings.beVerbose && mSettings.debugModeStep.has_value() &&
        mCurrentSimulationStep % mSettings.debugModeStep.value() == 0) {
      printCurrentState();
    }

    if ((mSettings.shouldCreatePlots && mSettings.debugModeStep.has_value()) &&
        mCurrentSimulationStep != 0) {
      history.stepPulse(*pMesh);
    }

    if (mSettings.shouldCreatePlots && mSettings.debugModeStep.has_value() &&
        mCurrentSimulationStep % mSettings.debugModeStep.value() == 0) {
      SimulationResultsReporter::createPlot(
          "Number of Steps", "Log of Residual Forces ",
          history.logResidualForces, {}, history.redMarks);
      //      SimulationResultsReporter::createPlot(
      //          "Number of Steps", "Log of Kinetic energy",
      //          history.kineticEnergy, {}, history.redMarks);
    }

#ifdef POLYSCOPE_DEFINED
    if (mSettings.shouldDraw && mSettings.debugModeStep.has_value() &&
        mCurrentSimulationStep % mSettings.debugModeStep.value() == 0) {
      //      std::string saveTo = std::filesystem::current_path()
      //                               .append("Debugging_files")
      //                               .append("Screenshots")
      //                               .string();
      draw(pJob);
    }
#endif

    const bool itIsNotTooEarly =
        numOfDampings > 0 &&
        (pJob->nodalForcedDisplacements.empty() ||
         mCurrentSimulationStep >
             mSettings.gradualForcedDisplacementSteps->endStep) &&
        (!mSettings.initialDistortion.has_value() ||
         mCurrentSimulationStep >
             mSettings.initialDistortion->duration.endStep) &&
        mCurrentSimulationStep > 20;
    if (itIsNotTooEarly) {
      const bool convergedUsingResidualForces =
          convergedUsingResidualForcesCriteria();
      if (convergedUsingResidualForces) {
        break;
      }
    }
    const bool isKineticEnergyPeak =
        pMesh->previousTotalKineticEnergy > pMesh->currentTotalKineticEnergy;
    if (isKineticEnergyPeak) {
      if (itIsNotTooEarly) {
        const bool fullfillsKineticEnergyTerminationCriterion =
            mSettings.threshold_totalTranslationalKineticEnergy.has_value() &&
            pMesh->currentTotalTranslationalKineticEnergy <
                mSettings.threshold_totalTranslationalKineticEnergy.value();
        const bool shouldTerminate = fullfillsKineticEnergyTerminationCriterion;
        if (shouldTerminate) {
          if (mSettings.beVerbose) {
            postConvergenceMessage = "Simulation converged. ";
            postConvergenceMessage +=
                "The kinetic energy of the system was "
                " used as a convergence criterion";
          }
          break;
        }
      }

      applyKineticDamping(pJob);
      Dt *= mSettings.xi;
      if (mSettings.shouldCreatePlots) {
        history.redMarks.push_back(mCurrentSimulationStep);
      }
    }
  }
  auto endTime = std::chrono::steady_clock::now();
  SimulationResults results = computeResults(pJob);
  results.executionTime =
      std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime)
          .count();

  if (mSettings.beVerbose && !postConvergenceMessage.empty()) {
    std::cout << std::endl;
    std::cout << postConvergenceMessage << std::endl;
    std::cout << "Converged in " << results.executionTime / (1000 * 60.0)
              << " minutes." << std::endl;
    printCurrentState();
  }
  if (!mSettings.debugModeStep.has_value() && mSettings.shouldCreatePlots) {
    SimulationResultsReporter reporter;
    reporter.reportResults({results}, "Results", pJob->pMesh->getLabel());
  }

#ifdef POLYSCOPE_DEFINED
  if (mSettings.shouldDraw) {
    pMesh->unregister();
    pJob->pMesh->unregister();
  }
#endif
  return results;
}

void DRMSimulationModel::Settings::save(
    const std::filesystem::path& jsonPath) const {
  const std::filesystem::path jsonFilePath = [&]() {
    if (std::filesystem::is_directory(jsonPath)) {
      return std::filesystem::path(jsonPath).append(jsonDefaultFileName);
    } else {
      return jsonPath;
    }
  }();
  nlohmann::json json;
  json = *this;
  std::filesystem::create_directories(jsonFilePath.parent_path());
  std::ofstream jsonFile(jsonFilePath);
  std::cout << "Saving DRM settings to:" << jsonFilePath << std::endl;
  jsonFile << json;
  jsonFile.close();
}

bool DRMSimulationModel::Settings::load(
    const std::filesystem::path& jsonFilePath) {
  if (!std::filesystem::exists(std::filesystem::path(jsonFilePath))) {
    std::cerr << "The json file does not exist. Json file provided:"
              << jsonFilePath.string() << std::endl;
    std::cout << "Default drm settings will be used" << std::endl;
    return true;
  }

  if (std::filesystem::path(jsonFilePath).extension() != ".json") {
    std::cerr << "A json file is expected as input. The given file has the "
                 "following extension:"
              << std::filesystem::path(jsonFilePath).extension() << std::endl;
    assert(false);
    return false;
  }

  nlohmann::json json;
  std::ifstream ifs(jsonFilePath.string());
  ifs >> json;
  *this = json;

  return true;
}

bool DRMSimulationModel::Settings::hasConvergenceCriterion() const {
  return threshold_residualForcesNorm.has_value() ||
         threshold_averageResidualToExternalForcesNorm.has_value() ||
         threshold_totalTranslationalKineticEnergy;
}

void DRMSimulationModel::Settings::printConvergenceCriteria() const {
  std::cout << "Using the following criteria for convergence:" << std::endl;
  if (threshold_residualForcesNorm.has_value()) {
    std::cout << "Residual forces norm:" << threshold_residualForcesNorm.value()
              << std::endl;
  }
  if (threshold_averageResidualToExternalForcesNorm.has_value()) {
    std::cout << "Average residual forces norm:"
              << threshold_averageResidualToExternalForcesNorm.value()
              << std::endl;
  }
  if (threshold_totalTranslationalKineticEnergy) {
    std::cout << "Total trasnaltional kinetic energy:"
              << threshold_totalTranslationalKineticEnergy.value() << std::endl;
  }
}
