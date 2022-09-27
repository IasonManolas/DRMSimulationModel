#ifndef BEAM_HPP
#define BEAM_HPP
#include <assert.h>
#include <cmath>
#include <iostream>
#include <string>
#include "nlohmann/json.hpp"

struct BeamDimensions {
  double dim1, dim2;
  double A{0};  // cross sectional area

  struct MomentsOfInertia {
    double I2{0};  // second moment of inertia
    double I3{0};  // third moment of inertia
    double J{0};   // torsional constant (polar moment of inertia)
  } inertia;

 public:
  double getDim1() const;
  double getDim2() const;

  virtual void to_json(nlohmann::json& j, const BeamDimensions& beamDim) const;

  virtual void from_json(const nlohmann::json& j,
                         BeamDimensions& beamDim) const;
  virtual float getDrawingRadius() const = 0;
};

inline double BeamDimensions::getDim2() const {
  return dim2;
}

inline void BeamDimensions::to_json(nlohmann::json& j,
                                    const BeamDimensions& beamDim) const {
  j = nlohmann::json{{"BeamDimension.dim1", beamDim.dim1},
                     {"BeamDimension.dim2", beamDim.dim2}};
}

inline void BeamDimensions::from_json(const nlohmann::json& j,
                                      BeamDimensions& beamDim) const {
  const std::string jsonKey_dim1 = "BeamDimension.dim1";
  const std::string jsonKey_dim2 = "BeamDimension.dim2";
  if (j.contains(jsonKey_dim1)) {
    j.at(jsonKey_dim1).get_to(beamDim.dim1);
  }
  if (j.contains(jsonKey_dim2)) {
    j.at(jsonKey_dim2).get_to(beamDim.dim2);
  }
}

inline double BeamDimensions::getDim1() const {
  return dim1;
}

struct RectangularBeamDimensions : public BeamDimensions {
  inline static std::string name{"Rectangular"};
  inline const static double defaultSize = 0.002;

  RectangularBeamDimensions(const double& width, const double& height) {
    assert(width > 0 && height > 0);
    dim1 = width;
    dim2 = height;
    updateProperties();
  }
  RectangularBeamDimensions() {
    dim1 = defaultSize;
    dim2 = defaultSize;
    updateProperties();
  }

  std::string toString() const {
    return std::string("w=") + std::to_string(dim1) + std::string(" h=") +
           std::to_string(dim2);
  }

  void updateProperties() {
    A = dim1 * dim2;
    inertia.I2 = dim1 * std::pow(dim2, 3) / 12;
    inertia.I3 = dim2 * std::pow(dim1, 3) / 12;
    //    inertia.J = inertia.I2 + inertia.I3;
    const double t = std::min(dim1, dim2);
    const double b = std::max(dim1, dim2);
    inertia.J = b * pow(t, 3) *
                ((1.0 / 3.0) -
                 0.210 * (t / b) * (1.0 - (1.0 / 12.0) * pow((t / b), 4)));
  }
  static void computeMomentsOfInertia(
      const RectangularBeamDimensions& dimensions,
      MomentsOfInertia& inertia);
  double getWidth() const { return dim1; }
  double getHeight() const { return dim2; }
  float getDrawingRadius() const override;
};

inline float RectangularBeamDimensions::getDrawingRadius() const {
  return std::sqrt(dim1 * dim1 + dim2 * dim2) / 2;
}

struct CylindricalBeamDimensions : public BeamDimensions {
  inline static std::string name{"Cylindrical"};
  // https://www.engineeringtoolbox.com/area-moment-inertia-d_1328.html
  CylindricalBeamDimensions(const double& outsideDiameter,
                            const double& insideDiameter) {
    assert(outsideDiameter > 0 && insideDiameter > 0 &&
           outsideDiameter > insideDiameter);
    dim1 = insideDiameter;
    dim2 = outsideDiameter;
    updateProperties();
  }
  CylindricalBeamDimensions() {
    dim1 = 0.026;
    dim2 = 0.03;
    updateProperties();
  }
  void updateProperties() {
    A = M_PI *
        (std::pow(getOutterDiameter(), 2) - std::pow(getInnerDiameter(), 2)) /
        4;
    inertia.I2 =
        M_PI *
        (std::pow(getOutterDiameter(), 4) - std::pow(getInnerDiameter(), 4)) /
        64;
    inertia.I3 = inertia.I2;
    //    inertia.J = inertia.I2 + inertia.I3;
    inertia.J = inertia.I2 + inertia.I3;  // this is not the torsion constant
                                          // for circular tube cross section
  }
  double getInnerDiameter() const { return dim1; }
  double getOutterDiameter() const { return dim2; }
  float getDrawingRadius() const override;
};

inline float CylindricalBeamDimensions::getDrawingRadius() const {
  return getOutterDiameter();
}

struct ElementMaterial {
  double poissonsRatio;
  double youngsModulus;
  double G;
  ElementMaterial(const double& poissonsRatio, const double& youngsModulus)
      : poissonsRatio(poissonsRatio), youngsModulus(youngsModulus) {
    assert(poissonsRatio <= 0.5 && poissonsRatio >= -1);
    updateShearModulus();
  }
  ElementMaterial() : poissonsRatio(0.3), youngsModulus(1e9) {
    updateShearModulus();
  }
  std::string toString() const {
    return std::string("Material:") + std::string("\nPoisson's ratio=") +
           std::to_string(poissonsRatio) +
           std::string("\nYoung's Modulus(GPa)=") +
           std::to_string(youngsModulus / 1e9);
  }
  void updateShearModulus() { G = youngsModulus / (2 * (1 + poissonsRatio)); }
};

#endif  // BEAM_HPP
