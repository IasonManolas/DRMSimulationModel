#ifndef UTILITIES_H
#define UTILITIES_H

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <nlohmann/json.hpp>
#include <numeric>
#include <regex>
#include <string_view>

#define GET_VARIABLE_NAME(Variable) (#Variable)
NLOHMANN_JSON_NAMESPACE_BEGIN
template <typename T>
struct adl_serializer<std::optional<T>> {
  static void to_json(json& j, const std::optional<T>& opt) {
    if (opt == std::nullopt) {
      j = nullptr;
    } else {
      j = *opt;  // this will call adl_serializer<T>::to_json which will
                 // find the free function to_json in T's namespace!
    }
  }

  static void from_json(const json& j, std::optional<T>& opt) {
    if (j.is_null()) {
      opt = std::nullopt;
    } else {
      opt = j.get<T>();  // same as above, but with
                         // adl_serializer<T>::from_json
    }
  }
};
NLOHMANN_JSON_NAMESPACE_END

namespace Utilities {

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 2) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

inline bool compareNat(const std::string& a, const std::string& b) {
  if (a.empty())
    return true;
  if (b.empty())
    return false;
  if (std::isdigit(a[0]) && !std::isdigit(b[0]))
    return true;
  if (!std::isdigit(a[0]) && std::isdigit(b[0]))
    return false;
  if (!std::isdigit(a[0]) && !std::isdigit(b[0])) {
    if (std::toupper(a[0]) == std::toupper(b[0]))
      return compareNat(a.substr(1), b.substr(1));
    return (std::toupper(a[0]) < std::toupper(b[0]));
  }

  // Both strings begin with digit --> parse both numbers
  std::istringstream issa(a);
  std::istringstream issb(b);
  int ia, ib;
  issa >> ia;
  issb >> ib;
  if (ia != ib)
    return ia < ib;

  // Numbers are the same --> remove numbers and recurse
  std::string anew, bnew;
  std::getline(issa, anew);
  std::getline(issb, bnew);
  return (compareNat(anew, bnew));
}

inline std::string_view leftTrimSpaces(const std::string_view& str) {
  std::string_view trimmedString = str;
  const auto pos(str.find_first_not_of(" \t\n\r\f\v"));
  trimmedString.remove_prefix(std::min(pos, trimmedString.length()));
  return trimmedString;
}

inline std::string_view rightTrimSpaces(const std::string_view& str) {
  std::string_view trimmedString = str;
  const auto pos(trimmedString.find_last_not_of(" \t\n\r\f\v"));
  trimmedString.remove_suffix(
      std::min(trimmedString.length() - pos - 1, trimmedString.length()));
  return trimmedString;
}

inline std::string_view trimLeftAndRightSpaces(std::string_view str) {
  std::string_view trimmedString = str;
  trimmedString = leftTrimSpaces(trimmedString);
  trimmedString = rightTrimSpaces(trimmedString);
  return trimmedString;
}

template <typename InputIt>
inline void normalize(InputIt itBegin, InputIt itEnd) {
  const auto squaredSumOfElements = std::accumulate(
      itBegin, itEnd, 0.0,
      [](const auto& sum, const auto& el) { return sum + el * el; });
  assert(squaredSumOfElements != 0);
  std::transform(itBegin, itEnd, itBegin, [&](auto& element) {
    return element / std::sqrt(squaredSumOfElements);
  });
}

inline std::vector<std::string> split(const std::string& text,
                                      std::string delim) {
  std::vector<std::string> vec;
  size_t pos = 0, prevPos = 0;
  while (1) {
    pos = text.find(delim, prevPos);
    if (pos == std::string::npos) {
      vec.push_back(text.substr(prevPos));
      return vec;
    }

    vec.push_back(text.substr(prevPos, pos - prevPos));
    prevPos = pos + delim.length();
  }
}

inline std::string toString(const std::vector<std::vector<int>>& vv) {
  std::string s;
  s.append("{");
  for (const std::vector<int>& v : vv) {
    s.append("{");
    for (const int& i : v) {
      s.append(std::to_string(i) + ",");
    }
    s.pop_back();
    s.append("}");
  }
  s.append("}");
  return s;
}
inline void parseIntegers(const std::string& str, std::vector<size_t>& result) {
  typedef std::regex_iterator<std::string::const_iterator> re_iterator;
  typedef re_iterator::value_type re_iterated;

  std::regex re("(\\d+)");

  re_iterator rit(str.begin(), str.end(), re);
  re_iterator rend;

  std::transform(rit, rend, std::back_inserter(result),
                 [](const re_iterated& it) { return std::stoi(it[1]); });
}

inline void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove) {
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
        matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

  matrix.conservativeResize(numRows, numCols);
}

inline void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove) {
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols() - 1;

  if (colToRemove < numCols)
    matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
        matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

  matrix.conservativeResize(numRows, numCols);
}

inline std::filesystem::path getFilepathWithExtension(
    const std::filesystem::path& folderPath,
    const std::string& extension) {
  assert(std::filesystem::exists(folderPath));
  for (const std::filesystem::directory_entry& dirEntry :
       std::filesystem::directory_iterator(folderPath)) {
    if (dirEntry.is_regular_file() &&
        std::filesystem::path(dirEntry).extension() == extension) {
      return std::filesystem::path(dirEntry);
    }
  }

  return "";
}

void createPlot(const std::string& xLabel,
                const std::string& yLabel,
                const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& markerSizes,
                const std::vector<double>& c,
                const std::string& saveTo = {});

}  // namespace Utilities

#ifdef POLYSCOPE_DEFINED
using RGBColor = std::array<float, 3>;

namespace PolyscopeInterface {
inline struct GlobalPolyscopeData {
  std::vector<std::function<void()>> userCallbacks;
} globalPolyscopeData;

void mainCallback();

void addUserCallback(const std::function<void()>& userCallback);

void deinitPolyscope();

void init();
using PolyscopeLabel = std::string;
std::pair<PolyscopeLabel, size_t> getSelection();

void registerWorldAxes();
}  // namespace PolyscopeInterface

#endif

template <typename T1, typename T2>
void constructInverseMap(const T1& map, T2& oppositeMap) {
  assert(!map.empty());
  oppositeMap.clear();
  for (const auto& mapIt : map) {
    oppositeMap[mapIt.second] = mapIt.first;
  }
}

template <typename T>
std::string toString(const T& v) {
  return "(" + std::to_string(v[0]) + "," + std::to_string(v[1]) + "," +
         std::to_string(v[2]) + ")";
}

template <typename T>
size_t computeHashUnordered(const std::vector<T>& v) {
  size_t hash = 0;
  for (const auto& el : v) {
    hash += std::hash<T>{}(el);
  }
  return hash;
}

inline size_t computeHashOrdered(const std::vector<int>& v) {
  std::string elementsString;
  for (const auto& el : v) {
    elementsString += std::to_string(el);
  }

  return std::hash<std::string>{}(elementsString);
}

struct Vector6d : public std::array<double, 6> {
  Vector6d() {
    for (size_t i = 0; i < 6; i++) {
      this->operator[](i) = 0;
    }
  }

  Vector6d(const std::vector<double>& v) {
    assert(v.size() == 6);
    std::copy(v.begin(), v.end(), this->begin());
  }

  Vector6d(const double& d) {
    for (size_t i = 0; i < 6; i++) {
      this->operator[](i) = d;
    }
  }

  Vector6d(const std::array<double, 6>& arr) : std::array<double, 6>(arr) {}

  Vector6d(const std::initializer_list<double>& initList) {
    std::copy(initList.begin(), initList.end(), std::begin(*this));
  }

  Vector6d operator*(const double& d) const {
    Vector6d result;
    for (size_t i = 0; i < 6; i++) {
      result[i] = this->operator[](i) * d;
    }
    return result;
  }

  Vector6d operator*(const Vector6d& v) const {
    Vector6d result;
    for (size_t i = 0; i < 6; i++) {
      result[i] = this->operator[](i) * v[i];
    }
    return result;
  }

  Vector6d operator/(const double& d) const {
    Vector6d result;
    for (size_t i = 0; i < 6; i++) {
      result[i] = this->operator[](i) / d;
    }
    return result;
  }

  Vector6d operator/(const Vector6d& v) const {
    Vector6d result;
    for (size_t i = 0; i < 6; i++) {
      result[i] = this->operator[](i) / v[i];
    }
    return result;
  }

  Vector6d operator+(const Vector6d& v) const {
    Vector6d result;
    for (size_t i = 0; i < 6; i++) {
      result[i] = this->operator[](i) + v[i];
    }
    return result;
  }

  Vector6d operator-(const Vector6d& v) const {
    Vector6d result;
    for (size_t i = 0; i < 6; i++) {
      result[i] = this->operator[](i) - v[i];
    }
    return result;
  }

  Vector6d inverted() const {
    Vector6d result;
    for (size_t i = 0; i < 6; i++) {
      assert(this->operator[](i) != 0);
      result[i] = 1 / this->operator[](i);
    }
    return result;
  }
  bool isZero() const {
    for (size_t i = 0; i < 6; i++) {
      if (this->operator[](i) != 0)
        return false;
    }
    return true;
  }

  double squaredNorm() const {
    double squaredNorm = 0;
    std::for_each(this->begin(), std::end(*this),
                  [&](const double& v) { squaredNorm += pow(v, 2); });
    return squaredNorm;
  }

  double norm() const { return sqrt(squaredNorm()); }

  bool isFinite() const {
    return std::any_of(std::begin(*this), std::end(*this), [](const double& v) {
      if (!std::isfinite(v)) {
        return false;
      }
      return true;
    });
  }

  Eigen::Vector3d getTranslation() const {
    return Eigen::Vector3d(this->operator[](0), this->operator[](1),
                           this->operator[](2));
  }

  void setTranslation(const Eigen::Vector3d& v) {
    this->operator[](0) = v(0);
    this->operator[](1) = v(1);
    this->operator[](2) = v(2);
  }

  void setRotation(const Eigen::Vector3d& v) {
    this->operator[](3) = v(0);
    this->operator[](4) = v(1);
    this->operator[](5) = v(2);
  }

  Eigen::Vector3d getRotation() const {
    return Eigen::Vector3d(this->operator[](3), this->operator[](4),
                           this->operator[](5));
  }

  std::string toString() const {
    std::string s;
    for (int i = 0; i < 6; i++) {
      s.append(Utilities::to_string_with_precision(this->operator[](i), 10) +
               ",");
    }
    s.pop_back();
    return s;
  }

  //  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Vector6d)
};

NLOHMANN_JSON_NAMESPACE_BEGIN
template <>
struct adl_serializer<Vector6d> {
  static void to_json(nlohmann::json& j, const Vector6d& v) {
    j = static_cast<std::array<double, 6>>(v);
  }
  static void from_json(const nlohmann::json& j, Vector6d& v) {
    v = j.get<std::array<double, 6>>();
  }
};
NLOHMANN_JSON_NAMESPACE_END

namespace Utilities {
inline Eigen::MatrixXd toEigenMatrix(const std::vector<Vector6d>& v) {
  Eigen::MatrixXd m(v.size(), 6);

  for (size_t vi = 0; vi < v.size(); vi++) {
    const Vector6d& vec = v[vi];
    for (size_t i = 0; i < 6; i++) {
      m(vi, i) = vec[i];
    }
  }

  return m;
}

inline std::vector<Vector6d> fromEigenMatrix(const Eigen::MatrixXd& m) {
  std::vector<Vector6d> v(m.rows());

  for (size_t vi = 0; vi < m.rows(); vi++) {
    const Eigen::RowVectorXd& row = m.row(vi);
    for (size_t i = 0; i < 6; i++) {
      v[vi][i] = row(i);
    }
  }

  return v;
}
}  // namespace Utilities

#endif  // UTILITIES_H
