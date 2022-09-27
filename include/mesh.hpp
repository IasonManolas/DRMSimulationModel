#ifndef MESH_HPP
#define MESH_HPP

#include <filesystem>
#include <string>

class Mesh
{
protected:
    std::string label{"empty_label"};

public:
  virtual ~Mesh() = default;
  virtual bool load(const std::filesystem::path &meshFilePath) { return false; }
  virtual bool save(const std::filesystem::path &meshFilePath) { return false; }

  std::string getLabel() const;
  void setLabel(const std::string &newLabel);
  void prependToLabel(const std::string &text);
  void appendToLabel(const std::string &text);
};

inline std::string Mesh::getLabel() const { return label; }

inline void Mesh::setLabel(const std::string &newLabel) { label = newLabel; }

inline void Mesh::prependToLabel(const std::string &text)
{
    label = text + label;
}

inline void Mesh::appendToLabel(const std::string &text)
{
    label = label + text;
}

#endif // MESH_HPP
