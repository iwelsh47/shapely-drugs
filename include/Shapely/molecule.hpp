#pragma once

#include <Shapely/coordinates.hpp>
#include <Shapely/serialisation.hpp>
#include <vector>
#include <stdx/string.hpp>
#include <filesystem>

namespace fs = std::filesystem;

namespace shapely {

class PrecalculatedGrid {
public:
  
  std::array<float, 3> origin;
  std::array<uint64_t, 3> delta;
  float spacing;
  std::vector<float> data;
  
  bool empty() const { return data.empty(); }
};

class Molecule {
public:
  Molecule() = default;
  Molecule(const Molecule&) = default;
  Molecule(Molecule&&) = default;
  Molecule& operator=(Molecule&&) = default;
  Molecule& operator=(const Molecule&) = default;
  Molecule(const fs::path& filename);
  
  void LoadFile(const fs::path& filename);
  void LoadBinaryFile(const fs::path& filename);
  
  void reserve(uint32_t num);
  
  template <typename Archive>
  void CEREAL_SAVE_FUNCTION_NAME(Archive& ar, const uint32_t version) const {
    ar(name, atoms, partial_charges, radii, elements);
  }
  
  template <typename Archive>
  void CEREAL_LOAD_FUNCTION_NAME(Archive& ar, const uint32_t version) {
    ar(name, atoms, partial_charges, radii, elements);
  }
  
  void PrintXYZ(std::ostream& os) const;
  
//  double CalculateOverlapVolume(const Molecule& mol) const {
//    
//  }
  
public:
  stdx::string name;
  Cartesian atoms;
  std::vector<float> partial_charges, radii;
  std::vector<stdx::string> elements;
  PrecalculatedGrid grid;
};

}
CEREAL_CLASS_VERSION(shapely::Molecule, 1);
