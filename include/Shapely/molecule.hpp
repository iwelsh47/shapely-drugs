#pragma once

#include <Shapely/coordinates.hpp>
#include <vector>
#include <stdx/string.hpp>
#include <fstream>
#include <iterator>
#include <map>

namespace shapely {

class Molecule {
public:
  Molecule() = default;
  
  const std::map<stdx::string, double> element_lookup = {
    {"H", 1.008}, {"C", 12.01}, {"O", 15.99}
  };
  
  void LoadFile(const stdx::string& filename) {
    // Load the entire file as raw data
    std::ifstream file(filename);
    file.seekg(0, std::ios::end);
    stdx::string raw;
    raw.resize(file.tellg());
    file.seekg(0, std::ios::beg);
    file.read(raw.data(), raw.size());
    
    // Split into lines
    std::vector<stdx::string> lines;
    raw.split(std::back_inserter(lines), '\n');
    
    // Go through the lines, skipping the first as it is a header
    for (auto it = lines.begin() + 1; it < lines.end(); ++it) {
      // Strip any whitespace
      it->strip_inplace();
      // Ignore blank lines
      if (!it->size()) { continue; }
      // Split into components. Will have element, x, y, z, charge, radius
      std::vector<stdx::string> components;
      it->split(std::back_inserter(components));
      if (components.size() != 6) { continue; } // Skip short/long lines
      atomic_weights.push_back(element_lookup.at(components[0]));
      atoms.NewPoint(components[1].convert<double>(),
                     components[2].convert<double>(),
                     components[3].convert<double>());
      partial_charges.push_back(components[4].convert<double>());
      radii.push_back(components[5].convert<double>());
    }
    
    std::cout << atoms << '\n';
    
    // Correctly lineup the atoms
    CentreInplace(atoms, atomic_weights.data());
    RotateInplace(atoms, PrincipalAxesAlignment(atoms, atomic_weights.data()));
    
    std::cout << atoms << '\n';
  }
  
  void reserve(uint32_t num) {
    atoms.Reserve(num);
    atomic_weights.reserve(num);
    partial_charges.reserve(num);
    radii.reserve(num);
  }
  
public:
  Cartesian atoms;
  std::vector<double> atomic_weights, partial_charges, radii;
};

}
