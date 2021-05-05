//
//  molecule.cpp
//  shapely_tester
//
//  Created by Ivan Welsh on 11/02/21.
//

#include <fstream>
#include <filesystem>

#include <Shapely/molecule.hpp>
#include <Shapely/serialisation.hpp>

static const std::vector<stdx::string> ElementSymbols = {
  "XX", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
  "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
  "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
  "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
  "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
  "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
  "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};

namespace shapely {
struct AtomSerialData {
  std::array<float, 3> coord;
  float partial_charge, radius;
  int16_t element;
  
  template <class Archive> void CEREAL_SERIALIZE_FUNCTION_NAME(Archive &ar) {
    ar(coord, partial_charge, radius, element);
  }
};

Molecule::Molecule(const fs::path& filename) {
  if (filename.extension() == ".isaq") { LoadFile(filename); }
  else { LoadBinaryFile(filename); }
}

void Molecule::LoadFile(const fs::path &filename) {
  // Must be an isaq file
  if (filename.extension() != ".isaq") {
    throw std::runtime_error("Incorrect file type");
  }
  
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
    atoms.NewPoint(components[1].convert<float>(),
                   components[2].convert<float>(),
                   components[3].convert<float>());
    partial_charges.push_back(components[4].convert<float>());
    radii.push_back(components[5].convert<float>());
    stdx::string::size_type counts = components[0].find_first_of(stdx::string::digits());
    elements.push_back(components[0].substr(0, counts));
  }
  
  // Correctly lineup the atoms
  // Cube the radii for weighting things
  std::vector<float> r3 = radii;
  for (float& r : r3) { r *= r * r; }
  
  CentreInplace(atoms, r3.data());
  RotateInplace(atoms, PrincipalAxesAlignment(atoms, r3.data()));
  
  name = filename.stem();
}

void Molecule::LoadBinaryFile(const fs::path &filename) {
  // Must be a density file
  if (filename.extension() != ".density") {
    throw std::runtime_error("Incorrect file type");
  }
  
  std::vector<AtomSerialData> atom_data;
  {
    char throwaway_byte;
    using Loader = cereal::BinaryInputArchive;
    std::ifstream is(filename, std::ios::binary);
    Loader load(is);
    load(throwaway_byte, atom_data);
    load(grid.delta, grid.origin, grid.spacing, grid.data);
  }
  
  // Generate the molecule from the atom data
  for (AtomSerialData& atom : atom_data) {
    atoms.NewPoint(atom.coord[0], atom.coord[1], atom.coord[2]);
    partial_charges.push_back(atom.partial_charge);
    radii.push_back(atom.radius);
    elements.push_back(ElementSymbols[atom.element]);
  }
  
  name = filename.stem();
}

void Molecule::reserve(uint32_t num)  {
  atoms.Reserve(num);
  partial_charges.reserve(num);
  radii.reserve(num);
  elements.reserve(num);
}

void Molecule::PrintXYZ(std::ostream &os) const {
  os << elements.size() << '\n'; // Number of atoms
  os << name << '\n'; // File name
  os << std::setprecision(5) << std::fixed;
  for (uint32_t i = 0; i < elements.size(); ++i) {
    os << std::setw(3) << elements[i] << " "
       << std::setw(15) << (atoms[i].x() * 0.5291772109)<< " "
       << std::setw(15) << (atoms[i].y() * 0.5291772109)<< " "
       << std::setw(15) << (atoms[i].z() * 0.5291772109)<< '\n';
  }
}

}
