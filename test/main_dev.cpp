#include <Eigen/Dense>
#include <iostream>
#include <Shapely/coordinates.hpp>
#include <Shapely/molecule.hpp>
#include <vector>

int main(int argc, const char * argv[]) {
  using namespace shapely;
  Molecule mol;
  mol.LoadFile("../data/ethanol_b3lyp.isaq");
  
  
  return 0;
}
