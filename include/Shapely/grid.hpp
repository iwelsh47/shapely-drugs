#pragma once

#include <vector>

#include <sul/dynamic_bitset.hpp>
#include <Shapely/compressed_bitset.hpp>
#include <Shapely/coordinates.hpp>
#include <Shapely/molecule.hpp>
#include <Shapely/utils.hpp>

namespace shapely {
struct Overlap {
  uint64_t AandB = 0, AplusB = 0, AminusB = 0, BminusA = 0, VA = 0, VB = 0;
  // ESP score with points not in either volume, in A not B, in B not A
  float esp_outerV = 0.f, esp_Va = 0.f, esp_Vb = 0.f;
};

class Grid {
public:
  static constexpr float esp_spacing = 1.5;
  static constexpr int64_t esp_dpoints = 35;
  static constexpr int64_t esp_apoints = 2 * esp_dpoints + 1;
  static constexpr int64_t esp_npoints = pow_n<3>(esp_apoints);
  
  static constexpr float vol_spacing = 0.25;
  static constexpr int64_t vol_dpoints = 200;
  static constexpr int64_t vol_apoints = 2 * vol_dpoints + 1;
  static constexpr int64_t vol_npoints = pow_n<3>(vol_apoints);
  
  using ESPGrid = std::vector<float>;
  using VGrid = sul::dynamic_bitset<>;
  using CompressedUnit = std::pair<VGrid::size_type, VGrid::size_type>;
  
  // IDX = (ix + dp) + ap ((iy + dp) + ap (iz + dp))
  template <bool is_vol = true>
  static uint64_t GetIndex(int64_t ix, int64_t iy, int64_t iz, uint8_t style = 0) {
    static constexpr int64_t dp = is_vol ? vol_dpoints : esp_dpoints;
    static constexpr int64_t ap = is_vol ? vol_apoints : esp_apoints;
    switch (style) {
      case 0: return ( ix + dp) + ap * (( iy + dp) + ap * ( iz + dp));
      case 1: return ( ix + dp) + ap * ((-iz + dp) + ap * ( iy + dp));
      case 2: return ( ix + dp) + ap * ((-iy + dp) + ap * (-iz + dp));
      case 3: return ( ix + dp) + ap * (( iz + dp) + ap * (-iy + dp));
      case 4: return (-ix + dp) + ap * (( iy + dp) + ap * (-iz + dp));
      case 5: return (-ix + dp) + ap * ((-iy + dp) + ap * ( iz + dp));
      case 6: return (-ix + dp) + ap * (( iz + dp) + ap * ( iy + dp));
      case 7: return (-ix + dp) + ap * ((-iz + dp) + ap * (-iy + dp));
        
      case 8:  return ( iy + dp) + ap * ((-ix + dp) + ap * ( iz + dp));
      case 9:  return ( iy + dp) + ap * ((-iz + dp) + ap * (-ix + dp));
      case 10: return ( iy + dp) + ap * (( ix + dp) + ap * (-iz + dp));
      case 11: return ( iy + dp) + ap * (( iz + dp) + ap * ( ix + dp));
      case 12: return (-iy + dp) + ap * (( ix + dp) + ap * ( iz + dp));
      case 13: return (-iy + dp) + ap * (( iz + dp) + ap * (-ix + dp));
      case 14: return (-iy + dp) + ap * ((-iz + dp) + ap * ( ix + dp));
      case 15: return (-iy + dp) + ap * ((-ix + dp) + ap * (-iz + dp));
        
      case 16: return ( iz + dp) + ap * (( iy + dp) + ap * (-ix + dp));
      case 17: return ( iz + dp) + ap * (( ix + dp) + ap * ( iy + dp));
      case 18: return ( iz + dp) + ap * ((-iy + dp) + ap * ( ix + dp));
      case 19: return ( iz + dp) + ap * ((-ix + dp) + ap * (-iy + dp));
      case 20: return (-iz + dp) + ap * (( iy + dp) + ap * ( ix + dp));
      case 21: return (-iz + dp) + ap * (( ix + dp) + ap * (-iy + dp));
      case 22: return (-iz + dp) + ap * ((-ix + dp) + ap * ( iy + dp));
      case 23: return (-iz + dp) + ap * ((-iy + dp) + ap * (-ix + dp));
        
      default: throw std::out_of_range("Unknown style");
    }
  }
  
  Grid() = default;
  Grid(Grid&& mol) = default;
  Grid(const Grid&) = default;
  Grid(const Molecule& mol) { PrepGrid(mol); }
  Grid& operator=(Grid&& mol) = default;
  Grid& operator=(const Grid&) = default;
  
  void PrepGrid(const Molecule& mol);
  
  void PopulateVolumes(const Molecule& mol);
  
  std::pair<Overlap, uint8_t> CalculateBestOverlap(const Grid& other) const;
  
  Overlap CalculateOverlap(const Grid& other, uint8_t style = 0) const;
  
  template <class Archive>
  void CEREAL_SAVE_FUNCTION_NAME(Archive& ar, const uint32_t version) const {
    ar(name, V, volume, compressed_volume, esp);
  }
  
  template <class Archive>
  void CEREAL_LOAD_FUNCTION_NAME(Archive& ar, const uint32_t version) {
    ar(name, V, volume, compressed_volume, esp);
  }
  
  void PrintXYZ(std::ostream& os, const CartesianLimits& limits, uint8_t style) const;
  
  void Compress();
  
  void Expand();
  
public:
  stdx::string name;
  uint64_t V = 0;
  std::array<VGrid, 24> volume;
  std::array<CompressedBitset, 24> compressed_volume;
  ESPGrid esp;
};
}
