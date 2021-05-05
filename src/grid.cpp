#include <Shapely/grid.hpp>
#include <Shapely/molecule.hpp>

namespace shapely {
void Grid::PrepGrid(const Molecule& mol) {
 // Resize everything
 esp.resize(esp_npoints, 0);
 for (VGrid& g : volume) { g.resize(vol_npoints); }
 name = mol.name;
 
 // ESP stuffs
 // X dimension
 for (int64_t ix = -esp_dpoints; ix <= esp_dpoints; ++ix) {
   Eigen::RowVector3f pt(0.0, 0.0, 0.0);
   pt.x() = esp_spacing * ix;
   
   // Y dimension
   for (int64_t iy = -esp_dpoints; iy <= esp_dpoints; ++iy) {
     pt.y() = esp_spacing * iy;
     
     // Z dimension
     for (int64_t iz = -esp_dpoints; iz <= esp_dpoints; ++iz) {
       pt.z() = esp_spacing * iz;
       
       // Atom testing
       uint64_t iv = GetIndex<false>(ix, iy, iz, 0);
       bool is_vol = false;
       std::vector<float> pt_contributions(mol.atoms.NumPoints(), 0);
       for (uint32_t atom = 0; atom < mol.atoms.NumPoints(); ++atom) {
         float dist = (mol.atoms[atom] - pt).norm();
         pt_contributions[atom] = mol.partial_charges[atom] / dist;
         if (!is_vol && dist < mol.radii[atom]) { is_vol = true; }
         if (is_vol) { break; }
       }
       if (!is_vol) { esp[iv] = accurate_sum(pt_contributions); }
     }
   }
 }
 
 PopulateVolumes(mol);
}

void Grid::PopulateVolumes(const Molecule& mol)  {
  CartesianLimits limits = Limits(mol.atoms, mol.radii.data());
  int64_t lowx, lowy, lowz, highx, highy, highz;
  if (limits.min_x < 0) { lowx = (int64_t)std::round((limits.min_x - vol_spacing) / vol_spacing); }
  else { lowx = (int64_t)std::round((limits.min_x + vol_spacing) / vol_spacing); }
  if (limits.max_x > 0) { highx = (int64_t)std::round((limits.max_x + vol_spacing) / vol_spacing); }
  else { highx = (int64_t)std::round((limits.max_x - vol_spacing) / vol_spacing); }
  
  if (limits.min_y < 0) { lowy = (int64_t)std::round((limits.min_y - vol_spacing) / vol_spacing); }
  else { lowy = (int64_t)std::round((limits.min_y + vol_spacing) / vol_spacing); }
  if (limits.max_y > 0) { highy = (int64_t)std::round((limits.max_y + vol_spacing) / vol_spacing); }
  else { highy = (int64_t)std::round((limits.max_y - vol_spacing) / vol_spacing); }
  
  if (limits.min_z < 0) { lowz = (int64_t)std::round((limits.min_z - vol_spacing) / vol_spacing); }
  else { lowz = (int64_t)std::round((limits.min_z + vol_spacing) / vol_spacing); }
  if (limits.max_z > 0) { highz = (int64_t)std::round((limits.max_z + vol_spacing) / vol_spacing); }
  else { highz = (int64_t)std::round((limits.max_z - vol_spacing) / vol_spacing); }
  
  for (int64_t ix = lowx - 1; ix <= highx + 1; ++ix) {
    Eigen::RowVector3f pt(0.0, 0.0, 0.0);
    pt.x() = vol_spacing * ix;
    for (int64_t iy = lowy - 1; iy <= highy + 1; ++iy) {
      pt.y() = vol_spacing * iy;
      for (int64_t iz = lowz - 1; iz <= highz + 1; ++iz) {
        pt.z() = vol_spacing * iz;
        for (uint32_t atom = 0; atom < mol.atoms.NumPoints(); ++atom) {
          float dist = (mol.atoms[atom] - pt).norm();
          if (dist < mol.radii[atom]) {
            for (uint8_t style = 0; style < volume.size(); ++style) {
              volume[style].set(GetIndex(ix, iy, iz, style));
            }
            break;
          }
        }
      }
    }
  }
  V = volume[0].count();
}

std::pair<Overlap, uint8_t> Grid::CalculateBestOverlap(const Grid& other) const  {
  double best_score = 0.0;
  uint8_t best_style = 0;
  for (uint8_t style = 0; style < volume.size(); ++style) {
    double score = std::min((double)V / (double)other.V, (double)other.V / (double)V);
    score *= (double)(volume[0] & other.volume[style]).count();
    if (score > best_score) {
      best_score = score;
      best_style = style;
    }
  }
  return std::make_pair(CalculateOverlap(other, best_style), best_style);
}

Overlap Grid::CalculateOverlap(const Grid& other, uint8_t style) const  {
  assert(style < volume.size() && "Unsupported style");
  Overlap res;
  res.VA = V;
  res.VB = other.V;
  VGrid intersection = (volume[0] & other.volume[style]);
  res.AandB = intersection.count();
  res.AplusB = (volume[0] | other.volume[style]).count();
  res.AminusB = (volume[0] - other.volume[style]).count();
  res.BminusA = (other.volume[style] - volume[0]).count();
  
  for (int64_t ix = -esp_dpoints; ix <= esp_dpoints; ++ix) {
    for (int64_t iy = -esp_dpoints; iy <= esp_dpoints; ++iy) {
      for (int64_t iz = -esp_dpoints; iz <= esp_dpoints; ++iz) {
        uint64_t ime = GetIndex<false>(ix, iy, iz, 0);
        float esp_me = esp[ime];
        float esp_ot = other.esp[GetIndex<false>(ix, iy, iz, style)];
        if (esp_me == 0.f && esp_ot != 0.f) {
          res.esp_Vb += std::fabs(esp_me - esp_ot);
        } else if (esp_me != 0.f && esp_ot == 0.f) {
          res.esp_Va += std::fabs(esp_me - esp_ot);
        } else if (esp_me != 0.f && esp_ot != 0.f) {
          res.esp_outerV += std::fabs(esp_me - esp_ot);
        }
      }
    }
  }
  
  return res;
}

void Grid::PrintXYZ(std::ostream& os, const CartesianLimits& limits, uint8_t style) const  {
  os << std::setprecision(5) << std::fixed;
  struct Position{ float x, y, z;
    Position(float ix, float iy, float iz) : x(ix * 0.5291772109), y(iy * 0.5291772109), z(iz * 0.5291772109) { }
  };
  std::vector<Position> pts; pts.reserve(volume[0].count());
  for (int64_t ix = -vol_dpoints; ix <= vol_dpoints; ++ix) {
    float x = vol_spacing * ix;
    if (x < limits.min_x - 1 || x > limits.max_x + 1) { continue; }
    
    for (int64_t iy = -vol_dpoints; iy <= vol_dpoints; ++iy) {
      float y = vol_spacing * iy;
      if (y < limits.min_y || y > limits.max_y) { continue; }
    
      for (int64_t iz = -vol_dpoints; iz <= vol_dpoints; ++iz) {
        float z = vol_spacing * iz;
        if (z < limits.min_z || z > limits.max_z) { continue; }
        if (volume[0].test(GetIndex<true>(ix, iy, iz, style))) {
          pts.emplace_back(x, y, z);
        }
      }
    }
  }
  os << pts.size() << '\n' << name << '\n';
  for (Position& pt : pts) {
    os << " He " << std::setw(15) << pt.x << " " << std::setw(15) << pt.y << " " << std::setw(15) << pt.z << '\n';
  }
}

void Grid::Compress() {
  using sz_t = VGrid::size_type;
  for (uint64_t idx = 0; idx < volume.size(); ++idx) {
    if (volume[idx].empty()) { continue; }
    compressed_volume[idx] = volume[idx];
    volume[idx].clear();
  }
}

void Grid::Expand() {
  for (uint64_t idx = 0; idx < volume.size(); ++idx) {
    if (volume[idx].empty()) {
      volume[idx] = (VGrid)compressed_volume[idx];
      compressed_volume[idx].clear();
    }
  }
}

}
