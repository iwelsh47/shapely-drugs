#include <algorithm>
#include <complex>
#include <filesystem>

#include <Shapely/FFTProcess.hpp>
#include <Shapely/utils.hpp>

namespace shapely {
#define NUM_FFT 243
static const fs::path wisdom_file("/Users/ivanwelsh/GitHub/shapely-drugs/data/fftw_wisdom.dat");

FFTProcess::FFTProcess(const fs::path& target_file, const std::vector<fs::path>& test_files)
: work_r(nullptr), work_c(nullptr), num_fft(NUM_FFT), grid_spacing(0.5), target(target_file) {
  tests.reserve(test_files.size());
  for (const fs::path& file : test_files) { tests.emplace_back(file); }
}

void GenerateWisdom(void* data_r, void* data_c, int32_t N) {
  fftwf_plan forward = fftwf_plan_dft_r2c_3d(N, N, N, (float*)data_r, (fftwf_complex*)data_c, FFTW_PATIENT);
  fftwf_plan reverse = fftwf_plan_dft_c2r_3d(N, N, N, (fftwf_complex*)data_c, (float*)data_r, FFTW_PATIENT);
  fftwf_export_wisdom_to_filename(wisdom_file.native().data());
  fftwf_destroy_plan(forward);
  fftwf_destroy_plan(reverse);
}

void FFTProcess::Initialise() {
  // Figure out what num_fft should be for this instance
  // Target will stay unmoved, so max of limits plus a bit
  CartesianLimits limits = Limits(target.atoms, target.radii.data());
  float max_dist = limits.MaxAbs();

  // Tests will rotate
  for (const Molecule& mol : tests) {
    limits = Limits(mol.atoms, mol.radii.data());
    float max_mol = limits.MaxAbs();
    if (max_mol > max_dist) { max_mol = max_dist; }
  }
  int32_t calcN = (int32_t)std::ceil(max_dist/grid_spacing) + 5;
  num_fft = 4 * calcN + 1; // *2 for pos, neg direction; *2 for grid -> FFT; +1 so is odd so proper centred
  
  std::cout << "Num FFT to use is: " << num_fft << std::endl;
  
  uint64_t nf3 = num_fft * num_fft * num_fft;
  
  // Use the fftw malloc for working data so is aligned correctly for SIMD instructions
  work_r = (float*) fftwf_malloc(sizeof(float) * nf3);
  work_c = (complex*) fftwf_malloc(sizeof(complex) * nf3);
  
  // Zero out working data
  ZeroWorkMemory();
  
  // Wisdom storeage
  if (!fs::exists(wisdom_file)) {
    fs::create_directories(wisdom_file.parent_path());
    GenerateWisdom(work_r, work_c, num_fft);
  }
  fftwf_import_wisdom_from_filename(wisdom_file.native().data());
  
  // Setup the FFTW plans
  using FC = fftwf_complex*;
  GenerateWisdom(work_r, work_c, num_fft);
  forward = fftwf_plan_dft_r2c_3d(num_fft, num_fft, num_fft, work_r, (FC)work_c, FFTW_PATIENT | FFTW_WISDOM_ONLY);
  reverse = fftwf_plan_dft_c2r_3d(num_fft, num_fft, num_fft, (FC)work_c, work_r, FFTW_PATIENT | FFTW_WISDOM_ONLY);
}

void FFTProcess::PerformFFT(bool do_forward) {
  if (do_forward) { fftwf_execute(this->forward); }
  else { fftwf_execute(this->reverse); }
}

void FFTProcess::SaveResult(std::vector<complex>& store) {
  store.resize(num_fft * num_fft * num_fft);
  // Copy result into storage
  memcpy(store.data(), work_c, sizeof(fftwf_complex) * store.size());
}

#define SWITCH_COORDS(dx, dy, dz, sx, sy, sz) \
test.atoms[atom].x() = sx * mol.atoms[atom].dx(); \
test.atoms[atom].y() = sy * mol.atoms[atom].dy(); \
test.atoms[atom].z() = sz * mol.atoms[atom].dz(); \
break

#define SWITCH_COORDS_RAW(output, input, dx, dy, dz, sx, sy, sz) \
output.x() = sx * input.dx(); \
output.y() = sy * input.dy(); \
output.z() = sz * input.dz(); \
break

Score::Score() {}
Score::Score(uint32_t Ni, float d) : N(Ni), delta(d) { }
Score& Score::operator=(const Score &other) {
  assert(other.N == N && other.delta == delta);
  score = other.score;
  index = other.index;
  rotation = other.rotation;
  return *this;
}
Score& Score::operator=(Score &&other)  {
  assert(other.N == N && other.delta == delta);
  score = std::move(other.score);
  index = std::move(other.index);
  rotation = std::move(other.rotation);
  return *this;
}

void FFTProcess::Run(int32_t num_results, std::map<uint32_t, std::vector<Score>>& results) {
  uint64_t nf3 = pow_n<3>(num_fft);
  // Setup the data to run on target
  Discretisation(target, 3.0f, 0.5f);
  // Run and save the FFT to get DFT of target
  PerformFFT(true);
  SaveResult(target_dft);
  ZeroWorkMemory();
  
  // Empty out the results
  results.clear();
  
  // Iterate through the test molecules
  uint32_t mol_idx = 0;
  for (const Molecule& mol : tests) {
    std::cout << "Running " << mol.name << '\n';
    results[mol_idx].reserve(24);
    // Iterate through the 24 simple rotations
    for (uint32_t i = 0; i < 24; ++i) {
      // Make a copy of the molecule and rotate the coordinates according to the simple rotation
      Molecule test = mol;
      for (uint32_t atom = 0; atom < mol.atoms.num_points; ++atom) {
        switch (i) {
          case  0: SWITCH_COORDS(x, y, z,  1,  1,  1);
          case  1: SWITCH_COORDS(x, z, y,  1, -1,  1);
          case  2: SWITCH_COORDS(x, y, z,  1, -1, -1);
          case  3: SWITCH_COORDS(x, z, y,  1,  1, -1);
          case  4: SWITCH_COORDS(x, y, z, -1,  1, -1);
          case  5: SWITCH_COORDS(x, y, z, -1, -1,  1);
          case  6: SWITCH_COORDS(x, z, y, -1,  1,  1);
          case  7: SWITCH_COORDS(x, z, y, -1, -1, -1);
            
          case  8: SWITCH_COORDS(y, x, z,  1, -1,  1);
          case  9: SWITCH_COORDS(y, z, x,  1, -1, -1);
          case 10: SWITCH_COORDS(y, x, z,  1,  1, -1);
          case 11: SWITCH_COORDS(y, z, x,  1,  1,  1);
          case 12: SWITCH_COORDS(y, x, z, -1,  1,  1);
          case 13: SWITCH_COORDS(y, z, x, -1,  1, -1);
          case 14: SWITCH_COORDS(y, z, x, -1, -1,  1);
          case 15: SWITCH_COORDS(y, x, z, -1, -1, -1);
            
          case 16: SWITCH_COORDS(z, y, x,  1,  1, -1);
          case 17: SWITCH_COORDS(z, x, y,  1,  1,  1);
          case 18: SWITCH_COORDS(z, y, x,  1, -1,  1);
          case 19: SWITCH_COORDS(z, x, y,  1, -1, -1);
          case 20: SWITCH_COORDS(z, y, x, -1,  1,  1);
          case 21: SWITCH_COORDS(z, x, y, -1,  1, -1);
          case 22: SWITCH_COORDS(z, x, y, -1, -1,  1);
          case 23: SWITCH_COORDS(z, y, x, -1, -1, -1);
          default: break;
        }
      }
      // Setup the data to run on test
      Discretisation(test, 1.5f, 0.5f);
      PerformFFT(true);
      
      // Perform the convolution: DFT_target* * DFT_test
      for (uint64_t j = 0; j < nf3; ++j) {
        work_out[j] *= std::conj(target_dft[j]);
      }
      
      // Perform the inverse DFT
      PerformFFT(false);
      
      // Get the score
      Score test_score(num_fft, grid_spacing);
      test_score.rotation = i;
      for (uint64_t j = 0; j < nf3; ++j) {
        double res = std::real(work_in[j]);
        if (res > test_score.score) {
          test_score.score = res;
          test_score.index = j;
        }
      }
      test_score.score /= nf3;
      results[mol_idx].emplace_back(std::move(test_score));
      
      // Zero out working data
      ZeroWorkMemory();
    }
    
    // Sort the scores
    std::sort(results[mol_idx].begin(), results[mol_idx].end(), std::greater<Score>());
    // Discard if more than asked for
    if (num_results > 0 && results[mol_idx].size() > num_results) {
      results[mol_idx].resize(num_results);
    }
    
    ++mol_idx;
  }
}

void FFTProcess::PrintTargetXYZ(std::ostream& os) const {
  target.PrintXYZ(os);
}
void FFTProcess::PrintTestXYZ(std::ostream& os, uint32_t idx) const {
  assert(idx < tests.size());
  tests[idx].PrintXYZ(os);
}
void FFTProcess::PrintScoredTestXYZ(std::ostream& os, uint32_t idx, const Score& score) const {
  assert(idx < tests.size());
  Molecule tmp = tests[idx];
  for (uint32_t atom = 0; atom < tmp.atoms.num_points; ++atom) {
    switch (score.GetRotationType()) {
      case  0: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], x, y, z,  1,  1,  1);
      case  1: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], x, z, y,  1, -1,  1);
      case  2: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], x, y, z,  1, -1, -1);
      case  3: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], x, z, y,  1,  1, -1);
      case  4: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], x, y, z, -1,  1, -1);
      case  5: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], x, y, z, -1, -1,  1);
      case  6: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], x, z, y, -1,  1,  1);
      case  7: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], x, z, y, -1, -1, -1);
        
      case  8: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], y, x, z,  1, -1,  1);
      case  9: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], y, z, x,  1, -1, -1);
      case 10: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], y, x, z,  1,  1, -1);
      case 11: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], y, z, x,  1,  1,  1);
      case 12: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], y, x, z, -1,  1,  1);
      case 13: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], y, z, x, -1,  1, -1);
      case 14: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], y, z, x, -1, -1,  1);
      case 15: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], y, x, z, -1, -1, -1);
        
      case 16: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], z, y, x,  1,  1, -1);
      case 17: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], z, x, y,  1,  1,  1);
      case 18: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], z, y, x,  1, -1,  1);
      case 19: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], z, x, y,  1, -1, -1);
      case 20: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], z, y, x, -1,  1,  1);
      case 21: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], z, x, y, -1,  1, -1);
      case 22: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], z, x, y, -1, -1,  1);
      case 23: SWITCH_COORDS_RAW(tmp.atoms[atom], tests[idx].atoms[atom], z, y, x, -1, -1, -1);
      default: break;
    }
    tmp.atoms[atom].x() -= score.GetDeltaX();
    tmp.atoms[atom].y() -= score.GetDeltaY();
    tmp.atoms[atom].z() -= score.GetDeltaZ();
  }
  tmp.PrintXYZ(os);
}

void FFTProcess::Discretisation(const Molecule &mol, float upper, float lower) {
  using namespace std::complex_literals;
  const float grid_length = grid_spacing * num_fft; // Length of the grid from side to side
  const float grid_min = -grid_length / 2; // Min and max coordinate value of the grid
  
  // Go through the atoms in the molecule
  for (uint32_t atom_idx = 0; atom_idx < mol.atoms.num_points; ++atom_idx) {
    const Eigen::RowVector3f pt = mol.atoms[atom_idx];
    const Eigen::RowVector3f mid_pt = (pt.array() - grid_min) / grid_spacing;
    const float r2 = pow_n<2>(mol.radii[atom_idx]); // square of atom radius for cutoff determination
    const float r2_lower = pow_n<2>(lower * mol.radii[atom_idx]);
    const float r2_upper = pow_n<2>(upper * mol.radii[atom_idx]);
    // Determine number of points to iterate in each direction around atom centre
    int32_t delta = std::ceil(upper * mol.radii[atom_idx] / grid_spacing) + 1;
    // Find the nearest index in X direction
    int32_t mid_x = std::round(mid_pt.x());
    // Iterate over X direction
    for (uint32_t ix = std::max(0, mid_x - delta); ix < std::min(num_fft, mid_x + delta); ++ix) {
      Eigen::RowVector3f grid_pt(grid_spacing * ix + grid_min, 0.0, 0.0);
      // Do the same for Y
      int32_t mid_y = std::round(mid_pt.y());
      for (uint32_t iy = std::max(0, mid_y - delta); iy < std::min(num_fft, mid_y + delta); ++iy) {
        grid_pt.y() = grid_spacing * iy + grid_min;
        // And Z
        int32_t mid_z = std::round(mid_pt.z());
        for (uint32_t iz = std::max(0, mid_z - delta); iz < std::min(num_fft, mid_z + delta); ++iz) {
          uint64_t i = num_fft * num_fft * ix + num_fft * iy + iz;
          // Don't do anything if inside an inner core
          if (work_in[i].imag() > 1.5f) { continue; }
          // Finish calculate the grid point position
          grid_pt.z() = grid_spacing * iz + grid_min;
          // Get the distance between the grid point and the atom
          float dist = (pt - grid_pt).squaredNorm();
          // Skip if outside radius
          if (dist > r2_upper) { continue; }
          
          // Check the three distances
          if (dist <= r2_lower) {
            // Inner core. No electrostatics. Set volume to 2 and electrostatics to 0
            work_in[i] = 0.f;// + 2if;
          } else if (dist <= r2) {
            // Core. Electrostatics and volume to 1
            float r = std::sqrt(dist);
            work_in[i] = work_in[i].real() + (float)mol.partial_charges[atom_idx] / r;// + 1if;
//            work_in[i] = 0.f + 1if;
          } else {
            // Outer. Only electrostatics
            float r = std::sqrt(dist);
            work_in[i] += mol.partial_charges[atom_idx] / r;
          }
        }
      }
    }
  }
  
  // Change in inner core values to i
  for (uint64_t i = 0; i < pow_n<3>(num_fft); ++i) {
    if (work_in[i].imag() > 0.f) { work_in[i].imag(5.e-2); }
  }
}

FFTProcess::~FFTProcess() {
  fftwf_free(work_in);
  fftwf_free(work_out);
  fftwf_destroy_plan(forward);
  fftwf_destroy_plan(reverse);
  fftwf_cleanup();
}

void FFTProcess::ZeroWorkMemory() {
  uint64_t nf3 = pow_n<3>(num_fft);
  for (uint64_t i = 0; i < nf3; ++i) {
    work_in[i] = 0.f;
    work_out[i] = 0.f;
  }
}

}
