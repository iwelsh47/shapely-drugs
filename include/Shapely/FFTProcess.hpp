#pragma once

#include <complex>
#include <filesystem>
#include <map>
#include <memory>
#include <vector>

#include <fftw3.h>
#include <Shapely/coordinates.hpp>
#include <Shapely/molecule.hpp>

namespace fs = std::filesystem;

namespace shapely {
class Score;

class FFTProcess {
private:
  using complex = std::complex<float>;
  fftwf_plan forward, reverse;
  int32_t num_fft;  // Number of grid points in the FFT in one axis. Twice the number as in the molecule grid?
  float grid_spacing; // Distance between grid points in bohr
  
  std::vector<complex> target_dft;
  float* work_r;
  complex* work_c;
  
  Molecule target;
  std::vector<Molecule> tests;
  
public:
  FFTProcess(const fs::path& target_file,
             const std::vector<fs::path>& test_files);
  
  void Initialise();
  inline std::map<uint32_t, std::vector<Score>> Run(int32_t num_results) {
    std::map<uint32_t, std::vector<Score>> results;
    Run(num_results, results);
    return results;
  }
  void Run(int32_t num_results, std::map<uint32_t, std::vector<Score>>& results);
  
  void PrintTargetXYZ(std::ostream& os) const;
  void PrintTestXYZ(std::ostream& os, uint32_t idx) const;
  void PrintScoredTestXYZ(std::ostream& os, uint32_t idx, const Score& score) const;
  
  ~FFTProcess();
private:
  void PerformFFT(bool forward);
  void SaveResult(std::vector<complex>& store);
  // upper_radius cutoff fraction for electrostatics maximum radius
  // lower radius cutoff fraction for electrostatics minimum radius
  void Discretisation(const Molecule& mol, float upper_radius, float lower_radius);
  void ZeroWorkMemory();
  
};

class Score {
private:
  const uint32_t N = 0;
  const float delta = 0.f;
  float score = 0.f;
  uint64_t index = 0;
  uint8_t rotation = 0;
  friend class FFTProcess;
  
public:
  Score();
  Score(uint32_t N, float delta);
  Score(const Score& other) = default;
  Score(Score&& other) = default;
  
  Score& operator=(const Score& other);
  Score& operator=(Score&& other);
  
  inline float GetScore() const { return score; }
  inline uint64_t GetRawIndex() const { return index; }
  inline int32_t GetXIndex() const {
    int32_t x = (int32_t)(index / (N * N));
    if (x > N/2) { x -= N; }
    return x;
  }
  inline float GetDeltaX() const { return GetXIndex() * delta; }
  inline int32_t GetYIndex() const {
    int32_t y = (index / N) % N;
    if (y > N/2) { y -= N; }
    return y;
  }
  inline float GetDeltaY() const { return GetYIndex() * delta; }
  inline int32_t GetZIndex() const {
    int32_t z = index % N;
    if (z > N/2) { z -= N; }
    return z;
  }
  inline float GetDeltaZ() const { return GetZIndex() * delta; }
  inline uint32_t GetRotationType() const { return (uint32_t)rotation; }
  
  inline bool operator<(const Score& b) const { return score < b.score; }
  inline bool operator>(const Score& b) const { return score > b.score; }
  inline bool operator<=(const Score& b) const { return score <= b.score; }
  inline bool operator>=(const Score& b) const { return score >= b.score; }
  
  
  inline void swap(Score& other) {
    float tmps = score; uint64_t tmpi = index; uint8_t tmpr = rotation;
    score = other.score;
    index = other.index;
    rotation = other.rotation;
    other.score = tmps;
    other.index = tmpi;
    other.rotation = tmpr;
  }
};

inline void swap(Score& a, Score& b) { a.swap(b); }

}
