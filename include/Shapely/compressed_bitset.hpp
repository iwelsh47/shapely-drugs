#pragma once

#include <Shapely/serialisation.hpp>
#include <sul/dynamic_bitset.hpp>
#include <vector>

namespace shapely {

class CompressedBitset {
  using sz_t = sul::dynamic_bitset<>::size_type;
  using CompressedUnit = std::pair<sz_t, uint32_t>;
  
public:
  CompressedBitset() = default;
  CompressedBitset(const sul::dynamic_bitset<>& bitset)  {
    SetBits(bitset);
  }
  
  CompressedBitset(const CompressedBitset&) = default;
  CompressedBitset(CompressedBitset&&) = default;
  
  CompressedBitset& operator=(const CompressedBitset&) = default;
  CompressedBitset& operator=(CompressedBitset&&) = default;
  CompressedBitset& operator=(const sul::dynamic_bitset<>& bitset) {
    bits.clear();
    SetBits(bitset);
    return *this;
  }
  
  void clear() {
    num_bits = 0;
    bits.clear();
  }
  
  CompressedBitset operator|(const CompressedBitset& other) const {
    if (num_bits != other.num_bits) {
      throw std::range_error("Bitsets must have same number of bits.");
    }
    
    CompressedBitset result;
    result.num_bits = num_bits;
    result.bits.reserve(bits.size() + other.bits.size());
    
    
    
    return result;
  }
  
  uint64_t count() const {
    uint64_t c = 0;
    for (const CompressedUnit& unit : bits) { c += unit.second; }
    return c;
  }
  
  template <class Archive>
  void CEREAL_SERIALIZE_FUNCTION_NAME(Archive& ar, const uint32_t version) {
    ar(num_bits, bits);
  }
  
  explicit operator sul::dynamic_bitset<>() const {
    sul::dynamic_bitset<> result(num_bits);
    for (const CompressedUnit& sc : bits) {
      result.set(sc.first, sc.second, true);
    }
    return result;
  }
  
private:
  void SetBits(const sul::dynamic_bitset<>& bitset) {
    num_bits = bitset.size();
    bits.reserve(bitset.num_blocks());
    sul::dynamic_bitset<> rbitset = ~bitset;
    
    sz_t start = bitset.find_first();
    sz_t end = rbitset.find_next(start);
    
    while (start != bitset.npos) {
      sz_t count = end - start;
      while (count > (sz_t)std::numeric_limits<uint32_t>::max()) {
        bits.emplace_back(start, (uint32_t)std::numeric_limits<uint32_t>::max());
        count -= (sz_t)std::numeric_limits<uint32_t>::max();
        start += (sz_t)std::numeric_limits<uint32_t>::max();
      }
      if (count) {
        bits.emplace_back(start, (uint32_t)count);
      }
      start = bitset.find_next(end - 1);
      end = rbitset.find_next(start);
    }
    bits.shrink_to_fit();
  }
  
private:
  uint64_t num_bits;
  std::vector<CompressedUnit> bits;
};

}
