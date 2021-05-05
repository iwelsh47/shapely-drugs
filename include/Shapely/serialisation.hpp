#pragma once

#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

#include <sul/dynamic_bitset.hpp>
#include <stdx/fixed_vector.hpp>
#include <Eigen/Dense>

namespace cereal {

// Serialising for sul::dynamic_bitset when BinaryData optimisation supported
template <class Archive, class Block, class Allocator > inline
std::enable_if_t<traits::is_output_serializable<BinaryData<Block>, Archive>::value, void>
CEREAL_SAVE_FUNCTION_NAME( Archive& ar, sul::dynamic_bitset<Block, Allocator> const & bits) {
  // Save in either raw dump or compressed format, depending on size.
  // Assume compressed first. Keep building compressed until larger than uncompressed.
  using sz_t = typename sul::dynamic_bitset<Block, Allocator>::size_type;
  using CompressedUnit = std::pair<sz_t, uint32_t>;
  
  sz_t max_compressed = sizeof(Block) * bits.num_blocks() / sizeof(CompressedUnit);
  
  std::vector<CompressedUnit> counts;
  counts.reserve(bits.num_blocks());
  
  sul::dynamic_bitset<Block, Allocator> rbits = ~bits;
  
  sz_t start = bits.find_first();
  sz_t end = rbits.find_next(start);
  while (start != bits.npos && counts.size() < max_compressed) {
    sz_t count = end - start;
    while (count > (sz_t)std::numeric_limits<uint32_t>::max()) {
      counts.emplace_back(start, (uint32_t)std::numeric_limits<uint32_t>::max());
      count -= (sz_t)std::numeric_limits<uint32_t>::max();
      start += (sz_t)std::numeric_limits<uint32_t>::max();
    }
    if (count) {
      counts.emplace_back(start, (uint32_t)count);
    }
    start = bits.find_next(end - 1);
    end = rbits.find_next(start);
  }
  
  ar( make_size_tag(static_cast<size_type>(bits.size()))); // Number of bits
  
  if (counts.size() < max_compressed) {
    ar((uint8_t)1, counts);
  } else {
    ar((uint8_t)0);
    ar( make_size_tag(static_cast<size_type>(bits.num_blocks()))); // Number of blocks
    ar( binary_data(bits.data(), bits.num_blocks() * sizeof(Block)) );
  }
}

// Deserialising for sul::dynamic_bitset when BinaryData optimisation supported
template <class Archive, class Block, class Allocator > inline
std::enable_if_t<traits::is_input_serializable<BinaryData<Block>, Archive>::value, void>
CEREAL_LOAD_FUNCTION_NAME( Archive& ar, sul::dynamic_bitset<Block, Allocator>& bits) {
  uint8_t is_compressed;
  size_type nbits;
  ar(make_size_tag(nbits), is_compressed);
  bits.clear();
  bits.resize(nbits);
  
  if (is_compressed) {
    using sz_t = typename sul::dynamic_bitset<Block, Allocator>::size_type;
    using CompressedUnit = std::pair<sz_t, uint32_t>;
    
    std::vector<CompressedUnit> counts;
    ar(counts);
    for (CompressedUnit& sc : counts) { bits.set(sc.first, sc.second, true); }
    
  } else {
    size_type nblocks;
    ar( make_size_tag(nblocks) );
    ar(binary_data(bits.data(), nblocks * sizeof(Block)));
  }
}

// Serialising for sul::dynamic_bitset when BinaryData optimisation is not supported
template <class Archive, class Block, class Allocator > inline
std::enable_if_t<!traits::is_output_serializable<BinaryData<Block>, Archive>::value, void>
CEREAL_SAVE_FUNCTION_NAME( Archive& ar, sul::dynamic_bitset<Block, Allocator> const & bits) {
  ar(make_size_tag(bits.num_blocks()), make_size_tag(bits.size()));
  std::vector<Block> data(bits.data(), bits.data() + bits.num_blocks());
  ar(CEREAL_NVP_("data", data));
}

// Deserialising for sul::dynamic_bitset when BinaryData optimisation is not supported
template <class Archive, class Block, class Allocator > inline
std::enable_if_t<!traits::is_input_serializable<BinaryData<Block>, Archive>::value, void>
CEREAL_LOAD_FUNCTION_NAME( Archive& ar, sul::dynamic_bitset<Block, Allocator>& bits) {
  std::string sbits;
  ar(CEREAL_NVP_("data", sbits));
  bits = sul::dynamic_bitset<Block, Allocator>(sbits);
}

// Serialising fixed_vector
template<class Archive, class T, class Allocator> inline
void CEREAL_SAVE_FUNCTION_NAME( Archive& ar, stdx::fixed_vector<T, Allocator> const& vector) {
  ar(CEREAL_NVP_("min_index", vector.min_index()),
     CEREAL_NVP_("max_index", vector.max_index()),
     CEREAL_NVP_("data", vector.get_elems()));
}

// Deserialising fixed_vector
template<class Archive, class T, class Allocator> inline
void CEREAL_LOAD_FUNCTION_NAME( Archive& ar, stdx::fixed_vector<T, Allocator>& vector) {
  int64_t min, max;
  ar(CEREAL_NVP_("min_index", min),
     CEREAL_NVP_("max_index", max),
     CEREAL_NVP_("data", vector.get_elems()));
  vector.resize(min, max);
}

// Serialising Eigen matrix()
template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
void
CEREAL_SAVE_FUNCTION_NAME(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> const & m)
    {
      uint64_t rows = m.rows();
      uint64_t cols = m.cols();
      ar(rows, cols);

      for (int i = 0; i < rows; i++)
          for (int j = 0; j < cols; j++)
              ar(m(i,j));
    }

  template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
void CEREAL_LOAD_FUNCTION_NAME(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & m)
    {
      uint64_t rows;
      uint64_t cols;
      ar(rows, cols);
      m.resize(rows, cols);

      for (int i = 0; i < rows; i++)
          for (int j = 0; j < cols; j++)
              ar(m(i,j));
    }

}
