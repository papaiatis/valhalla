#include "baldr/tilepatch.h"
#include "baldr/directededge.h"
#include "baldr/graphid.h"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using valhalla::baldr::TileWords;

constexpr size_t kEdgeWords = sizeof(valhalla::baldr::DirectedEdge) / 8;
static_assert(sizeof(valhalla::baldr::DirectedEdge) % 8 == 0, "A directed edge is whole words");
// not_thru is the top bit of a directed edge's second word
constexpr size_t kFlagsWord = 1;
constexpr uint64_t kNotThru = uint64_t{1} << 63;
constexpr uint64_t kInvalid = valhalla::baldr::kInvalidGraphId;

enum Class : size_t {
  kOther,
  kTransitions,
  kEdgeStart,
  kNotThruCleared,
  kBinInvalid,
  kBinXor,
  kClasses
};
// whether the words of a class carry an XOR value
constexpr std::array<bool, kClasses> kHasValues = {true, true, true, false, false, true};

uint64_t Word(const char* tile, size_t index) {
  uint64_t word;
  std::memcpy(&word, tile + 8 * index, 8);
  return word;
}

void SetWord(char* tile, size_t index, uint64_t word) {
  std::memcpy(tile + 8 * index, &word, 8);
}

uint64_t GetLE(const char* data, int bytes) {
  uint64_t value = 0;
  for (int i = bytes - 1; i >= 0; --i) {
    value = (value << 8) | static_cast<unsigned char>(data[i]);
  }
  return value;
}

void PutU32(std::string& out, uint32_t value) {
  for (int i = 0; i < 4; ++i) {
    out.push_back(static_cast<char>(value >> (8 * i)));
  }
}

void PutVarint(std::string& out, uint64_t value) {
  while (value >= 0x80) {
    out.push_back(static_cast<char>(value | 0x80));
    value >>= 7;
  }
  out.push_back(static_cast<char>(value));
}

bool GetVarint(const unsigned char*& at, const unsigned char* end, uint64_t& value) {
  value = 0;
  for (int shift = 0; shift < 64 && at < end; shift += 7) {
    const unsigned char byte = *at++;
    value |= static_cast<uint64_t>(byte & 0x7f) << shift;
    if (!(byte & 0x80)) {
      return true;
    }
  }
  return false;
}

// Runs of consecutive positions: how many, then every gap since the end of the previous run, then
// every length (the gaps together deflate better than gaps mixed with lengths).
void PutRuns(std::string& out, const std::vector<size_t>& positions) {
  std::vector<std::pair<size_t, size_t>> runs; // (first, length)
  for (const size_t position : positions) {
    if (!runs.empty() && runs.back().first + runs.back().second == position) {
      ++runs.back().second;
    } else {
      runs.emplace_back(position, 1);
    }
  }
  PutVarint(out, runs.size());
  size_t previous_end = 0;
  for (const auto& [first, length] : runs) {
    PutVarint(out, first - previous_end);
    previous_end = first + length;
  }
  for (const auto& run : runs) {
    PutVarint(out, run.second);
  }
}

struct Runs {
  std::vector<std::pair<size_t, size_t>> runs; // (first position, length)
  size_t words = 0;                            // positions in all runs
};

// Reads runs of positions below `limit`; false if the stream is malformed.
bool GetRuns(const unsigned char*& at, const unsigned char* end, size_t limit, Runs& out) {
  uint64_t count;
  if (!GetVarint(at, end, count) || count > static_cast<uint64_t>(end - at)) {
    return false;
  }
  std::vector<size_t> gaps(count);
  for (auto& gap : gaps) {
    uint64_t value;
    if (!GetVarint(at, end, value)) {
      return false;
    }
    gap = value;
  }
  size_t next = 0;
  out.runs.reserve(count);
  for (const auto gap : gaps) {
    uint64_t length;
    if (!GetVarint(at, end, length) || gap > limit - next || length == 0 ||
        length > limit - next - gap) {
      return false;
    }
    out.runs.emplace_back(next + gap, length);
    next += gap + length;
    out.words += length;
  }
  return true;
}

// Positions of the edges whose not_thru is set.
std::vector<size_t> NotThruEdges(const char* tile, const TileWords& layout) {
  std::vector<size_t> edges;
  for (size_t e = 0; e < layout.edges_count; ++e) {
    if (Word(tile, layout.edges_first + kEdgeWords * e + kFlagsWord) & kNotThru) {
      edges.push_back(e);
    }
  }
  return edges;
}

// Whether the region of `count` words from `first` lies inside a tile of `words` words.
bool Inside(size_t first, size_t count, size_t words, size_t stride = 1) {
  return count == 0 || (first <= words && count <= (words - first) / stride);
}

} // namespace

namespace valhalla {
namespace baldr {

TilePatch
MakeTilePatch(const char* base, const char* joined, size_t size, const TileWords& layout, int level) {
  const size_t words = size / 8;
  if (std::memcmp(base + 8 * words, joined + 8 * words, size - 8 * words) != 0) {
    throw std::logic_error("The join changed the bytes after the last 8-byte word of a tile");
  }
  if (!Inside(layout.transitions_first, layout.transitions_count, words) ||
      !Inside(layout.edges_first, layout.edges_count, words, kEdgeWords) ||
      !Inside(layout.bins_first, layout.bins_count, words)) {
    throw std::logic_error("The layout of a tile reaches beyond the tile");
  }

  // the rank of every edge that has not_thru set, for the class of edges that clear it
  std::vector<size_t> not_thru_rank(layout.edges_count, 0);
  {
    size_t rank = 0;
    for (size_t e = 0; e < layout.edges_count; ++e) {
      not_thru_rank[e] = rank;
      rank += (Word(base, layout.edges_first + kEdgeWords * e + kFlagsWord) & kNotThru) ? 1 : 0;
    }
  }

  std::array<std::vector<size_t>, kClasses> positions;
  std::array<std::vector<uint64_t>, kClasses> values;
  const auto add = [&](Class c, size_t position, uint64_t xor_value) {
    positions[c].push_back(position);
    if (kHasValues[c]) {
      values[c].push_back(xor_value);
    }
  };
  const auto inside = [](size_t i, size_t first, size_t count, size_t stride) {
    return i >= first && i - first < count * stride;
  };
  bool changed = false;
  for (size_t i = 0; i < words; ++i) {
    const uint64_t a = Word(base, i), b = Word(joined, i);
    if (a == b) {
      continue;
    }
    changed = true;
    if (inside(i, layout.transitions_first, layout.transitions_count, 1)) {
      add(kTransitions, i - layout.transitions_first, a ^ b);
    } else if (inside(i, layout.edges_first, layout.edges_count, kEdgeWords)) {
      const size_t edge = (i - layout.edges_first) / kEdgeWords;
      const size_t word = (i - layout.edges_first) % kEdgeWords;
      if (word == 0) {
        add(kEdgeStart, edge, a ^ b);
      } else if (word == kFlagsWord && (a ^ b) == kNotThru && (a & kNotThru)) {
        add(kNotThruCleared, not_thru_rank[edge], 0);
      } else {
        add(kOther, i, a ^ b);
      }
    } else if (inside(i, layout.bins_first, layout.bins_count, 1)) {
      if (b == kInvalid) {
        add(kBinInvalid, i - layout.bins_first, 0);
      } else {
        add(kBinXor, i - layout.bins_first, a ^ b);
      }
    } else {
      add(kOther, i, a ^ b);
    }
  }
  if (!changed) {
    return {};
  }
  TilePatch result;

  std::string stream;
  for (size_t c = 0; c < kClasses; ++c) {
    PutRuns(stream, positions[c]);
  }
  for (size_t c = 0; c < kClasses; ++c) {
    for (int plane = 0; plane < 8 && kHasValues[c]; ++plane) {
      for (const uint64_t value : values[c]) {
        stream.push_back(static_cast<char>(value >> (8 * plane)));
      }
    }
  }

  result.content_hash = 14695981039346656037ull; // FNV-1a
  for (const unsigned char byte : stream) {
    result.content_hash = (result.content_hash ^ byte) * 1099511628211ull;
  }

  uLongf deflated = compressBound(static_cast<uLong>(stream.size()));
  std::string& patch = result.bytes;
  PutU32(patch, static_cast<uint32_t>(size));
  PutU32(patch, static_cast<uint32_t>(stream.size()));
  patch.resize(patch.size() + deflated);
  if (compress2(reinterpret_cast<Bytef*>(&patch[8]), &deflated,
                reinterpret_cast<const Bytef*>(stream.data()), static_cast<uLong>(stream.size()),
                level) != Z_OK) {
    throw std::runtime_error("Could not compress a tile patch");
  }
  patch.resize(8 + deflated);
  return result;
}

void ApplyTilePatch(const char* patch,
                    size_t patch_size,
                    char* tile,
                    size_t size,
                    const TileWords& layout) {
  const auto damaged = [] { return std::runtime_error("Tile patch is damaged"); };
  if (patch_size < 8 || GetLE(patch, 4) != size) {
    throw std::runtime_error("Tile patch does not belong to a tile of " + std::to_string(size) +
                             " bytes");
  }
  const size_t words = size / 8;
  if (!Inside(layout.transitions_first, layout.transitions_count, words) ||
      !Inside(layout.edges_first, layout.edges_count, words, kEdgeWords) ||
      !Inside(layout.bins_first, layout.bins_count, words)) {
    throw std::runtime_error("The layout of a tile reaches beyond the tile");
  }

  const size_t stream_size = GetLE(patch + 4, 4);
  std::vector<unsigned char> stream(stream_size);
  uLongf inflated = static_cast<uLongf>(stream_size);
  if (uncompress(stream.data(), &inflated, reinterpret_cast<const Bytef*>(patch + 8),
                 static_cast<uLong>(patch_size - 8)) != Z_OK ||
      inflated != stream_size) {
    throw damaged();
  }

  // everything is read and checked before the first word changes
  const unsigned char* at = stream.data();
  const unsigned char* const end = at + stream_size;
  std::array<Runs, kClasses> runs;
  const std::array<size_t, kClasses> limits = {words,
                                               layout.transitions_count,
                                               layout.edges_count,
                                               layout.edges_count,
                                               layout.bins_count,
                                               layout.bins_count};
  for (size_t c = 0; c < kClasses; ++c) {
    if (!GetRuns(at, end, limits[c], runs[c])) {
      throw damaged();
    }
  }
  std::vector<size_t> not_thru_edges;
  if (runs[kNotThruCleared].words > 0) {
    not_thru_edges = NotThruEdges(tile, layout);
    if (runs[kNotThruCleared].runs.back().first + runs[kNotThruCleared].runs.back().second >
        not_thru_edges.size()) {
      throw damaged();
    }
  }
  size_t value_words = 0;
  std::array<size_t, kClasses> value_start = {};
  for (size_t c = 0; c < kClasses; ++c) {
    value_start[c] = value_words;
    value_words += kHasValues[c] ? runs[c].words : 0;
  }
  if (static_cast<size_t>(end - at) != 8 * value_words) {
    throw damaged();
  }

  // the n-th changed word of a class, with the XOR value of its class
  const auto value = [&](Class c, size_t n) {
    uint64_t result = 0;
    for (int plane = 0; plane < 8; ++plane) {
      result |= static_cast<uint64_t>(at[8 * value_start[c] + plane * runs[c].words + n])
                << (8 * plane);
    }
    return result;
  };
  for (size_t c = 0; c < kClasses; ++c) {
    size_t n = 0;
    for (const auto& [first, length] : runs[c].runs) {
      for (size_t position = first; position < first + length; ++position, ++n) {
        switch (static_cast<Class>(c)) {
          case kOther:
            SetWord(tile, position, Word(tile, position) ^ value(kOther, n));
            break;
          case kTransitions: {
            const size_t i = layout.transitions_first + position;
            SetWord(tile, i, Word(tile, i) ^ value(kTransitions, n));
            break;
          }
          case kEdgeStart: {
            const size_t i = layout.edges_first + kEdgeWords * position;
            SetWord(tile, i, Word(tile, i) ^ value(kEdgeStart, n));
            break;
          }
          case kNotThruCleared: {
            const size_t i = layout.edges_first + kEdgeWords * not_thru_edges[position] + kFlagsWord;
            SetWord(tile, i, Word(tile, i) ^ kNotThru);
            break;
          }
          case kBinInvalid:
            SetWord(tile, layout.bins_first + position, kInvalid);
            break;
          case kBinXor: {
            const size_t i = layout.bins_first + position;
            SetWord(tile, i, Word(tile, i) ^ value(kBinXor, n));
            break;
          }
          default:
            break;
        }
      }
    }
  }
}

} // namespace baldr
} // namespace valhalla
