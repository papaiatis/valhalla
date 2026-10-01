#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace valhalla {
namespace baldr {

/**
 * Where the parts of a tile that a package join changes lie, in 8-byte words from the start of the
 * tile. A region that is empty or not 8-byte aligned has count 0.
 */
struct TileWords {
  size_t transitions_first = 0, transitions_count = 0; // node transitions
  size_t edges_first = 0, edges_count = 0;             // directed edges, 6 words each
  size_t bins_first = 0, bins_count = 0;               // edge bin entries (graph ids)
};

/**
 * A tile patch, and a hash of its stream before deflating. Deflate output differs between zlib
 * builds, the stream doesn't: equal hashes mean equal changes on any platform.
 */
struct TilePatch {
  std::string bytes; // empty when the tiles are equal
  uint64_t content_hash = 0;
};

/**
 * The patch that turns a package tile into its joined tile. Both are `size` bytes, and the join
 * never changes bytes after the last whole 8-byte word (throws std::logic_error if it did).
 * `layout` describes the package tile.
 *
 * A join changes whole 8-byte words, mostly in a few regular ways, so each changed word goes into
 * the class that describes its change cheapest, and each class lists the positions of its words as
 * runs (the gap since the previous run's end, then the length) in a space of its own:
 *  - other: any word; the change is its XOR with the original; positions are word indexes
 *  - node transitions: XOR; positions are transition indexes
 *  - first word of a directed edge (end node and flags): XOR; positions are edge indexes
 *  - not_thru of a directed edge cleared: no value; positions are ranks among the edges that have
 *    not_thru set
 *  - bin entries set to the invalid graph id: no value; positions are bin entry indexes
 *  - other bin entry changes: XOR; positions are bin entry indexes
 * The XORs of a class follow, one byte plane at a time (the same byte of every word together).
 * Everything is deflated; the patch is the tile size and the deflated size (u32 each), then the
 * deflated stream.
 */
TilePatch
MakeTilePatch(const char* base, const char* joined, size_t size, const TileWords& layout, int level);

/**
 * Applies a patch to the package tile it was made for: `size` bytes at `tile`, described by
 * `layout`. The whole patch is checked before the first byte changes; applying a patch twice would
 * undo it, so a patch that fails must leave the tile alone. Throws std::runtime_error when the
 * patch is damaged or doesn't belong to the tile.
 */
void ApplyTilePatch(const char* patch,
                    size_t patch_size,
                    char* tile,
                    size_t size,
                    const TileWords& layout);

} // namespace baldr
} // namespace valhalla
