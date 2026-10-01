#include "baldr/tilepatch.h"
#include "baldr/directededge.h"
#include "baldr/graphid.h"
#include "baldr/nodetransition.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using namespace valhalla::baldr;

namespace {

constexpr size_t kWords = 4000;
constexpr int kLevel = 9;

// A tile-shaped buffer: transitions, directed edges (6 words each), and bin entries, with other
// words around them; the layout says where.
struct Fixture {
  std::vector<char> base;
  TileWords layout;

  explicit Fixture(uint32_t seed) : base(8 * kWords) {
    layout.transitions_first = 50;
    layout.transitions_count = 300;
    layout.edges_first = 400;
    layout.edges_count = 400; // 2400 words, to 2800
    layout.bins_first = 2900;
    layout.bins_count = 600;
    std::mt19937_64 random(seed);
    for (size_t i = 0; i < kWords; ++i) {
      uint64_t word = random();
      // not_thru is set on half of the edges, and a quarter of the bin entries are invalid
      if (i >= layout.edges_first && i < layout.edges_first + 6 * layout.edges_count &&
          (i - layout.edges_first) % 6 == 1) {
        word = (word & ~(uint64_t{1} << 63)) | (random() & 1) << 63;
      }
      if (i >= layout.bins_first && i < layout.bins_first + layout.bins_count) {
        word = random() % 4 == 0 ? kInvalidGraphId : word & kInvalidGraphId;
      }
      std::memcpy(&base[8 * i], &word, 8);
    }
  }

  uint64_t Word(const std::vector<char>& tile, size_t i) const {
    uint64_t word;
    std::memcpy(&word, &tile[8 * i], 8);
    return word;
  }
  void Set(std::vector<char>& tile, size_t i, uint64_t word) const {
    std::memcpy(&tile[8 * i], &word, 8);
  }
  size_t EdgeWord(size_t edge, size_t word) const {
    return layout.edges_first + 6 * edge + word;
  }
};

std::string Patch(const Fixture& f, const std::vector<char>& joined, int level = kLevel) {
  return MakeTilePatch(f.base.data(), joined.data(), f.base.size(), f.layout, level).bytes;
}

std::vector<char> Applied(const Fixture& f, const std::string& patch) {
  auto tile = f.base;
  ApplyTilePatch(patch.data(), patch.size(), tile.data(), tile.size(), f.layout);
  return tile;
}

} // namespace

// The classes of a patch assume these layouts of the tile structures.
TEST(TilePatch, TileLayoutAssumptions) {
  EXPECT_EQ(sizeof(DirectedEdge), 48u);
  EXPECT_EQ(sizeof(NodeTransition), 8u);
  EXPECT_EQ(sizeof(GraphId), 8u);
  EXPECT_EQ(GraphId().value, kInvalidGraphId);
  EXPECT_EQ(kInvalidGraphId, 0x3fffffffffffull);

  // not_thru is the top bit of the edge's second word, and no other bit of it
  DirectedEdge edge;
  uint64_t before[6], after[6];
  std::memcpy(before, &edge, sizeof(edge));
  edge.set_not_thru(true);
  std::memcpy(after, &edge, sizeof(edge));
  for (size_t w = 0; w < 6; ++w) {
    EXPECT_EQ(before[w] ^ after[w], w == 1 ? uint64_t{1} << 63 : 0u) << "word " << w;
  }
}

TEST(TilePatch, EqualTilesNeedNoPatch) {
  const Fixture f(1);
  EXPECT_TRUE(Patch(f, f.base).empty());
}

TEST(TilePatch, EveryKindOfChangeRoundTrips) {
  const Fixture f(2);
  auto joined = f.base;
  // other words (before the transitions, between regions, after the bins)
  for (const size_t i :
       {size_t{0}, size_t{1}, size_t{49}, size_t{360}, size_t{2850}, size_t{3600}, kWords - 1}) {
    f.Set(joined, i, ~f.Word(joined, i));
  }
  // transitions
  for (size_t t = 0; t < 300; t += 3) {
    f.Set(joined, f.layout.transitions_first + t,
          f.Word(joined, f.layout.transitions_first + t) + 77 * t);
  }
  size_t cleared = 0, other_flags = 0;
  for (size_t e = 0; e < f.layout.edges_count; ++e) {
    // end nodes
    if (e % 4 == 0) {
      f.Set(joined, f.EdgeWord(e, 0), f.Word(joined, f.EdgeWord(e, 0)) ^ (0xbf9000 + e));
    }
    const uint64_t flags = f.Word(joined, f.EdgeWord(e, 1));
    if (flags >> 63) {
      if (e % 3 != 0) {
        // not_thru cleared alone
        f.Set(joined, f.EdgeWord(e, 1), flags & ~(uint64_t{1} << 63));
        ++cleared;
      } else {
        // not_thru cleared along with another change: not that class
        f.Set(joined, f.EdgeWord(e, 1), (flags & ~(uint64_t{1} << 63)) ^ 5);
        ++other_flags;
      }
    } else if (e % 5 == 0) {
      // not_thru set where it was clear: not that class either
      f.Set(joined, f.EdgeWord(e, 1), flags | (uint64_t{1} << 63));
      ++other_flags;
    }
    // other words of an edge
    if (e % 7 == 0) {
      f.Set(joined, f.EdgeWord(e, 3), 12345 + e);
    }
  }
  EXPECT_GT(cleared, 20u);
  EXPECT_GT(other_flags, 20u);
  // bin entries: set invalid, remapped, and made valid again
  for (size_t b = 0; b < f.layout.bins_count; ++b) {
    const size_t i = f.layout.bins_first + b;
    if (b % 5 == 0) {
      f.Set(joined, i, kInvalidGraphId);
    } else if (b % 5 == 1) {
      f.Set(joined, i, f.Word(joined, i) ^ 0xff000);
    } else if (b % 5 == 2 && f.Word(joined, i) == kInvalidGraphId) {
      f.Set(joined, i, 424242);
    }
  }

  const auto patch = Patch(f, joined);
  ASSERT_FALSE(patch.empty());
  EXPECT_TRUE(Applied(f, patch) == joined);
  EXPECT_LT(patch.size(), joined.size());
}

// The hash covers what changes, not how the deflate of this zlib build wrote it.
TEST(TilePatch, ContentHashDoesNotDependOnTheDeflateLevel) {
  const Fixture f(6);
  std::mt19937_64 random(6);
  auto joined = f.base;
  for (size_t k = 0; k < 800; ++k) {
    f.Set(joined, random() % kWords, random());
  }
  const auto make = [&](int level) {
    return MakeTilePatch(f.base.data(), joined.data(), f.base.size(), f.layout, level);
  };
  const auto fast = make(1), best = make(9);
  EXPECT_NE(fast.bytes, best.bytes);
  EXPECT_EQ(fast.content_hash, best.content_hash);
  EXPECT_NE(fast.content_hash, 0u);
  // another change, another hash
  f.Set(joined, 7, f.Word(joined, 7) ^ 1);
  EXPECT_NE(make(9).content_hash, best.content_hash);
}

TEST(TilePatch, RandomChangesRoundTrip) {
  for (uint32_t seed = 10; seed < 30; ++seed) {
    const Fixture f(seed);
    std::mt19937_64 random(seed);
    auto joined = f.base;
    for (size_t k = 0; k < 1 + seed * 40; ++k) {
      f.Set(joined, random() % kWords, random());
    }
    const auto patch = Patch(f, joined);
    ASSERT_FALSE(patch.empty());
    EXPECT_TRUE(Applied(f, patch) == joined) << "seed " << seed;
  }
}

TEST(TilePatch, RegionsLeftOutGoIntoTheGeneralClass) {
  Fixture f(3);
  auto joined = f.base;
  for (size_t i = 0; i < kWords; i += 3) {
    f.Set(joined, i, ~f.Word(joined, i));
  }
  f.layout = {}; // a tile whose regions aren't 8-byte aligned
  const auto patch = Patch(f, joined);
  EXPECT_TRUE(Applied(f, patch) == joined);
}

TEST(TilePatch, TheLastBytesOfATileThatIsNotWholeWordsStay) {
  Fixture f(4);
  f.base.resize(8 * kWords + 5);
  auto joined = f.base;
  f.Set(joined, 100, 1);
  const auto patch = Patch(f, joined);
  EXPECT_TRUE(Applied(f, patch) == joined);
  joined.back() ^= 1;
  EXPECT_THROW(Patch(f, joined), std::logic_error);
}

TEST(TilePatch, DamagedPatchesLeaveTheTileAlone) {
  const Fixture f(5);
  std::mt19937_64 random(5);
  auto joined = f.base;
  for (size_t k = 0; k < 600; ++k) {
    f.Set(joined, random() % kWords, random());
  }
  const auto patch = Patch(f, joined);

  // wrong tile size, truncated, and with its sizes changed
  auto tile = f.base;
  EXPECT_THROW(ApplyTilePatch(patch.data(), patch.size(), tile.data(), tile.size() - 8, f.layout),
               std::runtime_error);
  EXPECT_THROW(ApplyTilePatch(patch.data(), 7, tile.data(), tile.size(), f.layout),
               std::runtime_error);
  for (const size_t keep : {size_t{8}, size_t{9}, patch.size() / 2, patch.size() - 1}) {
    EXPECT_THROW(ApplyTilePatch(patch.data(), keep, tile.data(), tile.size(), f.layout),
                 std::runtime_error)
        << keep;
    EXPECT_TRUE(tile == f.base) << keep;
  }
  // a layout reaching beyond the tile
  auto outside = f.layout;
  outside.edges_count = 100000;
  EXPECT_THROW(ApplyTilePatch(patch.data(), patch.size(), tile.data(), tile.size(), outside),
               std::runtime_error);
  EXPECT_TRUE(tile == f.base);

  // random damage: a patch either fails and leaves the tile as it was, or applies; it never reads
  // or writes outside the tile (checked in a tile followed by guard bytes)
  for (size_t round = 0; round < 3000; ++round) {
    auto damaged = patch;
    const size_t flips = 1 + random() % 3;
    for (size_t k = 0; k < flips; ++k) {
      damaged[random() % damaged.size()] ^= static_cast<char>(1 << (random() % 8));
    }
    std::vector<char> guarded = f.base;
    guarded.resize(f.base.size() + 64, 0x5a);
    try {
      ApplyTilePatch(damaged.data(), damaged.size(), guarded.data(), f.base.size(), f.layout);
    } catch (const std::runtime_error&) {
      EXPECT_TRUE(std::vector<char>(guarded.begin(), guarded.begin() + f.base.size()) == f.base)
          << "round " << round;
    }
    EXPECT_TRUE(std::vector<char>(guarded.begin() + f.base.size(), guarded.end()) ==
                std::vector<char>(64, 0x5a))
        << "round " << round;
  }
}
