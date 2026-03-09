#ifndef _CONTIG_TABLE_HPP_
#define _CONTIG_TABLE_HPP_

#include <vector>
#include <cstdint>
#include "Util.hpp"
#include "core/range.hpp"
#include "elias_fano.hpp"

namespace pufferfish {

class ContigTable {
public:
    using position_type = pufferfish::util::Position;
    using range_type = core::range<const position_type*>;

    ContigTable() = default;

    // Build from construction data: positions is the flat posting list,
    // offsets[i] is the start index of contig i's posting list (size = numContigs + 1)
    void build(std::vector<position_type>&& positions,
               const std::vector<uint64_t>& offsets) {
        positions_ = std::move(positions);
        num_contigs_ = offsets.size() > 0 ? offsets.size() - 1 : 0;
        offsets_.encode(offsets.begin(), offsets.size(),
                        offsets.empty() ? 0 : offsets.back());
    }

    // Get the posting list for contig `rank`
    range_type contigRange(uint64_t contigRank) const {
        uint64_t start = offsets_.access(contigRank);
        uint64_t end = offsets_.access(contigRank + 1);
        const position_type* base = positions_.data();
        return range_type(base + start, base + end);
    }

    uint64_t numContigs() const { return num_contigs_; }
    uint64_t numPositions() const { return positions_.size(); }

    // Access the raw positions vector (for legacy compatibility)
    const std::vector<position_type>& positions() const { return positions_; }
    std::vector<position_type>& positions() { return positions_; }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

private:
    std::vector<position_type> positions_;
    bits::elias_fano<> offsets_;
    uint64_t num_contigs_{0};

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.num_contigs_);
        visitor.visit(t.positions_);
        visitor.visit(t.offsets_);
    }
};

}  // namespace pufferfish

#endif  // _CONTIG_TABLE_HPP_
