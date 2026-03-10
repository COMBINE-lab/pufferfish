#include <fstream>
#include <iostream>
#include <future>
#include <stdlib.h>

#include "CLI/Timer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "PufferfishIndex.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"

#include "jellyfish/mer_dna.hpp"

// SSHash essentials for load/save
#include "essentials.hpp"

PufferfishIndex::PufferfishIndex() { }

PufferfishIndex::PufferfishIndex(const std::string& indexDir, pufferfish::util::IndexLoadingOpts opts) {
  if (!puffer::fs::DirExists(indexDir.c_str())) {
    std::cerr << "The index directory " << indexDir << " does not exist!\n";
    std::exit(1);
  }

  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("k", k_));
    infoArchive(cereal::make_nvp("num_kmers", numKmers_));
    infoArchive(cereal::make_nvp("have_edge_vec", haveEdges_));
    infoArchive(cereal::make_nvp("have_ref_seq", haveRefSeq_));
    infoArchive(cereal::make_nvp("num_decoys", numDecoys_));
    infoArchive(cereal::make_nvp("first_decoy_index", firstDecoyIndex_));
    infoStream.close();
    twok_ = 2 * k_;
  }
  haveEdges_ = opts.try_loading_edges and haveEdges_;
  haveRefSeq_ = opts.try_loading_ref_seqs and haveRefSeq_;
  haveEqClasses_ = opts.try_loading_eqclasses and haveEqClasses_;

  // Detect index format: SSHash vs legacy BooPHF
  std::string sshashFile = indexDir + "/sshash.bin";
  useSSHash_ = puffer::fs::FileExists(sshashFile.c_str());

  {
    CLI::AutoTimer timer{"Loading contig table", CLI::Timer::Big};
    std::ifstream contigTableStream(indexDir + "/" + pufferfish::util::CTABLE);
    cereal::BinaryInputArchive contigTableArchive(contigTableStream);
    contigTableArchive(refNames_);
    contigTableArchive(refExt_);
    // Load raw positions and offsets, then build ContigTable
    std::vector<pufferfish::util::Position> rawContigTable;
    contigTableArchive(rawContigTable);
    contigTableStream.close();

    // Load contig offsets and build the EF-encoded ContigTable
    std::string pfile = indexDir + "/" + pufferfish::util::CONTIG_OFFSETS;
    compact::vector<uint64_t> rawOffsets{16};
    rawOffsets.deserialize(pfile, false);
    std::vector<uint64_t> offsets(rawOffsets.size());
    for (size_t i = 0; i < rawOffsets.size(); ++i) {
      offsets[i] = rawOffsets[i];
    }
    contigTable_.build(std::move(rawContigTable), offsets);
  }
  numContigs_ = contigTable_.numContigs();

  {
    std::string rlPath = indexDir + "/" + pufferfish::util::REFLENGTH;
    std::ifstream refLengthStream(rlPath);
    if(refLengthStream.good()) {
      CLI::AutoTimer timer{"Loading reference lengths", CLI::Timer::Big};
      cereal::BinaryInputArchive refLengthArchive(refLengthStream);
      refLengthArchive(refLengths_);
    } else {
      refLengths_ = std::vector<uint32_t>(refNames_.size(), 1000);
    }
  }

  {
    std::string rlPath = indexDir + "/" + pufferfish::util::COMPLETEREFLENGTH;
    std::ifstream completeRefLengthStream(rlPath);
    if(completeRefLengthStream.good()) {
      cereal::BinaryInputArchive completeRefLengthArchive(completeRefLengthStream);
      completeRefLengthArchive(completeRefLengths_);
    } else {
      throw std::runtime_error("could not load complete reference lengths!");
    }
  }

  if (haveEqClasses_) {
    CLI::AutoTimer timer{"Loading eq table", CLI::Timer::Big};
    std::ifstream eqTableStream(indexDir + "/" + pufferfish::util::EQTABLE);
    cereal::BinaryInputArchive eqTableArchive(eqTableStream);
    eqTableArchive(eqClassIDs_);
    eqTableArchive(eqLabels_);
    eqTableStream.close();
  }

  if (useSSHash_) {
    // --- SSHash index loading path ---
    {
      CLI::AutoTimer timer{"Loading SSHash dictionary", CLI::Timer::Big};
      essentials::load(dict_, sshashFile.c_str());
      // Verify k matches
      if (dict_.k() != k_) {
        std::cerr << "SSHash dictionary k=" << dict_.k()
                  << " does not match info.json k=" << k_ << "\n";
        std::exit(1);
      }
      numKmers_ = dict_.num_kmers();
    }
  } else {
    // --- Legacy BooPHF index loading path ---
    {
      CLI::AutoTimer timer{"Loading mphf table", CLI::Timer::Big};
      std::string hfile = indexDir + "/" + pufferfish::util::MPH;
      std::ifstream hstream(hfile);
      hash_.reset(new boophf_t);
      hash_->load(hstream);
      hstream.close();
      hash_raw_ = hash_.get();
    }

    {
      CLI::AutoTimer timer{"Loading contig boundaries", CLI::Timer::Big};
      std::string bfile = indexDir + "/" + pufferfish::util::RANK;
      contigBoundary_.deserialize(bfile, false);
      rankSelDict = rank9sel(&contigBoundary_, (uint64_t)contigBoundary_.size());
    }

    {
      CLI::AutoTimer timer{"Loading sequence", CLI::Timer::Big};
      std::string sfile = indexDir + "/" + pufferfish::util::SEQ;
      seq_.deserialize(sfile, false);
      lastSeqPos_ = seq_.size() - k_;
    }

    {
      CLI::AutoTimer timer{"Loading positions", CLI::Timer::Big};
      std::string pfile = indexDir + "/" + pufferfish::util::POS;
      pos_.deserialize(pfile, false);
    }
  }

  if (haveRefSeq_) {
    CLI::AutoTimer timer{"Loading reference sequence", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::REFSEQ;
    refseq_.deserialize(pfile, false);
  }

  {
    std::string rlPath = indexDir + "/" + pufferfish::util::REFACCUMLENGTH;
    std::ifstream refLengthStream(rlPath);
    if(refLengthStream.good()) {
      CLI::AutoTimer timer{"Loading reference accumulative lengths", CLI::Timer::Big};
      cereal::BinaryInputArchive refLengthArchive(refLengthStream);
      refLengthArchive(refAccumLengths_);
    } else {
      refAccumLengths_ = std::vector<uint64_t>(refNames_.size(), 1000);
    }
  }

  if (haveEdges_) {
    CLI::AutoTimer timer{"Loading edges", CLI::Timer::Big};
    std::string pfile = indexDir + "/" + pufferfish::util::EDGE;
    edge_.deserialize(pfile, false);
  }

  firstDecoyEncodedIndex_ = (numDecoys_ > 0) ? getRefId(firstDecoyIndex_) : std::numeric_limits<decltype(firstDecoyEncodedIndex_)>::max();
}

// ============================================================================
// SSHash-backed getRefPos implementations
// ============================================================================

namespace {
inline pufferfish::util::ProjectedHits makeEmptyHit(uint32_t k) {
  using IterT = const pufferfish::util::Position*;
  return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k,
          core::range<IterT>{}};
}
} // anonymous namespace

auto PufferfishIndex::getRefPos(CanonicalKmer& mer, pufferfish::util::QueryCache& qc)
    -> pufferfish::util::ProjectedHits {
  using IterT = const pufferfish::util::Position*;

  if (useSSHash_) {
    // SSHash lookup path: use the canonical uint64 kmer representation
    sshash::default_kmer_t sshash_kmer;
    sshash_kmer.bits = mer.getCanonicalWord();
    auto res = dict_.lookup(sshash_kmer, /*check_reverse_complement=*/false);

    if (res.kmer_id == sshash::constants::invalid_uint64) {
      return makeEmptyHit(k_);
    }

    uint64_t string_id = res.string_id;
    uint32_t relPos = static_cast<uint32_t>(res.kmer_id_in_string);
    uint64_t globalPos = res.string_begin + relPos;
    uint64_t contigStart = res.string_begin;
    uint64_t contigEnd = res.string_end;  // exclusive
    uint64_t clen = contigEnd - contigStart;

    // Derive orientation from SSHash's lookup result directly,
    // avoiding a 2k-bit read from the internal bit_vector.
    // forward_orientation means the stored k-mer matches our query (the canonical word).
    // hitFW (IDENTITY_MATCH) means the stored k-mer equals fwWord().
    bool fwIsCanonical = (mer.fwWord() <= mer.rcWord());
    bool hitFW = (res.kmer_orientation == sshash::constants::forward_orientation) == fwIsCanonical;

    // Cache contig table range to avoid repeated Elias-Fano lookups
    // when consecutive k-mers land on the same contig.
    core::range<IterT> contigIterRange;
    if (string_id == qc.prevRank) {
      contigIterRange = core::range<IterT>(qc.rangeBegin, qc.rangeEnd);
    } else {
      contigIterRange = contigTable_.contigRange(string_id);
      qc.prevRank = string_id;
      qc.contigStart = contigStart;
      qc.contigEnd = contigEnd - 1;  // QueryCache uses inclusive end
      qc.rangeBegin = contigIterRange.begin();
      qc.rangeEnd = contigIterRange.end();
    }

    return {static_cast<uint32_t>(string_id),
            globalPos,
            relPos,
            hitFW,
            static_cast<uint32_t>(clen),
            k_,
            contigIterRange};
  }

  // Legacy BooPHF path
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = const_cast<const pos_vector_t&>(pos_)[res];
    uint64_t fk = seq_.get_int(2*pos, 2*k_);

    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      uint64_t rank = 0;
      uint64_t sp = 0;
      uint64_t contigEnd = 0;

      const bool same_contig = ((qc.contigStart <= pos) and (pos <= qc.contigEnd));
      if (same_contig) {
        rank = qc.prevRank;
        sp = qc.contigStart;
        contigEnd = qc.contigEnd;
      } else {
        rank = rankSelDict.rank(pos);
        sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank - 1)) + 1;
        contigEnd = rankSelDict.select(rank);
        qc.prevRank = rank;
        qc.contigStart = sp;
        qc.contigEnd = contigEnd;
      }
      auto contigIterRange = contigRange(rank);
      uint32_t relPos = static_cast<uint32_t>(pos - sp);
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);
      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      return {static_cast<uint32_t>(rank),
              pos,
              relPos,
              hitFW,
              static_cast<uint32_t>(clen),
              k_,
              contigIterRange};
    } else {
      return makeEmptyHit(k_);
    }
  }

  return makeEmptyHit(k_);
}

auto PufferfishIndex::getRefPos(CanonicalKmer& mer, pufferfish::util::QueryCache& qc,
                                sshash_streaming_query_t& sq, char const* kmer_str)
    -> pufferfish::util::ProjectedHits {
  using IterT = const pufferfish::util::Position*;

  // Use SSHash streaming query for sequential k-mer lookups
  auto res = sq.lookup(kmer_str);

  if (res.kmer_id == sshash::constants::invalid_uint64) {
    return makeEmptyHit(k_);
  }

  uint64_t string_id = res.string_id;
  // kmer_offset is only set during seed(), not updated during extensions.
  // Compute globalPos from string_begin + kmer_id_in_string which are always current.
  uint32_t relPos = static_cast<uint32_t>(res.kmer_id_in_string);
  uint64_t globalPos = res.string_begin + relPos;
  uint64_t clen = res.string_end - res.string_begin;

  // Derive orientation from SSHash's lookup result directly.
  bool fwIsCanonical = (mer.fwWord() <= mer.rcWord());
  bool hitFW = (res.kmer_orientation == sshash::constants::forward_orientation) == fwIsCanonical;

  // Cache contig table range to avoid repeated Elias-Fano lookups.
  core::range<IterT> contigIterRange;
  if (string_id == qc.prevRank) {
    contigIterRange = core::range<IterT>(qc.rangeBegin, qc.rangeEnd);
  } else {
    contigIterRange = contigTable_.contigRange(string_id);
    qc.prevRank = string_id;
    qc.contigStart = res.string_begin;
    qc.contigEnd = res.string_end - 1;
    qc.rangeBegin = contigIterRange.begin();
    qc.rangeEnd = contigIterRange.end();
  }

  return {static_cast<uint32_t>(string_id),
          globalPos,
          relPos,
          hitFW,
          static_cast<uint32_t>(clen),
          k_,
          contigIterRange};
}

auto PufferfishIndex::getRefPos(CanonicalKmer& mer) -> pufferfish::util::ProjectedHits {
  if (useSSHash_) {
    sshash::default_kmer_t sshash_kmer;
    sshash_kmer.bits = mer.getCanonicalWord();
    auto res = dict_.lookup(sshash_kmer, /*check_reverse_complement=*/false);

    if (res.kmer_id == sshash::constants::invalid_uint64) {
      return makeEmptyHit(k_);
    }

    uint64_t string_id = res.string_id;
    uint32_t relPos = static_cast<uint32_t>(res.kmer_id_in_string);
    uint64_t globalPos = res.string_begin + relPos;
    uint64_t clen = res.string_end - res.string_begin;

    bool fwIsCanonical = (mer.fwWord() <= mer.rcWord());
    bool hitFW = (res.kmer_orientation == sshash::constants::forward_orientation) == fwIsCanonical;

    auto contigIterRange = contigTable_.contigRange(string_id);

    return {static_cast<uint32_t>(string_id),
            globalPos,
            relPos,
            hitFW,
            static_cast<uint32_t>(clen),
            k_,
            contigIterRange};
  }

  // Legacy BooPHF path
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = const_cast<const pos_vector_t&>(pos_)[res];
    uint64_t fk = seq_.get_int(2*pos, 2*k_);
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      auto rank = rankSelDict.rank(pos);
      auto contigIterRange = contigRange(rank);
      uint64_t sp =
          (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank - 1)) + 1;
      uint64_t contigEnd = rankSelDict.select(rank);
      uint32_t relPos = static_cast<uint32_t>(pos - sp);
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);
      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      return {static_cast<uint32_t>(rank),
              pos,
              relPos,
              hitFW,
              static_cast<uint32_t>(clen),
              k_,
             contigIterRange};
    } else {
      return makeEmptyHit(k_);
    }
  }

  return makeEmptyHit(k_);
}

// ============================================================================
// SSHash-backed overrides for base class methods
// ============================================================================

CanonicalKmer PufferfishIndex::getStartKmer(uint64_t cid) {
  if (useSSHash_) {
    CanonicalKmer::k(k_);
    CanonicalKmer kb;
    auto [begin, end] = dict_.string_offsets(cid);
    uint64_t fk = dict_.strings().get_bits(2 * begin, twok_);
    kb.fromNum(fk);
    return kb;
  }
  // Legacy path: delegate to base class implementation
  CanonicalKmer::k(k_);
  CanonicalKmer kb;
  uint64_t sp = (cid == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(cid)) + 1;
  uint64_t fk = seq_.get_int(2*sp, 2*k_);
  kb.fromNum(fk);
  return kb;
}

CanonicalKmer PufferfishIndex::getEndKmer(uint64_t cid) {
  if (useSSHash_) {
    CanonicalKmer::k(k_);
    CanonicalKmer kb;
    auto [begin, end] = dict_.string_offsets(cid);
    // end is exclusive; last kmer starts at (end - k_)
    uint64_t fk = dict_.strings().get_bits(2 * (end - k_), twok_);
    kb.fromNum(fk);
    return kb;
  }
  // Legacy path
  CanonicalKmer::k(k_);
  CanonicalKmer kb;
  uint64_t contigEnd = rankSelDict.select(cid + 1);
  uint64_t fk = seq_.get_int(2*(contigEnd - k_ + 1), 2*k_);
  kb.fromNum(fk);
  return kb;
}

uint32_t PufferfishIndex::getContigLen(uint64_t rank) {
  if (useSSHash_) {
    auto [begin, end] = dict_.string_offsets(rank);
    return static_cast<uint32_t>(end - begin);
  }
  // Legacy path
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank)) + 1;
  uint64_t contigEnd = rankSelDict.select(rank + 1);
  return static_cast<uint32_t>(contigEnd - sp + 1);
}

uint64_t PufferfishIndex::getGlobalPos(uint64_t rank) {
  if (useSSHash_) {
    auto [begin, end] = dict_.string_offsets(rank);
    return begin;
  }
  // Legacy path
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict.select(rank)) + 1;
  return sp;
}

std::string PufferfishIndex::getSeqStr(size_t start, int64_t length, bool isFw) {
  if (length <= 0) { return ""; }
  std::string outstr(length, 'N');
  uint64_t c = isFw ? 0 : static_cast<uint64_t>(length - 1);

  uint64_t bucket_offset = start * 2;
  auto len_on_vector = length * 2;
  int32_t toFetch = len_on_vector;
  while (toFetch > 0) {
    uint32_t len = (toFetch >= 64) ? 64 : toFetch;
    toFetch -= len;
    uint64_t word = getSeqBits(bucket_offset, len);
    if (isFw) {
      for (uint32_t i = 0; i < len; i += 2) {
        uint8_t next_bits = ((word >> i) & 0x03);
        outstr[c++] = "ACGT"[next_bits];
      }
    } else {
      for (uint32_t i = 0; i < len; i += 2) {
        uint8_t next_bits = ((word >> i) & 0x03);
        outstr[c++] = "TGCA"[next_bits];
      }
    }
    bucket_offset += len;
  }
  if (!isFw) {
    std::reverse(outstr.begin(), outstr.end());
  }
  return outstr;
}

PufferfishIndex::~PufferfishIndex() { }
