#pragma once
//
// Created by anton on 7/20/20.
//

#include <deque>
#include <sequences/sequence.hpp>

template<typename T, typename U>
T pow_custom(T base, U p) {
  if (p == 0)
    return 1;
  T tmp = pow_custom(base, p / 2);
  if (p % 2 == 1)
    return base * tmp * tmp;
  else
    return tmp * tmp;
}

template<typename htype>
class RollingHash {
 public:
  const size_t k;
  const htype hbase;
  const htype kpow;
  const htype inv;

  RollingHash(size_t _k, htype _hbase) : k(_k),
                                         hbase(_hbase),
                                         kpow(pow_custom<htype, htype>(hbase, k - 1)),
                                         inv(pow_custom<htype, htype>(hbase, (htype(1u) << (sizeof(htype) * 8u - 1u)) - 1u)) {
    VERIFY(inv * hbase == htype(1));
  }

  RollingHash extensionHash() const {
    return RollingHash(k + 1, hbase);
  }

  htype hash(const Sequence &seq, size_t pos) const {
    htype hash = 0;
    for (size_t i = pos; i < pos + k; i++) {
      hash = hash * hbase + seq[i];
    }
    return hash;
  }

  htype extendRight(const Sequence &seq, size_t pos, htype hash, unsigned char c) const {
    return hash * hbase + c;
  }

  htype extendLeft(const Sequence &seq, size_t pos, htype hash, unsigned char c) const {
    return hash + c * kpow * hbase;
  }

  htype shiftRight(const Sequence &seq, size_t pos, htype hash, unsigned char c) const {
    return (hash - kpow * seq[pos]) * hbase + c;
  }

  htype shiftLeft(const Sequence &seq, size_t pos, htype hash, unsigned char c) const {
    return (hash - seq[pos + k - 1]) * inv + c * kpow;
  }

  htype next(const Sequence &seq, size_t pos, htype hash) const {
    return shiftRight(seq, pos, hash, seq[pos + k]);
  }

  htype prev(const Sequence &seq, size_t pos, htype hash) const {
    return shiftLeft(seq, pos, hash, seq[pos - 1]);
  }

  bool hasNext(const Sequence &seq, size_t pos) const {
    return pos + k < seq.size();
  }

  bool hasPrev(const Sequence &seq, size_t pos) const {
    return pos > 0;
  }
};

template<typename htype>
class KWH {
 private:
  KWH(const RollingHash<htype> &_hasher, const Sequence &_seq, size_t _pos, htype _fhash, htype _rhash) : hasher(_hasher),
                                                                                                          seq(_seq),
                                                                                                          pos(_pos),
                                                                                                          fhash(_fhash),
                                                                                                          rhash(_rhash) {
  }

  htype fhash;
  htype rhash;
  Sequence seq;

 public:
  const RollingHash<htype> &hasher;
  size_t pos;

  KWH(const RollingHash<htype> &_hasher, const Sequence &_seq, size_t _pos) : hasher(_hasher),
                                                                              seq(_seq),
                                                                              pos(_pos),
                                                                              fhash(_hasher.hash(_seq, _pos)),
                                                                              rhash(_hasher.hash(!_seq, _seq.size() - _pos - _hasher.k)) {
  }

  KWH(const KWH &other) = default;

  Sequence getSeq() const {
    return seq.Subseq(pos, pos + hasher.k);
  }

  KWH<htype> operator!() const {
    return KWH<htype>(hasher, !seq, seq.size() - pos - hasher.k, rhash, fhash);
  }

  htype hash() const {
    return std::min(fhash, rhash);
  }

  htype get_fhash() const {
    return fhash;
  }

  htype get_rhash() const {
    return rhash;
  }

  htype extendRight(unsigned char c) const {
    return std::min(hasher.extendRight(seq, pos, fhash, c), hasher.extendLeft(!seq, seq.size() - pos - hasher.k, rhash, c ^ 3u));
  }

  htype extendLeft(unsigned char c) const {
    return std::min(hasher.extendLeft(seq, pos, fhash, c), hasher.extendRight(!seq, seq.size() - pos - hasher.k, rhash, c ^ 3u));
  }

  KWH next() const {
    return {hasher, seq, pos + 1, hasher.next(seq, pos, fhash), hasher.prev(!seq, seq.size() - pos - hasher.k, rhash)};
  }

  KWH prev() const {
    return {hasher, seq, pos - 1, hasher.prev(seq, pos, fhash), hasher.next(!seq, seq.size() - pos - hasher.k, rhash)};
  }

  bool hasNext() const {
    return hasher.hasNext(seq, pos);
  }

  bool hasPrev() const {
    return hasher.hasPrev(seq, pos);
  }

  KWH &operator=(const KWH &other) {
    if (this == &other)
      return *this;
    seq = other.seq;
    pos = other.pos;
    fhash = other.fhash;
    rhash = other.rhash;
    return *this;
  }

  bool isCanonical() const {
    return fhash < rhash;
  }
};
