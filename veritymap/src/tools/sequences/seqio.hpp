#pragma once

#include <filesystem>
#include <functional>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "common/string_utils.hpp"
#include "contigs.hpp"
#include "stream.hpp"

namespace io {

typedef std::vector<std::filesystem::path> Library;

template<class Reader>
class ContigIterator {
 private:
  Reader &reader;
  bool isend;

 public:
  typedef StringContig value_type;

  ContigIterator(Reader &_reader, bool _isend) : reader(_reader), isend(_isend) {
    if (reader.eof()) {
      isend = true;
    }
  }

  void operator++() {
    reader.inner_read();
    if (reader.eof()) {
      isend = true;
    }
  }

  StringContig operator*() { return reader.get(); }

  bool operator==(const ContigIterator &other) const { return isend == other.isend; }

  bool operator!=(const ContigIterator &other) const { return isend != other.isend; }
};

template<class Reader>
class SeqIterator : public std::iterator<std::forward_iterator_tag, Sequence, size_t, Sequence *, Sequence &> {
 private:
  Reader &reader;
  bool isend;

 public:
  SeqIterator(Reader &_reader, bool _isend) : reader(_reader), isend(_isend) {
    if (reader.eof()) {
      isend = true;
    }
  }

  void operator++() {
    reader.inner_read();
    if (reader.eof()) {
      isend = true;
    }
  }

  const Sequence &operator*() { return reader.get().seq; }

  bool operator==(const SeqIterator &other) const { return isend == other.isend; }
  bool operator!=(const SeqIterator &other) const { return isend != other.isend; }
};

//    TODO: Deal with corrupted files, comments in read names
class SeqReader {
 public:
  typedef ContigIterator<SeqReader> Iterator;

 private:
  void choose_next_pos(size_t start) {
    cur_start = start;
    if (next.size() > start + 2 * min_read_size - overlap) {
      cur_end = start + min_read_size;
    } else {
      cur_end = next.size();
    }
  }

  void inner_read() {
    if (cur_end > 0 && cur_end < next.size()) {
      choose_next_pos(cur_end - overlap);
      return;
    }
    while (stream != nullptr) {
      std::string id, seq;
      std::getline(*stream, id);
      if (!id.empty()) {
        std::stringstream ss;
        while (stream->peek() != EOF && stream->peek() != '>' && stream->peek() != '+') {
          std::getline(*stream, seq);
          trim(seq);
          if (!seq.empty())
            ss << seq;
        }
        next = {ss.str(), std::move(trim(id.substr(1, id.size() - 1)))};
        choose_next_pos(0);
        cur_start = 0;
        if (fastq) {
          std::getline(*stream, id);
          size_t qlen = 0;
          while (!stream->eof() && qlen < next.size()) {
            std::getline(*stream, seq);
            trim(seq);
            qlen += seq.size();
            if (seq.empty())
              break;
          }
        }
        return;
      }
      nextFile();
    }
    next = StringContig();
    cur_start = 0;
  }

  void nextFile() {
    delete stream;
    if (file_it == lib.end()) {
      stream = nullptr;
    } else {
      std::filesystem::path file_name = *file_it;
      if (endsWith(file_name, ".gz")) {
        stream = new gzstream::igzstream(file_name.c_str());
        fastq = endsWith(file_name, "fastq.gz") or endsWith(file_name, "fq.gz");
      } else {
        stream = new std::ifstream(file_name);
        fastq = endsWith(file_name, "fastq") or endsWith(file_name, "fq");
      }
      ++file_it;
    }
  }

  const Library lib;
  Library::const_iterator file_it;
  std::istream *stream{};
  bool fastq{};
  size_t min_read_size;
  size_t overlap;
  StringContig next{};
  size_t cur_start = 0;
  size_t cur_end = 0;

 public:
  friend class ContigIterator<SeqReader>;
  friend class SeqIterator<SeqReader>;

  explicit SeqReader(Library _lib, size_t _min_read_size = size_t(-1) / 2, size_t _overlap = size_t(-1) / 8)
      : lib(std::move(_lib)),
        file_it(lib.begin()),
        min_read_size(_min_read_size),
        overlap(_overlap) {
    VERIFY(min_read_size >= overlap * 2);
    reset();
  }

  explicit SeqReader(const std::filesystem::path &file_name, size_t _min_read_size = size_t(-1) / 2,
                     size_t _overlap = size_t(-1) / 8)
      : SeqReader(Library({file_name}), _min_read_size, _overlap) {}

  void reset() {
    file_it = lib.begin();
    nextFile();
    cur_start = 0;
    cur_end = 0;
    inner_read();
  }

  SeqReader(SeqReader &&other) = default;

  //        static SeqReader CompressingReader(const std::string &file_name) {
  //            return SeqReader(file_name, [](std::string &s) {compress_inplace(s);});
  //        }
  //
  //        static SeqReader CompressingReader(const Library &lib_) {
  //            return SeqReader(lib_, [](std::string &s) {compress_inplace(s);});
  //        }

  ContigIterator<SeqReader> begin() { return {*this, false}; }

  ContigIterator<SeqReader> end() { return {*this, true}; }

  SeqIterator<SeqReader> seqbegin() { return {*this, false}; }

  SeqIterator<SeqReader> seqend() { return {*this, true}; }

  StringContig get() {
    StringContig tmp;
    if (cur_start != 0 || cur_end != next.size()) {
      return StringContig(next.seq.substr(cur_start, cur_end - cur_start), next.id + "_" + std::to_string(cur_start));
    } else {
      return std::move(next);
    }
  }

  StringContig read() {
    StringContig tmp = get();
    inner_read();
    return std::move(tmp);
  }

  std::vector<StringContig> readAll() {
    std::vector<StringContig> res;
    while (!eof()) { res.emplace_back(read()); }
    return std::move(res);
  }

  std::vector<Contig> readAllContigs() {
    std::vector<Contig> res;
    while (!eof()) { res.emplace_back(read().makeContig()); }
    return std::move(res);
  }

  bool eof() { return next.isNull(); }

  ~SeqReader() { delete stream; }
};

}// namespace io