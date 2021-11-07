//
// Created by Andrey Bzikadze on 2/28/21.
//

#pragma once

namespace tandem_mapper::dna_strand {

enum class Strand {
  forward,
  reverse
};

std::string strand2str(const Strand& strand) {
  return strand == Strand::forward ? "+" : "-";
}

}// namespace tandem_mapper::dna_strand
