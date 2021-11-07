//
// Created by Andrey Bzikadze on 10/19/21.
//

#pragma once

#include <vector>
#include "sequences/contigs.hpp"

namespace tools::common::coverage_utils {

    double get_coverage(const std::vector<Contig> & contigs, const std::vector<Contig> & readset);

} // End namespace tools::common::coverage_utils