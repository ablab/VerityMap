#pragma once

#include "cigar.hpp"

#include "ksw2/ksw2.h"

namespace tandem_mapper::ksw_align {

    cigar_utils::Cigar align(const Sequence &tseq, const int t_st, const int t_en,
                             const Sequence &qseq, const int q_st, const int q_en,
                             const int8_t match_score,
                             const int8_t mis_score,
                             const int8_t gapo,
                             const int8_t gape) {
        // Based on example at https://github.com/lh3/ksw2
        const int8_t a = match_score, b = mis_score < 0 ? mis_score : -mis_score; // a>0 and b<0
        const int8_t mat[25] { a,b,b,b,0,
                               b,a,b,b,0,
                               b,b,a,b,0,
                               b,b,b,a,0,
                               0,0,0,0,0 };
        const int32_t tl { t_en - t_st }, ql { q_en - q_st };
        uint8_t *ts, *qs, c[256];
        ksw_extz_t ez;

        memset(&ez, 0, sizeof(ksw_extz_t));
        memset(c, 4, 256);
        c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
        c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
        ts = (uint8_t*)malloc(tl);
        qs = (uint8_t*)malloc(ql);
        for (size_t i = 0; i < tl; ++i) ts[i] = (uint8_t)tseq[i]; // encode to 0/1/2/3
        for (size_t i = 0; i < ql; ++i) qs[i] = (uint8_t)qseq[i];
        ksw_extz2_sse(nullptr, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 10, 0, &ez);

        cigar_utils::Cigar cigar(ez);

        free(ez.cigar); free(ts); free(qs);

        return cigar;
    }

} // End namespace tandem_mapper::ksw_align
