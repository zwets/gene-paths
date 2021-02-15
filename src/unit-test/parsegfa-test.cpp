/* utils-test.cpp
 * 
 * Copyright (C) 2021  Marco van Zwetselaar <io@zwets.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <gtest/gtest.h>
#include <sstream>
#include <fstream>
#include "parsegfa.h"
#include "gfagraph.h"
#include "utils.h"

using namespace gene_paths;

namespace {

TEST(parsegfa_test, read_gfa) {

    std::ifstream gfa_file("data/with_seqs.gfa");
    ASSERT_TRUE(gfa_file);
    
    gfa::graph gfa = parse_gfa(gfa_file);

    ASSERT_EQ(gfa.segs.size(), 9);
}

TEST(parsegfa_test, read_gfa_and_fna) {

    std::ifstream gfa_file("data/without_seqs.gfa");
    ASSERT_TRUE(gfa_file);

    std::ifstream fna_file("data/seqs.fna");
    ASSERT_TRUE(fna_file);

    gfa::graph gfa = parse_gfa(gfa_file, fna_file);

    ASSERT_EQ(gfa.segs.size(), 9);
}

TEST(parsegfa_test, read_gfa_string) {

    std::istringstream s_gfa("H\tVN:Z:2.0\nS\t1\t4\t*\n");
    std::istringstream s_fna(">1\nACGT\n");

    gfa::graph gfa = parse_gfa(s_gfa, s_fna);
    ASSERT_EQ(gfa.segs.size(), 1);
}

TEST(parsegfa_test, read_gfa_mismatch_fna) {

    std::istringstream s_gfa("H\tVN:Z:2.0\nS\t1\t4\t*\n");
    std::istringstream s_fna(">1\nACG\n");

    ASSERT_EXIT( parse_gfa(s_gfa, s_fna);,
            testing::ExitedWithCode(1), 
            ": error: segment length in GFA \\(4\\) differs from FASTA \\(3\\) for seqid 1");
}

TEST(parsegfa_test, read_gfa_and_edge) {

    std::istringstream s_gfa("H\tVN:Z:2.0\n"
        "S\ts1\t4\tACGT\n"
        "S\ts2\t9\tTAGCATACG\n"
        "E\t*\ts1+\ts2-\t1\t4$\t5\t9\t*\n");

    gfa::graph gfa = parse_gfa(s_gfa);
    ASSERT_EQ(gfa.segs.size(), 2);
    ASSERT_EQ(gfa.arcs.size(), 8);

    auto s1 = gfa.segs[0];
    auto s2 = gfa.segs[1];
    ASSERT_EQ(s1.len, 4);
    ASSERT_EQ(s2.len, 9);

    // We have:
    // - v (v+) = 0, w (w-) = 3, v' (v-) = 0, w' (w+) = 2
    // - v_lv 0<<32|1, w_lw 3<<32|0, v'_lv 1<<32|0, w'_lw 2<<32|5
    // - v_lv+ov 0<<32|4, w_lw+ow 3<<32|4, v'_lv+ov 1<<32|3, w'_lw+ow 2<<32|9

    // 0: v_lv -> w_lw
    ASSERT_EQ(gfa.arcs[0].v_lv, 0L<<32|1);
    ASSERT_EQ(gfa.arcs[0].w_lw, 3L<<32|0);

    // 1: v_lv+ov -> w_lw+ow
    ASSERT_EQ(gfa.arcs[1].v_lv, 0L<<32|4);
    ASSERT_EQ(gfa.arcs[1].w_lw, 3L<<32|4);

    // 2: v'_lv -> w'_lw
    ASSERT_EQ(gfa.arcs[2].v_lv, 1L<<32|0);
    ASSERT_EQ(gfa.arcs[2].w_lw, 2L<<32|5);

    // 3: v'_lv+ov -> w'_lw+ow
    ASSERT_EQ(gfa.arcs[3].v_lv, 1L<<32|3);
    ASSERT_EQ(gfa.arcs[3].w_lw, 2L<<32|9);

    // 4: w'_lw -> v'_lv
    ASSERT_EQ(gfa.arcs[4].v_lv, gfa.arcs[2].w_lw); 
    ASSERT_EQ(gfa.arcs[4].w_lw, gfa.arcs[2].v_lv);

    // 5: w'_lw+ow -> v'_lv+ov
    ASSERT_EQ(gfa.arcs[5].v_lv, gfa.arcs[3].w_lw); 
    ASSERT_EQ(gfa.arcs[5].w_lw, gfa.arcs[3].v_lv);

    // 6: w_lw -> v_lv
    ASSERT_EQ(gfa.arcs[6].v_lv, gfa.arcs[0].w_lw);
    ASSERT_EQ(gfa.arcs[6].w_lw, gfa.arcs[0].v_lv);

    // 7: w_lw+ow -> v_lv+ov
    ASSERT_EQ(gfa.arcs[7].v_lv, gfa.arcs[1].w_lw);
    ASSERT_EQ(gfa.arcs[7].w_lw, gfa.arcs[1].v_lv);
}

} // namespace
  // vim: sts=4:sw=4:ai:si:et
