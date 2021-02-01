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
#include "gfagraph.h"

using namespace gfa;

namespace {

seg SEG1 = { 4, "s1", "ACGT" };
seg SEG2 = { 9, "s2", "TAGCATACG" };

TEST(gfagraph_test, empty_gfa) {
    graph gfa;
    ASSERT_EQ(gfa.segs.size(), 0);
    ASSERT_EQ(gfa.arcs.size(), 0);
}

TEST(gfagraph_test, add_1_seg) {
    graph gfa;
    seg s;
    s.len = 4;
    s.data = "ACGT";
    s.name = "s1";
    gfa.add_seg(s);
    ASSERT_EQ(gfa.segs.size(), 1);
    ASSERT_EQ(gfa.seg_ixs.size(), 1);
    ASSERT_EQ(gfa.seg_ixs[s.name], 0);
}

TEST(gfagraph_test, add_2_seg) {
    graph gfa;
    seg s1; s1.len = 4; s1.data = "ACGT"; s1.name = "s1";
    gfa.add_seg(s1);
    seg s2; s2.len = 5; s2.data = "GATCA"; s2.name = "s2";
    gfa.add_seg(s2);
    ASSERT_EQ(gfa.segs.size(), 2);
    ASSERT_EQ(gfa.seg_ixs.size(), 2);
    ASSERT_EQ(gfa.seg_ixs[s1.name], 0);
    ASSERT_EQ(gfa.seg_ixs[s2.name], 1);
}

TEST(gfagraph_test, add_dup_seg) {
    graph gfa;
    seg s1; s1.len = 4; s1.data = "ACGT"; s1.name = "s1";
    gfa.add_seg(s1);
    ASSERT_EXIT( gfa.add_seg(s1);,
            testing::ExitedWithCode(1), 
            ": error: duplicate segment name: s1");
}

TEST(gfagraph_test, add_len_wrong) {
    graph gfa;
    seg s1; s1.len = 4; s1.data = "ACG"; s1.name = "s1";
    ASSERT_EXIT( gfa.add_seg(s1);,
            testing::ExitedWithCode(1), 
            ": error: segment length in GFA \\(4\\) differs from FASTA \\(3\\) for seqid s1");
}

TEST(gfagraph_test, add_edge) {
    graph gfa;
    gfa.add_seg(SEG1);
    gfa.add_seg(SEG2);
    gfa.add_edge("s1+", 2, 4, "s2-", 6, 8);
    ASSERT_EQ(gfa.segs.size(), 2);
    ASSERT_EQ(gfa.arcs.size(), 2);

    auto s1 = gfa.segs[0];
    auto s2 = gfa.segs[1];
    ASSERT_EQ(s1.len, 4);
    ASSERT_EQ(s2.len, 9);

    auto a1 = gfa.arcs[0];
    ASSERT_EQ(a1.v_lv, 2);
    ASSERT_EQ(a1.ov, 2);
    ASSERT_EQ(a1.w_lw, std::uint64_t(3)<<32|1);
    ASSERT_EQ(a1.ow, 2);

    auto a2 = gfa.arcs[1];
    ASSERT_EQ(a2.v_lv, std::uint64_t(3)<<32|1);
    ASSERT_EQ(a2.ov, 2);
    ASSERT_EQ(a2.w_lw, 2);
    ASSERT_EQ(a2.ow, 2);
}

} // namespace
  // vim: sts=4:sw=4:ai:si:et
