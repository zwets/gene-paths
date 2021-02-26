/* graph-test.cpp
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
#include "graph.h"

using namespace gfa;

namespace {

TEST(graph_test, empty_gfa) {
    graph gfa;
    ASSERT_EQ(gfa.segs.size(), 0);
    ASSERT_EQ(gfa.arcs.size(), 0);
}

TEST(graph_test, add_1_seg) {
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

TEST(graph_test, add_2_seg) {
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

TEST(graph_test, add_dup_seg) {
    graph gfa;
    seg s1; s1.len = 4; s1.data = "ACGT"; s1.name = "s1";
    gfa.add_seg(s1);
    ASSERT_EXIT( gfa.add_seg(s1);,
            testing::ExitedWithCode(1), 
            ": error: duplicate segment name: s1");
}

TEST(graph_test, add_len_wrong) {
    graph gfa;
    seg s1; s1.len = 4; s1.data = "ACG"; s1.name = "s1";
    ASSERT_EXIT( gfa.add_seg(s1);,
            testing::ExitedWithCode(1), 
            ": error: segment length in GFA \\(4\\) differs from FASTA \\(3\\) for seqid s1");
}

static seg SEG1 = { 4, "s1", "ACGT" };
static seg SEG2 = { 9, "s2", "TAGCATACG" };
static seg SEG3 = { 5, "s3", "CATTA" };
static seg SEG4 = { 8, "s4", "GCGCAATT" };

TEST(graph_test, add_edge) {
    graph gfa;
    gfa.add_seg(SEG1);                         // ACGT
    gfa.add_seg(SEG2);                         //  CGTATGCTA
    gfa.add_edge("s1+", 1, 4, "s2-", 6, 9);
    ASSERT_EQ(gfa.segs.size(), 2);
    ASSERT_EQ(gfa.arcs.size(), 4);

    auto s1 = gfa.segs[0];
    auto s2 = gfa.segs[1];
    ASSERT_EQ(s1.len, 4);
    ASSERT_EQ(s2.len, 9);

    ASSERT_EQ(gfa.arcs[0].v_lv, 0L<<32 | 1);
    ASSERT_EQ(gfa.arcs[0].w_lw, 3L<<32 | 0);

    ASSERT_EQ(gfa.arcs[1].v_lv, 0L<<32 | 4);
    ASSERT_EQ(gfa.arcs[1].w_lw, 3L<<32 | 3);

    ASSERT_EQ(gfa.arcs[2].v_lv, 2L<<32 | 6);
    ASSERT_EQ(gfa.arcs[2].w_lw, 1L<<32 | 0);

    ASSERT_EQ(gfa.arcs[3].v_lv, 2L<<32 | 9);
    ASSERT_EQ(gfa.arcs[3].w_lw, 1L<<32 | 3);

}

TEST(graph_test, add_blunt_edge) {
    graph gfa;
    gfa.add_seg(SEG1);                         //          ACGT
    gfa.add_seg(SEG2);                         // CGTATGCTA
    gfa.add_edge("s1-", 4, 4, "s2+", 9, 9);
    ASSERT_EQ(gfa.segs.size(), 2);
    ASSERT_EQ(gfa.arcs.size(), 2);

    ASSERT_EQ(gfa.arcs[0].v_lv, 0L<<32 | 4);   // v is s1-
    ASSERT_EQ(gfa.arcs[0].w_lw, 3L<<32 | 0);

    ASSERT_EQ(gfa.arcs[1].v_lv, 2L<<32 | 9);
    ASSERT_EQ(gfa.arcs[1].w_lw, 1L<<32 | 0);
}

TEST(graph_test, vtx_iter) {
    graph gfa;
    gfa.add_seg(SEG1);
    gfa.add_seg(SEG2);
    gfa.add_seg(SEG3);
    gfa.add_edge("s1+", 1, 4, "s2-", 5, 9); // s1+ .[--) s2- (---].....   and s2+ .....[---) s1- (--].
    gfa.add_edge("s2-", 0, 0, "s3+", 0, 0); // s2- .........) s3+ (.....  and s3- .....) s2+ (.........
    gfa.add_edge("s3+", 4, 5, "s1+", 0, 1); // s3+ ....[) s1+ (]...       and s1- ...[) s3- (]...

    auto afv = gfa.arcs_from_vtx(0); // s1+
    ASSERT_EQ(std::distance(gfa.arcs.cbegin(), afv.first), 0);
    ASSERT_EQ(std::distance(afv.first, afv.second), 2);

    auto nxt = afv.first;
    ASSERT_EQ(nxt->v_lv, 1);
    ASSERT_EQ(nxt->w_lw, 3L<<32|0);

    nxt = ++afv.first;
    ASSERT_EQ(nxt->v_lv, 4);
    ASSERT_EQ(nxt->w_lw, 3L<<32|4);

    afv = gfa.arcs_from_vtx(1); // s1-
    ASSERT_EQ(std::distance(afv.first, afv.second), 2);
    ASSERT_EQ(afv.first->v_lv, 1L<<32|3);
    ASSERT_EQ(afv.first->w_lw, 5L<<32|0);

    afv = gfa.arcs_from_vtx(2); // s2+
    ASSERT_EQ(std::distance(afv.first, afv.second), 2);

    afv = gfa.arcs_from_vtx(3); // s2-
    ASSERT_EQ(std::distance(afv.first, afv.second), 1);

    afv = gfa.arcs_from_vtx(4); // s3+
    ASSERT_EQ(std::distance(afv.first, afv.second), 2);

    afv = gfa.arcs_from_vtx(5); // s3-
    ASSERT_EQ(std::distance(afv.first, afv.second), 1);
    ASSERT_EQ(afv.first->v_lv, 5L<<32|5);
    ASSERT_EQ(afv.first->w_lw, 2L<<32|0);
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
