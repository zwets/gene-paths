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

static seg SEG1 = { 4, "s1", "ACGT" };
static seg SEG2 = { 9, "s2", "TAGCATACG" };
static seg SEG3 = { 5, "s3", "CATTA" };
static seg SEG4 = { 8, "s4", "GCGCAATT" };

TEST(gfagraph_test, add_edge) {
    graph gfa;
    gfa.add_seg(SEG1);
    gfa.add_seg(SEG2);
    gfa.add_edge("s1+", 1, 4, "s2-", 5, 9);
    ASSERT_EQ(gfa.segs.size(), 2);
    ASSERT_EQ(gfa.arcs.size(), 2);

    auto s1 = gfa.segs[0];
    auto s2 = gfa.segs[1];
    ASSERT_EQ(s1.len, 4);
    ASSERT_EQ(s2.len, 9);

    auto a1 = gfa.arcs[0];
    ASSERT_EQ(a1.v_lv, 1);
    ASSERT_EQ(a1.ov, 3);
    ASSERT_EQ(a1.w, 3);
    ASSERT_EQ(a1.ow, 4);

    auto a2 = gfa.arcs[1];
    ASSERT_EQ(a2.v_lv, std::uint64_t(3^1)<<32|5);
    ASSERT_EQ(a2.w, 0^1);
    ASSERT_EQ(a2.ov, 4);
    ASSERT_EQ(a2.ow, 3);
}

TEST(gfagraph_test, arcs_from_vtx_empty) {
    graph gfa;
    ASSERT_EQ(gfa.segs.size(), 0);
    ASSERT_EQ(gfa.arcs.size(), 0);

    auto arcs_from_vtx = gfa.arcs_from_vtx(0);
    ASSERT_EQ(arcs_from_vtx.first, arcs_from_vtx.second);
    auto arcs_from_v_lv = gfa.arcs_from_v_lv(0);
    ASSERT_EQ(arcs_from_v_lv.first, arcs_from_v_lv.second);
}

TEST(gfagraph_test, arcs_order) {
    graph gfa;
    gfa.add_seg(SEG1);
    gfa.add_seg(SEG2);
    gfa.add_seg(SEG3);
    gfa.add_edge("s1+", 1, 4, "s2-", 5, 9); // s1+ .[--) s2- (---].....   and s2+ .....[---) s1- (--].
    gfa.add_edge("s2-", 0, 0, "s3+", 0, 0); // s2- .........) s3+ (.....  and s3- .....) s2+ (.........
    gfa.add_edge("s3+", 4, 5, "s1+", 0, 1); // s3+ ....[) s1+ (]...       and s1- ...[) s3- (]...

    ASSERT_EQ(gfa.segs.size(), 3);
    ASSERT_EQ(gfa.arcs.size(), 6);

    arc& o = gfa.arcs[0];

    o = gfa.arcs[0]; ASSERT_EQ(o.ov, 3); ASSERT_EQ(o.ow, 4); ASSERT_EQ(o.w, 3); ASSERT_EQ(o.v_lv, (std::uint64_t(0)<<32)|1); // s1+ s2-
    o = gfa.arcs[1]; ASSERT_EQ(o.ov, 1); ASSERT_EQ(o.ow, 1); ASSERT_EQ(o.w, 5); ASSERT_EQ(o.v_lv, (std::uint64_t(1)<<32)|3); // s1- s3-
    o = gfa.arcs[2]; ASSERT_EQ(o.ov, 4); ASSERT_EQ(o.ow, 3); ASSERT_EQ(o.w, 1); ASSERT_EQ(o.v_lv, (std::uint64_t(2)<<32)|5); // s2+ s1-
    o = gfa.arcs[3]; ASSERT_EQ(o.ov, 0); ASSERT_EQ(o.ow, 0); ASSERT_EQ(o.w, 4); ASSERT_EQ(o.v_lv, (std::uint64_t(3)<<32)|9); // s2- s3+
    o = gfa.arcs[4]; ASSERT_EQ(o.ov, 1); ASSERT_EQ(o.ow, 1); ASSERT_EQ(o.w, 0); ASSERT_EQ(o.v_lv, (std::uint64_t(4)<<32)|4); // s3+ s1+
    o = gfa.arcs[5]; ASSERT_EQ(o.ov, 0); ASSERT_EQ(o.ow, 0); ASSERT_EQ(o.w, 2); ASSERT_EQ(o.v_lv, (std::uint64_t(5)<<32)|5); // s3- s2+
}

TEST(gfagraph_test, vtx_iter) {
    graph gfa;
    gfa.add_seg(SEG1);
    gfa.add_seg(SEG2);
    gfa.add_seg(SEG3);
    gfa.add_edge("s1+", 1, 4, "s2-", 5, 9); // s1+ .[--) s2- (---].....   and s2+ .....[---) s1- (--].
    gfa.add_edge("s2-", 0, 0, "s3+", 0, 0); // s2- .........) s3+ (.....  and s3- .....) s2+ (.........
    gfa.add_edge("s3+", 4, 5, "s1+", 0, 1); // s3+ ....[) s1+ (]...       and s1- ...[) s3- (]...

    const auto afv0 = gfa.arcs_from_vtx(0);
    ASSERT_EQ(std::distance(gfa.arcs.cbegin(), afv0.first), 0);
    ASSERT_EQ(std::distance(afv0.first, afv0.second), 1);

    const auto afv1 = gfa.arcs_from_vtx(1);
    ASSERT_EQ(std::distance(afv1.first, afv1.second), 1);

    const auto afv2 = gfa.arcs_from_vtx(2);
    ASSERT_EQ(std::distance(afv2.first, afv2.second), 1);

    const auto afv3 = gfa.arcs_from_vtx(3);
    ASSERT_EQ(std::distance(gfa.arcs.cbegin(), afv3.first), 3);
    ASSERT_EQ(std::distance(afv3.first, afv3.second), 1);
}

//static seg SEG1 = { 4, "s1", "ACGT" };
//static seg SEG2 = { 9, "s2", "TAGCATACG" };
//static seg SEG3 = { 5, "s3", "CATTA" };
//static seg SEG4 = { 8, "s4", "CTATAATT" };

TEST(gfagraph_test, v_vl_iter) {
    graph gfa;
    gfa.add_seg(SEG1);
    gfa.add_seg(SEG2);
    gfa.add_seg(SEG3);
    gfa.add_seg(SEG4);
    gfa.add_edge("s1+", 1, 4, "s2-", 5, 9); // s1+ .[--) s2- (---].....     and s2+ .....[---) s1- (--].
    gfa.add_edge("s2-", 0, 0, "s3+", 0, 0); // s2- .........) s3+ (.....    and s3- .....) s2+ (.........
    gfa.add_edge("s2-", 0, 3, "s4+", 0, 3); // s2- ......[..) s4+ (..]..... and s4- .....[..) s2+ (..]......
    gfa.add_edge("s3+", 4, 5, "s1+", 0, 1); // s3+ ....[) s1+ (]...         and s1- ...[) s3- (]...

    arc& o = gfa.arcs[0];

    o = gfa.arcs[0]; ASSERT_EQ(o.ov, 3); ASSERT_EQ(o.ow, 4); ASSERT_EQ(o.w, 3); ASSERT_EQ(o.v_lv, (std::uint64_t(0)<<32)|1); // s1+ s2-
    o = gfa.arcs[1]; ASSERT_EQ(o.ov, 1); ASSERT_EQ(o.ow, 1); ASSERT_EQ(o.w, 5); ASSERT_EQ(o.v_lv, (std::uint64_t(1)<<32)|3); // s1- s3-
    o = gfa.arcs[2]; ASSERT_EQ(o.ov, 4); ASSERT_EQ(o.ow, 3); ASSERT_EQ(o.w, 1); ASSERT_EQ(o.v_lv, (std::uint64_t(2)<<32)|5); // s2+ s1-
    o = gfa.arcs[3]; ASSERT_EQ(o.ov, 3); ASSERT_EQ(o.ow, 3); ASSERT_EQ(o.w, 6); ASSERT_EQ(o.v_lv, (std::uint64_t(3)<<32)|6); // s2- s4+
    o = gfa.arcs[4]; ASSERT_EQ(o.ov, 0); ASSERT_EQ(o.ow, 0); ASSERT_EQ(o.w, 4); ASSERT_EQ(o.v_lv, (std::uint64_t(3)<<32)|9); // s2- s3+
    o = gfa.arcs[5]; ASSERT_EQ(o.ov, 1); ASSERT_EQ(o.ow, 1); ASSERT_EQ(o.w, 0); ASSERT_EQ(o.v_lv, (std::uint64_t(4)<<32)|4); // s3+ s1+
    o = gfa.arcs[6]; ASSERT_EQ(o.ov, 0); ASSERT_EQ(o.ow, 0); ASSERT_EQ(o.w, 2); ASSERT_EQ(o.v_lv, (std::uint64_t(5)<<32)|5); // s3- s2+
    o = gfa.arcs[7]; ASSERT_EQ(o.ov, 3); ASSERT_EQ(o.ow, 3); ASSERT_EQ(o.w, 2); ASSERT_EQ(o.v_lv, (std::uint64_t(7)<<32)|5); // s4- s2+

    const auto itp0 = gfa.arcs_from_vtx(3);
    ASSERT_EQ(std::distance(gfa.arcs.cbegin(), itp0.first), 3);
    ASSERT_EQ(std::distance(itp0.first, itp0.second), 2);

    const auto itp1 = gfa.arcs_from_v_lv(std::uint64_t(3)<<32|1);
    ASSERT_EQ(std::distance(gfa.arcs.cbegin(), itp1.first), 3);
    ASSERT_EQ(std::distance(itp1.first, itp1.second), 2);

    const auto itp2 = gfa.arcs_from_v_lv(std::uint64_t(3)<<32|6);
    ASSERT_EQ(std::distance(gfa.arcs.cbegin(), itp2.first), 3);
    ASSERT_EQ(std::distance(itp2.first, itp2.second), 2);

    const auto itp3 = gfa.arcs_from_v_lv(std::uint64_t(3)<<32|7);
    ASSERT_EQ(std::distance(gfa.arcs.cbegin(), itp3.first), 4);
    ASSERT_EQ(std::distance(itp3.first, itp3.second), 1);

    const auto itp4 = gfa.arcs_from_v_lv(std::uint64_t(3)<<32|9);
    ASSERT_EQ(std::distance(gfa.arcs.cbegin(), itp4.first), 4);
    ASSERT_EQ(std::distance(itp4.first, itp4.second), 1);

    const auto itp5 = gfa.arcs_from_v_lv(std::uint64_t(3)<<32|10);
    ASSERT_EQ(std::distance(itp5.first, itp5.second), 0);
}

} // namespace
  // vim: sts=4:sw=4:ai:si:et
