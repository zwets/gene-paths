/* targets-test.cpp
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
#include "targets.h"
#include "utils.h"

using namespace gfa;

namespace {

static seg SEG1 = { 10, "SEG1", "CATTAGTACT" };

static graph make_graph()
{
    graph g;
    g.segs.reserve(1 + 3 /* two targets and a terminator */);
    g.arcs.reserve(2 * 4 /* four for each target */);
    g.add_seg(SEG1); // 0  0  1 (seg_ix, v+, v-)
    // ter segment   // 1  2  3
    // target 1      // 2  4  5
    // target 2      // 3  6  7
    return g;
}

// arcs will be
//
// 0_b   to 4_0      SEG1+ to TGT1+ [END]
// 1_L-e to 5_0      SEG1- to TGT1- [END]
// 2_0   to 4_0      TER+  to TGT1+ [START]
// 2_0   to 5_0      TER+  to TGT1- [START]
// 4_e-b to 0_e      TGT1+ to SEG1+ [START]
// 4_e-b to 2_1      TGT1+ to TER+  [END]
// 5_e-b to 1_L-b    TGT1- to SEG1- [START]
// 5_e-b to 2_1      TGT1- to TER+  [END]


TEST(targets_test, start_pos_full) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1+", target::role_t::START);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 2L<<32|0);
    ASSERT_EQ(a.w_lw, 4L<<32|0);
    // 2_0   to 4_0      TER+  to TGT1+ [START]
    ASSERT_EQ(g.arcs.at(0).v_lv, 2L<<32|0);
    ASSERT_EQ(g.arcs.at(0).w_lw, 4L<<32|0);
    // 4_e-b to 0_e      TGT1+ to SEG1+ [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 4L<<32|10);
    ASSERT_EQ(g.arcs.at(1).w_lw, 0L<<32|10);
    ASSERT_EQ(g.segs.at(2).data, SEG1.data);
}

TEST(targets_test, start_pos_part) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1+:2:5", target::role_t::START);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 2L<<32|0);
    ASSERT_EQ(a.w_lw, 4L<<32|0);
    // 2_0   to 4_0      TER+  to TGT1+ [START]
    ASSERT_EQ(g.arcs.at(0).v_lv, 2L<<32|0);
    ASSERT_EQ(g.arcs.at(0).w_lw, 4L<<32|0);
    // 4_e-b to 0_e      TGT1+ to SEG1+ [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 4L<<32|(5-2));
    ASSERT_EQ(g.arcs.at(1).w_lw, 0L<<32|5);
    ASSERT_EQ(g.segs.at(2).data, "TTA");
}

TEST(targets_test, start_pos_point) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1+:7", target::role_t::START);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 2L<<32|0);
    ASSERT_EQ(a.w_lw, 4L<<32|0);
    // 2_0   to 4_0      TER+  to TGT1+ [START]
    ASSERT_EQ(g.arcs.at(0).v_lv, 2L<<32|0);
    ASSERT_EQ(g.arcs.at(0).w_lw, 4L<<32|0);
    // 4_e-b to 0_e      TGT1+ to SEG1+ [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 4L<<32|0);
    ASSERT_EQ(g.arcs.at(1).w_lw, 0L<<32|7);
    ASSERT_EQ(g.segs.at(2).data, "");
}


TEST(targets_test, start_neg_full) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1-", target::role_t::START);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 2L<<32|0);
    ASSERT_EQ(a.w_lw, 5L<<32|0);
    // 2_0   to 5_0      TER+  to TGT1- [START]
    ASSERT_EQ(g.arcs.at(0).v_lv, 2L<<32|0);
    ASSERT_EQ(g.arcs.at(0).w_lw, 5L<<32|0);
    // 5_e-b to 1_L-b    TGT1- to SEG1- [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 5L<<32|(10-0));
    ASSERT_EQ(g.arcs.at(1).w_lw, 1L<<32|(10-0));
    ASSERT_EQ(g.segs.at(2).data, SEG1.data);
}

TEST(targets_test, start_neg_part) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1-:2:5", target::role_t::START);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 2L<<32|0);
    ASSERT_EQ(a.w_lw, 5L<<32|0);
    // 2_0   to 5_0      TER+  to TGT1- [START]
    ASSERT_EQ(g.arcs.at(0).v_lv, 2L<<32|0);
    ASSERT_EQ(g.arcs.at(0).w_lw, 5L<<32|0);
    // 5_e-b to 1_L-b    TGT1- to SEG1- [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 5L<<32|(5-2));
    ASSERT_EQ(g.arcs.at(1).w_lw, 1L<<32|(10-2));
    ASSERT_EQ(g.segs.at(2).data, "TTA");
}

TEST(targets_test, start_neg_point) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1-:7", target::role_t::START);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 2L<<32|0);
    ASSERT_EQ(a.w_lw, 5L<<32|0);
    // 2_0   to 5_0      TER+  to TGT1- [START]
    ASSERT_EQ(g.arcs.at(0).v_lv, 2L<<32|0);
    ASSERT_EQ(g.arcs.at(0).w_lw, 5L<<32|0);
    // 5_e-b to 1_L-b    TGT1- to SEG1- [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 5L<<32|(7-7));
    ASSERT_EQ(g.arcs.at(1).w_lw, 1L<<32|(10-7));
    ASSERT_EQ(g.segs.at(2).data, "");
}


TEST(targets_test, end_pos_full) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1+", target::role_t::END);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 4L<<32|10);
    ASSERT_EQ(a.w_lw, 2L<<32|1);
    // 0_b   to 4_0      SEG1+ to TGT1+ [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 0L<<32|0);
    ASSERT_EQ(g.arcs.at(0).w_lw, 4L<<32|0);
    // 4_e-b to 2_1      TGT1+ to TER+  [END]
    ASSERT_EQ(g.arcs.at(1).v_lv, 4L<<32|10);
    ASSERT_EQ(g.arcs.at(1).w_lw, 2L<<32|1);
    ASSERT_EQ(g.segs.at(2).data, SEG1.data);
}

TEST(targets_test, end_pos_part) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1+:2:5", target::role_t::END);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 4L<<32|(5-2));
    ASSERT_EQ(a.w_lw, 2L<<32|1);
    // 0_b   to 4_0      SEG1+ to TGT1+ [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 0L<<32|2);
    ASSERT_EQ(g.arcs.at(0).w_lw, 4L<<32|0);
    // 4_e-b to 2_1      TGT1+ to TER+  [END]
    ASSERT_EQ(g.arcs.at(1).v_lv, 4L<<32|(5-2));
    ASSERT_EQ(g.arcs.at(1).w_lw, 2L<<32|1);
    ASSERT_EQ(g.segs.at(2).data, "TTA");
}

TEST(targets_test, end_pos_point) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1+:7", target::role_t::END);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 4L<<32|(7-7));
    ASSERT_EQ(a.w_lw, 2L<<32|1);
    // 0_b   to 4_0      SEG1+ to TGT1+ [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 0L<<32|7);
    ASSERT_EQ(g.arcs.at(0).w_lw, 4L<<32|0);
    // 4_e-b to 2_1      TGT1+ to TER+  [END]
    ASSERT_EQ(g.arcs.at(1).v_lv, 4L<<32|(7-7));
    ASSERT_EQ(g.arcs.at(1).w_lw, 2L<<32|1);
    ASSERT_EQ(g.segs.at(2).data, "");
}


TEST(targets_test, end_neg_full) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1-", target::role_t::END);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 5L<<32|10);
    ASSERT_EQ(a.w_lw, 2L<<32|1);
    // 1_L-e to 5_0      SEG1- to TGT1- [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 1L<<32|(10-10));
    ASSERT_EQ(g.arcs.at(0).w_lw, 5L<<32|0);
    // 5_e-b to 2_1      TGT1- to TER+  [END]
    ASSERT_EQ(g.arcs.at(1).v_lv, 5L<<32|(10-0));
    ASSERT_EQ(g.arcs.at(1).w_lw, 2L<<32|1);
    ASSERT_EQ(g.segs.at(2).data, SEG1.data);
}

TEST(targets_test, end_neg_part) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1-:2:5", target::role_t::END);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 5L<<32|(5-2));
    ASSERT_EQ(a.w_lw, 2L<<32|1);
    // 1_L-e to 5_0      SEG1- to TGT1- [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 1L<<32|(10-5));
    ASSERT_EQ(g.arcs.at(0).w_lw, 5L<<32|0);
    // 5_e-b to 2_1      TGT1- to TER+  [END]
    ASSERT_EQ(g.arcs.at(1).v_lv, 5L<<32|(5-2));
    ASSERT_EQ(g.arcs.at(1).w_lw, 2L<<32|1);
    ASSERT_EQ(g.segs.at(2).data, "TTA");
}

TEST(targets_test, end_neg_point) {
    graph g = make_graph();
    target t(g);
    t.set("SEG1-:7", target::role_t::END);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 5L<<32|(7-7));
    ASSERT_EQ(a.w_lw, 2L<<32|1);
    // 1_L-e to 5_0      SEG1- to TGT1- [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 1L<<32|(10-7));
    ASSERT_EQ(g.arcs.at(0).w_lw, 5L<<32|0);
    // 5_e-b to 2_1      TGT1- to TER+  [END]
    ASSERT_EQ(g.arcs.at(1).v_lv, 5L<<32|(7-7));
    ASSERT_EQ(g.arcs.at(1).w_lw, 2L<<32|1);
    ASSERT_EQ(g.segs.at(2).data, "");
}


TEST(targets_test, reinsert_test) {
    graph g = make_graph();
    target t(g);

    t.set("SEG1+:2:6", target::role_t::START);
    arc a = t.get_arc();
    ASSERT_EQ(a.v_lv, 2L<<32|0);
    ASSERT_EQ(a.w_lw, 4L<<32|0);
    // 2_0   to 4_0      TER+  to TGT1+ [START]
    ASSERT_EQ(g.arcs.at(0).v_lv, 2L<<32|0);
    ASSERT_EQ(g.arcs.at(0).w_lw, 4L<<32|0);
    // 4_e-b to 0_e      TGT1+ to SEG1+ [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 4L<<32|(6-2));
    ASSERT_EQ(g.arcs.at(1).w_lw, 0L<<32|6);
    ASSERT_EQ(g.segs.at(2).data, "TTAG");

    t.set("SEG1-:2:6", target::role_t::END);
    a = t.get_arc();
    ASSERT_EQ(a.v_lv, 5L<<32|(6-2));
    ASSERT_EQ(a.w_lw, 2L<<32|1);
    // 1_L-e to 5_0      SEG1- to TGT1- [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 1L<<32|(10-6));
    ASSERT_EQ(g.arcs.at(0).w_lw, 5L<<32|0);
    // 5_e-b to 2_1      TGT1- to TER+  [END]
    ASSERT_EQ(g.arcs.at(1).v_lv, 5L<<32|(6-2));
    ASSERT_EQ(g.arcs.at(1).w_lw, 2L<<32|1);
    ASSERT_EQ(g.segs.at(2).data, "TTAG");
}


TEST(targets_test, two_pos_tgts) {
    graph g = make_graph();
    target t1(g), t2(g);

    t1.set("SEG1+:1:3", target::role_t::START);
    t2.set("SEG1+:4:8", target::role_t::END);

    arc a1 = t1.get_arc();
    ASSERT_EQ(a1.v_lv, 2L<<32|0);
    ASSERT_EQ(a1.w_lw, 4L<<32|0);
    arc a2 = t2.get_arc();
    ASSERT_EQ(a2.v_lv, 6L<<32|4);
    ASSERT_EQ(a2.w_lw, 2L<<32|1);

    // 0_b   to 6_0      SEG1+ to TGT2+ [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 0L<<32|4);
    ASSERT_EQ(g.arcs.at(0).w_lw, 6L<<32|0);
    // 2_0   to 4_0      TER+  to TGT1+ [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 2L<<32|0);
    ASSERT_EQ(g.arcs.at(1).w_lw, 4L<<32|0);
    // 4_e-b to 0_e      TGT1+ to SEG1+ [START]
    ASSERT_EQ(g.arcs.at(2).v_lv, 4L<<32|2);
    ASSERT_EQ(g.arcs.at(2).w_lw, 0L<<32|3);
    // 6_e-b to 2_1      TGT2+ to TER+  [END]
    ASSERT_EQ(g.arcs.at(3).v_lv, 6L<<32|4);
    ASSERT_EQ(g.arcs.at(3).w_lw, 2L<<32|1);
}


// 0_b   to 6_0      SEG1+ to TGT1+ [END]
// 0_b   to 8_0      SEG1+ to TGT2+ [END]
// 1_L-e to 7_0      SEG1- to TGT1- [END]
// 1_L-e to 9_0      SEG1- to TGT2- [END]
// 2_b   to 6_0      SEG2+ to TGT1+ [END]
// 2_b   to 8_0      SEG2+ to TGT2+ [END]
// 3_L-e to 7_0      SEG2- to TGT1- [END]
// 3_L-e to 9_0      SEG2- to TGT2- [END]
// 4_0   to 6_0      TER+  to TGT1+ [START]
// 4_0   to 8_0      TER+  to TGT2+ [START]
// 4_0   to 7_0      TER+  to TGT1- [START]
// 4_0   to 9_0      TER+  to TGT2- [START]
// 6_e-b to 0_e      TGT1+ to SEG1+ [START]
// 6_e-b to 2_e      TGT1+ to SEG2+ [START]
// 6_e-b to 4_1      TGT1+ to TER+  [END]
// 7_e-b to 1_L-b    TGT1- to SEG1- [START]
// 7_e-b to 3_L-b    TGT1- to SEG2- [START]
// 7_e-b to 4_1      TGT1- to TER+  [END]
// 8_e-b to 0_e      TGT2+ to SEG1+ [START]
// 8_e-b to 2_e      TGT2+ to SEG2+ [START]
// 8_e-b to 4_1      TGT2+ to TER+  [END]
// 9_e-b to 1_L-b    TGT2- to SEG1- [START]
// 9_e-b to 3_L-b    TGT2- to SEG2- [START]
// 9_e-b to 4_1      TGT2- to TER+  [END]

TEST(targets_test, two_segs_two_tgts) {

    graph g;
    g.segs.reserve(2 + 3 /* two targets and a terminator */);
    g.arcs.reserve(2 * 4 /* four for each target */);
    g.add_seg({ 11, "SEG1", "GCTATGACAAT" });
    g.add_seg({ 9,  "SEG2", "TTGTATAGT" });

    target t1(g), t2(g);
    t1.set("SEG1-:4:9", target::role_t::START);
    t2.set("SEG2+:3:8", target::role_t::END);

    arc a1 = t1.get_arc();
    ASSERT_EQ(a1.v_lv, 4L<<32);
    ASSERT_EQ(a1.w_lw, 7L<<32);
    arc a2 = t2.get_arc();
    ASSERT_EQ(a2.v_lv, 8L<<32|(8-3));
    ASSERT_EQ(a2.w_lw, 4L<<32|1);

    // 2_b   to 8_0      SEG2+ to TGT2+ [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 2L<<32|3);
    ASSERT_EQ(g.arcs.at(0).w_lw, 8L<<32|0);
    // 4_0   to 7_0      TER+  to TGT1- [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 4L<<32|0);
    ASSERT_EQ(g.arcs.at(1).w_lw, 7L<<32|0);
    // 7_e-b to 1_L-b    TGT1- to SEG1- [START]
    ASSERT_EQ(g.arcs.at(2).v_lv, 7L<<32|(9-4));
    ASSERT_EQ(g.arcs.at(2).w_lw, 1L<<32|(11-4));
    // 8_e-b to 4_1      TGT2+ to TER+  [END]
    ASSERT_EQ(g.arcs.at(3).v_lv, 8L<<32|(8-3));
    ASSERT_EQ(g.arcs.at(3).w_lw, 4L<<32|1);

        // now flip them

    t1.set("SEG1-:4:9", target::role_t::END);
    t2.set("SEG2+:3:8", target::role_t::START);

    a1 = t1.get_arc();
    ASSERT_EQ(a1.v_lv, 7L<<32|(9-4));
    ASSERT_EQ(a1.w_lw, 4L<<32|1);
    a2 = t2.get_arc();
    ASSERT_EQ(a2.v_lv, 4L<<32|0);
    ASSERT_EQ(a2.w_lw, 8L<<32|0);

    // 1_L-e to 7_0      SEG1- to TGT1- [END]
    ASSERT_EQ(g.arcs.at(0).v_lv, 1L<<32|(11-9));
    ASSERT_EQ(g.arcs.at(0).w_lw, 7L<<32|0);
    // 4_0   to 8_0      TER+  to TGT2+ [START]
    ASSERT_EQ(g.arcs.at(1).v_lv, 4L<<32|0);
    ASSERT_EQ(g.arcs.at(1).w_lw, 8L<<32|0);
    // 7_e-b to 4_1      TGT1- to TER+  [END]
    ASSERT_EQ(g.arcs.at(2).v_lv, 7L<<32|(9-4));
    ASSERT_EQ(g.arcs.at(2).w_lw, 4L<<32|1);
    // 8_e-b to 2_e      TGT2+ to SEG2+ [START]
    ASSERT_EQ(g.arcs.at(3).v_lv, 8L<<32|(8-3));
    ASSERT_EQ(g.arcs.at(3).w_lw, 2L<<32|8);
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
