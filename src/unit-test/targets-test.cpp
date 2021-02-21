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

using namespace gfa;

namespace {

TEST(targets_test, no_pos) {
    gfa::target t = gfa::target::parse("C1+");
    ASSERT_EQ(t.ctg, "C1");
    ASSERT_FALSE(t.neg);
    ASSERT_EQ(t.beg, std::uint64_t(-1));
    ASSERT_EQ(t.end, std::uint64_t(-1));
}

TEST(targets_test, neg_ctg) {
    gfa::target t = gfa::target::parse("C2-");
    ASSERT_EQ(t.ctg, "C2");
    ASSERT_TRUE(t.neg);
    ASSERT_EQ(t.beg, std::uint64_t(-1));
    ASSERT_EQ(t.end, std::uint64_t(-1));
}

TEST(targets_test, beg_only) {
    gfa::target t = gfa::target::parse("C3-:0");
    ASSERT_EQ(t.ctg, "C3");
    ASSERT_TRUE(t.neg);
    ASSERT_EQ(t.beg, 0);
    ASSERT_EQ(t.end, std::uint64_t(-1));
}

TEST(targets_test, beg_and_end) {
    gfa::target t = gfa::target::parse("C4-:0:2");
    ASSERT_EQ(t.ctg, "C4");
    ASSERT_TRUE(t.neg);
    ASSERT_EQ(t.beg, 0);
    ASSERT_EQ(t.end, 2);
}

static seg SEG1 = { 10, "SEG1", "CATTAGTACT" };
static seg REV1 = { 10, "REV1", "AGTACTAATG" };

static graph make_graph()
{
    graph g;
    g.segs.reserve(3);
    g.arcs.reserve(4);
    g.add_seg(SEG1);
    g.add_seg(REV1);
    return g;
}

TEST(targets_test, in_graph) {
    graph g = make_graph();
    target t = target::parse("SEG1+");
    t.add_seg_to_graph(g, "WHOLE");
    ASSERT_EQ(t.name, "WHOLE");
    ASSERT_EQ(t.ctg, "SEG1");
    ASSERT_EQ(t.beg, 0);
    ASSERT_EQ(t.end, 10);
    const seg& s = g.get_seg("WHOLE");
    ASSERT_EQ(s.len, 10);
    ASSERT_EQ(s.data, SEG1.data);
    ASSERT_EQ(s.name, "WHOLE");
}

TEST(targets_test, in_graph_part) {
    graph g = make_graph();
    target t = target::parse("SEG1+:2:5");
    t.add_seg_to_graph(g, "PART");
    ASSERT_EQ(t.name, "PART");
    ASSERT_EQ(t.ctg, "SEG1");
    ASSERT_EQ(t.beg, 2);
    ASSERT_EQ(t.end, 5);
    const seg& s = g.get_seg("PART");
    ASSERT_EQ(s.len, 3);
    ASSERT_EQ(s.data, "TTA");
    ASSERT_EQ(s.name, "PART");
}

TEST(targets_test, in_graph_0) {
    graph g = make_graph();
    target t = target::parse("SEG1+:7");
    t.add_seg_to_graph(g, "ZERO");
    ASSERT_EQ(t.name, "ZERO");
    ASSERT_EQ(t.ctg, "SEG1");
    ASSERT_EQ(t.beg, 7);
    ASSERT_EQ(t.end, 7);
    const seg& s = g.get_seg("ZERO");
    ASSERT_EQ(s.len, 0);
    ASSERT_EQ(s.data, "");
}

TEST(targets_test, arc_into) {
    graph g = make_graph();
    target t = target::parse("SEG1+");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, true);
    ASSERT_EQ(p->v_lv, graph::v_lv(0,0));
    ASSERT_EQ(p->w_lw, graph::v_lv(4,0));
}

TEST(targets_test, arc_out_of) {
    graph g = make_graph();
    target t = target::parse("SEG1+");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, false);
    ASSERT_EQ(p->v_lv, graph::v_lv(4,10));
    ASSERT_EQ(p->w_lw, graph::v_lv(0,10));
}

TEST(targets_test, arc_neg_into) {
    graph g = make_graph();
    target t = target::parse("SEG1-");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, true);
    ASSERT_EQ(p->v_lv, graph::v_lv(1,0));
    ASSERT_EQ(p->w_lw, graph::v_lv(5,0));
}

TEST(targets_test, arc_part_into) {
    graph g = make_graph();
    target t = target::parse("SEG1+:2:5");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, true);
    ASSERT_EQ(p->v_lv, graph::v_lv(0,2));
    ASSERT_EQ(p->w_lw, graph::v_lv(4,0));
}

TEST(targets_test, arc_part_out_of) {
    graph g = make_graph();
    target t = target::parse("SEG1+:2:5");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, false);
    ASSERT_EQ(p->v_lv, graph::v_lv(4,3));
    ASSERT_EQ(p->w_lw, graph::v_lv(0,5));
}

TEST(targets_test, arc_part_neg_into) {
    graph g = make_graph();
    target t = target::parse("SEG1-:2:5");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, true);
    ASSERT_EQ(p->v_lv, graph::v_lv(1,5));
    ASSERT_EQ(p->w_lw, graph::v_lv(5,0));
}

TEST(targets_test, arc_part_neg_out_of) {
    graph g = make_graph();
    target t = target::parse("SEG1-:2:5");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, false);
    ASSERT_EQ(p->v_lv, graph::v_lv(5,3));
    ASSERT_EQ(p->w_lw, graph::v_lv(1,8));
}

TEST(targets_test, arc_1_into) {
    graph g = make_graph();
    target t = target::parse("SEG1+:7");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, true);
    ASSERT_EQ(p->v_lv, graph::v_lv(0,7));
    ASSERT_EQ(p->w_lw, graph::v_lv(4,0));
}

TEST(targets_test, arc_1_out_of) {
    graph g = make_graph();
    target t = target::parse("SEG1+:7");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, false);
    ASSERT_EQ(p->v_lv, graph::v_lv(4,0));
    ASSERT_EQ(p->w_lw, graph::v_lv(0,7));
}

TEST(targets_test, arc_1_neg_into) {
    graph g = make_graph();
    target t = target::parse("SEG1-:7");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, true);
    ASSERT_EQ(p->v_lv, graph::v_lv(1,3));
    ASSERT_EQ(p->w_lw, graph::v_lv(5,0));
}

TEST(targets_test, arc_1_neg_out_of) {
    graph g = make_graph();
    target t = target::parse("SEG1-:7");
    t.add_seg_to_graph(g, "TGT");
    arc *p = t.add_arc_to_graph(g, false);
    ASSERT_EQ(p->v_lv, graph::v_lv(5,0));
    ASSERT_EQ(p->w_lw, graph::v_lv(1,3));
}

} // namespace
  // vim: sts=4:sw=4:ai:si:et
