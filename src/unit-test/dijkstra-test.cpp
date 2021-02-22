/* dijkstra-test.cpp
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
#include <string>
#include "dijkstra.h"
#include "targets.h"

using namespace gfa;

namespace {

static seg SEG1 = { 3, "s1", "CAT" };
static seg SEG2 = { 4, "s2", "TAGT" };
static std::string FROM = "s1+:0:1";
static std::string TO = "s2+:2:3";

static graph simple_graph() {
    graph g;
    g.segs.reserve(4);
    g.add_seg(SEG1);
    g.add_seg(SEG2);
    g.arcs.reserve(1*8 + 2*2);
    g.add_edge("s1+", 2, 3, "s2+", 0, 1);
    return g;
}

static std::pair<arc*, arc*> add_targets(graph& g) {

    target t0 = target::parse(FROM);
    target t1 = target::parse(TO);

    t0.add_seg_to_graph(g, FROM);
    t0.add_arc_to_graph(g, false); // from tgt into ctg
    arc* p_start = t0.add_arc_to_graph(g, true);  // from ctg to tgt

    t1.add_seg_to_graph(g, TO);
    t1.add_arc_to_graph(g, true); // from ctg into tgt
    arc* p_end = t1.add_arc_to_graph(g, false);  // from tgt to ctg

    return std::make_pair(p_start, p_end);
}

TEST(dijkstra_test, make_graph) {
    graph g = simple_graph();
    ASSERT_EQ(g.segs.size(), 2);
    ASSERT_EQ(g.arcs.size(), 8);
}

TEST(dijkstra_test, add_targets) {
    graph g = simple_graph();
    std::pair<arc*, arc*> t = add_targets(g);
    ASSERT_EQ(g.segs.size(), 4);
    ASSERT_EQ(g.arcs.size(), 12);
    ASSERT_EQ(t.first->v_lv, 0L);
    ASSERT_EQ(t.first->w_lw, 2L<<33);
    ASSERT_EQ(t.second->v_lv, 3L<<33|1);
    ASSERT_EQ(t.second->w_lw, 1L<<33|3);
}

TEST(dijkstra_test, dijkstra_con) {
    graph g = simple_graph();
    add_targets(g);

    dijkstra dk(g);
    ASSERT_EQ(dk.ps.path_arcs.size(), 1);
    ASSERT_EQ(dk.ds.size(), 12);
    ASSERT_EQ(dk.vs.size(), 0);
    ASSERT_FALSE(dk.found);
}

TEST(dijkstra_test, dijkstra_restart) {
    graph g = simple_graph();
    std::pair<arc*, arc*> t = add_targets(g);

    dijkstra dk(g);
    dk.restart(t.first);
    ASSERT_EQ(dk.ps.path_arcs.size(), 2);
    ASSERT_EQ(dk.ds.size(), 12);
    ASSERT_EQ(dk.vs.size(), 1);
    ASSERT_FALSE(dk.found);
}

TEST(dijkstra_test, pop_visit) {
    graph g = simple_graph();
    std::pair<arc*, arc*> t = add_targets(g);
    dijkstra dk(g);
    dk.restart(t.first);

    dijkstra::dnode& d = dk.pop_visit();
    ASSERT_EQ(dk.vs.size(), 0);
    ASSERT_EQ(d.len, 0);
    ASSERT_EQ(d.p_ref, 1);
    ASSERT_FALSE(d.is_visited());
}

TEST(dijkstra_test, all_paths) {
    graph g = simple_graph();
    std::pair<arc*, arc*> t = add_targets(g);
    dijkstra dk(g);
    dk.all_paths(t.first);
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
