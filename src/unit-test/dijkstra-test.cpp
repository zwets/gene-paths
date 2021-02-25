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
#include "utils.h"

using namespace gfa;

namespace {

static seg SEG1 = { 3, "s1", "CAT" };
static seg SEG2 = { 4, "s2", "TAGT" };
static std::string FROM = "s1:0:1+";
static std::string TO = "s2:2:3+";

static graph simple_graph() {
    graph g;
    g.segs.reserve(2+3);  // 2 segs plus 2 targets and the TERminal
    g.add_seg(SEG1);  // seg_ix 0
    g.add_seg(SEG2);  // seg_ix 1
    // TER            // seg_ix 2
    // TGT1           // seg_ix 3
    // TGT2           // seg_ix 4
    g.arcs.reserve(1*8 + 4);  // 1 edge, 4 arcs for targets
    g.add_edge("s1+", 2, 3, "s2+", 0, 1);
    return g;
}

typedef std::pair<const arc*, const arc*> parc_pair;
static parc_pair add_targets(graph& g) {

    target t0(g);
    t0.set(FROM, target::role_t::START);
    target t1(g);
    t1.set(TO, target::role_t::END);

    return std::make_pair(t0.p_arc(), t1.p_arc());
}

TEST(dijkstra_test, make_graph) {
    graph g = simple_graph();
    ASSERT_EQ(g.segs.size(), 2);
    ASSERT_EQ(g.arcs.size(), 8);
}

TEST(dijkstra_test, add_targets) {
    graph g = simple_graph();
    parc_pair t = add_targets(g);
    ASSERT_EQ(g.segs.size(), 2+3);
    ASSERT_EQ(g.arcs.size(), 1*8+4);
    ASSERT_EQ(t.first->v_lv, ((2L<<1)<<32)|0);  // TER+:0
    ASSERT_EQ(t.first->w_lw, ((3L<<1)<<32)|0);  // to TGT1+:0
    ASSERT_EQ(t.second->v_lv,((4L<<1)<<32)|1);  // TGT2+:$
    ASSERT_EQ(t.second->w_lw,((2L<<1)<<32)|1);  // to TER+:1
}

TEST(dijkstra_test, dijkstra_construct) {
    graph g = simple_graph();
    add_targets(g);

    dijkstra dk(g);
    ASSERT_EQ(dk.ps.path_arcs.size(), 1);       // just the null path
    ASSERT_EQ(dk.ds.size(), 12);                // 8 for the edge, 4 for targets
    ASSERT_EQ(dk.vs.size(), 0);                 // nothing visitable
    ASSERT_FALSE(dk.found_pix);                 // nothing found
    ASSERT_EQ(dk.found_len, 0);
}

TEST(dijkstra_test, dijkstra_restart) {
    graph g = simple_graph();
    parc_pair t = add_targets(g);

    dijkstra dk(g);
    dk.restart(t.first);                        // restart with given start arc
    ASSERT_EQ(dk.ps.path_arcs.size(), 2);       // null and the start arc
    ASSERT_EQ(dk.ds.size(), 12);                // same as before
    ASSERT_EQ(dk.vs.size(), 1);                 // start node is single visitable
    ASSERT_FALSE(dk.found_pix);
    ASSERT_EQ(dk.found_len, 0);
}

TEST(dijkstra_test, pop_visit) {
    graph g = simple_graph();
    parc_pair t = add_targets(g);

    dijkstra dk(g);
    dk.restart(t.first);
    ASSERT_EQ(dk.vs.size(), 1);                 // start node is single visitable

    dijkstra::dnode& d = dk.pop_visit();
    ASSERT_EQ(dk.vs.size(), 0);                 // visitable was removed
    ASSERT_EQ(d.len, 0);                        // path to initial visit node is 0 length
    ASSERT_EQ(d.p_ref, 1);                      // its path index is 1
    ASSERT_FALSE(d.is_visited());               // and it is not yet marked visited
}

TEST(dijkstra_test, did_nay_run) {
    graph g = simple_graph();
    add_targets(g);
    dijkstra dk(g);

    ASSERT_EQ(dk.found_len, 0);
    ASSERT_EQ(dk.length(), 0);                  // we don't want it to crash when
    ASSERT_EQ(dk.route(), "");                  // no paths have yet been searched
    ASSERT_EQ(dk.sequence(), "");
}

TEST(dijkstra_test, shortest_path) {
    graph g = simple_graph();
    parc_pair t = add_targets(g);
    dijkstra dk(g);

    ASSERT_TRUE(dk.shortest_path(t.first, t.second));
    ASSERT_TRUE(dk.found_pix);
    ASSERT_EQ(dk.found_len, 5);
    ASSERT_EQ(dk.length(), 5);
    ASSERT_EQ(dk.route(), "s1:0:1+ s1:1:2+ s2:0:2+ s2:2:3+");
    ASSERT_EQ(dk.sequence(), "CATAG");
}

TEST(dijkstra_test, shortest_paths) {
    graph g = simple_graph();
    parc_pair t = add_targets(g);
    dijkstra dk(g);

    dk.shortest_paths(t.first);
    ASSERT_FALSE(dk.found_pix);                 // always unset when searching all paths
    ASSERT_EQ(dk.found_len, 0);
    ASSERT_EQ(dk.ps.path_arcs.size(), 4+4);     // half of the 8 from the edge, plus the target ones
    size_t last = dk.ps.path_arcs.size() - 1;
    ASSERT_EQ(dk.length(last), 5);
    ASSERT_EQ(dk.route(last), "s1:0:1+ s1:1:2+ s2:0:2+ s2:2:3+");
    ASSERT_EQ(dk.sequence(last), "CATAG");
}

TEST(dijkstra_test, furthest_path_from) {
    graph g = simple_graph();
    parc_pair t = add_targets(g);
    dijkstra dk(g);

    dk.furthest_path(t.first);
    ASSERT_TRUE(dk.found_pix);
    ASSERT_EQ(dk.found_len, 5);
    ASSERT_EQ(dk.length(), 5);
    ASSERT_EQ(dk.route(), "s1:0:1+ s1:1:2+ s2:0:2+ s2:2:3+");
    ASSERT_EQ(dk.sequence(), "CATAG");
}

/*
TEST(dijkstra_test, furthest_path) {
    gene_paths::set_verbose(true);
    graph g = simple_graph();
    add_targets(g);
    dijkstra dk(g);

    dk.furthest_path();
    ASSERT_TRUE(dk.found_pix);
    ASSERT_EQ(dk.found_len, 5);
    ASSERT_EQ(dk.length(), 5);
    ASSERT_EQ(dk.route(), "s1:0:1+ s1:1:2+ s2:0:2+ s2:2:3+");
    ASSERT_EQ(dk.sequence(), "CATAG");
}
*/

} // namespace
  // vim: sts=4:sw=4:ai:si:et
