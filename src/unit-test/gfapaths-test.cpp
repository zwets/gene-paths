/* gfapaths-test.cpp
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
#include "gfapaths.h"
#include "utils.h"

using namespace gfa;

namespace {

static seg SEG1 = { 4, "s1", "ACGT" };
static seg SEG2 = { 9, "s2", "TAGCATACG" };
static seg SEG3 = { 5, "s3", "CATTA" };
static seg SEG4 = { 8, "s4", "CTATAATT" };

static graph make_graph() {
    graph gfa;
    gfa.add_seg(SEG1);
    gfa.add_seg(SEG2);
    gfa.add_seg(SEG3);
    gfa.add_seg(SEG4);
    gfa.add_edge("s1+", 1, 4, "s2-", 5, 9); // s1+ .[--) s2- (---].....     and s2+ .....[---) s1- (--].
    gfa.add_edge("s2-", 0, 0, "s3+", 0, 0); // s2- .........) s3+ (.....    and s3- .....) s2+ (.........
    gfa.add_edge("s2-", 0, 3, "s4+", 0, 3); // s2- ......[..) s4+ (..]..... and s4- .....[..) s2+ (..]......
    gfa.add_edge("s3+", 4, 5, "s1+", 0, 1); // s3+ ....[) s1+ (]...         and s1- ...[) s3- (]...
    return gfa;
}

TEST(gfapaths_test, empty_path) {
    graph g = make_graph();
    paths p = paths(g, 1);
    std::size_t p_ix = p.start_path(graph::v_lv(0, 0));
    ASSERT_EQ(p_ix, 1);
    ASSERT_EQ(p.path_starts.size(), 1);
    ASSERT_EQ(p.path_arcs.size(), 2);
    ASSERT_EQ(p.path_arcs.at(1).pre_ix, 0);
    ASSERT_EQ(p.path_arcs.at(1).p_arc, &*p.path_starts.cbegin());
}

TEST(gfapaths_test, path_1) {
    graph g = make_graph();
    paths p = paths(g, 1);
    std::size_t p_ix = p.start_path(graph::v_lv(graph::seg_vtx_p(0), 0));
    p.extend(p_ix, g.arcs.cbegin());
    ASSERT_EQ(p.path_starts.size(), 1);
    ASSERT_EQ(p.path_arcs.size(), 3);
    ASSERT_EQ(p.path_arcs.at(2).pre_ix, p_ix);
    ASSERT_EQ(p.path_arcs.at(2).p_arc, &*g.arcs.cbegin());
}

static std::string get_seq(const paths& ps, std::size_t pix)
{
    std::stringstream ss;
    ps.write_path_seq(ss, pix);
    return ss.str();
}
    
TEST(gfapaths_test, write_empty) {
    graph g = make_graph();
    paths p = paths(g, 1);
    std::size_t p_ix = p.start_path(graph::v_lv(graph::seg_vtx_p(0), 0));
    ASSERT_EQ(p.ride_len(p.path_arcs.at(p_ix)), 0);
    ASSERT_EQ(p.length(p.path_arcs.at(p_ix)), 0);
    ASSERT_EQ(get_seq(p, p_ix), "");
}

TEST(gfapaths_test, write_1) {
    graph g = make_graph();
    paths p = paths(g, 1);
    std::size_t p_ix = p.start_path(graph::v_lv(graph::seg_vtx_p(2), 2)); // s3+ CA|TTA

    std::vector<arc>::const_iterator arc_it = g.arcs_from_v_lv(graph::v_lv(graph::seg_vtx_p(2), 2)).first;
    ASSERT_EQ(arc_it->v_lv, 2L<<33|4);  // s3+ CATT|A
    ASSERT_EQ(arc_it->w_lw, 0);         // s1+ |ACGT

    p.extend(p_ix, arc_it);
    ASSERT_EQ(p.path_starts.size(), 1);
    ASSERT_EQ(p.path_arcs.size(), 3);
    ASSERT_EQ(p.path_arcs.at(2).pre_ix, p_ix);
    ASSERT_EQ(p.path_arcs.at(2).p_arc, &*arc_it);
    ASSERT_EQ(p.ride_len(p.path_arcs.at(p_ix+1)), 2);
    ASSERT_EQ(p.length(p.path_arcs.at(p_ix+1)), 2);
    ASSERT_EQ(get_seq(p, p_ix+1), "TT");
}

TEST(gfapaths_test, write_2) {
    graph g = make_graph();
    paths p = paths(g, 1);
    std::size_t p_ix = p.start_path(graph::v_lv(graph::seg_vtx_p(g.get_seg_ix("s3")), 1));  // start at s3+: C|ATTA

    std::vector<arc>::const_iterator arc_it = g.arcs_from_v_lv(graph::v_lv(graph::seg_vtx_p(g.get_seg_ix("s3")),1)).first;
    ASSERT_EQ(arc_it->v_lv, 2L<<33|4);   // s3+:4 C|ATT|A
    ASSERT_EQ(arc_it->w_lw, 0);          // s1+:0 |A|CGT

    // add arc onto s1+ [A]CGT, so we have C|ATT|A
    p.extend(p_ix++, arc_it);
    ASSERT_EQ(p.path_starts.size(), 1);
    ASSERT_EQ(p.path_arcs.size(), p_ix+1);
    ASSERT_EQ(p.path_arcs.at(p_ix).pre_ix, p_ix-1);
    ASSERT_EQ(p.path_arcs.at(p_ix).p_arc, &*arc_it);
    ASSERT_EQ(p.ride_len(p.path_arcs.at(p_ix)), 3);
    ASSERT_EQ(p.length(p.path_arcs.at(p_ix)), 3);
    ASSERT_EQ(get_seq(p, p_ix), "ATT");

    // find first arc away from 1+, is return arc to where we came from
    arc_it = g.arcs_from_v_lv(arc_it->w_lw).first;
    ASSERT_EQ(arc_it->v_lv, 0);          // s1+:0
    ASSERT_EQ(arc_it->w_lw, 2L<<33|4);   // s3+:4

    // but next is A[CTG] overlap to s2- [CGTA]TGCTA
    arc_it = arc_it + 1;
    ASSERT_EQ(arc_it->v_lv, 1);          // s1+:1
    ASSERT_EQ(arc_it->w_lw, 3L<<32|0);   // s2-:0

    // add that so we have C|ATT|A| with CGTATGCTA coming
    p.extend(p_ix++, arc_it);
    ASSERT_EQ(p.ride_len(p.path_arcs.at(p_ix)), 1);
    ASSERT_EQ(p.length(p.path_arcs.at(p_ix)), 4);
    ASSERT_EQ(get_seq(p, p_ix), "ATTA");

    // find first arc away from 2- is return arc to where we came from
    arc_it = g.arcs_from_v_lv(arc_it->w_lw).first;
    // next is to s1+ end
    arc_it = arc_it + 1;
    // next is to s4+ start
    arc_it = arc_it + 1;
    ASSERT_EQ(arc_it->v_lv, 3L<<32|6);   // s2-:6
    ASSERT_EQ(arc_it->w_lw, 6L<<32|0);   // s4+:0

    // add arc to s1+ A[CGT] and we have C|ATT|A|CGTATG|CTA
    p.extend(p_ix++, arc_it);
    ASSERT_EQ(p.ride_len(p.path_arcs.at(p_ix)), 6);
    ASSERT_EQ(p.length(p.path_arcs.at(p_ix)), 10);
    ASSERT_EQ(get_seq(p, p_ix), "ATTACGTATG");

    // find first arc away from 4+, is return arc to where we came from
    arc_it = g.arcs_from_v_lv(arc_it->w_lw).first;
    // and next is to s2- end
    arc_it = arc_it + 1;
    ASSERT_EQ(arc_it->v_lv, 6L<<32|3);   // s4+:3   // adds CTA
    ASSERT_EQ(arc_it->w_lw, 3L<<32|9);   // s2-:9

    // we take it and have |ATT|A|CGTATG|CTA|
    p.extend(p_ix++, arc_it);
    ASSERT_EQ(p.ride_len(p.path_arcs.at(p_ix)), 3);
    ASSERT_EQ(p.length(p.path_arcs.at(p_ix)), 13);
    ASSERT_EQ(get_seq(p, p_ix), "ATTACGTATGCTA");
}

} // namespace
  // vim: sts=4:sw=4:ai:si:et
