/* paths-test.cpp
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
#include "paths.h"
#include "targets.h"
#include "utils.h"

using namespace gfa;

namespace {

static seg SEG1 = { 4, "s1", "ACGT" };
static seg SEG2 = { 9, "s2", "TAGCATACG" }; // rc: CGTATGCTA
static seg SEG3 = { 5, "s3", "CATTA" };
static seg SEG4 = { 8, "s4", "CTATAATT" };

static graph make_graph() {
    graph g;
    g.segs.reserve(5);
    g.add_seg(SEG1);
    g.add_seg(SEG2);
    g.add_seg(SEG3);
    g.add_seg(SEG4);
    g.arcs.reserve(4*8+1);
    g.add_edge("s1+", 1, 4, "s2-", 5, 9); // s1+ .[--) s2- (---].....     and s2+ .....[---) s1- (--].
    g.add_edge("s2-", 0, 0, "s3+", 0, 0); // s2- .........) s3+ (.....    and s3- .....) s2+ (.........
    g.add_edge("s2-", 0, 3, "s4+", 0, 3); // s2- ......[..) s4+ (..]..... and s4- .....[..) s2+ (..]......
    g.add_edge("s3+", 4, 5, "s1+", 0, 1); // s3+ ....[) s1+ (]...         and s1- ...[) s3- (]...
    return g;
}

static const arc* add_start(graph& g, std::string ref) {
    target t(g);
    t.set(ref, target::role_t::START);
    return &*g.arcs_from_v_lv(t.get_arc().w_lw).first;
}

TEST(paths_test, empty_path) {
    graph g = make_graph();
    ASSERT_EQ(g.get_seg("s1").len, 4);
    paths p = paths(g);
    ASSERT_EQ(p.path_arcs.size(), 1);
    ASSERT_EQ(p.path_arcs.at(0).pre_ix, 0);
}

TEST(paths_test, path_1) {
    graph g = make_graph();
    const arc  *a = add_start(g, "s1:0+");
    paths p(g);
    std::size_t i = p.extend(0, a);
    ASSERT_EQ(p.path_arcs.size(), 2);
    ASSERT_EQ(p.path_arcs.at(i).pre_ix, 0);
    ASSERT_EQ(p.path_arcs.at(i).p_arc, a);
}

TEST(paths_test, write_empty) {
    graph g = make_graph();
    const arc  *a = add_start(g, "s1:0+");
    paths p(g);
    std::size_t i = p.extend(0, a);
    ASSERT_EQ(p.path_arcs.size(), 2);

    const path_arc& pa = p.path_arcs.at(i);
    ASSERT_EQ(pa.pre_ix, 0);
    ASSERT_EQ(pa.p_arc, a);
    ASSERT_EQ(p.ride_len(pa), 0);
    ASSERT_EQ(p.length(pa), 0);
    ASSERT_EQ(p.route(pa), "");
    ASSERT_EQ(p.sequence(pa), "");
}

TEST(paths_test, write_1) {
    graph g = make_graph();
    const arc* a = add_start(g, "s3:2+"); // s3+ CA|TTA
    paths p(g);
    std::size_t i = p.extend(0, a);
    ASSERT_EQ(i, i);

    std::vector<arc>::const_iterator arc_it = g.arcs_from_v_lv(graph::v_lv(graph::seg_vtx_p(2), 2)).first;
    ASSERT_EQ(arc_it->v_lv, 2L<<33|4);  // s3+ CATT|A
    ASSERT_EQ(arc_it->w_lw, 0);         // s1+ |ACGT

    i = p.extend(i, &*arc_it);

    const path_arc& pa = p.path_arcs.at(i);
    ASSERT_EQ(pa.pre_ix, 1);
    ASSERT_EQ(pa.p_arc, &*arc_it);
    ASSERT_EQ(p.ride_len(pa), 2);
    ASSERT_EQ(p.length(pa), 2);
    ASSERT_EQ(p.route(pa), "s3:2:4+");
    ASSERT_EQ(p.sequence(pa), "TT");
}

TEST(paths_test, write_2) {
    graph g = make_graph();
    const arc* a = add_start(g, "s3:1+"); // s3+ C|ATTA
    paths p = paths(g);
    std::size_t i = p.extend(0, a);

    std::vector<arc>::const_iterator arc_it = g.arcs_from_v_lv(graph::v_lv(graph::seg_vtx_p(g.get_seg_ix("s3")),1)).first;
    ASSERT_EQ(arc_it->v_lv, 2L<<33|4);   // s3+:4 C|ATT|A
    ASSERT_EQ(arc_it->w_lw, 0);          // s1+:0 |A|CGT

    // add arc onto s1+ [A]CGT, so we have C|ATT|A
    i = p.extend(i, &*arc_it);

    const path_arc* pa = &p.path_arcs.at(i);

    ASSERT_EQ(pa->pre_ix, 1);
    ASSERT_EQ(pa->p_arc, &*arc_it);
    ASSERT_EQ(p.ride_len(*pa), 3);
    ASSERT_EQ(p.length(*pa), 3);
    ASSERT_EQ(p.route(*pa), "s3:1:4+");
    ASSERT_EQ(p.sequence(*pa), "ATT");

    // find first arc away from 1+, is return arc to where we came from
    arc_it = g.arcs_from_v_lv(arc_it->w_lw).first;
    ASSERT_EQ(arc_it->v_lv, 0);          // s1+:0
    ASSERT_EQ(arc_it->w_lw, 2L<<33|4);   // s3+:4

    // but next is A[CTG] overlap to s2- [CGTA]TGCTA
    arc_it = arc_it + 1;
    ASSERT_EQ(arc_it->v_lv, 1);          // s1+:1
    ASSERT_EQ(arc_it->w_lw, 3L<<32|0);   // s2-:0

    // add that so we have C|ATT|A| with CGTATGCTA coming
    i = p.extend(i, &*arc_it);
    pa = &p.path_arcs.at(i);
    ASSERT_EQ(p.ride_len(*pa), 1);
    ASSERT_EQ(p.length(*pa), 4);
    ASSERT_EQ(p.route(*pa), "s3:1:4+ s1:0:1+");
    ASSERT_EQ(p.sequence(*pa), "ATTA");

    // find first arc away from 2- is return arc to where we came from
    arc_it = g.arcs_from_v_lv(arc_it->w_lw).first;
    // next is to s1+ end
    arc_it = arc_it + 1;
    // next is to s4+ start
    arc_it = arc_it + 1;
    ASSERT_EQ(arc_it->v_lv, 3L<<32|6);   // s2-:6
    ASSERT_EQ(arc_it->w_lw, 6L<<32|0);   // s4+:0

    // add arc to s1+ A[CGT] and we have C|ATT|A|CGTATG|CTA
    i = p.extend(i, &*arc_it);
    pa = &p.path_arcs.at(i);
    ASSERT_EQ(p.ride_len(*pa), 6);
    ASSERT_EQ(p.length(*pa), 10);
    ASSERT_EQ(p.route(*pa), "s3:1:4+ s1:0:1+ s2:3:9-");
    ASSERT_EQ(p.sequence(*pa), "ATTACGTATG");

    // find first arc away from 4+, is return arc to where we came from
    arc_it = g.arcs_from_v_lv(arc_it->w_lw).first;
    // and next is to s2- end
    arc_it = arc_it + 1;
    ASSERT_EQ(arc_it->v_lv, 6L<<32|3);   // s4+:3   // adds CTA
    ASSERT_EQ(arc_it->w_lw, 3L<<32|9);   // s2-:9

    // we take it and have |ATT|A|CGTATG|CTA|
    i = p.extend(i, &*arc_it);
    pa = &p.path_arcs.at(i);
    ASSERT_EQ(p.ride_len(*pa), 3);
    ASSERT_EQ(p.length(*pa), 13);
    ASSERT_EQ(p.route(*pa), "s3:1:4+ s1:0:1+ s2:3:9- s4:0:3+");
    ASSERT_EQ(p.sequence(*pa), "ATTACGTATGCTA");
}

} // namespace
  // vim: sts=4:sw=4:ai:si:et
