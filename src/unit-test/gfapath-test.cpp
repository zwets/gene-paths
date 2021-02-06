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
#include "gfapath.h"
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

TEST(gfapath_tes, empty_path) {
    path p(1, 100);
    ASSERT_EQ(p.len, 0);
    ASSERT_EQ(p.start, std::uint64_t(1)<<32|100);
    ASSERT_EQ(p.now_at, p.start);
    ASSERT_TRUE(p.arc_ixs.empty());
}

TEST(gfapath_test, test) {
    graph g = make_graph();
    path p(graph::seg_vtx_p(2), 2);
    ASSERT_EQ(p.len, 0);
    ASSERT_EQ(p.start, graph::seg_vtx_p(2)<<32|2);
    ASSERT_EQ(p.now_at, p.start);
    p.add_arc(g, 5);
    ASSERT_EQ(p.now_at, graph::seg_vtx_p(0)<<32|1);
    ASSERT_EQ(p.len, 3);
    ASSERT_EQ(p.arc_ixs.size(), 1);
    p.add_arc(g, 0);
    ASSERT_EQ(p.now_at, graph::seg_vtx_n(1)<<32|4);
    ASSERT_EQ(p.len, 8);
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
