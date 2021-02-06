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
#include "gfa2logic.h"

using namespace gfa2;

namespace {

std::string s("s+");
std::string sn("s-");
std::string s1p("s1+");
std::string s1n("s1-");
std::string s2p("s2+");
std::string s2n("s2-");

TEST(gfa2logic_test, empty_seq) {
    vtx v = { s, 0, 0, 0, true };
    ASSERT_EXIT( v.validate(),
            testing::ExitedWithCode(1), 
            ": error: segment length is 0 for vertex s\\+");
}

TEST(gfa2logic_test, wrong_sign) {
    vtx v = { s, 1, 0, 1, false };
    ASSERT_EXIT( v.validate(),
            testing::ExitedWithCode(1), 
            ": error: inconsistent name and orientation: s\\+ defined neg");
}

TEST(gfa2logic_test, begin_past_length) {
    vtx v = { s, 1, 2, 3, true };
    ASSERT_EXIT( v.validate(),
            testing::ExitedWithCode(1), 
            ": error: begin or end beyond segment length on vertex s\\+");
}

TEST(gfa2logic_test, end_past_length) {
    vtx v = { s, 1, 0, 3, true };
    ASSERT_EXIT( v.validate(),
            testing::ExitedWithCode(1), 
            ": error: begin or end beyond segment length on vertex s\\+");
}

TEST(gfa2logic_test, begin_after_end) {
    vtx v = { s, 1, 1, 0, true };
    ASSERT_EXIT( v.validate(),
            testing::ExitedWithCode(1), 
            ": error: begin past end on vertex s\\+");
}

TEST(gfa2logic_test, contained_vtx) {
    vtx v = { s, 3, 0, 3, true };
    v.validate();
    ASSERT_TRUE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), v.l);
    ASSERT_EQ(v.overhang_l(), 0);
    ASSERT_EQ(v.overhang_r(), 0);
}

TEST(gfa2logic_test, contained_vtx_inv) {
    vtx v = { sn, 3, 0, 3, false };
    v.validate();
    ASSERT_TRUE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), v.l);
    ASSERT_EQ(v.overhang_l(), 0);
    ASSERT_EQ(v.overhang_r(), 0);
}

TEST(gfa2logic_test, right_blunt) {
    vtx v = { s, 5, 5, 5, true };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_TRUE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 0);
    ASSERT_EQ(v.overhang_l(), v.l);
    ASSERT_EQ(v.overhang_r(), 0);
}

TEST(gfa2logic_test, right_blunt_inv) {
    vtx v = { sn, 5, 0, 0, false };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_TRUE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 0);
    ASSERT_EQ(v.overhang_l(), v.l);
    ASSERT_EQ(v.overhang_r(), 0);
}

TEST(gfa2logic_test, left_blunt) {
    vtx v = { s, 5, 0, 0, true };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_TRUE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 0);
    ASSERT_EQ(v.overhang_l(), 0);
    ASSERT_EQ(v.overhang_r(), v.l);
}

TEST(gfa2logic_test, left_blunt_inv) {
    vtx v = { sn, 5, 5, 5, false };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_TRUE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 0);
    ASSERT_EQ(v.overhang_l(), 0);
    ASSERT_EQ(v.overhang_r(), v.l);
}

TEST(gfa2logic_test, right_dovetail) {
    vtx v = { s, 5, 3, 5, true };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_TRUE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 2);
    ASSERT_EQ(v.overhang_l(), 3);
    ASSERT_EQ(v.overhang_r(), 0);
}

TEST(gfa2logic_test, right_dovetail_inv) {
    vtx v = { sn, 5, 0, 2, false };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_TRUE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 2);
    ASSERT_EQ(v.overhang_l(), 3);
    ASSERT_EQ(v.overhang_r(), 0);
}

TEST(gfa2logic_test, left_dovetail) {
    vtx v = { s, 5, 0, 2, true };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_TRUE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 2);
    ASSERT_EQ(v.overhang_l(), 0);
    ASSERT_EQ(v.overhang_r(), 3);
}

TEST(gfa2logic_test, left_dovetail_inv) {
    vtx v = { sn, 5, 3, 5, false };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_FALSE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_TRUE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 2);
    ASSERT_EQ(v.overhang_l(), 0);
    ASSERT_EQ(v.overhang_r(), 3);
}

TEST(gfa2logic_test, containing_vtx) {
    vtx v = { s, 6, 1, 3, true };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_TRUE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 2);
    ASSERT_EQ(v.overhang_l(), 1);
    ASSERT_EQ(v.overhang_r(), 3);
}

TEST(gfa2logic_test, containing_vtx_inv) {
    vtx v = { sn, 6, 1, 3, false };
    v.validate();
    ASSERT_FALSE(v.is_contained());
    ASSERT_TRUE(v.is_container());
    ASSERT_FALSE(v.is_blunt_r());
    ASSERT_FALSE(v.is_blunt_l());
    ASSERT_FALSE(v.dovetails_r());
    ASSERT_FALSE(v.dovetails_l());
    ASSERT_EQ(v.overlap(), 2);
    ASSERT_EQ(v.overhang_l(), 3);
    ASSERT_EQ(v.overhang_r(), 1);
}

TEST(gfa2logic_test, impossible_edge1) {
    edge e = {
     { s1p, 3, 2, 3, true },
     { s2p, 5, 3, 5, true }};
    ASSERT_EXIT( e.validate(),
            testing::ExitedWithCode(1), 
            ": error: impossible edge: s1\\+ to s2+");
}

TEST(gfa2logic_test, possible_edge1_inv) {
    edge e = {
     { s1p, 3, 2, 3, true },
     { s2n, 5, 3, 5, false }};
    e.validate();
}

TEST(gfa2logic_test, flippable_possible) {
    edge e = {
     { s2n, 5, 3, 5, false },
     { s1p, 3, 2, 3, true }};
    e.validate();
}

TEST(gfa2logic_test, edge_no_flip) {
    edge e = {
     { s1p, 3, 3, 3, true },
     { s2p, 5, 0, 2, true }};
    e.validate();
    ASSERT_FALSE(e.needs_flip());
}

TEST(gfa2logic_test, edge_inv_no_flip) {
    edge e = {
     { s1p, 3, 3, 3, true },
     { s2n, 5, 3, 5, false }};
    e.validate();
    ASSERT_FALSE(e.needs_flip());
}

TEST(gfa2logic_test, edge_needs_flip) {
    edge e = {
     { s1p, 3, 0, 1, true },
     { s2p, 5, 3, 5, true }};
    e.validate();
    ASSERT_FALSE(e.s.goes_left());
    ASSERT_FALSE(e.d.goes_right());
    ASSERT_TRUE(e.d.goes_left());
    ASSERT_TRUE(e.s.goes_right());
    ASSERT_TRUE(e.needs_flip());
}

TEST(gfa2logic_test, edge_inv_needs_flip) {
    edge e = {
     { s1p, 3, 0, 1, true },
     { s2n, 5, 0, 2, false }};
    e.validate();
    ASSERT_FALSE(e.s.goes_left());
    ASSERT_FALSE(e.d.goes_right());
    ASSERT_TRUE(e.d.goes_left());
    ASSERT_TRUE(e.s.goes_right());
    ASSERT_TRUE(e.needs_flip());
}

/*
TEST(gfa2logic_test, container_flip) {
    edge e = {
     { s1p, 3, 0, 3, true },
     { s2p, 5, 3, 5, true }};
    e.validate();
    ASSERT_TRUE(e.s.goes_left());
    ASSERT_TRUE(e.s.goes_right());
    ASSERT_TRUE(e.d.goes_left());
    ASSERT_FALSE(e.d.goes_right());
    ASSERT_TRUE(e.needs_flip());
    ASSERT_EQ(e.vtx_l().id, s2p);
    ASSERT_EQ(e.vtx_r().id, s1p);
}

TEST(gfa2logic_test, container_no_flip) {
    edge e = {
     { s1p, 3, 0, 3, true },
     { s2p, 5, 0, 2, true }};
    e.validate();
    ASSERT_TRUE(e.s.goes_left());
    ASSERT_TRUE(e.s.goes_right());
    ASSERT_TRUE(e.d.goes_right());
    ASSERT_FALSE(e.d.goes_left());
    ASSERT_FALSE(e.needs_flip());
    ASSERT_EQ(e.vtx_l().id, s1p);
    ASSERT_EQ(e.vtx_r().id, s2p);
}
*/

} // namespace
  // vim: sts=5:sw=4:ai:si:et
