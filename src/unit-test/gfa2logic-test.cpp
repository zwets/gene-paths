/* gfa2logic-test.cpp
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
    ASSERT_EQ(v.o(), v.l);
    ASSERT_EQ(v.l1(), 0);
    ASSERT_EQ(v.l2(), v.l);
    ASSERT_EQ(v.r1(), 0);
    ASSERT_EQ(v.r2(), v.l);
    ASSERT_EQ(v.l1i(), 0);
    ASSERT_EQ(v.l2i(), v.l);
    ASSERT_EQ(v.r1i(), 0);
    ASSERT_EQ(v.r2i(), v.l);
}

TEST(gfa2logic_test, contained_vtx_inv) {
    vtx v = { sn, 3, 0, 3, false };
    v.validate();
    ASSERT_EQ(v.o(), v.l);
    ASSERT_EQ(v.l1(), 0);
    ASSERT_EQ(v.l2(), v.l);
    ASSERT_EQ(v.r1(), 0);
    ASSERT_EQ(v.r2(), v.l);
    ASSERT_EQ(v.l1i(), 0);
    ASSERT_EQ(v.l2i(), v.l);
    ASSERT_EQ(v.r1i(), 0);
    ASSERT_EQ(v.r2i(), v.l);
}

TEST(gfa2logic_test, right_blunt) {
    vtx v = { s, 5, 5, 5, true };
    v.validate();
    ASSERT_EQ(v.o(), 0);
    ASSERT_EQ(v.l1(), v.l);
    ASSERT_EQ(v.l2(), v.l);
    ASSERT_EQ(v.r1(), 0);
    ASSERT_EQ(v.r2(), 0);
    ASSERT_EQ(v.l1i(), 0);
    ASSERT_EQ(v.l2i(), 0);
    ASSERT_EQ(v.r1i(), v.l);
    ASSERT_EQ(v.r2i(), v.l);
}

TEST(gfa2logic_test, right_blunt_inv) {
    vtx v = { sn, 5, 0, 0, false };
    v.validate();
    ASSERT_EQ(v.o(), 0);
    ASSERT_EQ(v.l1(), v.l);
    ASSERT_EQ(v.l2(), v.l);
    ASSERT_EQ(v.r1(), 0);
    ASSERT_EQ(v.r2(), 0);
    ASSERT_EQ(v.l1i(), 0);
    ASSERT_EQ(v.l2i(), 0);
    ASSERT_EQ(v.r1i(), v.l);
    ASSERT_EQ(v.r2i(), v.l);
}

TEST(gfa2logic_test, left_blunt) {
    vtx v = { s, 5, 0, 0, true };
    v.validate();
    ASSERT_EQ(v.o(), 0);
    ASSERT_EQ(v.l1(), 0);
    ASSERT_EQ(v.l2(), 0);
    ASSERT_EQ(v.r1(), v.l);
    ASSERT_EQ(v.r2(), v.l);
    ASSERT_EQ(v.l1i(), v.l);
    ASSERT_EQ(v.l2i(), v.l);
    ASSERT_EQ(v.r1i(), 0);
    ASSERT_EQ(v.r2i(), 0);
}

TEST(gfa2logic_test, left_blunt_inv) {
    vtx v = { sn, 5, 5, 5, false };
    v.validate();
    ASSERT_EQ(v.o(), 0);
    ASSERT_EQ(v.l1(), 0);
    ASSERT_EQ(v.l2(), 0);
    ASSERT_EQ(v.r1(), v.l);
    ASSERT_EQ(v.r2(), v.l);
    ASSERT_EQ(v.l1i(), v.l);
    ASSERT_EQ(v.l2i(), v.l);
    ASSERT_EQ(v.r1i(), 0);
    ASSERT_EQ(v.r2i(), 0);
}

TEST(gfa2logic_test, right_dovetail) {
    vtx v = { s, 5, 3, 5, true };  // ---==
    v.validate();
    ASSERT_EQ(v.o(), 2);
    ASSERT_EQ(v.l1(), 3);  // ---==
    ASSERT_EQ(v.l2(), 5);
    ASSERT_EQ(v.r1(), 0);
    ASSERT_EQ(v.r2(), 2);
    ASSERT_EQ(v.l1i(), 0); // ==----
    ASSERT_EQ(v.l2i(), 2);
    ASSERT_EQ(v.r1i(), 3);
    ASSERT_EQ(v.r2i(), 5);
}

TEST(gfa2logic_test, right_dovetail_inv) {
    vtx v = { sn, 5, 0, 2, false }; // ---== (==---)
    v.validate();
    ASSERT_EQ(v.o(), 2);
    ASSERT_EQ(v.l1(), 3);  // ---==
    ASSERT_EQ(v.l2(), 5);
    ASSERT_EQ(v.r1(), 0);
    ASSERT_EQ(v.r2(), 2);
    ASSERT_EQ(v.l1i(), 0); // ==---
    ASSERT_EQ(v.l2i(), 2);
    ASSERT_EQ(v.r1i(), 3);
    ASSERT_EQ(v.r2i(), 5);
}

TEST(gfa2logic_test, left_dovetail) {
    vtx v = { s, 5, 0, 2, true }; // ==----
    v.validate();
    ASSERT_EQ(v.o(), 2);
    ASSERT_EQ(v.l1(), 0);  // ==---
    ASSERT_EQ(v.l2(), 2);
    ASSERT_EQ(v.r1(), 3);
    ASSERT_EQ(v.r2(), 5);
    ASSERT_EQ(v.l1i(), 3); // ---==
    ASSERT_EQ(v.l2i(), 5);
    ASSERT_EQ(v.r1i(), 0);
    ASSERT_EQ(v.r2i(), 2);
}

TEST(gfa2logic_test, left_dovetail_inv) {
    vtx v = { sn, 5, 3, 5, false };
    v.validate();
    ASSERT_EQ(v.o(), 2);
    ASSERT_EQ(v.l1(), 0);
    ASSERT_EQ(v.l2(), 2);
    ASSERT_EQ(v.r1(), 3);
    ASSERT_EQ(v.r2(), 5);
    ASSERT_EQ(v.l1i(), 3);
    ASSERT_EQ(v.l2i(), 5);
    ASSERT_EQ(v.r1i(), 0);
    ASSERT_EQ(v.r2i(), 2);
}

TEST(gfa2logic_test, containing_vtx) {
    vtx v = { s, 6, 1, 3, true };
    v.validate();
    ASSERT_EQ(v.o(), 2);
    ASSERT_EQ(v.l1(), 1);
    ASSERT_EQ(v.l2(), 3);
    ASSERT_EQ(v.r1(), 3);
    ASSERT_EQ(v.r2(), 5);
    ASSERT_EQ(v.l1i(), 3);
    ASSERT_EQ(v.l2i(), 5);
    ASSERT_EQ(v.r1i(), 1);
    ASSERT_EQ(v.r2i(), 3);
}

TEST(gfa2logic_test, containing_vtx_inv) {
    vtx v = { sn, 6, 1, 3, false };  // ---==-  (-==---)
    v.validate();
    ASSERT_EQ(v.o(), 2);
    ASSERT_EQ(v.l1(), 3);  // ---==-
    ASSERT_EQ(v.l2(), 5);
    ASSERT_EQ(v.r1(), 1);
    ASSERT_EQ(v.r2(), 3);
    ASSERT_EQ(v.l1i(), 1); // -==---
    ASSERT_EQ(v.l2i(), 3);
    ASSERT_EQ(v.r1i(), 3);
    ASSERT_EQ(v.r2i(), 5);
}

TEST(gfa2logic_test, dovetail_edge) {
    edge e = {
     { s1p, 3, 2, 3, true },
     { s2n, 5, 3, 5, false }};
    e.validate();
    ASSERT_EQ(e.ov(), 1);
    ASSERT_EQ(e.lv(), 2);
    ASSERT_EQ(e.lv2(), 3);
    ASSERT_EQ(e.lvi(), 0);
    ASSERT_EQ(e.lv2i(), 1);
    ASSERT_EQ(e.ow(), 2);
    ASSERT_EQ(e.lw(), 0);
    ASSERT_EQ(e.lw2(), 2);
    ASSERT_EQ(e.lwi(), 3);
    ASSERT_EQ(e.lw2i(), 5);
}

TEST(gfa2logic_test, general_edge) {
    edge e = {
     { s2n, 6, 3, 5, false },   // neg -==---    (---==-)
     { s1p, 9, 2, 5, true }};   // pos --===----
    e.validate();
    ASSERT_EQ(e.ov(), 2);
    ASSERT_EQ(e.lv(), 1);
    ASSERT_EQ(e.lv2(), 3);
    ASSERT_EQ(e.lvi(), 3);
    ASSERT_EQ(e.lv2i(), 5);
    ASSERT_EQ(e.ow(), 3);
    ASSERT_EQ(e.lw(), 2);
    ASSERT_EQ(e.lw2(), 5);
    ASSERT_EQ(e.lwi(), 4);
    ASSERT_EQ(e.lw2i(), 7);
}


} // namespace
  // vim: sts=5:sw=4:ai:si:et
