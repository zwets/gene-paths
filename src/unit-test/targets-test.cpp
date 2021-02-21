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
#include "targets.h"

using namespace gfa;

namespace {

TEST(targets_test, no_pos) {
    gfa::target t = gfa::target::parse("CONTIG+");
    ASSERT_EQ(t.ctg, "CONTIG");
    ASSERT_FALSE(t.neg);
    ASSERT_EQ(t.beg, std::uint64_t(-1));
    ASSERT_EQ(t.end, std::uint64_t(-1));
}

TEST(targets_test, beg_only) {
    gfa::target t = gfa::target::parse("CONTIG-:0");
    ASSERT_EQ(t.ctg, "CONTIG");
    ASSERT_TRUE(t.neg);
    ASSERT_EQ(t.beg, 0);
    ASSERT_EQ(t.end, std::uint64_t(-1));
}

TEST(targets_test, both_pos) {
    gfa::target t = gfa::target::parse("CONTIG-:0:2");
    ASSERT_EQ(t.ctg, "CONTIG");
    ASSERT_TRUE(t.neg);
    ASSERT_EQ(t.beg, 0);
    ASSERT_EQ(t.end, 2);
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
