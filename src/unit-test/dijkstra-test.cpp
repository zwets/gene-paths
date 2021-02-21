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
#include "dijkstra.h"

using namespace gfa;

namespace {

TEST(dijkstra_test, dummy) {
    graph g;
    ASSERT_EQ(g.segs.size(), 0);
    ASSERT_EQ(g.arcs.size(), 0);
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
