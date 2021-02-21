/* dijkstra.cpp
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

#include "dijkstra.h"

#include "gfapaths.h"
#include "utils.h"

namespace gfa {

void
dijkstra::all_paths(const arc* start)
{
    ps.path_arcs.clear();
    ps.extend(0, start);
    found = 1;
}

void
dijkstra::shortest_path(const arc* start, const arc* end)
{
    ps.path_arcs.clear();
    ps.extend(0, start);
    found = start == end;
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
