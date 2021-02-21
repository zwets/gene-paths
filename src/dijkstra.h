/* dijkstra.h
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
#ifndef dijkstra_h_INCLUDED
#define dijkstra_h_INCLUDED

#include <vector>
#include "gfapaths.h"

namespace gfa {

// Return the path_ix of the shortest path between v_lv0 and v_lv1
std::size_t shortest_path(paths&, std::uint64_t v_lv0, std::uint64_t v_lv1);

} // namespace gfa

#endif // dijkstra_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
