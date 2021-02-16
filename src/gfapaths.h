/* gfapaths.h
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
#ifndef gfapaths_h_INCLUDED
#define gfapaths_h_INCLUDED

#include <vector>
#include "gfagraph.h"

namespace gfa {

/* path_arc
 *
 * Recursively defines a path as an existing path extended with an arc.
 *
 * The existing path is referenced by its index in path_arcs, or the
 * reserved value 0 for the "null" path at the start.
 *
 * The extension arc is a pointer to an arc in a gfa::graph.
 */
struct path_arc {
    std::size_t pre_ix;     // index of preceding path in paths or 0
    const arc* p_arc;       // points to the arc added to make this path
};

/* paths
 *
 * Holds paths defined over a graph.
 */
struct paths {

    const graph& g;
    std::vector<arc> path_starts;
    std::vector<path_arc> path_arcs;

    paths(const graph& gr, std::size_t max)
        : g(gr) { path_starts.reserve(max); path_arcs.push_back({0,0}); }

    // starts path at pos on vtx_ix, returns path_ix
    std::size_t start_path(std::uint32_t vtx_ix, std::uint32_t pos);

    // creates new path that extends path_ix with the arc at arc_it
    inline void grow_path(std::size_t path_ix, std::vector<arc>::const_iterator arc_it) {
        path_arc p = { path_ix, reinterpret_cast<const arc*>(&*arc_it) };
        path_arcs.push_back(p);
    }

    // write the path sequence with id path_ix to an ostream
    std::ostream& write_path_seq(std::ostream& os, std::size_t path_ix) const;
};


} // namespace gfa

#endif // gfapaths_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
