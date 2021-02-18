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

/* This unit defines the data structures to hold paths over a graph.
 *
 * We define a path recursively as the "null" path or a path followed by
 * an arc.  An arc, as defined in gfagraph, holds a source (v_lv) and
 * destination location (w_lv).  Both contain a vertex identifier (v, w),
 * and a distance from the start of the vertex (lv, lw).
 *
 * A way to think about paths is to view a gfa::graph as a metro network,
 * where the segments are lines (that travel in one direction, and where
 * you can get off at any point), and each arc is a station where you can
 * hop onto another line.
 *
 * A path then is a sequence of rides (on a segment) and hops (between
 * segments).  We do not store the rides, only the hops.  The rides are
 * implicit - they go from the w_lw of one hop to the v_lv of the next.
 *
 * This recursive definition of paths allows for very concise storage:
 * we keep an array of 'path_arc' objects, each representing a new path
 * by having a pointer to an existing path, and a pointer to an arc that
 * that extends the path to be a new path.
 */

struct path_arc {
    std::size_t pre_ix;     // index of preceding path in paths or 0
    const arc* p_arc;       // points to the arc added to make this path
};

/* The paths struct holds any number of paths defined over a graph.
 *
 * To start a new set of paths from some location (position on a vertex),
 * call start_path(vtx, pos).  To create a path that extends an existing
 * path, call extend(path_ix, arc).
 */
struct paths {

    const graph& g;
    std::vector<arc> path_starts;
    std::vector<path_arc> path_arcs;

    paths(const graph& gr, std::size_t max)
        : g(gr) { 
        path_starts.reserve(max);
        path_arcs.push_back( {0,0} /* the 'null' path_arc at path_ix 0 */ );
    }

    // starts a path at v_lv (vtx<<32|pos)
    std::size_t start_path(std::uint64_t v_lv);

    // creates new path that extends path_ix with the arc at p_arc
    inline void extend(std::size_t path_ix, const arc *p_arc) {
        path_arcs.push_back( { path_ix, p_arc } );
    }

    // creates new path that extends path_ix with the arc at it
    inline void extend(std::size_t path_ix, std::vector<arc>::const_iterator it) {
        extend(path_ix, reinterpret_cast<const arc*>(&*it));
    }

    // returns the length of the ride from previous hop to current hop
    inline std::size_t ride_len(const path_arc& p) const {
        return p.pre_ix == 0 ? 0 : p.p_arc->v_lv - path_arcs.at(p.pre_ix).p_arc->w_lw;
    }

    // return the length of the path
    std::size_t length(const path_arc& p) const;

    // write the path route for path path_ix to an ostream
    std::ostream& write_route(std::ostream& os, std::size_t path_ix) const;

    // write the path sequence with id path_ix to an ostream
    std::ostream& write_seq(std::ostream& os, std::size_t path_ix) const;
};


} // namespace gfa

#endif // gfapaths_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
