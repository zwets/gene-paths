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
#include <string>
#include "gfagraph.h"

#ifndef NDEBUG
#include "utils.h"
using gene_paths::raise_error;
#endif // NDEBUG

namespace gfa {

/* This unit defines the data structures to hold paths over a graph.
 *
 * We define a path recursively as the "null" path or a path followed by
 * an arc.  An arc, as defined in gfagraph, holds a source (v_lv) and
 * destination location (w_lv).  Locations consist of a vertex id (v, w),
 * and an offset from the start of the vertex (lv, lw).
 *
 * A way to think about paths is to view a gfa::graph as a metro network,
 * where the segments are lines (that travel in one direction, and where
 * you can get off at any point), and each arc is a station where you can
 * switch onto another line.
 *
 * A path then is a sequence of rides (on a segment) and hops (between
 * segments).  We do not store the rides, only the hops.  The rides are
 * implicit - they go from the w_lw of one hop to the v_lv of the next.
 *
 * This recursive definition of paths allows for very concise storage:
 * we keep an array of 'path_arc' objects, each representing a new path
 * by pointing back to an existing path, and adding an arc that extends
 * the path.
 *
 * When we later run Dijkstra's algorithm to find the shortest path, the
 * vertices are the w_lw (ends) of the arcs, and the path to a vertex
 * can be updated to a shorter path by changing its back pointer.
 */

struct path_arc {
    std::size_t pre_ix;     // index of preceding path in paths or 0
    const arc* p_arc;       // pointer to the arc its that path

        // convenience selectors of the src (v) and dst (w) fields

    inline std::uint64_t src_v() const { return p_arc->v_lv >> 32; }
    inline std::uint64_t src_lv() const { return p_arc->v_lv & 0xFFFFFFFFL; }
    inline std::uint64_t dst_v() const { return p_arc->w_lw >> 32; }
    inline std::uint64_t dst_lv() const { return p_arc->w_lw & 0xFFFFFFFFL; }

        // convenience aliases for the lengthy this->p_arc->v_lw

    inline std::uint64_t v_lv() const { return p_arc->v_lv; }
    inline std::uint64_t w_lw() const { return p_arc->w_lw; }
};

/* The paths struct holds any number of paths defined over a graph.
 *
 * Its core operation is extend(path_ix, p_arc), which adds a path_arc
 * that extends the existing path at path_ix with the arc pointed at by
 * p_arc, producing a new path.
 *
 * The "null" path at path_ix 0 signifies the start of a path.  Therefore
 * to create a new path starting with some arc* pa use extend(0, pa).
 */
struct paths {

    const graph& g;
    std::vector<path_arc> path_arcs;

    paths(const graph& gr)
        : g(gr) {
        path_arcs.push_back( {0,0} /* the 'null' path_arc at path_ix 0 */ );
    }

    // resets to empty
    inline void clear() { path_arcs.clear(); path_arcs.push_back({0,0}); }

    // selector for the path_arc at p_ix, just forwards
    inline const path_arc& at(std::size_t ix) const { return path_arcs.at(ix); }
    inline path_arc& at(std::size_t ix) { return path_arcs.at(ix); }

    // creates new path that extends path_ix with the arc at p_arc
    inline std::size_t extend(std::size_t path_ix, const arc *p_arc) {
#ifndef NDEBUG
        if (path_ix && p_arc->v() != path_arcs.at(path_ix).p_arc->w() )
            raise_error("programmer error: invalid path extension");
#endif
        path_arcs.push_back({ path_ix, p_arc });
        return path_arcs.size() - 1;
    }

    // returns the length of the 'ride' from previous arc to current arc
    inline std::size_t ride_len(const path_arc& p) const {
        return p.pre_ix ? p.p_arc->v_lv - path_arcs.at(p.pre_ix).p_arc->w_lw : 0;
    }

    // return the length of the path
    std::size_t length(const path_arc& p) const;

    // write the path route for p to an ostream or string
    std::ostream& write_route(std::ostream& os, const path_arc& p) const;
    std::string route(const path_arc& p) const;

    // write the path sequence for p to an ostream
    std::ostream& write_seq(std::ostream& os, const path_arc& p) const;
    std::string sequence(const path_arc& p) const;
};


} // namespace gfa

#endif // gfapaths_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
