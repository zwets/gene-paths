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
#include <map>
#include <set>
#include "graph.h"
#include "paths.h"

namespace gfa {


// dijkstra - structure to perform shortest path searches on a graph

struct dijkstra
{
    const graph& g;
    paths ps;
    std::size_t found_pix;  // holds the index into ps when path is found
    std::size_t found_len;  // holds the length of the path that was found

    dijkstra(const graph& gr)
        : g(gr), ps(g) { restart(); }

        // finder functions

    // shortest path from start to end arc, false if no path, sets found to its index
    inline bool shortest_path(const arc* start, const arc* end) { return find_paths(start, end); }

    // find the shortest paths from start to every destination in the graph, put their indices in ps
    inline void shortest_paths(const arc* start) { find_paths(start); }  // to every destination

    // find the shortest path to the destination arc that is furthest from start
    // NOTE: we do not currently detect or flag circular paths
    void furthest_path(const arc* start);

//  // find the shortest path between the two arcs that are furthest apart, by iterating
//  // (NON-OPTIMISED!) over all possible start arcs
//  void furthest_path();

        // retrieval of route, sequence and length of a path

    // return the length of the path with index p_ix, or the found path by default
    std::size_t length(std::size_t p_ix = std::size_t(-1)) const
    { return p_ix == std::size_t(-1) ? found_len : ps.length(ps.path_arcs.at(p_ix)); }

    // return the route string of the path with index p_ix, or of the found path by default
    std::string route(std::size_t p_ix = std::size_t(-1)) const
    { return ps.route(ps.path_arcs.at(p_ix == std::size_t(-1) ? found_pix : p_ix)); }

    // write the route string of the path with index p_ix, or the found path by default, to os
    std::ostream& write_route(std::ostream& os, std::size_t p_ix = std::size_t(-1)) const
    { return ps.write_route(os, ps.path_arcs.at(p_ix == std::size_t(-1) ? found_pix : p_ix)); }

    // return the sequence of the path with index p_ix, or of the found path by default
    std::string sequence(std::size_t p_ix = std::size_t(-1)) const
    { return ps.sequence(ps.path_arcs.at(p_ix == std::size_t(-1) ? found_pix : p_ix)); }

    // write the sequence of the path with index p_ix, or the found path by default, to os
    std::ostream& write_sequence(std::ostream& os, std::size_t p_ix = std::size_t(-1)) const
    { return ps.write_seq(os, ps.path_arcs.at(p_ix == std::size_t(-1) ? found_pix : p_ix)); }

#ifdef NDEBUG
    private:    // hide implementation detail as private unless debugging
#endif
        // clears all data structures for another search
        void restart(const arc* start = 0);

        // the core finder function
        bool find_paths(const arc* start, const arc* end = 0);

        // dnode - pointer to current shortest path to a destination
        struct dnode {

            std::size_t len;        // total path length upto path p_ref
            std::uint64_t p_ref;    // high bit marks visited, rest indexes into ps

            inline std::uint64_t p_ix() const { return p_ref & 0x7FFFFFFFFFFFFFFFL; }
            inline bool is_visited() const { return p_ref>>63; }
            inline void mark_visited() { p_ref |= 0x8000000000000000L; }
        };

        // ds - map to look up the dnode for each destination (w_lw)
        typedef std::map<std::uint64_t, dnode> dmap_t;
        dmap_t ds;

        // comparator for keeping the vs ordered on length
        struct nearest_dnode {
            bool operator()(dmap_t::iterator const& i1, dmap_t::iterator const& i2) const noexcept {
                return i1->second.len < i2->second.len || (i1->second.len == i2->second.len && i1->first < i2->first);
            }
        };

        // vs - ordered set of visitable nodes sorted on increasing path length
        std::set<dmap_t::iterator, nearest_dnode> vs;

        // pops the nearest visitable off the vs
        dnode& pop_visit();
};


} // namespace gfa

#endif // dijkstra_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
