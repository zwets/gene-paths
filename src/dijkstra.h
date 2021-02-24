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
    std::size_t found;  // holds the index into ps when path is found

    dijkstra(const graph& gr)
        : g(gr), ps(g), found(0) { restart(); }

        // finder functions

    inline void all_paths(const arc* start) { find_paths(start); }
    inline bool shortest_path(const arc* start, const arc* end) { return find_paths(start, end); }

        // convenience retrieval of route, sequence and length of a/the found path

    std::size_t length(std::size_t path_ix = std::size_t(-1)) const
        { return ps.length(ps.path_arcs.at(path_ix == std::size_t(-1) ? found : path_ix)); }

    std::string route(std::size_t path_ix = std::size_t(-1)) const
        { return ps.route(ps.path_arcs.at(path_ix == std::size_t(-1) ? found : path_ix)); }

    std::ostream& write_route(std::ostream& os, std::size_t path_ix = std::size_t(-1)) const
        { return ps.write_route(os, ps.path_arcs.at(path_ix == std::size_t(-1) ? found : path_ix)); }

    std::string sequence(std::size_t path_ix = std::size_t(-1)) const
        { return ps.sequence(ps.path_arcs.at(path_ix == std::size_t(-1) ? found : path_ix)); }

    std::ostream& write_sequence(std::ostream& os, std::size_t path_ix = std::size_t(-1)) const
        { return ps.write_seq(os, ps.path_arcs.at(path_ix == std::size_t(-1) ? found : path_ix)); }

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
