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
#include <functional>
#include "gfagraph.h"
#include "gfapaths.h"

namespace gfa {


// dijkstra - structure to perform shortest path searches on a graph

struct dijkstra
{
    const graph& g;
    paths ps;
    std::size_t found;

    dijkstra(const graph& gr)
        : g(gr), ps(g), found(0) { restart(); }

    void all_paths(const arc* start);
    void shortest_path(const arc* start, const arc* end);

        // convenience retrieval of route and sequence for the found path

    std::string route(std::size_t path_ix = std::size_t(-1)) const
        { return ps.route(ps.path_arcs.at(path_ix == std::size_t(-1) ? found : path_ix)); }

    std::ostream& write_route(std::ostream& os, std::size_t path_ix = std::size_t(-1)) const
        { return ps.write_route(os, ps.path_arcs.at(path_ix == std::size_t(-1) ? found : path_ix)); }

    std::string sequence(std::size_t path_ix = std::size_t(-1)) const
        { return ps.sequence(ps.path_arcs.at(path_ix == std::size_t(-1) ? found : path_ix)); }

    std::ostream& write_sequence(std::ostream& os, std::size_t path_ix = std::size_t(-1)) const
        { return ps.write_seq(os, ps.path_arcs.at(path_ix == std::size_t(-1) ? found : path_ix)); }

        // implementation detail - private unless testing

#ifdef NDEBUG
    private:    // implementation detail, in production have private
#endif
        // clears all data structures for another search
        void restart(const arc* start = 0);

        // dnode, dmap_t, diter_t - support data structures

        struct dnode {
            std::size_t len;
            std::uint64_t p_ref;  // high bit marks visited, rest is p_ix into ps

            inline std::uint64_t p_ix() const { return p_ref & 0x7FFFFFFFFFFFFFFFL; }
            inline bool is_visited() const { return p_ref>>63; }
            inline void mark_visited() { p_ref |= 0x8000000000000000L; }
        };

        typedef std::map<std::uint64_t, dnode> dmap_t;

        // destination map indexed on w_lw (destination)
        dmap_t ds;

        // comparator for keeping the vs set ordered on length (and destination second)
        struct nearest_dnode { // for keeping the vs
            bool operator()(dmap_t::iterator const& i1, dmap_t::iterator const& i2) const noexcept {
                return i1->second.len < i2->second.len || (i1->second.len == i2->second.len && i1->first < i2->first);
            }
        };

        // visitables list sorted on shortest path
        std::set<dmap_t::iterator, nearest_dnode> vs;

        // pops nearest visitable off the vs
        dnode& pop_visit();
};


} // namespace gfa

#endif // dijkstra_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
