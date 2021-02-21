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
dijkstra::restart(const arc* start)
{
    ps.clear();
    ls.clear();
    found = 0;
    setup_ds();

    if (start) {
        ps.extend(0, start);
        auto d_it = ds.find(start->w_lw);
        if (d_it == ds.end())
            gene_paths::raise_error("start arc not found in graph");
        d_it->second = { 1, 0 };
        ls.insert({ 0, d_it });
    }
}

void
dijkstra::setup_ds()
{
    ds.clear();
    std::pair<std::uint64_t, dnode> val = { 0, { 0, std::size_t(-1) } };
    for (auto it = g.arcs.cbegin(); it != g.arcs.cend(); ++it) {
        val.first = it->w_lw;
        ds.insert(val);
    }
}

dijkstra::dnode&
dijkstra::pop_visit()
{
    auto top = ls.begin();
    dnode& d = top->second->second;
#ifndef NDEBUG
    if (top->first != d.len)
        gene_paths::raise_error("dijkstra: inconsistent len (programmer error)");
    if (top->second->first != ps.path_arcs.at(d.path_ix).p_arc->w_lw)
        gene_paths::raise_error("dijkstra: inconsistent w_lw (programmer error)");
#endif
    ls.erase(top);
    return d;
}

void
dijkstra::all_paths(const arc* start)
{
    restart(start);
    found = 1;

    while (!ls.empty()) {
        dnode& dn = pop_visit();
        // TODO more
    }
    // put start in the tentative nodes with len 0
    // while there are tentative nodes
    //  pick the one with the shortest path
    //    for each of its tentative destinations
    //      if it was visited, next
    //      else if it is tentative
    //        and this has a shorter path, update that
    //    mark this node as visited
}

void
dijkstra::shortest_path(const arc* start, const arc* end)
{
    restart(start);
    found = start == end ? 1 : 0;

    while (!found && !ls.empty()) {
        dnode& dn = pop_visit();
        // TODO more
        break;
    }
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
