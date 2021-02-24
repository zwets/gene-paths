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

#include "paths.h"
#include "utils.h"

namespace gfa {

using gene_paths::raise_error;
using gene_paths::verbose_emit;

void
dijkstra::restart(const arc* start)
{
    // clear the paths, dnodes, and visitables

    ps.clear();
    ds.clear();
    vs.clear();
    found = 0;

    // set up ds to have all w_lw (destinations), initialising them
    // with infinite length and a null path reference.

    std::pair<std::uint64_t, dnode> val = { 0, { std::size_t(-1), 0 } };

    for (auto it = g.arcs.cbegin(); it != g.arcs.cend(); ++it)
    {
        val.first = it->w_lw;
        ds.insert(val);
    }

    // if we have a start arc, add it to visitables
    if (start) {

        // add start arc to path_arcs in ps, it will have p_ix 1
        std::size_t p_ix = ps.extend(0, start);

        // look up the start arc destination in ds
        dmap_t::iterator d_it = ds.find(start->w_lw);
        if (d_it == ds.end())
            raise_error("start arc not found in graph");

        // update its dnode to have len 0 and p_ix 1
        d_it->second = { 0, p_ix };

        // add a visitable for the start arc
        vs.insert(d_it);
    }
}


dijkstra::dnode&
dijkstra::pop_visit()
{
    dmap_t::iterator top = *vs.begin();
    dnode& d = top->second;
#ifndef NDEBUG
    if (d.is_visited())
        raise_error("dijkstra: visitable dnode already visited (programmer error)");
    if (!d.p_ref)
        raise_error("dijkstra: visitable dnode without a p_ref (programmer error)");
    if (top->first != ps.path_arcs.at(d.p_ref).p_arc->w_lw)
        raise_error("dijkstra: visitable dnode indexed at wrong w_lw (programmer error)");
#endif
    vs.erase(top);
    return d;
}


bool
dijkstra::find_paths(const arc* start, const arc* end)
{
    restart(start);
    found = 0;

    while (!found && !vs.empty()) {

        // pick the next node to visit
        dnode& vn = pop_visit();    // has .len and .p_ref

        // retrieve the path index, path arc and len to arrive at vn
        std::size_t cur_pix = vn.p_ref;
        std::size_t cur_len = vn.len;
        const path_arc& cur_pa = ps.at(cur_pix);
#ifndef NDEBUG
        verbose_emit("start visit of p_ref %lu at %lu", cur_pix, cur_len);
#endif
        // the dest (w_lw) of that arc is the new start (v_lv)
        std::uint64_t v_lv = cur_pa.p_arc->w_lw;

        // get all arcs leaving from vn's vertex at lv or later
        const auto iters = g.arcs_from_v_lv(v_lv);

        // look at the arcs to each tentative destination in turn
        for (auto a_it = iters.first; a_it != iters.second; ++a_it) {

            // ignore any arc that would take us right back
            if (a_it->w_lw == cur_pa.p_arc->v_lv)
                continue;

            // Note how we iterate over outbound arcs, where added length
            // lies on vn's contig, and then a (zero-length) jump is made:
            //
            //            w: --1------2---o---->
            //                 |     /
            //    v: ---x======1----2------>
            //
            // We are at vn=x (the w_lw of vn's arc) and iterate the arcs
            // downstream on v (v1-w1 and v2-w2).  We compute the length of
            // the '===' segments as the added distance.

            // compute distance to the departing arc
            std::uint64_t add_len = a_it->v_lv - v_lv;

            // locate the dnode for the tentative destination, which will
            // have the shortest path found to it so far (INF if not seen)
            const dmap_t::iterator d_it = ds.find(a_it->w_lw);
            dnode& dn = d_it->second;

            // if the new path is shorter, update the tentative dest
            if (cur_len + add_len < dn.len) {
#ifndef NDEBUG
                if (dn.is_visited())
                    raise_error("dijkstra: visited node found with shorter path (programmer error)");
#endif
                // if we haven't seen this destination yet
                if (!dn.p_ref) {
                    // add a path_arc from us to it to ps
                    dn.p_ref = ps.extend(cur_pix, &*a_it);
#ifndef NDEBUG
                    verbose_emit("- extended with new p_ref %lu (+%lu)", dn.p_ref, add_len);
#endif
                }
                else {
                    // remove old record from the visitables list
                    vs.erase(d_it);
                    // repoint its pre-path to the vn, and set new arc
                    // note: it can't be the pre_ix of anything yet
                    path_arc& d_pa = ps.at(dn.p_ref);
                    d_pa.pre_ix = cur_pix;
                    d_pa.p_arc = &*a_it;
#ifndef NDEBUG
                    verbose_emit("- updated existing p_ref %lu (-%lu)", dn.p_ref, dn.len - (cur_len + add_len));
#endif
                }

                // update the dnode with the new shortest length
                dn.len = cur_len + add_len;

                // and (re)add it to the visitables
                vs.insert(d_it);

            } // end if shorter path

        } // end iterate over tentatives

        // the vn is now visited and the shortest path to its arc
        vn.mark_visited();

        // check if we are done, i.e. the vn took the end arc
        if (end && cur_pa.p_arc == end) {
            found = vn.p_ix();
            verbose_emit("shortest path of length %lu found (index %lu)", vn.len, found);
        }

    } // end while !found and !vs.empty()

    verbose_emit("done exploring %lu (potential) paths", ps.path_arcs.size());

    // return true if we found path or not end was specified (find all)
    return !end || found;
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
