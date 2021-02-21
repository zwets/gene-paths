/* gfapaths.cpp
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

#include "gfapaths.h"

#include <vector>
#include <ostream>
#include "utils.h"

namespace gfa {

using gene_paths::raise_error;
using gene_paths::verbose_emit;

std::size_t
paths::length(const path_arc& p) const
{
    return p.pre_ix ? length(path_arcs.at(p.pre_ix)) + ride_len(p) : 0;
}

std::ostream&
paths::write_seq(std::ostream& os, std::size_t path_ix) const
{
    const path_arc& p = path_arcs.at(path_ix);

    // unless we are the start arc
    if (p.pre_ix)
    {
        // write the path leading up to path_ix
        write_seq(os, p.pre_ix);

        // retrieve the vertex (metro line) that path_ix hops off from
        std::uint64_t v = graph::get_v(p.p_arc->v_lv);

        // write the sequence of vertex from where we got on it (w_lw of
        // the previous arc) up to where we hop off it (v_lv of its arc)
        g.get_seg(graph::vtx_seg(v))
            .write_vtx(os, graph::is_neg(v),
                    graph::get_lv(path_arcs.at(p.pre_ix).p_arc->w_lw),
                    graph::get_lv(p.p_arc->v_lv));
    }

    return os;
}

std::ostream&
paths::write_route(std::ostream& os, std::size_t path_ix) const
{
    const path_arc& p = path_arcs.at(path_ix);

    // unless we are the start arc
    if (p.pre_ix)
    {
        // write the route path leading up to path_ix
        write_route(os, p.pre_ix);

        // retrieve the vertex (metro line) and where we got on and off
        std::uint64_t v = graph::get_v(p.p_arc->v_lv);
        std::uint64_t e = graph::get_lv(p.p_arc->v_lv);
        std::uint64_t b = graph::get_lv(path_arcs.at(p.pre_ix).p_arc->w_lw);
        const seg& s = g.get_seg(graph::vtx_seg(v));

        // write separator if we are not the first route element
        if (path_arcs.at(p.pre_ix).pre_ix != 0)
            os << ' ';

        // write segment name and orientation
        os << s.name << (graph::is_pos(v) ? '+' : '-');

        // unless the whole segment is traversed, write the ride
        if (b != 0 || e != s.len)
            os  << ':' << (graph::is_pos(v) ? b : s.len-e)
                << ':' << (graph::is_pos(v) ? e : s.len-b);
    }

    return os;
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
