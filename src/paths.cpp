/* paths.cpp
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

#include "paths.h"

#include <vector>
#include <sstream>
#include <ostream>
#include "utils.h"

namespace gfa {

using gene_paths::raise_error;
using gene_paths::verbose_emit;

std::size_t
paths::length(const path_arc& tip) const
{
    // recursive is nicer, but iterate friendlier
    //return p.pre_ix ? ride_len(p) + length(path_arcs.at(p.pre_ix)) : 0;

    std::size_t len = 0L;
    const path_arc* p = &tip;

    while (p->pre_ix) {
        len += ride_len(*p);
        p = &path_arcs.at(p->pre_ix);
    }

    return len;
}

std::ostream&
paths::write_seq(std::ostream& os, const path_arc& p) const
{
    if (p.pre_ix) // unless we are the start arc
    {
            // recurse over the pre-path

        const path_arc& pp = path_arcs.at(p.pre_ix);
        write_seq(os, pp);

            // write the seq of the ride, is v from pp.dst to p.src

        std::uint64_t v = pp.dst_v(); // same as p.src_v()
        g.get_seg(graph::vtx_seg(v))
            .write_vtx(os, graph::is_neg(v), pp.dst_lv(), p.src_lv());
    }

    return os;
}

std::string
paths::sequence(const path_arc& p) const
{
    std::stringstream ss;
    write_seq(ss, p);
    return ss.str();
}

std::ostream&
paths::write_route(std::ostream& os, const path_arc& p) const
{
    if (p.pre_ix) { // unless we are the start arc

            // recurse into the pre-path (const to inline)

        const path_arc& pp = path_arcs.at(p.pre_ix);
        write_route(os, pp);

            // append the seg name of final ride on v

        const std::uint64_t v = p.src_v();
        const seg& s = g.get_seg(graph::vtx_seg(v));

        if (pp.pre_ix) os << ' ';
        os << s.name;

            // append section unless v was traversed all the way

        const std::uint64_t b = pp.dst_lv(), e = p.src_lv();

        if (b != 0 || e != s.len)
            os  << ':' << (graph::is_pos(v) ? b : s.len-e)
                << ':' << (graph::is_pos(v) ? e : s.len-b);

            // append orientation of v

        os << (graph::is_pos(v) ? '+' : '-');
    }

    return os;
}

std::string
paths::route(const path_arc& p) const
{
    std::stringstream ss;
    write_route(ss, p);
    return ss.str();
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
