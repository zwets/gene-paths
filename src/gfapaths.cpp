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

std::size_t
paths::start_path(std::uint32_t vtx_id, std::uint32_t pos)
{
    // We allocate the array path_starts for a fixed size, arbitrarily
    // chosen so that one path could start from every segment.
    // We want it fixed size because we use pointers to the arcs in it,
    // and growing it will reallocate and invalidate those.
    // Not ideal but works for now, and can be optimised later.
    if (path_starts.size() == path_starts.capacity())
        raise_error("sorry, start_path array exhausted");

    // The path start is a 'pseudo-arc' arriving on w=vtx_id at lw=pos,
    // coming from v_lv at that same location.
    arc arc0;
    arc0.v_lv = std::uint64_t(vtx_id)<<32 | pos;
    arc0.w_lw = std::uint64_t(vtx_id)<<32 | pos;

    // Add the pseudo arc to the path_starts array
    path_starts.push_back(arc0);

    // Create the new path and add to the paths vector
    path_arc p;
    p.pre_ix = 0;
    p.p_arc = reinterpret_cast<const arc*>(&*path_starts.crbegin());

    path_arcs.push_back(p);
    return path_arcs.size() - 1;
}

std::ostream&
paths::write_path_seq(std::ostream& os, std::size_t path_ix) const
{
    const path_arc& p = path_arcs.at(path_ix);

    if (p.pre_ix)
    {
        write_path_seq(os, p.pre_ix);

        const path_arc& p0 = path_arcs.at(p.pre_ix);

        std::uint64_t vtx_ix = p.p_arc->v_lv >> 32;

        g.segs.at(graph::vtx_seg(vtx_ix)).write_vtx(os, (vtx_ix & 1) == 1, std::uint32_t(p0.p_arc->w_lw), std::uint32_t(p.p_arc->v_lv));
    }

    return os;
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
