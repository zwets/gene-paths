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
paths::start_path(std::uint32_t vtx_id, std::uint32_t pos)
{
    // The array path_starts is fixed at construction so we can store
    // pointers to the arcs in it.  Growing it would invalidate these.
    // Not ideal but works for now, and can be optimised later.
    if (path_starts.size() == path_starts.capacity())
        raise_error("sorry, start_path array exhausted");

    // Add a 'pseudo arc' pointing at the start location to the
    // path starts array.
    path_starts.push_back({
        std::uint64_t(vtx_id)<<32 | pos,
        std::uint64_t(vtx_id)<<32 | pos
    });

    // Add the first path_arc of the path, having a null lead path, and
    // the pseudo arc pointing at its start location as its extension.
    extend(0, path_starts.cend() - 1);

    // Return the index of the new path
    return path_arcs.size() - 1;
}

std::size_t
paths::length(const path_arc& p) const
{
    return p.pre_ix ? length(path_arcs.at(p.pre_ix)) + ride_len(p) : 0;
}

std::ostream&
paths::write_path_seq(std::ostream& os, std::size_t path_ix) const
{
    const path_arc& p = path_arcs.at(path_ix);

    // Unless we are the start arc
    if (p.pre_ix)
    {
        // Write the path leading up to path_ix
        write_path_seq(os, p.pre_ix);

        // Retrieve the vertex (metro line) that path_ix hops off from
        std::uint64_t vtx_ix = p.p_arc->v_lv >> 32;

        // Write the sequence of vertex from where we got on it (w_lw of
        // the previous arc) up to where we hop off it (v_lv of its arc)
        g.segs.at(graph::vtx_seg(vtx_ix))
            .write_vtx(os, vtx_ix & 1, 
                std::uint32_t(path_arcs.at(p.pre_ix).p_arc->w_lw), 
                std::uint32_t(p.p_arc->v_lv));
    }

    return os;
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
