/* gfagraph.cpp
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

#include "gfagraph.h"

#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <ostream>
#include "gfa2logic.h"
#include "utils.h"

namespace gfa {

using gene_paths::raise_error;

const std::size_t path::START = std::uint64_t(-1);

static bool // for lower_bound - returns true when it is before v_lv
v_lv_less_l(const arc& it, std::uint64_t v_lv)
{
    return it.v_lv < v_lv;
}

static bool // for upper_bound - returns true when v_lv is before it
arc_less_u(const arc& a, const arc& it)
{
    return a.v_lv < it.v_lv || (a.v_lv == it.v_lv && a.w_lw < it.w_lw);
}

static const char RC_MAP[256] = {
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 , 'T', 'V', 'G', 'H',  0 ,  0 , 'C', 'D',  0 ,  0 , 'M',  0 , 'K', 'N',  0 ,
  0 ,  0 , 'Y', 'W', 'A',  0 , 'B', 'S',  0 , 'R',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 , 't', 'v', 'g', 'h',  0 ,  0 , 'c', 'd',  0 ,  0 , 'm',  0 , 'k', 'n',  0 ,
  0 ,  0 , 'y', 'w', 'a',  0 , 'b', 's',  0 , 'r',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0
};

std::ostream& 
seg::write_seq(std::ostream& os, bool rc, std::uint32_t beg, std::uint32_t end) const
{
    if (end == std::uint32_t(-1))
        end = len;

    const char *p0 = data.data() + beg;
    const char *p1 = data.data() + end;

    if (rc)
        while (p1 != p0) os << RC_MAP[std::size_t(*--p1)]; 
    else
        while (p0 != p1) os << *p0++;

    return os;
}

const seg*
graph::find_seg(const std::string& name) const
{
    const auto it = seg_ixs.find(name);
    return it == seg_ixs.cend() ? 0 : &(segs[it->second]);
}

std::size_t
graph::find_seg_ix(const std::string& name) const
{
    const auto it = seg_ixs.find(name);
    return it == seg_ixs.cend() ? std::size_t(-1) : it->second;
}

std::size_t
graph::get_seg_ix(const std::string& name) const
{
    const auto it = seg_ixs.find(name);
    if (it == seg_ixs.cend())
        raise_error("unknown segment: %s", name.c_str());
    return it->second;
}

void
graph::add_seg(const seg& s)
{
    if (s.name.empty())
        raise_error("segment name is empty");

    if (s.len != s.data.length())
        raise_error("segment length in GFA (%d) differs from FASTA (%d) for seqid %s", s.len, s.data.length(), s.name.c_str());

    if (seg_ixs.find(s.name) != seg_ixs.end())
        raise_error("duplicate segment name: %s", s.name.c_str());

    seg_ixs[s.name] = segs.size();
    segs.push_back(s);
}

void
graph::add_edge(const std::string& sref, std::uint32_t sbeg, std::uint32_t send,
                const std::string& dref, std::uint32_t dbeg, std::uint32_t dend)
{
        // determine orientations and segment names

    auto ps = sref.cend() - 1;
    auto pd = dref.cend() - 1;

    if (*ps != '-' && *ps != '+')
        raise_error("sequence reference without sign: %s", sref.c_str());

    if (*pd != '-' && *pd != '+')
        raise_error("sequence reference without sign: %s", dref.c_str());

    std::uint64_t s_neg = *ps == '-' ? 1 : 0;
    std::uint64_t d_neg = *pd == '-' ? 1 : 0;

        // look up segments

    std::string s_name(sref.cbegin(), ps);
    std::string d_name(dref.cbegin(), pd);

    std::size_t s_ix = get_seg_ix(s_name);
    std::size_t d_ix = get_seg_ix(d_name);

    const seg& s_seg = segs[s_ix];
    const seg& d_seg = segs[d_ix];

        // create the GFA2 representation having the logic

    gfa2::edge edge = {
        { sref, std::uint32_t(s_seg.len), sbeg, send, bool(!s_neg) },
        { dref, std::uint32_t(d_seg.len), dbeg, dend, bool(!d_neg) }
    };

    edge.validate();

        // add the paths from v to w

    std::uint64_t v_id = (s_ix<<1)|s_neg;
    std::uint64_t w_id = (d_ix<<1)|d_neg;

    arc as[4];

        // the left-to-right first turn-off from v and w and v.v.

    as[0].v_lv = as[1].w_lw = v_id<<32 | edge.lv();
    as[0].w_lw = as[1].v_lv = w_id<<32 | edge.lw();

        // same for the complement arcs

    as[2].v_lv = as[3].w_lw = (v_id^1)<<32 | edge.rv();
    as[2].w_lw = as[3].v_lv = (w_id^1)<<32 | edge.rw();

        // add them to the arcs array

    for (arc a : as) {
        arcs.emplace(std::upper_bound(arcs.cbegin(), arcs.cend(), a, arc_less_u), a);
    }

        // unless zero overlap, add the second turn-off for all

    if (edge.ov() != 0 || edge.ow() != 0) {
        for (arc a : as) {
            a.v_lv += edge.ov();
            a.w_lw += edge.ow();
            arcs.emplace(std::upper_bound(arcs.cbegin(), arcs.cend(), a, arc_less_u), a);
        }
    }
}

std::pair<std::vector<arc>::const_iterator, std::vector<arc>::const_iterator>
graph::arcs_from_v_lv(std::uint64_t v_lv) const
{
    arc next = { ((v_lv>>32)+1)<<32, 0 };
    std::vector<arc>::const_iterator lo = std::lower_bound(arcs.cbegin(), arcs.cend(), v_lv, v_lv_less_l);

    return std::make_pair(lo, std::upper_bound(lo, arcs.cend(), next, arc_less_u));
}

std::size_t
graph::start_path(std::uint32_t vtx_id, std::uint32_t pos)
{
    // We allocate the array path_starts for a fixed size, arbitrarily
    // chosen so that one path could start from every segment.
    // We want it fixed size because we use pointers to the arcs in it,
    // and growing it will reallocate and invalidate those.
    // Not ideal but works for now, and can be optimised later.
    if (path_starts.empty())
        path_starts.reserve(segs.size());
    else if (path_starts.size() == path_starts.capacity())
        raise_error("sorry, start_path array exhausted");

    // The path start is a 'pseudo-arc' arriving on w=vtx_id at lw=pos,
    // coming from v_lv at that same location.
    arc arc0;
    arc0.v_lv = std::uint64_t(vtx_id)<<32 | pos;
    arc0.w_lw = std::uint64_t(vtx_id)<<32 | pos;

    // Add the pseudo arc to the path_starts array
    path_starts.push_back(arc0);

    // Create the new path and add to the paths vector
    path p;
    p.pre_ix = path::START;
    p.p_arc = reinterpret_cast<const arc*>(&*path_starts.crbegin());

    paths.push_back(p);
    return paths.size() - 1;
}

std::ostream&
graph::write_path_seq(std::ostream& os, std::size_t path_ix) const
{
    const path& p = paths.at(path_ix);

    if (p.pre_ix != path::START)
    {
        write_path_seq(os, p.pre_ix);

        const path& p0 = paths.at(p.pre_ix);

        std::uint64_t vtx_ix = p.p_arc->v_lv >> 32;

        segs.at(graph::vtx_seg(vtx_ix)).write_vtx(os, (vtx_ix & 1) == 1, std::uint32_t(p0.p_arc->w_lw), std::uint32_t(p.p_arc->v_lv));
    }

    return os;
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
