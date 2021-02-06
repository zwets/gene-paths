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
#include "gfa2logic.h"
#include "utils.h"

namespace gfa {

using gene_paths::raise_error;

static bool // for lower_bound - returns true when it is before v_lv
v_lv_less_l(const arc& it, std::uint64_t v_lv)
{
    return it.v_lv < v_lv;
}

static bool // for upper_bound - returns true when v_lv is before it
v_lv_less_u(std::uint64_t v_lv, const arc& it)
{
    return v_lv < it.v_lv;
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

    gfa2::edge g2e = {
        { sref, std::uint32_t(s_seg.len), sbeg, send, bool(!s_neg) },
        { dref, std::uint32_t(d_seg.len), dbeg, dend, bool(!d_neg) }
    };

    g2e.validate();

        // We ignore containment and containers (for now)

    if (g2e.s.is_contained())
        raise_error("we do not handle containing/ed segments (yet): %s", g2e.s.id.c_str());
    if (g2e.d.is_contained())
        raise_error("we do not handle containing/ed segments (yet): %s", g2e.d.id.c_str());

        // Now we are left with dovetailing or blunt

    std::uint32_t v = g2e.needs_flip() ? (d_ix<<1)|d_neg : (s_ix<<1)|s_neg;
    std::uint32_t w = g2e.needs_flip() ? (s_ix<<1)|s_neg : (d_ix<<1)|d_neg;

    arc a1;
    a1.v_lv = (std::uint64_t(v)<<32) | g2e.vtx_l().overhang_l();
    a1.w = w;
    a1.ov = g2e.vtx_l().overlap();
    a1.ow = g2e.vtx_r().overlap();

    const auto it1 = std::upper_bound(arcs.cbegin(), arcs.cend(), a1.v_lv, v_lv_less_u);
    arcs.emplace(it1, a1);

    arc a2;
    a2.v_lv = (std::uint64_t(w^1)<<32) | g2e.vtx_r().overhang_r();
    a2.w = v^1;
    a2.ov = a1.ow;
    a2.ow = a1.ov;

    const auto it2 = std::upper_bound(arcs.cbegin(), arcs.cend(), a2.v_lv, v_lv_less_u);
    arcs.emplace(it2, a2);
}

std::pair<std::vector<arc>::const_iterator, std::vector<arc>::const_iterator>
graph::arcs_from_v_lv(std::uint64_t v_lv) const
{
    std::uint64_t next = ((v_lv>>32)+1)<<32;
    std::vector<arc>::const_iterator lo = std::lower_bound(arcs.cbegin(), arcs.cend(), v_lv, v_lv_less_l);

    return std::make_pair(lo, std::upper_bound(lo, arcs.cend(), next, v_lv_less_u));
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
