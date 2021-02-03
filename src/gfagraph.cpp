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
#include <map>
#include "gfa2logic.h"
#include "utils.h"

namespace gfa {

using gene_paths::raise_error;

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

    auto ps = sref.end() - 1;
    auto pd = dref.end() - 1;

    if (*ps != '-' && *ps != '+')
        raise_error("sequence reference without sign: %s", sref.c_str());

    if (*pd != '-' && *pd != '+')
        raise_error("sequence reference without sign: %s", dref.c_str());

    std::uint64_t s_neg = *ps == '-' ? 1 : 0;
    std::uint64_t d_neg = *pd == '-' ? 1 : 0;

        // look up segments

    std::string s_name(sref.begin(), ps);
    std::string d_name(dref.begin(), pd);

    const auto s_iter = seg_ixs.find(s_name);
    if (s_iter == seg_ixs.end())
        raise_error("unknown sequence in edge: %s", s_name.c_str());

    const auto d_iter = seg_ixs.find(d_name);
    if (d_iter == seg_ixs.end())
        raise_error("unknown sequence in edge: %s", d_name.c_str());

    std::size_t s_ix = s_iter->second;
    std::size_t d_ix = d_iter->second;

    const seg& s_seg = segs[s_ix];
    const seg& d_seg = segs[d_ix];

        // create the GFA2 representation with the logic

    gfa2::edge gfa2edge = {
        { sref, std::uint32_t(s_seg.len), sbeg, send, bool(!s_neg) },
        { dref, std::uint32_t(d_seg.len), dbeg, dend, bool(!d_neg) }
    };

    gfa2edge.validate();
}

} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
