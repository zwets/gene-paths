/* targets.cpp
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

#include "targets.h"

#include <vector>
#include <string>
#include <regex>
#include "utils.h"

namespace gfa {

using gene_paths::raise_error;
using gene_paths::verbose_emit;

target
target::parse(const std::string& s)
{
    target tgt;

    std::regex re("([^[:space:]]+)(\\+|-)(:([[:digit:]]+)(:([[:digit:]]+))?)?");
    std::smatch m;

    if (!std::regex_match(s, m, re))
        raise_error("invalid target syntax: %s", s.c_str());

    tgt.ctg = m[1].str();
    tgt.neg = m[2].str().at(0) == '-';

    std::string n;

    n = m[4].str();
    tgt.beg = n.empty() ? std::uint64_t(-1) : std::stoul(n);

    n = m[6].str();
    tgt.end = n.empty() ? std::uint64_t(-1) : std::stoul(n);

    verbose_emit("parsed target: %s%c:%lu:%ld", tgt.ctg.c_str(), tgt.neg ? '-' : '+', tgt.beg, tgt.end);
    return tgt;
}

std::size_t
target::add_seg_to_graph(graph& g, std::string n)
{
    if (g.find_seg_ix(name) != std::uint64_t(-1))
        raise_error("contig with target name already in graph: %s", name.c_str());

    std::size_t ref_seg_ix = g.find_seg_ix(ctg);
    if (ref_seg_ix == std::uint64_t(-1))
        raise_error("contig not in graph: %s", ctg.c_str());

    const seg& ref_seg = g.get_seg(ref_seg_ix);
    if (beg == std::uint64_t(-1))
        { beg = 0; end = ref_seg.len; }
    else if (beg > ref_seg.len)
        raise_error("start pos %lu exceeds segment length %lu for target: %s", beg, ref_seg.len, ctg.c_str());
    else if (end == std::uint64_t(-1))
        end = beg;
    else if (beg > end)
        raise_error("begin position beyond end position on target: %s", ctg.c_str());

    std::stringstream ss;
    ref_seg.write_seq(ss, neg, beg, end);

    name = n;
    g.add_seg({ end-beg, name, ss.str() });
    std::size_t seg_ix = g.get_seg_ix(name);

    verbose_emit("added segment %lu: %s", seg_ix, name.c_str());
    return seg_ix;
}

arc*
target::add_arc_to_graph(graph& g, bool to_tgt) const
{
    if (g.arcs.size() == g.arcs.capacity())
        raise_error("arcs vector exhausted (programmer error)");

    std::size_t tgt_ix = g.get_seg_ix(name);
    std::size_t ref_ix = g.get_seg_ix(ctg);
    const seg& ref_seg = g.get_seg(ref_ix);

    std::uint64_t v, w, lv, lw;

    if (to_tgt) { // from referenced contig to target
        v = neg ? graph::seg_vtx_n(ref_ix) : graph::seg_vtx_p(ref_ix);
        w = neg ? graph::seg_vtx_n(tgt_ix) : graph::seg_vtx_p(tgt_ix);
        lv = neg ? ref_seg.len - end : beg;
        lw = 0;
    }
    else { // from tgt to ctg
        v = neg ? graph::seg_vtx_n(tgt_ix) : graph::seg_vtx_p(tgt_ix);
        w = neg ? graph::seg_vtx_n(ref_ix) : graph::seg_vtx_p(ref_ix);
        lv = end - beg;
        lw = neg ? ref_seg.len - beg : end;
    }

    return reinterpret_cast<arc*>(&*g.add_arc({ graph::v_lv(v, lv), graph::v_lv(w, lw) }));
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
