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

void
target::set(const std::string& ref, role_t role)
{
        // parse the reference

    std::regex re("([^:[:space:]]+)(:([[:digit:]]+)(:([[:digit:]]+))?)?(\\+|-)");
    std::smatch m;

    if (!std::regex_match(ref, m, re))
        raise_error("invalid target syntax: %s", ref.c_str());

    std::string ctg = m[1].str();

    std::string _b = m[3].str();
    std::size_t beg = _b.empty() ? std::uint64_t(-1) : std::stoul(_b);

    std::string _e = m[5].str();
    std::size_t end = _e.empty() ? beg : std::stoul(_e);

    bool neg = m[6].str().at(0) == '-';

    verbose_emit("parsed target: %s:%ld:%ld%c", ctg.c_str(), beg, end, neg ? '-' : '+');

        // locate or create the terminator

    static const std::string TER("__T__");

    std::size_t ter_ix = g.find_seg_ix(TER);

    if (ter_ix == std::uint64_t(-1)) {
        g.add_seg( { 1, TER, "X" } );
        ter_ix = g.get_seg_ix(TER);

        verbose_emit("added terminal segment %lu: %s", ter_ix, TER.c_str());
    }
    else {
        verbose_emit("terminal segment %lu: %s", ter_ix, TER.c_str());
    }

        // locate the referenced contig in graph

    std::size_t ref_ix = g.find_seg_ix(ctg);
    if (ref_ix == std::uint64_t(-1))
        raise_error("contig not in graph: %s", ctg.c_str());

    const seg& ref_seg = g.get_seg(ref_ix);
    if (beg == std::uint64_t(-1))
        beg = end = role == START ? 0 : ref_seg.len;
    else if (beg > ref_seg.len)
        raise_error("start pos %lu exceeds segment length %lu for target: %s", beg, ref_seg.len, ctg.c_str());
    else if (beg > end)
        raise_error("begin position beyond end position on target: %s", ctg.c_str());

        // locate or create the target segment

    std::uint64_t seg_ix = std::uint64_t(-1);

    if (beg != end) { // see targets.h, no need if ref is 0 length

        std::stringstream ss;
        ss << ctg << ':' << beg << ':' << end;
        std::string seg_name = ss.str();

        seg_ix= g.find_seg_ix(seg_name);
        if (seg_ix == std::uint64_t(-1)) {

            std::stringstream ss;
            ref_seg.write_seq(ss, false, beg, end);  // Store the + segment always

            g.add_seg({ end-beg, seg_name, ss.str() });
            seg_ix = g.get_seg_ix(seg_name);

            verbose_emit("added target segment %lu: %s", seg_ix, seg_name.c_str());
        }
        else {
            verbose_emit("found target segment %lu: %s", seg_ix, seg_name.c_str());
        }
    }
    else {
        seg_ix = ref_ix;
        verbose_emit("target segment is contig %lu: %s", seg_ix, ctg.c_str());
    }

        // remove existing arcs

    if (ter_arc.v_lv != std::uint64_t(-1)) {
        g.arcs.erase(g.arcs_from_v_lv(ter_arc.v_lv).first);
    }

    if (ctg_arc.v_lv != std::uint64_t(-1)) {
        g.arcs.erase(g.arcs_from_v_lv(ctg_arc.v_lv).first);
    }

        // create the new ctg_arc

    std::uint64_t v = 0L, w = 0L, lv = 0L, lw = 0L;

    if (seg_ix != ref_ix) {

        if (g.arcs.size() == g.arcs.capacity())
            raise_error("programmer error: arcs vector exhausted (cap %lu)", g.arcs.size());

        if (role == START) { // from seg_$ to ctg_end
            v = graph::seg_vtx(seg_ix, neg);
            w = graph::seg_vtx(ref_ix, neg);
            lv = end - beg;
            lw = end;
        }
        else if (role == END) { // from ctg_beg to seg_0
            v = graph::seg_vtx(ref_ix, neg);
            w = graph::seg_vtx(seg_ix, neg);
            lv = beg;
            lw = 0;
        }
        else {
            raise_error("programmer error: target role not START or END");
        }

        ctg_arc = { graph::v_lv(v, lv), graph::v_lv(w, lw) };
        g.add_arc(ctg_arc);

        verbose_emit("added ctg arc %lu: %lu_%lu to %lu_%lu", g.arcs.size(), v, lv, w, lw);
    }
    else {
        ctg_arc = NO_ARC;
    }

        // create the new ter_arc

    if (g.arcs.size() == g.arcs.capacity())
        raise_error("programmer error: arcs vector exhausted (cap %lu)", g.arcs.size());

    if (role == START) { // from ter_0 to seg_0 or ctg_b
        v = graph::seg_vtx(ter_ix, false); // pos
        w = graph::seg_vtx(seg_ix, neg);
        lv = 0;
        lw = seg_ix == ref_ix ? end : 0;
    }
    else if (role == END) { // from seg_$ or ctg_e to ter_1
        v = graph::seg_vtx(seg_ix, neg);
        w = graph::seg_vtx(ter_ix, false); // pos
        lv = seg_ix == ref_ix ? beg : end - beg;
        lw = 1;
    }
    else {
        raise_error("programmer error: target role not START or END");
    }

    ter_arc = { graph::v_lv(v, lv), graph::v_lv(w, lw) };
    g.add_arc(ter_arc);

    verbose_emit("added terminal arc %lu: %lu_%lu to %lu_%lu", g.arcs.size(), v, lv, w, lw);
}


} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
