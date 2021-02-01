/* gfagraph.h
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
#ifndef gfagraph_h_INCLUDED
#define gfagraph_h_INCLUDED

#include <string>
#include <vector>
#include <map>

namespace gfa {

/* Data structures inspired by https://github.com/lh3/gfatools.
 *
 * A segment is a sequence with length and data.  A vertex is one side
 * of a segment.  All objects are identified by their index, and vtx_id
 * is seg_ix<<1|orientation, thus seg_ix = vtx_ix>>2.
 *
 * An arc is a directed edge between two vertices, which in GFA2 is like:
 *
 *      |<--- lv --->|<-- ov -->|<- xv ->
 *   v: -------------=========== - - - ->
 *         
 *        w: - - - - ===========-------------->
 *           <- xw ->|<-- ow -->|<---- lw ---->
 *
 * An edge connects a pair of vertices in the GFA2 way, i.e. with an
 * overlap [beg,end] on each vertex.  We ignore the cigar for now.
 */

struct seg {
    std::uint64_t len;
    std::string name;
    std::string data;
};

struct vtx {
};

struct arc {
    std::uint64_t v_lv;     // 32 bits index, lower 32 bits lv
    std::uint64_t w_lw;     // 32 bits index, lower 32 bits lw
    std::uint32_t ov, ow;
};

struct graph {
    std::vector<seg> segs;
    std::vector<arc> arcs;

    std::map<std::string, std::size_t> seg_ixs;

    void add_seg(const seg&);

    void add_edge(const std::string& sref, std::uint32_t sbeg, std::uint32_t send,
                  const std::string& dref, std::uint32_t dbeg, std::uint32_t dend);
};

} // namespace gfa

#endif // gfagraph_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
