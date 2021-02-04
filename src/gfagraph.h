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
 * A segment is a sequence with length and data.  A vertex is one side of
 * the segment, and corresponds to the + or - orientation of the segment.
 * The sequence data in the + and - orientations are reverse complements.
 * We store the data for the + orientation.
 *
 * Segments and vertices are identified by indices.  The two vertices of
 * segment seg_ix are given by seg_ix<<1|ori, thus seg_ix = vtx_ix>>1.
 *
 * GFA2 allows arbitrary overlaps between segments in an edge, allowing:
 *
 *      \         / s1
 *       \_______/       Segments s1 and s2 overlap in the middle
 *       /       \
 *      /         \ s2
 *
 * We do not cater for this type of overlap, nor its special case of
 * containment where the length of s1 or s2 is equal to the overlap.  
 * Overlaps must be at the ends of the segments (dovetailing), or be
 * absent altogether (typical in a final assembly graph).
 *
 * This means we can use Heng Li's arc model from gfatools.  An arc is a
 * directed edge between two vertices v and w:
 *
 *      |<--- lv --->|<-- ov -->|
 *   v: -------------============
 *                     overlap         
 *                w: ============------------->
 *                   |<-- ow -->|<---- lw ---->
 *
 * Values lv are lw are lengths that contribute to the sort order of the
 * forward and reverse arcs respectively.  For each edge we store both
 * the forward and reverse arc.
 *
 * Arcs are stored in a vector that is sorted on v_vl = vtx_id<<32|vl.
 * Hence the outbound arcs from any v are contiguous and sorted on how
 * "early" they leaves the vertex.
 */

struct seg {
    std::uint64_t len;
    std::string name;
    std::string data;
};

struct vtx {
};

struct arc {
    std::uint64_t v_lv;     // 32 bits vtx_ix = seg_ix<<1|ori, lower 32 bits lv
    std::uint32_t w;        // 32 bits vtx_ix = seg_ix<<1|ori
    std::uint32_t ov, ow;
    std::uint32_t dummy;    // reserved (align struct to 64 bits)
};

struct graph {

    std::vector<seg> segs;
    std::map<std::string, std::size_t> seg_ixs;

    std::vector<arc> arcs;

    std::pair<std::vector<arc>::const_iterator, std::vector<arc>::const_iterator> 
        arcs_from_vtx(std::uint64_t) const;

    std::pair<std::vector<arc>::const_iterator, std::vector<arc>::const_iterator> 
        arcs_from_v_lv(std::uint64_t) const;

    void add_seg(const seg&);

    void add_edge(const std::string& sref, std::uint32_t sbeg, std::uint32_t send,
                  const std::string& dref, std::uint32_t dbeg, std::uint32_t dend);
};

} // namespace gfa

#endif // gfagraph_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
