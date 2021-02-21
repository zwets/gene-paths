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
 * Segments store the data for the + orientation.
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
 * Though this H-shape should not occur in final assembly graphs, we can
 * cater for it by slightly extending Heng Li's gfatools data model. That
 * model assumes that edges dovetail or do not overlap at all:
 *
 *      |<--- lv --->|<-- ov -->|
 *   v: -------------============
 *                     overlap
 *                w: ============--------------
 *                   |<-- ow -->|<---- lw ---->
 *
 * The gfatools model encodes the graph d as an array of directed arcs,
 * two per edge.  One leaves v at length lv, its complement leaves w' at
 * lw.  The arcs always arrive at 0 (on +) or $ (on -).
 *
 * Our extension is to additionally encode the lengths rv and rw, where
 * l{v,w} are upstream of the overlap, and r{v,w} downstream:
 *
 *      |<--- lv --->|<-- ov -->|<- rv ->|
 *   v: -------------============---------
 *                     overlap
 *       w: ---------============---------------
 *          |<- lw ->|<-- ow -->|<---- rw ---->|
 *
 * We then encode each edge as eight arcs, corresponding to going from
 * (1) v to w at (a) lv, (b) lv+ov, (2) w to v at (a) lw, (b) lw+ow, and
 * their four complements.  Iff ov=ow=0, there are no (b) arcs.
 *
 * Note that v and w above may be either + or -, as lv and lw are always
 * upstream of the overlap.  They are measured from 0 for a + segment,
 * and from $ for a - segment.
 *
 * As in gfatools, arcs are stored in an array sorted on v_lv, which is
 * vtx_ix<<32|lv, so that the outbound arcs from every vertex vtx_ix are
 * contiguous and sorted on how "early" they leave the vertex.
 *
 * -- Note 1: Dovetailing is special case rv=0 && lw=0 (or lv=0 && rw=0)
 *
 * It would seem that in this case there is no point in storing w_0→v_lv
 * and v'_0→w'_rw, but note that adding these arcs adds a path to the
 * overlap on the other segment from some other edge Z linked to w_0:
 *
 *   v: ------------->>>>>>>>>>>>
 *                w: ^===========---------------
 *        z: >>>>>>>>^
 *
 * -- Note 2: Non-overlap additionally has lv=$ && rw=$ (or rv=$ && lw=$)
 *
 *      |<--- lv --->|
 *   v: --------------              (or v and w reversed)
 *                 w: ----------
 *                    |<- rw ->|
 *
 * Again, we do not eliminate w_0→v_lv or its complement v'_0→w'_rw, as
 * other arcs may continue (niche case ...) the path to other segments:
 *
 *                u: >>>>>>>>>>>>
 *   v: -------------^
 *                w: ^-------
 *
 * We can however leave out v_(lv+ov)→w_(lw+ow) and w_(lw+ow)→v_(lv+ov),
 * and their complements, as ov=ow=0 so they are equivalent to the first
 * arcs.
 *
 * -- Note 3: Zero-overlap can theoretically happen in middle too
 *
 * But we treat it the same as at the end, including the ov=ow=0 case.
 *
 * -- Note 4: Containment can now be expressed too
 *
 * And this will be convenient for defining paths, as the start and end
 * regions can be defined as containments:
 *
 *      |<--- lv --->|<-- ov -->|<- rv ->|
 *   v: -------------============---------
 *                w: ============             (lw=rw=0)
 */

struct seg {
    std::uint64_t len;
    std::string name;
    std::string data;

    // writes the sequence content in [beg,end) to os, optionally reverse complementing
    // note that beg and end are forward positions as in GFA2, i.e. applied before rc
    std::ostream& write_seq(std::ostream& os, bool rc = false, std::uint32_t beg = 0, std::uint32_t end = std::uint32_t(-1)) const;

    // write the sequence content in [beg,end) on the pos or neg vertex of the segment
    // where beg and end are interpreted on the vertex, i.e. from end when neg is true
    std::ostream& write_vtx(std::ostream& os, bool neg, std::uint32_t beg = 0, std::uint32_t end = std::uint32_t(-1)) const {
        return neg ? write_seq(os, true, end == std::uint32_t(-1) ? 0 : len-end, len-beg) : write_seq(os, false, beg, end);
    }
};

struct arc {
    std::uint64_t v_lv;     // vtx_ix<<32|lv packed for sorting
    std::uint64_t w_lw;     // vtx_ix<<32|lw

        // convenience selectors

    inline std::uint64_t v() const { return v_lv>>32; }
    inline std::uint64_t lv() const { return v_lv & 0xFFFFFFFFL; }
    inline std::uint64_t w() const { return w_lw>>32; }
    inline std::uint64_t lw() const { return w_lw & 0xFFFFFFFFL; }
};

struct graph {

        // building the graph

    void add_seg(const seg&);

    void add_edge(const std::string& sref, std::uint32_t sbeg, std::uint32_t send,
                  const std::string& dref, std::uint32_t dbeg, std::uint32_t dend);

    std::vector<arc>::iterator add_arc(const arc&);

        // segment storage and lookup

    std::vector<seg> segs;
    std::map<std::string, std::size_t> seg_ixs;

    std::size_t find_seg_ix(const std::string& name) const;     // ix or size_t(-1)
    std::size_t get_seg_ix(const std::string& name) const;      // ix or error out

    const seg* find_seg(const std::string& name) const;         // pointer or null
    inline const seg& get_seg(const std::string& name) const    // ref or error out
        { return segs.at(get_seg_ix(name)); }
    inline const seg& get_seg(const std::size_t seg_ix) const   // ref or exception
        { return segs.at(seg_ix); }

        // segment to vertex mapping

    inline static std::uint64_t seg_vtx_p(std::size_t seg_ix) { return seg_ix<<1; }
    inline static std::uint64_t seg_vtx_n(std::size_t seg_ix) { return seg_ix<<1|1; }
    inline static std::uint64_t inv_vtx(std::uint64_t vtx_ix) { return vtx_ix^1; }
    inline static std::uint64_t vtx_seg(std::uint64_t vtx_ix) { return vtx_ix>>1; }

        // query vertex properties

    inline static bool is_pos(std::uint64_t v) { return !(v & 1); }
    inline static bool is_neg(std::uint64_t v) { return v & 1; }

        // convenience functions on v_lv 

    inline static std::uint64_t v_lv(std::uint64_t v, std::uint64_t lv) { return v<<32|lv; }
    inline static std::uint64_t vlv_v(std::uint64_t v_lv) { return v_lv>>32; }
    inline static std::uint64_t vlv_lv(std::uint64_t v_lv) { return v_lv & 0xFFFFFFFFL; }
    inline static std::uint64_t vlv_seg(std::uint64_t v_lv) { return v_lv>>33; }

        // arc storage and lookup

    std::vector<arc> arcs;

    // begin and past-the-end iterator for all arcs leaving v at lv or further downstream
    std::pair<std::vector<arc>::const_iterator, std::vector<arc>::const_iterator>
        arcs_from_v_lv(std::uint64_t) const;

    // begin and past-the-end iterator for all arcs leaving vtx
    inline std::pair<std::vector<arc>::const_iterator, std::vector<arc>::const_iterator>
        arcs_from_vtx(std::uint64_t vtx_ix) const
        { return arcs_from_v_lv(vtx_ix<<32); }
};


} // namespace gfa

#endif // gfagraph_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
