/* graph.h
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
#ifndef graph_h_INCLUDED
#define graph_h_INCLUDED

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
 * model assumes that edges abut or dovetail:
 *
 *      |<--- lv --->|<-- ov -->|
 *   v: -------------============
 *                     overlap
 *                w: ============--------------
 *                   |<-- ow -->|<---- lw ---->
 *
 * The gfatools model encodes the graph d as an array of directed arcs,
 * two per edge.  One leaves v at length lv, its complement leaves w' at
 * lw.  The arcs always arrive at the head of the other segment.
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
 * To completely implement this model, each edge requires eight arcs:
 * - v:lv → w:lw, v:lv+ov → w:lw+ow, w:lw → v:lv, w:lw+ow → v:lv+ov
 * and their four complements:
 * - v':rv → w':rw, v':rv+ov → w':rw+ow, w':rw → v':rv, w':rw+ow → v':rv+ov
 *
 * As in gfatools, arcs are stored in an array sorted on v_lv, which is
 * vtx_ix<<32|lv, so that the outbound arcs from every vertex vtx_ix are
 * contiguous and sorted on how "early" they leave the vertex.
 *
 * -- Special case 1: ov = ow = 0 (X-shape)
 *
 * If there is no overlap on either vertex (ov=ow=0), then half the arcs
 * coincide, so we are left with:
 *
 *   v:lv → w:lw, w:lw → v:lv, v':rv → w':rw, w':rw → v':rv
 *
 * Note that we don't need all four if the "overlap" is at the start/end
 * of the vertices (connecting them tail-to-head), as shown below.
 *
 * -- Special case 2: dovetailing (rv = lw = 0 || lv = rw = 0)
 *
 * This is as in the gfatools model above, except that in GFA2 arcs can
 * "land" anywhere on the destination.  So, in addition to v:lv → w:0 and
 * w':rw → v':0 , we also have v:lv+ov → w:ow and w':rw+ow → v':ov.
 *
 *      |<--- lv --->|<-- ov -->|
 *   v: -------------============
 *                   ^          v 
 *                w: ============--------------
 *                   |<-- ow -->|<---- rw ---->
 *
 * Logically we also have their inverses, but they are odd: w:0 → v:lv
 * takes arrivals on w straight off to ov, and v':0 → w':rw does this
 * for arrivals on v': they never traverse the segment they came for.
 *
 * Also awkward are w:ow → v:$ and v':ov → w':$.  They take arrivals
 * over the overlap, but then to the end of the other vertex (without
 * including any of its bases), then hop onto any arc.
 *
 * -- Special case 3 = 1+2: abutment (ov=ow=0 && (rv=lw=0 || lv=rw=0))
 *
 * Combining special cases 1 & 2, we have non-overlapping links, with
 * arcs v:$ → w:0 and its complement w':$ → v':$, and their inverses:
 *
 *      |<--- lv --->|
 *   v: --------------              (or v and w reversed)
 *                 w: ----------
 *                    |<- rw ->|
 *
 * However the inverses (w:0 → v:$ and v':0 →w':$) are even weirder, as
 * they take arrivals on either vertex straight off the end of the other,
 * serving stricly as "transit station".
 *
 * == What to do with the "niche cases"?
 *
 * The links in SPAdes GFA output are mostly #3 or #2.  Unicycler seems
 * to produce only #3: abutting links with no (or little) overlap.
 * Looking at their typical L lines in GFA (simplified for clarity):
 *
 *     L * s1 s2 *
 *     L * s1 s3 *
 *     L * s4 s2 *
 *
 * Then should we make path s4-s3?  It is generated by s4→s2:0→s1:$→s3,
 * whereas no base on s2 or s1 is actually involved.  If it were real,
 * then it would have its own "L s4 s3". 
 *
 * Adding "L s3 s2" to the above implies the path s3-s3, which certainly
 * doesn't look like it was intended, yet by the letter it is allowed.
 *
 * CONCLUSION: some subtle semantic context is lost after translating
 * GFA1 L nodes to GFA2's generic E nodes.  Safer to assume that the
 * "strange arcs" are unintentional.  This minimises the number of
 * possibly "false positive" paths.
 *
 * WHAT THIS MEANS: we do not create arcs that terminate at $ or start
 * at 0.
 */

struct seg {
    std::uint64_t len;
    std::string name;
    std::string data;

    // writes the sequence content in [beg,end) to os, optionally reverse complementing
    // note that beg and end are positions as in GFA2, i.e. before orienting the segment
    std::ostream& write_seq(std::ostream& os, bool rc = false, std::uint32_t beg = 0, std::uint32_t end = std::uint32_t(-1)) const;

    // write the sequence content in [beg,end) on the pos or neg vertex of the segment,
    // here beg and end are interpreted on the vertex, that is after reversing
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

    inline static std::uint64_t seg_vtx(std::size_t seg_ix, bool neg) { return seg_ix<<1|!!neg; }
    inline static std::uint64_t vtx_seg(std::uint64_t vtx_ix) { return vtx_ix>>1; }
    inline static std::uint64_t inv_vtx(std::uint64_t vtx_ix) { return vtx_ix^1; }

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

#endif // graph_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
