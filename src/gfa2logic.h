/* gfa2logic.h
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
#ifndef gfa2logic_h_INCLUDED
#define gfa2logic_h_INCLUDED

#include <string>

namespace gfa2 {

/* Data structures from GFA2, used while parsing and converting to arcs. */

struct vtx {

    const std::string& id;  // seg name, including + or -
    std::uint32_t l;        // length of the segment
    std::uint32_t b;        // begin overlap
    std::uint32_t e;        // end overlap
    bool p;                 // is positive

    // validates consistency of the fields, errors out if invalid
    void validate() const;

        // The following methods return the indicated lengths,
        // after taking into account orientation.
        // So if !p, they are measured from the end.
        //
        //      <--- l1 ---><----- r1 ------>
        //      <------ l2 -------><-- r2 --> 
        //                  <- o ->
        //   v: ------------=======----------
        //                  <- o ->
        //      <------ r1i ------><-- l1i -->
        //      <-- r2i ---><------ l2i ----->

    inline std::uint32_t o() const { return e - b; }

    inline std::uint32_t l1() const { return b; }
    inline std::uint32_t l2() const { return e; }
    inline std::uint32_t r1() const { return l-b; }
    inline std::uint32_t r2() const { return l-e; }

        // same for inverse of the vertex

    inline std::uint32_t l1i() const { return l-e; }
    inline std::uint32_t l2i() const { return l-b; }
    inline std::uint32_t r1i() const { return e; }
    inline std::uint32_t r2i() const { return b; }

};

struct edge {

    vtx v;   // first vertex
    vtx w;   // second vertex

    // validates consistency of the edge and its vertices
    inline void validate() const { v.validate(); w.validate(); }

        // Return the lv and lw lengths as discussed in graph.h,
        // in the orientation of the segments.  lv2 and lw2 point
        // at the end of the overlap.  The ~i methods return the
        // parameters for the inverse orientation of the segment.

    inline std::uint32_t ov() const   { return v.o(); }
    inline std::uint32_t lv() const   { return v.l1(); }
    inline std::uint32_t lv2() const  { return v.l2(); }
    inline std::uint32_t lvi() const  { return v.l1i(); }
    inline std::uint32_t lv2i() const { return v.l2i(); }
    inline std::uint32_t rv() const   { return v.r1(); }
    inline std::uint32_t rv2() const  { return v.r2(); }
    inline std::uint32_t rvi() const  { return v.r1i(); }
    inline std::uint32_t rv2i() const { return v.r2i(); }

    inline std::uint32_t ow() const   { return w.o(); }
    inline std::uint32_t lw() const   { return w.l1(); }
    inline std::uint32_t lw2() const  { return w.l2(); }
    inline std::uint32_t lwi() const  { return w.l1i(); }
    inline std::uint32_t lw2i() const { return w.l2i(); }
    inline std::uint32_t rw() const   { return w.r1(); }
    inline std::uint32_t rw2() const  { return w.r2(); }
    inline std::uint32_t rwi() const  { return w.r1i(); }
    inline std::uint32_t rw2i() const { return w.r2i(); }

};

} // namespace gfa2

#endif // gfa2logic_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
