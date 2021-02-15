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

/* Data structures from GFA2, used while parsing. */

struct vtx {
    const std::string& id;  // seg name, including + or -
    std::uint32_t l;        // length of the segment
    std::uint32_t b;        // begin overlap
    std::uint32_t e;        // end overlap
    bool p;                 // positive (false is neg)

    void validate() const;

    inline bool is_contained() const { return b == 0 && e == l; }
    inline bool is_container() const { return b != 0 && e != l; }
    inline bool is_blunt_r() const { return p ? b == l : e == 0; }
    inline bool is_blunt_l() const { return p ? e == 0 : b == l; }
    inline bool dovetails_r() const { return p ? e == l && b < l && b > 0 : b == 0 && e > 0 && e < l; }
    inline bool dovetails_l() const { return p ? b == 0 && e > 0 && e < l : e == l && b < l && b > 0; }
    inline std::uint32_t overlap() const { return e - b; }
    inline std::uint32_t overhang_l() const { return p ? b : l - e; }
    inline std::uint32_t overhang_r() const { return p ? l - e : b; }
};

struct edge {
    vtx v;   // source vertex
    vtx w;   // dest vertex vertex

    void validate() const;

    inline std::uint32_t lv() const { return v.overhang_l(); }
    inline std::uint32_t ov() const { return v.overlap(); }
    inline std::uint32_t rv() const { return v.overhang_r(); }

    inline std::uint32_t lw() const { return w.overhang_l(); }
    inline std::uint32_t ow() const { return w.overlap(); }
    inline std::uint32_t rw() const { return w.overhang_r(); }
};

} // namespace gfa2

#endif // gfa2logic_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
