/* gfapath.h
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
#ifndef gfapath_h_INCLUDED
#define gfapath_h_INCLUDED

#include <vector>
#include "gfagraph.h"

namespace gfa {

struct path {
    std::uint64_t len;
    std::uint64_t start;
    std::uint64_t now_at;
    std::vector<std::size_t> arc_ixs;

    path(std::uint32_t vtx, std::uint32_t pos) 
	    : len(0), start(std::uint64_t(vtx)<<32|pos), now_at(start) { }

    void add_arc(const graph& g, std::size_t arc_ix);
};

} // namespace gfa

#endif // gfapath_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
