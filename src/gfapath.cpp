/* gfapath.cpp
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

#include "gfapath.h"

namespace gfa {

void
path::add_arc(const graph& g, std::size_t arc_ix)
{
    const arc& a = g.arcs.at(arc_ix);

    arc_ixs.push_back(arc_ix);
    len += a.v_lv - now_at + a.ow;
    now_at = std::uint64_t(a.w)<<32 | a.ow;
}

} // namespace gfa

// vim: sts=4:sw=4:ai:si:et
