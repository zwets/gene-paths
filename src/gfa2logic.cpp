/* gfa2logic.cpp
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

#include "gfa2logic.h"

#include "utils.h"

namespace gfa2 {

using gene_paths::raise_error;

void
vtx::validate() const
{
    if (l == 0)
        raise_error("segment length is 0 for vertex %s", id.c_str());
    if (b > l || e > l)
        raise_error("begin or end beyond segment length on vertex %s", id.c_str());
    if (b > e)
        raise_error("begin past end on vertex %s", id.c_str());
    if ((p && *(id.end()-1) != '+') || (!p && *(id.end()-1) != '-'))
        raise_error("inconsistent name and orientation: %s defined %s", id.c_str(), (p ? "pos" : "neg"));
}

void
edge::validate() const
{
    s.validate();
    d.validate();
}

} // namespace gfa2

// vim: sts=4:sw=4:ai:si:et
