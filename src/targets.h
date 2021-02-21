/* targets.h
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
#ifndef targets_h_INCLUDED
#define targets_h_INCLUDED

#include <string>
#include "gfagraph.h"

namespace gfa {


// target - helper structure capturing a target on a graph
//
struct target
{
    std::string name;   // name of target seg in graph
    std::string ctg;    // name of referent contig in graph
    bool neg;           // true iff target is on neg vertex
    std::uint64_t beg;  // begin of target on contig
    std::uint64_t end;  // end of target on contig

    // parses a "CONTIG[+-][:BEG[:END]]" referent into target
    static target parse(const std::string&);

    // add target segment to g with name n, returning its seg_ix
    std::size_t add_seg_to_graph(graph& g, std::string n);

    // add arc from ctg to target (default) or vice versa
    // returns pointer to the inserted arc (for path start/end)
    arc* add_arc_to_graph(graph& g, bool to_tgt=true) const;
};


} // namespace gfa

#endif // targets_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
