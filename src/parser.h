/* parser.h
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
#ifndef parser_h_INCLUDED
#define parser_h_INCLUDED

#include <iostream>
#include "graph.h"

namespace gfa {

// parse a GFA file with embedded sequences into a gfa::graph
// reserving extra room for spare_segs and spare_arcs
extern graph parse(std::istream& gfa, int spare_segs = 0, int spare_arcs = 0);

// parse a GFA file with sequences in a FASTA file into a gfa::graph
// reserving extra room for spare_segs and spare_arcs
extern graph parse(std::istream& gfa, std::istream& fna, int spare_segs = 0, int spare_arcs = 0);

} // namespace gfa

#endif // parser_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
