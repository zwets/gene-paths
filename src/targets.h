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
#include "graph.h"

namespace gfa {


// target - helper structure to capture start and end targets on a graph
//
// In order to find a path between the start and end regions, we add them
// as segments to the graph, but with arcs such that they are traversed
// only at the start and end of the path:
//
//      FROM seg     i------o           TO seg     i------x
//       arc s_o            v          arc e_i     ^
//        contig  ---b======e--->       contig  ---b======e--->
//
// The s_o (start-out) arc goes from the end of the FROM segment to the
// end of the target region on the contig.  Thus if we start the path at
// i, we first cover the target and then leave it (for good) over s_o.
//
// The e_i (end-in) arc is trickier.  We need to get to x (rather than i)
// so need an arc at x.  The solution is this:
//
//       "S" seg     S                 "E" seg            E
//       arc s_i     v                 arc e_o            ^
//      FROM seg     i------o           TO seg     i------x
//       arc s_o            v          arc e_i     ^
//        contig  ---b======e--->       contig  ---b======e--->
//
// We add two dummy zero-length segments that are the "terminals" for the
// the path search.  The path starts with arc s_i and ends with e_o.
//
// In fact we can optimise and do with a single "terminal" of length 1,
// with s_i leaving it at 0, and e_o arriving at 1 (to prevent traversal)
// We refer to this terminal segment as "ter" below.
//
// We can optimise further by noting that when the target has 0 length,
// then we just need the ter can be omitted altogether, and we just need
// the ter:
//
//       ter seg     TER               ter seg          TER
//       arc s_i     v                 arc e_o            ^
//        contig  ---i---------->       contig  ----------o--->
//
// When the graph has only links (non-overlapping tail-to-head edges),
// as is the case with Unicycler, then any path to the end (or from the
// start) of a contig must traverse the whole contig anyway.
//
struct target
{
    enum role_t { START, END };

    // construct a target on graph g
    target(graph& gr)
        : g(gr), ter_arc(NO_ARC), ctg_arc(NO_ARC) { }

    // set the target at ref and give it START or END role
    // the ref must have format "CONTIG[+-][:BEG[:END]]"
    void set(const std::string&, role_t);

    // get the arc that is start/end of the path (depending on role)
    arc get_arc() const { return ter_arc; }

    // get pointer to the arc in graph.arcs (convenience method)
    // NOTE: this invalidates as soon as you add other arcs/targets
    inline const arc* p_arc() const {
        return reinterpret_cast<const arc*>(
                &*(g.arcs_from_v_lv(get_arc().v_lv).first));
    }

#ifdef NDEBUG
    private:   // implementation detail private except when debug/test
#endif
        static constexpr arc NO_ARC = { std::uint64_t(-1), std::uint64_t(-1) };

        graph& g;           // the graph on which target (will) sit
        arc ter_arc;        // the arc between terminal and target
        arc ctg_arc;        // the arc between target and contig
};


} // namespace gfa

#endif // targets_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
