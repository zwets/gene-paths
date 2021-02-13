/* parsegfa.cpp
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

#include "gfagraph.h"

#include <fstream>
#include <iostream>
#include <string>
#include "gfakluge.hpp"
#include "utils.h"

namespace gene_paths {

static void
gfak_to_graph(gfak::GFAKluge& gfak, gfa::graph& graph)
{
    auto n2s = gfak.get_name_to_seq();
    std::size_t n_segs = n2s.size() << 1;

    verbose_emit("graph has %ld segs", n_segs);
    graph.segs.reserve(n_segs);

    for (auto p : n2s) {
        gfa::seg seg;
        seg.name = p.first;
        seg.len = p.second.length;
        seg.data = p.second.sequence;
        graph.add_seg(seg);
    }

    auto s2e = gfak.get_seq_to_edges();
    std::size_t n_arcs = 0;
    for (auto p : n2s)
        n_arcs += 2 * s2e[p.first].size();

    verbose_emit("graphs has %ld arcs", n_arcs);
    graph.arcs.reserve(n_arcs);

    for (auto p : n2s) {
        auto edges = s2e[p.first];
        for (auto e = edges.cbegin(); e != edges.cend(); ++e) {
            std::string sname = e->source_name + (e->source_orientation_forward ? '+' : '-');
            std::string dname = e->sink_name + (e->sink_orientation_forward ? '+' : '-');
            graph.add_edge(
                    sname, e->source_begin, e->source_end,
                    dname, e->sink_begin, e->sink_end);
        }
    }
}

static void
add_fasta_to_gfak(gfak::GFAKluge& gfak, std::istream& fasta)
{
    std::string line;
    std::string seqid;
    std::string data;

    while (line.empty() && std::getline(fasta, line))
        /* be lenient about empty lines at start */;

    while (!line.empty()) {

        if (line[0] != '>')
            raise_error("invalid FASTA header: %s", line.c_str());

        data.clear();

        std::string::const_iterator p = line.cbegin() + 1;
        while (p != line.cend() && !std::isspace(*p)) ++p;
        seqid = std::string(line, 1, p - line.cbegin() - 1);

        while (std::getline(fasta, line)) {
            if (line.empty())
                continue;
            if (line[0] == '>')
                break;
            data.append(line);
            line.clear();
        }

        gfak.set_sequence_data(seqid, data);
        //gfak::sequence_elem *psel = gfak.find_sequence_elem(seqid);
        //if (psel)
        //    psel->sequence = data;
    }
}

gfa::graph
parse_gfa(std::istream& file)
{
    gfak::GFAKluge gfak;
    gfa::graph graph;

    if (!gfak.parse_gfa_file(file))
        raise_error("failed to parse GFA");

    gfak_to_graph(gfak, graph);

    return graph;
}

gfa::graph
parse_gfa(std::istream& gfa, std::istream& fasta)
{
    gfak::GFAKluge gfak;
    gfa::graph graph;

    if (!gfak.parse_gfa_file(gfa))
        raise_error("failed to parse GFA");

    add_fasta_to_gfak(gfak, fasta);
    gfak_to_graph(gfak, graph);

    return graph;
}

} // namespace gene_paths

// vim: sts=4:sw=4:ai:si:et
