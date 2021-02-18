/* gene-paths.cpp
 *
 * gene-paths - determine gene order and orientation in assemblies
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

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <cstring>
#include <cstdlib>
#include "parsegfa.h"
#include "gfagraph.h"
#include "utils.h"


using namespace gene_paths;


static const std::string USAGE(
"Usage: gene-paths [OPTIONS] GFA_FILE FROM TO\n"
"\n"
"  Output the shortest path between regions FROM and TO in the assembly\n"
"  graph in GFA_FILE.\n"
"\n"
"  FROM and TO are specified as SEG[+-][:BEG[:END]]], where SEG is the\n"
"  segment identifier, + or - the STRAND, and BEG and END the start and\n"
"  end positions of the target region on SEG.\n"
"\n"
"  When END is omitted, BEG specifies a zero-length position.  When BEG\n"
"  and END are both omitted they default to the start and end of SEG.\n"
"\n"
"  OPTIONS\n"
"   -f, --fasta FILE  read sequences for GFA_FILE from FILE\n"
"   -u, --undirected  also search for TO upstream of FROM\n"
"   -v, --verbose     write detailed information to stderr\n"
"   -h, --help        output this information and exit\n"
"\n"
"  The default is to search for a path that has FROM upstream of TO.\n"
"  Option -u/--undirected searches in both directions.\n"
"\n"
"  STRAND and POSITION are interpreted as in GFA2:\n"
"  - We call the sequence data in GFA_FILE (or FASTA FILE when external)\n"
"    the + strand.  The - strand is its reverse complement.\n"
"  - Positions are inbetween bases, starting at 0 left of the sequence,\n"
"    and ending at L at its right, where L is the sequence length.\n"
"  - Positions are interpreted on the segment BEFORE it is oriented, so\n"
"    e.g. C-:0:5 refers to the final 5 bases on C's reverse complement.\n"
"\n");

static void usage_exit()
{
    std::cerr << USAGE;
    std::exit(1);
}

// parses the 
static void parse_tgt_ref(std::string s, std::string* seg, bool* neg, std::uint64_t* beg, std::uint64_t* end)
{
    std::regex re("([^[:space:]]+)(\\+|-)(:([[:digit:]]+)(:([[:digit:]]+))?)?");
    std::smatch m;
    std::string n;
    
    if (!std::regex_match(s, m, re))
        raise_error("invalid target syntax (see --help): %s", s.c_str());

    *seg = m[1].str();
    *neg = m[2].str().at(0) == '-';
    if (!(n = m[4].str()).empty())
        *beg = std::stoul(n);
    if (!(n = m[6].str()).empty())
        *end = std::stoul(n);

    verbose_emit("parsed target: %s%c:%lu:%ld", seg->c_str(), *neg ? '-' : '+', *beg, *end);
}

static std::size_t add_target(gfa::graph& g, std::string ref, std::string name, int dir /* -1 from, 0 both, 1 to */)
{
    std::string ctg;
    bool neg = false;
    std::uint64_t beg = std::uint64_t(-1);
    std::uint64_t end = std::uint64_t(-1);

        // add segment to the graph

    parse_tgt_ref(ref, &ctg, &neg, &beg, &end);

    std::size_t ref_seg_ix = g.find_seg_ix(ctg);
    if (ref_seg_ix == std::uint64_t(-1))
        raise_error("contig not in graph: %s", ctg.c_str());

    const gfa::seg& ref_seg = g.get_seg(ref_seg_ix);
    if (beg == std::uint64_t(-1))
        { beg = 0; end = ref_seg.len; }
    else if (beg > ref_seg.len)
        raise_error("start pos %lu exceeds segment length %lu for target: %s", beg, ref_seg.len, ref.c_str());
    else if (end == std::uint64_t(-1))
        end = beg;
    else if (beg > end)
        raise_error("begin position beyond end position on target: %s", ref.c_str());

    std::stringstream ss;
    ref_seg.write_seq(ss, neg, beg, end);

    g.add_seg( { end - beg, name, ss.str() });
    std::size_t seg_ix = g.get_seg_ix(name);

    verbose_emit("added segment %lu: %s", seg_ix, name.c_str());

        // add arcs to the graph

    if (dir != 1) { // arcs when FROM or BOTH
        // TODO
    }

    if (dir != -1) { // arcs when TO or BOTH
        // TODO
    }

    return seg_ix;
}

static void add_targets(gfa::graph& g, std::string from_ref, std::string to_ref, bool undirected)
{
    std::size_t seg_ix0 = add_target(g, from_ref, "__TGT0__", undirected ? 0 : -1);
    std::size_t seg_ix1 = add_target(g, to_ref, "__TGT1__", undirected ? 0 : 1);
}

int main (int /*argc*/, char *argv[])
{
    set_progname("gene-paths");

    std::string gfa_fname;
    std::string fna_fname;
    std::string from_tgt;
    std::string to_tgt;
    bool undirected = false;

        // parse options

    while (*++argv && **argv == '-')
    {
        if (!std::strcmp("-v", *argv) || !std::strcmp("--verbose", *argv)) {
            set_verbose(true);
        }
        else if (!std::strcmp("-h", *argv) || !std::strcmp("--help", *argv)) {
            std::cout << USAGE;
            return 0;
        }
        else if (!std::strcmp("-u", *argv) || !std::strcmp("--undirected", *argv)) {
            undirected = true;
        }
        else if ((!std::strcmp("-f", *argv) || !std::strcmp("--fasta", *argv)) && *++argv) {
            fna_fname = *argv;
        }
        else {
            usage_exit();
        }
    }

        // parse arguments

    if (!*argv) usage_exit();
    gfa_fname = *argv++;

    std::ifstream gfa_file(gfa_fname);
    if (!gfa_file)
        raise_error("failed to open file: %s", gfa_fname.c_str());

    if (!*argv) usage_exit();
    from_tgt = *argv++;

    if (!*argv) usage_exit();
    to_tgt = *argv++;

    if (*argv) usage_exit();

        // read GFA (and FASTA) into graph

    gfa::graph graph;

    verbose_emit("reading GFA file: %s", gfa_fname.c_str());
    if (!fna_fname.empty()) {

        std::ifstream fna_file(fna_fname);
        if (!fna_file)
            raise_error("failed to open file: %s", fna_fname.c_str());
        verbose_emit("reading FASTA from file: %s", fna_fname.c_str());

        graph = parse_gfa(gfa_file, fna_file, 2, 2*4);  // reserve 2 segments and 8 arcs for start and end
    }
    else {
        graph = parse_gfa(gfa_file, 2, 2*4);            // reserve 2 segments and 8 arcs for start and end
    }

        // add targets to the graph

    add_targets(graph, from_tgt, to_tgt, undirected);

        // call dijkstra

    if (!get_verbose())
        std::cerr << "nothing happening, try --verbose" << std::endl;

    return 0;
}

// vim: sts=4:sw=4:et:si:ai
