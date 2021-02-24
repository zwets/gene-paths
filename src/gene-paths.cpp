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
#include <cstring>
#include <cstdlib>
#include "parsegfa.h"
#include "gfagraph.h"
#include "dijkstra.h"
#include "targets.h"
#include "utils.h"


using namespace gene_paths;


static const std::string USAGE(
"Usage: gene-paths [OPTIONS] GFA_FILE FROM TO\n"
"\n"
"  Output the shortest path between regions FROM and TO in the assembly\n"
"  graph in GFA_FILE.\n"
"\n"
"  FROM and TO are specified as SEG[+-][:BEG[:END]]], where SEG is the\n"
"  contig identifier, + or - the STRAND, and BEG and END the start and\n"
"  end positions of the target region on SEG.\n"
"\n"
"  When END is omitted, BEG specifies a zero-length position.  When BEG\n"
"  and END are both omitted they default to the start and end of SEG.\n"
"\n"
"  OPTIONS\n"
"   -f, --fasta FILE  read sequences for GFA_FILE from FILE\n"
"   -b, --both-dirs   search for FROM both up- and downstream of TO\n"
"   -v, --verbose     write detailed informationo stderr\n"
"   -h, --help        output this information and exit\n"
"\n"
"  The default is to search for a path going downstream from FROM.  Use\n"
"  option -b/--both to also search for a path that has TO upstream.\n"
"\n"
"  STRAND and POSITION are interpreted as in GFA2:\n"
"  - We call the sequence data in GFA_FILE (or FASTA FILE) the + strand.\n"
"    The - strand is its reverse complement.\n"
"  - Positions are inbetween bases, starting at 0 left of the sequence,\n"
"    and ending at L at its right, where L is the sequence length.\n"
"  - Positions are interpreted on the segment BEFORE it is oriented, so\n"
"    e.g. C-:0:5 refers to the FINAL 5 bases on C's reverse complement.\n"
"\n");

static void usage_exit()
{
    std::cerr << USAGE;
    std::exit(1);
}


int main (int /*argc*/, char *argv[])
{
    set_progname("gene-paths");

    std::string gfa_fname;
    std::string fna_fname;
    bool both_dirs = false;

        // parse options

    while (*++argv && **argv == '-')
    {
        if (!std::strcmp("-v", *argv) || !std::strcmp("--verbose", *argv)) {
            set_verbose(true);
        }
        else if (!std::strcmp("-h", *argv) || !std::strncmp("--help", *argv, 6)) {
            std::cout << USAGE;
            return 0;
        }
        else if (!std::strcmp("-b", *argv) || !std::strncmp("--both", *argv, 6)) {
            both_dirs = true;
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
    std::string from_ref = *argv++;

    if (!*argv) usage_exit();
    std::string to_ref = *argv++;

    if (*argv) usage_exit();

        // read GFA (and FASTA) into graph

    gfa::graph g;

    verbose_emit("reading GFA file: %s", gfa_fname.c_str());
    if (!fna_fname.empty()) {

        std::ifstream fna_file(fna_fname);
        if (!fna_file)
            raise_error("failed to open file: %s", fna_fname.c_str());
        verbose_emit("reading FASTA from file: %s", fna_fname.c_str());

        g = parse_gfa(gfa_file, fna_file, 4, 2*2);  // reserve 2 segs and 4 arcs
    }
    else {
        g = parse_gfa(gfa_file, 2, 2*2);            // reserve 2 segs and 4 arcs
    }

        // add targets to the graph

    gfa::target from(g), to(g);

    from.set(from_ref, gfa::target::START);
    to.set(to_ref, gfa::target::END);

        // call dijkstra

    gfa::dijkstra dijkstra(g);

    verbose_emit("searching shortest path: %s -> %s", from_ref.c_str(), to_ref.c_str());

    bool fwd = dijkstra.shortest_path(from.p_arc(), to.p_arc());
    if (fwd)
    {
        std::cout << ">PATH ";
        dijkstra.write_route(std::cout);
        std::cout << " (length " << dijkstra.length() << ")";
        std::cout << std::endl;

        dijkstra.write_sequence(std::cout);
        std::cout << std::endl;
    }

        // if bidirectional, call again with from and to swapped

    bool rev = false;
    if (both_dirs) {

        from.set(to_ref, gfa::target::START);
        to.set(from_ref, gfa::target::END);

        verbose_emit("searching inverse path: %s -> %s", to_ref.c_str(), from_ref.c_str());

        rev = dijkstra.shortest_path(from.p_arc(), to.p_arc());
        if (rev)
        {
            std::cout << ">PATH_REV ";
            dijkstra.write_route(std::cout);
            std::cout << " (length " << dijkstra.length() << ")";
            std::cout << std::endl;
            dijkstra.write_sequence(std::cout);
            std::cout << std::endl;
        }
    }

    if (!(fwd || rev)) {
        std::cerr << "No path was found\n";
    }

    return fwd || rev;
}

// vim: sts=4:sw=4:et:si:ai
