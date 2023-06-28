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
#include "graph.h"
#include "parser.h"
#include "targets.h"
#include "dijkstra.h"
#include "utils.h"

using namespace gene_paths;

static const std::string USAGE(
"Usage: gene-paths [OPTIONS] GFA_FILE FROM TO\n"
"\n"
"  Find the shortest path between locations FROM and TO in the genome\n"
"  assembly graph in GFA_FILE.\n"
"\n"
"  OPTIONS\n"
"   -b, --bidir       search for TO both upstream and downstream of FROM\n"
"   -f, --fasta FILE  read sequences for GFA_FILE from FILE\n"
"   -v, --verbose     write detailed progress information to stderr\n"
"   -h, --help        print this information and exit\n"
"\n"
"  The path search looks for TO downstream of FROM.  Use option -b/--bidir\n"
"  to also search for a path that has TO upstream of FROM.  Both paths (if\n"
"  any exist) will be reported.\n"
"\n"
"  FROM and TO are specified as CTG[:BEG[:END]]S, where CTG is the name of\n"
"  the contig, BEG and END are the optional start and end positions on CTG,\n"
"  and S is the mandatory strand identifier (+ or -).\n"
"\n"
"  BEG and END can be '$' to signify the end of CTG.  When END is omitted\n"
"  it defaults to BEG, so the reference is a (zero length) position.\n"
"  When BEG and END are both omitted, they default to 0 and $ respectively,\n"
"  so the reference is the whole contig.\n"
"\n"
"  STRAND and POSITION are interpreted as in GFA2:\n"
"  - We define the data in the GFA or FASTA file to be the + strand,\n"
"    and its reverse reverse complement the - strand;\n"
"  - Positions are in between bases, with 0 to the left of the sequence,\n"
"    and $ to the right, $ being the sequence length;\n"
"  - Positions are interpreted before orienting the segment, so pos 0 is\n"
"    at the upstream end of a segment, regardless of sign.\n"
"\n"
"  TIP: 'gene-paths 1+ 1:0+' finds the shortest CYCLICAL path that contains\n"
"  contig 1.  It starts out going across all of 1, then searches for a path\n"
"  to where it started.  Note how 'gene-paths 1:$+ 1:0+' is similar, but\n"
"  excludes contig 1 itself.\n"
"\n"
"  Note: if you use '$' in FROM or TO you will likely need to quote it,\n"
"  to prevent interpretation by your command shell.\n"
"\n");

static void usage_exit(int err = 1)
{
    (err ? std::cerr : std::cout) << USAGE;
    std::exit(err);
}

static void write_path(const gfa::dijkstra& d)
{
    if (d.found_pix) {
        std::cout << ">PATH ";
        d.write_route(std::cout);
        std::cout << " (length " << d.length() << ")";
        std::cout << std::endl;

        d.write_sequence(std::cout);
        std::cout << std::endl;
    }
}

int main (int /*argc*/, char *argv[])
{
    set_progname("gene-paths");

    std::string gfa_fname;
    std::string fna_fname;
    bool bidirectional = false;
    bool furthest = false;

        // parse options

    while (*++argv && **argv == '-')
    {
        if (!std::strcmp("-v", *argv) || !std::strcmp("--verbose", *argv)) {
            set_verbose(true);
        }
        else if (!std::strcmp("-h", *argv) || !std::strcmp("--help", *argv)) {
            usage_exit(0);
        }
        // option is parsed and executed, just not documented in --help yet
        else if (!std::strcmp("-u", *argv) || !std::strcmp("--furthest", *argv)) {
            furthest = true;
        }
        else if (!std::strcmp("-b", *argv) || !std::strncmp("--bidir", *argv, 7)) {
            bidirectional = true;
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

    std::string to_ref;
    if (!furthest && *argv)
        to_ref = *argv++;

    if (*argv) usage_exit();

        // read GFA (and FASTA) into graph,
        // reserving 3 segs and 4 arcs for the targets (see paths.h)

    gfa::graph g;

    verbose_emit("reading GFA file: %s", gfa_fname.c_str());
    if (!fna_fname.empty()) {

        std::ifstream fna_file(fna_fname);
        if (!fna_file)
            raise_error("failed to open file: %s", fna_fname.c_str());
        verbose_emit("reading FASTA from file: %s", fna_fname.c_str());

        g = gfa::parse(gfa_file, fna_file, 3, 4);
    }
    else {
        g = gfa::parse(gfa_file, 3, 4);
    }

        // set up for program return value

    bool success = true;

        // create targets on the graph, and set the from

    gfa::target from(g), to(g);

    from.set(from_ref, gfa::target::START);

        // create dijkstra algorithm and point it at the graph

    gfa::dijkstra dijkstra(g);

        // @TODO@ document the -u/--furthest option
        // it finds the shortest path from FROM to every possible TO,
        // then returns the longest of these shortest paths

    if (furthest) // find longest of all shortest paths from FROM
    {
        verbose_emit("searching furthest path from: %s", from_ref.c_str());

        dijkstra.furthest_path(from.p_arc());
        write_path(dijkstra);
    }
    else // find shortest path from FROM to TO
    {
        verbose_emit("searching shortest path: %s -> %s", from_ref.c_str(), to_ref.c_str());

        to.set(to_ref, gfa::target::END);

        success = dijkstra.shortest_path(from.p_arc(), to.p_arc());
        write_path(dijkstra);

        if (bidirectional) // also find shortest path with TO upstream of FROM
        {
            verbose_emit("searching inverse path: %s -> %s", to_ref.c_str(), from_ref.c_str());

            from.set(to_ref, gfa::target::START);
            to.set(from_ref, gfa::target::END);

            success |= dijkstra.shortest_path(from.p_arc(), to.p_arc());
            write_path(dijkstra);
        }

        if (!success) {
            std::cerr << "No path was found\n";
        }
    }

    return success ? 0 : 1;
}

// vim: sts=4:sw=4:et:si:ai
