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
#include <string>
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
"   -u, --undirected  search for TO both up- and downstream from FROM\n"
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

int main (int /*argc*/, char *argv[])
{
    set_progname("gene-paths");

    std::string gfa_fname;
    std::string fna_fname;
    std::string arg_from;
    std::string arg_to;
    bool arg_undirected = false;

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
            arg_undirected = true;
        }
        else if ((!std::strcmp("-f", *argv) || !std::strcmp("--fasta", *argv)) && *++argv) {
            fna_fname = *argv;
        }
        else {
            usage_exit();
        }
    }

    if (!*argv) usage_exit();
    gfa_fname = *argv++;

    std::ifstream gfa_file(gfa_fname);
    if (!gfa_file)
        raise_error("failed to open file: %s", gfa_fname.c_str());

    if (!*argv) usage_exit();
    arg_from = *argv++;

    if (!*argv) usage_exit();
    arg_to = *argv++;

    if (*argv) usage_exit();

    gfa::graph graph;

    verbose_emit("reading GFA file: %s", gfa_fname.c_str());
    if (!fna_fname.empty()) {

        std::ifstream fna_file(fna_fname);
        if (!fna_file)
            raise_error("failed to open file: %s", fna_fname.c_str());

        verbose_emit("reading FASTA from file: %s", fna_fname.c_str());
        graph = parse_gfa(gfa_file, fna_file);
    }
    else {
        graph = parse_gfa(gfa_file);
    }

    std::cout << graph.segs.size() <<std::endl;

    return 0;
}

// vim: sts=4:sw=4:et:si:ai
