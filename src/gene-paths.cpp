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
#include "gfakluge.hpp"
#include "parsegfa.h"
#include "utils.h"


using namespace gene_paths;


static const std::string USAGE(
"Usage: gene-paths [OPTIONS] GFA_FILE\n"
"\n"
"  Placeholder usage information.\n"
"\n"
"  OPTIONS\n"
"   -f, --fasta FILE  read sequences for GFA_FILE from FILE\n"
"   -v, --verbose     write progress information to stderr\n"
"   -h, --help        output this information and exit\n"
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

    while (*++argv && **argv == '-') 
    {
        if (!std::strcmp("-v", *argv) || !std::strcmp("--verbose", *argv)) {
            set_verbose(true);
        }
        else if (!std::strcmp("-h", *argv) || !std::strcmp("--help", *argv)) {
            std::cout << USAGE;
            return 0;
        }
        else if ((!std::strcmp("-f", *argv) || !std::strcmp("--fasta", *argv)) && *++argv) {
            fna_fname = *argv;
        }
        else {
            usage_exit();
        }
    }

    if (!*argv)
        usage_exit();

    gfa_fname = *argv++;

    if (*argv)
        usage_exit();

    gfak::GFAKluge gfa;

    std::ifstream gfa_file(gfa_fname);
    if (!gfa_file)
        raise_error("failed to open file: %s", gfa_fname.c_str());

    if (!fna_fname.empty()) {

        std::ifstream fna_file(fna_fname);
        if (!fna_file)
            raise_error("failed to open file: %s", fna_fname.c_str());

        verbose_emit("reading GFA file: %s", gfa_fname.c_str());
        verbose_emit("reading FASTA file: %s", fna_fname.c_str());
        parse_gfa(gfa, gfa_file, fna_file);
    }
    else {
        verbose_emit("reading GFA file: %s", gfa_fname.c_str());
        parse_gfa(gfa, gfa_file);
    }

    return 0;
}

// vim: sts=4:sw=4:et:si:ai
