/* utils-test.cpp
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

#include <gtest/gtest.h>
#include <sstream>
#include <fstream>
#include "gfakluge.hpp"
#include "parsegfa.h"

using namespace gene_paths;

namespace {

TEST(parsegfa_test, read_gfa) {

    std::ifstream gfa_file("data/with_seqs.gfa");
    ASSERT_TRUE(gfa_file);
    
    gfak::GFAKluge gfak;
    parse_gfa(gfak, gfa_file);

    auto name2seq = gfak.get_name_to_seq();
    ASSERT_EQ(name2seq.size(), 9);
}

TEST(parsegfa_test, read_gfa_and_fna) {

    std::ifstream gfa_file("data/without_seqs.gfa");
    ASSERT_TRUE(gfa_file);

    std::ifstream fna_file("data/seqs.fna");
    ASSERT_TRUE(fna_file);

    gfak::GFAKluge gfak;
    parse_gfa(gfak, gfa_file, fna_file);

    auto name2seq = gfak.get_name_to_seq();
    ASSERT_EQ(name2seq.size(), 9);
}

TEST(parsegfa_test, read_gfa_string) {

    std::istringstream gfa("H\tVN:Z:2.0\nS\t1\t4\t*\n");
    std::istringstream fna(">1\nACGT\n");

    gfak::GFAKluge gfak;
    parse_gfa(gfak, gfa, fna);
    ASSERT_EQ(gfak.get_name_to_seq().size(), 1);
}

TEST(parsegfa_test, read_gfa_mismatch_fna) {

    std::istringstream gfa("H\tVN:Z:2.0\nS\t1\t4\t*\n");
    std::istringstream fna(">1\nACG\n");

    gfak::GFAKluge gfak;
    ASSERT_EXIT( parse_gfa(gfak, gfa, fna);,
            testing::ExitedWithCode(1), 
            ": error: FASTA length \\(3\\) does not match length in GFA \\(4\\) for seqid: 1");
}

} // namespace
// vim: sts=4:sw=4:ai:si:et
