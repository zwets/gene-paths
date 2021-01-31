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

#include <iostream>
#include <fstream>
#include <string>
#include "gfakluge.hpp"
#include "utils.h"

namespace gene_paths {

void
parse_gfa(gfak::GFAKluge& gfak, std::istream& file)
{
    if (!gfak.parse_gfa_file(file))
        raise_error("failed to parse GFA");

    verbose_emit("GFA file parsed successfully");
}

void
parse_gfa(gfak::GFAKluge& gfak, std::istream& gfa, std::istream& fasta)
{
    parse_gfa(gfak, gfa);

    auto name2seq = gfak.get_name_to_seq(); // maps name to seq_elem
    int n_found = 0, n_skipped = 0;
    std::string line;
    std::string seqid;
    std::string data;

    std::getline(fasta, line);

    while (line.length()) {

        if (line[0] != '>')
            raise_error("invalid FASTA header: %s", line.c_str());

        std::string::const_iterator p = line.begin() + 1;
        while (p != line.end() && !std::isspace(*p)) ++p;
        seqid = std::string(line, 1, p - line.begin() - 1);

        auto found = name2seq.find(seqid);
        if (found != name2seq.end()) {

            ++n_found;
            data.clear();

            while (std::getline(fasta, line)) {
                if (!line.length())
                    continue;
                if (line[0] == '>')
                    break;
                data.append(line);
                line.clear();
            }

            if (data.length() == found->second.length)
                found->second.sequence = data;
            else
                raise_error("FASTA length (%d) does not match length in GFA (%d) for seqid: %s",
                        data.length(), found->second.length, seqid.c_str());
        }
        else 
            ++n_skipped;
    }

    verbose_emit("attached %d sequences from FASTA; not found: %d", n_found, n_skipped);
}

} // namespace gene_paths

// vim: sts=4:sw=4:ai:si:et
