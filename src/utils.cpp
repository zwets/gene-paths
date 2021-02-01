/* utils.cpp
 * 
 * Copyright (C) 2018,2021  Marco van Zwetselaar <io@zwets.it>
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

#include "utils.h"

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>

namespace gene_paths {

static bool verbose = false;
static const char* progname = "";

void
set_progname(const char *p)
{
    progname = p;
}

void
set_verbose(bool v)
{
    verbose = v;
}

void
raise_error(const char *fmt, ...)
{
    char buf[2048];

    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);

    std::cerr << progname << ": error: " << buf << std::endl;
    std::exit(1);
}

void
verbose_emit(const char *fmt, ...)
{
    if (verbose)
    {
        char buf[2048];

        va_list ap;
        va_start(ap, fmt);
        vsnprintf(buf, sizeof(buf), fmt, ap);
        va_end(ap);

        std::cerr << progname << ": " << buf << std::endl;
    }
}

} // namespace gene_paths

// vim: sts=4:sw=4:ai:si:et
