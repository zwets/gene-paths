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
#include "utils.h"

using namespace gene_paths;

namespace {

TEST(utils_test, dummy) {
    for (int i=0; i < 10; ++i) {
        EXPECT_EQ(i, -i * -1);
    }
}

TEST(utils_test, raise_error) {
    ASSERT_EXIT( raise_error("test raise"), 
            testing::ExitedWithCode(1), ": error: test raise");
}

TEST(utils_test, set_verbose) {
    ASSERT_EXIT({ set_verbose(true); verbose_emit("test verbose"); std::exit(0); },
            testing::ExitedWithCode(0), ": test verbose");
}
TEST(utils_test, no_set_verbose) {
    ASSERT_EXIT({ verbose_emit("test no verbose"); std::exit(0); },
            testing::ExitedWithCode(0), "");
}

TEST(utils_test, set_progname_exit) {
    ASSERT_EXIT({ set_progname("test name"); raise_error("test raise"); }, 
            testing::ExitedWithCode(1), "test name: error: test raise");
}

TEST(utils_test, set_error_args) {
    ASSERT_EXIT({ raise_error("test message: %d", 42); }, 
            testing::ExitedWithCode(1), ": error: test message: 42");
}

} // namespace
// vim: sts=4:sw=4:ai:si:et
