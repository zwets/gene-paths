# Adapted from the sample Makefile that comes with libgtest-dev on Ubuntu.
#
# SYNOPSIS
#   make clean
#   make [all]
#   make test

# Points to the root of Google Test, relative to where this file is.
GTEST_DIR = ./gtest

# Where to find the code under test.
USER_DIR = ..

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include

# Flags passed to the C++ compiler.
CXXFLAGS += -std=c++14 -g -Wall -Wextra -pedantic -Wno-unknown-pragmas -march=native -pthread -I$(USER_DIR) -I$(USER_DIR)/gfakluge

OBJS =
# All Google Test headers.  Normally you shouldn't change this definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# USER SECTION

TARGET = run-all-tests

USER_HEADERS = $(USER_DIR)/*.h

USER_OBJS = dijkstra.o paths.o targets.o graph.o gfa2logic.o parser.o utils.o

TEST_OBJS = utils-test.o parser-test.o gfa2logic-test.o graph-test.o targets-test.o paths-test.o dijkstra-test.o 

# Build targets.

all : $(TARGET)

clean :
	rm -f $(TARGET) $(TEST_OBJS) $(USER_OBJS) gtest.a gtest_main.a gtest-all.o gtest_main.o

test : $(TARGET)
	./$(TARGET)

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# USER RULES

%.o : $(USER_DIR)/%.cpp $(USER_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

%-test.o : %-test.cpp $(USER_HEADERS) $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

$(TARGET): $(TEST_OBJS) $(USER_OBJS) gtest_main.a $(USER_LIBS)
	$(CXX) -pthread $^ -o $@

