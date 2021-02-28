CXXFLAGS += -std=c++14 -O3 -DNDEBUG -Wall -Wextra -pedantic -Wno-unknown-pragmas -march=native
# For debug:
#CXXFLAGS += -std=c++14 -g -Wall -Wextra -pedantic -Wno-unknown-pragmas -march=native

OBJS = gene-paths.o dijkstra.o paths.o targets.o graph.o gfa2logic.o parser.o utils.o

LIBS =

HDRS = *.h gfakluge/*.hpp

TARGET = gene-paths

$(TARGET): $(OBJS) $(HDRS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all: $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
	$(MAKE) -C unit-test clean

test:
	$(MAKE) -C unit-test test

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -I gfakluge -c $<

