CXXFLAGS += -std=c++14 -O3 -DNDEBUG -Wall -Wextra -pedantic -Wno-unknown-pragmas -march=native

OBJS = gene-paths.o gfagraph.o parsegfa.o utils.o 

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

