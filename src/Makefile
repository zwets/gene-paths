CXXFLAGS += -std=c++14 -O3 -DNDEBUG -Wall -Wextra -pedantic #-march=native

OBJS = gene-paths.o utils.o 

LIBS =

HDRS = utils.h

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
	$(CXX) $(CXXFLAGS) -c $<

