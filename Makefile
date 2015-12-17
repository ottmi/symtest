CXX = g++
CXXFLAGS = -O3 -Wall -pedantic
LDFLAGS = 
LIBS = -lm

VERSION := $(shell awk '/VERSION/ {print $$3}' src/globals.h| sed 's/\"\(.*\)\"/\1/')
MODULES := Alignment.cpp AlignmentReader.cpp Matrix.cpp Sequence.cpp chi.cpp heatmap.cpp helper.cpp symtest.cpp
SRC := $(addprefix src/,$(MODULES))
OBJ := $(patsubst src/%.cpp,src/%.o,$(SRC))
INCLUDES := $(addprefix -I,$(SRC_DIR))
BIN = symtest

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(LDFLAGS) $(LIBS) -o $(BIN) $(OBJ)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

tarball:
	tar --transform "s,^,symtest-$(VERSION)/," -cjf symtest-$(VERSION).tar.bz2 Makefile src/*.cpp src/*.h

clean:
	$(RM) $(BIN) $(OBJ)
