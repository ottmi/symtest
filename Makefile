CXX = g++
CXXFLAGS = -g3 -gdwarf-2 -O3 
LDFLAGS = 
LIBS = -lm

VERSION := $(shell awk '/VERSION/ {print $$3}' src/globals.h| sed 's/\"\(.*\)\"/\1/')
MODULES := Alignment.cpp AlignmentReader.cpp Sequence.cpp helper.cpp symtest.cpp
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
	tar --transform "s,^,shuffle-$(VERSION)/," -cjf shuffle-$(VERSION).tar.bz2 configure Makefile.in src/*.cpp src/*.h


clean:
	$(RM) $(BIN) $(OBJ)
