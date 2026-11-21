GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Wno-unused-parameter -Icxxopts/include -Izstr/src `pkg-config --cflags zlib`

# Pick one of these two
# -----
# Use this by default when bamtools lib is available.
# silly workaround: bamtools does not have pkg-config cflags for finding the include directory
# instead assume it's a folder at the same location as zlib
CPPFLAGS+=`pkg-config --cflags zlib`/bamtools
LIBS=-lbamtools
# -----
# Use this if bamtools lib is not available and cannot be installed. This will disable all functionality that uses bam files and only allows sam files
#CPPFLAGS+=-DNOBAMTOOLSLIB
# -----

ODIR=obj
BINDIR=bin
SRCDIR=src


LINKFLAGS= $(CPPFLAGS) `pkg-config --libs zlib` -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

_DEPS = Logger.h Common.h FileHandler.h EM.h AlleleSpecificExpression.h SNPCounter.h PseudobulkAnalysis.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = Logger.o Common.o FileHandler.o EM.o AlleleSpecificExpression.o SNPCounter.o PseudobulkAnalysis.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

VERSION := Git branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

TEST_SRC = $(wildcard test/*.cpp)
TEST_DEPS = $(wildcard test/*.h)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

all: $(BINDIR)/XCIE-EM $(BINDIR)/CountSNPs

$(BINDIR)/XCIE-EM: $(DEPS) $(OBJ) $(SRCDIR)/main.cpp
	$(GPP) -o $@ $^ $(LINKFLAGS) -DVERSION="\"$(VERSION)\""

$(BINDIR)/CountSNPs: $(DEPS) $(OBJ) $(SRCDIR)/CountSNPs.cpp
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(BINDIR)/tests: $(DEPS) $(OBJ) $(TEST_SRC) $(TEST_DEPS)
	$(GPP) -o $@ $(OBJ) $(TEST_SRC) $(LINKFLAGS) -Idoctest/doctest

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*

test: $(BINDIR)/tests
	$(BINDIR)/tests
