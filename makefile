GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Wno-unused-parameter -Icxxopts/include

ODIR=obj
BINDIR=bin
SRCDIR=src

_DEPS = Logger.h Common.h FileHandler.h EM.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = main.o Logger.o Common.o FileHandler.o EM.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

VERSION := Git branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

all: $(BINDIR)/XCIE-EM

$(BINDIR)/XCIE-EM: $(DEPS) $(OBJ)
	$(GPP) -o $@ $^ $(CPPFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(ODIR)/main.o: $(SRCDIR)/main.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
