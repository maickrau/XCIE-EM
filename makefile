GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Wno-unused-parameter -Icxxopts/include

ODIR=obj
BINDIR=bin
SRCDIR=src

VERSION := Git branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

all: $(BINDIR)/XCIE-EM

$(BINDIR)/XCIE-EM: src/EM.cpp
	$(GPP) -o $@ $^ $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
