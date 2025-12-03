GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Wno-unused-parameter

ODIR=obj
BINDIR=bin
SRCDIR=src

$(shell mkdir -p bin)
$(shell mkdir -p obj)

all: $(BINDIR)/EM

$(BINDIR)/EM: src/EM.cpp
	$(GPP) -o $@ $^ $(CPPFLAGS)

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
