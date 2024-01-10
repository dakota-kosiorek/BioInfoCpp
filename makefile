IDIR =include
CC=g++
CFLAGS=-I$(IDIR)

# Adjusted paths for src directory
SRCDIR = src

ODIR=$(SRCDIR)/obj
LDIR =lib

LIBS=-lm

_DEPS = analysis.hpp biomath.hpp fundamentals.hpp genetics.hpp query.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o analysis.o biomath.o fundamentals.o genetics.o query.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

testing: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 