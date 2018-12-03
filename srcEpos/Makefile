CC=gcc
export VERSION = $(shell git describe)
export DATE    = $(shell git log | grep Date | head -n 1 | awk '{printf "%s\\ %s\\ %s\\ %s\\ %s\n", $$2, $$3, $$4, $$5, $$6}')
CFLAGS = -DVERSION=\"$(VERSION)\" -DDATE=\"$(DATE)\" -O3 -Wall -Wshadow -pedantic -std=gnu99 -g #-pg
# The source files, object files, libraries and executable name.
SRCFILES= epos.c interface.c eprintf.c sfs.c tab.c util.c popSizes.c \
newton.c greedy.c exhaustive.c test.c search.c 
OBJFILES= epos.o interface.o eprintf.o sfs.o tab.o util.o popSizes.o \
newton.o greedy.o exhaustive.o test.o search.o 
LIBS= -lm -lgsl -lblas #-pg
EXECFILE=epos
DIRECTORY=Epos
# The make rule for the executable
.PHONY : all
all : $(EXECFILE)
$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)
interface.o: interface.h
eprintf.o: eprintf.h
# Other Standard make rules
lint : 
	lint $(SRCFILES) | more
clean:
	rm -f *.o *~
test:
	bash ../scripts/testEpos.sh 
