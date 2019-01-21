CC=gcc
export VERSION = $(shell bash ../scripts/getVersion.sh) 
export DATE    = $(shell bash ../scripts/getDate.sh) 
CFLAGS = -DVERSION=\"$(VERSION)\" -DDATE=\"$(DATE)\" -O3 -Wall -Wshadow -pedantic -std=gnu99 -g -I../common/ -L../common/ -m64 #-pg
# The source files, object files, libraries and executable name.
SRCFILES= epos.c interface.c sfs.c util.c popSizes.c  \
newton.c greedy.c exhaustive.c test.c search.c xval.c \
gsl_rng.c
OBJFILES= epos.o interface.o sfs.o util.o popSizes.o  \
newton.o greedy.o exhaustive.o test.o search.o xval.o \
gsl_rng.o
LIBS= -lm -lgsl -lblas -lcommon #-pg
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
