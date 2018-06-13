CC=gcc
CFLAGS= -O3 -Wall -Wshadow -pedantic -std=gnu99 -g #-fopenmp # -pthread -fsanitize=thread #-pg
# The source files, object files, libraries and executable name.
SRCFILES= epos.c interface.c eprintf.c sfs.c tab.c popSizes.c printMatrix.c util.c foldedE.c unfolded.c gsl_rng.c xval.c
OBJFILES= epos.o interface.o eprintf.o sfs.o tab.o popSizes.o printMatrix.o util.o foldedE.o unfolded.c gsl_rng.o xval.o
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
