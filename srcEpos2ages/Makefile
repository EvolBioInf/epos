CC=gcc
export VERSION = $(shell bash ../scripts/getVersion.sh) 
export DATE    = $(shell bash ../scripts/getDate.sh) 
CFLAGS= -DVERSION=\"$(VERSION)\" -DDATE=\"$(DATE)\" -O3 -Wall -Wshadow -pedantic -std=gnu99 -I../common -L../common -g
# The source files, object files, libraries and executable name.
SRCFILES= epos2ages.c interface.c util.c popSizes.c ages.c
OBJFILES= epos2ages.o interface.o util.o popSizes.o ages.o
LIBS= -lm -lgsl -lblas -lcommon
EXECFILE=epos2ages
DIRECTORY=Epos2ages
# The make rule for the executable
.PHONY : all
all : $(EXECFILE)
$(EXECFILE) : $(OBJFILES)

	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)
interface.o: interface.h
eprintf.o: eprintf.h

doc:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)

# Other Standard make rules
lint : 
	lint $(SRCFILES) | more
clean:
	rm -f *.o *~
tarfile:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)
	cp -rf $(SRCFILES) *.h  README COPYING Makefile ../Doc/epos2agesDoc.pdf $(DIRECTORY)_$(VERSION)
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	/bin/rm -r $(DIRECTORY)_$(VERSION)

test:
	bash ../scripts/testEpos2ages.sh 
