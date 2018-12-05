TAG := $(shell git describe | sed 's/-.*//')
all:
	make -C common
	make -C srcEpos
	mkdir -p build
	cp srcEpos/epos build
	make -C srcEpos2ages
	cp srcEpos2ages/epos2ages build
test:
	make -C srcEpos test
	make -C srcEpos2ages test
clean:
	make -C common clean
	make -C srcEpos clean
	make -C srcEpos2ages clean
	make -C doc clean
.PHONY:	doc
doc:	
	echo $(TAG) > doc/version.tex
	make -C doc
