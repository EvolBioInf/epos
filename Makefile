TAG := $(shell git describe | sed 's/-.*//')
DATE := $(shell git log -1 --format=%ai $(TAG) | sed 's/-/ /g' | awk '{printf "\\\\DTMdisplaydate{%s}{%s}{%s}{-1}\n", $$1, $$2, $$3}')
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
	echo $(TAG)  > doc/version.tex
	echo $(DATE) > doc/date.tex
	make -C doc
