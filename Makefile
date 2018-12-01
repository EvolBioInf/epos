export VERSION = $(shell git describe)
all:
	make -C srcEpos
	mkdir -p build
	mv srcEpos/epos build
	make -C srcEpos2ages
	mv srcEpos2ages/epos2ages build
test:
	make -C srcEpos test
	make -C srcEpos2ages test
clean:
	make -C srcEpos clean
	make -C srcEpos2ages clean
