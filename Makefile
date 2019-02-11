all:
	make -C common
	make -C srcEpos
	mkdir -p build
	cp srcEpos/epos build
	make -C srcEpos2ages
	cp srcEpos2ages/epos2ages build
	make -C srcEpos2plot
test:
	make -s -C srcEpos test
	make -s -C srcEpos2ages test
	make -s -C srcEpos2plot test
clean:
	make -C common clean
	make -C srcEpos clean
	make -C srcEpos2ages clean
	make -C srcEpos2plot clean
	make -C doc clean
.PHONY:	doc
doc:	
	make -C doc
