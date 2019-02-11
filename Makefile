all:
	make -C common
	make -C srcEpos
	mkdir -p build
	cp srcEpos/epos build
	make -C srcEpos2ages
	cp srcEpos2ages/epos2ages build
test:
	make -s -C srcEpos test
	make -s -C srcEpos2ages test
clean:
	make -C common clean
	make -C srcEpos clean
	make -C srcEpos2ages clean
.PHONY:	doc
doc:	
	make -C doc
