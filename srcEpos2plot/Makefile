export VERSION = $(shell bash ./scripts/version.sh)
export DATE    = $(shell bash ./scripts/date.sh)
all:
	go build -ldflags "-X main.VERSION=${VERSION} -X main.DATE=${DATE}"
install:
	go install -ldflags "-X main.VERSION=${VERSION} -X main.DATE=${DATE}"
test:
	bash ./scripts/test.sh
.PHONY:	doc
doc:	
	make -C doc
clean:
	rm -f epos2plot

