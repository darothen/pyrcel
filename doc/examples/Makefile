
NOTEBOOKS := $(wildcard *.ipynb)
RSTS := $(NOTEBOOKS:.ipynb=.rst)

.PHONY: all
all: $(RSTS)

%.rst: %.ipynb
	jupyter-nbconvert $^ --to rst

clean:
	rm *.rst
	rm -rf *_files/
	rm -rf _build/