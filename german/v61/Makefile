PYTHON ?= $(if $(wildcard ../../.venv/bin/python),../../.venv/bin/python,python3)
DATA := $(wildcard data/*.txt* data/*.HEIC)
LATEX_ENV = TEXINPUTS=..: BIBINPUTS=..: max_print_line=1048576 \
	TEXMFCACHE="$(CURDIR)/build/tex-cache" TEXMFVAR="$(CURDIR)/build/tex-cache"

all: build/v61.pdf

# plot.py creates every plot, measured-data table, and correlated fit result.
build/results.tex: plot.py $(DATA) ../matplotlibrc ../header-matplotlib.tex ../../pyproject.toml ../../uv.lock | build
	$(PYTHON) plot.py > build/analysis.log

build/v61.pdf: v61.tex build/results.tex FORCE | build
	$(LATEX_ENV) latexmk --lualatex --output-directory=build \
		--interaction=nonstopmode --halt-on-error $<

preview: build/results.tex | build
	$(LATEX_ENV) latexmk -pvc --lualatex --output-directory=build \
		--interaction=nonstopmode --halt-on-error v61.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean preview FORCE
