PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: install

build:
	cd ..;\
	R CMD build $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGSRC)

check:
	Rscript -e "devtools::check(document = FALSE, args = '--as-cran')"

precompile:
	Rscript -e "source('vignettes/precompile.R')"

test:
	Rscript -e "library('testthat',quietly=TRUE);library('Spower',quietly=TRUE);options(warn=2);test_dir('tests/testthat')"

kniterrors:
	grep -Hrn '## Error'
	grep -Hrn '## Warning'

pkgdown:
	sed -i 's/# opts$$verbose/opts$$verbose/g' R/03-estimation.R
	sed -i s/dontrun/donttest/g man/*.Rd
	make install
	Rscript -e "library('pkgdown',quietly=TRUE);build_site()"
	git checkout -- .
	make install
	cd docs/reference
	make kniterrors

clean:
	$(RM) src/*.o
	$(RM) src/*.so
	$(RM) ../$(PKGNAME)_$(PKGVERS).tar.gz
	$(RM) -r ../$(PKGNAME).Rcheck/


