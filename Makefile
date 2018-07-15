build: attrib
	R CMD build .

attrib:
	R -e 'Rcpp::compileAttributes()'

install:
	R CMD INSTALL --build .

doc:
	R CMD Rd2pdf -o TDA.pdf .

test:
	R -e 'devtools::test()'

check:
	R CMD check --as-cran $(TGZ)

.PHONY: veryclean clean

veryclean: clean
	rm -f src/TDA.so  TDA_*.tar.gz TDA.pdf
	rm -rf TDA.Rcheck

clean:
	find . -name '*.o' -exec rm {} \;
	rm -f src/RcppExports.cpp R/RcppExports.R

