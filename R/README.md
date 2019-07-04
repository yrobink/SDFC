
# SDFC (Statistical Distribution Fit with Covariates)


## Features
- python3 and R version
- c++ independent files for QuantileRegression
- Fit with (or without) covariates parametric laws (Normal, Exp, Gamma, GEV, GPD)
- Fit of non parametric distribution (QuantileRegression)
- Support for fit with fix parameters
- Support for user defined custom link function


## R instruction

Requires:
- R
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- devtools
- roxygen2
- Rcpp
- RcppEigen

I'm not really sure of how to build, but this sequence works:
```
roxygen2::roxygenize("SDFC") ## Return an error, but necessary to generate NAMESPACE file...
devtools::load_all("SDFC")
roxygen2::roxygenize("SDFC") ## Now, no errors
devtools::build("SDFC")      ## Generate the package
install.packages( "SDFC_version.tar.gz" )
```


## Note for quantile regression

The quantile regression is solved with the Frish-Newton algorithm, written in c++,
depending of the Eigen c++ library. It is a re-written of the Fortran code of
Koenker, available in the R package [quantreg](https://cran.r-project.org/web/packages/quantreg/index.html)


## License

Copyright Yoann Robin, 2019, [CeCILL-C](https://cecill.info/licences.en.html) license (open source)

This software is a computer program that is part of the SDFC (Statistical
Distribution Fit with Covariates) library. This library makes it possible
to regress the parameters of some statistical law with co-variates.

This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-C 
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-C license and that you accept its terms.



