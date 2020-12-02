# SDFC (Statistical Distribution Fit with Covariates)


## Features
- python3 and R version
- c++ independent files for QuantileRegression
- Fit with (or without) covariates parametric laws (Normal, Exp, Gamma, GEV, GPD)
- MLE and Bayesian fit available for all laws
- Fit of non parametric distribution (QuantileRegression)
- Support for fit with fix parameters
- Support for user defined custom link function

## Example

![Alt](/figures/Example_py.png)



## Note for quantile regression

The quantile regression is solved with the Frish-Newton algorithm, written in c++,
depending of the Eigen c++ library. It is a re-written of the Fortran code of
Koenker, available in the R package [quantreg](https://cran.r-project.org/web/packages/quantreg/index.html)


## License

Copyright(c) 2020 Yoann Robin

This file is part of SDFC.

SDFC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SDFC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SDFC.  If not, see <https://www.gnu.org/licenses/>.

