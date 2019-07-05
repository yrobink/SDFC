
library(roxygen2)
library(devtools)
roxygenize("SDFC")
load_all("SDFC")
roxygenize("SDFC")
build("SDFC")


