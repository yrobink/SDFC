#!/bin/sh

## Delete python temporary files
rm -rf python/SDFC.egg*
rm -rf python/build
rm -rf python/dist
rm -rf python/tmp
rm -rf python/var


## Delete R temporary files
rm -f R/SDFC/NAMESPACE
rm -f R/SDFC/man/*.Rd
rm -f R/*.tar.gz
rm -f R/SDFC/src/RcppExports.cpp
rm -f R/SDFC/R/RcppExports.R
rm -f R/SDFC/src/*.o
rm -f R/SDFC/src/*.so

