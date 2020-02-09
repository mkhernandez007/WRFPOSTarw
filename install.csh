#!/bin/csh

#==============================================================================
# This code is the installation code to install the GRIB I library only, and to
# copile the NetCDF2GRIB_V1.f coverter code.
#------------------------------------------------------------------------------
# @creator: Michael Kevin Hernandez
# @email  : mkh182@psu.edu
# @usage  : csh install.csh
#==============================================================================

#===================================
# User defined variables
#===================================

set CURRENTdir      = `pwd`
set LIBNETCDF       = /usr/local/netcdf-3.6.1/lib/libnetcdf.a


#
## DON'T edit the variables below.
#

set LIBGRIBdir      = $CURRENTdir/w3lib
set LIBGRIB         = $CURRENTdir/w3lib/lib/libw3.a


#
## Unpacking the GRIB 1 library and creating the Library
#

echo " Creating your Grib 1 Library ..."

tar -xf w3lib.tar
cp Makefile $LIBGRIBdir/src/.
cd $LIBGRIBdir/src/

rm -rf *.o
make
cp libw3.a $LIBGRIBdir/lib/.

#
## Now to compile the NetCDF2GRIB_V1.f converter code.
### NOTE:  Use the FORTRAN90 compiler that your system has.  Ask your 
###        System Administrator for this information if your unfamiliar
###        with this information.
#### NOTE: If you wish to compile version one which does grib files in 
####       in sigma levels then switch the names on line 60.  The default
####       setting is for version two which does grib files in pressure
####       levels
#

echo " Compiling NetCDF2GRIB_V4.f ..."

cd $CURRENTdir

echo " pgf90 -o NetCDF2GRIB_V4 NetCDF2GRIB_V4.f $LIBGRIB $LIBNETCDF"


pgf90 -g -o NetCDF2GRIB_V4 NetCDF2GRIB_V4.f $LIBGRIB $LIBNETCDF

echo " |============================| "
echo " | Completed the installation | "
echo " |============================| "







