#!/bin/csh

#==============================================================================
# This code is the will search for all netcdf data files in d01 and d02 domain
# assign it to $file and run the NetCDF2GRIB_V1 coverter to get grib ouputs in
# both domains.
#------------------------------------------------------------------------------
# @creator: Michael Kevin Hernandez
# @email  : mkh182@psu.edu
# @usage  : csh runme.csh
#==============================================================================

#===================================
# User defined variables
#===================================

set CURRENTdir      = `pwd`
set DATAdir         = /echidna/s0/mkh182/mk/ 

#
## This section of the code is specifically set for HRD
#

if ($#argv > 0) then
   set NAME  = $1
   set YEAR  = $2
   set MONTH = $3
   set DAY   = $4
   set HOUR  = $5
   set DATAdir         = /ptmp/$NAME.$YEAR$MONTH$DAY$HOUR/ATMOS
   echo $DATAdir
endif

cd $DATAdir
mkdir GRIBfiles

#
## Search for files with the prefex wrfout_d01_*:00 and assign it to $file_d01.
## Do the same for wrfout_d02_*:00.
#


set file_d01 = `ls wrfout_d01_*:00`
set file_d02 = `ls wrfout_d02_*:00`

@ count = 1
foreach processfile($file_d01)
        echo "  "
        $CURRENTdir/NetCDF2GRIB_V4 $processfile
        mv $processfile.gb GRIBfiles
        @ count ++
end 

echo " |===============================| "
echo " | Completed the Gribbing of d01 | "
echo " |===============================| "

@ count = 1
foreach processfile($file_d02)
        echo " "
        $CURRENTdir/NetCDF2GRIB_V4 $processfile
        mv $processfile.gb GRIBfiles
        @ count ++
end

cd $CURRENTdir

echo " |===============================| "
echo " | Completed the Gribbing of d02 | "
echo " |===============================| "







