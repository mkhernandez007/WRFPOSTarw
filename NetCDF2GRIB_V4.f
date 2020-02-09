C========================================================================
C                        NetCDF to Grib 1 writer
C========================================================================
C
C  This program is to write data from the WRF-ARW netcdf output into
C  Grib 1 format.  The output file created by this program will have
C  the same name as the input netcdf file, but with the extension of .gb.
C  This will take will convert Pressure, u- and v- component wind,
C  temperature, specific humidity, sea-land mask, albedo, ground
C  roughness, skin temperature, and mean sea level pressure.  With this
C  information this code calculates geopotential height.
C
C  Updates:
C    10/31/08: Michael Kevin Hernandez wrote the code to read in an ascii
C              file and make a grib file as a result.
C    11/07/08: Erica Collura and Shelly Guare wrote the code to write out
C              the NetCDF data into the ascii file.
C    04/07/09: Michael Kevin Hernandez changed the code to write out
C              grib data from ARW. A subroutine was written to calculate
C              temperature using geopotential height.
C    04/13/09: Michael Kevin Hernandez fixed the U and V variables becuase
C              they were on C-grid data.  Now everything is in A-grid.
C    08/26/09: Michael Kevin Hernandez made the code to calculate tempe-
C              rature much easier (from potential temp, and not geopoten-
C              tial hieght).  The addition of Latent Heat flux, W compo-
C              nent of the wind, Relative humidity, geopotenial, vapor
C              pressure, and density, allows for more detailed analysis.
C    08/27/09: The addition of equivalent potential temperature, virtual
C              temperature. 
C    03/01/10: Discontinuities in wind measuremtns were found and delt with.
C
C  Note:
C    Without the aid of Erica's and Shelly's NetCDF to ascii, this
C    program will not have been accomplished in a timely matter.
C
C    This was written for the WRF-ARW model with a lambert projection
C    of the data points.
C
C------------------------------------------------------------------------
C List of argument
C------------------------------------------------------------------------
C
C PSFC     = Mean Sea Level Pressure
C P        = Pressure
C TEMP     = Temperature
C U        = U-wind component
C V        = V-wind component
C QVAPOR   = Specific Humidity * 1000
C LANDMASK = Sea and Land Mask
C TSK      = Skin Temperature
C GEOPOT   = PHI = Geopotential Hieght
C GEOPOTT  = Geopotential
C W        = W-wind component
C RAINNC   = Total precipitation
C RAINC    = Convective total precipitation
C T        = Potential Temperature
C TV       = Virtual Temperature
C TE       = Equivalent Potential Temperature
C LH       = Latent Heat Flux
C RH       = Relative Humidity
C VP       = Vapor Pressure
C RHO      = Density * 1000
C
C========================================================================
C @ program: NetCDF2GRIB_V4.f
C @ author : Michael Kevin Hernandez, Erica Collura, and Shelly Guare
C @ email  : mkh182@psu.edu, elc5046@psu.edu, sag5085@psu.edu
C @ running: ./NetCDF2GRIB_V4 netcdfouputfile
C @ output : netcdfoutputfile.gb
C @ version: 4.0.0
C========================================================================
      PROGRAM NetCDF2GRIB

        IMPLICIT NONE
        
C
CC Must include the NetCDF Include file so that the program recognizes 
CC the NetCDF function and subroutine calls.
C

       include '/usr/global/netcdf/netcdf-4.0.1-pgi9/include/netcdf.inc'

C
CC NetCDF usage varaibles
C

        INTEGER*4 ncid, ncrcode, ndims, nvars, ngatts, recdim, varid 
        INTEGER*4 start(9), vartype, nvdims, nvatts, vdims(9), ivc, ill
        INTEGER*4 dimsize(9), count(9), status, dimid

        CHARACTER(len=16),ALLOCATABLE, dimension(:) :: varskeep
        CHARACTER(len = 256) dimname(9)

C
CC GRIB usage variables
C

        INTEGER*4, dimension(200) :: KPDS, KGDS
        LOGICAL*1, allocatable :: LB (:)
        REAL*4,    allocatable :: F(:)
        REAL*4,    allocatable :: Ftmp(:)
        INTEGER*4 LUGB, KF, iret, fid, iexists
        
C
CC For reading in variable
C

        INTEGER*4 nargs,IARGC, file1, file2, numbvar, modelhieghtphi
        INTEGER*4 numvars, londim, latdim, modelhieght, time, latdimV
        INTEGER*4 dim1, dim2, dim3, i, j, k, l, m, zed, numptsphi, minpt
        INTEGER*4 numpts, maxpt, minpt2, maxpt2, onelevel, onelevel2
        CHARACTER(LEN=256) ncfile, var, newname
        REAL*4, ALLOCATABLE, dimension(:) :: datum, temp, pressure, ww
        REAL*4, ALLOCATABLE, dimension(:) :: u_ave, v_ave, theta
        REAL*4, ALLOCATABLE, dimension(:) :: presperturb, geopotperturb
        REAL*4, ALLOCATABLE, dimension(:,:,:):: prestran
        REAL*4, ALLOCATABLE, dimension(:) :: glat, glon
        REAL*4, ALLOCATABLE, dimension(:) :: geopot 
        REAL*4    glatmax, glatmin, glonmax, glonmin

C
CC In order to get the date of the file
C

        CHARACTER*4 yy, mo, dd, hh, mm
        INTEGER*4 iyy, imo, idd, ihh, imm

C
CC Variable you want to be gribbed.
C
      numbvar = 15
      ALLOCATE(varskeep(numbvar))

      varskeep(1)  = 'LANDMASK'
      varskeep(2)  = 'PSFC'
      varskeep(3)  = 'PB'
      varskeep(4)  = 'U'
      varskeep(5)  = 'V'
      varskeep(6)  = 'QVAPOR'
      varskeep(7)  = 'TSK'
      varskeep(8)  = 'P'
      varskeep(9)  = 'PHB'
      varskeep(10) = 'PH'
      varskeep(11) = 'W'
      varskeep(12) = 'RAINNC'
      varskeep(13) = 'RAINC'
      varskeep(14) = 'T'
      varskeep(15) = 'LH'
  
C
CC Retrieving a file
C 
        nargs = IARGC()
        if (nargs < 1) then
          write(*,*) 'You must enter an NetCDF file:'
          stop
        end if

        CALL GETARG(1, ncfile)

C
CC Retrieving the date of the forecast for which it is valid for.
C

        open(unit=130,file='date')
        yy = ncfile(12:15)
        mo = ncfile(17:18)
        dd = ncfile(20:21)
        hh = ncfile(23:24)
        mm = ncfile(26:27)
        write (130,'(5A7)') yy, mo, dd, hh, mm
        write (130,'(A30,A3)') ncfile,".gb"
        close (130)
        open(unit=130,file='date')
        read(130,'(5I7)')  iyy, imo, idd, ihh, imm
        read(130,'(A50)', end = 10) newname
10      write(*,'(A17,I3,A1,I2,A1,I4,A1,I2,A3)') " Gribbing date : ",  
     &       imo,"/",idd,"/",iyy,"/",ihh,":00"
        close(130)

C
CC Opening the NetCDF file.
C

        ncid = NCOPN(ncfile, NCNOWRIT, ncrcode)

C
CC Handle return if there is problems with NCOPN.
C

        if (ncrcode .ne. 0) then
           write(*,*) 'NCOPN return code:', ncrcode
           stop
        end if

C
CC NCINQ returns information about the open NetCDF file.
C

        CALL NCINQ(ncid, ndims, nvars, ngatts, recdim, ncrcode)

C
CC Handle return if problems occur from NCINQ.
C

        if (ncrcode .ne. 0) then
           stop
        else 
           continue
        end if

C
CC NCQINQ returns the name and size of a dimension.
CCC Loop though the seven dimenions of each variable extracting the
CCC dimension ID, name, and size. 
C
       
        do i = 1, ndims, 1
           dimid = i
           CALL NCDINQ(ncid, dimid, dimname(dimid), dimsize(dimid),
     &              ncrcode)
        end do

C
CC NCVINQ returns information about a NetCDF variable.  The information
CC returned is the name, type, number of dimensions, a list of dimension
CC IDs describing the shape of the variable, and the number of variable
CC attributes that have been assigned to the variable.
C

       do i=1, nvars, 1
          varid = i
          CALL NCVINQ (ncid, varid, var, vartype, nvdims, vdims,
     &                 nvatts, ncrcode)
       end do

C
CC Looping to set start equal to one each time program loops through for
CC each number of dimensions.
C

       do i=1, ndims, 1
          start(i) = 1
       end do
    
C
CC Open Grib1 file and setting LUGB
C

       LUGB = 50
       call baopenw(LUGB,newname,iret)

C
CC Looping to collect information to form the KGDS array.
C
     
       do i=1, nvars, 1
          varid = i
          CALL NCVINQ (ncid, varid, var, vartype, nvdims, vdims,
     &               nvatts, ncrcode)

C
CC Looping through count, it is the number of dimensions of the specified 
CC variable.  
C 

          count = 1
         
          do j=1, nvdims, 1
             count(j) = dimsize(vdims(j))
          end do
C
CC Setting count to their normal dimensions
C

          londim = count(1)
          latdim = count(2)

          numpts = londim * latdim 

          if (var == 'XLAT') then
              ALLOCATE(glat(numpts))
              CALL NCVGT(ncid, varid, start, count, glat, ncrcode)
              glat = glat * 1000        ! in millidegrees
              glatmin  = minval(glat)!glat(1)
              glatmax  = glat(londim*(latdim-1))!glat((latdim-1)*londim +1)!glat(londim*(latdim-1)+ latdim/2)

          else if (var == 'XLONG') then
              ALLOCATE(glon(numpts))
              CALL NCVGT(ncid, varid, start, count, glon, ncrcode)
              glon = glon * 1000        ! in millidegrees
              glonmax  = glon(londim)
              glonmin  = glon(1)
              
C
CC             Creating a simple near-neighbor interpolation map.
C

               DEALLOCATE(glat)
               DEALLOCATE(glon)

              call define_kgds(londim, latdim, modelhieght,
     1                         glonmax, glatmax, glonmin, glatmin,
     2                         KGDS)
              goto 200
           endif 
       end do
      

     

200    write(*,'(A34)') '**********************************'
       write(*,'(A24)') 'Grid has been defined   '
       write(*,'(A24)') '** Commence Gribbing    '
       write(*,'(A34)') '**********************************'
 
      call define_kgds(londim, latdim, modelhieght, glonmax, glatmax,
     &                 glonmin, glatmin, KGDS)

C
CC Looping to collect information about the pressure inorder to
CC begin the pressure interpolation.
C

        do i=1, nvars, 1

          varid = i

          CALL NCVINQ (ncid, varid, var, vartype, nvdims, vdims,
     &               nvatts, ncrcode)

          count = 1

          do j=1, nvdims, 1
             count(j) = dimsize(vdims(j))
          end do

C
CC Setting count to thier normal dimensions.
C
          londim = count(1)
          latdim = count(2)
          modelhieght = count(3)
          time =  count(4)
          dim1 =  count(5)
          dim2 =  count(6)
          dim3 =  count(7)

          numpts = londim * latdim * time * modelhieght

C
CCC NCVGT reads an array of values from NetCDF variable from the file.
CCC Looping to recognize each variable. The if statement recongnizes
CCC each of the variabels listed in varkeep(i).
C
          do k = 1, numbvar, 1
             if (var .eq. varskeep(k)) then
                ALLOCATE(datum(numpts))
                CALL NCVGT(ncid, varid, start, count, datum, ncrcode)
                if (var == 'P') then
                   ALLOCATE(presperturb(numpts))
                   presperturb = datum
                   DEALLOCATE(datum)
                else if(var == 'PB') then
                   var = 'PRES'
                   ALLOCATE(pressure(numpts))
                   numptsphi = numpts
                   modelhieghtphi = modelhieght
                   latdimV = latdim
                   pressure = (datum + presperturb)
                   call ToGrib(numpts, modelhieght, datum, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         latdim)
                   DEALLOCATE(datum)
                end if
             end if
          end do
       end do

C
CC Looping to collect information about each variable.
C

       do i=1, nvars, 1
          varid = i
          CALL NCVINQ (ncid, varid, var, vartype, nvdims, vdims,
     &               nvatts, ncrcode)

          count = 1

          do j=1, nvdims, 1
             count(j) = dimsize(vdims(j))
          end do

C
CC Setting count to their normal dimensions
C

          londim = count(1)
          latdim = count(2)
          modelhieght = count(3)
          time =  count(4)
          dim1 =  count(5)
          dim2 =  count(6)
          dim3 =  count(7)

          numpts = londim * latdim * time * modelhieght

C
CC NCVGT reads an array of values from NetCDF variable from the file.
CCC Looping to recognize each variable. The if statement recongnizes
CCC each of the variabels listed in varkeep(i).
C

          do k = 1, numbvar, 1
             if (var .eq. varskeep(k)) then
                ALLOCATE(datum(numpts))
                CALL NCVGT(ncid, varid, start, count, datum, ncrcode)
                if (var == 'PH') then
                   ALLOCATE(geopotperturb(numpts))
                   geopotperturb = datum
                   DEALLOCATE(datum)
                else if(var == 'PHB') then
                   var = 'GEOPOT'
                   ALLOCATE(geopot(numpts))
                   geopot = (datum + geopotperturb)/9.81
                   call ToGrib(numpts, modelhieght, geopot, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         latdim)
                   DEALLOCATE(datum)
                else if (var == 'QVAPOR') then
                   ALLOCATE(ww(numpts))
                   ww = datum
                   datum = ((datum)/(1+datum) * 1000)
                   call ToGrib(numpts, modelhieght, datum, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         latdim)
                   DEALLOCATE(datum)
                else if (var == 'W') then
                   datum = (datum * 10000)
                   call ToGrib(numpts, modelhieght, datum, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         latdim)
                   DEALLOCATE(datum)
                else if (var == 'T') then
                   datum = (datum + 300)
                   call ToGrib(numpts, modelhieght, datum, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         latdim)
                   Allocate(theta(numpts))
                   theta = datum
                   DEALLOCATE(datum)
                else if (var == 'P') then
                else if (var == 'PB') then
                else if (var == 'U') then
                   ALLOCATE(u_ave(numptsphi))
C                   ALLOCATE(u_ave(numpts))
                   onelevel = (londim-1)*latdim
                   onelevel2 = londim * latdim 
                   do m = 1, modelhieght, 1
                      do l = 1, latdim, 1
                         do j = 1, londim-1, 1
                          u_ave(((m-1)*onelevel)+((l-1)*(londim-1)+j))=
     &                   (datum(((m-1)*onelevel2)+((l-1)*(londim)+j))+
     &                    datum(((m-1)*onelevel2)+((l-1)*(londim)+(j+1)
     &                      )))/2
                          enddo
                       enddo
                   enddo
                   call ToGrib(numptsphi, modelhieght, u_ave, iyy, imo,
     &                         idd, ihh, imm, var, KGDS,
     &                         (londim - 1), latdim)
                   DEALLOCATE(u_ave)
                   DEALLOCATE(datum)
                else if (var == 'V') then
                   ALLOCATE(v_ave(numptsphi))
                   onelevel = (latdim-1)*londim 
                   onelevel2 = latdim*londim
                   do m = 1, modelhieght, 1
                      do l = 1, latdim -1, 1
                         do j = 1, londim, 1
                           v_ave(((m-1)*onelevel)+((l-1)*(londim)+j))=
     &                    (datum(((m-1)*onelevel2)+((l-1)*(londim)+j))+
     &                     datum(((m-1)*onelevel2)+((l-1)*(londim))
     &                     +(j+londim)))/2
                          enddo
                       enddo
                    enddo
                    
                   call ToGrib(numptsphi, modelhieght, v_ave, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         (latdim - 1))
                   DEALLOCATE(v_ave)
                   DEALLOCATE(datum)
                else
                   call ToGrib(numpts, modelhieght, datum, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         latdim)
                   DEALLOCATE(datum)
                end if
             end if
          end do
       end do

       DEALLOCATE(geopotperturb)

C
CC The subroutine to calculate both geopotential and potential 
CC temperature.
C

      DEALLOCATE(geopot)

      var = 'RH'
      call rhcalc(pressure, theta, numptsphi, modelhieghtphi, iyy,
     &            imo, idd, ihh, imm, KGDS, var, londim, latdim,
     &            ww)
      var = 'VP'
      call vpcalc(pressure, theta, numptsphi, modelhieghtphi, iyy,
     &            imo, idd, ihh, imm, KGDS, var, londim, latdim,
     &            ww)
      var = 'TV'
      call Tvcalc(pressure, theta, numptsphi, modelhieghtphi, iyy,
     &            imo, idd, ihh, imm, KGDS, var, londim, latdim,
     &            ww)
      var = 'TE'
      call Tecalc(pressure, theta, numptsphi, modelhieghtphi, iyy,
     &            imo, idd, ihh, imm, KGDS, var, londim, latdim,
     &            ww)

      DEALLOCATE(ww)

      var = 'RHO'
      call density(pressure, theta, numptsphi, modelhieghtphi, iyy, imo,
     &             idd, ihh, imm, KGDS, var, londim, latdim)

      var = 'TEMP'
      call temps(pressure, theta, numptsphi, modelhieghtphi, iyy, imo,
     &           idd, ihh, imm, KGDS, var, londim, latdim)

      call baclose(LUGB,iret)

      END

C========================================================================
C Subroutine Tecalc
C========================================================================
C This subroutine calculates Equivalent Potential Temperature which is
C not given by ARW.
C========================================================================
      subroutine Tecalc(pres, theta, numpts, modelhieght, iyy, imo, idd,
     &                 ihh, imm, KGDS, var, londim, latdim, ww)


        implicit none

        INTEGER*4  numpts, onelevel, modelhieght, g, Rd, minpt, maxpt
        INTEGER*4  iyy, imo, idd, ihh, imm, cp, j, lv
        INTEGER*4  londim, latdim
        CHARACTER(len=16) var
        INTEGER*4, dimension(200) :: KGDS
        REAL*4, dimension(numpts) :: pres, ww, theta, temp
        REAL*4, ALLOCATABLE :: Te(:)

        ALLOCATE(Te(numpts))

        g = 9.81
        Rd = 287.15
        cp = 1005.7
        lv = 2501000

        onelevel = numpts/modelhieght
        do j = 1, modelhieght, 1                    ! modelhieght, 1
           minpt = ((j - 1) * onelevel) + 1
           maxpt = (j * onelevel)
           temp(minpt:maxpt)= theta(minpt:maxpt)*(
     &                        (pres(minpt:maxpt)/100000)**(.286))
           Te(minpt:maxpt) = theta(minpt:maxpt)*
     1                       exp((lv*ww(minpt:maxpt))/
     2                       (cp*(temp(minpt:maxpt))))
        end do

        call ToGrib(numpts, modelhieght, Te, iyy, imo, idd, ihh, imm,
     &              var, KGDS, londim, latdim)

        DEALLOCATE(Te)

      end subroutine
C========================================================================
C Subroutine Tvcalc
C========================================================================
C This subroutine calculates virtual Temperature which is not given by
C ARW.
C========================================================================
      subroutine Tvcalc(pres, theta, numpts, modelhieght, iyy, imo, idd,
     &                 ihh, imm, KGDS, var, londim, latdim, ww)


        implicit none

        INTEGER*4  numpts, onelevel, modelhieght, g, Rd, minpt, maxpt
        INTEGER*4  iyy, imo, idd, ihh, imm, j
        INTEGER*4  londim, latdim
        CHARACTER(len=16) var
        INTEGER*4, dimension(200) :: KGDS
        REAL*4, dimension(numpts) :: pres, ww, theta
        REAL*4, ALLOCATABLE :: temp(:), Tv(:)

        ALLOCATE(Tv(numpts))
        ALLOCATE(temp(numpts))


        g = 9.81
        Rd = 287.15

        onelevel = numpts/modelhieght
        do j = 1, modelhieght, 1                    ! modelhieght, 1
           minpt = ((j - 1) * onelevel) + 1
           maxpt = (j * onelevel)
           temp(minpt:maxpt)= theta(minpt:maxpt)*(
     &                        (pres(minpt:maxpt)/100000)**(.286))
           Tv(minpt:maxpt) = (temp(minpt:maxpt))*(1+((461.5/287.15)-1)*
     &                           (ww(minpt:maxpt)))
        end do

        DEALLOCATE(temp)

        call ToGrib(numpts, modelhieght, Tv, iyy, imo, idd, ihh, imm,
     &              var, KGDS, londim, latdim)

        DEALLOCATE(Tv)

      end subroutine
C========================================================================
C Subroutine rhcalc
C========================================================================
C This subroutine calculates Relative Humidity which is not given by ARW.
C========================================================================
      subroutine rhcalc(pres, theta, numpts, modelhieght, iyy, imo, idd,
     &                 ihh, imm, KGDS, var, londim, latdim, ww)


        implicit none

        INTEGER*4  numpts, onelevel, modelhieght, g, Rd, minpt, maxpt
        INTEGER*4  iyy, imo, idd, ihh, imm, j
        INTEGER*4  londim, latdim
        CHARACTER(len=16) var
        INTEGER*4, dimension(200) :: KGDS
        REAL*4, dimension(numpts) :: pres, ww, theta
        REAL*4, ALLOCATABLE :: temp(:), rh(:), es(:), ws(:)

        ALLOCATE(temp(numpts))
        ALLOCATE(ws(numpts))
        ALLOCATE(es(numpts))
        ALLOCATE(rh(numpts))

        g = 9.81
        Rd = 287.15

        onelevel = numpts/modelhieght

        do j = 1, modelhieght, 1                    ! modelhieght, 1
           minpt = ((j - 1) * onelevel) + 1
           maxpt = (j * onelevel)
           temp(minpt:maxpt)= theta(minpt:maxpt)*(
     &                        (pres(minpt:maxpt)/100000)**(.286))-273.15
           es(minpt:maxpt)  = 6.112*exp((17.67 * temp(minpt:maxpt))/
     &                                  (243.5 + temp(minpt:maxpt)))
           ws(minpt:maxpt)   = (0.622*es(minpt:maxpt))/
     1                         (pres(minpt:maxpt)/100 - es(minpt:maxpt))
C     2                          0.378*temp(minpt:maxpt))
           rh(minpt:maxpt) = (ww(minpt:maxpt))/(ws(minpt:maxpt))*100
           end do

        DEALLOCATE(temp)
        DEALLOCATE(ws)
        DEALLOCATE(es)

        call ToGrib(numpts, modelhieght, rh, iyy, imo, idd, ihh, imm,
     &              var, KGDS, londim, latdim)

        DEALLOCATE(rh)

      end subroutine
C========================================================================
C Subroutine Density
C========================================================================
C This subroutine calculates density which is not given by ARW.
C========================================================================
      subroutine density(pres, theta,  numpts, modelhieght, iyy, imo,
     &                  idd, ihh, imm, KGDS, var, londim, latdim)


        implicit none

        INTEGER*4  numpts, onelevel, modelhieght, g, Rd, minpt, maxpt
        INTEGER*4  iyy, imo, idd, ihh, imm, j
        INTEGER*4  londim, latdim
        CHARACTER(len=16) var
        INTEGER*4, dimension(200) :: KGDS
        REAL*4, dimension(numpts) :: pres, theta, temp
        REAL*4, ALLOCATABLE :: rho(:)

        ALLOCATE(rho(numpts))

        g = 9.81
        Rd = 287.15

        onelevel = numpts/modelhieght

        do j = 1, modelhieght, 1                    ! modelhieght, 1
           minpt = ((j - 1) * onelevel) + 1
           maxpt = (j * onelevel)
           temp(minpt:maxpt)= theta(minpt:maxpt)*(
     &                        (pres(minpt:maxpt)/100000)**(.286))

           rho(minpt:maxpt)=  (pres(minpt:maxpt)/100)/(Rd*
     &                         temp(minpt:maxpt))*100000
        end do

        call ToGrib(numpts, modelhieght, rho, iyy, imo, idd, ihh, imm,
     &              var, KGDS, londim, latdim)

        DEALLOCATE(rho)

      end subroutine
C========================================================================
C Subroutine temps
C========================================================================
C This subroutine calculates temperature which is not given by ARW.
C========================================================================
      subroutine temps(pres, theta, numpts, modelhieght, iyy, imo, idd,
     &                 ihh, imm, KGDS, var, londim, latdim)


        implicit none

        INTEGER*4  numpts, onelevel, modelhieght, minpt, maxpt
        INTEGER*4  iyy, imo, idd, ihh, imm, j
        INTEGER*4  londim, latdim
        CHARACTER(len=16) var
        INTEGER*4, dimension(200) :: KGDS
        REAL*4, dimension(numpts) :: pres, theta
        REAL*4, ALLOCATABLE :: temp(:)

        ALLOCATE(temp(numpts))

        onelevel = numpts/modelhieght

        do j = 1, modelhieght, 1                    ! modelhieght, 1
           minpt = ((j - 1) * onelevel) + 1
           maxpt = (j * onelevel)
           temp(minpt:maxpt)= theta(minpt:maxpt)*(
     &                        (pres(minpt:maxpt)/100000)**(.286))
        end do

        call ToGrib(numpts, modelhieght, temp, iyy, imo, idd, ihh, imm,
     &              var, KGDS, londim, latdim)

        DEALLOCATE(temp)

      end subroutine
C========================================================================
C Subroutine vpcalc
C========================================================================
C This subroutine calculates vapor pressure which is not given by ARW.
C========================================================================
      subroutine vpcalc(pres, theta, numpts, modelhieght, iyy, imo, idd,
     &                 ihh, imm, KGDS, var, londim, latdim, ww)


        implicit none

        INTEGER*4  numpts, onelevel, modelhieght, g, Rd, minpt, maxpt
        INTEGER*4  iyy, imo, idd, ihh, imm, j
        INTEGER*4  londim, latdim
        CHARACTER(len=16) var
        INTEGER*4, dimension(200) :: KGDS
        REAL*4, dimension(numpts) :: pres, ww, theta
        REAL*4, ALLOCATABLE :: e(:), es(:), ws(:), rh(:), temp(:)

        ALLOCATE(e(numpts))
        ALLOCATE(temp(numpts))
        ALLOCATE(ws(numpts))
        ALLOCATE(es(numpts))
        ALLOCATE(rh(numpts))


        g = 9.81
        Rd = 287.15

        onelevel = numpts/modelhieght

        do j = 1, modelhieght, 1                    ! modelhieght, 1
           minpt = ((j - 1) * onelevel) + 1
           maxpt = (j * onelevel)
           temp(minpt:maxpt)= theta(minpt:maxpt)*(
     &                        (pres(minpt:maxpt)/100000)**(.286))-273.15
           es(minpt:maxpt)  = 6.112*exp((17.67 * temp(minpt:maxpt))/
     &                                  (243.5 + temp(minpt:maxpt)))
           ws(minpt:maxpt)   = (0.622*es(minpt:maxpt))/
     &                         (pres(minpt:maxpt)/100 - es(minpt:maxpt))
           rh(minpt:maxpt) = (ww(minpt:maxpt))/(ws(minpt:maxpt))
           e(minpt:maxpt)  = rh(minpt:maxpt) * es(minpt:maxpt) * 100
        end do

        DEALLOCATE(temp)
        DEALLOCATE(ws)
        DEALLOCATE(es)
        DEALLOCATE(rh)

        call ToGrib(numpts, modelhieght, e, iyy, imo, idd, ihh, imm,
     &              var, KGDS, londim, latdim)

        DEALLOCATE(e)

      end subroutine
C========================================================================
C Subroutine ToGrib
C========================================================================
      subroutine ToGrib(numpts, modelhieght, datum, iyy, imo, idd, ihh, 
     &                  imm, var, KGDS, londim, latdim)

        implicit none

        INTEGER*4  onelevel, numpts, modelhieght, j, zed, KF, minpt, i
        INTEGER*4  maxpt, LUGB, iyy, imo, idd, ihh, imm, iret
        INTEGER*4  londim, latdim
        CHARACTER(len=16) var

        INTEGER*4, dimension(200)  :: KPDS, KGDS
        LOGICAL*1, allocatable     :: LB (:)
        REAL*4,    allocatable     :: F(:),Finterp(:)
        REAL*4,    dimension(numpts):: datum

        LUGB = 50
        onelevel = numpts/modelhieght
        KF = onelevel

        ALLOCATE (F(onelevel))
        ALLOCATE (LB(onelevel))
        ALLOCATE (Finterp(onelevel))

        write(*,*) '** Gribbing variable:  ',var

        if (modelhieght .eq. 1) then

           LB = .false.
           zed = 1

C           write(*,'(2A24)')    ' variable             :', var
C           write(*,'(A24,I10)') ' Current Height       :', zed

           F = datum

           call define_kpds(zed, iyy, imo, idd, ihh, imm, var, KPDS)
   
           CALL putgb(LUGB,KF,KPDS,KGDS,LB,F,iret)
C           CALL irets(iret)

        else 
          LB = .false.

          do j = 1, modelhieght, 1 
           zed = j
           if (zed == 57) then
               goto 210
           end if
 
C           write(*,'(2A24)')    ' variable             :', var
C           write(*,'(A24,I10)') ' Current Height       :', zed

            call define_kpds(zed, iyy, imo, idd, ihh, imm, var, KPDS)

            minpt = ((j - 1) * onelevel) + 1
            maxpt = (j * onelevel)

            F = datum(minpt:maxpt)

 
            CALL putgb(LUGB,KF,KPDS,KGDS,LB,F,iret)
C            CALL irets(iret)

          enddo

        end if

210     continue
        DEALLOCATE(LB)
        DEALLOCATE(F)
        DEALLOCATE(Finterp)


      end subroutine

C========================================================================
C Subroutine define_kpds
C========================================================================
      subroutine define_kpds(zed, iyy, imo, idd, ihh, imm, var, kpds)

        implicit none      

        CHARACTER*256 var
        INTEGER*4, dimension(200) :: KPDS
        INTEGER*4 iyy, imo, idd, ihh, imm, zed, kpds5 

C
CC  The below code is a manual hard coding of the KPDS array needed for 
CC  this program to save the ascii file into a Grib file.
C

        KPDS(1)  = 51   ! Originating Center is Miami
        KPDS(2)  = 112  ! WRF-NMM, generic resolution (HRS is a derivative of
                        ! of WRF-NMM)
        KPDS(3)  = 255  ! non-standard grid - defined in the GDS 
        KPDS(4)  = 00000001   ! Section 2 included, Section 3 omitted
        KPDS(5)  = kpds5(var) ! function that assings units to variables
        KPDS(6)  = 107  ! 107 for sigma levels
        KPDS(7)  = zed  ! if KPDS(6) = 100 then zed which is the pressure 
                        ! at the isobaric pressure level
        KPDS(8)  = iyy  ! Year
        KPDS(9)  = imo  ! Month 
        KPDS(10) = idd  ! Day
        KPDS(11) = ihh  ! Hour 
        KPDS(12) = imm  ! Minute
        KPDS(13) = 10   ! Unit of time is 3 hour
        KPDS(14) = 0    ! Time Range 1
        KPDS(15) = 0    ! Time Range 2
        KPDS(16) = 0    ! Forecast product valid for reference time 
        KPDS(17) = 0    ! Number included in an average (no average)
        KPDS(18) = 1    ! Version of Grib  
        KPDS(19) = 2    ! GRIB parameter table version
        KPDS(20) = 0    ! Number missing from averages or accumulations
        KPDS(21) = 21   ! Current Century of data 
        KPDS(22) = 0    ! Unit scaling power of 10
        KPDS(23) = 0    ! Hurricane Research Division doesn't have a subcenter 

        return

      end subroutine
C========================================================================
C Subroutine irets
C========================================================================
      subroutine irets(iret)

        implicit none

        integer*4 iret

        if (iret == 0) then
           write(*,*) "COMPLETED MAKING GRIB FIELD WITHOUT ERROR"
        else if (iret == 1) then
           write(*,*) "IPFLAG NOT 0 OR 1"
        else if (iret == 2) then
           write(*,*) "IGFLAG NOT 0 OR 1"
        else if (iret == 3) then
           write(*,*) "ERROR CONVERTING IEEE F.P. NUMBER TO IBM370 F.P."
        else if (iret == 4) then
           write(*,*) "W3FI71 ERROR/IGRID NOT DEFINED"
        else if (iret == 5) then
           write(*,*) "W3FK74 ERROR/GRID REPRESENTATION TYPE NOT VALID"
        else if (iret == 6) then
           write(*,*) "GRID TOO LARGE FOR PACKER DIMENSION ARRAYS"
        else if (iret == 7) then
           write(*,*) "LENGTH OF BIT MAP NOT EQUAL TO SIZE OF FLD/IFLD"
        else if (iret == 8) then
           write(*,*) "W3FI73 ERROR, ALL VALUES IN IBMAP ARE ZERO"
        end if
  
        return

      end subroutine
C========================================================================
C Subroutine define_kgds
C========================================================================
      subroutine define_kgds(londim, latdim, modelhieght, glonmax,
     &                       glatmax, glonmin, glatmin, kgds)

        implicit none

        CHARACTER*40 var, asciifilename
        INTEGER*4, dimension(200) :: KGDS
        INTEGER*4 londim, latdim, modelhieght, total
        REAL*4 glonmax, glatmax, glonmin, glatmin, dlat, dlon, cenlat, 
     &         cenlon, latone, lattwo, lontrue, glongmin, glongmax

C
CC Total number of actual points
C   

        total = londim * latdim

C
CC Calculating the central lat and lon in m^o
C

        cenlon = (-72.751+360) * 1000
        cenlat = 21.0 * 1000

C
CC Resolution in m^o degrees
C

        dlat  = abs((glatmax-glatmin)/latdim * 100)
        dlon  = abs((glonmax-glonmin)/londim * 100)
C
CC Lat from pole of secant
C
        latone = glatmax!abs(glatmax - glatmin)/2 + glatmin + 15000
        lattwo = glatmin!abs(glatmax - glatmin)/2 + glatmin - 15000

C
CC True longitude
C

        lontrue = glonmin + (londim/2 * dlon)/100
        

C
CC  The below code is a manual hard coding of the KGDS array needed for
CC  this program to save the ascii file into a Grib file.
C

        KGDS(1)  = 3        ! lambert projection
        KGDS(2)  = londim   ! Number of points on Latitude    
        KGDS(3)  = latdim   ! Number of points on Longitude
        KGDS(4)  = glatmin  ! latitude of first grid point
        KGDS(5)  = glonmin  ! longitude of first grid point
        KGDS(6)  = 136      ! Direction increments given, Earth assumed
                            ! spherical, reserved, reserved, u and v are
                            ! resolved relative to the defined direction 
                            ! of increasing x and y coordinates 
                            ! respectively
        KGDS(7)  = lontrue  ! LOV - Orientation of Grid
        KGDS(8)  = dlon     ! longitudinal direction of increment
        KGDS(9)  = dlat     ! latitudinal direction increment
        KGDS(10) = 0        ! projection center (North pole is on the 
                            ! plane, if it were the South pole this
                            ! would equal 1)
        KGDS(11) = 64       ! points scanned in +i direction, +j
                            ! direction, Adjacent points in i are
                            ! consecutive, reserved
        KGDS(12) = latone   ! First lat from pole of secant cone
        KGDS(13) = lattwo   ! Second lat from pole of secant cone
        KGDS(19) = 0        ! number for vertical coordinate parameters
        KGDS(20) = 255      ! no octet number of the list of vertical 
                            ! coordinate parameters
        KGDS(21) = 0        ! no PL
        KGDS(22) = 0        ! number of words in each row

        return

      end subroutine
C========================================================================
C Subroutine kpds5
C========================================================================
      integer function kpds5(var)
      
         implicit none 

         character*256 var
   
         if (var == 'PRES') then
            kpds5 = 1           ! Pressure [Pa]
         else if (var == 'PSFC') then
            kpds5 = 2           ! Pressure reduce to MSL [Pa]
         else if (var == 'GEOPOTT') then
            kpds5 = 6           ! Geopotential [m&2/s^2]
         else if (var == 'GEOPOT') then
            kpds5 = 7           ! Geopotential height [gpm]
         else if (var == 'TEMP') then
            kpds5 = 11          ! Temperature[K]
         else if (var == 'TV') then
            kpds5 = 12          !Virtual Temperature [K]
         else if (var == 'T') then
            kpds5 = 13          ! Potential tempertuare [K]
         else if (var == 'TE') then
            kpds5 = 14          ! Equivalent Potential tempertuare [K]
         else if (var == 'TD') then
            kpds5 = 17          ! Dew Point Temperature [K]
         else if (var == 'U') then
            kpds5 = 33          ! u-wind component [m/s]
         else if (var == 'V') then
            kpds5 = 34          ! v-wind component [m/s]
         else if (var == 'W') then
            kpds5 = 40          ! w-wind component [m/s]
         else if (var == 'QVAPOR') then
            kpds5 = 51          ! specific humidity[g/kg] 
         else if (var == 'RH') then
            kpds5 = 52          ! Relative Humidity[%]
         else if (var == 'VP') then
            kpds5 = 55          ! vapor pressure [pa]
         else if (var == 'RAINNC') then
            kpds5 = 61          ! total precip [kg/m^2]
         else if (var == 'RAINC') then
            kpds5 = 63          ! convective total precip [kg/m^2]
         else if (var == 'SNOW') then
            kpds5 = 65          ! Snow Depth [m]
         else if (var == 'SNOWH') then
            kpds5 = 66          ! Water equiv. of accum. snow depth [kg/m^2]
         else if (var == 'LANDMASK') then
            kpds5 = 81          ! Land [=0] Sea [=1] Mask
         else if (var == 'ALBEDO') then
            kpds5 = 84          ! Base albedo [%]
         else if (var == 'TSK') then
            kpds5 = 85          ! Deep ground soil temperature [K]
         else if (var == 'VEGFRA') then
            kpds5 = 87          ! Vegetation [%]
         else if (var == 'RHO') then
            kpds5 = 89          ! Density [kg/m^2]
         else if (var == 'LH') then
            kpds5 = 121         ! Latent heat Flux [W/m^2]
         else
            kpds5 = 255         ! Missing value
         end if

         return

      endfunction
C========================================================================



