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
C              temperature.
C    04/08/09: Michael Kevin Hernandez changed the code to write out
C              grib data from ARW, and interpolated it to pressure levels.
C              The interpolation scheme that was implemented is similar to
C              the bilinear vertical interpolation by Kevin Yeh from
C              the Hurricane Research Division. Potential temperature
C              calculations were deleted.
C    04/13/09: Michael kevin Hernandez fixed the U and V variables becuase
C              they were on C-grid data.  Now everything is in A-grid.
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
C ALBEDO   = Base Albedo
C PSFC     = Mean Sea Level Pressure
C P        = Pressure
C TEMP     = Temperature
C U        = U-wind component
C V        = V-wind component
C QVAPOR   = Specific Humidity
C LANDMASK = Sea and Land Mask
C TSK      = Skin Temperature
C GEOPOT  = PHI = Geopotential Hieght
C
C========================================================================
C @ program: NetCDF2GRIB.f
C @ author : Michael Kevin Hernandez, Erica Collura, and Shelly Guare
C @ email  : mkh182@psu.edu, elc5046@psu.edu, sag5085@psu.edu
C @ running: ./NetCDF2GRIB_V2.1 netcdfouputfile
C @ output : netcdfoutputfile.gb
C @ version: 1.0.0
C========================================================================
      PROGRAM NetCDF2GRIB

        IMPLICIT NONE
        
C
CC Must include the NetCDF Include file so that the program recognizes 
CC the NetCDF function and subroutine calls.
C

        include 'netcdf.inc'

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
        INTEGER*4 LUGB, KF, iret
        
C
CC For reading in variable
C

        INTEGER*4 nargs,IARGC, file1, file2, numbvar, modelhieghtphi
        INTEGER*4 numvars, londim, latdim, modelhieght, time, latdimV
        INTEGER*4 dim1, dim2, dim3, i, j, k, l, zed, numptsphi
        INTEGER*4 numpts, minpt, maxpt, minpt2, maxpt2, onelevel
        CHARACTER(LEN=256) ncfile, var, newname
        REAL*4, ALLOCATABLE, dimension(:) :: datum, pres, temp, pressure
        REAL*4, ALLOCATABLE, dimension(:) :: u_ave, v_ave
        REAL*4, ALLOCATABLE, dimension(:) :: presperturb, geopotperturb
        REAL*4, ALLOCATABLE, dimension(:,:,:):: prestran
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
      numbvar = 11
      ALLOCATE(varskeep(numbvar))

      varskeep(1)  = 'LANDMASK'
      varskeep(2)  = 'PSFC'
      varskeep(3)  = 'PB'
      varskeep(4)  = 'U'
      varskeep(5)  = 'V'
      varskeep(6)  = 'QVAPOR'
      varskeep(7)  = 'TSK'
      varskeep(8)  = 'ALBEDO'
      varskeep(9)  = 'P'
      varskeep(10) = 'PHB'
      varskeep(11) = 'PH'
  
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
              ALLOCATE(datum(numpts))
              CALL NCVGT(ncid, varid, start, count, datum, ncrcode)
C              datum = datum * 57.2957795  ! 1 radian = 57.2957795^o
              datum = datum * 1000        ! in millidegrees
              glatmax  = maxval(datum)
              glatmin  = minval(datum)
              DEALLOCATE(datum)

          else if (var == 'XLONG') then
              ALLOCATE(datum(numpts))
              CALL NCVGT(ncid, varid, start, count, datum, ncrcode)
C              datum    = datum * 57.2957795  ! 1 radian = 57.2958^o
              datum    = datum * 1000 ! positive lon in m^o
              glonmax  = maxval(datum)
              glonmin  = minval(datum)
              DEALLOCATE(datum)
              call define_kgds(londim, latdim, modelhieght,
     1                         glonmax, glatmax, glonmin, glatmin,
     2                         KGDS)
              goto 200
           endif 
       end do

200    write(*,'(A34)') '**********************************'
       write(*,'(A24)') 'Grid has been defined   '
       write(*,'(A34)') '** Commence bilinear p interp     '

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
                ALLOCATE(pres(numpts))
                ALLOCATE(prestran(londim,latdim,60))
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
                   DEALLOCATE(datum)
                end if
             end if
          end do
       end do

       call prestransform(numpts, modelhieght, pressure, londim, latdim,
     &                    prestran)

       write(*,'(A24)') '** Commence Gribbing    '
       write(*,'(A34)') '**********************************'

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
     &                         latdim, prestran)

                   DEALLOCATE(datum)
                else if (var == 'QVAPOR') then
                   datum = ((datum)/(1+datum) * 1000)
                   call ToGrib(numpts, modelhieght, datum, iyy, imo, 
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         latdim, prestran)

                   DEALLOCATE(datum)
                else if (var == 'P') then
                else if (var == 'PB') then
                else if (var == 'U') then
                   ALLOCATE(u_ave(numptsphi))
                   do j = 1, numptsphi, 1
                      u_ave(j) = (datum(j) + datum(j+1))/2
                   enddo
                   call ToGrib(numptsphi, modelhieght, u_ave, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, (londim - 1),
     &                         latdim, prestran)
                else if (var == 'V') then
                   ALLOCATE(v_ave(numptsphi))
                   do j = 1, numptsphi, 1
                      v_ave(j) = (datum(j) + datum(latdimV + j))/2
                   enddo
                   call ToGrib(numptsphi, modelhieght, v_ave, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         (latdim - 1), prestran)
                else
                   call ToGrib(numpts, modelhieght, datum, iyy, imo,
     &                         idd, ihh, imm, var, KGDS, londim,
     &                         latdim, prestran)


                   DEALLOCATE(datum)
                end if
             end if
          end do
       end do

C
CC The subroutine to calculate both geopotential and potential 
CC temperature.
C
      var = 'TEMP'
     
      call temps(pressure, geopot, numptsphi, modelhieghtphi, iyy, imo,
     &           idd, ihh, imm, KGDS, var, londim, latdim, prestran)

      call baclose(LUGB,iret)

      STOP      
      END

C========================================================================
C Subroutine temps
C========================================================================
C This subroutine calculates temperature which is not given by ARW.
C========================================================================
      subroutine temps(pres, phi,  numpts, modelhieght, iyy, imo, idd,
     &                  ihh, imm, KGDS, var, londim, latdim, prestran)


        implicit none

        INTEGER*4  numpts, onelevel, modelhieght, g, Rd, minpt, maxpt
        INTEGER*4  iyy, imo, idd, ihh, imm, minpt2, maxpt2, j
        INTEGER*4  londim, latdim
        CHARACTER(len=16) var
        INTEGER*4, dimension(200) :: KGDS
        REAL*4, dimension(numpts) :: temp, pres 
        REAL*4, dimension(numpts/modelhieght*61) ::  phi
        REAL*4, dimension(londim,latdim,60) :: prestran

        g = 9.81
        Rd = 287.15

     
        onelevel = numpts/modelhieght
 
C
CC The hydrostatic equation provides the algoritm to calculate this variable.
C

        do j = 1, modelhieght, 1                    ! modelhieght, 1
           minpt = ((j - 1) * onelevel) + 1
           minpt2 = (j*onelevel) + 1
           maxpt = (j * onelevel)
           maxpt2 = (j+1)*onelevel
           temp(minpt:maxpt)= ((phi(minpt2:maxpt2)-phi(minpt:maxpt))*g)/
     1                        (alog(pres(minpt:maxpt)/
     2                              pres(minpt2:maxpt2))*Rd) 
        end do



        call ToGrib(numpts, modelhieght, temp, iyy, imo, idd, ihh, imm,
     &              var, KGDS, londim, latdim, prestran)


      end subroutine

C========================================================================
C Subroutine prestransform
C========================================================================
C This subroutine recieves in the pres file, and creates a new pres file
C which will aid in vertically interpolating all other meteorological
C variables.
C========================================================================
      subroutine  prestransform(numpts, nz, pres, nx, ny, prestran)

        implicit none

        INTEGER*4  i, j, k, z, x, numpts, nx, ny, nz
        REAL*4,    dimension(nz*nx*ny):: pres
        REAL*4,    dimension(nx,ny,nz):: prestran

C
CC Here we seperate the one dimensional datum array into a 3-D array.
CC This will enable me to use Kevin Yeh's vertical interpolation code
CC to place the data into the much needed pressure level instead of
CC the current sigma level.
C

        do k = 0, nz - 1, 1
           do j = 0, ny - 1, 1
              do i = 1, nx, 1
                 x = i + (nx)*j + (nx*(ny))*k
                 prestran(i,(j + 1),(k + 1))=pres(x)
              enddo
           enddo
        enddo

      return

      end subroutine

C========================================================================
C Subroutine datatransform
C========================================================================
C This subroutine recieves in the data file, and creates a new data file
C which has been vertically interpolated into the common pressure levels,
C in the meteorological field.
C========================================================================
      subroutine  datatransform(numpts, nz, oneDarray, nx, ny, pm,
     &                          datumtransformed)

        implicit none

        INTEGER*4  i, j, k, z, x, numpts, nx, ny, nz
        INTEGER*4  nzp, ko, m, mb
        REAL*4,    deltaP, P
        REAL*4,    dimension(nz*nx*ny):: oneDarray, datumtransformed
        REAL*4,    dimension(nx*ny*60):: pres
        REAL*4,    dimension(nx,ny,nz):: threeDarray
        REAL*4,    dimension(nx,ny,41):: datumtwo, datumpres
        REAL*4,    dimension(nx,ny,60):: pm
        REAL*4,    dimension(41):: plev

C
CC Here we seperate the one dimensional datum array into a 3-D array.
CC This will enable me to use Kevin Yeh's vertical interpolation code
CC to place the data into the much needed pressure level instead of
CC the current sigma level.
C

        if (nz > 1) then
           do k = 0, nz - 1, 1
              do j = 0, ny - 1, 1
                 do i = 1, nx, 1
                    x = i + (nx)*j + (nx*(ny))*k
                    threeDarray(i,(j + 1),(k + 1))=oneDarray(x)
                 enddo
              enddo
           enddo

C
CC This is the p-interpolation code.
C

           nzp = 41

           do i = 1, nzp, 1
              plev(i) = 100000-2500*(i-1)
           end do
         
           do j = 1, ny, 1
              do i = 1, nx, 1
                 ko = 1
                 do m = 1, nzp, 1
                    mb = plev(m)
                    if (pm(i,j,1) >= mb .and. mb > pm(i,j,nz)) then
                       do k = ko, 60-1, 1
                          IF ( mb > pm(i,j,k+1) ) EXIT
                       enddo
                       deltaP= (mb-pm(i,j,k+1))/(pm(i,j,k)-pm(i,j,k+1))
                       P = 1 - (deltaP)
                       datumtwo(i,j,m) = deltaP * threeDarray(i,j,k) +
     &                                    P * threeDarray(i,j,k+1)
                       ko = k
                    endif
                 enddo
              enddo
           enddo

C
CC This is the part of the code that places the data back into a one
CC dimensional array to be fed into the ToGrib subroutine.
C

           do k = 0, nzp - 1, 1
              do j = 0, ny - 1, 1
                 do i = 1, nx, 1
                    x = i + (nx)*j + (nx*(ny))*k
                    datumtransformed(x) = datumtwo(i,(j + 1),(k + 1))
                 enddo
              enddo
           enddo
        end if
      
      return

      end subroutine

C========================================================================
C Subroutine ToGrib
C========================================================================
      subroutine ToGrib(numpts, modelhieght, datum, iyy, imo, idd, ihh, 
     &                  imm, var, KGDS, londim, latdim, prestran)

        implicit none

        INTEGER*4  onelevel, numpts, modelhieght, j, zed, KF, minpt
        INTEGER*4  maxpt, LUGB, iyy, imo, idd, ihh, imm, iret
        INTEGER*4  londim, latdim
        CHARACTER(len=16) var

        INTEGER*4, dimension(200)   :: KPDS, KGDS
        LOGICAL*1, allocatable      :: LB (:)
        REAL*4,    allocatable      :: F(:), datumtranformed(:)
        REAL*4,    dimension(numpts):: datum
        REAL*4,    dimension(londim,latdim,60):: prestran
        REAL*4,    dimension(41):: plev

        do j = 1, 41, 1
           plev(j) = 1000-25*(j-1)
        end do

        LUGB = 50
        onelevel = numpts/modelhieght
        KF = onelevel

        ALLOCATE(datumtranformed(onelevel*41))
        ALLOCATE (F(onelevel))
        ALLOCATE (LB(onelevel))

        write(*,*) '** Gribbing variable:  ',var


        if (modelhieght == 1) then

           LB = .false.
           zed = 1000

C           write(*,'(2A24)')    ' variable             :', var
C           write(*,'(A24,I10)') ' Current Height       :', zed

           call define_kpds(zed, iyy, imo, idd, ihh, imm, var, KPDS)

           F = datum

           CALL putgb(LUGB,KF,KPDS,KGDS,LB,F,iret)
C          CALL irets(iret)

        else
          call datatransform(numpts, modelhieght, datum, londim, latdim,
     &                       prestran, datumtranformed)

          LB = .false.

          do j = 1, 41, 1 
             zed = plev(j)


C           write(*,'(2A24)')    ' variable             :', var
C           write(*,'(A24,I10)') ' Current Height       :', zed

             call define_kpds(zed, iyy, imo, idd, ihh, imm, var, KPDS)
             minpt = ((j - 1) * onelevel) + 1
             maxpt = (j * onelevel)

             F = datumtranformed(minpt:maxpt)
             CALL putgb(LUGB,KF,KPDS,KGDS,LB,F,iret)
C            CALL irets(iret)

          enddo

        end if

210     continue
        DEALLOCATE(LB)
        DEALLOCATE(F)


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
        KPDS(6)  = 100  ! 100 for isobaric pressure levels
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
     &         cenlon

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

        dlat  = abs((glatmax-glatmin)/latdim)
        dlon  = abs((glonmax-glonmin)/londim)

C
CC  The below code is a manual hard coding of the KGDS array needed for
CC  this program to save the ascii file into a Grib file.
C

        KGDS(1)  = 0        ! regular lat lon
        KGDS(2)  = londim   ! Number of points on Latitude    
        KGDS(3)  = latdim   ! Number of points on Longitude
        KGDS(4)  = glatmin  ! latitude of first grid point
        KGDS(5)  = glonmin  ! longitude of first grid point
        KGDS(6)  = 136      ! Direction increments given, Earth assumed
                            ! spherical, reserved, reserved, u and v are
                            ! resolved relative to the defined direction 
                            ! of increasing x and y coordinates 
                            ! respectively
        KGDS(7)  = glatmax  ! latitude of extreme point
        KGDS(8)  = glonmax  ! longitude of extreme point
        KGDS(9)  = dlon     ! longitudinal direction of increment
        KGDS(10) = dlat     ! latitudinal direction increment
        KGDS(11) = 64       ! points scanned in +i direction, +j
                            ! direction, Adjacent points in i are
                            ! consecutive, reserved
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
         else if (var == 'GEOPOT') then
            kpds5 = 7           ! Geopotential [m]
         else if (var == 'TEMP') then
            kpds5 = 11          ! Temperature[K]
         else if (var == 'THETA') then
            kpds5 = 13          ! Potential tempertuare [k]
         else if (var == 'U') then
            kpds5 = 33          ! u-wind component [m/s]
         else if (var == 'V') then
            kpds5 = 34          ! v-wind component [m/s]
         else if (var == 'W') then
            kpds5 = 40          ! w-wind component [m/s]
         else if (var == 'QVAPOR') then
            kpds5 = 51          ! specific humidity[g/kg] 
         else if (var == 'ACPREC') then
            kpds5 = 61          ! total precip
         else if (var == 'LANDMASK') then
            kpds5 = 81          ! Land [=0] Sea [=1] Mask
         else if (var == 'ALBEDO') then
            kpds5 = 84          ! Base albedo
         else if (var == 'TSK') then
            kpds5 = 85          ! Deep ground soil temperature
         else
            kpds5 = 255         ! Missing value
         end if

         return

      endfunction
C========================================================================

