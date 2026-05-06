subroutine SwashBCspecfile ( filenm, nbps, bpix, bpiy, xcgrid, ycgrid, kgrpnt, tcycl, lnest, btype, blayk, tsmo )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWASH team                               |
!   --|-----------------------------------------------------------|--
!
!
!     SWASH (Simulating WAves till SHore); a non-hydrostatic wave-flow model
!     Copyright (C) 2010-2026  Delft University of Technology
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program. If not, see <http://www.gnu.org/licenses/>.
!
!
!   Authors
!
!    1.00: Marcel Zijlema
!
!   Updates
!
!    1.00, July 2011: New subroutine
!
!   Purpose
!
!   Reads 1D or 2D spectrum from file, and makes it appropriate for synthesizing time series along open boundaries
!
!   Modules used
!
    use ocpcomm1
    use ocpcomm4
    use SwashCommdata3
    use SwashCommdata4
    use m_bndspec
    use m_parall
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                          :: blayk    ! indicate how boundary condition is specified for each layer
    integer, intent(in)                          :: btype    ! boundary type
    integer, intent(in)                          :: nbps     ! number of boundary grid points per segment / side
    !
    real   , intent(in)                          :: tcycl    ! cyclic period of time record
    real   , intent(in)                          :: tsmo     ! period for smoothing boundary values during cold start
    !
    logical, intent(in)                          :: lnest    ! indicates whether boundary conditions are obtained from nesting in SWAN
    !
    character(80), intent(in)                    :: filenm   ! file name
    !
    integer, dimension(mxc,myc), intent(in)      :: kgrpnt   ! index table containing the address of each (active) grid point
                                                             ! =1; not active grid point
                                                             ! >1; active grid point
    integer, dimension(nbps)   , intent(in)      :: bpix     ! indices of boundary points in x-direction
    integer, dimension(nbps)   , intent(in)      :: bpiy     ! indices of boundary points in y-direction
    !
    real   , dimension(mxc,myc), intent(in)      :: xcgrid   ! coordinates of computational grid in x-direction
    real   , dimension(mxc,myc), intent(in)      :: ycgrid   ! coordinates of computational grid in y-direction
!
!   Local variables
!
    integer                                      :: dorder   ! indicate the order of reading spectral directions
                                                             ! > 0: ascending order
                                                             ! < 0: descending order
    integer                                      :: i        ! loop counter
    integer, save                                :: ient = 0 ! number of entries in this subroutine
    integer                                      :: ierr     ! error indicator: ierr=0: no error, otherwise error
    integer                                      :: iostat   ! I/O status in call FOR
    integer                                      :: ival     ! auxiliary integer
    integer                                      :: ix       ! index of point in x-direction
    integer                                      :: iy       ! index of point in y-direction
    integer                                      :: j        ! loop counter
    integer                                      :: k        ! loop counter
    integer                                      :: mdc      ! number of SWAN directions in spectrum file
    integer                                      :: msc      ! number of SWAN frequencies in spectrum file
    integer                                      :: nbloc    ! number of boundary point locations
    integer, save                                :: nbv1 = 1 ! auxiliary integer to store nbv2 from previous variable part of segment / side
    integer                                      :: nbv2     ! user-defined number of boundary values
    integer                                      :: ndir     ! number of directions to be used for synthesizing time series
    integer                                      :: ndsd     ! unit reference number of spectrum file
    integer                                      :: nfreq    ! number of frequencies to be used for synthesizing time series
    integer                                      :: nhedt    ! number of heading lines per time step
    integer                                      :: nqua     ! number of quantities on spectrum file
    integer, dimension(3)                        :: switch   ! to switch to several options
    !
    real                                         :: ampl     ! amplitude of a wave component
    real                                         :: bfac     ! basic multiplication factor for spectral density
    real                                         :: cnorm    ! normalisation coefficient
    real                                         :: cosdir   ! cosine of directions
    real                                         :: dd       ! coefficient of directional spreading (standard deviation or power)
    real                                         :: ddir     ! increment in directional space
    real                                         :: df       ! increment in frequency space
    real                                         :: dir      ! spectral direction
    real                                         :: dsum     ! integral over normalised directional distribution
    real                                         :: dwidth   ! directional width
    real                                         :: edir     ! energy in a directional bin
    real                                         :: efac     ! multiplication factor to get proper variance density in m2/Hz/rad
    real                                         :: excval   ! exception value for spectral density
    real                                         :: f        ! frequency
    real                                         :: fac      ! auxiliary factor
    real                                         :: fmax     ! maximum frequency
    real                                         :: fmin     ! minimum frequency
    real                                         :: hs       ! significant wave height
    real                                         :: m0       ! moment of zeroth order
    real                                         :: ms       ! spreading coefficient for cosine-power
    real                                         :: rdist    ! relative distance from point on boundary to line between two consecutive locations on spectrum file
    real                                         :: rval     ! auxiliary real
    real                                         :: rval2    ! another auxiliary real
    real                                         :: rx       ! difference in x-coordinates between two consecutive locations on spectrum file
    real                                         :: ry       ! difference in y-coordinates between two consecutive locations on spectrum file
    real                                         :: sector   ! directional sector to be employed for 1D SWAN spectrum
    real                                         :: spec     ! spectral density
    real                                         :: specmax  ! maximum of 1D spectral density
    real                                         :: thp      ! peak wave direction
    real                                         :: x1       ! x-coordinate of previous location on spectrum file
    real                                         :: x2       ! x-coordinate of current location on spectrum file
    real                                         :: xp       ! x-coordinate of boundary point of segment or side
    real                                         :: y1       ! y-coordinate of previous location on spectrum file
    real                                         :: y2       ! y-coordinate of current location on spectrum file
    real                                         :: yp       ! y-coordinate of boundary point of segment or side
    !
    real, dimension(:,:), allocatable            :: bspden   ! spectral density on spectrum file
    real, dimension(:)  , allocatable            :: bspdir   ! spectral directions on spectrum file
    real, dimension(:)  , allocatable            :: bspfrq   ! spectral frequencies on spectrum file
    real, dimension(:)  , allocatable            :: cdf      ! cumulative distribution function
    real, dimension(:)  , allocatable            :: ddstr    ! non-normalised directional distribution
    real, dimension(:)  , allocatable            :: freq     ! uniformly distributed frequencies
    real, dimension(:)  , allocatable            :: spec1    ! auxiliary array to store 1D SWAN wave spectrum
    real, dimension(:)  , allocatable            :: spec2    ! auxiliary array to store 1D wave spectrum to be used for synthesizing time series
    real, dimension(:)  , allocatable            :: tarr     ! temporary array to store random values
    real, dimension(:,:), allocatable            :: varden   ! variance density to be used for synthesizing time series
    real, dimension(:)  , allocatable            :: xbp      ! x-coordinates of locations on spectrum file
    real, dimension(:)  , allocatable            :: ybp      ! y-coordinates of locations on spectrum file
    !
    real                                         :: GAMMAF   ! the gamma function
    !
    logical                                      :: EQCSTR   ! compares two strings
    logical                                      :: EQREAL   ! compares two reals
    logical                                      :: found    ! grid points of boundary found in between consecutvie locations on spectrum file
    logical                                      :: lloc     ! indicates whether there are locations on spectrum file or not
    logical                                      :: STPNOW   ! indicates that program must stop
    !
    character(80)                                :: hedlin   ! heading line
    character(120)                               :: msgstr   ! string to pass message
    !
    type(bfsdat), pointer                        :: bfstmp   ! list containing parameters for Fourier series
    type(bgpdat), pointer                        :: bgptmp   ! list containing parameters for boundary grid points
    !
    type fspt                                                ! linked list for Fourier components
       real                  :: ampl, omega
       type(fspt), pointer   :: nextfs
    end type fspt
    type(fspt), target       :: frstf
    type(fspt), pointer      :: currf, tmpf
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    call strace (ient,'SwashBCspecfile')
    !
    ! open file
    !
    ndsd   = 0
    iostat = 0
    call FOR (ndsd, filenm, 'OF', iostat)
    if (STPNOW()) return
    !
    ! imposed 1D or 2D spectrum
    !
    read (ndsd, '(a)', err=920) hedlin
    if ( EQCSTR(hedlin,'SPEC1D') ) then
       !
       ! file consisting of 2 columns: frequency, variance density
       !
       if ( lnest ) call msgerr (3, 'this spectrum file is not appropriate in case of command SPECSWAN')
       !
       df          = 0.
       nfreq       = 0
       frstf%ampl  = 0.
       frstf%omega = 0.
       nullify(frstf%nextfs)
       currf => frstf
       do
          read (ndsd, *, end=110, err=910) f, spec
          if ( spec < 0. ) then
             call msgerr (3, 'error in reading 1D spectrum file: density is negative')
             exit
          endif
          if ( nfreq < 2 ) df = f - df
          nfreq = nfreq + 1
          allocate(tmpf)
          tmpf%ampl  = sqrt(2.*spec)
          tmpf%omega = 2.*pi*f
          nullify(tmpf%nextfs)
          currf%nextfs => tmpf
          currf => tmpf
       enddo
       !
       ! close spectrum file
       !
 110   close (ndsd)
       !
       if ( tcycl < 9.e9 ) then
          write (msgstr,'(a,f9.2,a)') 'user-defined cyclic period is ignored. Instead, a cyclic period of ',1./df,' sec will be used.'
          call msgerr (1, trim(msgstr) )
       endif
       !
       ! store Fourier components
       !
       nbval = nbval + 1
       allocate(bfstmp)
       bfstmp%nbfs      = nbval
       bfstmp%nfreq     = nfreq
       bfstmp%azero     = -1.e10
       bfstmp%spparm(3) = -999.    ! indicent angle is assumed to be normal on the boundary
       !
       if ( nfreq > 0 ) then
          !
          currf => frstf%nextfs
          allocate(bfstmp%ampl (nfreq))
          allocate(bfstmp%omega(nfreq))
          allocate(bfstmp%phase(nfreq))
          allocate(bfstmp%theta(nfreq))
          !
          ! generate random values for phase
          allocate(tarr(nfreq))
          if ( INODE == MASTER ) then
             call random_seed(put=seedf)
             do j = 1, nfreq
                call random_number(rval)
                tarr(j) = rval
             enddo
          endif
          !
          ! scatter array tarr to all nodes
          !
          call SWBROADC ( tarr, nfreq, SWREAL )
          if (STPNOW()) return
          !
          do j = 1, nfreq
             bfstmp%ampl (j) = currf%ampl * sqrt(df)
             bfstmp%omega(j) = currf%omega
             bfstmp%phase(j) = 2.*pi*tarr(j)
             bfstmp%theta(j) = 0.            ! no directional spreading
             currf => currf%nextfs
          enddo
          deallocate(tmpf)
          deallocate(tarr)
          !
       endif
       !
       nullify(bfstmp%nextbfs)
       if ( .not.lbfs ) then
          fbfs = bfstmp
          cubfs => fbfs
          lbfs = .true.
       else
          cubfs%nextbfs => bfstmp
          cubfs => bfstmp
       endif
       !
    else if ( EQCSTR(hedlin,'SWAN') ) then
       !
       ! SWAN spectrum file
       !
 120   read (ndsd, '(a)', end=130, err=910) hedlin
       if ( ITEST >= 60 ) write (PRTEST, 201) hedlin
       ! skip heading lines starting with comment sign
       if ( hedlin(1:1) == COMID .or. hedlin(1:1) == '!' ) goto 120
       if ( EQCSTR(hedlin,'TIME') ) then
          call msgerr (3, 'nonstationary spectrum not allowed as boundary condition')
          read (ndsd, *, end=130, err=910) ival
          read (ndsd, '(a)', end=130, err=910) hedlin
          if ( ITEST >= 60 ) write (PRTEST, 201) hedlin
          nhedt = 1
       else
          nhedt = 0
       endif
       !
       ! read boundary locations
       !
       lloc = .true.
       if ( EQCSTR(hedlin,'LOC') ) then
          if ( lnest .and. kspher == 1 ) call msgerr (3, 'boundary locations in spectrum file are Cartesian, while computational ones are spherical')
       else if ( EQCSTR(hedlin,'LONLAT') ) then
          if ( lnest .and. kspher == 0 ) call msgerr (3, 'boundary locations in spectrum file are spherical, while computational ones are Cartesian')
       else
          lloc = .false.
       endif
       !
       if ( lloc ) then
          read (ndsd, *, end=130, err=910) nbloc
          if ( .not.lnest .and. nbloc >  1 ) call msgerr (2, 'spectrum file contains more than 1 boundary location, which is not appropriate in case of command SPECFILE')
          if (      lnest .and. nbloc == 1 ) call msgerr (3, 'spectrum file contains just 1 boundary location, which is not appropriate in case of command SPECSWAN')
          if (.not.allocated(xbp)) allocate(xbp(nbloc))
          if (.not.allocated(ybp)) allocate(ybp(nbloc))
          do j = 1, nbloc
             ierr = 0
             call REFIXY ( ndsd, xbp(j), ybp(j), ierr )
             if ( ierr == -1 ) goto 130
             if ( ierr == -2 ) goto 910
             if ( ITEST >= 80 ) write (PRTEST, *) ' coordinates of boundary locations in spectrum file ', j, xbp(j)+xoffs, ybp(j)+yoffs, ierr
          enddo
          !
          if ( .not.lnest ) nbloc = 1
          !
          read (ndsd, '(a)', end=130, err=910) hedlin
          if ( ITEST >= 60 ) write (PRTEST, 201) hedlin
       else
          if ( lnest ) call msgerr (3, 'spectrum file does not contain boundary locations, which is not appropriate in case of command SPECSWAN')
          nbloc = 1
          if (.not.allocated(xbp)) allocate(xbp(nbloc))
          if (.not.allocated(ybp)) allocate(ybp(nbloc))
          xbp = 0.
          ybp = 0.
       endif
       if ( nbloc > nbps ) call msgerr (1, 'there are more boundary locations than grid points on segment or side')
       if ( ITEST >= 60 ) write (PRTEST, 202) nbloc
       !
       ! read spectral resolution
       !
       if ( EQCSTR(hedlin(2:5),'FREQ') ) then
          read (ndsd, *, end=130, err=910) msc
          if (.not.allocated(bspfrq)) allocate(bspfrq(msc))
          do j = 1, msc
             ! read frequency in Hz and convert to rad/s
             read (ndsd, *, end=130, err=910) f
             bspfrq(j) = 2.*pi*f
          enddo
          read (ndsd, '(a)', end=130, err=910) hedlin
          if ( ITEST >= 60 ) write (PRTEST, 201) hedlin
       else
          msc = 0
          if (.not.allocated(bspfrq)) allocate(bspfrq(msc))
       endif
       if ( msc < 3 ) call msgerr (3, 'number of frequencies on SWAN spectrum file too small')
       if ( ITEST >= 60 ) write (PRTEST, 203) msc
       !
       if ( EQCSTR(hedlin(2:4),'DIR') ) then
          read (ndsd, *, end=130, err=910) mdc
          if (.not.allocated(bspdir)) allocate(bspdir(mdc))
          do j = 1, mdc
             ! read direction in degrees and convert to radians
             read (ndsd, *, end=130, err=910) dir
             if ( EQCSTR(hedlin,'N') ) dir = 180. + dnorth - dir
             dir = dir * degrad
             ! reverse order if second direction is smaller than first
             if ( j == 1 ) then
                dorder = 1
             else if ( j == 2 ) then
                if ( dir < rval ) then
                   dorder = -1
                   bspdir(mdc) = rval
                endif
             else
                if ( dorder < 0 ) then
                   if ( dir > rval ) call msgerr (3, 'spectral directions in SWAN spectrum file not in right order')
                else
                   if ( dir < rval ) call msgerr (3, 'spectral directions in SWAN spectrum file not in right order')
                endif
             endif
             rval = dir
             if ( dorder < 0 ) then
                bspdir(mdc+1-j) = dir
             else
                bspdir(j) = dir
             endif
          enddo
          ddir = abs(bspdir(2) - bspdir(1))
          ndir = mdc
          read (ndsd, '(a)', end=130, err=910) hedlin
          if ( ITEST >= 60 ) write (PRTEST, 201) hedlin
       else
          mdc    = 0
          ndir   = 51
          dorder = 0
          if (.not.allocated(bspdir)) allocate(bspdir(ndir))
       endif
       if ( oned .and. mdc /= 0 ) call msgerr (3, 'directional wave spectrum not allowed in 1D mode')
       if ( ITEST >= 60 ) write (PRTEST, 204) mdc
       !
       ! read characteristics of a quantity (name, unit, exception value)
       !
       if ( EQCSTR(hedlin,'QUANT') ) then
          read (ndsd, *, end=130, err=910) nqua
          if ( .not. ( (nqua == 1 .and. mdc > 0) .or. (nqua == 3 .and. mdc == 0) ) ) then
             call msgerr (2, 'incompatible data on SWAN spectrum file')
             write (PRINTF, 205) nqua, mdc
          endif
          do j = 1, nqua
             read (ndsd, '(a)', end=130, err=910) hedlin
             ! if first quantity is 'ENDENS' divide by rho*grav
             if ( j == 1 ) then
                if ( EQCSTR(hedlin,'ENDENS') ) then
                   ! quantity on file is energy density
                   switch(1) = 1
                else if ( EQCSTR(hedlin,'VADENS') ) then
                   ! quantity on file is variance density
                   switch(1) = 2
                else
                   call msgerr (2, 'incorrect quantity in SWAN spectrum file: '//hedlin(1:10))
                   switch(1) = 2
                endif
             else if ( j == 2 ) then
                ! if second quantity is 'NDIR' transform from Nautical to Cartesian direction
                if ( EQCSTR(hedlin,'NDIR') ) then
                   ! quantity on file is Nautical direction
                   switch(2) = 2
                else if ( EQCSTR(hedlin,'CDIR') ) then
                   ! quantity on file is Cartesian direction
                   switch(2) = 1
                else
                   call msgerr (2, 'incorrect quantity in SWAN spectrum file: '//hedlin(1:10))
                   switch(2) = 1
                endif
             else if ( j == 3 ) then
                ! if third quantity is 'DSPRP' or 'POWER' power is given, otherwise calculate power from directional spread in degrees
                if ( EQCSTR(hedlin,'DSPRP') .or. EQCSTR(hedlin,'POWER') ) then
                   ! quantity on file is power of cosine
                   switch(3) = 2
                else if ( EQCSTR(hedlin,'DSPR') .or. EQCSTR(hedlin,'DEGR') ) then
                   ! quantity on file is directional spread in degrees
                   switch(3) = 1
                else
                   call msgerr (2, 'incorrect quantity in SWAN spectrum file: '//hedlin(1:10))
                   switch(3) = 1
                endif
             endif
             ! check unit and exception value
             read (ndsd, '(a)', end=130, err=910) hedlin
             if ( j == 3 .and. EQCSTR(hedlin,'DEGR') ) then
                if ( switch(3) /= 1 ) then
                   call msgerr (2, 'incompatible options in SWAN spectrum file')
                   switch(3) = 1
                endif
             endif
             if ( j == 1 ) then
                read (ndsd, *, end=130, err=910) excval
             else
                read (ndsd, '(a)', end=130, err=910) hedlin
             endif
          enddo
       endif
       if ( ITEST >= 60 ) write (PRTEST, 206) nqua
       !
       ! determine frequency step
       !
       if ( tcycl < 9.e9 ) then
          df = 2.*pi/tcycl
       else
          df = bspfrq(2) - bspfrq(1)
       endif
       !
       ! determine frequency range and number of frequencies (uniformly distributed)
       !
       fmin  = bspfrq(1)
       fmax  = bspfrq(msc)
       nfreq = nint((fmax-fmin)/df) + 1
       !
       if (.not.allocated(freq)) allocate(freq(nfreq))
       !
       freq(1) = fmin
       do j = 2, nfreq
          freq(j) = freq(j-1) + df
       enddo
       !
       if (.not.allocated(bspden)) allocate(bspden(ndir,msc))
       if (.not.allocated(varden)) allocate(varden(ndir,nfreq))
       if (.not.allocated(ddstr )) allocate(ddstr (ndir))
       if (.not.allocated(cdf   )) allocate(cdf   (ndir))
       if (.not.allocated(spec1 )) allocate(spec1 (msc  ))
       if (.not.allocated(spec2 )) allocate(spec2 (nfreq))
       !
       ! multiply with Jacobian (= 1/2pi) to transform density E(f) with f in Hz into density E(omega) with omega in rad/s
       bfac = 1./(2.*pi)
       !
       ! multiply with 180/pi to account for directions in radians instead of degrees in case of 2D spectrum
       if ( mdc > 0 ) bfac = bfac * 180. / pi
       !
       ! divide by rhow*grav if quantity on file is energy density
       if ( switch(1) == 1 )  bfac = bfac / (rhow*grav)
       !
       ! read heading lines per time step
       !
       do j = 1, nhedt
          read (ndsd, '(a)', end=130, err=910) hedlin
          if ( ITEST >= 90 ) write (PRINTF, 201) hedlin
       enddo
       !
       ! read additional information from the header and spectrum for each boundary location
       !
       bsploop: do k = 1, nbloc
          !
          if ( k > nbps ) goto 130
          !
          ! read location number or ZERO/NODATA/FACTOR from heading line
          read (ndsd, '(a)', end=130, err=910) hedlin
          if ( ITEST >= 100 ) write (PRINTF, 201) hedlin
          !
          ! take proper action if heading line contains ZERO or NODATA
          if ( hedlin(1:6) == 'NODATA' .or. hedlin(1:4) == 'ZERO' ) then
             if ( nbloc == 1 ) then
                call msgerr (3, 'no data found on SWAN spectrum file')
             endif
             ! go to next location
             cycle bsploop
          else if ( hedlin(1:6) == 'FACTOR' ) then
             ! multiply spectral density with factor read from file
             read (ndsd, *, end=130, err=910) rval
             efac = bfac * rval
          else
             ! in case of 1D spectrum heading line can be ignored
             efac = bfac
             if ( mdc > 0 ) call msgerr (3, 'incorrect code in SWAN spectrum file: '//hedlin(1:20))
          endif
          !
          ! read spectral densities from SWAN spectrum file
          !
          if ( mdc == 0 ) then
             !
             ! 1D spectral input
             !
             specmax = -999.
             !
             do j = 1, msc
                !
                if ( j == 1 ) then
                   read (ndsd,*, end=130, err=910) spec1(j), rval, rval2
                else
                   read (ndsd,*, end=930, err=910) spec1(j), rval, rval2
                endif
                !
                ! in case of exception value, no energy
                if ( EQREAL(spec1(j),excval) ) then
                   spec1(j) = 0.
                   rval     = 0.
                   rval2    = 180./pi
                endif
                !
                ! search for peak wave energy
                if ( spec1(j) > specmax ) then
                   specmax = spec1(j)
                   thp     = rval
                   dd      = rval2
                endif
                !
             enddo
             !
             if ( switch(2) == 1 ) then
                ! conversion from degrees to radians
                thp = pi * thp / 180.
             else
                ! conversion from Nautical to Cartesian convention
                thp = pi * ( 180. + dnorth - thp ) / 180.
             endif
             !
             ! calculate spreading coefficient
             !
             if ( switch(3) == 1 ) then
                if ( dd > 23. ) then
                   fac = 1.2
                elseif ( dd > 17. ) then
                   fac = 1.096
                else
                   fac = 1.01
                endif
                dwidth = pi * dd / 180.
                if ( dwidth /= 0. ) then
                   ms = max(fac*dwidth**(-2) - 2., 1.)
                else
                   ms = 1000.
                endif
             else
                ms = dd
             endif
             !
             if ( ITEST >= 80 ) write (PRTEST, 207) specmax, 180.*thp/pi, ms
             !
             ! determine sector of directional space depending on spreading
             !
             if ( ms <= 6. ) then
                sector = pi
             else if ( ms <= 30. ) then
                sector = 0.5*pi
             else if ( ms <= 100. ) then
                sector = 0.25*pi
             else
                sector = 0.125*pi
             endif
             !
             ! determine directional space
             !
             ddir = sector / real(ndir-1)
             !
             do i = 1, ndir
                bspdir(i) = sector*(real(i-1)/real(ndir-1) - 0.5)
             enddo
             !
             ! compute normalisation coefficient of directional distribution (based on cos m-model)
             !
             if ( ms > 10. ) then
                cnorm = sqrt(ms/(2.*pi)) * (1. + 0.25/ms)
             else
                cnorm = GAMMAF(0.5*ms+1.) / (sqrt(pi) * GAMMAF(0.5*ms+0.5))
             endif
             !
             ! multiply non-directional wave spectrum with normalised directional distribution
             !
             dsum = 0.
             do i = 1, ndir
                cosdir = cos(bspdir(i))
                if ( cosdir > 0. ) then
                   edir = cnorm * max(cosdir**ms, 1.e-10)
                else
                   edir = 1.e-10
                endif
                if ( ITEST >= 20 ) dsum = dsum + edir * ddir
                do j = 1, msc
                   bspden(i,j) = edir * spec1(j)
                enddo
             enddo
             if ( ITEST >= 20 ) then
                if ( abs(dsum-1.) > 0.1 ) write (PRTEST, 208) dsum, cnorm, ms
             endif
             !
          else
             !
             ! 2D spectral input
             !
             if ( dorder < 0 ) then
                read (ndsd, *, end=930 ,err=910) ((bspden(i,j),i=mdc,1,-1),j=1,msc)
             else
                read (ndsd, *, end=930 ,err=910) ((bspden(i,j),i=1,mdc),j=1,msc)
             endif
             !
             if ( ITEST >= 120 ) then
                write (PRINTF,*)' spectrum from file'
                do j = 1, msc
                   write (PRINTF,*) j, (bspden(i,j),i=1,mdc)
                enddo
             endif
             !
             ! in case of exception value, no energy
             !
             do i = 1, mdc
                do j = 1, msc
                   if ( EQREAL(bspden(i,j),excval) ) bspden(i,j) = 0.
               enddo
             enddo
             !
             ! search for peak wave direction (use ddstr as a temporary array)
             ddstr = sum(bspden, dim=2)
             i     = maxval(maxloc(ddstr))
             thp   = bspdir(i)
             !
             bspdir = bspdir - thp
             !
          endif
          !
          ! compute two-dimensional variance density to be used for synthesizing time series
          !
          do i = 1, ndir
             !
             ! convert to variance density in m2/Hz/rad as function of radian frequency
             do j = 1, msc
                spec1(j) = efac * bspden(i,j)
             enddo
             !
             ! interpolate over uniformly distributed frequencies while keeping energy constant
             call CHGBAS (bspfrq, freq, 0., spec1, spec2, msc, nfreq, ITEST, PRTEST)
             !
             do j = 1, nfreq
                varden(i,j) = spec2(j)
             enddo
             !
          enddo
          !
          ! compute significant wave height
          !
          hs = 4.*sqrt(sum(varden)*df*ddir)
          !
          if ( ITEST >= 20 ) write(PRTEST, 209) xbp(k)+xoffs, ybp(k)+yoffs, hs
          !
          ! compute non-normalised directional distribution and its cumulative distribution function
          !
          ddstr = sum(varden, dim=2)
          do i = 1, ndir
             cdf(i) = sum(ddstr(1:i))
          enddo
          cdf = cdf/sum(ddstr)
          !
          ! store Fourier components
          !
          nbval = nbval + 1
          nbv2  = nbval
          allocate(bfstmp)
          bfstmp%nbfs      = nbval
          bfstmp%nfreq     = nfreq
          bfstmp%azero     = -1.e10
          if ( .not.oned ) then
             bfstmp%spparm(3) = thp
          else
             bfstmp%spparm(3) = -999.
          endif
          !
          if ( nfreq > 0 ) then
             !
             allocate(bfstmp%ampl (nfreq))
             allocate(bfstmp%omega(nfreq))
             allocate(bfstmp%phase(nfreq))
             allocate(bfstmp%theta(nfreq))
             !
             ! phase at each frequency is randomly chosen between 0 and 2pi
             !
             allocate(tarr(nfreq))
             if ( INODE == MASTER ) then
                call random_seed(put=seedf)
                do j = 1, nfreq
                   call random_number(rval)
                   tarr(j) = rval
                enddo
             endif
             call SWBROADC ( tarr, nfreq, SWREAL )
             if (STPNOW()) return
             !
             do j = 1, nfreq
                bfstmp%omega(j) = freq(j)
                bfstmp%phase(j) = 2.*pi*tarr(j)
             enddo
             !
             ! pick a random number in distribution function and assign corresponding (interpolated) direction to each frequency
             ! next, find the corresponding variance density and compute the amplitude
             !
             if ( INODE == MASTER ) then
                call random_seed(put=seedt)
                do j = 1, nfreq
                   call random_number(rval)
                   tarr(j) = rval
                enddo
             endif
             call SWBROADC ( tarr, nfreq, SWREAL )
             if (STPNOW()) return
             !
             m0 = 0.
             do j = 1, nfreq
                rval = tarr(j)
                if ( rval <= cdf(1) ) then
                   bfstmp%theta(j) = bspdir(1)
                   rval2           = varden(1,j)
                   ampl            = sqrt(2.*rval2*df*ddir)
                   bfstmp%ampl (j) = ampl
                else
                   do i = 2, ndir
                      if ( rval > cdf(i-1) .and. rval <= cdf(i) ) then
                         fac = (rval-cdf(i-1))/(cdf(i)-cdf(i-1))
                         bfstmp%theta(j) = (1.-fac)*bspdir(i-1)   + fac*bspdir(i)
                         rval2           = (1.-fac)*varden(i-1,j) + fac*varden(i,j)
                         ampl            = sqrt(2.*rval2*df*ddir)
                         bfstmp%ampl (j) = ampl
                      endif
                   enddo
                endif
                m0 = m0 + 0.5*ampl*ampl
             enddo
             !
             ! correct amplitude spectrum for wave height
             bfstmp%ampl = hs * bfstmp%ampl / ( 4.*sqrt(m0) )
             !
             deallocate(tarr)
             !
          endif
          !
          nullify(bfstmp%nextbfs)
          if ( .not.lbfs ) then
             fbfs = bfstmp
             cubfs => fbfs
             lbfs = .true.
          else
             cubfs%nextbfs => bfstmp
             cubfs => bfstmp
          endif
          !
          if ( mdc > 0 ) bspdir = bspdir + thp
          !
          if ( lnest ) then
             !
             found = .false.
             !
             x2 = xbp(k)
             y2 = ybp(k)
             !
             if ( k /= 1 ) then
                !
                rx   = x2 - x1
                ry   = y2 - y1
                rval = rx**2 + ry**2
                !
                if ( rval > 0. ) then
                   rx = rx/rval
                   ry = ry/rval
                   !
                   ! loop over boundary and select points between (x1,y1) and (x2,y2)
                   do i = 1, nbps
                      ix = bpix(i)
                      if ( optg /= 5 ) then
                         iy = bpiy(i)
                         xp = xcgrid(ix,iy)
                         yp = ycgrid(ix,iy)
                      else
                         xp = xcugrd(ix)
                         yp = ycugrd(ix)
                      endif
                      !
                      rdist = abs( rx*(yp-y1) - ry*(xp-x1) )
                      ! this is relative distance from current point to line between (x1,y1) and (x2,y2) with respect to the length of that line
                      !
                      if ( rdist < 0.1 ) then
                         !
                         fac = rx*(xp-x1) + ry*(yp-y1)
                         ! this is relative length of projection on line between (x1,y1) and (x2,y2)
                         !
                         if ( fac > -0.001 .and. fac < 1.001 ) then
                            !
                            found = .true.
                            !
                            if ( fac < 0.01 ) fac = 0.
                            if ( fac > 0.99 ) fac = 1.
                            !
                            allocate(bgptmp)
                            if ( optg /= 5 ) then
                                bgptmp%bgp(1) = kgrpnt(ix,iy)
                            else
                               bgptmp%bgp(1) = ix
                            endif
                            bgptmp%bgp(2) = btype
                            bgptmp%bgp(3) = nint(1000.*fac)
                            bgptmp%bgp(4) = nbv2
                            bgptmp%bgp(5) = nint(1000.*(1.-fac))
                            bgptmp%bgp(6) = nbv1
                            bgptmp%bgp(7) = nint(100000.*tsmo)
                            bgptmp%bgp(8) = blayk
                            bgptmp%bgp(9) = ibound
                            nullify(bgptmp%nextbgp)
                            if ( .not.lbgp ) then
                               fbgp = bgptmp
                               cubgp => fbgp
                               lbgp = .true.
                            else
                               cubgp%nextbgp => bgptmp
                               cubgp => bgptmp
                            endif
                            !
                         endif
                      endif
                   enddo
                endif
                !
                if ( .not.found ) then
                   write (msgstr,'(a,2f12.4,a,2f12.4)') 'no grid points on interval from ',x1+xoffs, y1+yoffs,' to ',x2+xoffs, y2+yoffs
                   call msgerr (1, trim(msgstr) )
                endif
             endif
             !
             x1   = x2
             y1   = y2
             nbv1 = nbv2
             !
          endif
          !
       enddo bsploop
       !
       ! close spectrum file
       !
 130   close(ndsd)
       !
       ! deallocate all temporary arrays
       !
       if ( allocated( xbp    ) ) deallocate ( xbp    )
       if ( allocated( ybp    ) ) deallocate ( ybp    )
       if ( allocated( bspfrq ) ) deallocate ( bspfrq )
       if ( allocated( bspdir ) ) deallocate ( bspdir )
       if ( allocated( freq   ) ) deallocate ( freq   )
       if ( allocated( bspden ) ) deallocate ( bspden )
       if ( allocated( varden ) ) deallocate ( varden )
       if ( allocated( ddstr  ) ) deallocate ( ddstr  )
       if ( allocated( cdf    ) ) deallocate ( cdf    )
       if ( allocated( spec1  ) ) deallocate ( spec1  )
       if ( allocated( spec2  ) ) deallocate ( spec2  )
       !
    else
       call msgerr (3, 'unsupported boundary data file')
    endif
    !
    return
    !
 201 format (' heading line: ', a)
 202 format (i6, ' boundary locations')
 203 format (i6, ' boundary frequencies')
 204 format (i6, ' boundary directions')
 205 format (i3, ' quantities; ', i5, ' directions')
 206 format (i6, ' boundary quantities')
 207 format (' SWAN peak density ', e10.3, '; cart dir ',f7.1, '; cos power ', f7.2)
 208 format (' integral over directions is ', f9.4,' with norm. coeff=', f10.3,'; power=', f8.2)
 209 format (' Hs of incident spectrum at location ', 2f12.4, ' equals ', f8.3)
    !
 910 call msgerr (4, 'error reading data from spectrum file '//filenm)
    return
 920 call msgerr (4, 'error opening spectrum file '//filenm)
    return
 930 call msgerr (2, 'insufficient data on spectrum file '//filenm)
    return
    !
end subroutine SwashBCspecfile
