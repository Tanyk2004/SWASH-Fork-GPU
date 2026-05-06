subroutine SwashInitCompGrid ( logcom )
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
!    1.00, February 2010: New subroutine
!
!   Purpose
!
!   Initialises arrays for description of computational grid in case of regular grid
!   These arrays need to be partitioned, if appropriate
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use m_genarr
    use m_parall
    use SwashFlowdata
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    logical, dimension(6), intent(inout) :: logcom   ! indicates which commands have been given to know if all the
                                                     ! information for a certain command is available. Meaning:
                                                     ! (1) no meaning
                                                     ! (2) command CGRID has been carried out
                                                     ! (3) command READINP BOTTOM has been carried out
                                                     ! (4) command READ COOR has been carried out
                                                     ! (5) command READ UNSTRUC has been carried out
                                                     ! (6) arrays s1, u1 and v1 have been allocated
!
!   Local variables
!
    integer                              :: i        ! loop counter
    integer, save                        :: ient = 0 ! number of entries in this subroutine
    integer                              :: istat    ! indicate status of allocation
    integer                              :: j        ! loop counter
    !
    real                                 :: dep      ! local bottom level
    real                                 :: SVALQI   ! interpolated value of an input array for point given in problem coordinates
    real                                 :: xp       ! x-coordinate of grid point
    real                                 :: yp       ! y-coordinate of grid point
    !
    logical                              :: EQREAL   ! compares two reals
    logical                              :: STPNOW   ! indicates that program must stop
    !
    character(80)                        :: msgstr   ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    call strace (ient,'SwashInitCompGrid')
    !
    ! fill index table of indirect addressing for computational grid in global domain
    !
    if ( .not.allocated(KGRPGL) ) allocate(KGRPGL(mxc,myc))
    KGRPGL = 1
    !
    do j = 1, myc-1
       do i = 1, mxc-1
          !
          if ( EQREAL(xcgrid(i,j), excfld(8)) ) then
             !
             if ( .not.EQREAL(ycgrid(i,j), excfld(9)) ) then
                call msgerr (2, 'incorrect grid coordinates')
                write (PRINTF, 201) xcgrid(i,j), ycgrid(i,j)
             endif
             !
          else
             !
             if ( EQREAL(ycgrid(i,j), excfld(9)) ) then
                call msgerr (2, 'incorrect grid coordinates')
                write (PRINTF, 201) xcgrid(i,j), ycgrid(i,j)
             else
                !
                ! compute coordinate offsets and reset grid coordinates in case of curvilinear grid
                !
                if ( optg == 3 ) then
                   !
                   if ( .not.lxoffs ) then
                      xoffs = xcgrid(i,j)
                      yoffs = ycgrid(i,j)
                      lxoffs = .true.
                      xcgrid(i,j) = 0.
                      ycgrid(i,j) = 0.
                   else
                      xcgrid(i,j) = real(xcgrid(i,j) - dble(xoffs))
                      ycgrid(i,j) = real(ycgrid(i,j) - dble(yoffs))
                   endif
                   !
                endif
                !
                xp = xcgrid(i,j)
                yp = ycgrid(i,j)
                !
                ! compute bottom level in corner point (xp,yp)
                !
                if ( lflgrd(1) ) then          ! bottom grid equals computational grid, no interpolation
                   !
                   dep = depth((j-1)*mxg(1)+i)
                   !
                else
                   !
                   dep = SVALQI (xp, yp, 1, depth, NEAREST)
                   !
                endif
                !
                if ( .not.EQREAL(dep, excfld(1)) ) then
                   !
                   mcgrd       = mcgrd + 1
                   KGRPGL(i,j) = mcgrd
                   !
                endif
                !
                if ( ITEST >= 250 .or. intes >= 30 ) write (PRINTF,202) i, j, xp+xoffs, yp+yoffs, dep, KGRPGL(i,j)
                !
             endif
             !
          endif
          !
       enddo
       !
       if ( KGRPGL(mxc-1,j) /= 1 ) then
          !
          mcgrd         = mcgrd + 1
          KGRPGL(mxc,j) = mcgrd
          !
       endif
       !
    enddo
    !
    if ( .not.oned ) then
       !
       do i = 1, mxc
          !
          if ( KGRPGL(i,myc-1) /= 1 ) then
             !
             mcgrd         = mcgrd + 1
             KGRPGL(i,myc) = mcgrd
             !
          endif
          !
       enddo
       !
    endif
    !
    excfld(8) = real(excfld(8) - dble(xoffs))
    excfld(9) = real(excfld(9) - dble(yoffs))
    !
    if ( mcgrd <= 1 ) call msgerr (3, 'no valid grid points found')
    !
    ! extrapolate virtual grid coordinates in case of curvilinear grid
    !
    if ( optg == 3 ) then
       !
       xcgrid(mxc,:) = 2.*xcgrid(mxc-1,:) - xcgrid(mxc-2,:)
       ycgrid(mxc,:) = 2.*ycgrid(mxc-1,:) - ycgrid(mxc-2,:)
       !
       xcgrid(:,myc) = 2.*xcgrid(:,myc-1) - xcgrid(:,myc-2)
       ycgrid(:,myc) = 2.*ycgrid(:,myc-1) - ycgrid(:,myc-2)
       !
    endif
    !
    ! check the curvilinear grid
    !
    if ( optg == 3 ) call CVCHEK ( KGRPGL, xcgrid, ycgrid )
    !
    ! compute xcgmin, xcgmax, ycgmin, ycgmax (of global domain)
    !
    xcgmin =  1.e9
    ycgmin =  1.e9
    xcgmax = -1.e9
    ycgmax = -1.e9
    do j = 1, myc-1
       do i = 1, mxc-1
          if ( KGRPGL(i,j) > 1 ) then
             if ( xcgrid(i,j) < xcgmin ) xcgmin = xcgrid(i,j)
             if ( ycgrid(i,j) < ycgmin ) ycgmin = ycgrid(i,j)
             if ( xcgrid(i,j) > xcgmax ) xcgmax = xcgrid(i,j)
             if ( ycgrid(i,j) > ycgmax ) ycgmax = ycgrid(i,j)
          endif
       enddo
    enddo
    !
    ! compute lengths and orientation of computational (global) domain in case of curvilinear grid
    ! note: it is assumed that the domain is rectangular with (non)uniform grid
    !       these parameters are needed for imposing wave spectrum on incoming boundaries
    !
    if ( optg == 3 ) then
       !
       xclen = sqrt( (xcgrid(mxc-1,1)-xcgrid(1,1))*(xcgrid(mxc-1,1)-xcgrid(1,1)) + (ycgrid(mxc-1,1)-ycgrid(1,1))*(ycgrid(mxc-1,1)-ycgrid(1,1)) )
       yclen = sqrt( (xcgrid(1,myc-1)-xcgrid(1,1))*(xcgrid(1,myc-1)-xcgrid(1,1)) + (ycgrid(1,myc-1)-ycgrid(1,1))*(ycgrid(1,myc-1)-ycgrid(1,1)) )
       !
       ! compute direction of the positive x-axis of the domain
       alpc = atan2(ycgrid(mxc-1,1)-ycgrid(1,1),xcgrid(mxc-1,1)-xcgrid(1,1))
       !
    endif
    !
    ! carry out domain decomposition meant for distributed-memory approach
    !
!TIMG    call SWTSTA(211)
    if ( PARLL ) then
       if ( kpart < 0 .or. kpart > 4 ) then
          if ( kpart > 4 ) call msgerr (1, 'incorrect grid partitioning option')
          if ( .not.oned .and. kmax > 1 ) then
             LORB = .true.
          else
             LORB = .false.
          endif
       else if ( kpart == 1 ) then
          LORB = .false.
       else if ( kpart == 2 ) then
          LORB = .true.
       else if ( kpart > 2 ) then
          LORB = .false.
       else
          LORB = .false.
          call msgerr (1, 'incorrect grid partitioning option')
       endif
    else
       LORB = .false.
    endif
    !
    allocate(INCORN(mxc,myc))
    INCORN = .false.
    !
    call SWDECOMP
!TIMG    call SWTSTO(211)
    if (STPNOW()) return
    !
!TIMG    call SWTSTA(212)
    !
    ! create copy of parts of KGRPGL for each subdomain -> kgrpnt
    !
    if ( .not.allocated(kgrpnt) ) allocate(kgrpnt(mxc,myc))
    kgrpnt = 1
    mcgrd  = 1
    do i = MXF, MXL
       do j = MYF, MYL
          if ( KGRPGL(i,j) /= 1 .and. .not.INCORN(i,j) ) then
             mcgrd                   = mcgrd + 1
             kgrpnt(i-MXF+1,j-MYF+1) = mcgrd
          endif
       enddo
    enddo
    deallocate(INCORN)
    !
    ! create copies of parts of xcgrid and ycgrid for each subdomain
    !
    if ( .not.allocated(XGRDGL) ) allocate(XGRDGL(MXCGL,MYCGL))
    if ( .not.allocated(YGRDGL) ) allocate(YGRDGL(MXCGL,MYCGL))
    XGRDGL = xcgrid
    YGRDGL = ycgrid
    deallocate(xcgrid,ycgrid)
    allocate(xcgrid(mxc,myc))
    allocate(ycgrid(mxc,myc))
    do i = MXF, MXL
       do j = MYF, MYL
          xcgrid(i-MXF+1,j-MYF+1) = XGRDGL(i,j)
          ycgrid(i-MXF+1,j-MYF+1) = YGRDGL(i,j)
       enddo
    enddo
    !
    ! compute XCLMIN, XCLMAX, YCLMIN, YCLMAX (of own subdomain)
    !
    XCLMIN =  1.e09
    YCLMIN =  1.e09
    XCLMAX = -1.e09
    YCLMAX = -1.e09
    do j = 1, myc
       do i = 1, mxc
          if ( kgrpnt(i,j) > 1 ) then
             if ( xcgrid(i,j) < XCLMIN ) XCLMIN = xcgrid(i,j)
             if ( ycgrid(i,j) < YCLMIN ) YCLMIN = ycgrid(i,j)
             if ( xcgrid(i,j) > XCLMAX ) XCLMAX = xcgrid(i,j)
             if ( ycgrid(i,j) > YCLMAX ) YCLMAX = ycgrid(i,j)
          endif
       enddo
    enddo
    !
!TIMG    call SWTSTO(212)
    !
    istat = 0
    if ( .not.allocated(s1) ) allocate(s1(mcgrd), stat = istat)
    if ( istat == 0 ) then
       if ( .not.allocated(u1) ) allocate (u1(mcgrd,kmax), stat = istat)
    endif
    if ( istat == 0 ) then
       if ( .not.allocated(v1) ) then
          if ( .not.oned ) then
             allocate (v1(mcgrd,kmax), stat = istat)
          else
             allocate (v1(0,0))
          endif
       endif
    endif
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: array s1, u1 or v1 and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    s1 = 0.
    u1 = 0.
    if ( .not.oned ) v1 = 0.
    logcom(6) = .true.
    !
    if ( kmax > 1 ) then
       !
       if ( .not.allocated(udep) ) allocate(udep(mcgrd), stat = istat)
       if ( istat == 0 ) then
          if ( .not.allocated(vdep) ) then
             if ( .not.oned ) then
                allocate (vdep(mcgrd), stat = istat)
             else
                allocate (vdep(0))
             endif
          endif
       endif
       if ( istat /= 0 ) then
          write (msgstr, '(a,i6)') 'allocation problem: array udep or vdep and return code is ',istat
          call msgerr ( 4, trim(msgstr) )
          return
       endif
       !
       udep = 0.
       if ( .not.oned ) vdep = 0.
       !
    endif
    !
    ! the following arrays for unstructured grids are allocated as empty ones
    !
    if ( .not.allocated(xcugrd) ) allocate(xcugrd(0))
    if ( .not.allocated(ycugrd) ) allocate(ycugrd(0))
    if ( .not.allocated( vmark) ) allocate( vmark(0))
    !
 201 format (' X= ', e12.4, '  Y= ', e12.4)
 202 format (2(i3,1x),1x,3(f10.1,1x),5x,i5)
    !
end subroutine SwashInitCompGrid
