subroutine SwashReadTestpnts ( lpntst, pntst )
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
!   Reads test points and generates point set TESTPNTS
!
!   Modules used
!
    use ocpcomm2
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata3
    use SwashCommdata4
    use outp_data
    use m_genarr
    use m_parall
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                     :: lpntst   ! maximum length of array pntst
    integer, dimension(lpntst), intent(out) :: pntst    ! grid point indices of test points
!
!   Local variables
!
    integer                                 :: i        ! loop counter
    integer, save                           :: ient = 0 ! number of entries in this subroutine
    integer                                 :: ilpos    ! actual length of error message filename
    integer                                 :: iostat   ! I/O status in call FOR
    integer                                 :: ix       ! index of point closest to test point in x-direction
    integer                                 :: iy       ! index of point closest to test point in y-direction
    integer                                 :: k        ! vertex index of point closest to test point
    integer                                 :: mptst    ! maximum number of test points
    !
    real                                    :: xc       ! broken index of test point in x-direction
    real                                    :: xp       ! x-coordinate of grid point
    real                                    :: yc       ! broken index of test point in y-direction
    real                                    :: yp       ! y-coordinate of grid point
    !
    logical                                 :: KEYWIS   ! indicates whether keyword in user manual is found or not
    logical                                 :: locgri   ! indicates whether test point is given as pair of indices or location
    logical                                 :: STPNOW   ! indicates that program must stop
    !
    type(opsdat), pointer                   :: opstmp   ! list containing parameters for grid
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    call strace (ient,'SwashReadTestpnts')
    !
    if ( optg /= 5 ) then
       mptst = lpntst/2
    else
       mptst = lpntst
    endif
    !
    call INKEYW ('STA','IJ')
    if ( mcgrd > 1 .or. nverts > 0 ) then
       if ( KEYWIS('XY') ) then
          locgri = .true.
       else if ( KEYWIS('IJ') ) then
          locgri = .false.
       else
          call WRNKEY
       endif
    else
       call msgerr (3, 'command READ BOT or READ UNSTRUC must precede this command')
    endif
    !
 100 if ( locgri ) then
       !
       call READXY ('X','Y', xp, yp, 'REP', -1.e10, -1.e10)
       if ( xp < -.9e10 ) then
          goto 130
       else
          if ( optg /= 5 ) then
             call CVMESH ( xp, yp, xc, yc, kgrpnt, xcgrid, ycgrid )
             if ( xc < 0. ) then
                if ( xp >= xcgmin .and. xp <= xcgmax .and. yp >= ycgmin .and. yp <= ycgmax ) then
                   goto 120
                else
                   goto 110
                endif
             endif
             ix = nint(xc) + MXF
             iy = nint(yc) + MYF
             if ( ITEST >= 30 ) write (PRTEST, 201) xp + xoffs, yp + yoffs, ix, iy
          else
             call SwanFindPoint ( xp, yp, k )
             if ( k < 0 ) then
                if ( xp >= xcgmin .and. xp <= xcgmax .and. yp >= ycgmin .and. yp <= ycgmax ) then
                   goto 120
                else
                   goto 110
                endif
             endif
             ix = k
             if ( ITEST >= 30 ) write (PRTEST, 202) xp + xoffs, yp + yoffs, ix
          endif
       endif
       !
    else
       !
       call ININTG ('I' , ix, 'REP', -1)
       if ( ix < 0 ) goto 130
       if ( optg /= 5 ) call ININTG ('J' , iy, 'REQ', 1)
       !
    endif
    !
    if ( optg /= 5 ) then
       !
       if ( ix >= 1 .and. ix <= MXCGL-1 .and. iy >= 1 .and. iy <= MYCGL-1 ) then
          ix = ix - MXF + 1
          iy = iy - MYF + 1
          if ( ix >= 1 .and. ix <= mxc-1 .and. iy >= 1 .and. iy <= myc-1 ) then
             if ( kgrpnt(ix,iy) > 1 ) then
                nptst = nptst + 1
                pntst(2*nptst-1) = ix
                pntst(2*nptst  ) = iy
                goto 120
             endif
          else
             goto 120
          endif
       endif
       !
    else
       !
       if ( ix >= 1 .and. ix <= nverts ) then
          nptst = nptst + 1
          pntst(nptst) = ix
          goto 120
       endif
       !
    endif
    !
 110 call msgerr (1, 'test point is not active')
    write (PRINTF, *) xp + xoffs, yp + yoffs
 120 if ( nptst <= mptst ) goto 100
    call msgerr (2, 'too many test points')
    !
    ! generate point set 'TESTPNTS' for test output
    !
 130 allocate(opstmp)
    opstmp%psname = 'TESTPNTS'
    opstmp%pstype = 'P'
    opstmp%mip    = nptst
    allocate(opstmp%xp(nptst))
    allocate(opstmp%yp(nptst))
    if ( optg /= 5 ) then
       do i = 1, nptst
          ix = pntst(2*i-1)
          iy = pntst(2*i)
          opstmp%xp(i) = xcgrid(ix,iy)
          opstmp%yp(i) = ycgrid(ix,iy)
       enddo
    else
       do i = 1, nptst
          k = pntst(i)
          opstmp%xp(i) = xcugrd(k)
          opstmp%yp(i) = ycugrd(k)
       enddo
    endif
    nullify(opstmp%nextops)
    if ( .not.lops ) then
       fops = opstmp
       cops => fops
       lops = .true.
    else
       cops%nextops => opstmp
       cops => opstmp
    endif
    !
    ! open file for test output of flow variables
    !
    call INCSTR ('FNAME', FILENM, 'STA', 'DIAGNOSTIC')
    ! append node number to FILENM in case of parallel computing
    if ( PARLL ) then
       ilpos = index ( FILENM, ' ' )-1
       write (FILENM(ilpos+1:ilpos+4), 203) INODE
    endif
    iostat = 0
    call FOR (ifpar, FILENM, 'UF', iostat)
    if (STPNOW()) return
    write (ifpar, 204) 1
    write (ifpar, 205) VERTXT
    write (ifpar, 206) PROJID, PROJNR
    if ( nstatm == 1 ) then
       write (ifpar, 207) 'TIME', 'time-dependent data'
       write (ifpar, 208) ITMOPT, 'time coding option'
    else
       write (ifpar, 207) 'ITER', 'iteration-dependent data'
       write (ifpar, 208) 0
    endif
    if ( kspher == 0 ) then
       write (ifpar, 207) 'LOCATIONS', 'locations in Cartesian coordinates (x,y)'
    else
       write (ifpar, 207) 'LONLAT', 'locations in spherical coordinates (lon,lat)'
    endif
    write (ifpar, 208) nptst, 'number of locations'
    if ( optg /= 5 ) then
       do i = 1, nptst
          ix = pntst(2*i-1)
          iy = pntst(2*i)
          write (ifpar, 209) xcgrid(ix,iy) + xoffs, ycgrid(ix,iy) + yoffs
       enddo
    else
       do i = 1, nptst
          k = pntst(i)
          write (ifpar, 209) xcugrd(k) + xoffs, ycugrd(k) + yoffs
       enddo
    endif
    !
    if ( oned .and. kmax == 1 ) then
       write (ifpar, 210) 2
    else if ( .not.oned .and. kmax > 1 ) then
       write (ifpar, 210) 4
    else
       write (ifpar, 210) 3
    endif
    !
    write (ifpar, 207) ovsnam(6), ovlnam(6)
    write (ifpar, 207) ovunit(6), 'unit'
    write (ifpar, 211) ovexcv(6), 'exception value'
    !
    write (ifpar, 207) ovsnam(11), ovlnam(11)
    write (ifpar, 207) ovunit(11), 'unit'
    write (ifpar, 211) ovexcv(11), 'exception value'
    !
    if ( .not.oned ) then
       !
       write (ifpar, 207) ovsnam(12), ovlnam(12)
       write (ifpar, 207) ovunit(12), 'unit'
       write (ifpar, 211) ovexcv(12), 'exception value'
       !
    endif
    !
    if ( kmax > 1 ) then
       !
       write (ifpar, 207) ovsnam(14), ovlnam(14)
       write (ifpar, 207) ovunit(14), 'unit'
       write (ifpar, 211) ovexcv(14), 'exception value'
       !
    endif
    !
 201 format (' test point ', 2f12.2, ' to grid point ', 2i4)
 202 format (' test point ', 2f12.2, ' to vertex ', i6)
 203 format ('-',i3.3)
 204 format ('SWASH', i4, t41,' diagnostic file')
 205 format ('$   Data produced by SWASH version ', a)
 206 format ('$   Project: ', a, ';  run number: ', a)
 207 format (a, t41, a)
 208 format (i6, t41, a)
 209 format (2(1x,f12.2))
 210 format ('QUANT', /, i6, t41, 'number of quantities in table')
 211 format (f14.6, t41, a)
    !
end subroutine SwashReadTestpnts
