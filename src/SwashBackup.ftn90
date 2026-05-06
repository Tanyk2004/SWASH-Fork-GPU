subroutine SwashBackup
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
!   Backups current state of the flow variables to a file
!
!   Modules used
!
    use ocpcomm2
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use SwashTimecomm
    use SwashFlowdata
    use m_genarr
    use m_parall
    use SwanGriddata
!
    implicit none
!
!   Local variables
!
    integer                 :: i                   ! loop counter
    integer, save           :: ient = 0            ! number of entries in this subroutine
    integer                 :: ilpos               ! actual length of error message filename
    integer                 :: indx                ! point index
    integer                 :: iostat              ! I/O status in call FOR
    integer                 :: j                   ! loop counter
    integer                 :: k                   ! loop counter
    integer                 :: l                   ! counter
    integer                 :: nref                ! reference number
    !
    real                    :: efac                ! reference factor for flow variables
    !
    logical                 :: EQREAL              ! compares two reals
    logical                 :: STPNOW              ! indicates that program must stop
    !
    character(8)            :: crform = '(2f14.4)' ! string format
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashBackup')
    !
    ! ==========================================================================
    !
    ! HOTFile 'fname'
    !
    ! ==========================================================================
    !
    call INCSTR ('FNAME', FILENM, 'REQ', ' ')
    ! append node number to FILENM in case of parallel computing
    if ( PARLL ) then
       ilpos = index ( FILENM, ' ' )-1
       write (FILENM(ilpos+1:ilpos+4), 201) INODE
    endif
    nref   = 0
    iostat = 0
    call FOR (nref, FILENM, 'UF', iostat)
    if (STPNOW()) return
    write (nref, 202) 'SWASH   1', 'SWASH standard file, version'
    if ( nstatm == 1 ) then
       write (nref, 202) 'TIME', 'time-dependent data'
       write (nref, 203) ITMOPT, 'time coding option'
    endif
    if ( kspher == 0 ) then
       write (nref, 202) 'LOCATIONS', 'locations in Cartesian coordinates (x,y)'
       crform = '(2f14.4)'
    else
       write (nref, 202) 'LONLAT', 'locations in spherical coordinates (lon,lat)'
       crform = '(2f12.6)'
    endif
    if ( optg /= 5 ) then
       write (nref, 203) (mxc-1)*(myc-1), 'number of locations'
       do i = 1, mxc-1
          do j = 1, myc-1
             if ( EQREAL(xcgrid(i,j), excfld(8)) .and. EQREAL(ycgrid(i,j), excfld(9)) ) then
                write (nref, fmt=crform) dble(excfld(8)) + dble(xoffs), dble(excfld(9)) + dble(yoffs)
             else
                write (nref, fmt=crform) dble(xcgrid(i,j)) + dble(xoffs), dble(ycgrid(i,j)) + dble(yoffs)
             endif
          enddo
       enddo
    else
       write (nref, 203) nverts, 'number of locations'
       do k = 1, nverts
          write (nref, fmt=crform) dble(xcugrd(k)) + dble(xoffs), dble(ycugrd(k)) + dble(yoffs)
       enddo
    endif
    !
    if ( kmax > 1 ) then
       !
       write (nref, 202) 'ZK', 'vertical layer schematization'
       write (nref, 203) kmax, 'number of layers'
       do k = 1, kmax
          write (nref, '(i2,1pe13.6e2)') indlay(k), hlay(k)
       enddo
       !
    endif
    !
    if ( oned ) then
       if ( kmax == 1 ) then
          write (nref, 204) 3
       else
          write (nref, 204) 5
       endif
    else
       if ( kmax == 1 ) then
          write (nref, 204) 4
       else
          write (nref, 204) 7
       endif
    endif
    !
    write (nref, 202) ovsnam(5), ovlnam(5)
    write (nref, 202) ovunit(5), 'unit'
    write (nref, 205) ovexcv(5), 'exception value'
    !
    write (nref, 202) ovsnam(6), ovlnam(6)
    write (nref, 202) ovunit(6), 'unit'
    write (nref, 205) ovexcv(6), 'exception value'
    !
    write (nref, 202) ovsnam(10), ovlnam(10)
    write (nref, 202) ovunit(10), 'unit'
    write (nref, 205) ovexcv(10), 'exception value'
    !
    if ( kmax > 1 ) then
       !
       write (nref, 202) ovsnam(10), trim(ovlnam(10))//' per layer'
       write (nref, 202) ovunit(10), 'unit'
       write (nref, 205) ovexcv(10), 'exception value'
       !
    endif
    !
    if ( .not.oned ) then
       !
       write (nref, 202) ovsnam(11), ovlnam(11)
       write (nref, 202) ovunit(11), 'unit'
       write (nref, 205) ovexcv(11), 'exception value'
       !
       if ( kmax > 1 ) then
          !
          write (nref, 202) ovsnam(11), trim(ovlnam(11))//' per layer'
          write (nref, 202) ovunit(11), 'unit'
          write (nref, 205) ovexcv(11), 'exception value'
          !
       endif
       !
    endif
    !
    if ( kmax > 1 ) then
       !
       write (nref, 202) ovsnam(14), ovlnam(14)
       write (nref, 202) ovunit(14), 'unit'
       write (nref, 205) ovexcv(14), 'exception value'
       !
    endif
    !
    ! writing of heading is completed, write time if nonstationary
    !
    if ( nstatm == 1 ) then
       if ( ITMOPT /= 7 ) then
          write (nref, 206) chtime
       else
          write (nref, 207) chtime
       endif
    endif
    !
    ! write flow variables
    !
    if ( optg /= 5 ) then
       !
       l = 0
       do i = 1, mxc-1
          do j = 1, myc-1
             !
             l = l + 1
             indx = kgrpnt(i,j)
             !
             if ( indx == 1 ) then
                !
                write (nref, 208) 'NODATA'
                !
             else
                !
                write (nref, 209) l
                !
                ! bottom level
                !
                write (nref, 210) dps(indx)
                !
                ! water level
                !
                write (nref, 210) s1(indx)
                !
                ! depth-averaged u-velocity
                !
                if ( kmax == 1 ) then
                   !
                   write (nref, 210) u1(indx,1)
                   !
                else
                   !
                   write (nref, 210) udep(indx)
                   !
                endif
                !
                ! u-velocity per layer
                !
                if ( kmax > 1 ) then
                   !
                   efac = maxval(abs(u1(indx,:)))
                   if ( efac < 1.e-8 ) then
                      write (nref, 211) 'ZERO'
                   else
                      efac = efac * 10.**(-5)
                      write (nref, 212) efac
                      write (nref, '(10(1x,i7))') (nint(u1(indx,k)/efac), k=1,kmax)
                   endif
                   !
                endif
                !
                if ( .not.oned ) then
                   !
                   ! depth-averaged v-velocity
                   !
                   if ( kmax == 1 ) then
                      !
                      write (nref, 210) v1(indx,1)
                      !
                   else
                      !
                      write (nref, 210) vdep(indx)
                      !
                   endif
                   !
                   ! v-velocity per layer
                   !
                   if ( kmax > 1 ) then
                      !
                      efac = maxval(abs(v1(indx,:)))
                      if ( efac < 1.e-8 ) then
                         write (nref, 211) 'ZERO'
                      else
                         efac = efac * 10.**(-5)
                         write (nref, 212) efac
                         write (nref, '(10(1x,i7))') (nint(v1(indx,k)/efac), k=1,kmax)
                      endif
                      !
                   endif
                   !
                endif
                !
             endif
          enddo
       enddo
       !
    else
       !
       do i = 1, nverts
          !
          write (nref, 209) i
          !
          ! bottom level
          !
          write (nref, 210) dps(i)
          !
          ! water level
          !
          write (nref, 210) s1(i)
          !
          ! depth-averaged u-velocity
          !
          if ( kmax == 1 ) then
             !
             write (nref, 210) u1(i,1)
             !
          else
             !
             write (nref, 210) udep(i)
             !
          endif
          !
          ! u-velocity per layer
          !
          if ( kmax > 1 ) then
             !
             efac = maxval(abs(u1(i,:)))
             if ( efac < 1.e-8 ) then
                write (nref, 211) 'ZERO'
             else
                efac = efac * 10.**(-5)
                write (nref, 212) efac
                write (nref, '(10(1x,i7))') (nint(u1(i,k)/efac), k=1,kmax)
             endif
             !
          endif
          !
          if ( .not.oned ) then
             !
             ! depth-averaged v-velocity
             !
             if ( kmax == 1 ) then
                !
                write (nref, 210) v1(i,1)
                !
             else
                !
                write (nref, 210) vdep(i)
                !
             endif
             !
             ! v-velocity per layer
             !
             if ( kmax > 1 ) then
                !
                efac = maxval(abs(v1(i,:)))
                if ( efac < 1.e-8 ) then
                   write (nref, 211) 'ZERO'
                else
                   efac = efac * 10.**(-5)
                   write (nref, 212) efac
                   write (nref, '(10(1x,i7))') (nint(v1(i,k)/efac), k=1,kmax)
                endif
                !
             endif
             !
          endif
          !
       enddo
       !
    endif
    !
    close (nref)
    !
 201 format ('-',i3.3)
 202 format (a, t41, a)
 203 format (i6, t41, a)
 204 format ('QUANT', /, i6, t41, 'number of quantities in table')
 205 format (f14.6, t41, a)
 206 format (a18, t41, 'date and time')
 207 format (a18, t41, 'time')
 208 format (a6)
 209 format ('LOCATION', i6)
 210 format (1x,1pe14.7e2)
 211 format (a4)
 212 format ('FACTOR', /, e18.8)
    !
end subroutine SwashBackup
