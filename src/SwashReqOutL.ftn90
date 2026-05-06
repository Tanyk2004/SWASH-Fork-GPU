subroutine SwashReqOutL ( found )
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
!   Reading and processing of specifying output locations
!
!   Modules used
!
    use ocpcomm1
    use ocpcomm4
    use SwashCommdata2, only: igtype
    use SwashCommdata3
    use SwashCommdata4
    use outp_data
    use m_genarr, only: depth, wlevl
    use m_parall
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    logical, intent(inout)             :: found          ! keyword found
!
!   Local variables
!
    integer                            :: i              ! loop counter
    integer, save                      :: ient = 0       ! number of entries in this subroutine
    integer                            :: ierr           ! error indicator: ierr=0: no error, otherwise error
    integer                            :: ii             ! indicator of function SIRAY
                                                         ! =0 no point found
                                                         ! =1 point found
    integer                            :: intd           ! interval containing number of rays
    integer                            :: intv           ! interval containing number of given locations on curve
    integer                            :: iostat         ! I/O status in call FOR
    integer                            :: ix1            ! lowest grid index of (sub)grid in x-direction
    integer                            :: ix2            ! highest grid index of (sub)grid in x-direction
    integer                            :: iy1            ! lowest grid index of (sub)grid in y-direction
    integer                            :: iy2            ! highest grid index of (sub)grid in y-direction
    integer                            :: mip            ! number of output points
    integer                            :: mipr           ! number of ray points
    integer                            :: mxfr           ! number of cells in x-direction in frame
    integer                            :: myfr           ! number of cells in y-direction in frame
    integer                            :: nds            ! file reference number
    integer                            :: SIRAY          ! searching for point on a ray
    !
    real                               :: alpfr          ! direction of x-axis of frame w.r.t. user coordinates
    real                               :: dp             ! local depth
    real                               :: xlenfr         ! length of frame in x-direction
    real                               :: xp             ! x-coordinate of subsequent points / start point of ray
    real                               :: xp1            ! x-coordinate of first point
    real                               :: xpfr           ! x-coordinate of origin of frame
    real                               :: xq             ! x-coordinate end point of subsequent rays
    real                               :: xq1            ! x-coordinate of end point of first ray
    real                               :: xx             ! x-coordinate point of ray with depth dp
    real                               :: ylenfr         ! length of frame in y-direction
    real                               :: yp             ! y-coordinate of subsequent points / start point of ray
    real                               :: yp1            ! y-coordinate of first point
    real                               :: ypfr           ! y-coordinate of origin of frame
    real                               :: yq             ! y-coordinate end point of subsequent rays
    real                               :: yq1            ! y-coordinate of end point of first ray
    real                               :: yy             ! y-coordinate point of ray with depth dp
    !
    character(16)                      :: prname         ! name of ray
    character(16)                      :: psname         ! name of point set
    !
    logical                            :: botdep         ! indicates bottom depth or water depth
    logical                            :: EQREAL         ! compares two reals
    logical                            :: KEYWIS         ! indicates whether keyword in user manual is found or not
    logical                            :: lfile          ! indicates reading from file or not
    logical                            :: STPNOW         ! indicates that program must stop
    !
    type(opsdat), pointer              :: opstmp         ! list containing parameters for grid
    type(opsdat), pointer              :: rops           ! list containing parameters for ray
    !
    type xypt                                            ! linked list for points of curve or ray (see command CURV / RAY)
       real                :: x, y, xq, yq
       type(xypt), pointer :: nextxy
    end type xypt
    type(xypt), target     :: frst
    type(xypt), pointer    :: curr, tmp
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashReqOutL')
    !
    ! ==========================================================================
    !
    ! FRAme   'sname'  [xpfr] [ypfr] [alpfr] [xlenfr] [ylenfr] [mxfr] [myfr]
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('FRA') ) then
       !
       allocate(opstmp)
       call INCSTR ('SNAME',psname,'REQ',' ')
       if ( LENCST > 8 ) call msgerr (2, 'sname is too long')
       opstmp%psname = psname
       call READXY ('XPFR', 'YPFR', xpfr, ypfr, 'REQ', 0., 0.)
       opstmp%opr(1) = xpfr
       opstmp%opr(2) = ypfr
       call INREAL('ALPFR',alpfr,'REQ',0.)
       if ( kspher > 0 .and. .not.EQREAL(alpfr,0.) ) call msgerr (2, '[alpfr] must be 0 in case of spherical coordinates')
       call INREAL('XLENFR', xlenfr,'REQ',0.)
       call INREAL('YLENFR', ylenfr,'REQ',0.)
       opstmp%opr(3) = xlenfr
       opstmp%opr(4) = ylenfr
       opstmp%opr(5) = pi2 * ( alpfr/360.-nint(alpfr/360.) )
       opstmp%opr(6) = 0.
       !
       ! the user gives number of cells along each side
       !
       call ININTG ('MXFR',mxfr,'STA',20)
       call ININTG ('MYFR',myfr,'STA',20)
       !
       opstmp%pstype = 'F'
       opstmp%opi(1) = mxfr + 1
       opstmp%opi(2) = myfr + 1
       allocate(opstmp%xp(0))
       allocate(opstmp%yp(0))
       nullify(opstmp%nextops)
       if ( .not.lops ) then
          fops = opstmp
          cops => fops
          lops = .true.
       else
          cops%nextops => opstmp
          cops => opstmp
       endif
       goto 800
       !
    endif
    !
    ! ==========================================================================
    !
    ! GROUP   'sname'  SUBGRID [ix1] [ix2] [iy1] [iy2]
    !
    ! ==========================================================================
    !
    if (KEYWIS('GROUP') .or. KEYWIS ('SUBG')) then
       !
       if ( optg == 5 ) then
          call msgerr (2, 'keyword GROUP not supported in unstructured grid')
          goto 800
       else
          allocate(opstmp)
          call INCSTR ('SNAME',psname,'REQ',' ')
          if ( LENCST > 8 ) call msgerr (2, 'sname is too long')
          opstmp%psname = psname
          call INKEYW ('STA', ' ')
          call IGNORE ('SUBG')
          call ININTG ('IX1', ix1, 'REQ', 1)
          call ININTG ('IX2', ix2, 'REQ', 1)
          call ININTG ('IY1', iy1, 'REQ', 1)
          call ININTG ('IY2', iy2, 'REQ', 1)
          if ( ix1 < 1 .or. ix2 > MXCGL-1 .or. ix1 > ix2 .or. iy1 < 1 .or. iy2 > MYCGL-1 .or. iy1 > iy2 ) then
             call msgerr (3, 'Check corners of GROUP command')
             call msgerr (3, ' .........the values should be.........')
             call msgerr (3, 'ix1 <= ix2 and ix1, ix2 in [1,MXC+1]')
             call msgerr (3, 'iy1 <= iy2 and iy1, iy2 in [1,MYC+1]')
          endif
          !
          ix1 = ix1 - 1
          ix2 = ix2 - 1
          iy1 = iy1 - 1
          iy2 = iy2 - 1
          !
          if ( optg == 3 ) then
             opstmp%pstype = 'H'
             opstmp%opr(1) = float(ix2)
             opstmp%opr(2) = float(iy2)
             opstmp%opr(3) = float(ix1)
             opstmp%opr(4) = float(iy1)
             opstmp%opr(5) = 0.
             opstmp%opr(6) = 1.
             opstmp%opi(1) = ix2 - ix1 + 1
             opstmp%opi(2) = iy2 - iy1 + 1
          else if ( optg == 1 ) then
             opstmp%pstype = 'F'
             if ( ix1 /= ix2 ) then
               opstmp%opr(3) = ( ix2 - ix1 ) * dx
             else
               opstmp%opr(3) = 0.01
             endif
             if ( iy1 /= iy2 ) then
               opstmp%opr(4) = ( iy2 - iy1 ) * dy
             else
               opstmp%opr(4) = 0.01
             endif
             opstmp%opr(1) = xpc + ix1 *dx * cospc - iy1 * dy * sinpc
             opstmp%opr(2) = ypc + ix1 *dx * sinpc + iy1 * dy * cospc
             opstmp%opr(5) = alpc
             opstmp%opr(6) = 1.
             opstmp%opi(1) = ix2 - ix1 + 1
             opstmp%opi(2) = iy2 - iy1 + 1
             if ( ITEST >= 20 .or. intes >= 10 ) write (PRINTF, 201) (opstmp%opr(i), i=1,6)
          endif
          allocate(opstmp%xp(0))
          allocate(opstmp%yp(0))
          nullify(opstmp%nextops)
          if ( .not.lops ) then
             fops = opstmp
             cops => fops
             lops = .true.
          else
             cops%nextops => opstmp
             cops => opstmp
          endif
          goto 800
       endif
       !
    endif
    !
    ! ==========================================================================
    !
    ! CURVe   'sname'  [xp1] [yp1]   < [int]  [xp]  [yp] >
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('CURV') ) then
       !
       allocate(opstmp)
       call INCSTR ('SNAME',psname,'REQ',' ')
       if ( LENCST > 8 ) call msgerr (2, 'sname is too long')
       opstmp%psname = psname
       opstmp%pstype = 'C'
       opstmp%mip    = 0
       !
       ! first point of a curve
       !
       call NWLINE
       if (STPNOW()) return
       call READXY ('XP1', 'YP1', xp, yp, 'REQ', 0., 0.)
       frst%x = xp
       frst%y = yp
       nullify(frst%nextxy)
       curr => frst
       mip = 1
       !
       ! interval and next corner point
       !
 100   call ININTG ('INT',intv,'REP',-1)
       if ( intv /= -1 ) then
          if ( intv <= 0 ) then
             call msgerr (2, '[int] is negative or zero')
             intv = 1
          endif
          xp1 = xp
          yp1 = yp
          call READXY ('XP', 'YP', xp, yp, 'REQ', 0., 0.)
          if ( ITEST >= 200 .or. intes >= 20 ) write (PRINTF, 202) psname
          do i = 1, intv
             mip = mip + 1
             allocate(tmp)
             tmp%x = xp1 + real(i) * ( xp - xp1 )/ real(intv)
             tmp%y = yp1 + real(i) * ( yp - yp1 )/ real(intv)
             nullify(tmp%nextxy)
             curr%nextxy => tmp
             curr => tmp
          enddo
          goto 100
       endif
       allocate(opstmp%xp(mip))
       allocate(opstmp%yp(mip))
       curr => frst
       do i = 1, mip
          opstmp%xp(i) = curr%x
          opstmp%yp(i) = curr%y
          if ( ITEST >= 200 .or. intes >= 50 ) then
             write (PRINTF,203) i, curr%x, curr%y
          endif
          curr => curr%nextxy
       enddo
       deallocate(tmp)
       !
       ! store number of points of the curve
       !
       opstmp%mip = mip
       if ( mip == 0 ) call msgerr (1, 'no output points on curve found')
       nullify(opstmp%nextops)
       if ( .not.lops ) then
          fops = opstmp
          cops => fops
          lops = .true.
       else
          cops%nextops => opstmp
          cops => opstmp
       endif
       goto 800
       !
    endif
    !
    ! ==========================================================================
    !
    ! POINts  'sname'  < [xp]  [yp]  >   |  FILE 'fname'
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('POIN') ) then
       !
       allocate(opstmp)
       call INCSTR ('SNAME',psname,'REQ',' ')
       if ( LENCST > 8 ) call msgerr (2, 'sname is too long')
       opstmp%psname = psname
       opstmp%pstype = 'P'
       mip = 0
       opstmp%mip = mip
       call INKEYW ('STA', ' ')
       if ( KEYWIS('FILE') ) then
          iostat = 0
          nds    = 0
          lfile  = .true.
          call INCSTR ('FNAME', FILENM, 'REQ', ' ')
          call FOR (nds, FILENM, 'OF', iostat)
          if (STPNOW()) return
       else
          lfile = .false.
       endif
       frst%x = 0.
       frst%y = 0.
       nullify(frst%nextxy)
       curr => frst
       do
         if ( lfile ) then
           ierr = 0
           call REFIXY (nds, xp, yp, ierr)
           if ( ierr == -1 ) goto 110
           if ( ierr == -2 ) then
              call msgerr (2, 'error reading point coordinates from file')
              goto 800
           endif
         else
           call READXY ('XP', 'YP', xp, yp, 'REP', -1.e10, -1.e10)
           if ( xp < -0.9e10 ) goto 110
         endif
         mip = mip + 1
         allocate(tmp)
         tmp%x = xp
         tmp%y = yp
         nullify(tmp%nextxy)
         curr%nextxy => tmp
         curr => tmp
       enddo
 110   allocate(opstmp%xp(mip))
       allocate(opstmp%yp(mip))
       curr => frst%nextxy
       do i = 1, mip
          opstmp%xp(i) = curr%x
          opstmp%yp(i) = curr%y
          curr => curr%nextxy
       enddo
       deallocate(tmp)
       !
       ! store number of output points
       !
       opstmp%mip = mip
       if ( mip == 0 ) call msgerr (2, 'no output points found')
       nullify(opstmp%nextops)
       if ( .not.lops ) then
          fops = opstmp
          cops => fops
          lops = .true.
       else
          cops%nextops => opstmp
          cops => opstmp
       endif
       goto 800
       !
    endif
    !
    ! ==========================================================================
    !
    ! RAY     'rname'  [xp1] [yp1] [xq1] [yq1]                       &
    !        <  [int]  [xp]  [yp]  [xq]  [yq]  >
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('RAY') ) then
       !
       if ( oned ) then
          call msgerr (2,'illegal keyword (RAY) in combination with 1D computation')
          goto 800
       else if ( igtype(1) == 2 .or. igtype(7) == 2 ) then
          call msgerr (2, 'keyword RAY not supported in case of curvilinear bottom grid')
          goto 800
       else
          allocate(opstmp)
          call INCSTR ('RNAME',psname,'REQ',' ')
          if ( LENCST > 8 ) call msgerr (2, 'rname is too long')
          opstmp%psname = psname
          opstmp%pstype = 'R'
          mip = 1
          opstmp%mip = mip
          !
          ! first ray
          !
          call NWLINE
          if (STPNOW()) return
          call READXY ('XP1', 'YP1', xp, yp, 'REQ', 0., 0.)
          call READXY ('XQ1', 'YQ1', xq, yq, 'REQ', 0., 0.)
          frst%x  = xp
          frst%y  = yp
          frst%xq = xq
          frst%yq = yq
          nullify(frst%nextxy)
          curr => frst
          !
          ! following rays
          !
  120     call ININTG ('INT',intd,'REP',-1)
          if ( intd /= -1 ) then
             if ( intd <= 0 ) then
                call msgerr (2, '[int] negative or zero')
                intd = 1
             endif
             xp1 = xp
             yp1 = yp
             xq1 = xq
             yq1 = yq
             call READXY ('XP', 'YP', xp, yp, 'REQ', 0., 0.)
             call READXY ('XQ', 'YQ', xq, yq, 'REQ', 0., 0.)
             do i = 1, intd
                mip = mip + 1
                allocate(tmp)
                tmp%x  = xp1 + real(i) * ( xp - xp1 ) / real(intd)
                tmp%y  = yp1 + real(i) * ( yp - yp1 ) / real(intd)
                tmp%xq = xq1 + real(i) * ( xq - xq1 ) / real(intd)
                tmp%yq = yq1 + real(i) * ( yq - yq1 ) / real(intd)
                nullify(tmp%nextxy)
                curr%nextxy => tmp
                curr => tmp
             enddo
             goto 120
          endif
          allocate(opstmp%xp(mip))
          allocate(opstmp%yp(mip))
          allocate(opstmp%xq(mip))
          allocate(opstmp%yq(mip))
          curr => frst
          do i = 1, mip
             opstmp%xp(i) = curr%x
             opstmp%yp(i) = curr%y
             opstmp%xq(i) = curr%xq
             opstmp%yq(i) = curr%yq
             curr => curr%nextxy
          enddo
          deallocate(tmp)
          !
          opstmp%mip = mip
          if ( mip == 1 ) call msgerr (1,'only one ray is defined')
          nullify(opstmp%nextops)
          if ( .not.lops ) then
             fops = opstmp
             cops => fops
             lops = .true.
          else
             cops%nextops => opstmp
             cops => opstmp
          endif
          goto 800
       endif
       !
    endif
    !
    ! ==========================================================================
    !
    ! ISOline 'sname'  'rname'  DEPTH / BOTTOM [dep]
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('ISO') ) then
       !
       if ( oned ) then
          call msgerr (2, 'illegal keyword (ISO) in combination with 1D computation')
          goto 800
       else if ( igtype(1) == 2 .or. igtype(7) == 2 ) then
          call msgerr (2, 'keyword ISO not supported in case of curvilinear bottom grid')
          goto 800
       else
          allocate(opstmp)
          call INCSTR ('SNAME',psname,'REQ',' ')
          if ( LENCST > 8 ) call msgerr (2, 'sname is too long')
          opstmp%psname = psname
          call INCSTR ('RNAME', prname, 'REQ', ' ')
          if ( LENCST > 8 ) call msgerr (2, 'rname is too long')
          call INKEYW ('STA', 'DEP')
          if ( KEYWIS ('BOT') ) then
             botdep = .true.
          else
             call IGNORE ('DEP')
             botdep = .false.
          endif
          call INREAL ('DEP',dp,'REQ',0.)
          rops => fops
          do
            if ( rops%psname == prname ) exit
            if (.not.associated(rops%nextops) ) then
               call msgerr (2, 'set of rays not defined')
               goto 800
            endif
            rops => rops%nextops
          enddo
          if ( rops%pstype /= 'R' ) then
             call msgerr (2, 'ray name set assigned to set of output locations')
             goto 800
          endif
          mipr = rops%mip
          opstmp%pstype = 'C'
          mip = 0
          opstmp%mip = mip
          frst%x = 0.
          frst%y = 0.
          nullify(frst%nextxy)
          curr => frst
          do i = 1, mipr
             xp = rops%xp(i)
             yp = rops%yp(i)
             xq = rops%xq(i)
             yq = rops%yq(i)
             ii = SIRAY (dp, xp, yp, xq, yq, xx, yy, botdep, depth, wlevl)
             if ( ii == 0 ) then
                write (PRINTF, 204) dp, xp+xoffs, yp+yoffs, xq+xoffs, yq+yoffs
             else
                mip = mip + 1
                allocate(tmp)
                tmp%x = xx
                tmp%y = yy
                nullify(tmp%nextxy)
                curr%nextxy => tmp
                curr => tmp
             endif
          enddo
          allocate(opstmp%xp(mip))
          allocate(opstmp%yp(mip))
          curr => frst%nextxy
          do i = 1, mip
             opstmp%xp(i) = curr%x
             opstmp%yp(i) = curr%y
             curr => curr%nextxy
          enddo
          deallocate(tmp)
          if ( mip == 0 ) call msgerr (2, 'no points with valid depth found')
          !
          ! store number of points of the curve
          !
          opstmp%mip = mip
          nullify(opstmp%nextops)
          if ( .not.lops ) then
             fops = opstmp
             cops => fops
             lops = .true.
          else
             cops%nextops => opstmp
             cops => opstmp
          endif
          goto 800
       endif
       !
    endif
    !
    ! command not found
    !
    return
    !
 800 found = .true.
    return
    !
 201 format (' Subgrid parms.', 6(1x,e12.4))
 202 format (' COORDINATES OF OUTPUT POINTS FOR CURVE  : ', a)
 203 format (' POINT(',i4,')','  (IX,IY) -> ',2f10.2)
 204 format (' No point with depth ',f5.2,' is found in ray :',4f10.2)
 205 format ((6x,2e22.15,i3))
    !
end subroutine SwashReqOutL
