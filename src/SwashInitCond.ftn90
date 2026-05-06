subroutine SwashInitCond
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
!   Sets initial state of the flow and turbulence fields
!
!   Modules used
!
    use ocpcomm1
    use ocpcomm2
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use SwashTimecomm
    use m_genarr
    use m_parall
    use SwanGriddata
    use SwanGridobjects
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer       :: i         ! loop counter
    integer, save :: ient = 0  ! number of entries in this subroutine
    integer       :: iface     ! loop counter of faces
    integer       :: iiopt     ! time coding option at hot-start file
    integer       :: ilpos     ! actual length of error message filename
    integer       :: indx      ! point index
    integer       :: iostat    ! I/O status in call FOR
    integer       :: ix        ! grid index in x-direction of (sub)domain in current run
    integer       :: iy        ! grid index in y-direction of (sub)domain in current run
    integer       :: j         ! loop counter
    integer       :: k         ! loop counter
    integer       :: mxcur     ! number of grid points in x-direction of (sub)domain in current run
    integer       :: mycur     ! number of grid points in y-direction of (sub)domain in current run
    integer       :: npcur     ! number of grid points in current run
    integer       :: nref      ! reference number
    integer       :: numlay    ! number of vertical layers at hot-start file
    integer       :: numpts    ! number of grid points at hot-start file
    integer       :: numqua    ! number of quantities at hot-start file
    !
    real          :: efac      ! reference factor for flow variables
    real          :: nx        ! x-component of normal to face
    real          :: ny        ! y-component of normal to face
    real          :: un        ! velocity normal to face
    real          :: vx        ! initial constant value for u-velocity
    real          :: vy        ! initial constant value for v-velocity
    real          :: wlev      ! initial constant value for water level
    !
    logical       :: EQCSTR    ! compares two strings
    logical       :: EQREAL    ! compares two reals
    logical       :: KEYWIS    ! indicates whether keyword in user manual is found or not
    logical       :: ptnsubgrd ! indicates whether grid point is in subdomain or not
    logical       :: singlehot ! indicates reading single or multiple hot-start file(s)
    logical       :: STPNOW    ! indicates that program must stop
    !
    character(80) :: rline     ! a line to be read
    !
    type(facetype), dimension(:), pointer :: face ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    call strace (ient,'SwashInitCond')
    !
    ! point to face object
    !
    face => gridobject%face_grid
    !
    ! ==========================================================================
    !
    !           | -> CONstant [wlev] [vx] [vy] [tke] [epsilon]
    !           |
    ! INITial  <  ZERO
    !           |
    !           | STEAdy
    !           |
    !           |           | -> MULTiple |
    !           | HOTStart <               >  'fname'
    !           |           |    SINGle   |
    !
    ! ==========================================================================
    !
    call INKEYW ('STA', 'CON')
    !
    if ( KEYWIS('ZERO') ) then
       !
       ! zero flow variables
       !
       s1 = 0.
       u1 = 0.
       if ( optg /= 5 ) then
          if ( .not.oned ) v1 = 0.
          if ( kmax > 1 ) then
             udep = 0.
             if ( .not.oned ) vdep = 0.
          endif
       endif
       !
    else if ( KEYWIS('STEA') ) then
       !
       ! (quasi-)steady flow condition for river flow simulation
       ! will be computed later by subroutine SwashInitSteady
       !
       instead = .true.
       !
    else if ( KEYWIS('HOTS') .or. KEYWIS('REST') ) then
       !
       ! hot-starting
       !
       if ( mxc <= 1 .and. optg /= 5 ) call msgerr (2, 'command INIT should follow CGRID')
       if ( mcgrd <= 1 .and. nverts <= 0 ) call msgerr (2, 'command INIT should follow read BOT or read UNSTRUC')
       !
       restrt = .true.
       !
       call INKEYW ('STA', 'MULT')
       if ( KEYWIS ('SING') ) then
          singlehot = .true.
          mxcur = MXCGL
          mycur = MYCGL
          call IGNORE ('SING')
       else
          singlehot = .false.
          mxcur = mxc - 1
          mycur = myc - 1
          call IGNORE ('MULT')
       endif
       npcur = mxcur*mycur
       if ( optg == 5 ) npcur = nverts
       call INCSTR ('FNAME', FILENM, 'REQ', ' ')
       ! append node number to FILENM in case of parallel computing
       if ( PARLL .and. .not.singlehot ) then
          ilpos = index ( FILENM, ' ' )-1
          write (FILENM(ilpos+1:ilpos+4),201) INODE
       endif
       nref   = 0
       iostat = 0
       call FOR (nref, FILENM, 'OF', iostat)
       if (STPNOW()) return
       read (nref, 202) rline
       if (rline(1:5) /= 'SWASH') call msgerr (3, FILENM//' is not a correct hotstart file')
 100   read (nref, 202) rline
       if (rline(1:1) == COMID .or. rline(1:1) == '!') goto 100
       if ( EQCSTR(rline,'TIME') ) then
          read (nref, *) iiopt
          if ( ITEST >= 50 ) write (PRTEST, 203) iiopt
          read (nref, 202) rline
          if ( nstatm == 0 ) call msgerr (1, 'time info in hotstart file ignored')
       else
          iiopt = -1
          if ( nstatm == 1 ) call msgerr (1, 'no time info in hotstart file')
       endif
       if ( EQCSTR(rline,'LOCA') .or. EQCSTR(rline,'LONLAT') ) then
          read (nref, *) numpts
          if ( numpts /= npcur ) then
             call msgerr (2, 'grid on hotstart file differs from one in CGRID command')
             write (PRINTF, 204) npcur, numpts
          endif
          if ( ITEST >= 50 ) write (PRTEST, 205) numpts
          do i = 1, numpts
             read (nref, *)
          enddo
          read (nref, 202) rline
       endif
       if ( EQCSTR(rline(1:2),'ZK') ) then
          read (nref, *) numlay
          if ( numlay /= kmax ) call msgerr (2, 'vertical grid on hotstart file differs from one in VERT command')
          if ( ITEST >= 50 ) write (PRTEST, 206) numlay
          do i = 1, numlay
             read (nref, *)
          enddo
          read (nref, 202) rline
       else
          numlay = 1
          if ( kmax > 1 ) call msgerr (1, 'no vertical grid schematisation on hotstart file')
       endif
       !
       read (nref, *) numqua
       if ( numqua < 3 .or. numqua > 7 .or. numqua ==6 ) call msgerr (2,'number of quantities on hotstart file incorrect')
       if ( oned ) then
          if ( kmax == 1 ) then
             if ( numqua /= 3 ) call msgerr (2,'number of quantities on hotstart file incorrect')
          else
             if ( numqua /= 5 ) call msgerr (2,'number of quantities on hotstart file incorrect')
          endif
       else
          if ( kmax == 1 ) then
             if ( numqua /= 4 ) call msgerr (2,'number of quantities on hotstart file incorrect')
          else
             if ( numqua /= 7 ) call msgerr (2,'number of quantities on hotstart file incorrect')
          endif
       endif
       !
       do i = 1, numqua
          read (nref, 202) rline
          if ( ITEST >= 50 ) write (PRTEST, 207) rline
          read (nref, 202) rline
          read (nref, 202) rline
       enddo
       !
       ! reading of heading is completed, read time if nonstationary
       !
       if ( .not. iiopt < 0 ) then
          read (nref, 202) rline
          call DTRETI (rline(1:18), iiopt, timco)
          write (PRINTF, 208) rline(1:18)
       endif
       !
       if ( optg /= 5 ) then
          !
          ! structured grid
          !
          do i = 1, mxcur
             !
             if ( singlehot ) then
                ix = i -MXF + 1
             else
                ix = i
             endif
             !
             do j = 1, mycur
                !
                if ( singlehot ) then
                   iy = j -MYF + 1
                else
                   iy = j
                endif
                indx = kgrpnt(ix,iy)
                !
                ptnsubgrd = .true.
                if ( singlehot .and. (MXF > i .or. MXL < i .or.MYF > j .or. MYL < j) ) ptnsubgrd = .false.
                !
                read (nref, 202) rline
                if ( indx == 1 .and. ptnsubgrd ) then
                   !
                   if ( rline(1:6) /= 'NODATA' ) then
                      call msgerr (2, 'valid flow variables for non-existing grid point')
                      write (PRINTF, *) ix, iy
                   endif
                   !
                else
                   !
                   read (nref, 202) rline
                   !
                   ! read water level
                   !
                   if ( ptnsubgrd ) then
                      read (nref, *) s1(indx)
                   else
                      read (nref, *) rline
                   endif
                   !
                   ! read depth-averaged u-velocity
                   !
                   if ( ptnsubgrd ) then
                      if ( kmax == 1 ) then
                         read (nref, *) u1(indx,1)
                      else
                         read (nref, *) udep(indx)
                      endif
                   else
                      read (nref, *) rline
                   endif
                   !
                   ! read u-velocity per layer
                   !
                   if ( kmax > 1 ) then
                      !
                      read (nref, 202) rline
                      if ( EQCSTR(rline,'NODATA') .or. EQCSTR(rline,'ZERO') ) then
                         !
                         if ( ptnsubgrd ) then
                            u1(indx,:) = 0.
                            if ( ITEST >= 150 ) write (PRTEST, 209) 'u-velocity', ix, iy
                         endif
                         !
                      else
                         !
                         ! first read factor
                         read (nref, *) efac
                         !
                         ! read u-velocity for all layers
                         if ( ptnsubgrd ) then
                            read (nref, *) (u1(indx,k), k=1,kmax)
                            u1(indx,:) = efac * u1(indx,:)
                            if ( ITEST >= 150 ) write (PRTEST, 210) 'u-velocity', ix, iy, efac
                         else
                            read (nref, *) rline
                         endif
                         !
                      endif
                      !
                   endif
                   !
                   if ( .not.oned ) then
                      !
                      ! read depth-averaged v-velocity
                      !
                      if ( ptnsubgrd ) then
                         if ( kmax == 1 ) then
                            read (nref, *) v1(indx,1)
                         else
                            read (nref, *) vdep(indx)
                         endif
                      else
                         read (nref, *) rline
                      endif
                      !
                      ! read v-velocity per layer
                      !
                      if ( kmax > 1 ) then
                         !
                         read (nref, 202) rline
                         if ( EQCSTR(rline,'NODATA') .or. EQCSTR(rline,'ZERO') ) then
                            !
                            if ( ptnsubgrd ) then
                               v1(indx,:) = 0.
                               if ( ITEST >= 150 ) write (PRTEST, 209) 'v-velocity', ix, iy
                            endif
                            !
                         else
                            !
                            ! first read factor
                            read (nref, *) efac
                            !
                            ! read v-velocity for all layers
                            if ( ptnsubgrd ) then
                               read (nref, *) (v1(indx,k), k=1,kmax)
                               v1(indx,:) = efac * v1(indx,:)
                               if ( ITEST >= 150 ) write (PRTEST, 210) 'v-velocity', ix, iy, efac
                            else
                               read (nref, *) rline
                            endif
                            !
                         endif
                         !
                      endif
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          enddo
          !
       else
          !
          ! unstructured grids
          !
          do i = 1, nverts
             !
             read (nref, 202) rline
             read (nref, 202) rline
             !
             ! read water level
             !
             read (nref, *) s1(i)
             !
             ! read depth-averaged u-velocity
             !
             if ( kmax == 1 ) then
                read (nref, *) u1(i,1)
             else
                read (nref, *) udep(i)
             endif
             !
             ! read u-velocity per layer
             !
             if ( kmax > 1 ) then
                !
                read (nref, 202) rline
                if ( EQCSTR(rline,'NODATA') .or. EQCSTR(rline,'ZERO') ) then
                   !
                   u1(i,:) = 0.
                   if ( ITEST >= 150 ) write (PRTEST, 211) 'u-velocity', i
                   !
                else
                   !
                   ! first read scaling factor
                   read (nref, *) efac
                   !
                   ! read u-velocity for all layers
                   read (nref, *) (u1(i,k), k=1,kmax)
                   u1(i,:) = efac * u1(i,:)
                   if ( ITEST >= 150 ) write (PRTEST, 212) 'u-velocity', i, efac
                   !
                endif
                !
             endif
             !
             if ( .not.oned ) then
                !
                ! read depth-averaged v-velocity
                !
                if ( kmax == 1 ) then
                   read (nref, *) v1(i,1)
                else
                   read (nref, *) vdep(i)
                endif
                !
                ! read v-velocity per layer
                !
                if ( kmax > 1 ) then
                   !
                   read (nref, 202) rline
                   if ( EQCSTR(rline,'NODATA') .or. EQCSTR(rline,'ZERO') ) then
                      !
                      v1(i,:) = 0.
                      if ( ITEST >= 150 ) write (PRTEST, 211) 'v-velocity', i
                      !
                   else
                      !
                      ! first read scaling factor
                      read (nref, *) efac
                      !
                      ! read v-velocity for all layers
                      read (nref, *) (v1(i,k), k=1,kmax)
                      v1(i,:) = efac * v1(i,:)
                      if ( ITEST >= 150 ) write (PRTEST, 212) 'v-velocity', i, efac
                      !
                   endif
                   !
                endif
                !
             endif
             !
          enddo
          !
       endif
       close (nref)
       !
    else
       !
       ! default initial variables
       !
       call IGNORE ('CON')
       call INREAL('WLEV',    wlev,   'STA',0.)
       call INREAL('VX',      vx,     'STA',0.)
       call INREAL('VY',      vy,     'STA',0.)
       call INREAL('TKE',     tkeini, 'UNC',0.)
       call INREAL('EPSILON', epsini, 'UNC',0.)
       !
       s1 = wlev
       if ( optg /= 5 ) then
          u1 = vx
          if ( .not.oned ) v1 = vy
          if ( kmax > 1 ) then
             udep = vx
             if ( .not.oned ) vdep = vy
          endif
          !
          ! set to zero for permanently dry points
          !
          s1(1  ) = 0.
          u1(1,:) = 0.
          if ( .not.oned ) v1(1,:) = 0.
          if ( kmax > 1 ) then
             udep(1) = 0.
             if ( .not.oned ) vdep(1) = 0.
          endif
       else
          do iface = 1, nfaces
             nx = face(iface)%attr(FACENORMX)
             ny = face(iface)%attr(FACENORMY)
             un = nx * vx + ny * vy
             u1(iface,:) = un
             if ( kmax > 1 ) udep(iface) = un
          enddo
       endif
       !
    endif
    !
 201 format ('-',i3.3)
 202 format (a)
 203 format (' time coding option:', i2)
 204 format (1x, i6, ' points in computational grid, on hotstart file:', i6)
 205 format (1x, i6, '  output locations')
 206 format (1x, i6, '  vertical layers')
 207 format (1x, 'quantity: ', a)
 208 format (' initial condition read for time: ', a)
 209 format (' zero ', a10, ' or no data for point:', 2i4)
 210 format (1x, a10, ' in point:', 2i4,'  factor=', e12.4)
 211 format (' zero ', a10, 'or no data for vertex:', i6)
 212 format (1x, a10, ' in vertex:', i6, '  factor=', e12.4)
    !
end subroutine SwashInitCond
