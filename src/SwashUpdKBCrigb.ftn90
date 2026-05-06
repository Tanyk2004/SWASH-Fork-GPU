subroutine SwashUpdKBCrigb
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
!    9.01: Dirk Rijnsdorp, Marcel Zijlema
!
!   Updates
!
!    9.01, October 2022: New subroutine
!
!   Purpose
!
!   Updates kinematic boundary conditions at rigid body surface
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr, only: kgrpnt, guu, gvv, xcgrid, ycgrid
    use SwashFlowdata
    use SwashRigBoddata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: l        ! label of body
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mu       ! index of point m+1
    integer       :: nd       ! index of point n-1
    integer       :: nu       ! index of point n+1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: num      ! pointer to m,n+1
    !
    real          :: dsdx     ! slope of hull surface in x-direction
    real          :: dsdy     ! slope of hull surface in y-direction
    real          :: dxl      ! local mesh size in x-direction
    real          :: dyl      ! local mesh size in y-direction
    real          :: rx       ! x-component of the position vector on wetted hull surface relative to COG
    real          :: ry       ! y-component of the position vector on wetted hull surface relative to COG
    real          :: rz       ! z-component of the position vector on wetted hull surface relative to COG
    real          :: theta    ! implicitness factor
    real          :: vrb      ! velocity component of rigid body (either linear or angular)
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdKBCrigb')
    !
    ! if not moving rigid bodies, return
    !
    if ( ifloat /= 2 ) return
    !
    theta = pship(11)
    !
    skc = 0.
    !
    if ( oned ) then
       !
       ! loop over wl-points in 1D computational grid
       !
       do m = mfu, ml
          !
          md = m - 1
          mu = m + 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          nmu = kgrpnt(mu,1)
          !
          ! surge - contribution underneath body
          !
          if ( presp(nm) == 1 ) then
             !
             l = lfbs(nm)
             if ( l == 0 .or. .not.bdof(l,1,1) ) goto 10
             !
             ! surge velocity
             !
             vrb = theta * vfot1(l,1) + (1.-theta) * vfot0(l,1)
             !
             ! slope of hull
             ! note: flou is positive downwards
             !
             dsdx = ( flou(nm) - flou(nmd) ) / dx
             !
             skc(nm) = skc(nm) + vrb * dsdx
             !
          endif
          !
          ! surge - contribution next to body (from free surface to pressurized in x-direction)
          !
  10      if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
             !
             l = lfbs(nmu)
             if ( l == 0 .or. .not.bdof(l,1,1) ) goto 20
             !
             ! surge velocity
             !
             vrb = theta * vfot1(l,1) + (1.-theta) * vfot0(l,1)
             !
             ! slope of hull on the port side
             ! note: flos is positive downwards
             !
             dsdx = -( flos(nmu) + s1(nm) ) / dx
             !
             skc(nm) = skc(nm) - vrb * dsdx
             !
          endif
          !
          ! surge - contribution next to body (from pressurized to free surface in x-direction)
          !
  20      if ( presp(nm) == 0 .and. presp(nmd) == 1 ) then
             !
             l = lfbs(nmd)
             if ( l == 0 .or. .not.bdof(l,1,1) ) goto 30
             !
             ! surge velocity
             !
             vrb = theta * vfot1(l,1) + (1.-theta) * vfot0(l,1)
             !
             ! slope of hull on the starboard side
             ! note: flos is positive downwards
             !
             dsdx = ( s1(nm) + flos(nmd) ) / dx
             !
             skc(nm) = skc(nm) - vrb * dsdx
             !
          endif
          !
          ! heave - contribution underneath body
          !
  30      if ( presp(nm) == 1 ) then
             !
             l = lfbs(nm)
             if ( l == 0 .or. .not.bdof(l,3,1) ) goto 40
             !
             ! heave velocity
             !
             vrb = theta * vfot1(l,3) + (1.-theta) * vfot0(l,3)
             !
             skc(nm) = skc(nm) + vrb
             !
          endif
          !
          ! pitch - contribution underneath body
          !
  40      if ( presp(nm) == 1 ) then
             !
             l = lfbs(nm)
             if ( l == 0 .or. .not.bdof(l,2,2) ) goto 50
             !
             ! pitch velocity
             !
             vrb = theta * vfor1(l,2) + (1.-theta) * vfor0(l,2)
             !
             ! position of wetted hull surface
             ! note: flou is positive downwards
             !
             rx = 0.5*( xcgrid(m,1) + xcgrid(md,1) ) - bcog(l,1)
             !
             rz = -0.5*( flou(nm) + flou(nmd) ) - bcog(l,3)
             !
             ! slope of hull
             !
             dsdx = ( flou(nm) - flou(nmd) ) / dx
             !
             skc(nm) = skc(nm) + vrb * rz * dsdx - vrb * rx
             !
          endif
          !
          ! pitch - contribution next to body (from free surface to pressurized in x-direction)
          !
  50      if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
             !
             l = lfbs(nmu)
             if ( l == 0 .or. .not.bdof(l,2,2) ) goto 60
             !
             ! pitch velocity
             !
             vrb = theta * vfor1(l,2) + (1.-theta) * vfor0(l,2)
             !
             ! position of wetted hull surface
             ! note: flos is positive downwards
             !
             rz = 0.5*( s1(nm) - flos(nmu) ) - bcog(l,3)
             !
             ! slope of hull on the port side
             !
             dsdx = -( flos(nmu) + s1(nm) ) / dx
             !
             skc(nm) = skc(nm) - vrb * rz * dsdx
             !
          endif
          !
          ! pitch - contribution next to body (from pressurized to free surface in x-direction)
          !
  60      if ( presp(nm) == 0 .and. presp(nmd) == 1 ) then
             !
             l = lfbs(nmd)
             if ( l == 0 .or. .not.bdof(l,2,2) ) goto 70
             !
             ! pitch velocity
             !
             vrb = theta * vfor1(l,2) + (1.-theta) * vfor0(l,2)
             !
             ! position of wetted hull surface
             ! note: flos is positive downwards
             !
             rz = 0.5*( s1(nm) - flos(nmd) ) - bcog(l,3)
             !
             ! slope of hull on the starboard side
             !
             dsdx = ( s1(nm) + flos(nmd) ) / dx
             !
             skc(nm) = skc(nm) - vrb * rz * dsdx
             !
          endif
          !
  70      continue
          !
       enddo
       !
    else
       !
       ! loop over wl-points in 2D computational grid
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             nmu = kgrpnt(mu,n )
             ndm = kgrpnt(m ,nd)
             num = kgrpnt(m ,nu)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmd == 1 ) nmd = nm
             if ( nmu == 1 ) nmu = nm
             if ( ndm == 1 ) ndm = nm
             if ( num == 1 ) num = nm
             !
             ! surge - contribution underneath body
             !
             if ( presp(nm) == 1 ) then
                !
                l = lfbs(nm)
                if ( l == 0 .or. .not.bdof(l,1,1) ) goto 80
                !
                ! surge velocity
                !
                vrb = theta * vfot1(l,1) + (1.-theta) * vfot0(l,1)
                !
                dxl = 0.5*( gvv(nm) + gvv(ndm) )
                !
                ! slope of hull
                ! note: flou is positive downwards
                !
                dsdx = ( flou(nm) - flou(nmd) ) / dxl
                !
                skc(nm) = skc(nm) + vrb * dsdx
                !
             endif
             !
             ! surge - contribution next to body (from free surface to pressurized in x-direction)
             !
  80         if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
                !
                l = lfbs(nmu)
                if ( l == 0 .or. .not.bdof(l,1,1) ) goto 90
                !
                ! surge velocity
                !
                vrb = theta * vfot1(l,1) + (1.-theta) * vfot0(l,1)
                !
                dxl = 0.5*( gvv(nm) + gvv(ndm) )
                !
                ! slope of hull on the port side
                ! note: flos is positive downwards
                !
                dsdx = -( flos(nmu) + s1(nm) ) / dxl
                !
                skc(nm) = skc(nm) - vrb * dsdx
                !
             endif
             !
             ! surge - contribution next to body (from pressurized to free surface in x-direction)
             !
  90         if ( presp(nm) == 0 .and. presp(nmd) == 1 ) then
                !
                l = lfbs(nmd)
                if ( l == 0 .or. .not.bdof(l,1,1) ) goto 100
                !
                ! surge velocity
                !
                vrb = theta * vfot1(l,1) + (1.-theta) * vfot0(l,1)
                !
                dxl = 0.5*( gvv(nm) + gvv(ndm) )
                !
                ! slope of hull on the starboard side
                ! note: flos is positive downwards
                !
                dsdx = ( s1(nm) + flos(nmd) ) / dxl
                !
                skc(nm) = skc(nm) - vrb * dsdx
                !
             endif
             !
             ! sway - contribution underneath body
             !
 100         if ( presp(nm) == 1 ) then
                !
                l = lfbs(nm)
                if ( l == 0 .or. .not.bdof(l,2,1) ) goto 110
                !
                ! sway velocity
                !
                vrb = theta * vfot1(l,2) + (1.-theta) * vfot0(l,2)
                !
                dyl = 0.5*( guu(nm) + guu(nmd) )
                !
                ! slope of hull
                ! note: flov is positive downwards
                !
                dsdy = ( flov(nm) - flov(ndm) ) / dyl
                !
                skc(nm) = skc(nm) + vrb * dsdy
                !
             endif
             !
             ! sway - contribution next to body (from free surface to pressurized in y-direction)
             !
 110         if ( presp(nm) == 0 .and. presp(num) == 1 ) then
                !
                l = lfbs(num)
                if ( l == 0 .or. .not.bdof(l,2,1) ) goto 120
                !
                ! sway velocity
                !
                vrb = theta * vfot1(l,2) + (1.-theta) * vfot0(l,2)
                !
                dyl = 0.5*( guu(nm) + guu(nmd) )
                !
                ! slope of hull on the port side
                ! note: flos is positive downwards
                !
                dsdy = -( flos(num) + s1(nm) ) / dyl
                !
                skc(nm) = skc(nm) - vrb * dsdy
                !
             endif
             !
             ! sway - contribution next to body (from pressurized to free surface in y-direction)
             !
 120         if ( presp(nm) == 0 .and. presp(ndm) == 1 ) then
                !
                l = lfbs(ndm)
                if ( l == 0 .or. .not.bdof(l,2,1) ) goto 130
                !
                ! sway velocity
                !
                vrb = theta * vfot1(l,2) + (1.-theta) * vfot0(l,2)
                !
                dyl = 0.5*( guu(nm) + guu(nmd) )
                !
                ! slope of hull on the starboard side
                ! note: flos is positive downwards
                !
                dsdy = ( s1(nm) + flos(ndm) ) / dyl
                !
                skc(nm) = skc(nm) - vrb * dsdy
                !
             endif
             !
             ! heave - contribution underneath body
             !
 130         if ( presp(nm) == 1 ) then
                !
                l = lfbs(nm)
                if ( l == 0 .or. .not.bdof(l,3,1) ) goto 140
                !
                ! heave velocity
                !
                vrb = theta * vfot1(l,3) + (1.-theta) * vfot0(l,3)
                !
                skc(nm) = skc(nm) + vrb
                !
             endif
             !
             ! roll - contribution underneath body
             !
 140         if ( presp(nm) == 1 ) then
                !
                l = lfbs(nm)
                if ( l == 0 .or. .not.bdof(l,1,2) ) goto 150
                !
                ! roll velocity
                !
                vrb = theta * vfor1(l,1) + (1.-theta) * vfor0(l,1)
                !
                ! position of wetted hull surface
                ! note: flov is positive downwards
                !
                ry = 0.25*( ycgrid(m,n) + ycgrid(m,nd) + ycgrid(md,n) + ycgrid(md,nd) ) - bcog(l,2)
                !
                rz = -0.5*( flov(nm) + flov(ndm) ) - bcog(l,3)
                !
                dyl = 0.5*( guu(nm) + guu(nmd) )
                !
                ! slope of hull
                !
                dsdy = ( flov(nm) - flov(ndm) ) / dyl
                !
                skc(nm) = skc(nm) - vrb * rz * dsdy + vrb * ry
                !
             endif
             !
             ! roll - contribution next to body (from free surface to pressurized in y-direction)
             !
 150         if ( presp(nm) == 0 .and. presp(num) == 1 ) then
                !
                l = lfbs(num)
                if ( l == 0 .or. .not.bdof(l,1,2) ) goto 160
                !
                ! roll velocity
                !
                vrb = theta * vfor1(l,1) + (1.-theta) * vfor0(l,1)
                !
                ! position of wetted hull surface
                ! note: flos is positive downwards
                !
                rz = 0.5*( s1(nm) - flos(num) ) - bcog(l,3)
                !
                dyl = 0.5*( guu(nm) + guu(nmd) )
                !
                ! slope of hull on the port side
                !
                dsdy = -( flos(num) + s1(nm) ) / dyl
                !
                skc(nm) = skc(nm) + vrb * rz * dsdy
                !
             endif
             !
             ! roll - contribution next to body (from pressurized to free surface in y-direction)
             !
 160         if ( presp(nm) == 0 .and. presp(ndm) == 1 ) then
                !
                l = lfbs(ndm)
                if ( l == 0 .or. .not.bdof(l,1,2) ) goto 170
                !
                ! roll velocity
                !
                vrb = theta * vfor1(l,1) + (1.-theta) * vfor0(l,1)
                !
                ! position of wetted hull surface
                ! note: flos is positive downwards
                !
                rz = 0.5*( s1(nm) - flos(ndm) ) - bcog(l,3)
                !
                dyl = 0.5*( guu(nm) + guu(nmd) )
                !
                ! slope of hull on the starboard side
                !
                dsdy = ( s1(nm) + flos(ndm) ) / dyl
                !
                skc(nm) = skc(nm) + vrb * rz * dsdy
                !
             endif
             !
             ! pitch - contribution underneath body
             !
 170         if ( presp(nm) == 1 ) then
                !
                l = lfbs(nm)
                if ( l == 0 .or. .not.bdof(l,2,2) ) goto 180
                !
                ! pitch velocity
                !
                vrb = theta * vfor1(l,2) + (1.-theta) * vfor0(l,2)
                !
                ! position of wetted hull surface
                ! note: flou is positive downwards
                !
                rx = 0.25*( xcgrid(m,n) + xcgrid(md,n) + xcgrid(m,nd) + xcgrid(md,nd) ) - bcog(l,1)
                !
                rz = -0.5*( flou(nm) + flou(nmd) ) - bcog(l,3)
                !
                dxl = 0.5*( gvv(nm) + gvv(ndm) )
                !
                ! slope of hull
                !
                dsdx = ( flou(nm) - flou(nmd) ) / dxl
                !
                skc(nm) = skc(nm) + vrb * rz * dsdx - vrb * rx
                !
             endif
             !
             ! pitch - contribution next to body (from free surface to pressurized in x-direction)
             !
 180         if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
                !
                l = lfbs(nmu)
                if ( l == 0 .or. .not.bdof(l,2,2) ) goto 190
                !
                ! pitch velocity
                !
                vrb = theta * vfor1(l,2) + (1.-theta) * vfor0(l,2)
                !
                ! position of wetted hull surface
                ! note: flos is positive downwards
                !
                rz = 0.5*( s1(nm) - flos(nmu) ) - bcog(l,3)
                !
                dxl = 0.5*( gvv(nm) + gvv(ndm) )
                !
                ! slope of hull on the port side
                !
                dsdx = -( flos(nmu) + s1(nm) ) / dxl
                !
                skc(nm) = skc(nm) - vrb * rz * dsdx
                !
             endif
             !
             ! pitch - contribution next to body (from pressurized to free surface in x-direction)
             !
 190         if ( presp(nm) == 0 .and. presp(nmd) == 1 ) then
                !
                l = lfbs(nmd)
                if ( l == 0 .or. .not.bdof(l,2,2) ) goto 200
                !
                ! pitch velocity
                !
                vrb = theta * vfor1(l,2) + (1.-theta) * vfor0(l,2)
                !
                ! position of wetted hull surface
                ! note: flos is positive downwards
                !
                rz = 0.5*( s1(nm) - flos(nmd) ) - bcog(l,3)
                !
                dxl = 0.5*( gvv(nm) + gvv(ndm) )
                !
                ! slope of hull on the starboard side
                !
                dsdx = ( s1(nm) + flos(nmd) ) / dxl
                !
                skc(nm) = skc(nm) - vrb * rz * dsdx
                !
             endif
             !
             ! yaw - contribution underneath body
             !
 200         if ( presp(nm) == 1 ) then
                !
                l = lfbs(nm)
                if ( l == 0 .or. .not.bdof(l,3,2) ) goto 210
                !
                ! yaw velocity
                !
                vrb = theta * vfor1(l,3) + (1.-theta) * vfor0(l,3)
                !
                ! position of wetted hull surface
                !
                rx = 0.25*( xcgrid(m,n) + xcgrid(md,n) + xcgrid(m,nd) + xcgrid(md,nd) ) - bcog(l,1)
                !
                ry = 0.25*( ycgrid(m,n) + ycgrid(m,nd) + ycgrid(md,n) + ycgrid(md,nd) ) - bcog(l,2)
                !
                dxl = 0.5*( gvv(nm) + gvv(ndm) )
                dyl = 0.5*( guu(nm) + guu(nmd) )
                !
                ! slope of hull
                ! note: flou/flov is positive downwards
                !
                dsdx = ( flou(nm) - flou(nmd) ) / dxl
                !
                dsdy = ( flov(nm) - flov(ndm) ) / dyl
                !
                skc(nm) = skc(nm) - vrb * ry * dsdx + vrb * rx * dsdy
                !
             endif
             !
             ! yaw - contribution next to body (from free surface to pressurized in x-direction)
             !
 210         if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
                !
                l = lfbs(nmu)
                if ( l == 0 .or. .not.bdof(l,3,2) ) goto 220
                !
                ! yaw velocity
                !
                vrb = theta * vfor1(l,3) + (1.-theta) * vfor0(l,3)
                !
                ! position of wetted hull surface
                !
                ry = 0.5*( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                !
                dxl = 0.5*( gvv(nm) + gvv(ndm) )
                !
                ! slope of hull on the port side
                ! note: flos is positive downwards
                !
                dsdx = -( flos(nmu) + s1(nm) ) / dxl
                !
                skc(nm) = skc(nm) + vrb * ry * dsdx
                !
             endif
             !
             ! yaw - contribution next to body (from pressurized to free surface in x-direction)
             !
 220         if ( presp(nm) == 0 .and. presp(nmd) == 1 ) then
                !
                l = lfbs(nmd)
                if ( l == 0 .or. .not.bdof(l,3,2) ) goto 230
                !
                ! yaw velocity
                !
                vrb = theta * vfor1(l,3) + (1.-theta) * vfor0(l,3)
                !
                ! position of wetted hull surface
                !
                ry = 0.5*( ycgrid(md,n) + ycgrid(md,nd) ) - bcog(l,2)
                !
                dxl = 0.5*( gvv(nm) + gvv(ndm) )
                !
                ! slope of hull on the starboard side
                ! note: flos is positive downwards
                !
                dsdx = ( s1(nm) + flos(nmd) ) / dxl
                !
                skc(nm) = skc(nm) + vrb * ry * dsdx
                !
             endif
             !
             ! yaw - contribution next to body (from free surface to pressurized in y-direction)
             !
 230         if ( presp(nm) == 0 .and. presp(num) == 1 ) then
                !
                l = lfbs(num)
                if ( l == 0 .or. .not.bdof(l,3,2) ) goto 240
                !
                ! yaw velocity
                !
                vrb = theta * vfor1(l,3) + (1.-theta) * vfor0(l,3)
                !
                ! position of wetted hull surface
                !
                rx = 0.5*( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                !
                dyl = 0.5*( guu(nm) + guu(nmd) )
                !
                ! slope of hull on the port side
                ! note: flos is positive downwards
                !
                dsdy = -( flos(num) + s1(nm) ) / dyl
                !
                skc(nm) = skc(nm) - vrb * rx * dsdy
                !
             endif
             !
             ! yaw - contribution next to body (from pressurized to free surface in y-direction)
             !
 240         if ( presp(nm) == 0 .and. presp(ndm) == 1 ) then
                !
                l = lfbs(ndm)
                if ( l == 0 .or. .not.bdof(l,3,2) ) goto 250
                !
                ! yaw velocity
                !
                vrb = theta * vfor1(l,3) + (1.-theta) * vfor0(l,3)
                !
                ! position of wetted hull surface
                !
                rx = 0.5*( xcgrid(m,nd) + xcgrid(md,nd) ) - bcog(l,1)
                !
                dyl = 0.5*( guu(nm) + guu(nmd) )
                !
                ! slope of hull on the starboard side
                ! note: flos is positive downwards
                !
                dsdy = ( s1(nm) + flos(ndm) ) / dyl
                !
                skc(nm) = skc(nm) - vrb * rx * dsdy
                !
             endif
             !
 250         continue
             !
          enddo
       enddo
       !
    endif
    !
end subroutine SwashUpdKBCrigb
