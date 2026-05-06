subroutine SwashUpdateUDepths ( u )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
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
!   10.00: Marcel Zijlema
!
!   Updates
!
!    1.00, February 2020: New subroutine
!   10.00,  January 2023: add high order correction
!
!   Purpose
!
!   Initialize / update water depths in both cells and faces
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata
    use outp_data, only: hrunp
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    real, dimension(nfaces), intent(in) :: u ! flow velocity at faces
!
!   Local variables
!
    integer               :: icell    ! cell index / loop counter over cells
    integer               :: icelll   ! left cell of present face
    integer               :: icellr   ! right cell of present face
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer               :: iface    ! loop counter over faces
    integer               :: jc       ! loop counter
    integer               :: k        ! loop counter
    integer, dimension(3) :: v        ! vertices of present cell
    integer               :: vf1      ! first vertex of present face
    integer               :: vf2      ! second vertex of present face
    integer               :: vu       ! upwind vertex
    !
    real                  :: area     ! area of present cell
    real                  :: depmin   ! local minimum of bottom depth
    real                  :: du       ! averaged bed level in upwind vertex
    real                  :: finp     ! interpolation factor
    real                  :: fluxlim  ! flux limiter
    real                  :: grad1    ! solution gradient
    real                  :: grad2    ! another solution gradient
    real                  :: h        ! local water depth
    real                  :: s1min    ! local minimum of water level
    real                  :: totarea  ! total area of all cells around vertex
    real                  :: wlu      ! averaged water level in upwind vertex
    !
    logical               :: adapted  ! true if value of epsdry has been changed
    !
    character(80)         :: msgstr   ! string to pass message
    !
    type(verttype), dimension(:), pointer :: vert ! datastructure for vertices with their attributes
    type(celltype), dimension(:), pointer :: cell ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdateUDepths')
    !
    ! point to vertex, cell and face objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    ! compute and check water depth in cells
    !
    adapted = .false.
    !
    do icell = 1, ncells
       !
       hs(icell) = s1(icell) + dps(icell)
       !
       if ( hs(icell) < 0. .and. inewt == 0 ) then
          !
          if ( hs(icell) < -epsdry ) then
             !
             write (msgstr,'(a,i8,a,f8.2,a)') 'water depth is negative in cell=',icell,'; water depth = ',1000.*hs(icell),' mm'
             call msgerr (1, trim(msgstr) )
             !
             s1(icell) = 0.99*epsdry - dps(icell)
             epsdry  = -1.01*hs(icell)
             adapted = .true.
             !
          else
             s1(icell) = 0.99*epsdry - dps(icell)
          endif
          !
          hs(icell) = s1(icell) + dps(icell)
          !
       endif
       !
       ! in case of Newton iteration, water depth must be bounded below by bottom level
       !
       if ( inewt /= 0 ) hs(icell) = max( epsdry, hs(icell) )
       !
    enddo
    !
    ! check minimal depth in cells
    !
    if ( inewt == 0 ) then
       !
       if ( adapted ) then
          write (msgstr,'(a,e12.4)') 'new minimal depth for checking drying and flooding: DEPMIN = ',epsdry
          call msgerr (1, trim(msgstr) )
       endif
       !
       ! check minimal depth for drying and flooding (at most 1 cm)
       !
       if ( epsdry > 0.01 ) then
          !
          call msgerr ( 4, 'INSTABLE: water level is too far below the bottom level!' )
          call msgerr ( 0, '          Please reduce the time step!' )
          return
          !
       endif
       !
       ! also reset minimal depth at faces
       !
       epshu = epsdry
       !
    endif
    !
    ! compute the water depth at faces based on averaging
    !
    humo = hum
    !
    do iface = 1, nfaces
       !
       finp = face(iface)%attr(FACELINPF)
       !
       ! consider left and right cells of current face
       !
       icelll = face(iface)%atti(FACECL)
       icellr = face(iface)%atti(FACECR)
       !
       if ( icelll /= 0 .and. icellr /= 0 ) then
          !
          h = finp * hs(icelll) + (1.-finp) * hs(icellr)
          !
       else if ( icelll == 0 ) then
          !
          h = hs(icellr)
          !
       else
          !
          h = hs(icelll)
          !
       endif
       !
       hum(iface) = h
       !
    enddo
    !
    ! extrapolate water depth at midface in time to improve accuracy of momentum-conservative time integration
    !
    humn = 1.5*hum - 0.5*humo
    !
    ! compute the water depth at faces based on upwinding
    !
    if ( .not.depcds ) then
       !
       do iface = 1, nfaces
          !
          ! consider left and right cells of current face
          !
          icelll = face(iface)%atti(FACECL)
          icellr = face(iface)%atti(FACECR)
          !
          if ( icelll == 0 .or. icellr == 0 ) then   ! boundary cell
             !
             if ( icelll /= 0 ) then
                h = s1(icelll) + dpu(iface)
             else
                h = s1(icellr) + dpu(iface)
             endif
             !
          else if ( u(iface) > epsuf ) then
             !
             h = s1(icelll) + dpu(iface)
             !
          else if ( u(iface) < -epsuf ) then
             !
             h = s1(icellr) + dpu(iface)
             !
          else
             !
             h = max( s1(icelll), s1(icellr) ) + dpu(iface)
             !
          endif
          !
          ! in case of Newton iteration, water depth must be bounded below by bottom level
          !
          if ( inewt /= 0 ) h = max( epsdry, h )
          !
          hu(iface) = h
          !
       enddo
       !
    else
       !
       hu = hum
       !
    endif
    !
    ! compute higher order correction to the water depth in cell circumcenter (if appropriate)
    ! Note: for the r-ratio the most upwave vertex of upwind cell is used
    !
    if ( corrdep ) then
       !
       propsc = nint(pnums(11))
       kappa  = pnums(12)
       mbound = pnums(13)
       phieby = pnums(14)
       !
       floop: do iface = 1, nfaces
          !
          ! get vertices of current face
          !
          vf1 = face(iface)%atti(FACEV1)
          vf2 = face(iface)%atti(FACEV2)
          !
          ! consider left and right cells of current face
          !
          icelll = face(iface)%atti(FACECL)
          icellr = face(iface)%atti(FACECR)
          !
          if ( icelll == 0 .or. icellr == 0 ) then
             !
             ! no correction in boundary cell
             !
          else if ( u(iface) > epsuf ) then
             !
             ! get vertices of upwind cell
             !
             v(1) = cell(icelll)%atti(CELLV1)
             v(2) = cell(icelll)%atti(CELLV2)
             v(3) = cell(icelll)%atti(CELLV3)
             !
             ! search for most upwave vertex
             !
             do k = 1, 3
                if ( v(k) /= vf1 .and. v(k) /= vf2 ) then
                   vu = v(k)
                   exit
                endif
             enddo
             !
             ! compute area-weighted averaged water and bed levels at upwave vertex
             !
             wlu     = 0.
             du      = 0.
             totarea = 0.
             !
             do jc = 1, vert(vu)%noc
                !
                icell = vert(vu)%cell(jc)%atti(CELLID)
                !
                area = cell(icell)%attr(CELLAREA)
                !
                wlu = wlu + area * s1 (icell)
                du  = du  + area * dps(icell)
                !
                totarea = totarea + area
                !
             enddo
             !
             wlu = wlu / totarea
             du  = du  / totarea
             !
             depmin = min( du , dps(icelll), dps(icellr) )
             s1min  = min( wlu, s1 (icelll), s1 (icellr) )
             !
             if ( s1min + depmin < 0. ) cycle floop
             !
             ! compute solution gradients
             !
             grad1 = s1(icellr) - s1(icelll)
             grad2 = s1(icelll) - wlu
             !
             ! update water depth
             !
             hu(iface) = hu(iface) + 0.5 * fluxlim(grad1,grad2)
             !
          else if ( u(iface) < -epsuf ) then
             !
             ! get vertices of upwind cell
             !
             v(1) = cell(icellr)%atti(CELLV1)
             v(2) = cell(icellr)%atti(CELLV2)
             v(3) = cell(icellr)%atti(CELLV3)
             !
             ! search for most upwave vertex
             !
             do k = 1, 3
                if ( v(k) /= vf1 .and. v(k) /= vf2 ) then
                   vu = v(k)
                   exit
                endif
             enddo
             !
             ! compute area-weighted averaged water and bed levels at upwave vertex
             !
             wlu     = 0.
             du      = 0.
             totarea = 0.
             !
             do jc = 1, vert(vu)%noc
                !
                icell = vert(vu)%cell(jc)%atti(CELLID)
                !
                area = cell(icell)%attr(CELLAREA)
                !
                wlu = wlu + area * s1 (icell)
                du  = du  + area * dps(icell)
                !
                totarea = totarea + area
                !
             enddo
             !
             wlu = wlu / totarea
             du  = du  / totarea
             !
             depmin = min( dps(icelll), dps(icellr), du  )
             s1min  = min( s1 (icelll), s1 (icellr), wlu )
             !
             if ( s1min + depmin < 0. ) cycle floop
             !
             ! compute solution gradients
             !
             grad1 = s1(icelll) - s1(icellr)
             grad2 = s1(icellr) - wlu
             !
             ! update water depth
             !
             hu(iface) = hu(iface) + 0.5 * fluxlim(grad1,grad2)
             !
          endif
          !
       enddo floop
       !
    endif
    !
    ! compute inundation depth
    !
    do icell = 1, ncells
       if ( hindun(icell) /= 1. .and. hs(icell) > hrunp ) hindun(icell) = 1.
    enddo
    !
end subroutine SwashUpdateUDepths
