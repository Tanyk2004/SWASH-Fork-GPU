subroutine SwashUpdUDepu ( u )
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
!    1.00, October 2021: New subroutine
!   10.00, January 2023: add high order correction
!
!   Purpose
!
!   Update water depths at cell faces
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata
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
    integer               :: icell    ! cell index
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
    real                  :: htotl    ! water depth in left cell
    real                  :: htotr    ! water depth in right cell
    real                  :: s1min    ! local minimum of water level
    real                  :: totarea  ! total area of all cells around vertex
    real                  :: wlu      ! averaged water level in upwind vertex
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
    if (ltrace) call strace (ient,'SwashUpdUDepu')
    !
    ! point to vertex, cell and face objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    if ( .not.depcds ) then
       !
       ! compute the water depth at faces based on upwinding
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
          hu(iface) = h
          !
       enddo
       !
    else
       !
       ! compute the water depth at faces based on averaging
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
             htotl = s1(icelll) + dps(icelll)
             htotr = s1(icellr) + dps(icellr)
             !
             h = finp * htotl + (1.-finp) * htotr
             !
          else if ( icelll == 0 ) then
             !
             h = s1(icellr) + dps(icellr)
             !
          else
             !
             h = s1(icelll) + dps(icelll)
             !
          endif
          !
          hu(iface) = h
          !
       enddo
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
end subroutine SwashUpdUDepu
