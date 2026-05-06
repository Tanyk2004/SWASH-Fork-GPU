subroutine SwashImpDepUflow ( u1, u0, uvc, qn, quf, q, dq, gmat, rho )
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
!   10.01: Marcel Zijlema
!
!   Updates
!
!    1.00, January 2020: New subroutine
!   10.01, January 2023: add higher order advection
!   10.01,   April 2023: Newton-type iteration added
!
!   Purpose
!
!   Performs the time integration for the non-hydrostatic, depth-averaged shallow water equations on triangular mesh
!
!   Method
!
!   The time integration with respect to the continuity equation and the water level gradient of the
!   u-momentum equation is based on a theta-scheme. Only a value of 0.5 <= theta <= 1 will be taken.
!
!   The time integration with respect to the advective and viscous terms and bottom friction is based
!   on Euler implicit.
!
!   The space discretization of the advective term is momentum conservative and is approximated by
!   first order upwind or higher order (flux-limited) scheme (CDS, Fromm, BDF, QUICK, MUSCL, Koren, etc.).
!   The r-ratio formulation based on most upwave vertex is employed.
!
!   The w-momentum equation only contains the z-gradient of the non-hydrostatic pressure and is
!   discretized by means of the Keller-box scheme.
!
!   The non-hydrostatic pressure is obtained by means of the second order accurate pressure correction technique.
!
!   Modules used
!
    use ocpcomm4
    use m_genarr, only: wlimp
    use SwashCommdata3
    use SwashTimecomm
    use SwashFlowdata, rhotmp => rho, &
                       u1tmp  => u1 , &
                       u0tmp  => u0 , &
                       uvctmp => uvc, &
                       qtmp   => q  , &
                       qntmp  => qn , &
                       quftmp => quf, &
                       dqtmp  => dq
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    real, dimension(ncells)  , intent(out)   :: dq     ! non-hydrostatic pressure correction
    real, dimension(nfaces,2), intent(out)   :: gmat   ! gradient matrix for pressure at face
    real, dimension(ncells)  , intent(inout) :: q      ! non-hydrostatic pressure
    real, dimension(nfaces)  , intent(out)   :: qn     ! mass flux at triangular face
    real, dimension(ncells,3), intent(out)   :: quf    ! mass flux at cell face in flow direction
    real, dimension(ncells)  , intent(in)    :: rho    ! density of water
    real, dimension(nfaces)  , intent(in)    :: u0     ! face velocity at previous time level
    real, dimension(nfaces)  , intent(inout) :: u1     ! face velocity at current time level
    real, dimension(ncells,2), intent(out)   :: uvc    ! velocity components in cell circumcenter
!
!   Local variables
!
    integer                                :: bcell    ! boundary cell index for sponge layer
    integer                                :: bface    ! boundary face index for sponge layer
    integer                                :: icell    ! cell index / loop counter over cells
    integer                                :: icella   ! cell index for computing area-weighted averaged quantity
    integer                                :: icellb   ! boundary cell
    integer                                :: icelld   ! cell index of downwind cell
    integer                                :: icelll   ! left cell of present face
    integer                                :: icellr   ! right cell of present face
    integer                                :: icellu   ! cell index of upwind cell
    integer, save                          :: ient = 0 ! number of entries in this subroutine
    integer                                :: iface    ! face index / loop counter over faces
    integer                                :: ifacel   ! local face
    integer                                :: j        ! loop counter
    integer                                :: jc       ! loop counter
    integer                                :: jf       ! loop counter
    integer                                :: l        ! loop counter
    integer, dimension(3)                  :: vc       ! vertices of present cell
    integer                                :: vf1      ! first vertex of present face
    integer                                :: vf2      ! second vertex of present face
    integer                                :: vu       ! upwind vertex
    !
    real                                   :: area     ! area of present cell
    real                                   :: cadv     ! contribution to advection term
    real                                   :: cfl      ! CFL number
    real                                   :: contrib  ! contribution to various parts of momentum equation
    real                                   :: ctrs     ! contribution to system of water level equations
    real                                   :: cvisc    ! contribution to viscosity term
    real                                   :: dbs      ! change in water level at boundary face over time
    real                                   :: denom    ! a denominator
    real                                   :: dpx      ! x-component of bottom gradient
    real                                   :: dpy      ! y-component of bottom gradient
    real                                   :: ener     ! total energy of closed system
    real                                   :: epsab2   ! small amount to modify weights of AB2 scheme
    real                                   :: fac      ! a factor
    real                                   :: faca     ! contribution weight to left/right cells for advection/viscosity term
    real                                   :: finp     ! interpolation factor
    real                                   :: fluxlim  ! flux limiter
    real                                   :: gamma    ! relaxation parameter for sponge layer
    real                                   :: grad1    ! solution gradient
    real                                   :: grad1x   ! x-component of solution gradient
    real                                   :: grad1y   ! y-component of solution gradient
    real                                   :: grad2    ! another solution gradient
    real                                   :: grad2x   ! x-component of another solution gradient
    real                                   :: grad2y   ! y-component of another solution gradient
    real                                   :: hf       ! arithmetic average of water depth at face
    real                                   :: hox      ! x-component of high order limited contribution to flow velocity
    real                                   :: hoy      ! y-component of high order limited contribution to flow velocity
    real                                   :: lf       ! length of face
    real                                   :: lwfac    ! Lax-Wendroff factor
    real                                   :: moutf    ! net mass outflow
    real                                   :: nx       ! x-component of normal to face
    real                                   :: ny       ! y-component of normal to face
    real                                   :: qf       ! (total) mass flux
    real                                   :: rdx      ! reciprocal of distance between circumcenters adjacent to face
    real                                   :: rsgn     ! sign for indicating face orientation
    real                                   :: sumqf    ! sum of outgoing mass fluxes over present cell
    real                                   :: theta    ! implicitness factor for time integration of continuity equation
    real                                   :: theta2   ! implicitness factor for water level gradient
    real                                   :: theta3   ! implicitness factor for non-hydrostatic pressure gradient
    real                                   :: totarea  ! total area of all cells around vertex
    real                                   :: uf       ! updated flow velocity at present face
    real                                   :: ut       ! tangential velocity component
    real                                   :: utot     ! velocity magnitude
    real                                   :: ux       ! x-component of flow velocity in cell circumenter
    real                                   :: uy       ! y-component of flow velocity in cell circumenter
    real                                   :: vol      ! total displaced volume of water
    real                                   :: w0u      ! averaged w-velocity in upwind vertex
    real                                   :: xc       ! x-coordinate of cell circumcenter
    real                                   :: xf       ! x-coordinate of face center
    real                                   :: yc       ! y-coordinate of cell circumcenter
    real                                   :: yf       ! y-coordinate of face center
    real                                   :: zgrad    ! water level gradient term
    !
    logical                                :: STPNOW   ! indicates that program must stop
    !
    type(verttype), dimension(:), pointer  :: vert     ! datastructure for vertices with their attributes
    type(celltype), dimension(:), pointer  :: cell     ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer  :: face     ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashImpDepUflow')
    !
    ! point to vertex, cell and face objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    theta  = pnums(1)
    theta2 = pnums(4)
    theta3 = pnums(5)
    !
    ! compute the mass flux
    !
    do icell = 1, ncells
       !
       ! loop over faces of the cell
       !
       do jf = 1, cell(icell)%nof
          !
          ! face identifier
          !
          iface = cell(icell)%face(jf)%atti(FACEID)
          !
          ! get length of current face
          !
          lf = face(iface)%attr(FACELEN)
          !
          ! consider left and right cells of current face
          !
          icelll = face(iface)%atti(FACECL)
          icellr = face(iface)%atti(FACECR)
          !
          ! take into account orientation of the current face
          !
          if ( icell == icelll ) then
             rsgn =  1.
          else if ( icell == icellr ) then
             rsgn = -1.
          endif
          !
          ! compute mass flux at current face
          !
          quf(icell,jf) = rsgn * lf * hu(iface) * u0(iface)
          !
       enddo
       !
    enddo
    !
    ! build gradient matrix for non-hydrostatic pressure
    !
    if ( ihydro == 1 ) then
       !
       do iface = 1, nfaces
          !
          rdx = face(iface)%attr(FACEDISTC)
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             gmat(iface,1) = -rdx * ( (0.5-finp)*dps(icelll) + 0.5*s0(icelll) +      finp * dps(icellr) )
             gmat(iface,2) =  rdx * ( (finp-0.5)*dps(icellr) + 0.5*s0(icellr) + (1.-finp) * dps(icelll) )
             !
          else if ( wlimp(iface) .and. wetu(iface) == 1 ) then   ! described water level at boundary face
             !
             ! consider boundary cell of current face
             !
             icellb = face(iface)%atti(FACEC1)
             !
             gmat(iface,1) = -rdx * hs(icellb)
             gmat(iface,2) = 0.
             !
          else
             !
             gmat(iface,1) = 0.
             gmat(iface,2) = 0.
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute gradient of non-hydrostatic pressure in momentum equation
    !
    if ( ihydro == 1 ) then
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             qgrad(iface) = gmat(iface,1) * q(icelll) + gmat(iface,2) * q(icellr)
             !
          else if ( wlimp(iface) .and. wetu(iface) == 1 ) then   ! described water level at boundary face
             !
             ! consider boundary cell of current face
             !
             icellb = face(iface)%atti(FACEC1)
             !
             qgrad(iface) = gmat(iface,1) * q(icellb)
             !
          endif
          !
       enddo
       !
       if ( iproj == 2 ) qgrad = (1.-theta3) * qgrad
       !
    else
       !
       qgrad = 0.
       !
    endif
    !
    ! compute cell-based velocity vector
    !
    call perot ( u0, 1, 1 )
    !
    ! compute bottom friction term, if appropriate
    !
    if ( irough == 11 ) then
       !
       ! linear bottom friction
       !
       do iface = 1, nfaces
          !
          ! consider left and right cells of current face
          !
          icelll = face(iface)%atti(FACECL)
          icellr = face(iface)%atti(FACECR)
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             cbot(iface) = ( finp * cfricu(icelll) + (1.-finp) * cfricu(icellr) ) / hum(iface)
             !
          elseif ( wetu(iface) == 1 ) then   ! wet boundary face
             !
             if ( icelll == 0 ) then
                !
                cbot(iface) = cfricu(icellr) / hum(iface)
                !
             elseif ( icellr == 0 ) then
                !
                cbot(iface) = cfricu(icelll) / hum(iface)
                !
             endif
             !
          else
             !
             cbot = 0.
             !
          endif
          !
       enddo
       !
    else if ( irough /= 0 ) then
       !
       ! quadratic bottom friction
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             cbot(iface) = 0.
             !
             do l = 1, 2
                !
                if ( l == 1 ) then
                   !
                   ! left cell of current face
                   !
                   icell = face(iface)%atti(FACECL)
                   faca  = finp
                   !
                else
                   !
                   ! right cell of current face
                   !
                   icell = face(iface)%atti(FACECR)
                   faca  = 1. - finp
                   !
                endif
                !
                ! compute velocity magnitude
                !
                utot = sqrt( uvc(icell,1)*uvc(icell,1) + uvc(icell,2)*uvc(icell,2) )
                !
                cbot(iface) = cbot(iface) + faca * cfricu(icell) * utot / hum(iface)
                !
             enddo
             !
          elseif ( wetu(iface) == 1 ) then   ! wet boundary face
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             if ( icelll == 0 ) then
                !
                ! compute velocity magnitude
                !
                utot = sqrt( uvc(icellr,1)*uvc(icellr,1) + uvc(icellr,2)*uvc(icellr,2) )
                !
                cbot(iface) = cfricu(icellr) * utot / hum(iface)
                !
             elseif ( icellr == 0 ) then
                !
                ! compute velocity magnitude
                !
                utot = sqrt( uvc(icelll,1)*uvc(icelll,1) + uvc(icelll,2)*uvc(icelll,2) )
                !
                cbot(iface) = cfricu(icelll) * utot / hum(iface)
                !
             endif
             !
          else
             !
             cbot = 0.
             !
          endif
          !
       enddo
       !
    else
       !
       cbot = 0.
       !
    endif
    !
    ! compute advection term using first order upwind (momentum conservative)
    !
    do iface = 1, nfaces
       !
       if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
          !
          finp = face(iface)%attr(FACELINPF)
          !
          nx = face(iface)%attr(FACENORMX)
          ny = face(iface)%attr(FACENORMY)
          !
          advec(iface) = 0.
          !
          do l = 1, 2
             !
             if ( l == 1 ) then
                !
                ! left cell of current face
                !
                icell = face(iface)%atti(FACECL)
                faca  = finp
                !
             else
                !
                ! right cell of current face
                !
                icell = face(iface)%atti(FACECR)
                faca  = 1. - finp
                !
             endif
             !
             cadv = 0.
             !
             ! loop over faces of the cell under consideration
             !
             do jf = 1, cell(icell)%nof
                !
                ! local face
                !
                ifacel = cell(icell)%face(jf)%atti(FACEID)
                !
                if ( face(ifacel)%atti(FMARKER) == 0 .and. wetu(ifacel) == 1 ) then
                   !
                   ! consider left and right cells of local face
                   !
                   icelll = face(ifacel)%atti(FACECL)
                   icellr = face(ifacel)%atti(FACECR)
                   !
                   ! take upwind cell of local face
                   !
                   if ( u0(ifacel) > 0. ) then
                      icellu = icelll
                   else
                      icellu = icellr
                   endif
                   !
                   ! compute contribution to advection term - only update for ingoing flux
                   !
                   if ( quf(icell,jf) < 0. ) cadv = cadv + quf(icell,jf) * ( nx * uvc(icellu,1) + ny * uvc(icellu,2) - u0(iface) )
                   !
                endif
                !
             enddo
             !
             area = cell(icell)%attr(CELLAREA)
             !
             advec(iface) = advec(iface) + faca * cadv / area
             !
          enddo
          !
       else if ( face(iface)%atti(FBTYPE) == 2 ) then
          !
          nx = face(iface)%attr(FACENORMX)
          ny = face(iface)%attr(FACENORMY)
          !
          ! consider left and right cells of boundary face and ...
          !
          icelll = face(iface)%atti(FACECL)
          icellr = face(iface)%atti(FACECR)
          !
          ! ... boundary cell
          !
          icellb = face(iface)%atti(FACEC1)
          !
          ! get orientation at boundary face
          !
          if ( icellb == icelll ) then
             rsgn =  1.
          else if ( icellb == icellr ) then
             rsgn = -1.
          endif
          !
          if ( rsgn * u0(iface) > 0. ) then
             !
             ! outflow
             !
             advec(iface) = 0.
             !
             ! loop over faces of the boundary cell
             !
             do jf = 1, cell(icellb)%nof
                !
                ! local face
                !
                ifacel = cell(icellb)%face(jf)%atti(FACEID)
                !
                if ( face(ifacel)%atti(FMARKER) == 0 .and. wetu(ifacel) == 1 ) then
                   !
                   ! take upwind cell of local face
                   !
                   if ( u0(ifacel) > 0. ) then
                      icellu = face(ifacel)%atti(FACECL)
                   else
                      icellu = face(ifacel)%atti(FACECR)
                   endif
                   !
                   ! update for ingoing flux
                   !
                   if ( quf(icellb,jf) < 0. ) advec(iface) = advec(iface) + quf(icellb,jf) * ( nx * uvc(icellu,1) + ny * uvc(icellu,2) - u0(iface) )
                   !
                endif
                !
             enddo
             !
             area = cell(icellb)%attr(CELLAREA)
             !
             advec(iface) = advec(iface) / area
             !
          endif
          !
       endif
       !
    enddo
    !
    ! add second order approximation based on Lax-Wendroff method (if appropriate)
    !
    propsc = nint(pnums(6))
    !
    if ( propsc /= 1 ) then
       !
       kappa  = pnums(7)
       mbound = pnums(8)
       phieby = pnums(9)
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             nx = face(iface)%attr(FACENORMX)
             ny = face(iface)%attr(FACENORMY)
             !
             do l = 1, 2
                !
                if ( l == 1 ) then
                   !
                   ! left cell of current face
                   !
                   icell = face(iface)%atti(FACECL)
                   faca  = finp
                   !
                else
                   !
                   ! right cell of current face
                   !
                   icell = face(iface)%atti(FACECR)
                   faca  = 1. - finp
                   !
                endif
                !
                cadv = 0.
                !
                ! loop over faces of the cell under consideration
                !
                do jf = 1, cell(icell)%nof
                   !
                   ! local face
                   !
                   ifacel = cell(icell)%face(jf)%atti(FACEID)
                   !
                   if ( face(ifacel)%atti(FMARKER) == 0 .and. wetu(ifacel) == 1 ) then
                      !
                      ! compute local Lax-Wendroff factor
                      !
                      rdx   = face(ifacel)%attr(FACEDISTC)
                      lwfac = 1. - rdx * dt * abs(u0(ifacel))
                      !
                      ! get vertices of face under consideration
                      !
                      vf1 = face(ifacel)%atti(FACEV1)
                      vf2 = face(ifacel)%atti(FACEV2)
                      !
                      ! consider up- and downwind cells of local face
                      !
                      if ( u0(ifacel) > 0. ) then
                         icellu = face(ifacel)%atti(FACECL)
                         icelld = face(ifacel)%atti(FACECR)
                      else
                         icellu = face(ifacel)%atti(FACECR)
                         icelld = face(ifacel)%atti(FACECL)
                      endif
                      !
                      ! get vertices of upwind cell
                      !
                      vc(1) = cell(icellu)%atti(CELLV1)
                      vc(2) = cell(icellu)%atti(CELLV2)
                      vc(3) = cell(icellu)%atti(CELLV3)
                      !
                      ! search for most upwave vertex
                      !
                      do j = 1, 3
                         if ( vc(j) /= vf1 .and. vc(j) /= vf2 ) then
                            vu = vc(j)
                            exit
                         endif
                      enddo
                      !
                      ! compute area-weighted averaged velocity components at upwave vertex
                      !
                      ux      = 0.
                      uy      = 0.
                      totarea = 0.
                      !
                      do jc = 1, vert(vu)%noc
                         !
                         icella = vert(vu)%cell(jc)%atti(CELLID)
                         !
                         area = cell(icella)%attr(CELLAREA)
                         !
                         ux = ux + area * uvc(icella,1)
                         uy = uy + area * uvc(icella,2)
                         !
                         totarea = totarea + area
                         !
                      enddo
                      !
                      ux = ux / totarea
                      uy = uy / totarea
                      !
                      ! compute solution gradients
                      !
                      grad1x = uvc(icelld,1) - uvc(icellu,1)
                      grad2x = uvc(icellu,1) - ux
                      !
                      grad1y = uvc(icelld,2) - uvc(icellu,2)
                      grad2y = uvc(icellu,2) - uy
                      !
                      ! compute high order limited part
                      !
                      hox = 0.5 * max(0.,lwfac) * fluxlim(grad1x,grad2x)
                      hoy = 0.5 * max(0.,lwfac) * fluxlim(grad1y,grad2y)
                      !
                      cadv = cadv + quf(icell,jf) * ( nx * hox + ny * hoy )
                      !
                   endif
                   !
                enddo
                !
                area = cell(icell)%attr(CELLAREA)
                !
                advec(iface) = advec(iface) + faca * cadv / area
                !
             enddo
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute viscosity term
    !
    if ( ihvisc /= 0 ) then
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             nx = face(iface)%attr(FACENORMX)
             ny = face(iface)%attr(FACENORMY)
             !
             visc(iface) = 0.
             !
             do l = 1, 2
                !
                if ( l == 1 ) then
                   !
                   ! left cell of current face
                   !
                   icell = face(iface)%atti(FACECL)
                   faca  = finp
                   !
                else
                   !
                   ! right cell of current face
                   !
                   icell = face(iface)%atti(FACECR)
                   faca  = 1. - finp
                   !
                endif
                !
                cvisc = 0.
                !
                ! loop over faces of the cell under consideration
                !
                do jf = 1, cell(icell)%nof
                   !
                   ! local face
                   !
                   ifacel = cell(icell)%face(jf)%atti(FACEID)
                   !
                   if ( face(ifacel)%atti(FMARKER) == 0 .and. wetu(ifacel) == 1 ) then
                      !
                      ! get length and normal distance of local face
                      !
                      lf  = face(ifacel)%attr(FACELEN)
                      rdx = face(ifacel)%attr(FACEDISTC)
                      !
                      ! consider left and right cells of local face
                      !
                      icelll = face(ifacel)%atti(FACECL)
                      icellr = face(ifacel)%atti(FACECR)
                      !
                      ! take into account orientation of local face
                      !
                      if ( icell == icelll ) then
                         rsgn =  1.
                      else if ( icell == icellr ) then
                         rsgn = -1.
                      endif
                      !
                      ! compute contribution to viscosity term
                      !
                      cvisc = cvisc + rsgn * vnu2d(ifacel) * hum(ifacel) * lf * rdx * ( nx * ( uvc(icellr,1) - uvc(icelll,1) ) + ny * ( uvc(icellr,2) - uvc(icelll,2) ) )
                      !
                   endif
                   !
                enddo
                !
                area = cell(icell)%attr(CELLAREA)
                !
                visc(iface) = visc(iface) + faca * cvisc / area
                !
             enddo
             !
          endif
          !
       enddo
       !
    else
       !
       visc = 0.
       !
    endif
    !
    ! compute implicit part of wind stress term, if appropriate
    !
    if ( relwnd ) then
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             wndimp(iface) = cwndu(iface) / hum(iface)
             !
          endif
          !
       enddo
       !
    else
       !
       wndimp = 0.
       !
    endif
    !
    ! compute baroclinic forcing
    !
    if ( idens /= 0 ) then
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
             !
             if ( wetu(iface) == 1 ) then
                !
                rdx = face(iface)%attr(FACEDISTC)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                pgrad(iface) = 0.5 * grav * rdx * hum(iface) * hum(iface) * (rho(icellr) - rho(icelll)) / rhow
                !
             else
                !
                pgrad(iface) = 0.
                !
             endif
             !
          endif
          !
       enddo
       !
    else
       !
       pgrad = 0.
       !
    endif
    !
    ! compute atmospheric pressure gradient
    !
    if ( svwp ) then
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
             !
             if ( wetu(iface) == 1 ) then
                !
                rdx = face(iface)%attr(FACEDISTC)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                pgrad(iface) = pgrad(iface) + rdx * hum(iface) * (patm(icellr) - patm(icelll)) / rhow
                !
             endif
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute Coriolis force
    !
    if ( coriolis ) then
       !
       epsab2 = pcor(2)
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             nx = face(iface)%attr(FACENORMX)
             ny = face(iface)%attr(FACENORMY)
             !
             ! at time level n-1
             !
             ! --- modified quasi second order Adams-Bashforth scheme which is conditionally stable for the inviscid case (see Marshall et al., 1997)
             !
             corf(iface) = -(0.5+epsab2) * cfu(iface,1)
             !
             ! at time level n
             !
             cfu(iface,1) = 0.
             !
             do l = 1, 2
                !
                if ( l == 1 ) then
                   !
                   ! left cell of current face
                   !
                   icell = face(iface)%atti(FACECL)
                   faca  = finp
                   !
                else
                   !
                   ! right cell of current face
                   !
                   icell = face(iface)%atti(FACECR)
                   faca  = 1. - finp
                   !
                endif
                !
                ! get coordinates of the cell circumcenter and cell area
                !
                xc = cell(icell)%attr(CELLCCX)
                yc = cell(icell)%attr(CELLCCY)
                !
                area = cell(icell)%attr(CELLAREA)
                !
                ! compute cell-centered depth-integrated velocity components
                !
                ux = 0.
                uy = 0.
                !
                do jf = 1, cell(icell)%nof
                   !
                   ! local face
                   !
                   ifacel = cell(icell)%face(jf)%atti(FACEID)
                   !
                   ! get length of local face
                   !
                   lf = face(ifacel)%attr(FACELEN)
                   !
                   ! get coordinates of the midface
                   !
                   xf = face(ifacel)%attr(FACEMX)
                   yf = face(ifacel)%attr(FACEMY)
                   !
                   ! consider left and right cells of local face
                   !
                   icelll = face(ifacel)%atti(FACECL)
                   icellr = face(ifacel)%atti(FACECR)
                   !
                   ! take into account orientation of local face
                   !
                   if ( icell == icelll ) then
                      rsgn =  1.
                   else if ( icell == icellr ) then
                      rsgn = -1.
                   endif
                   !
                   ! compute the depth-integrated velocity vector at current cell using Perot's formula
                   !
                   ux = ux + rsgn * lf * ( xf - xc ) * hu(ifacel) * u0(ifacel)
                   uy = uy + rsgn * lf * ( yf - yc ) * hu(ifacel) * u0(ifacel)
                   !
                enddo
                !
                cfu(iface,1) = cfu(iface,1) + faca * fcor(icell,1) * ( ny * ux - nx * uy ) / area
                !
             enddo
             !
             corf(iface) = corf(iface) + (1.5+epsab2) * cfu(iface,1)
             !
          endif
          !
       enddo
       !
    else
       !
       corf = 0.
       !
    endif
    !
    ! compute the flow velocity
    !
    ! - loop over boundary faces where the water level is prescribed, if appropriate
    !
    do iface = 1, nfaces
       !
       if ( wlimp(iface) ) then   ! described water level at boundary face
          !
          if ( wetu(iface) == 1 ) then
             !
             rdx = face(iface)%attr(FACEDISTC)
             !
             ! consider boundary cell of current face
             !
             icellb = face(iface)%atti(FACEC1)
             !
             ! compute water depth at boundary face
             !
             hf = 0.5 * ( s0(icellb) + bcso(iface) ) + dps(icellb)
             !
             ! compute water level gradient
             !
             zgrad = pfac(iface) * grav * rdx * hf * ( bcso(iface) - s0(icellb) )
             !
             ! compute total contributions of the momentum equation
             !
             contrib = advec(iface) / humn(iface) + ( zgrad + qgrad(iface) ) / hum(iface)
             !
             denom = 1. + dt * cbot(iface)
             !
             ! update normal component of flow velocity at boundary face
             !
             u1(iface) = ( u0(iface) - dt * contrib ) / denom
             !
          else
             !
             u1(iface) = 0.
             !
          endif
          !
       endif
       !
    enddo
    !
    ! - loop over internal faces
    !
    do iface = 1, nfaces
       !
       if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
          !
          if ( wetu(iface) == 1 ) then
             !
             rdx = face(iface)%attr(FACEDISTC)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             ! compute arithmetic average water depth at present face
             !
             hf = 0.5 * ( hs(icelll) + hs(icellr) )
             !
             ! compute water level gradient
             !
             zgrad = grav * rdx * hf * ( s0(icellr) - s0(icelll) )
             !
             ! compute total contributions of the momentum equation
             !
             contrib = advec(iface) / humn(iface) + ( zgrad + pgrad(iface) + qgrad(iface) - windu(iface) - visc(iface) + corf(iface) ) / hum(iface)
             !
             denom = 1. + dt * ( cbot(iface) + wndimp(iface) )
             !
             ! update flow velocity at face
             !
             u1(iface) = ( u0(iface) - dt * contrib ) / denom
             !
          else
             !
             u1(iface) = 0.
             !
          endif
          !
       endif
       !
    enddo
    !
    ! compute the water level
    !
    fac = grav * theta * theta2 * dt * dt
    !
    ! first, build the equation for water level correction
    !
    if ( inewt == 0 ) then
       !
       do icell = 1, ncells
          !
          amat(icell,0) = cell(icell)%attr(CELLAREA)
          rhs (icell  ) = 0.
          !
          ! loop over faces of the cell
          !
          do jf = 1, cell(icell)%nof
             !
             ! face identifier
             !
             iface = cell(icell)%face(jf)%atti(FACEID)
             !
             lf  = face(iface)%attr(FACELEN)
             rdx = face(iface)%attr(FACEDISTC)
             !
             ! get the necessary cells
             !
             icellb = face(iface)%atti(FACEC1)
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             ! build matrix and right-hand side
             !
             if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
                !
                hf             = 0.5 * ( hs(icelll) + hs(icellr) )
                ctrs           = -fac * lf * hu(iface) * rdx * hf / hum(iface)
                amat(icell,jf) = ctrs
                !
             else if ( wlimp(iface) ) then   ! described water level at boundary face
                !
                hf         = 0.5 * ( s0(icellb) + bcso(iface) ) + dps(icellb)
                ctrs       = -pfac(iface) * fac * lf * hu(iface) * abs(rdx) * hf / hum(iface)
                dbs        = bcs(iface) - bcso(iface)
                rhs(icell) = rhs(icell) - ctrs * dbs
                !
             else
                !
                ctrs = 0.
                !
             endif
             !
             amat(icell,0) = amat(icell,0) - ctrs
             !
             ! take into account orientation of the current face
             !
             if ( icell == icelll ) then
                rsgn =  1.
             else if ( icell == icellr ) then
                rsgn = -1.
             endif
             !
             rhs(icell) = rhs(icell) - dt * rsgn * lf * hu(iface) * ( theta*u1(iface) + (1.-theta)*u0(iface) )
             !
          enddo
          !
       enddo
       !
       ! add mass source due to internal wave generation
       !
       if ( iwvgen /= 0 ) then
          !
          do icell = 1, ncells
             !
             area = cell(icell)%attr(CELLAREA)
             !
             rhs(icell) = rhs(icell) + dt * area * srcm(icell)
             !
          enddo
          !
       endif
       !
       ! next, solve the equation for water level correction
       !
       call pcgu( amat, rhs, ds )
       if (STPNOW()) return
       !
    else
       !
       do icell = 1, ncells
          !
          amatn(icell,0) = 0.
          rhsn (icell  ) = 0.
          !
          ! loop over faces of the cell
          !
          do jf = 1, cell(icell)%nof
             !
             ! face identifier
             !
             iface = cell(icell)%face(jf)%atti(FACEID)
             !
             lf  = face(iface)%attr(FACELEN)
             rdx = face(iface)%attr(FACEDISTC)
             !
             ! get the necessary cells
             !
             icellb = face(iface)%atti(FACEC1)
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             ! build matrix and right-hand side
             !
             if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
                !
                hf              = 0.5 * ( hs(icelll) + hs(icellr) )
                ctrs            = -fac * lf * hu(iface) * rdx * hf / hum(iface)
                amatn(icell,jf) = ctrs
                !
             else if ( wlimp(iface) ) then   ! described water level at boundary face
                !
                hf          = 0.5 * ( s0(icellb) + bcso(iface) ) + dps(icellb)
                ctrs        = -pfac(iface) * fac * lf * hu(iface) * abs(rdx) * hf / hum(iface)
                dbs         = bcs(iface) - bcso(iface)
                rhsn(icell) = rhsn(icell) - ctrs * dbs
                !
             else
                !
                ctrs = 0.
                !
             endif
             !
             amatn(icell,0) = amatn(icell,0) - ctrs
             !
             ! take into account orientation of the current face
             !
             if ( icell == icelll ) then
                rsgn =  1.
             else if ( icell == icellr ) then
                rsgn = -1.
             endif
             !
             rhsn(icell) = rhsn(icell) - dt * rsgn * lf * hu(iface) * ( theta*u1(iface) + (1.-theta)*u0(iface) )
             !
          enddo
          !
          area = cell(icell)%attr(CELLAREA)
          !
          rhsn(icell) = rhsn(icell) + area * hs(icell)
          !
       enddo
       !
       ! check sum(rhs) > 0, important for convergence of Newton-type iteration
       !
       if ( .not. sum(rhsn) > 0d0 ) then
          !
          call msgerr ( 2, 'total water volume is negative!')
          call msgerr ( 0, 'It is advised to reduce the time step!' )
          !
       endif
       !
       ! add mass source due to internal wave generation
       !
       if ( iwvgen /= 0 ) then
          !
          do icell = 1, ncells
             !
             area = cell(icell)%attr(CELLAREA)
             !
             rhsn(icell) = rhsn(icell) + dt * area * srcm(icell)
             !
          enddo
          !
       endif
       !
       ! next, solve the equation for water level correction
       !
       call newtonU ( amatn, rhsn, ds )
       if (STPNOW()) return
       !
    endif
    !
    ! correct flow velocity
    !
    do iface = 1, nfaces
       !
       if ( wlimp(iface) ) then   ! described water level at boundary face
          !
          if ( wetu(iface) == 1 ) then
             !
             rdx = face(iface)%attr(FACEDISTC)
             !
             ! consider boundary cell of current face
             !
             icellb = face(iface)%atti(FACEC1)
             !
             ! compute water depth at boundary face
             !
             hf = 0.5 * ( s0(icellb) + bcso(iface) ) + dps(icellb)
             !
             ! get change in boundary water level over time
             !
             dbs = bcs(iface) - bcso(iface)
             !
             u1(iface) = u1(iface) - pfac(iface) * grav * theta2 * dt * rdx * hf * ( dbs - ds(icellb) ) / hum(iface)
             !
          endif
          !
       endif
       !
    enddo
    !
    do iface = 1, nfaces
       !
       if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
          !
          if ( wetu(iface) == 1 ) then
             !
             rdx = face(iface)%attr(FACEDISTC)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             ! compute arithmetic average water depth at present face
             !
             hf = 0.5 * ( hs(icelll) + hs(icellr) )
             !
             u1(iface) = u1(iface) - grav * theta2 * dt * rdx * hf * ( ds(icellr) - ds(icelll) ) / hum(iface)
             !
          endif
          !
       endif
       !
    enddo
    !
    ! update water level
    !
    s1 = s0 + ds
    !
    ! check positivity of the water depth in each cell of the domain, if appropriate
    !
    if ( inewt == 0 ) then
       !
       cflmax = -999.
       !
       do icell = 1, ncells
          !
          ! compute the sum of mass fluxes leaving the cell
          !
          sumqf = 0.
          !
          do jf = 1, cell(icell)%nof
             !
             ! face identifier
             !
             iface = cell(icell)%face(jf)%atti(FACEID)
             !
             ! get length of current face
             !
             lf = face(iface)%attr(FACELEN)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             ! take into account orientation of the current face
             !
             if ( icell == icelll ) then
                rsgn =  1.
             else if ( icell == icellr ) then
                rsgn = -1.
             endif
             !
             ! compute mass flux at current face
             !
             uf = theta*u1(iface) + (1.-theta)*u0(iface)
             qf = rsgn * lf * hu(iface) * uf
             !
             if ( qf > 0. ) sumqf = sumqf + qf
             !
          enddo
          !
          ! compute the "flow" Courant number
          !
          if ( hs(icell) > epsdry ) then
             !
             area = cell(icell)%attr(CELLAREA)
             !
             cfl = sumqf * dt / hs(icell) / area
             if ( cfl > cflmax ) cflmax = cfl
             !
          endif
          !
       enddo
       !
       ! give warning in case of CFL > 1
       !
       if ( .not. cflmax < 1. ) then
          !
          call msgerr ( 1, 'positivity of the water depth cannot be guaranteed')
          call msgerr ( 0, 'It is advised to reduce the time step' )
          !
       endif
       !
    endif
    !
    ! compute intermediate w-velocity
    !
    if ( ihydro == 1 ) then
       !
       if ( iproj == 1 ) then
          fac = 1.
       else if ( iproj == 2 ) then
          fac = 1. - theta3
       endif
       !
       ! compute cell-based velocity vector
       !
       call perot ( u1, 1, 1 )
       !
       do icell = 1, ncells
          !
          if ( wets(icell) == 1 ) then
             !
             w1bot(icell) = 0.
             !
             ! loop over faces of the cell
             !
             do jf = 1, cell(icell)%nof
                !
                ! face identifier
                !
                iface = cell(icell)%face(jf)%atti(FACEID)
                !
                ! get length of current face
                !
                lf = face(iface)%attr(FACELEN)
                !
                nx = face(iface)%attr(FACENORMX)
                ny = face(iface)%attr(FACENORMY)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                ! take into account orientation of the current face
                !
                if ( icell == icelll ) then
                   rsgn =  1.
                else if ( icell == icellr ) then
                   rsgn = -1.
                endif
                !
                ! compute kinematic condition
                !
                w1bot(icell) = w1bot(icell) - rsgn * lf * dpu(iface) * ( nx * uvc(icell,1) + ny * uvc(icell,2) )
                !
             enddo
             !
             area = cell(icell)%attr(CELLAREA)
             !
             w1bot(icell) = w1bot(icell) / area
             !
             w1top(icell) = w0top(icell) + w0bot(icell) - w1bot(icell) + 2.*dt*fac*q(icell)/hs(icell)
             !
          else
             !
             w1bot(icell) = 0.
             w1top(icell) = 0.
             !
          endif
          !
       enddo
       !
       if ( horwinc ) then
          !
          ! compute horizontal advection term (momentum conservative and explicit)
          ! Note: for the r-ratio the most upwave vertex of upwind cell is used
          !
          propsc = nint(pnums(16))
          kappa  = pnums(17)
          mbound = pnums(18)
          phieby = pnums(19)
          !
          do icell = 1, ncells
             !
             if ( wets(icell) == 1 ) then
                !
                area = cell(icell)%attr(CELLAREA)
                !
                fac = 2. * dt / ( area * hs(icell) )
                !
                ! loop over faces of the cell
                !
                floop: do jf = 1, cell(icell)%nof
                   !
                   ! face identifier
                   !
                   iface = cell(icell)%face(jf)%atti(FACEID)
                   !
                   if ( face(iface)%atti(FMARKER) /= 0 ) cycle floop   ! skip boundary face
                   !
                   ! get vertices of current face
                   !
                   vf1 = face(iface)%atti(FACEV1)
                   vf2 = face(iface)%atti(FACEV2)
                   !
                   ! get length of current face
                   !
                   lf = face(iface)%attr(FACELEN)
                   !
                   ! consider left and right cells of current face
                   !
                   icelll = face(iface)%atti(FACECL)
                   icellr = face(iface)%atti(FACECR)
                   !
                   ! take into account orientation of the current face
                   !
                   if ( icell == icelll ) then
                      rsgn =  1.
                   else if ( icell == icellr ) then
                      rsgn = -1.
                   endif
                   !
                   ! compute mass flux at current face
                   !
                   qf = fac * lf * hu(iface) * u0(iface)
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      w1top(icell) = w1top(icell) - rsgn * qf * ( 0.5 * (w0top(icellr) + w0top(icelll)) - w0top(icell) )
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      if ( qf > 0. ) then
                         !
                         w1top(icell) = w1top(icell) - rsgn * qf * ( w0top(icelll) - w0top(icell) )
                         !
                      else
                         !
                         w1top(icell) = w1top(icell) - rsgn * qf * ( w0top(icellr) - w0top(icell) )
                         !
                      endif
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         rdx = face(iface)%attr(FACEDISTC)
                         !
                         lwfac = 1. - rdx * dt * abs(u0(iface))
                         !
                         if ( qf > 0. ) then
                            !
                            ! get vertices of upwind cell
                            !
                            vc(1) = cell(icelll)%atti(CELLV1)
                            vc(2) = cell(icelll)%atti(CELLV2)
                            vc(3) = cell(icelll)%atti(CELLV3)
                            !
                            ! search for most upwave vertex
                            !
                            do j = 1, 3
                               if ( vc(j) /= vf1 .and. vc(j) /= vf2 ) then
                                  vu = vc(j)
                                  exit
                               endif
                            enddo
                            !
                            ! compute area-weighted averaged concentration at upwave vertex
                            !
                            w0u     = 0.
                            totarea = 0.
                            !
                            do jc = 1, vert(vu)%noc
                               !
                               icella = vert(vu)%cell(jc)%atti(CELLID)
                               !
                               area = cell(icella)%attr(CELLAREA)
                               !
                               w0u = w0u + area * w0top(icella)
                               !
                               totarea = totarea + area
                               !
                            enddo
                            !
                            w0u = w0u / totarea
                            !
                            ! compute solution gradients
                            !
                            grad1 = w0top(icellr) - w0top(icelll)
                            grad2 = w0top(icelll) - w0u
                            !
                            w1top(icell) = w1top(icell) - 0.5 * max(0.,lwfac) * rsgn * qf * fluxlim(grad1,grad2)
                            !
                         else
                            !
                            ! get vertices of upwind cell
                            !
                            vc(1) = cell(icellr)%atti(CELLV1)
                            vc(2) = cell(icellr)%atti(CELLV2)
                            vc(3) = cell(icellr)%atti(CELLV3)
                            !
                            ! search for most upwave vertex
                            !
                            do j = 1, 3
                               if ( vc(j) /= vf1 .and. vc(j) /= vf2 ) then
                                  vu = vc(j)
                                  exit
                               endif
                            enddo
                            !
                            ! compute area-weighted averaged concentration at upwave vertex
                            !
                            w0u     = 0.
                            totarea = 0.
                            !
                            do jc = 1, vert(vu)%noc
                               !
                               icella = vert(vu)%cell(jc)%atti(CELLID)
                               !
                               area = cell(icella)%attr(CELLAREA)
                               !
                               w0u = w0u + area * w0top(icella)
                               !
                               totarea = totarea + area
                               !
                            enddo
                            !
                            w0u = w0u / totarea
                            !
                            ! compute solution gradients
                            !
                            grad1 = w0top(icelll) - w0top(icellr)
                            grad2 = w0top(icellr) - w0u
                            !
                            w1top(icell) = w1top(icell) - 0.5 * max(0.,lwfac) * rsgn * qf * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                      endif
                      !
                   endif
                   !
                enddo floop
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
    ! compute the non-hydrostatic pressure correction
    !
    if ( ihydro == 1 ) then
       !
       ! build the Poisson equation
       !
       do icell = 1, ncells
          !
          if ( wets(icell) == 1 ) then
             !
             ! get coordinates of the cell circumcenter and cell area
             !
             xc = cell(icell)%attr(CELLCCX)
             yc = cell(icell)%attr(CELLCCY)
             !
             area = cell(icell)%attr(CELLAREA)
             !
             ! compute bottom gradient
             !
             dpx = 0.
             dpy = 0.
             !
             do jf = 1, cell(icell)%nof
                !
                ! face identifier
                !
                iface = cell(icell)%face(jf)%atti(FACEID)
                !
                ! get length of the face
                !
                lf = face(iface)%attr(FACELEN)
                !
                nx = face(iface)%attr(FACENORMX)
                ny = face(iface)%attr(FACENORMY)
                !
                ! consider left and right cells of the face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                ! take into account orientation of the face
                !
                if ( icell == icelll ) then
                   rsgn =  1.
                else if ( icell == icellr ) then
                   rsgn = -1.
                endif
                !
                dpx = dpx + rsgn * lf * dpu(iface) * nx
                dpy = dpy + rsgn * lf * dpu(iface) * ny
                !
             enddo
             !
             dpx = 2. * dpx / hs(icell) / area
             dpy = 2. * dpy / hs(icell) / area
             !
             amat(icell,0) = -2.*area / (hs(icell)*hs(icell))
             rhs (icell  ) = area * ( w1top(icell) + w1bot(icell) ) / hs(icell)
             !
             ! loop over faces of the cell
             !
             do jf = 1, cell(icell)%nof
                !
                ! face identifier
                !
                iface = cell(icell)%face(jf)%atti(FACEID)
                !
                lf = face(iface)%attr(FACELEN)
                !
                ! get coordinates of the midface
                !
                xf = face(iface)%attr(FACEMX)
                yf = face(iface)%attr(FACEMY)
                !
                ! consider left and right cells of the face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                fac = lf * ( 1. + dpx * ( xf - xc ) + dpy * ( yf - yc ) )
                !
                ! build matrix and right-hand side
                !
                if ( icell == icelll ) then
                   !
                   if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
                      !
                      amat(icell, 0) = amat(icell,0) + fac * gmat(iface,1) / hum(iface)
                      amat(icell,jf) =                 fac * gmat(iface,2) / hum(iface)
                      !
                   else if ( wlimp(iface) ) then   ! described water level at boundary face
                      !
                      amat(icell,0) = amat(icell,0) + fac * gmat(iface,1) / hum(iface)
                      !
                   endif
                   !
                   rhs(icell) = rhs(icell) + fac * u1(iface)
                   !
                else if ( icell == icellr ) then
                   !
                   if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
                      !
                      amat(icell,0 ) = amat(icell,0) - fac * gmat(iface,2) / hum(iface)
                      amat(icell,jf) =               - fac * gmat(iface,1) / hum(iface)
                      !
                   else if ( wlimp(iface) ) then   ! described water level at boundary face
                      !
                      amat(icell,0) = amat(icell,0) - fac * gmat(iface,1) / hum(iface)
                      !
                   endif
                   !
                   rhs(icell) = rhs(icell) - fac * u1(iface)
                   !
                endif
                !
             enddo
             !
             rhs(icell) = rhs(icell) / (dt*theta3)
             !
          else
             !
             amat(icell,0  ) = 1.
             amat(icell,1:3) = 0.
             rhs (icell    ) = 0.
             !
          endif
          !
       enddo
       !
       ! solve the Poisson equation
       !
       call bicgstabu ( amat, rhs, dq )
       if (STPNOW()) return
       !
    endif
    !
    ! update the non-hydrostatic pressure
    !
    if ( ihydro == 1 ) then
       if ( iproj == 1 ) then
          q = q + dq
       else if ( iproj == 2 ) then
          q = dq
       endif
    endif
    !
    ! correct the flow velocities
    !
    if ( ihydro == 1 ) then
       !
       ! u-velocity
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             u1(iface) = u1(iface) - dt*theta3*( gmat(iface,1) * dq(icelll) + gmat(iface,2) * dq(icellr) ) / hum(iface)
             !
          else if ( wlimp(iface) .and. wetu(iface) == 1 ) then   ! described water level at boundary face
             !
             ! consider boundary cell of current face
             !
             icellb = face(iface)%atti(FACEC1)
             !
             u1(iface) = u1(iface) - dt*theta3* gmat(iface,1) * dq(icellb) / hum(iface)
             !
          endif
          !
       enddo
       !
       ! compute cell-based velocity vector
       !
       call perot ( u1, 1, 1 )
       !
       ! w-velocity
       !
       do icell = 1, ncells
          !
          if ( wets(icell) == 1 ) then
             !
             w1top(icell) = w1top(icell) + w1bot(icell) + 2.*dt*theta3*dq(icell)/hs(icell)
             !
             w1bot(icell) = 0.
             !
             ! loop over faces of the cell
             !
             do jf = 1, cell(icell)%nof
                !
                ! face identifier
                !
                iface = cell(icell)%face(jf)%atti(FACEID)
                !
                ! get length of current face
                !
                lf = face(iface)%attr(FACELEN)
                !
                nx = face(iface)%attr(FACENORMX)
                ny = face(iface)%attr(FACENORMY)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                ! take into account orientation of the current face
                !
                if ( icell == icelll ) then
                   rsgn =  1.
                else if ( icell == icellr ) then
                   rsgn = -1.
                endif
                !
                ! compute kinematic condition
                !
                w1bot(icell) = w1bot(icell) - rsgn * lf * dpu(iface) * ( nx * uvc(icell,1) + ny * uvc(icell,2) )
                !
             enddo
             !
             area = cell(icell)%attr(CELLAREA)
             !
             w1bot(icell) = w1bot(icell) / area
             !
             w1top(icell) = w1top(icell) - w1bot(icell)
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute the mass flux at faces
    !
    do iface = 1, nfaces
       !
       lf = face(iface)%attr(FACELEN)
       !
       uf = theta*u1(iface) + (1.-theta)*u0(iface)
       !
       qn(iface) = lf * hu(iface) * uf
       !
    enddo
    !
    ! calculate net mass outflow based on local continuity equation
    !
    if ( ITEST >= 30 .and. ihydro == 1 ) then
       !
       moutf = 0.
       !
       do icell = 1, ncells
          !
          ! loop over faces of the cell
          !
          do jf = 1, cell(icell)%nof
             !
             ! face identifier
             !
             iface = cell(icell)%face(jf)%atti(FACEID)
             !
             ! get length of current face
             !
             lf = face(iface)%attr(FACELEN)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             ! take into account orientation of the current face
             !
             if ( icell == icelll ) then
                rsgn =  1.
             else if ( icell == icellr ) then
                rsgn = -1.
             endif
             !
             moutf = moutf + rsgn * lf * u1(iface)
             !
          enddo
          !
          area = cell(icell)%attr(CELLAREA)
          !
          moutf = moutf + area * ( w1top(icell) - w1bot(icell) ) / hs(icell)
          !
       enddo
       !
       write(PRINTF,101) moutf
       !
    endif
    !
    ! calculate total displaced volume of water (reference is at z = 0) and energy
    ! Note: this testing is appropriate for a closed domain with reflective boundaries
    !
    if ( ITEST >= 30 ) then
       !
       vol  = 0.
       ener = 0.
       !
       do icell = 1, ncells
          !
          ! compute displaced volume in the cell center
          vol = vol + s1(icell)
          !
          area = cell(icell)%attr(CELLAREA)
          !
          ! compute potential energy in the cell
          ener = ener + 0.5 * grav * area * s1(icell) * s1(icell)
          !
       enddo
       !
       do iface = 1, nfaces
          !
          lf  = face(iface)%attr(FACELEN)
          rdx = face(iface)%attr(FACEDISTC)
          !
          ! compute kinetic energy in the control volume of the face
          ! note: the total erea of the control volumes is twice the domain area
          !       which compensates the amount of velocity components (only half)
          ener = ener + 0.5 * lf * hu(iface) * u1(iface) * u1(iface) / rdx
          !
       enddo
       !
       write(PRINTF,102) vol
       write(PRINTF,103) ener
       !
    endif
    !
    ! apply wave absorption by means of sponge layer, if appropriate
    !
    do j = 1, nspl
       !
       do iface = 1, nfaces
          !
          gamma = sponu(iface,j)%gamma
          bface = sponu(iface,j)%bface
          !
          if ( wetu(iface) == 1 ) then
             !
             u1(iface) = (1.-gamma) * u1(iface) + gamma * u1(bface)
             !
          endif
          !
       enddo
       !
       do icell = 1, ncells
          !
          gamma = spons(icell,j)%gamma
          bcell = spons(icell,j)%bcell
          !
          if ( wets(icell) == 1 ) then
             !
             s1(icell) = (1.-gamma) * s1(icell) + gamma * s1(bcell)
             !
          endif
          !
       enddo
       !
    enddo
    !
    ! compute the maximum CFL number
    !
    cflmax = -999.
    !
    do iface = 1, nfaces
       !
       if ( wetu(iface) == 1 ) then                    ! only consider wet face
          !
          rdx = face(iface)%attr(FACEDISTG)
          !
          if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
             !
             cfl = rdx * dt * abs(u1(iface))
             if ( cfl > cflmax ) cflmax = cfl
             !
          endif
          !
       endif
       !
    enddo
    !
    ! give warning in case of CFL > 1
    !
    if ( .not. cflmax < 1. ) then
       !
       call msgerr ( 1, 'CFL condition is violated!')
       call msgerr ( 0, 'It is advised to reduce the time step!' )
       !
    endif
    !
 101 format (2x,'the net mass outflow is ',e14.8e2)
 102 format (2x,'the total displaced volume of water is ',e14.8e2)
 103 format (2x,'the total energy is ',e14.8e2)
    !
end subroutine SwashImpDepUflow
