subroutine SwashImpLayUflow
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
!    1.00,   April 2020: New subroutine
!   10.01, January 2023: add higher order advection
!   10.01,   April 2023: Newton-type iteration added
!
!   Purpose
!
!   Performs the time integration for the non-hydrostatic, layer-averaged shallow water equations on triangular mesh
!
!   Method
!
!   The time integration with respect to the continuity equation and the water level gradient of the
!   u-momentum equation is based on a theta-scheme. Only a value of 0.5 <= theta <= 1 will be taken.
!
!   The time integration with respect to the advective term is based on Euler explicit,
!   while that for the bottom friction is based on Euler implicit and for the
!   non-hydrostatic pressure gradient a semi-implicit approach is employed (theta-scheme). Both vertical
!   advective and viscosity terms are treated semi-implicit as well. This results in a tri-diagonal system.
!
!   The space discretization of the advective term is momentum conservative and is approximated by
!   first order upwind or higher order (flux-limited) scheme (CDS, Fromm, BDF, QUICK, MUSCL, Koren, etc.).
!   The r-ratio formulation based on most upwave vertex is employed.
!
!   The space discretization of the vertical advective and viscosity terms is based on central differences
!   in a finite volume fashion.
!
!   The vertical grid schematization gives rise to the definition of the vertical velocity with respect to
!   the moving layer interfaces. This relative velocity, stored in array wom, is defined as the difference
!   between the vertical velocity along the streamline and the vertical velocity along the interface.
!   However, in case of hydrostatic flows, the relative vertical velocity will be derived from the
!   layer-averaged continuity equation.
!
!   The vertical velocity in z-direction, stored in array w1, is obtained from the solution of the
!   w-momentum equation which contains the z-gradient of the non-hydrostatic pressure that is discretized
!   by means of either (explicit) central differences or the (implicit) Keller-box scheme. Optionally,
!   horizontal terms are treated explicit, while vertical terms are treated semi-implicit. This results in
!   a tri-diagonal system. The space discretization of the vertical advective and viscosity terms is based
!   on higher order (flux-limited) schemes and central differences, respectively, in a finite volume fashion.
!
!   The non-hydrostatic pressure is obtained by means of the second order accurate pressure correction technique.
!
!   Modules used
!
    use ocpcomm4
    use m_genarr, only: hlay, wlimp, work
    use SwashCommdata3
    use SwashTimecomm
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Parameter variables
!
    real   , parameter :: epswom  = 0.005 ! tolerance for relative vertical velocity at surface
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
    integer                                :: k        ! loop counter over vertical layers
    integer                                :: kd       ! index of layer k-1
    integer                                :: kdd      ! index of layer k-2
    integer                                :: ku       ! index of layer k+1
    integer                                :: kuu      ! index of layer k+2
    integer                                :: l        ! loop counter
    integer, dimension(3)                  :: vc       ! vertices of present cell
    integer                                :: vf1      ! first vertex of present face
    integer                                :: vf2      ! second vertex of present face
    integer                                :: vu       ! upwind vertex
    !
    real                                   :: area     ! area of present cell
    real                                   :: bi       ! inverse of main diagonal of the matrix
    real                                   :: cadv     ! contribution to advection term
    real                                   :: cfl      ! CFL number
    real                                   :: ctrkb    ! contribution of vertical terms below considered point
    real                                   :: ctrkt    ! contribution of vertical terms above considered point
    real                                   :: ctrs     ! contribution to system of water level equations
    real                                   :: cvisc    ! contribution to viscosity term
    real                                   :: dbs      ! change in water level at boundary face over time
    real                                   :: ener     ! total energy of closed system
    real                                   :: epsab2   ! small amount to modify weights of AB2 scheme
    real                                   :: fac      ! a factor
    real                                   :: fac1     ! another factor
    real                                   :: fac2     ! some other factor
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
    real                                   :: kwd      ! =1. if layer k-1 exists otherwise 0.
    real                                   :: kwu      ! =1. if layer k+1 exists otherwise 0.
    real                                   :: lf       ! length of face
    real                                   :: lwfac    ! Lax-Wendroff factor
    real                                   :: moutf    ! net mass outflow
    real                                   :: nx       ! x-component of normal to face
    real                                   :: ny       ! y-component of normal to face
    real                                   :: qf       ! (total) mass flux
    real                                   :: rdx      ! reciprocal of distance between circumcenters adjacent to face
    real                                   :: rhou     ! density of water in velocity point
    real                                   :: rsgn     ! sign for indicating face orientation
    real                                   :: theta    ! implicitness factor for time integration of continuity equation
    real                                   :: theta2   ! implicitness factor for water level gradient
    real                                   :: theta3   ! implicitness factor for non-hydrostatic pressure gradient
    real                                   :: thetau   ! implicitness factor for vertical terms in u-momentum equation
    real                                   :: thetaw   ! implicitness factor for vertical terms in w-momentum equation
    real                                   :: totarea  ! total area of all cells around vertex
    real                                   :: u        ! layer-averaged u-velocity component
    real                                   :: uf       ! updated depth-integrated flow velocity at present face
    real                                   :: ut       ! tangential velocity component
    real                                   :: utot     ! velocity magnitude
    real                                   :: ux       ! x-component of flow velocity in cell circumenter
    real                                   :: uy       ! y-component of flow velocity in cell circumenter
    real                                   :: v        ! layer-averaged v-velocity component
    real                                   :: vol      ! total displaced volume of water
    real                                   :: w0u      ! averaged w-velocity in upwind vertex
    real                                   :: xc       ! x-coordinate of cell circumcenter
    real                                   :: xf       ! x-coordinate of face center
    real                                   :: yc       ! y-coordinate of cell circumcenter
    real                                   :: yf       ! y-coordinate of face center
    real                                   :: zkxd     ! x-component of gradient in zku at (k-1)-th interface
    real                                   :: zkxu     ! x-component of gradient in zku at k-th interface
    real                                   :: zkyd     ! y-component of gradient in zku at (k-1)-th interface
    real                                   :: zkyu     ! y-component of gradient in zku at k-th interface
    !
    logical                                :: STPNOW   ! indicates that program must stop
    !
    character(85)                          :: msgstr   ! string to pass message
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
    if (ltrace) call strace (ient,'SwashImpLayUflow')
    !
    ! point to vertex, cell and face objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    theta  = pnums( 1)
    theta2 = pnums( 4)
    theta3 = pnums( 5)
    thetau = pnums(31)
    thetaw = pnums(32)
    !
    ! build the u-momentum equation
    !
    ! initialize system of equations on dry faces
    !
    do iface = 1, nfaces
       !
       if ( face(iface)%atti(FMARKER) == 0 .or. wlimp(iface) ) then   !  face is either internal or boundary face with described water level
          !
          if ( wetu(iface) /= 1 ) then
             !
             amatu(iface,:,1) = 1.
             amatu(iface,:,2) = 0.
             amatu(iface,:,3) = 0.
             rhsu (iface,:  ) = 0.
             !
          endif
          !
       endif
       !
    enddo
    !
    ! compute the time derivative
    !
    do iface = 1, nfaces
       !
       if ( face(iface)%atti(FMARKER) == 0 .or. wlimp(iface) ) then   !  face is either internal or boundary face with described water level
          !
          if ( wetu(iface) == 1 ) then
             !
             fac = dt * thetau
             !
             do k = 1, kmax
                amatu(iface,k,1) = 1. / fac
                rhsu (iface,k  ) = u0(iface,k) / fac
             enddo
             !
          endif
          !
       endif
       !
    enddo
    !
    ! compute the mass flux
    !
    do icell = 1, ncells
       !
       do k = 1, kmax
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
             quf(icell,k,jf) = rsgn * lf * hku(iface,k) * u0(iface,k)
             !
          enddo
          !
       enddo
       !
    enddo
    !
    ! compute cell-based velocity vector
    !
    if ( irough /= 0 .and. irough /= 4 .and. irough /= 11 ) then
       call perot ( udep, 1, 1 )
       work(:,1) = uvc(:,1,1)
       work(:,2) = uvc(:,1,2)
    endif
    call perot ( u0, 1, kmax )
    !
    ! compute horizontal advection term using first order upwind (momentum conservative)
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
             area = cell(icell)%attr(CELLAREA)
             !
             do k = 1, kmax
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
                      if ( u0(ifacel,k) > 0. ) then
                         icellu = icelll
                      else
                         icellu = icellr
                      endif
                      !
                      ! compute contribution to advection term - only update for ingoing flux
                      !
                      if ( quf(icell,k,jf) < 0. ) cadv = cadv + quf(icell,k,jf) * ( nx * uvc(icellu,k,1) + ny * uvc(icellu,k,2) - u0(iface,k) )
                      !
                   endif
                   !
                enddo
                !
                rhsu(iface,k) = rhsu(iface,k) - faca * cadv / area / hkumn(iface,k)
                !
             enddo
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
          area = cell(icellb)%attr(CELLAREA)
          !
          ! get orientation at boundary face
          !
          if ( icellb == icelll ) then
             rsgn =  1.
          else if ( icellb == icellr ) then
             rsgn = -1.
          endif
          !
          do k = 1, kmax
             !
             if ( rsgn * u0(iface,k) > 0. ) then
                !
                ! outflow
                !
                cadv = 0.
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
                      if ( u0(ifacel,k) > 0. ) then
                         icellu = face(ifacel)%atti(FACECL)
                      else
                         icellu = face(ifacel)%atti(FACECR)
                      endif
                      !
                      ! update for ingoing flux
                      !
                      if ( quf(icellb,k,jf) < 0. ) cadv = cadv + quf(icellb,k,jf) * ( nx * uvc(icellu,k,1) + ny * uvc(icellu,k,2) - u0(iface,k) )
                      !
                   endif
                   !
                enddo
                !
                rhsu(iface,k) = rhsu(iface,k) - cadv / area / hkumn(iface,k)
                !
             endif
             !
          enddo
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
                do k = 1, kmax
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
                         lwfac = 1. - rdx * dt * abs(u0(ifacel,k))
                         !
                         ! get vertices of face under consideration
                         !
                         vf1 = face(ifacel)%atti(FACEV1)
                         vf2 = face(ifacel)%atti(FACEV2)
                         !
                         ! consider up- and downwind cells of local face
                         !
                         if ( u0(ifacel,k) > 0. ) then
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
                            ux = ux + area * uvc(icella,k,1)
                            uy = uy + area * uvc(icella,k,2)
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
                         grad1x = uvc(icelld,k,1) - uvc(icellu,k,1)
                         grad2x = uvc(icellu,k,1) - ux
                         !
                         grad1y = uvc(icelld,k,2) - uvc(icellu,k,2)
                         grad2y = uvc(icellu,k,2) - uy
                         !
                         ! compute high order limited part
                         !
                         hox = 0.5 * max(0.,lwfac) * fluxlim(grad1x,grad2x)
                         hoy = 0.5 * max(0.,lwfac) * fluxlim(grad1y,grad2y)
                         !
                         cadv = cadv + quf(icell,k,jf) * ( nx * hox + ny * hoy )
                         !
                      endif
                      !
                   enddo
                   !
                   area = cell(icell)%attr(CELLAREA)
                   !
                   rhsu(iface,k) = rhsu(iface,k) - faca * cadv / area / hkumn(iface,k)
                   !
                enddo
                !
             enddo
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute water level gradient
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
             do k = 1, kmax
                !
                ! compute arithmetic average layer thickness at present face
                !
                hf = 0.5 * ( hks(icelll,k) + hks(icellr,k) )
                !
                ! store water level gradient
                !
                rhsu(iface,k) = rhsu(iface,k) - grav * rdx * hf * ( s0(icellr) - s0(icelll) ) / hkum(iface,k)
                !
             enddo
             !
          endif
          !
       endif
       !
    enddo
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
             do k = 1, kmax
                !
                ! store water level gradient
                !
                rhsu(iface,k) = rhsu(iface,k) - pfac(iface) * grav * rdx * hf * hlay(k) * ( bcso(iface) - s0(icellb) ) / hkum(iface,k)
                !
             enddo
             !
          endif
          !
       endif
       !
    enddo
    !
    ! build gradient matrix for non-hydrostatic pressure
    !
    if ( ihydro == 1 ) then
       !
       ! Keller-box scheme, so non-hydrostatic pressure is located at the centers of layer interfaces
       !
       do k = 1, kmax
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
                fac = rdx / hkum(iface,k)
                !
                fac1 = zks(icellr,k-1) - zks(icelll,k-1)
                fac2 = zks(icellr,k  ) - zks(icelll,k  )
                !
                gmatu(iface,k,1) = (-0.5 * hks(icelll,k) -     finp  * fac1) * fac
                gmatu(iface,k,2) = (-0.5 * hks(icelll,k) +     finp  * fac2) * fac
                gmatu(iface,k,3) = ( 0.5 * hks(icellr,k) - (1.-finp) * fac1) * fac
                gmatu(iface,k,4) = ( 0.5 * hks(icellr,k) + (1.-finp) * fac2) * fac
                !
             else if ( wlimp(iface) .and. wetu(iface) == 1 ) then   ! described water level at boundary face
                !
                ! consider boundary cell of current face
                !
                icellb = face(iface)%atti(FACEC1)
                !
                fac = -rdx * hks(icellb,k) / hkum(iface,k)
                !
                gmatu(iface,k,1) = fac
                gmatu(iface,k,2) = fac
                gmatu(iface,k,3) = 0.
                gmatu(iface,k,4) = 0.
                !
             else
                !
                gmatu(iface,k,1) = 0.
                gmatu(iface,k,2) = 0.
                gmatu(iface,k,3) = 0.
                gmatu(iface,k,4) = 0.
                !
             endif
             !
          enddo
          !
       enddo
       !
       do iface = 1, nfaces
          !
          gmatu(iface,1,1) = 0.
          gmatu(iface,1,3) = 0.
          !
       enddo
       !
       ! to reduce the pressure Poisson equation set pressure of bottom face to that of top face for a number of layers
       !
       do l = 1, qlay
          !
          do iface = 1, nfaces
             !
             gmatu(iface,kmax-l+1,1) = gmatu(iface,kmax-l+1,1) + gmatu(iface,kmax-l+1,2)
             gmatu(iface,kmax-l+1,2) = 0.
             gmatu(iface,kmax-l+1,3) = gmatu(iface,kmax-l+1,3) + gmatu(iface,kmax-l+1,4)
             gmatu(iface,kmax-l+1,4) = 0.
             !
          enddo
          !
       enddo
       !
    else if ( ihydro == 2 ) then
       !
       ! central differences, so non-hydrostatic pressure is located at the cell centers
       !
       do k = 1, kmax
          !
          kd = max(k-1,1   )
          ku = min(k+1,kmax)
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
                fac1 = rdx * ( zks(icellr,k-1) - zks(icelll,k-1) ) / ( hkum(iface,kd) + hkum(iface,k ) )
                fac2 = rdx * ( zks(icellr,k  ) - zks(icelll,k  ) ) / ( hkum(iface,k ) + hkum(iface,ku) )
                !
                gmatu(iface,k,1) = -    finp  * fac1
                gmatu(iface,k,3) = -(1.-finp) * fac1
                gmatu(iface,k,2) = ( -rdx*hks(icelll,k) -     finp  * hkum(iface,kd)*fac1 +     finp  * hkum(iface,ku)*fac2 ) / hkum(iface,k)
                gmatu(iface,k,4) = (  rdx*hks(icellr,k) - (1.-finp) * hkum(iface,kd)*fac1 + (1.-finp) * hkum(iface,ku)*fac2 ) / hkum(iface,k)
                gmatu(iface,k,5) =      finp  * fac2
                gmatu(iface,k,6) =  (1.-finp) * fac2
                !
             else if ( wlimp(iface) .and. wetu(iface) == 1 ) then   ! described water level at boundary face
                !
                ! consider boundary cell of current face
                !
                icellb = face(iface)%atti(FACEC1)
                !
                fac = -2. * rdx * hks(icellb,k) / hkum(iface,k)
                !
                gmatu(iface,k,1) = 0.
                gmatu(iface,k,2) = fac
                gmatu(iface,k,3) = 0.
                gmatu(iface,k,4) = 0.
                gmatu(iface,k,5) = 0.
                gmatu(iface,k,6) = 0.
                !
             else
                !
                gmatu(iface,k,1) = 0.
                gmatu(iface,k,2) = 0.
                gmatu(iface,k,3) = 0.
                gmatu(iface,k,4) = 0.
                gmatu(iface,k,5) = 0.
                gmatu(iface,k,6) = 0.
                !
             endif
             !
          enddo
          !
       enddo
       !
       do iface = 1, nfaces
          !
          gmatu(iface,kmax,2) = gmatu(iface,kmax,2) + 2.*gmatu(iface,kmax,5)
          gmatu(iface,kmax,1) = gmatu(iface,kmax,1) -    gmatu(iface,kmax,5)
          gmatu(iface,kmax,4) = gmatu(iface,kmax,4) + 2.*gmatu(iface,kmax,6)
          gmatu(iface,kmax,3) = gmatu(iface,kmax,3) -    gmatu(iface,kmax,6)
          gmatu(iface,kmax,5) = 0.
          gmatu(iface,kmax,6) = 0.
          !
          gmatu(iface,1,2) = gmatu(iface,1,2) - gmatu(iface,1,1)
          gmatu(iface,1,4) = gmatu(iface,1,4) - gmatu(iface,1,3)
          gmatu(iface,1,1) = 0.
          gmatu(iface,1,3) = 0.
          !
       enddo
       !
    endif
    !
    ! compute gradient of non-hydrostatic pressure
    !
    if ( ihydro /= 0 ) then
       !
       if ( iproj == 1 ) then
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
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
                   rhsu(iface,k) = rhsu(iface,k) - gmatu(iface,k,1)*q(icelll,kd) - gmatu(iface,k,2)*q(icelll,k ) - gmatu(iface,k,3)*q(icellr,kd)  &
                                                 - gmatu(iface,k,4)*q(icellr,k ) - gmatu(iface,k,5)*q(icelll,ku) - gmatu(iface,k,6)*q(icellr,ku)
                   !
                else if ( wlimp(iface) .and. wetu(iface) == 1 ) then   ! described water level at boundary face
                   !
                   ! consider boundary cell of current face
                   !
                   icellb = face(iface)%atti(FACEC1)
                   !
                   rhsu(iface,k) = rhsu(iface,k) - gmatu(iface,k,1)*q(icellb,kd) - gmatu(iface,k,2)*q(icellb,k) - gmatu(iface,k,5)*q(icellb,ku)
                   !
                endif
                !
             enddo
             !
          enddo
          !
       else if ( iproj == 2 .and. theta3 /= 1. ) then
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
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
                   rhsu(iface,k) = rhsu(iface,k) - (1.-theta3) * ( gmatu(iface,k,1)*q(icelll,kd) + gmatu(iface,k,2)*q(icelll,k ) + gmatu(iface,k,3)*q(icellr,kd)  &
                                                                 + gmatu(iface,k,4)*q(icellr,k ) + gmatu(iface,k,5)*q(icelll,ku) + gmatu(iface,k,6)*q(icellr,ku) )
                   !
                else if ( wlimp(iface) .and. wetu(iface) == 1 ) then   ! described water level at boundary face
                   !
                   ! consider boundary cell of current face
                   !
                   icellb = face(iface)%atti(FACEC1)
                   !
                   rhsu(iface,k) = rhsu(iface,k) - (1.-theta3) * ( gmatu(iface,k,1)*q(icellb,kd) + gmatu(iface,k,2)*q(icellb,k) + gmatu(iface,k,5)*q(icellb,ku) )
                   !
                endif
                !
             enddo
             !
          enddo
          !
       endif
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
                do k = 1, kmax
                   !
                   rhsu(iface,k) = rhsu(iface,k) - 0.5 * grav * rdx * hkum(iface,k) * (rho(icellr,k) - rho(icelll,k)) / rhow
                   !
                enddo
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
                finp = face(iface)%attr(FACELINPF)
                !
                rdx = face(iface)%attr(FACEDISTC)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                fac1 = 0.
                fac2 = 0.
                !
                do k = 2, kmax
                   !
                   rhou = rhow + finp * rho(icelll,k-1) + (1.-finp) * rho(icellr,k-1)
                   !
                   fac1 = fac1 + rhou * (hks(icellr,k-1) - hks(icelll,k-1)) + hkum(iface,k-1) * (rho(icellr,k-1) - rho(icelll,k-1))
                   !
                   fac2 = fac2 + hks(icellr,k-1) - hks(icelll,k-1)
                   !
                   rhou = rhow + finp * rho(icelll,k) + (1.-finp) * rho(icellr,k)
                   !
                   rhsu(iface,k) = rhsu(iface,k) - grav * rdx * ( fac1 - rhou*fac2 ) / rhow
                   !
                enddo
                !
             endif
             !
          endif
          !
       enddo
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
                rhsu(iface,:) = rhsu(iface,:) - rdx * ( patm(icellr) - patm(icelll) ) / rhow
                !
             endif
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute horizontal viscosity term
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
                area = cell(icell)%attr(CELLAREA)
                !
                do k = 1, kmax
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
                         cvisc = cvisc + rsgn * vnu2d(ifacel) * hkum(ifacel,k) * lf * rdx * ( nx * ( uvc(icellr,k,1)- uvc(icelll,k,1) ) + &
                                                                                              ny * ( uvc(icellr,k,2)- uvc(icelll,k,2) ) )
                         !
                      endif
                      !
                   enddo
                   !
                   rhsu(iface,k) = rhsu(iface,k) + faca * cvisc / area / hkum(iface,k)
                   !
                enddo
                !
             enddo
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute explicit part of wind stress term (top layer only), if appropriate
    !
    if ( iwind /= 0 ) then
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             rhsu(iface,1) = rhsu(iface,1) + windu(iface) / hkum(iface,1)
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute implicit part of wind stress term (top layer only), if appropriate
    !
    if ( relwnd ) then
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             amatu(iface,1,1) = amatu(iface,1,1) + cwndu(iface) / hkum(iface,1)
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute bottom friction term (bottom layer only), if appropriate
    !
    if ( irough == 4 ) then
       !
       ! logarithmic wall-law
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .or. wlimp(iface) ) then   !  face is either internal or boundary face with described water level
             !
             if ( wetu(iface) == 1 ) then
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
                   fac = finp * logfrc(icelll,1) + (1.-finp) * logfrc(icellr,1)
                   !
                else if ( icelll == 0 ) then
                   !
                   fac = logfrc(icellr,1)
                   !
                else
                   !
                   fac = logfrc(icelll,1)
                   !
                endif
                !
                amatu(iface,kmax,1) = amatu(iface,kmax,1) + fac / hkum(iface,kmax)
                !
             endif
             !
          endif
          !
       enddo
       !
    else if ( irough == 11 ) then
       !
       ! linear bottom friction
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .or. wlimp(iface) ) then   !  face is either internal or boundary face with described water level
             !
             if ( wetu(iface) == 1 ) then
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
                   fac = finp * cfricu(icelll) + (1.-finp) * cfricu(icellr)
                   !
                else if ( icelll == 0 ) then
                   !
                   fac = cfricu(icellr)
                   !
                else
                   !
                   fac = cfricu(icelll)
                   !
                endif
                !
                amatu(iface,kmax,1) = amatu(iface,kmax,1) + fac / hkum(iface,kmax)
                !
             endif
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
          if ( face(iface)%atti(FMARKER) == 0 .or. wlimp(iface) ) then   !  face is either internal or boundary face with described water level
             !
             if ( wetu(iface) == 1 ) then
                !
                finp = face(iface)%attr(FACELINPF)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                do l = 1, 2
                   !
                   if ( icelll /= 0 .and. icellr /= 0 ) then
                      !
                      if ( l == 1 ) then
                         !
                         ! left cell of current face
                         !
                         icell = icelll
                         faca  = finp
                         !
                      else
                         !
                         ! right cell of current face
                         !
                         icell = icellr
                         faca  = 1. - finp
                         !
                      endif
                      !
                   else if ( icelll == 0 ) then
                      !
                      ! right cell of current face
                      !
                      icell = icellr
                      !
                      if ( l == 1 ) then
                         !
                         faca  = 0.
                         !
                      else
                         faca  = 1.
                         !
                      endif
                      !
                   else
                      !
                      ! left cell of current face
                      !
                      icell = icelll
                      !
                      if ( l == 1 ) then
                         !
                         faca  = 1.
                         !
                      else
                         !
                         faca  = 0.
                         !
                      endif
                      !
                   endif
                   !
                   ! compute velocity magnitude
                   !
                   utot = sqrt( uvc(icell,kmax,1)*uvc(icell,kmax,1) + uvc(icell,kmax,2)*uvc(icell,kmax,2) )
                   !
                   if ( utot > 1.e-8 ) then
                      !
                      fac = work(icell,1)*work(icell,1) + work(icell,2)*work(icell,2)
                      !
                      amatu(iface,kmax,1) = amatu(iface,kmax,1) + faca * cfricu(icell) * fac / ( utot * hkum(iface,kmax) )
                      !
                   endif
                   !
                enddo
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
             rhsu(iface,:) = rhsu(iface,:) + (0.5+epsab2) * cfu(iface,:) / hkum(iface,:)
             !
             ! at time level n
             !
             cfu(iface,:) = 0.
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
                do k = 1, kmax
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
                      ux = ux + rsgn * lf * ( xf - xc ) * hku(ifacel,k) * u0(ifacel,k)
                      uy = uy + rsgn * lf * ( yf - yc ) * hku(ifacel,k) * u0(ifacel,k)
                      !
                   enddo
                   !
                   cfu(iface,k) = cfu(iface,k) + faca * fcor(icell,1) * ( ny * ux - nx * uy ) / area
                   !
                enddo
                !
             enddo
             !
             rhsu(iface,:) = rhsu(iface,:) - (1.5+epsab2) * cfu(iface,:) / hkum(iface,:)
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute vertical terms (implicit)
    !
    propsc = nint(pnums(36))
    kappa  = pnums(37)
    mbound = pnums(38)
    phieby = pnums(39)
    !
    do k = 1, kmax
       !
       kd  = max(k-1,1   )
       kdd = max(k-2,1   )
       ku  = min(k+1,kmax)
       kuu = min(k+2,kmax)
       !
       kwd = 1.
       kwu = 1.
       if ( k == 1    ) kwd = 0.
       if ( k == kmax ) kwu = 0.
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .or. wlimp(iface) ) then   !  face is either internal or boundary face with described water level
             !
             if ( wetu(iface) == 1 ) then
                !
                finp = face(iface)%attr(FACELINPF)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                ! advection term
                !
                if ( propsc == 3 .and. kappa == 1. ) then
                   !
                   ! central differences
                   !
                   if ( icelll /= 0 .and. icellr /= 0 ) then
                      !
                      ctrkt = finp * wom(icelll,k-1) + (1.-finp) * wom(icellr,k-1)
                      ctrkb = finp * wom(icelll,k  ) + (1.-finp) * wom(icellr,k  )
                      !
                   else if ( icelll == 0 ) then
                      !
                      ctrkt = wom(icellr,k-1)
                      ctrkb = wom(icellr,k  )
                      !
                   else
                      !
                      ctrkt = wom(icelll,k-1)
                      ctrkb = wom(icelll,k  )
                      !
                   endif
                   !
                   ctrkt = kwd * ctrkt / ( hkum(iface,kd) + hkum(iface,k ) )
                   ctrkb = kwu * ctrkb / ( hkum(iface,k ) + hkum(iface,ku) )
                   !
                   amatu(iface,k,1) = amatu(iface,k,1) - ctrkt + ctrkb
                   amatu(iface,k,2) =  ctrkt
                   amatu(iface,k,3) = -ctrkb
                   !
                else
                   !
                   ! first order upwind scheme
                   !
                   if ( icelll /= 0 .and. icellr /= 0 ) then
                      !
                      ctrkt = finp * wom(icelll,k-1) + (1.-finp) * wom(icellr,k-1)
                      ctrkb = finp * wom(icelll,k  ) + (1.-finp) * wom(icellr,k  )
                      !
                   else if ( icelll == 0 ) then
                      !
                      ctrkt = wom(icellr,k-1)
                      ctrkb = wom(icellr,k  )
                      !
                   else
                      !
                      ctrkt = wom(icelll,k-1)
                      ctrkb = wom(icelll,k  )
                      !
                   endif
                   !
                   ctrkt = kwd * ctrkt / hkum(iface,k)
                   ctrkb = kwu * ctrkb / hkum(iface,k)
                   !
                   amatu(iface,k,1) = amatu(iface,k,1) - ctrkt + ctrkb + max(ctrkt,0.) - min(ctrkb,0.)
                   amatu(iface,k,2) =  min(ctrkt,0.)
                   amatu(iface,k,3) = -max(ctrkb,0.)
                   !
                   ! add higher order (flux-limited) correction, if appropriate
                   !
                   if ( propsc /= 1 ) then
                      !
                      if ( ctrkt > 0. ) then
                         !
                         grad1 = u0(iface,kd) - u0(iface,k )
                         grad2 = u0(iface,k ) - u0(iface,ku)
                         !
                         rhsu(iface,k) = rhsu(iface,k) - 0.5 * ctrkt * fluxlim(grad1,grad2)
                         !
                      else if ( ctrkt < 0. ) then
                         !
                         grad1 = u0(iface,kd ) - u0(iface,k )
                         grad2 = u0(iface,kdd) - u0(iface,kd)
                         !
                         rhsu(iface,k) = rhsu(iface,k) + 0.5 * ctrkt * fluxlim(grad1,grad2)
                         !
                      endif
                      !
                      if ( ctrkb > 0. ) then
                         !
                         grad1 = u0(iface,k ) - u0(iface,ku )
                         grad2 = u0(iface,ku) - u0(iface,kuu)
                         !
                         rhsu(iface,k) = rhsu(iface,k) + 0.5 * ctrkb * fluxlim(grad1,grad2)
                         !
                      else if ( ctrkb < 0. ) then
                         !
                         grad1 = u0(iface,k ) - u0(iface,ku)
                         grad2 = u0(iface,kd) - u0(iface,k )
                         !
                         rhsu(iface,k) = rhsu(iface,k) - 0.5 * ctrkb * fluxlim(grad1,grad2)
                         !
                      endif
                      !
                   endif
                   !
                endif
                !
                ! viscosity term
                !
                if ( iturb < 2 ) then
                   !
                   if ( icelll /= 0 .and. icellr /= 0 ) then
                      !
                      ctrkt = finp * vnu3d(icelll,k-1) + (1.-finp) * vnu3d(icellr,k-1)
                      ctrkb = finp * vnu3d(icelll,k  ) + (1.-finp) * vnu3d(icellr,k  )
                      !
                   else if ( icelll == 0 ) then
                      !
                      ctrkt = vnu3d(icellr,k-1)
                      ctrkb = vnu3d(icellr,k  )
                      !
                   else
                      !
                      ctrkt = vnu3d(icelll,k-1)
                      ctrkb = vnu3d(icelll,k  )
                      !
                   endif
                   !
                   ctrkt = 2.* ctrkt / ( hkum(iface,k)*( hkum(iface,kd) + hkum(iface,k ) ) )
                   ctrkb = 2.* ctrkb / ( hkum(iface,k)*( hkum(iface,k ) + hkum(iface,ku) ) )
                   !
                   amatu(iface,k,1) = amatu(iface,k,1) + kwd*ctrkt + (2.*kwu-1.)*ctrkb
                   amatu(iface,k,2) = amatu(iface,k,2) - kwd*ctrkt +    (1.-kwu)*ctrkb
                   amatu(iface,k,3) = amatu(iface,k,3) -                     kwu*ctrkb
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
    ! solve the u-momentum equation
    !
    do iface = 1, nfaces
       !
       if ( face(iface)%atti(FMARKER) == 0 .or. wlimp(iface) ) then   !  face is either internal or boundary face with described water level
          !
          bi = 1./amatu(iface,1,1)
          !
          amatu(iface,1,1) = bi
          amatu(iface,1,3) = amatu(iface,1,3)*bi
          rhsu (iface,1  ) = rhsu (iface,1  )*bi
          !
          do k = 2, kmax
             !
             bi = 1./(amatu(iface,k,1) - amatu(iface,k,2)*amatu(iface,k-1,3))
             amatu(iface,k,1) = bi
             amatu(iface,k,3) = amatu(iface,k,3)*bi
             rhsu (iface,k  ) = (rhsu(iface,k) - amatu(iface,k,2)*rhsu(iface,k-1))*bi
             !
          enddo
          !
          u1(iface,kmax) = rhsu(iface,kmax)
          do k = kmax-1, 1, -1
             u1(iface,k) = rhsu(iface,k) - amatu(iface,k,3)*u1(iface,k+1)
          enddo
          !
       endif
       !
    enddo
    !
    ! re-update the solution in case of thetau <> 1
    !
    if ( thetau /= 1. ) then
       !
       do k = 1, kmax
          !
          do iface = 1, nfaces
             !
             if ( face(iface)%atti(FMARKER) == 0 .or. wlimp(iface) ) then   !  face is either internal or boundary face with described water level
                !
                if ( wetu(iface) == 1 ) then
                   !
                   u1(iface,k) = ( u1(iface,k) - (1.-thetau) * u0(iface,k) ) / thetau
                   !
                endif
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
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
             uf = 0.
             !
             do k = 1, kmax
                !
                uf = uf + hku(iface,k) * ( theta*u1(iface,k) + (1.-theta)*u0(iface,k) )
                !
             enddo
             !
             rhs(icell) = rhs(icell) - dt * rsgn * lf * uf
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
             uf = 0.
             !
             do k = 1, kmax
                !
                uf = uf + hku(iface,k) * ( theta*u1(iface,k) + (1.-theta)*u0(iface,k) )
                !
             enddo
             !
             rhsn(icell) = rhsn(icell) - dt * rsgn * lf * uf
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
             do k = 1, kmax
                !
                u1(iface,k) = u1(iface,k) - pfac(iface) * grav * theta2 * dt * rdx * hf * ( dbs - ds(icellb) ) / hum(iface)
                !
             enddo
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
             do k = 1, kmax
                !
                u1(iface,k) = u1(iface,k) - grav * theta2 * dt * rdx * hf * ( ds(icellr) - ds(icelll) ) / hum(iface)
                !
             enddo
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
          qf = 0.
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
             ! compute total mass flux at current face
             !
             uf = 0.
             !
             do k = 1, kmax
                !
                uf = uf + hku(iface,k) * ( theta*u1(iface,k) + (1.-theta)*u0(iface,k) )
                !
             enddo
             !
             uf = rsgn * lf * uf
             !
             if ( uf > 0. ) qf = qf + uf
             !
          enddo
          !
          ! compute the "flow" Courant number
          !
          if ( hs(icell) > epsdry ) then
             !
             area = cell(icell)%attr(CELLAREA)
             !
             cfl = qf * dt / hs(icell) / area
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
    ! build the w-momentum equation
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       ! initialize system of equations in dry points
       !
       do icell = 1, ncells
          !
          if ( wets(icell) /= 1 ) then
             !
             amatw(icell,:,1) = 1.
             amatw(icell,:,2) = 0.
             amatw(icell,:,3) = 0.
             rhsw (icell,:  ) = 0.
             !
          endif
          !
       enddo
       !
       ! bottom:
       !
       ! compute cell-based velocity vector
       !
       call perot ( u1, 1, kmax )
       !
       ! the kinematic condition is imposed
       !
       do icell = 1, ncells
          !
          if ( wets(icell) == 1 ) then
             !
             amatw(icell,kmax,1) = 1.
             amatw(icell,kmax,2) = 0.
             amatw(icell,kmax,3) = 0.
             !
             w1(icell,kmax) = 0.
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
                w1(icell,kmax) = w1(icell,kmax) + rsgn * lf * zku(iface,kmax) * ( nx * (2.*uvc(icell,kmax,1)-uvc(icell,kmax-1,1)) + &
                                                                                  ny * (2.*uvc(icell,kmax,2)-uvc(icell,kmax-1,2)) )
                !
             enddo
             !
             area = cell(icell)%attr(CELLAREA)
             !
             w1(icell,kmax) = w1(icell,kmax) / area
             !
             rhsw(icell,kmax) = w1(icell,kmax)
             !
          endif
          !
       enddo
       !
       ! free surface and interior part:
       !
       ! compute the time derivative
       !
       do k = 0, kmax-1
          !
          do icell = 1, ncells
             !
             if ( wets(icell) == 1 ) then
                !
                amatw(icell,k,1) = 1. / (dt*thetaw)
                rhsw (icell,k  ) = w0(icell,k) / (dt*thetaw)
                !
             endif
             !
          enddo
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
          do k = 0, kmax-1
             !
             kd = max(k,1)
             !
             do icell = 1, ncells
                !
                if ( wets(icell) == 1 ) then
                   !
                   area = cell(icell)%attr(CELLAREA)
                   !
                   fac = 0.5 * area * ( hks(icell,kd) + hks(icell,k+1) )
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
                      fac1 = lf * hku(iface,kd ) * u0(iface,kd )
                      fac2 = lf * hku(iface,k+1) * u0(iface,k+1)
                      !
                      qf = ( fac2*hku(iface,kd) + fac1*hku(iface,k+1) ) / ( fac*( hku(iface,kd) + hku(iface,k+1) ) )
                      !
                      if ( propsc == 3 .and. kappa == 1. ) then
                         !
                         ! central differences
                         !
                         rhsw(icell,k) = rhsw(icell,k) - rsgn * qf * ( 0.5 * (w0(icellr,k) + w0(icelll,k)) - w0(icell,k) )
                         !
                      else
                         !
                         ! first order upwind scheme
                         !
                         if ( qf > 0. ) then
                            !
                            rhsw(icell,k) = rhsw(icell,k) - rsgn * qf * ( w0(icelll,k) - w0(icell,k) )
                            !
                         else
                            !
                            rhsw(icell,k) = rhsw(icell,k) - rsgn * qf * ( w0(icellr,k) - w0(icell,k) )
                            !
                         endif
                         !
                         ! add higher order (flux-limited) correction, if appropriate
                         !
                         if ( propsc /= 1 ) then
                            !
                            rdx = face(iface)%attr(FACEDISTC)
                            !
                            lwfac = 1. - rdx * dt * abs(u0(iface,k))
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
                                  w0u = w0u + area * w0(icella,k)
                                  !
                                  totarea = totarea + area
                                  !
                               enddo
                               !
                               w0u = w0u / totarea
                               !
                               ! compute solution gradients
                               !
                               grad1 = w0(icellr,k) - w0(icelll,k)
                               grad2 = w0(icelll,k) - w0u
                               !
                               rhsw(icell,k) = rhsw(icell,k) - 0.5 * max(0.,lwfac) * rsgn * qf * fluxlim(grad1,grad2)
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
                                  w0u = w0u + area * w0(icella,k)
                                  !
                                  totarea = totarea + area
                                  !
                               enddo
                               !
                               w0u = w0u / totarea
                               !
                               ! compute solution gradients
                               !
                               grad1 = w0(icelll,k) - w0(icellr,k)
                               grad2 = w0(icellr,k) - w0u
                               !
                               rhsw(icell,k) = rhsw(icell,k) - 0.5 * max(0.,lwfac) * rsgn * qf * fluxlim(grad1,grad2)
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
          enddo
          !
       endif
       !
       if ( verwinc ) then
          !
          ! compute vertical terms (implicit)
          !
          propsc = nint(pnums(41))
          kappa  = pnums(42)
          mbound = pnums(43)
          phieby = pnums(44)
          !
          do k = 0, kmax-1
             !
             kd  = max(k  ,1   )
             kdd = max(k-1,1   )
             kuu = min(k+2,kmax)
             !
             do icell = 1, ncells
                !
                if ( wets(icell) == 1 ) then
                   !
                   ! advection term
                   !
                   ctrkt = ( wom(icell,kd-1) + wom(icell,kd  ) ) / ( hks(icell,kd) + hks(icell,k+1) )
                   ctrkb = ( wom(icell,k   ) + wom(icell,k +1) ) / ( hks(icell,kd) + hks(icell,k+1) )
                   !
                   if ( k == 0 ) ctrkt = 0.
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      amatw(icell,k,1) = amatw(icell,k,1) - 0.5*ctrkt + 0.5*ctrkb
                      amatw(icell,k,2) =  0.5*ctrkt
                      amatw(icell,k,3) = -0.5*ctrkb
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      amatw(icell,k,1) = amatw(icell,k,1) - ctrkt + ctrkb + max(ctrkt,0.) - min(ctrkb,0.)
                      amatw(icell,k,2) =  min(ctrkt,0.)
                      amatw(icell,k,3) = -max(ctrkb,0.)
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         if ( ctrkt > 0. ) then
                            !
                            grad1 = w0(icell,kd-1) - w0(icell,k  )
                            grad2 = w0(icell,k   ) - w0(icell,k+1)
                            !
                            rhsw(icell,k) = rhsw(icell,k) - 0.5 * ctrkt * fluxlim(grad1,grad2)
                            !
                         else if ( ctrkt < 0. ) then
                            !
                            grad1 = w0(icell,kd -1) - w0(icell,k   )
                            grad2 = w0(icell,kdd-1) - w0(icell,kd-1)
                            !
                            rhsw(icell,k) = rhsw(icell,k) + 0.5 * ctrkt * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                         if ( ctrkb > 0. ) then
                            !
                            grad1 = w0(icell,k  ) - w0(icell,k  +1)
                            grad2 = w0(icell,k+1) - w0(icell,kuu  )
                            !
                            rhsw(icell,k) = rhsw(icell,k) + 0.5 * ctrkb * fluxlim(grad1,grad2)
                            !
                         else if ( ctrkb < 0. ) then
                            !
                            grad1 = w0(icell,k   ) - w0(icell,k+1)
                            grad2 = w0(icell,kd-1) - w0(icell,k  )
                            !
                            rhsw(icell,k) = rhsw(icell,k) - 0.5 * ctrkb * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                      endif
                      !
                   endif
                   !
                   ! viscosity term
                   !
                   if ( iturb < 2 ) then
                      !
                      ctrkt = ( vnu3d(icell,kd-1) + vnu3d(icell,kd  ) ) / ( hks(icell,kd  )*( hks(icell,kd) + hks(icell,k+1) ) )
                      ctrkb = ( vnu3d(icell,k   ) + vnu3d(icell,k +1) ) / ( hks(icell,k +1)*( hks(icell,kd) + hks(icell,k+1) ) )
                      !
                      if ( k == 0 ) ctrkt = 0.
                      !
                      amatw(icell,k,1) = amatw(icell,k,1) + ctrkt + ctrkb
                      amatw(icell,k,2) = amatw(icell,k,2) - ctrkt
                      amatw(icell,k,3) = amatw(icell,k,3) - ctrkb
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          enddo
          !
       endif
       !
       ! build gradient matrix for non-hydrostatic pressure
       !
       if ( ihydro == 1 ) then
          !
          ! Keller-box scheme, so non-hydrostatic pressure is located at the centers of layer interfaces
          !
          do icell = 1, ncells
             !
             if ( wets(icell) == 1 ) then
                !
                gmatw(icell,:,1) =  2./hks(icell,:)
                gmatw(icell,:,2) = -gmatw(icell,:,1)
                !
                do k = 2, kmax
                   !
                   fac = 1.
                   !
                   do j = 1, kmax-k+1
                      !
                      fac = -fac
                      !
                      gmatw(icell,k-1,2*j+1) =  2.*fac/hks(icell,j+k-1)
                      gmatw(icell,k-1,2*j+2) = -gmatw(icell,k-1,2*j+1)
                      !
                   enddo
                   !
                enddo
                !
             else
                !
                gmatw(icell,:,:) = 0.
                !
             endif
             !
             gmatw(icell,1,1) = 0.
             !
          enddo
          !
          ! to reduce the pressure Poisson equation set pressure of bottom face to that of top face for a number of layers
          !
          do l = 1, qlay
             !
             do icell = 1, ncells
                !
                do k = 1, kmax
                   !
                   j = kmax +1 - k - l
                   if ( j < 0 ) cycle
                   !
                   gmatw(icell,k,2*j+1) = gmatw(icell,k,2*j+1) + gmatw(icell,k,2*j+2)
                   gmatw(icell,k,2*j+2) = 0.
                   !
                enddo
                !
             enddo
             !
          enddo
          !
       else if ( ihydro == 2 ) then
          !
          ! central differences, so non-hydrostatic pressure is located at the cell centers
          !
          do icell = 1, ncells
             !
             if ( wets(icell) == 1 ) then
                !
                do k = 1, kmax
                   !
                   kd = max(k-1,1)
                   !
                   gmatw(icell,k,1) =  2./(hks(icell,kd) + hks(icell,k))
                   gmatw(icell,k,2) = -gmatw(icell,k,1)
                   !
                enddo
                !
             else
                !
                gmatw(icell,:,1) = 0.
                gmatw(icell,:,2) = 0.
                !
             endif
             !
             gmatw(icell,1,2) = gmatw(icell,1,2) - gmatw(icell,1,1)
             gmatw(icell,1,1) = 0.
             !
          enddo
          !
       endif
       !
       ! compute gradient of non-hydrostatic pressure
       !
       fac = (-1.)**kmax
       !
       if ( iproj == 1 ) then
          fac1 = 1.
       else if ( iproj == 2 ) then
          fac1 = 1. - theta3
       endif
       !
       do k = 0, kmax-1
          !
          kd = max(k,1)
          !
          fac = -fac
          !
          do icell = 1, ncells
             !
             if ( wets(icell) == 1 ) then
                !
                if ( ihydro == 1 ) then
                   !
                   do j = 0, kmax-1
                      !
                      kd = max(k+j  ,1   )
                      kd = min(kd   ,kmax)
                      ku = min(k+j+1,kmax)
                      !
                      rhsw(icell,k) = rhsw(icell,k) - fac1 * ( gmatw(icell,k+1,2*j+1)*q(icell,kd) + gmatw(icell,k+1,2*j+2)*q(icell,ku) )
                      !
                   enddo
                   !
                   !rhsw(icell,k) = rhsw(icell,k) - fac * ( w1(icell,kmax) - w0(icell,kmax) ) / (dt*thetaw)
                   !
                else
                   !
                   rhsw(icell,k) = rhsw(icell,k) - fac1 * ( gmatw(icell,k+1,1)*q(icell,kd) + gmatw(icell,k+1,2)*q(icell,k+1) )
                   !
                endif
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! solve the w-momentum equation
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       do icell = 1, ncells
          !
          bi = 1./amatw(icell,0,1)
          !
          amatw(icell,0,1) = bi
          amatw(icell,0,3) = amatw(icell,0,3)*bi
          rhsw (icell,0  ) = rhsw (icell,0  )*bi
          !
          do k = 1, kmax
             !
             bi = 1./(amatw(icell,k,1) - amatw(icell,k,2)*amatw(icell,k-1,3))
             amatw(icell,k,1) = bi
             amatw(icell,k,3) = amatw(icell,k,3)*bi
             rhsw (icell,k  ) = (rhsw(icell,k) - amatw(icell,k,2)*rhsw(icell,k-1))*bi
             !
          enddo
          !
          w1(icell,kmax) = rhsw(icell,kmax)
          do k = kmax-1, 0, -1
             w1(icell,k) = rhsw(icell,k) - amatw(icell,k,3)*w1(icell,k+1)
          enddo
          !
       enddo
       !
       ! re-update the solution in case of thetaw <> 1
       !
       if ( thetaw /= 1. ) then
          !
          do k = 0, kmax-1
             !
             do icell = 1, ncells
                !
                if ( wets(icell) == 1 ) then
                   !
                   w1(icell,k) = ( w1(icell,k) - (1.-thetaw) * w0(icell,k) ) / thetaw
                   !
                endif
                !
             enddo
             !
          enddo
          !
       endif
       !
    endif
    !
    ! compute the non-hydrostatic pressure correction
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       ! build divergence matrix
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
             do k = 1, kmax
                !
                kd = max(k-1,1   )
                ku = min(k+1,kmax)
                !
                ! compute gradient of layer interfaces of current layer
                !
                zkxd = 0.
                zkyd = 0.
                zkxu = 0.
                zkyu = 0.
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
                   zkxd = zkxd + rsgn * lf * zku(iface,k-1) * nx
                   zkyd = zkyd + rsgn * lf * zku(iface,k-1) * ny
                   !
                   zkxu = zkxu + rsgn * lf * zku(iface,k  ) * nx
                   zkyu = zkyu + rsgn * lf * zku(iface,k  ) * ny
                   !
                enddo
                !
                zkxd = zkxd / area
                zkyd = zkyd / area
                zkxu = zkxu / area
                zkyu = zkyu / area
                !
                if ( k == kmax ) then
                   zkxu = 0.
                   zkyu = 0.
                endif
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
                   fac1 = lf * ( zkxd * ( xf - xc ) + zkyd * ( yf - yc ) ) / ( hkum(iface,k ) + hkum(iface,kd) )
                   fac2 = lf * ( zkxu * ( xf - xc ) + zkyu * ( yf - yc ) ) / ( hkum(iface,ku) + hkum(iface,k ) )
                   !
                   dmatu(icell,k,jf,1) = lf * hku(iface,k) - fac1 * hkum(iface,kd) + fac2 * hkum(iface,ku)
                   dmatu(icell,k,jf,2) =                   - fac1 * hkum(iface,k )
                   dmatu(icell,k,jf,3) =                                           + fac2 * hkum(iface,k )
                   !
                enddo
                !
             enddo
             !
          else
             !
             dmatu(icell,:,:,:) = 0.
             !
          endif
          !
       enddo
       !
       do icell = 1, ncells
          !
          dmatu(icell,1,:,1) = dmatu(icell,1,:,1) + 2.*dmatu(icell,1,:,2)
          dmatu(icell,1,:,3) = dmatu(icell,1,:,3) -    dmatu(icell,1,:,2)
          dmatu(icell,1,:,2) = 0.
          !
       enddo
       !
       ! build the Poisson equation
       !
       do icell = 1, ncells
          !
          area = cell(icell)%attr(CELLAREA)
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
             !
             kwu = 1.
             if ( k == kmax ) kwu = 0.
             !
             if ( ihydro == 1 ) then
                !
                amatp(icell,k,0) = area * ( gmatw(icell,k,2) + gmatw(icell,k,3) - kwu * gmatw(icell,ku,1) )
                !
                do j = 1, kmax-2
                   !
                   amatp(icell,k,ishif(j)) = area * ( gmatw(icell,k,2*j+2) + gmatw(icell,k,2*j+3) - kwu * gmatw(icell,ku,2*j) - kwu * gmatw(icell,ku,2*j+1) )
                   !
                enddo
                !
                amatp(icell,k,ishif(kmax-1)) = area * ( gmatw(icell,k,2*kmax) - kwu * gmatw(icell,ku,2*kmax-2) )
                !
             else
                !
                amatp(icell,k, 0) = area * ( gmatw(icell,k,2) - kwu * gmatw(icell,ku,1) )
                amatp(icell,k,12) = area * (                  - kwu * gmatw(icell,ku,2) )
                amatp(icell,k,16) = 0.
                !
             endif
             !
             amatp(icell,k,4) = 0.
             amatp(icell,k,8) = area * gmatw(icell,k,1)
             rhsp (icell,k  ) = area * ( w1(icell,k-1) - kwu*w1(icell,k) )
             !
             ! loop over faces of the cell
             !
             do jf = 1, cell(icell)%nof
                !
                ! face identifier
                !
                iface = cell(icell)%face(jf)%atti(FACEID)
                !
                ! consider left and right cells of the face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                ! build matrix and right-hand side
                !
                if ( icell == icelll ) then
                   !
                   if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
                      !
                      amatp(icell,k,    0) = amatp(icell,k, 0) + dmatu(icell,k,jf,1) * gmatu(iface,k ,2) + dmatu(icell,k,jf,2) * gmatu(iface,kd,5) +  &
                                                                 dmatu(icell,k,jf,3) * gmatu(iface,ku,1)
                      amatp(icell,k,   jf) =                     dmatu(icell,k,jf,1) * gmatu(iface,k ,4) + dmatu(icell,k,jf,2) * gmatu(iface,kd,6) +  &
                                                                 dmatu(icell,k,jf,3) * gmatu(iface,ku,3)
                      amatp(icell,k,    4) = amatp(icell,k, 4) + dmatu(icell,k,jf,2) * gmatu(iface,kd,1)
                      amatp(icell,k, 4+jf) =                     dmatu(icell,k,jf,2) * gmatu(iface,kd,3)
                      amatp(icell,k,    8) = amatp(icell,k, 8) + dmatu(icell,k,jf,1) * gmatu(iface,k ,1) + dmatu(icell,k,jf,2) * gmatu(iface,kd,2)
                      amatp(icell,k, 8+jf) =                     dmatu(icell,k,jf,1) * gmatu(iface,k ,3) + dmatu(icell,k,jf,2) * gmatu(iface,kd,4)
                      amatp(icell,k,   12) = amatp(icell,k,12) + dmatu(icell,k,jf,1) * gmatu(iface,k ,5) + dmatu(icell,k,jf,3) * gmatu(iface,ku,2)
                      amatp(icell,k,12+jf) =                     dmatu(icell,k,jf,1) * gmatu(iface,k ,6) + dmatu(icell,k,jf,3) * gmatu(iface,ku,4)
                      amatp(icell,k,   16) = amatp(icell,k,16) + dmatu(icell,k,jf,3) * gmatu(iface,ku,5)
                      amatp(icell,k,16+jf) =                     dmatu(icell,k,jf,3) * gmatu(iface,ku,6)
                      !
                   else if ( wlimp(iface) ) then   ! described water level at boundary face
                      !
                      amatp(icell,k, 0) = amatp(icell,k, 0) + dmatu(icell,k,jf,1) * gmatu(iface,k ,2) + dmatu(icell,k,jf,3) * gmatu(iface,ku,1)
                      amatp(icell,k, 4) = amatp(icell,k, 4) + dmatu(icell,k,jf,2) * gmatu(iface,kd,1)
                      amatp(icell,k, 8) = amatp(icell,k, 8) + dmatu(icell,k,jf,1) * gmatu(iface,k ,1) + dmatu(icell,k,jf,2) * gmatu(iface,kd,2)
                      amatp(icell,k,12) = amatp(icell,k,12) + dmatu(icell,k,jf,3) * gmatu(iface,ku,2)
                      !
                   endif
                   !
                   rhsp(icell,k) = rhsp(icell,k) + dmatu(icell,k,jf,1) * u1(iface,k) + dmatu(icell,k,jf,2) * u1(iface,kd) + dmatu(icell,k,jf,3) * u1(iface,ku)
                   !
                else if ( icell == icellr ) then
                   !
                   if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
                      !
                      amatp(icell,k,    0) = amatp(icell,k, 0) - dmatu(icell,k,jf,1) * gmatu(iface,k ,4) - dmatu(icell,k,jf,2) * gmatu(iface,kd,6) +  &
                                                               - dmatu(icell,k,jf,3) * gmatu(iface,ku,3)
                      amatp(icell,k,   jf) =                   - dmatu(icell,k,jf,1) * gmatu(iface,k ,2) - dmatu(icell,k,jf,2) * gmatu(iface,kd,5) +  &
                                                               - dmatu(icell,k,jf,3) * gmatu(iface,ku,1)
                      amatp(icell,k,    4) = amatp(icell,k, 4) - dmatu(icell,k,jf,2) * gmatu(iface,kd,3)
                      amatp(icell,k, 4+jf) =                   - dmatu(icell,k,jf,2) * gmatu(iface,kd,1)
                      amatp(icell,k,    8) = amatp(icell,k, 8) - dmatu(icell,k,jf,1) * gmatu(iface,k ,3) - dmatu(icell,k,jf,2) * gmatu(iface,kd,4)
                      amatp(icell,k, 8+jf) =                   - dmatu(icell,k,jf,1) * gmatu(iface,k ,1) - dmatu(icell,k,jf,2) * gmatu(iface,kd,2)
                      amatp(icell,k,   12) = amatp(icell,k,12) - dmatu(icell,k,jf,1) * gmatu(iface,k ,6) - dmatu(icell,k,jf,3) * gmatu(iface,ku,4)
                      amatp(icell,k,12+jf) =                   - dmatu(icell,k,jf,1) * gmatu(iface,k ,5) - dmatu(icell,k,jf,3) * gmatu(iface,ku,2)
                      amatp(icell,k,   16) = amatp(icell,k,16) - dmatu(icell,k,jf,3) * gmatu(iface,ku,6)
                      amatp(icell,k,16+jf) =                   - dmatu(icell,k,jf,3) * gmatu(iface,ku,5)
                      !
                   else if ( wlimp(iface) ) then   ! described water level at boundary face
                      !
                      amatp(icell,k, 0) = amatp(icell,k, 0) - dmatu(icell,k,jf,1) * gmatu(iface,k ,2) - dmatu(icell,k,jf,3) * gmatu(iface,ku,1)
                      amatp(icell,k, 4) = amatp(icell,k, 4) - dmatu(icell,k,jf,2) * gmatu(iface,kd,1)
                      amatp(icell,k, 8) = amatp(icell,k, 8) - dmatu(icell,k,jf,1) * gmatu(iface,k ,1) - dmatu(icell,k,jf,2) * gmatu(iface,kd,2)
                      amatp(icell,k,12) = amatp(icell,k,12) - dmatu(icell,k,jf,3) * gmatu(iface,ku,2)
                      !
                   endif
                   !
                   rhsp(icell,k) = rhsp(icell,k) - dmatu(icell,k,jf,1) * u1(iface,k) - dmatu(icell,k,jf,2) * u1(iface,kd) - dmatu(icell,k,jf,3) * u1(iface,ku)
                   !
                endif
                !
             enddo
             !
             rhsp(icell,k) = rhsp(icell,k) / (dt*theta3)
             !
          enddo
          !
       enddo
       !
       ! reduce the pressure Poisson equation
       !
       do l = 1, qlay
          !
          do icell = 1, ncells
             !
             amatp(icell,qmax, 0) = amatp(icell,qmax, 0) + real(qlay+2-l)*amatp(icell,kmax-l+1, 8)
             amatp(icell,qmax, 1) = amatp(icell,qmax, 1) + real(qlay+2-l)*amatp(icell,kmax-l+1, 9)
             amatp(icell,qmax, 2) = amatp(icell,qmax, 2) + real(qlay+2-l)*amatp(icell,kmax-l+1,10)
             amatp(icell,qmax, 3) = amatp(icell,qmax, 3) + real(qlay+2-l)*amatp(icell,kmax-l+1,11)
             amatp(icell,qmax, 8) = amatp(icell,qmax, 8) + real(qlay+2-l)*amatp(icell,kmax-l+1, 4)
             amatp(icell,qmax, 9) = amatp(icell,qmax, 9) + real(qlay+2-l)*amatp(icell,kmax-l+1, 5)
             amatp(icell,qmax,10) = amatp(icell,qmax,10) + real(qlay+2-l)*amatp(icell,kmax-l+1, 6)
             amatp(icell,qmax,11) = amatp(icell,qmax,11) + real(qlay+2-l)*amatp(icell,kmax-l+1, 7)
             rhsp (icell,qmax   ) = rhsp (icell,qmax   ) + real(qlay+2-l)*rhsp (icell,kmax-l+1   )
             !
          enddo
          !
       enddo
       !
       do icell = 1, ncells
          !
          do k = 1, kmax
             !
             if ( .not. amatp(icell,k,0) /= 0. ) then
                amatp(icell,k,:) =  0.
                amatp(icell,k,0) = -1.
                rhsp (icell,k  ) =  0.
             endif
             !
          enddo
          !
       enddo
       !
       ! solve the Poisson equation
       !
       call bicgstab3 ( amatp(1:ncells,1:qmax,0:nconct), rhsp(1:ncells,1:qmax), dq(1:ncells,1:qmax) )
       if (STPNOW()) return
       !
       do k = qmax+1, kmax
          !
          dq(:,k) = dq(:,qmax)
          !
       enddo
       !
    endif
    !
    ! update the non-hydrostatic pressure
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       if ( iproj == 1 ) then
          q = q + dq
       else if ( iproj == 2 ) then
          q = dq
       endif
       !
    endif
    !
    ! correct the flow velocities
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       ! u-velocity
       !
       do k = 1, kmax
          !
          kd = max(k-1,1   )
          ku = min(k+1,kmax)
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
                u1(iface,k) = u1(iface,k) - dt*theta3 * ( gmatu(iface,k,1)*dq(icelll,kd) + gmatu(iface,k,2)*dq(icelll,k ) + gmatu(iface,k,3)*dq(icellr,kd)  &
                                                        + gmatu(iface,k,4)*dq(icellr,k ) + gmatu(iface,k,5)*dq(icelll,ku) + gmatu(iface,k,6)*dq(icellr,ku) )
                !
             else if ( wlimp(iface) .and. wetu(iface) == 1 ) then   ! described water level at boundary face
                !
                ! consider boundary cell of current face
                !
                icellb = face(iface)%atti(FACEC1)
                !
                u1(iface,k) = u1(iface,k) - dt*theta3 * ( gmatu(iface,k,1)*dq(icellb,kd) + gmatu(iface,k,2)*dq(icellb,k) )
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! w-velocity
       !
       do k = 0, kmax-1
          !
          kd = max(k,1)
          !
          do icell = 1, ncells
             !
             if ( wets(icell) == 1 ) then
                !
                if ( ihydro == 1 ) then
                   !
                   do j = 0, kmax-1
                      !
                      kd = max(k+j  ,1   )
                      kd = min(kd   ,kmax)
                      ku = min(k+j+1,kmax)
                      !
                      w1(icell,k) = w1(icell,k) - dt*theta3*( gmatw(icell,k+1,2*j+1)*dq(icell,kd) + gmatw(icell,k+1,2*j+2)*dq(icell,ku) )
                      !
                   enddo
                   !
                else
                   !
                   w1(icell,k) = w1(icell,k) - dt*theta3*( gmatw(icell,k+1,1)*dq(icell,kd) + gmatw(icell,k+1,2)*dq(icell,k+1) )
                   !
                endif
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! compute cell-based velocity vector
       !
       call perot ( u1, 1, kmax )
       !
       do icell = 1, ncells
          !
          if ( wets(icell) == 1 ) then
             !
             w1(icell,kmax) = 0.
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
                w1(icell,kmax) = w1(icell,kmax) + rsgn * lf * zku(iface,kmax) * ( nx * (2.*uvc(icell,kmax,1)-uvc(icell,kmax-1,1)) + &
                                                                                  ny * (2.*uvc(icell,kmax,2)-uvc(icell,kmax-1,2)) )
                !
             enddo
             !
             area = cell(icell)%attr(CELLAREA)
             !
             w1(icell,kmax) = w1(icell,kmax) / area
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
       do k = 1, kmax
          !
          u = theta*u1(iface,k) + (1.-theta)*u0(iface,k)
          !
          qn(iface,k) = lf * hku(iface,k) * u
          !
       enddo
       !
    enddo
    !
    ! compute the depth-averaged flow velocity
    !
    udep = 0.
    !
    do iface = 1, nfaces
       !
       if ( wetu(iface) == 1 ) then
          !
          do k = 1, kmax
             !
             udep(iface) = udep(iface) + hku(iface,k) * u1(iface,k)
             !
          enddo
          !
          udep(iface) = udep(iface) / hu(iface)
          !
       endif
       !
    enddo
    !
    ! calculate net mass outflow based on local continuity equation
    !
    if ( ITEST >= 30 ) then
       !
       moutf = 0.
       !
       if ( ihydro == 1 .or. ihydro == 2 ) then
          !
          do icell = 1, ncells
             !
             area = cell(icell)%attr(CELLAREA)
             !
             do k = 1, kmax
                !
                kd = max(k-1,1   )
                ku = min(k+1,kmax)
                !
                kwu = 1.
                if ( k == kmax ) kwu = 0.
                !
                moutf = area * ( w1(icell,k-1) - kwu*w1(icell,k) )
                !
                ! loop over faces of the cell
                !
                do jf = 1, cell(icell)%nof
                   !
                   ! face identifier
                   !
                   iface = cell(icell)%face(jf)%atti(FACEID)
                   !
                   ! consider left and right cells of the face
                   !
                   icelll = face(iface)%atti(FACECL)
                   icellr = face(iface)%atti(FACECR)
                   !
                   if ( icell == icelll ) then
                      !
                      moutf = moutf + dmatu(icell,k,jf,1) * u1(iface,k) + dmatu(icell,k,jf,2) * u1(iface,kd) + dmatu(icell,k,jf,3) * u1(iface,ku)
                      !
                   else if ( icell == icellr ) then
                      !
                      moutf = moutf - dmatu(icell,k,jf,1) * u1(iface,k) - dmatu(icell,k,jf,2) * u1(iface,kd) - dmatu(icell,k,jf,3) * u1(iface,ku)
                      !
                   endif
                   !
                enddo
                !
             enddo
             !
          enddo
          !
       endif
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
          ener = ener + 0.5 * lf * hu(iface) * udep(iface) * udep(iface) / rdx
          !
       enddo
       !
       write(PRINTF,102) vol
       write(PRINTF,103) ener
       !
    endif
    !
    ! determine layer interfaces and layer thicknesses at new time level
    !
    work(:,1) = s1(:) + dps(:)
    !
    zksnew(:,0) = -dps(:) + work(:,1)
    !
    call sigmacoor ( zksnew, work, ncells )
    !
    do k = 1, kmax
       !
       hksnew(:,k) = zksnew(:,k-1) - zksnew(:,k)
       !
    enddo
    !
    ! compute the relative vertical velocity
    !
    if ( ihydro == 0 ) then
       !
       wom(:,kmax) = 0.
       !
       do k = kmax-1, 0, -1
          !
          do icell = 1, ncells
             !
             if ( wets(icell) == 1 ) then
                !
                qf = 0.
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
                   qf = qf + rsgn * lf * hku(iface,k+1) * ( theta * u1(iface,k+1) + (1.-theta) * u0(iface,k+1) )
                   !
                enddo
                !
                ! area of cell
                !
                area = cell(icell)%attr(CELLAREA)
                !
                wom(icell,k) = wom(icell,k+1) - qf / area - ( hksnew(icell,k+1) - hks(icell,k+1) ) / dt
                !
             else
                !
                wom(icell,k) = 0.
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! check if relative vertical velocity at surface is zero
       !
       do icell = 1, ncells
          !
          if ( abs(wom(icell,0)) > epswom ) then
             !
             write (msgstr,'(a,i8,a,e9.3,a)') 'nonzero relative vertical velocity at surface in cell=',icell,'; omega = ',wom(icell,0),' m/s'
             call msgerr (2, trim(msgstr) )
             !
             wom(icell,0) = 0.
             !
          endif
          !
       enddo
       !
    else
       !
       call perot ( u1, 1, kmax )
       !
       do icell = 1, ncells
          !
          if ( wets(icell) == 1 ) then
             !
             do k = 1, kmax-1
                !
                u = ( uvc(icell,k+1,1) * hks(icell,k) + uvc(icell,k,1) * hks(icell,k+1) ) / ( hks(icell,k) + hks(icell,k+1) )
                v = ( uvc(icell,k+1,2) * hks(icell,k) + uvc(icell,k,2) * hks(icell,k+1) ) / ( hks(icell,k) + hks(icell,k+1) )
                !
                qf = 0.
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
                   qf = qf + rsgn * lf * zku(iface,k) * ( nx * u + ny * v )
                   !
                enddo
                !
                area = cell(icell)%attr(CELLAREA)
                !
                wom(icell,k) = w1(icell,k) - qf / area - ( zksnew(icell,k) - zks(icell,k) ) / dt
                !
             enddo
             !
             wom(icell,0   ) = 0.
             wom(icell,kmax) = 0.
             !
          else if ( brks(icell) == 1 ) then
             !
             ! hydrostatic pressure is assumed at steep front of breaking wave, so relative vertical velocity is derived from local continuity equation
             !
             wom(icell,kmax) = 0.
             !
             do k = kmax-1, 0, -1
                !
                qf = 0.
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
                   qf = qf + rsgn * lf * hku(iface,k+1) * ( theta * u1(iface,k+1) + (1.-theta) * u0(iface,k+1) )
                   !
                enddo
                !
                ! area of cell
                !
                area = cell(icell)%attr(CELLAREA)
                !
                wom(icell,k) = wom(icell,k+1) - qf / area - ( hksnew(icell,k+1) - hks(icell,k+1) ) / dt
                !
             enddo
             !
             if ( abs(wom(icell,0)) > epswom ) then
                !
                write (msgstr,'(a,i8,a,e9.3,a)') 'nonzero relative vertical velocity at surface in cell=',icell,'; omega = ',wom(icell,0),' m/s'
                if ( ITEST >= 50 ) call msgerr (1, trim(msgstr) )
                !
                wom(icell,0) = 0.
                !
             endif
             !
          else
             !
             wom(icell,:) = 0.
             !
          endif
          !
       enddo
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
             u1(iface,1:kmax) = (1.-gamma) * u1(iface,1:kmax) + gamma * u1(bface,1:kmax)
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
    do k = 1, kmax
       !
       do iface = 1, nfaces
          !
          if ( wetu(iface) == 1 ) then                    ! only consider wet face
             !
             rdx = face(iface)%attr(FACEDISTG)
             !
             if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
                !
                cfl = rdx * dt * abs(u1(iface,k))
                if ( cfl > cflmax ) cflmax = cfl
                !
             endif
             !
          endif
          !
       enddo
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
end subroutine SwashImpLayUflow
