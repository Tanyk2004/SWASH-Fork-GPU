subroutine SwashExpLay1DHflow ( ibl, ibr, kgrpnt )
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
!    6.01: Panagiotis Vasarmidis
!    6.02: Tom Bogaard
!    8.01: Marcel Zijlema
!   10.05: Dirk Rijnsdorp
!
!   Updates
!
!    1.00,    July 2010: New subroutine
!    6.01,    June 2019: extension internal wave generation
!    6.02,    July 2019: extension 3D (non)linear k-eps model
!    8.01, October 2021: Hancock scheme for global continuity equation
!   10.05,    July 2023: extension ambient current
!
!   Purpose
!
!   Performs the time integration for the non-hydrostatic, layer-averaged 1D shallow water equations
!
!   Method
!
!   The time integration with respect to the continuity equation and the water level gradient of the
!   u-momentum equation is based on the leap-frog technique (or the Hansen scheme).
!
!   The time integration with respect to the horizontal advective term is based on the predictor-corrector
!   scheme of MacCormack, while that for the bottom friction is based on Euler implicit and for the
!   non-hydrostatic pressure gradient a semi-implicit approach is employed (theta-scheme). Both vertical
!   advective and viscosity terms are treated semi-implicit as well. This results in a tri-diagonal system.
!
!   The space discretization of the horizontal advective term is momentum conservative (or energy head
!   conservative in case of flow contraction) and is approximated by either first order upwind or higher
!   order (flux-limited) scheme (CDS, Fromm, BDF, QUICK, MUSCL, Koren, etc.). The higher order scheme is
!   treated by defect correction consistent with the MacCormack scheme.
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
    use SwashCommdata3
    use SwashTimecomm
    use m_genarr, only: work
    use m_parall
    use SwashSolvedata, only: icond
    use SwashFlowdata, ibltmp => ibl, &
                       ibrtmp => ibr
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                 :: ibl    ! boundary condition type at left boundary
    integer, intent(in)                 :: ibr    ! boundary condition type at right boundary
    !
    integer, dimension(mxc), intent(in) :: kgrpnt ! index table containing the address of each (active) grid point
                                                  ! =1: not active grid point
                                                  ! >1: active grid point
!
!   Parameter variables
!
    real   , parameter :: epswom  = 0.001 ! tolerance for relative vertical velocity at surface
!
!   Local variables
!
    integer, save                   :: ient = 0 ! number of entries in this subroutine
    integer                         :: j        ! loop counter
    integer                         :: k        ! loop counter over vertical layers
    integer                         :: kd       ! index of layer k-1
    integer                         :: kdd      ! index of layer k-2
    integer                         :: ku       ! index of layer k+1
    integer                         :: kuu      ! index of layer k+2
    integer                         :: l        ! loop counter
    integer                         :: m        ! loop counter over horizontal grid points
    integer                         :: md       ! index of point m-1
    integer                         :: mdd      ! index of point m-2
    integer                         :: mend     ! end index of loop over u-points
    integer                         :: msta     ! start index of loop over u-points
    integer                         :: mu       ! index of point m+1
    integer                         :: muu      ! index of point m+2
    integer                         :: nm       ! pointer to m
    integer                         :: nmd      ! pointer to m-1
    integer                         :: nmdd     ! pointer to m-2
    integer                         :: nmf      ! pointer to mf
    integer                         :: nmfu     ! pointer to mfu
    integer                         :: nmfuu    ! pointer to mfuu
    integer                         :: nml      ! pointer to ml
    integer                         :: nmld     ! pointer to mld
    integer                         :: nmlu     ! pointer to mlu
    integer                         :: nmu      ! pointer to m+1
    integer                         :: nmuu     ! pointer to m+2
    !
    real                            :: bi       ! inverse of main diagonal of the matrix
    real                            :: cfl      ! CFL number
    real                            :: ctrkb    ! contribution of vertical terms below considered point
    real                            :: ctrkt    ! contribution of vertical terms above considered point
    real                            :: dep      ! local water depth
    real                            :: ener     ! total energy of closed system
    real                            :: fac      ! a factor
    real                            :: fac1     ! another factor
    real                            :: fac2     ! some other factor
    real                            :: fh       ! factor for Hancock/RK3 scheme
    real                            :: fluxlim  ! flux limiter
    real                            :: grad1    ! solution gradient
    real                            :: grad2    ! another solution gradient
    real                            :: kwd      ! =1. if layer k-1 exists otherwise 0.
    real                            :: kwu      ! =1. if layer k+1 exists otherwise 0.
    real                            :: moutf    ! net mass outflow
    real                            :: rdx      ! reciprocal of mesh size
    real                            :: rhou     ! density of water in velocity point
    real                            :: theta    ! implicitness factor for non-hydrostatic pressure gradient
    real                            :: thetamb  ! implicitness factor for ambient current in global continuity eq
    real                            :: thetau   ! implicitness factor for vertical terms in u-momentum equation
    real                            :: thetaw   ! implicitness factor for vertical terms in w-momentum equation
    real                            :: u        ! u-velocity at point different from its point of definition
    real                            :: utot     ! velocity magnitude
    real                            :: vol      ! total displaced volume of water
    real                            :: w        ! w-velocity at point different from its point of definition
    !
    logical                         :: slimp    ! water level prescribed at left boundary
    logical                         :: srimp    ! water level prescribed at right boundary
    logical                         :: STPNOW   ! indicates that program must stop
    !
    character(80)                   :: msgstr   ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashExpLay1DHflow')
    !
    nmf   = kgrpnt(mf   )
    nmfu  = kgrpnt(mfu  )
    nmfuu = kgrpnt(mfu+1)
    nml   = kgrpnt(ml   )
    nmld  = kgrpnt(ml -1)
    nmlu  = kgrpnt(mlu  )
    !
    rdx = 1./dx
    !
    theta   = pnums( 5)
    thetau  = pnums(31)
    thetaw  = pnums(32)
    thetamb = pnums(71)
    !
    ! check if water level is prescribed at the boundaries
    !
    if ( ibl == 2 .or. ibl == 10 ) then
       !
       slimp = .true.
       !
    else
       !
       slimp = .false.
       !
    endif
    !
    if ( ibr == 2 .or. ibr == 10 ) then
       !
       srimp = .true.
       !
    else
       !
       srimp = .false.
       !
    endif
    !
    ! build the u-momentum equation
    !
    if ( slimp ) then    ! described water level at left boundary point
       !
       msta = mf
       !
    else
       !
       msta = mf + 1
       !
    endif
    !
    if ( srimp ) then    ! described water level at right boundary point
       !
       mend = ml
       !
    else
       !
       mend = ml - 1
       !
    endif
    !
    ! initialize system of equations in dry points
    !
    do m = msta, mend
       !
       nm = kgrpnt(m)
       !
       if ( wetu(nm) /= 1 ) then
          !
          amatu(nm,:,1) = 1.
          amatu(nm,:,2) = 0.
          amatu(nm,:,3) = 0.
          rhsu (nm,:  ) = 0.
          !
       endif
       !
    enddo
    !
    ! compute the time derivative
    !
    do k = 1, kmax
       !
       do m = msta, mend
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             fac = dt * thetau
             !
             amatu(nm,k,1) = 1. / fac
             rhsu (nm,k  ) = u0(nm,k) / fac
             !
          endif
          !
       enddo
       !
    enddo
    !
    ! compute the mass flux
    !
    do k = 1, kmax
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          qx(nm,k) = hku(nm,k)*u0(nm,k)
          !
       enddo
       !
    enddo
    !
    ! compute the mass flux in wl-point based on averaging
    !
    do k = 1, kmax
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          qm(nm,k) = 0.5 * ( qx(nm,k) + qx(nmd,k) )
          !
       enddo
       !
    enddo
    !
    ! compute the advective velocity in wl-point based on first order upwind
    !
    do k = 1, kmax
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac = 0.5 * ( u0(nm,k) + u0(nmd,k) )
             !
          else
             !
             fac = qm(nm,k)
             !
          endif
          !
          if ( fac > 1.0e-5 ) then
             !
             ua(nm,k) = u0(nmd,k)
             !
          else if ( fac < -1.0e-5 ) then
             !
             ua(nm,k) = u0(nm ,k)
             !
          endif
          !
       enddo
       !
    enddo
    !
    ! compute horizontal advection term (momentum/energy head conservative) at internal u-point (prediction)
    !
    do k = 1, kmax
       !
       do m = mf+1, ml-1
          !
          md = m - 1
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          nmu = kgrpnt(mu)
          !
          if ( linswe ) then
             !
             fac1 = 0.
             fac2 = 0.
             !
          else if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac1 = 0.5 * ( ua(nmu,k) + ua(nm,k) )
             fac2 = fac1
             !
          else if ( hkumn(nm,k) > 0. ) then
             !
             fac1 = qm(nmu,k) / hkumn(nm,k)
             fac2 = qm(nm ,k) / hkumn(nm,k)
             !
          else
             !
             fac1 = 0.
             fac2 = 0.
             !
          endif
          !
          if ( wetu(nm) == 1 ) then
             !
             rhsu(nm,k) = rhsu(nm,k) - rdx * ( fac1 * (ua(nmu,k) - u0(nm,k)) - fac2 * (ua(nm,k) - u0(nm,k)) )
             !
          endif
          !
       enddo
       !
    enddo
    !
    if ( slimp .and. .not.linswe ) then    ! described water level at left boundary point
       !
       do k = 1, kmax
          !
          if ( qm(nmfu,k) < -1.0e-5 ) then
             !
             rhsu(nmf,k) = rhsu(nmf,k) - rdx * qm(nmfu,k) * (u0(nmfu,k) - u0(nmf,k)) / hkum(nmf,k)
             !
          endif
          !
       enddo
       !
    endif
    !
    if ( srimp .and. .not.linswe ) then    ! described water level at right boundary point
       !
       do k = 1, kmax
          !
          if ( qm(nml,k) > 1.0e-5 ) then
             !
             rhsu(nml,k) = rhsu(nml,k) - rdx * qm(nml,k) * (u0(nml,k) - u0(nmld,k)) / hkum(nml,k)
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute flow resistance inside porous medium, if appropriate
    !
    if ( iporos == 1 ) then
       !
       do k = 1, kmax
          !
          do m = msta, mend
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) == 1 ) then
                !
                w = 0.25 * ( w0(nm,k-1) + w0(nmu,k-1) + w0(nm,k) + w0(nmu,k) )
                !
                utot = sqrt( u0(nm,k)*u0(nm,k) + w*w )
                !
                amatu(nm,k,1) = amatu(nm,k,1) + apomu(nm,k) + bpomu(nm,k) * utot
                !
                fac = dt * thetau
                !
                amatu(nm,k,1) = amatu(nm,k,1) + cpomu(nm,k) / fac
                rhsu (nm,k  ) = rhsu (nm,k  ) + cpomu(nm,k) * u0(nm,k) / fac
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! compute friction due to vegetation, if appropriate
    !
    if ( iveg /= 0 ) then
       !
       do k = 1, kmax
          !
          do m = msta, mend
             !
             nm = kgrpnt(m)
             !
             if ( wetu(nm) == 1 ) then
                !
                amatu(nm,k,1) = amatu(nm,k,1) + cvegu(nm,k,1) * abs(u0(nm,k))
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! add mass due to inertia, if appropriate
       !
       if ( cvm > 0. ) then
          !
          do k = 1, kmax
             !
             do m = msta, mend
                !
                nm = kgrpnt(m)
                !
                if ( wetu(nm) == 1 ) then
                   !
                   fac = dt * thetau
                   !
                   amatu(nm,k,1) = amatu(nm,k,1) + cvegu(nm,k,2) / fac
                   rhsu (nm,k  ) = rhsu (nm,k  ) + cvegu(nm,k,2) * u0(nm,k) / fac
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
    ! compute water level gradient
    !
    do m = msta, mend
       !
       mu = m + 1
       !
       nm  = kgrpnt(m )
       nmu = kgrpnt(mu)
       !
       if ( wetu(nm) == 1 ) then
          !
          rhsu(nm,:) = rhsu(nm,:) - grav * rdx * (s0(nmu) - s0(nm))
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
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) * wetu(nmu) == 1 .and. (m /= mf .or. slimp .or. .not.LMXF) .and. (m /= ml .or. srimp.or. .not.LMXL) ) then
                !
                fac = 0.5 * rdx / hkum(nm,k)
                !
                fac1 = zks(nmu,k-1) - zks(nm,k-1)
                fac2 = zks(nmu,k  ) - zks(nm,k  )
                !
                gmatu(nm,k,1) = (-hks(nm ,k) - fac1) * fac
                gmatu(nm,k,2) = (-hks(nm ,k) + fac2) * fac
                gmatu(nm,k,3) = ( hks(nmu,k) - fac1) * fac
                gmatu(nm,k,4) = ( hks(nmu,k) + fac2) * fac
                !
             else
                !
                gmatu(nm,k,1) = 0.
                gmatu(nm,k,2) = 0.
                gmatu(nm,k,3) = 0.
                gmatu(nm,k,4) = 0.
                !
             endif
             !
          enddo
          !
          if ( LMXF ) then
             !
             gmatu(nmf,k,3) = gmatu(nmf,k,3) - gmatu(nmf,k,1)
             gmatu(nmf,k,1) = 0.
             gmatu(nmf,k,4) = gmatu(nmf,k,4) - gmatu(nmf,k,2)
             gmatu(nmf,k,2) = 0.
             !
          endif
          !
          if ( LMXL ) then
             !
             gmatu(nml,k,1) = gmatu(nml,k,1) - gmatu(nml,k,3)
             gmatu(nml,k,3) = 0.
             gmatu(nml,k,2) = gmatu(nml,k,2) - gmatu(nml,k,4)
             gmatu(nml,k,4) = 0.
             !
          endif
          !
       enddo
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          gmatu(nm,1,1) = 0.
          gmatu(nm,1,3) = 0.
          !
       enddo
       !
       ! to reduce the pressure Poisson equation set pressure of bottom face to that of top face for a number of layers
       !
       do l = 1, qlay
          !
          do m = mf, ml
             !
             nm = kgrpnt(m)
             !
             gmatu(nm,kmax-l+1,1) = gmatu(nm,kmax-l+1,1) + gmatu(nm,kmax-l+1,2)
             gmatu(nm,kmax-l+1,2) = 0.
             gmatu(nm,kmax-l+1,3) = gmatu(nm,kmax-l+1,3) + gmatu(nm,kmax-l+1,4)
             gmatu(nm,kmax-l+1,4) = 0.
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
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) * wetu(nmu) == 1 .and. (m /= mf .or. slimp .or. .not.LMXF) .and. (m /= ml .or. srimp .or. .not.LMXL) ) then
                !
                fac1 = 0.5 * rdx * ( zks(nmu,k-1) - zks(nm,k-1) ) / ( hkum(nm,kd) + hkum(nm,k ) )
                fac2 = 0.5 * rdx * ( zks(nmu,k  ) - zks(nm,k  ) ) / ( hkum(nm,k ) + hkum(nm,ku) )
                !
                gmatu(nm,k,1) = -fac1
                gmatu(nm,k,2) = ( -rdx*hks(nm ,k) - hkum(nm,kd)*fac1 + hkum(nm,ku)*fac2 ) / hkum(nm,k)
                gmatu(nm,k,3) = -fac1
                gmatu(nm,k,4) = (  rdx*hks(nmu,k) - hkum(nm,kd)*fac1 + hkum(nm,ku)*fac2 ) / hkum(nm,k)
                gmatu(nm,k,5) =  fac2
                gmatu(nm,k,6) =  fac2
                !
             else
                !
                gmatu(nm,k,1) = 0.
                gmatu(nm,k,2) = 0.
                gmatu(nm,k,3) = 0.
                gmatu(nm,k,4) = 0.
                gmatu(nm,k,5) = 0.
                gmatu(nm,k,6) = 0.
                !
             endif
             !
          enddo
          !
          if ( LMXF ) then
             !
             gmatu(nmf,k,3) = gmatu(nmf,k,3) - gmatu(nmf,k,1)
             gmatu(nmf,k,1) = 0.
             gmatu(nmf,k,4) = gmatu(nmf,k,4) - gmatu(nmf,k,2)
             gmatu(nmf,k,2) = 0.
             gmatu(nmf,k,6) = gmatu(nmf,k,6) - gmatu(nmf,k,5)
             gmatu(nmf,k,5) = 0.
             !
          endif
          !
          if ( LMXL ) then
             !
             gmatu(nml,k,1) = gmatu(nml,k,1) - gmatu(nml,k,3)
             gmatu(nml,k,3) = 0.
             gmatu(nml,k,2) = gmatu(nml,k,2) - gmatu(nml,k,4)
             gmatu(nml,k,4) = 0.
             gmatu(nml,k,5) = gmatu(nml,k,5) - gmatu(nml,k,6)
             gmatu(nml,k,6) = 0.
             !
          endif
          !
       enddo
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          gmatu(nm,kmax,2) = gmatu(nm,kmax,2) + 2.*gmatu(nm,kmax,5)
          gmatu(nm,kmax,1) = gmatu(nm,kmax,1) -    gmatu(nm,kmax,5)
          gmatu(nm,kmax,4) = gmatu(nm,kmax,4) + 2.*gmatu(nm,kmax,6)
          gmatu(nm,kmax,3) = gmatu(nm,kmax,3) -    gmatu(nm,kmax,6)
          gmatu(nm,kmax,5) = 0.
          gmatu(nm,kmax,6) = 0.
          !
          gmatu(nm,1,2) = gmatu(nm,1,2) - gmatu(nm,1,1)
          gmatu(nm,1,4) = gmatu(nm,1,4) - gmatu(nm,1,3)
          gmatu(nm,1,1) = 0.
          gmatu(nm,1,3) = 0.
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
             do m = msta, mend
                !
                mu = m + 1
                !
                nm  = kgrpnt(m )
                nmu = kgrpnt(mu)
                !
                if ( wetu(nm) == 1 ) then
                   !
                   rhsu(nm,k) = rhsu(nm,k) - gmatu(nm,k,1)*q(nm ,kd) - gmatu(nm,k,2)*q(nm,k ) - gmatu(nm,k,3)*q(nmu,kd)  &
                                           - gmatu(nm,k,4)*q(nmu,k ) - gmatu(nm,k,5)*q(nm,ku) - gmatu(nm,k,6)*q(nmu,ku)
                   !
                endif
                !
             enddo
             !
          enddo
          !
       else if ( iproj == 2 .and. theta /= 1. ) then
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
             !
             do m = msta, mend
                !
                mu = m + 1
                !
                nm  = kgrpnt(m )
                nmu = kgrpnt(mu)
                !
                if ( wetu(nm) == 1 ) then
                   !
                   rhsu(nm,k) = rhsu(nm,k) - (1.-theta) * ( gmatu(nm,k,1)*q(nm ,kd) + gmatu(nm,k,2)*q(nm,k ) + gmatu(nm,k,3)*q(nmu,kd)  &
                                                          + gmatu(nm,k,4)*q(nmu,k ) + gmatu(nm,k,5)*q(nm,ku) + gmatu(nm,k,6)*q(nmu,ku) )
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
    ! compute baroclinic forcing at internal u-point
    !
    if ( idens /= 0 ) then
       !
       work(:,1) = 0.
       work(:,2) = 0.
       !
       do k = 1, kmax
          !
          do m = mf+1, ml-1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) - 0.5 * grav * hkum(nm,k) * rdx * (rho(nmu,k) - rho(nm,k)) / rhow
                !
             endif
             !
          enddo
          !
       enddo
       !
       do k = 2, kmax
          !
          do m = mf+1, ml-1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) == 1 ) then
                !
                rhou = rhow + 0.5 * (rho(nmu,k-1) + rho(nm,k-1))
                !
                work(nm,1) = work(nm,1) + rhou * (hks(nmu,k-1) - hks(nm,k-1)) + hkum(nm,k-1) * (rho(nmu,k-1) - rho(nm,k-1))
                !
                work(nm,2) = work(nm,2) + hks(nmu,k-1) - hks(nm,k-1)
                !
                rhou = rhow + 0.5 * (rho(nmu,k) + rho(nm,k))
                !
                rhsu(nm,k) = rhsu(nm,k) - grav * rdx * ( work(nm,1) - rhou*work(nm,2) ) / rhow
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! compute atmospheric pressure gradient at internal u-point
    !
    if ( svwp ) then
       !
       do m = mf+1, ml-1
          !
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmu = kgrpnt(mu)
          !
          if ( wetu(nm) == 1 ) then
             !
             rhsu(nm,:) = rhsu(nm,:) - rdx * (patm(nmu) - patm(nm)) / rhow
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute horizontal viscosity term at internal u-point
    !
    if ( ihvisc == 1 .and. hvisc > 0. ) then
       !
       do k = 1, kmax
          !
          do m = mf+1, ml-1
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) * wetu(nmd) * wetu(nmu) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) + hvisc * rdx * rdx * ( hks(nmu,k)*(u0(nmu,k) - u0(nm,k)) - hks(nm,k)*(u0(nm,k) - u0(nmd,k)) ) / hkum(nm,k)
                !
             endif
             !
          enddo
          !
       enddo
       !
    else if ( ihvisc > 1 ) then
       !
       do k = 1, kmax
          !
          do m = mf+1, ml-1
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) * wetu(nmd) * wetu(nmu) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) + rdx * rdx * ( vnu2d(nmu) * hks(nmu,k)*(u0(nmu,k) - u0(nm,k)) - vnu2d(nm) * hks(nm,k)*(u0(nm,k) - u0(nmd,k)) ) / hkum(nm,k)
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! compute divergence of Reynolds stress tensor at internal u-point
    !
    if ( iturb > 1 ) then
       !
       do k = 1, kmax
          !
          do m = mf+1, ml-1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) + rdx * ( rsuu(nmu,k) - rsuu(nm,k) ) + ( rsuw(nm,k-1) - rsuw(nm,k) ) / hkum(nm,k)
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! compute explicit part of wind stress term at internal u-point (top layer only), if appropriate
    !
    if ( iwind /= 0 ) then
       !
       do m = mf+1, ml-1
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             rhsu(nm,1) = rhsu(nm,1) + windu(nm)/max(1.e-3,hkum(nm,1))
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute implicit part of wind stress term at internal u-point (top layer only), if appropriate
    !
    if ( relwnd ) then
       !
       do m = mf+1, ml-1
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             amatu(nm,1,1) = amatu(nm,1,1) + cwndu(nm) / max(1.e-3,hkum(nm,1))
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute bottom friction (bottom layer only), if appropriate
    !
    if ( irough == 4 ) then
       !
       ! logarithmic wall-law
       !
       do m = msta, mend
          !
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmu = kgrpnt(mu)
          !
          if ( wetu(nm) == 1 ) then
             !
             amatu(nm,kmax,1) = amatu(nm,kmax,1) + 0.5 * ( logfrc(nm,1) + logfrc(nmu,1) ) / hkum(nm,kmax)
             !
          endif
          !
       enddo
       !
    else if ( irough == 11 ) then
       !
       ! linear bottom friction
       !
       do m = msta, mend
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             amatu(nm,kmax,1) = amatu(nm,kmax,1) + cfricu(nm) / hkum(nm,kmax)
             !
          endif
          !
       enddo
       !
    else if ( irough /= 0 ) then
       !
       do m = msta, mend
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 .and. abs(u0(nm,kmax)) > 1.e-8 ) then
             !
             amatu(nm,kmax,1) = amatu(nm,kmax,1) + cfricu(nm) * udep(nm) * udep(nm) / ( abs(u0(nm,kmax)) * hkum(nm,kmax) )
             !
          endif
          !
       enddo
       !
    endif
    !
    ! include effect of ambient current - first order approximation
    !
    if ( iamb == 1 ) then
       !
       ! compute the advective velocity in wl-point based on first order upwind
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             if ( Ubck(nm,k) > 1.0e-5 ) then
                !
                ua(nm,k) = u0(nmd,k)
                !
             else if ( Ubck(nm,k) < -1.0e-5 ) then
                !
                ua(nm,k) = u0(nm ,k)
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! compute horizontal advection term (momentum conservative) at internal u-point (prediction)
       !
       do k = 1, kmax
          !
          do m = mf+1, ml-1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) - rdx * ( Ubck(nmu,k) * ua(nmu,k) - Ubck(nm,k) * ua(nm,k) )
                !
             endif
             !
          enddo
          !
       enddo
       !
    else if ( iamb == 2 ) then
       !
       ! compute the advective velocity in wl-point based on first order upwind
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             fac = 0.5 * ( Ubck(nm,k) + Ubck(nmd,k) )
             !
             if ( fac > 1.0e-5 ) then
                !
                ua(nm,k) = u0(nmd,k)
                !
             else if ( fac < -1.0e-5 ) then
                !
                ua(nm,k) = u0(nm ,k)
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! compute horizontal advection term (momentum conservative) at internal u-point (prediction)
       !
       do k = 1, kmax
          !
          do m = mf+1, ml-1
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             nmu = kgrpnt(mu)
             !
             fac1 = 0.5 * ( Ubck(nmu,k) + Ubck(nm ,k) )
             fac2 = 0.5 * ( Ubck(nm ,k) + Ubck(nmd,k) )
             !
             if ( wetu(nm) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) - rdx * ( fac1 * ua(nmu,k) - fac2 * ua(nm,k) )
                !
             endif
             !
          enddo
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
       do m = msta, mend
          !
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmu = kgrpnt(mu)
          !
          if ( wetu(nm) == 1 ) then
             !
             ! advection term
             !
             if ( propsc == 3 .and. kappa == 1. ) then
                !
                ! central differences
                !
                ctrkt = 0.5 * kwd * ( wom(nm,k-1) + wom(nmu,k-1) ) / ( hkum(nm,kd) + hkum(nm,k ) )
                ctrkb = 0.5 * kwu * ( wom(nm,k  ) + wom(nmu,k  ) ) / ( hkum(nm,k ) + hkum(nm,ku) )
                !
                amatu(nm,k,1) = amatu(nm,k,1) - ctrkt + ctrkb
                amatu(nm,k,2) =  ctrkt
                amatu(nm,k,3) = -ctrkb
                !
             else
                !
                ! first order upwind scheme
                !
                ctrkt = 0.5 * kwd * ( wom(nm,k-1) + wom(nmu,k-1) ) / hkum(nm,k)
                ctrkb = 0.5 * kwu * ( wom(nm,k  ) + wom(nmu,k  ) ) / hkum(nm,k)
                !
                amatu(nm,k,1) = amatu(nm,k,1) - ctrkt + ctrkb + max(ctrkt,0.) - min(ctrkb,0.)
                amatu(nm,k,2) =  min(ctrkt,0.)
                amatu(nm,k,3) = -max(ctrkb,0.)
                !
                ! add higher order (flux-limited) correction, if appropriate
                !
                if ( propsc /= 1 ) then
                   !
                   if ( ctrkt > 0. ) then
                      !
                      grad1 = u0(nm,kd) - u0(nm,k )
                      grad2 = u0(nm,k ) - u0(nm,ku)
                      !
                      rhsu(nm,k) = rhsu(nm,k) - 0.5 * ctrkt * fluxlim(grad1,grad2)
                      !
                   else if ( ctrkt < 0. ) then
                      !
                      grad1 = u0(nm,kd ) - u0(nm,k )
                      grad2 = u0(nm,kdd) - u0(nm,kd)
                      !
                      rhsu(nm,k) = rhsu(nm,k) + 0.5 * ctrkt * fluxlim(grad1,grad2)
                      !
                   endif
                   !
                   if ( ctrkb > 0. ) then
                      !
                      grad1 = u0(nm,k ) - u0(nm,ku )
                      grad2 = u0(nm,ku) - u0(nm,kuu)
                      !
                      rhsu(nm,k) = rhsu(nm,k) + 0.5 * ctrkb * fluxlim(grad1,grad2)
                      !
                   else if ( ctrkb < 0. ) then
                      !
                      grad1 = u0(nm,k ) - u0(nm,ku)
                      grad2 = u0(nm,kd) - u0(nm,k )
                      !
                      rhsu(nm,k) = rhsu(nm,k) - 0.5 * ctrkb * fluxlim(grad1,grad2)
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
                ctrkt = ( vnu3d(nm,k-1) + vnu3d(nmu,k-1) ) / ( hkum(nm,k)*( hkum(nm,kd) + hkum(nm,k ) ) )
                ctrkb = ( vnu3d(nm,k  ) + vnu3d(nmu,k  ) ) / ( hkum(nm,k)*( hkum(nm,k ) + hkum(nm,ku) ) )
                !
                amatu(nm,k,1) = amatu(nm,k,1) + kwd*ctrkt + (2.*kwu-1.)*ctrkb
                amatu(nm,k,2) = amatu(nm,k,2) - kwd*ctrkt +    (1.-kwu)*ctrkb
                amatu(nm,k,3) = amatu(nm,k,3) -                     kwu*ctrkb
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
    do m = msta, mend
       !
       nm = kgrpnt(m)
       !
       bi = 1./amatu(nm,1,1)
       !
       amatu(nm,1,1) = bi
       amatu(nm,1,3) = amatu(nm,1,3)*bi
       rhsu (nm,1  ) = rhsu (nm,1  )*bi
       !
       do k = 2, kmax
          !
          bi = 1./(amatu(nm,k,1) - amatu(nm,k,2)*amatu(nm,k-1,3))
          amatu(nm,k,1) = bi
          amatu(nm,k,3) = amatu(nm,k,3)*bi
          rhsu (nm,k  ) = (rhsu(nm,k) - amatu(nm,k,2)*rhsu(nm,k-1))*bi
          !
       enddo
       !
       u1(nm,kmax) = rhsu(nm,kmax)
       do k = kmax-1, 1, -1
          u1(nm,k) = rhsu(nm,k) - amatu(nm,k,3)*u1(nm,k+1)
       enddo
       !
    enddo
    !
    ! re-update the solution in case of thetau <> 1
    !
    if ( thetau /= 1. ) then
       !
       do k = 1, kmax
          !
          do m = msta, mend
             !
             nm = kgrpnt(m)
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm,k) = ( u1(nm,k) - (1.-thetau) * u0(nm,k) ) / thetau
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! exchange u-velocities with neighbouring subdomains
    !
    call SWEXCHG ( u1, kgrpnt, 1, kmax )
    if (STPNOW()) return
    !
    ! compute higher order correction to the advective velocity in wl-point (if appropriate)
    !
    propsc = nint(pnums(6))
    !
    if ( propsc /= 1 ) then
       !
       kappa  = pnums(7)
       mbound = pnums(8)
       phieby = pnums(9)
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             !
             md  = m  - 1
             mdd = md - 1
             mu  = m  + 1
             if ( LMXF .and. mdd < mf ) mdd = mf
             if ( LMXL .and. mu  > ml ) mu  = ml
             !
             nm   = kgrpnt(m  )
             nmd  = kgrpnt(md )
             nmdd = kgrpnt(mdd)
             nmu  = kgrpnt(mu )
             !
             if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac = 0.5 * ( u0(nm,k) + u0(nmd,k) )
                !
             else
                !
                fac = qm(nm,k)
                !
             endif
             !
             if ( fac > 1.0e-5 ) then
                !
                grad1 = u1(nm ,k) - u0(nmd ,k)
                grad2 = u1(nmd,k) - u0(nmdd,k)
                !
                up(nm,k) = u0(nmd,k)
                ua(nm,k) = 0.5 * fluxlim(grad1,grad2)
                !
             else if ( fac < -1.0e-5 ) then
                !
                grad1 = u1(nmd,k) - u0(nm ,k)
                grad2 = u1(nm ,k) - u0(nmu,k)
                !
                up(nm,k) = u0(nm,k)
                ua(nm,k) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! correct u-velocity associated with higher order approximation of advection term at internal u-point (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do k = 1, kmax
          !
          do m = mf+1, ml-1
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             nmu = kgrpnt(mu)
             !
             if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                !fac1 = up(nmu,k) + 0.5*ua(nmu,k)
                !fac2 = up(nm ,k) + 0.5*ua(nm ,k)
                !jumps back to old implementation (020117)
                fac1 = 0.
                fac2 = 0.
                !
             else if ( hkumn(nm,k) > 0. ) then
                !
                fac1 = qm(nmu,k) / hkumn(nm,k)
                fac2 = qm(nm ,k) / hkumn(nm,k)
                !
             else
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm,k) = u1(nm,k) - dt * rdx * ( fac1 * ua(nmu,k)  - fac2 * ua(nm,k) )
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! exchange u-velocities with neighbouring subdomains (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       call SWEXCHG ( u1, kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
    endif
    !
    ! include effect of ambient current - high order correction
    !
    propsc = nint(pnums(61))
    !
    if ( iamb == 1 .and. propsc /= 1 ) then
       !
       ! compute higher order correction to the advective velocity in wl-point (if appropriate)
       !
       kappa  = pnums(62)
       mbound = pnums(63)
       phieby = pnums(64)
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             !
             md  = m  - 1
             mdd = md - 1
             mu  = m  + 1
             if ( LMXF .and. mdd < mf ) mdd = mf
             if ( LMXL .and. mu  > ml ) mu  = ml
             !
             nm   = kgrpnt(m  )
             nmd  = kgrpnt(md )
             nmdd = kgrpnt(mdd)
             nmu  = kgrpnt(mu )
             !
             if ( Ubck(nm,k) > 1.0e-5 ) then
                !
                grad1 = u1(nm ,k) - u0(nmd ,k)
                grad2 = u1(nmd,k) - u0(nmdd,k)
                !
                ua(nm,k) = 0.5 * fluxlim(grad1,grad2)
                !
             else if ( Ubck(nm,k) < -1.0e-5 ) then
                !
                grad1 = u1(nmd,k) - u0(nm ,k)
                grad2 = u1(nm ,k) - u0(nmu,k)
                !
                ua(nm,k) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! correct u-velocity associated with higher order approximation of advection term at internal u-point (if appropriate)
       !
       do k = 1, kmax
          !
          do m = mf+1, ml-1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm,k) = u1(nm,k) - dt * rdx * ( Ubck(nmu,k) * ua(nmu,k)  - Ubck(nm,k) * ua(nm,k) )
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! exchange u-velocities with neighbouring subdomains (if appropriate)
       !
       call SWEXCHG ( u1, kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
    else if ( iamb == 2 .and. propsc /= 1 ) then
       !
       ! compute higher order correction to the advective velocity in wl-point (if appropriate)
       !
       kappa  = pnums(62)
       mbound = pnums(63)
       phieby = pnums(64)
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             !
             md  = m  - 1
             mdd = md - 1
             mu  = m  + 1
             if ( LMXF .and. mdd < mf ) mdd = mf
             if ( LMXL .and. mu  > ml ) mu  = ml
             !
             nm   = kgrpnt(m  )
             nmd  = kgrpnt(md )
             nmdd = kgrpnt(mdd)
             nmu  = kgrpnt(mu )
             !
             fac = 0.5 * ( Ubck(nm,k) + Ubck(nmd,k) )
             !
             if ( fac > 1.0e-5 ) then
                !
                grad1 = u1(nm ,k) - u0(nmd ,k)
                grad2 = u1(nmd,k) - u0(nmdd,k)
                !
                ua(nm,k) = 0.5 * fluxlim(grad1,grad2)
                !
             else if ( fac < -1.0e-5 ) then
                !
                grad1 = u1(nmd,k) - u0(nm ,k)
                grad2 = u1(nm ,k) - u0(nmu,k)
                !
                ua(nm,k) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! correct u-velocity associated with higher order approximation of advection term at internal u-point (if appropriate)
       !
       do k = 1, kmax
          !
          do m = mf+1, ml-1
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             nmu = kgrpnt(mu)
             !
             fac1 = 0.5 * ( Ubck(nmu,k) + Ubck(nm ,k) )
             fac2 = 0.5 * ( Ubck(nm ,k) + Ubck(nmd,k) )
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm,k) = u1(nm,k) - dt * rdx * ( fac1 * ua(nmu,k)  - fac2 * ua(nm,k) )
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! exchange u-velocities with neighbouring subdomains (if appropriate)
       !
       call SWEXCHG ( u1, kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
    endif
    !
    ! build the w-momentum equation
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       ! initialize system of equations in dry points
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m)
          !
          if ( wets(nm) /= 1 ) then
             !
             amatw(nm,:,1) = 1.
             amatw(nm,:,2) = 0.
             amatw(nm,:,3) = 0.
             rhsw (nm,:  ) = 0.
             !
          endif
          !
       enddo
       !
       ! bottom:
       !
       ! the kinematic condition is imposed
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          if ( wets(nm) == 1 ) then
             !
             amatw(nm,kmax,1) = 1.
             amatw(nm,kmax,2) = 0.
             amatw(nm,kmax,3) = 0.
             w1   (nm,kmax  ) = 0.25 * rdx * ( 3.*(u1(nm,kmax)+u1(nmd,kmax)) - u1(nm,kmax-1) - u1(nmd,kmax-1) ) * ( zku(nm,kmax) - zku(nmd,kmax) )
             rhsw (nm,kmax  ) = w1(nm,kmax)
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
          do m = mfu, ml
             !
             nm = kgrpnt(m)
             !
             if ( wets(nm) == 1 ) then
                !
                amatw(nm,k,1) = 1. / (dt*thetaw)
                rhsw (nm,k  ) = w0(nm,k) / (dt*thetaw)
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
             do m = mfu, ml
                !
                md  = m  - 1
                mdd = md - 1
                mu  = m  + 1
                muu = mu + 1
                if ( LMXF .and. mdd < mf ) mdd = mf
                if ( LMXL .and. muu > ml ) muu = ml
                !
                nm   = kgrpnt(m  )
                nmd  = kgrpnt(md )
                nmdd = kgrpnt(mdd)
                nmu  = kgrpnt(mu )
                nmuu = kgrpnt(muu)
                !
                if ( wets(nm) == 1 ) then
                   !
                   fac = 0.5 * dx * ( hks(nm,kd) + hks(nm,k+1) )
                   !
                   if ( wetu(nmd) == 1 ) then
                      fac1 = ( qx(nmd,k+1)*hku(nmd,kd) + qx(nmd,kd)*hku(nmd,k+1) ) / ( fac*( hku(nmd,kd) + hku(nmd,k+1) ) )
                   else
                      fac1 = 0.
                   endif
                   if ( wetu(nm) == 1 ) then
                      fac2 = ( qx(nm ,k+1)*hku(nm ,kd) + qx(nm ,kd)*hku(nm ,k+1) ) / ( fac*( hku(nm ,kd) + hku(nm ,k+1) ) )
                   else
                      fac2 = 0.
                   endif
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      rhsw(nm,k) = rhsw(nm,k) - 0.5 * ( fac1 * (w0(nm,k) - w0(nmd,k)) + fac2 * (w0(nmu,k) - w0(nm,k)) )
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      rhsw(nm,k) = rhsw(nm,k) - ( max(fac1,0.) * (w0(nm,k) - w0(nmd,k)) + min(fac2,0.) * (w0(nmu,k) - w0(nm,k)) )
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         if ( fac2 > 0. ) then
                            !
                            grad1 = w0(nmu,k) - w0(nm ,k)
                            grad2 = w0(nm ,k) - w0(nmd,k)
                            !
                            rhsw(nm,k) = rhsw(nm,k) - 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         else if ( fac2 < 0. ) then
                            !
                            grad1 = w0(nmu ,k) - w0(nm ,k)
                            grad2 = w0(nmuu,k) - w0(nmu,k)
                            !
                            rhsw(nm,k) = rhsw(nm,k) + 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                         if ( fac1 > 0. ) then
                            !
                            grad1 = w0(nm ,k) - w0(nmd ,k)
                            grad2 = w0(nmd,k) - w0(nmdd,k)
                            !
                            rhsw(nm,k) = rhsw(nm,k) + 0.5 * fac1 * fluxlim(grad1,grad2)
                            !
                         else if ( fac1 < 0. ) then
                            !
                            grad1 = w0(nm ,k) - w0(nmd,k)
                            grad2 = w0(nmu,k) - w0(nm ,k)
                            !
                            rhsw(nm,k) = rhsw(nm,k) - 0.5 * fac1 * fluxlim(grad1,grad2)
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
          ! compute horizontal viscosity term
          !
          if ( ihvisc == 1 .and. hvisc > 0. ) then
             !
             do k = 0, kmax-1
                !
                kd = max(k,1)
                !
                do m = mfu, ml
                   !
                   md = m - 1
                   mu = m + 1
                   !
                   nm  = kgrpnt(m )
                   nmd = kgrpnt(md)
                   nmu = kgrpnt(mu)
                   !
                   if ( wets(nm) * wets(nmd) * wets(nmu) == 1 ) then
                      !
                      fac1 = hvisc * ( hkum(nm ,kd) + hkum(nm ,k+1) )
                      !
                      fac2 = hvisc * ( hkum(nmd,kd) + hkum(nmd,k+1) )
                      !
                      rhsw(nm,k) = rhsw(nm,k) + rdx * rdx * ( fac1 * (w0(nmu,k) - w0(nm,k)) - fac2 * (w0(nm,k) - w0(nmd,k)) ) / ( hks(nm,kd) + hks(nm,k+1) )
                      !
                   endif
                   !
                enddo
                !
             enddo
             !
          else if ( ihvisc > 1 ) then
             !
             do k = 0, kmax-1
                !
                kd = max(k,1)
                !
                do m = mfu, ml
                   !
                   md = m - 1
                   mu = m + 1
                   !
                   nm  = kgrpnt(m )
                   nmd = kgrpnt(md)
                   nmu = kgrpnt(mu)
                   !
                   if ( wets(nm) * wets(nmd) * wets(nmu) == 1 ) then
                      !
                      fac1 = 0.5 * ( hkum(nm ,kd) + hkum(nm ,k+1) ) * ( vnu2d(nm ) + vnu2d(nmu) )
                      !
                      fac2 = 0.5 * ( hkum(nmd,kd) + hkum(nmd,k+1) ) * ( vnu2d(nmd) + vnu2d(nm ) )
                      !
                      rhsw(nm,k) = rhsw(nm,k) + rdx * rdx * ( fac1 * (w0(nmu,k) - w0(nm,k)) - fac2 * (w0(nm,k) - w0(nmd,k)) ) / ( hks(nm,kd) + hks(nm,k+1) )
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
       if ( iamb /= 0 ) then
          !
          ! compute horizontal advection term due to ambient current (momentum conservative and explicit)
          !
          propsc = nint(pnums(66))
          kappa  = pnums(67)
          mbound = pnums(68)
          phieby = pnums(69)
          !
          do k = 0, kmax-1
             !
             kd = max(k,1)
             !
             do m = mfu, ml
                !
                md  = m  - 1
                mdd = md - 1
                mu  = m  + 1
                muu = mu + 1
                if ( LMXF .and. mdd < mf ) mdd = mf
                if ( LMXL .and. muu > ml ) muu = ml
                !
                nm   = kgrpnt(m  )
                nmd  = kgrpnt(md )
                nmdd = kgrpnt(mdd)
                nmu  = kgrpnt(mu )
                nmuu = kgrpnt(muu)
                !
                if ( wets(nm) == 1 ) then
                   !
                   fac = 0.5 * dx * ( hks(nm,kd) + hks(nm,k+1) )
                   !
                   if ( iamb == 1 ) then
                      fac1 = 0.25 * ( hku(nmd,kd)*(Ubck(nm ,kd) + Ubck(nmd,kd)) + hku(nmd,k+1)*(Ubck(nm ,k+1) + Ubck(nmd,k+1)) )
                      fac2 = 0.25 * ( hku(nm ,kd)*(Ubck(nmu,kd) + Ubck(nm ,kd)) + hku(nm ,k+1)*(Ubck(nmu,k+1) + Ubck(nm ,k+1)) )
                   else if ( iamb == 2 ) then
                      fac1 = 0.5 * ( hku(nmd,kd)*Ubck(nmd,kd) + hku(nmd,k+1)*Ubck(nmd,k+1) )
                      fac2 = 0.5 * ( hku(nm ,kd)*Ubck(nm ,kd) + hku(nm ,k+1)*Ubck(nm ,k+1) )
                   endif
                   fac1 = fac1 / fac
                   fac2 = fac2 / fac
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      rhsw(nm,k) = rhsw(nm,k) - 0.5 * ( fac1 * (w0(nm,k) - w0(nmd,k)) + fac2 * (w0(nmu,k) - w0(nm,k)) )
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      rhsw(nm,k) = rhsw(nm,k) - ( max(fac1,0.) * (w0(nm,k) - w0(nmd,k)) + min(fac2,0.) * (w0(nmu,k) - w0(nm,k)) )
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         if ( fac2 > 0. ) then
                            !
                            grad1 = w0(nmu,k) - w0(nm ,k)
                            grad2 = w0(nm ,k) - w0(nmd,k)
                            !
                            rhsw(nm,k) = rhsw(nm,k) - 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         else if ( fac2 < 0. ) then
                            !
                            grad1 = w0(nmu ,k) - w0(nm ,k)
                            grad2 = w0(nmuu,k) - w0(nmu,k)
                            !
                            rhsw(nm,k) = rhsw(nm,k) + 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                         if ( fac1 > 0. ) then
                            !
                            grad1 = w0(nm ,k) - w0(nmd ,k)
                            grad2 = w0(nmd,k) - w0(nmdd,k)
                            !
                            rhsw(nm,k) = rhsw(nm,k) + 0.5 * fac1 * fluxlim(grad1,grad2)
                            !
                         else if ( fac1 < 0. ) then
                            !
                            grad1 = w0(nm ,k) - w0(nmd,k)
                            grad2 = w0(nmu,k) - w0(nm ,k)
                            !
                            rhsw(nm,k) = rhsw(nm,k) - 0.5 * fac1 * fluxlim(grad1,grad2)
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
             do m = mfu, ml
                !
                nm = kgrpnt(m)
                !
                if ( wets(nm) == 1 ) then
                   !
                   ! advection term
                   !
                   ctrkt = ( wom(nm,kd-1) + wom(nm,kd  ) ) / ( hks(nm,kd) + hks(nm,k+1) )
                   ctrkb = ( wom(nm,k   ) + wom(nm,k +1) ) / ( hks(nm,kd) + hks(nm,k+1) )
                   !
                   if ( k == 0 ) ctrkt = 0.
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      amatw(nm,k,1) = amatw(nm,k,1) - 0.5*ctrkt + 0.5*ctrkb
                      amatw(nm,k,2) =  0.5*ctrkt
                      amatw(nm,k,3) = -0.5*ctrkb
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      amatw(nm,k,1) = amatw(nm,k,1) - ctrkt + ctrkb + max(ctrkt,0.) - min(ctrkb,0.)
                      amatw(nm,k,2) =  min(ctrkt,0.)
                      amatw(nm,k,3) = -max(ctrkb,0.)
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         if ( ctrkt > 0. ) then
                            !
                            grad1 = w0(nm,kd-1) - w0(nm,k  )
                            grad2 = w0(nm,k   ) - w0(nm,k+1)
                            !
                            rhsw(nm,k) = rhsw(nm,k) - 0.5 * ctrkt * fluxlim(grad1,grad2)
                            !
                         else if ( ctrkt < 0. ) then
                            !
                            grad1 = w0(nm,kd -1) - w0(nm,k   )
                            grad2 = w0(nm,kdd-1) - w0(nm,kd-1)
                            !
                            rhsw(nm,k) = rhsw(nm,k) + 0.5 * ctrkt * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                         if ( ctrkb > 0. ) then
                            !
                            grad1 = w0(nm,k  ) - w0(nm,k  +1)
                            grad2 = w0(nm,k+1) - w0(nm,kuu  )
                            !
                            rhsw(nm,k) = rhsw(nm,k) + 0.5 * ctrkb * fluxlim(grad1,grad2)
                            !
                         else if ( ctrkb < 0. ) then
                            !
                            grad1 = w0(nm,k   ) - w0(nm,k+1)
                            grad2 = w0(nm,kd-1) - w0(nm,k  )
                            !
                            rhsw(nm,k) = rhsw(nm,k) - 0.5 * ctrkb * fluxlim(grad1,grad2)
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
                      ctrkt = ( vnu3d(nm,kd-1) + vnu3d(nm,kd  ) ) / ( hks(nm,kd  )*( hks(nm,kd) + hks(nm,k+1) ) )
                      ctrkb = ( vnu3d(nm,k   ) + vnu3d(nm,k +1) ) / ( hks(nm,k +1)*( hks(nm,kd) + hks(nm,k+1) ) )
                      !
                      if ( k == 0 ) ctrkt = 0.
                      !
                      amatw(nm,k,1) = amatw(nm,k,1) + ctrkt + ctrkb
                      amatw(nm,k,2) = amatw(nm,k,2) - ctrkt
                      amatw(nm,k,3) = amatw(nm,k,3) - ctrkb
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
          do m = mfu, ml
             !
             nm = kgrpnt(m)
             !
             if ( wets(nm) == 1 ) then
                !
                gmatw(nm,:,1) =  2./hks(nm,:)
                gmatw(nm,:,2) = -gmatw(nm,:,1)
                !
                do k = 2, kmax
                   !
                   fac = 1.
                   !
                   do j = 1, kmax-k+1
                      !
                      fac = -fac
                      !
                      gmatw(nm,k-1,2*j+1) =  2.*fac/hks(nm,j+k-1)
                      gmatw(nm,k-1,2*j+2) = -gmatw(nm,k-1,2*j+1)
                      !
                   enddo
                   !
                enddo
                !
             else
                !
                gmatw(nm,:,:) = 0.
                !
             endif
             !
             gmatw(nm,1,1) = 0.
             !
          enddo
          !
          ! to reduce the pressure Poisson equation set pressure of bottom face to that of top face for a number of layers
          !
          do l = 1, qlay
             !
             do m = mfu, ml
                !
                nm = kgrpnt(m)
                !
                do k = 1, kmax
                   !
                   j = kmax +1 - k - l
                   if ( j < 0 ) cycle
                   !
                   gmatw(nm,k,2*j+1) = gmatw(nm,k,2*j+1) + gmatw(nm,k,2*j+2)
                   gmatw(nm,k,2*j+2) = 0.
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
          do m = mfu, ml
             !
             nm = kgrpnt(m)
             !
             if ( wets(nm) == 1 ) then
                !
                do k = 1, kmax
                   !
                   kd = max(k-1,1)
                   !
                   gmatw(nm,k,1) =  2./(hks(nm,kd) + hks(nm,k))
                   gmatw(nm,k,2) = -gmatw(nm,k,1)
                   !
                enddo
                !
             else
                !
                gmatw(nm,:,1) = 0.
                gmatw(nm,:,2) = 0.
                !
             endif
             !
             gmatw(nm,1,2) = gmatw(nm,1,2) - gmatw(nm,1,1)
             gmatw(nm,1,1) = 0.
             !
          enddo
          !
       endif
       !
       ! compute gradient of non-hydrostatic pressure
       !
       if ( linswe ) then
          fac = (-1.)**kmax
       else
          fac = 0.
       endif
       !
       if ( iproj == 1 ) then
          fac1 = 1.
       else if ( iproj == 2 ) then
          fac1 = 1. - theta
       endif
       !
       do k = 0, kmax-1
          !
          kd = max(k,1)
          !
          fac = -fac
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m)
             !
             if ( wets(nm) == 1 ) then
                !
                if ( ihydro == 1 ) then
                   !
                   do j = 0, kmax-1
                      !
                      kd = max(k+j  ,1   )
                      kd = min(kd   ,kmax)
                      ku = min(k+j+1,kmax)
                      !
                      rhsw(nm,k) = rhsw(nm,k) - fac1 * ( gmatw(nm,k+1,2*j+1)*q(nm,kd) + gmatw(nm,k+1,2*j+2)*q(nm,ku) )
                      !
                   enddo
                   !
                   rhsw(nm,k) = rhsw(nm,k) - fac * ( w1(nm,kmax) - w0(nm,kmax) ) / (dt*thetaw)
                   !
                else
                   !
                   rhsw(nm,k) = rhsw(nm,k) - fac1 * ( gmatw(nm,k+1,1)*q(nm,kd) + gmatw(nm,k+1,2)*q(nm,k+1) )
                   !
                endif
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! compute flow resistance inside porous medium, if appropriate
       !
       if ( iporos == 1 ) then
          !
          do k = 0, kmax-1
             !
             kd = max(k,1)
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                !
                if ( wets(nm) == 1 ) then
                   !
                   u = 0.25 * ( u0(nm,kd) + u0(nmd,kd) + u0(nm,k+1) + u0(nmd,k+1) )
                   !
                   utot = sqrt( u*u + w0(nm,k)*w0(nm,k) )
                   !
                   amatw(nm,k,1) = amatw(nm,k,1) + apoks(nm,k) + bpoks(nm,k) * utot
                   !
                   amatw(nm,k,1) = amatw(nm,k,1) + cpoks(nm,k) / (dt*thetaw)
                   rhsw (nm,k  ) = rhsw (nm,k  ) + cpoks(nm,k) * w0(nm,k) / (dt*thetaw)
                   !
                endif
                !
             enddo
             !
          enddo
          !
       endif
       !
       ! compute friction due to vegetation modelled as horizontal cylinders, if appropriate
       !
       if ( ivegw /= 0 ) then
          !
          do k = 0, kmax-1
             !
             kd = max(k,1)
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                !
                if ( wets(nm) == 1 ) then
                   !
                   amatw(nm,k,1) = amatw(nm,k,1) + 0.25 * ( cvegu(nm,kd,1) + cvegu(nmd,kd,1) + cvegu(nm,k+1,1) + cvegu(nmd,k+1,1) ) * abs(w0(nm,k))
                   !
                endif
                !
             enddo
             !
          enddo
          !
       endif
       !
       ! compute divergence of Reynolds stress tensor, if appropriate
       !
       if ( iturb > 1 ) then
          !
          do k = 1, kmax-1
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                !
                if ( wets(nm) == 1 ) then
                   !
                   rhsw(nm,k) = rhsw(nm,k) + rdx * ( rswu(nm,k) - rswu(nmd,k) ) + 2. * ( rsww(nm,k) - rsww(nm,k+1) ) / ( hks(nm,k) + hks(nm,k+1) )
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
    ! solve the w-momentum equation
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m)
          !
          bi = 1./amatw(nm,0,1)
          !
          amatw(nm,0,1) = bi
          amatw(nm,0,3) = amatw(nm,0,3)*bi
          rhsw (nm,0  ) = rhsw (nm,0  )*bi
          !
          do k = 1, kmax
             !
             bi = 1./(amatw(nm,k,1) - amatw(nm,k,2)*amatw(nm,k-1,3))
             amatw(nm,k,1) = bi
             amatw(nm,k,3) = amatw(nm,k,3)*bi
             rhsw (nm,k  ) = (rhsw(nm,k) - amatw(nm,k,2)*rhsw(nm,k-1))*bi
             !
          enddo
          !
          w1(nm,kmax) = rhsw(nm,kmax)
          do k = kmax-1, 0, -1
             w1(nm,k) = rhsw(nm,k) - amatw(nm,k,3)*w1(nm,k+1)
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
             do m = mfu, ml
                !
                nm = kgrpnt(m)
                !
                if ( wets(nm) == 1 ) then
                   !
                   w1(nm,k) = ( w1(nm,k) - (1.-thetaw) * w0(nm,k) ) / thetaw
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
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          if ( wets(nm) == 1 ) then
             !
             do k = 1, kmax
                !
                kd = max(k-1,1   )
                ku = min(k+1,kmax)
                !
                fac1 = 0.5 * rdx * ( zku(nm,k-1) - zku(nmd,k-1) ) / ( hks(nm,k ) + hks(nm,kd) )
                fac2 = 0.5 * rdx * ( zku(nm,k  ) - zku(nmd,k  ) ) / ( hks(nm,ku) + hks(nm,k ) )
                !
                if ( k == kmax ) fac2 = 0.
                !
                dmat(nm,k,1) = -fac1 * hks(nm ,k)
                dmat(nm,k,2) =   rdx * hku(nm ,k) - fac1 * hks(nm,kd) + fac2 * hks(nm,ku)
                dmat(nm,k,3) = -fac1 * hks(nm ,k)
                dmat(nm,k,4) =  -rdx * hku(nmd,k) - fac1 * hks(nm,kd) + fac2 * hks(nm,ku)
                dmat(nm,k,5) =  fac2 * hks(nm ,k)
                dmat(nm,k,6) =  fac2 * hks(nm ,k)
                !
             enddo
             !
          else
             !
             dmat(nm,:,:) = 0.
             !
          endif
          !
       enddo
       !
       do m = mfu, ml
          !
          nm  = kgrpnt(m)
          !
          dmat(nm,1,2) = dmat(nm,1,2) + 2.*dmat(nm,1,1)
          dmat(nm,1,5) = dmat(nm,1,5) -    dmat(nm,1,1)
          dmat(nm,1,4) = dmat(nm,1,4) + 2.*dmat(nm,1,3)
          dmat(nm,1,6) = dmat(nm,1,6) -    dmat(nm,1,3)
          dmat(nm,1,1) = 0.
          dmat(nm,1,3) = 0.
          !
       enddo
       !
       ! build the Poisson equation
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
             !
             kwu = 1.
             if ( k == kmax ) kwu = 0.
             !
             amatp(nm,k, 1) = dmat(nm,k,2) * gmatu(nm ,k ,2) + dmat(nm,k,4) * gmatu(nmd,k ,4) + dmat(nm,k,1) * gmatu(nm ,kd,5) +  &
                              dmat(nm,k,3) * gmatu(nmd,kd,6) + dmat(nm,k,5) * gmatu(nm ,ku,1) + dmat(nm,k,6) * gmatu(nmd,ku,3)
             amatp(nm,k, 2) = dmat(nm,k,4) * gmatu(nmd,k ,2) + dmat(nm,k,3) * gmatu(nmd,kd,5) + dmat(nm,k,6) * gmatu(nmd,ku,1)
             amatp(nm,k, 3) = dmat(nm,k,2) * gmatu(nm ,k ,4) + dmat(nm,k,1) * gmatu(nm ,kd,6) + dmat(nm,k,5) * gmatu(nm ,ku,3)
             amatp(nm,k, 4) = dmat(nm,k,2) * gmatu(nm ,k ,1) + dmat(nm,k,4) * gmatu(nmd,k ,3) +  &
                              dmat(nm,k,1) * gmatu(nm ,kd,2) + dmat(nm,k,3) * gmatu(nmd,kd,4) + gmatw(nm,k,1)
             amatp(nm,k, 5) = dmat(nm,k,2) * gmatu(nm ,k ,5) + dmat(nm,k,4) * gmatu(nmd,k ,6) +  &
                              dmat(nm,k,5) * gmatu(nm ,ku,2) + dmat(nm,k,6) * gmatu(nmd,ku,4)
             amatp(nm,k, 8) = dmat(nm,k,2) * gmatu(nm ,k ,3) + dmat(nm,k,1) * gmatu(nm ,kd,4)
             amatp(nm,k,10) = dmat(nm,k,4) * gmatu(nmd,k ,1) + dmat(nm,k,3) * gmatu(nmd,kd,2)
             amatp(nm,k,12) = dmat(nm,k,2) * gmatu(nm ,k ,6) + dmat(nm,k,5) * gmatu(nm ,ku,4)
             amatp(nm,k,14) = dmat(nm,k,4) * gmatu(nmd,k ,5) + dmat(nm,k,6) * gmatu(nmd,ku,2)
             amatp(nm,k,16) = dmat(nm,k,1) * gmatu(nm ,kd,1) + dmat(nm,k,3) * gmatu(nmd,kd,3)
             amatp(nm,k,17) = dmat(nm,k,5) * gmatu(nm ,ku,5) + dmat(nm,k,6) * gmatu(nmd,ku,6)
             amatp(nm,k,18) = dmat(nm,k,1) * gmatu(nm ,kd,3)
             amatp(nm,k,20) = dmat(nm,k,3) * gmatu(nmd,kd,1)
             amatp(nm,k,22) = dmat(nm,k,5) * gmatu(nm ,ku,6)
             amatp(nm,k,24) = dmat(nm,k,6) * gmatu(nmd,ku,5)
             !
             if ( ihydro == 1 ) then
                !
                do j = 3, kmax-1
                   !
                   amatp(nm,k,ishif(j)) = 0.
                   !
                enddo
                !
                amatp(nm,k,1) = amatp(nm,k,1) + gmatw(nm,k,2) + gmatw(nm,k,3) - kwu * gmatw(nm,ku,1)
                !
                do j = 1, kmax-2
                   !
                   amatp(nm,k,ishif(j)) = amatp(nm,k,ishif(j)) + gmatw(nm,k,2*j+2) + gmatw(nm,k,2*j+3) - kwu * gmatw(nm,ku,2*j) - kwu * gmatw(nm,ku,2*j+1)
                   !
                enddo
                !
                amatp(nm,k,ishif(kmax-1)) = amatp(nm,k,ishif(kmax-1)) + gmatw(nm,k,2*kmax) - kwu * gmatw(nm,ku,2*kmax-2)
                !
             else
                !
                amatp(nm,k,1) = amatp(nm,k,1) + gmatw(nm,k,2) - kwu * gmatw(nm,ku,1)
                amatp(nm,k,5) = amatp(nm,k,5)                 - kwu * gmatw(nm,ku,2)
                !
             endif
             !
             rhsp(nm,k) = ( dmat(nm,k,2) * u1(nm,k ) + dmat(nm,k,4) * u1(nmd,k ) +  &
                            dmat(nm,k,1) * u1(nm,kd) + dmat(nm,k,3) * u1(nmd,kd) +  &
                            dmat(nm,k,5) * u1(nm,ku) + dmat(nm,k,6) * u1(nmd,ku) +  &
                            w1(nm,k-1) - kwu*w1(nm,k) ) / (dt*theta)
             !
          enddo
          !
       enddo
       !
       ! reduce the pressure Poisson equation
       !
       do l = 1, qlay
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m)
             !
             amatp(nm,qmax, 1) = amatp(nm,qmax, 1) + real(qlay+2-l)*amatp(nm,kmax-l+1, 4)
             amatp(nm,qmax, 2) = amatp(nm,qmax, 2) + real(qlay+2-l)*amatp(nm,kmax-l+1,10)
             amatp(nm,qmax, 3) = amatp(nm,qmax, 3) + real(qlay+2-l)*amatp(nm,kmax-l+1, 8)
             amatp(nm,qmax, 4) = amatp(nm,qmax, 4) + real(qlay+2-l)*amatp(nm,kmax-l+1,16)
             amatp(nm,qmax, 8) = amatp(nm,qmax, 8) + real(qlay+2-l)*amatp(nm,kmax-l+1,18)
             amatp(nm,qmax,10) = amatp(nm,qmax,10) + real(qlay+2-l)*amatp(nm,kmax-l+1,20)
             rhsp (nm,qmax   ) = rhsp (nm,qmax   ) + real(qlay+2-l)*rhsp (nm,kmax-l+1   )
             !
          enddo
          !
       enddo
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m)
          !
          do k = 1, kmax
             !
             if ( .not. amatp(nm,k,1) /= 0. ) then
                amatp(nm,k,:) =  0.
                amatp(nm,k,1) = -1.
                rhsp (nm,k  ) =  0.
             endif
             !
          enddo
          !
       enddo
       !
       ! solve the Poisson equation
       !
       if ( qmax == 1 ) then
          !
          call tridiag ( amatp(1,1,2), amatp(1,1,1), amatp(1,1,3), rhsp(1,1), dq(1,1), kgrpnt )
          if (STPNOW()) return
          !
       else
          !
          if ( lprecon ) then
             !
!TIMG             call SWTSTA(296)
             if ( icond == 1 ) then
                !
                ! compute incomplete LU factorization restricted to diagonal used as split preconditioner
                !
                call iluds ( amatp(1:mcgrd,1:qmax,1:nconct) )
                !
             else if ( icond == 2 ) then
                !
                ! compute incomplete LU factorization restricted to diagonal used as right preconditioner
                !
                call iludr ( amatp(1:mcgrd,1:qmax,1:nconct) )
                !
             else if ( icond == 3 ) then
                !
                ! compute classical incomplete LU factorization
                !
                call ilu ( amatp(1:mcgrd,1:qmax,1:nconct) )
                !
             endif
!TIMG             call SWTSTO(296)
             !
          endif
          !
!TIMG          call SWTSTA(297)
          call bicgstab ( amatp(1:mcgrd,1:qmax,1:nconct), rhsp(1:mcgrd,1:qmax), dq(1:mcgrd,1:qmax) )
!TIMG          call SWTSTO(297)
          if (STPNOW()) return
          !
       endif
       !
       ! exchange pressure corrections with neighbouring subdomains
       !
       call SWEXCHG ( dq(1:mcgrd,1:qmax), kgrpnt, 1, qmax )
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
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm,k) = u1(nm,k) - dt*theta*( gmatu(nm,k,1)*dq(nm ,kd) + gmatu(nm,k,2)*dq(nm,k ) + gmatu(nm,k,3)*dq(nmu,kd) + &
                                                 gmatu(nm,k,4)*dq(nmu,k ) + gmatu(nm,k,5)*dq(nm,ku) + gmatu(nm,k,6)*dq(nmu,ku) )
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
          do m = mfu, ml
             !
             nm = kgrpnt(m)
             !
             if ( wets(nm) == 1 ) then
                !
                if ( ihydro == 1 ) then
                   !
                   do j = 0, kmax-1
                      !
                      kd = max(k+j  ,1   )
                      kd = min(kd   ,kmax)
                      ku = min(k+j+1,kmax)
                      !
                      w1(nm,k) = w1(nm,k) - dt*theta*( gmatw(nm,k+1,2*j+1)*dq(nm,kd) + gmatw(nm,k+1,2*j+2)*dq(nm,ku) )
                      !
                   enddo
                   !
                else
                   !
                   w1(nm,k) = w1(nm,k) - dt*theta*( gmatw(nm,k+1,1)*dq(nm,kd) + gmatw(nm,k+1,2)*dq(nm,k+1) )
                   !
                endif
                !
             endif
             !
          enddo
          !
       enddo
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          if ( wets(nm) == 1 ) then
             !
             w1(nm,kmax) = 0.25 * rdx * ( 3.*(u1(nm,kmax)+u1(nmd,kmax)) - u1(nm,kmax-1) - u1(nmd,kmax-1) ) * ( zku(nm,kmax) - zku(nmd,kmax) )
             !
          endif
          !
       enddo
       !
    endif
    !
    ! impose Neumann condition for w-velocity at boundaries, if appropriate
    !
    if ( horwinc ) then
       !
       if ( ihydro == 1 .or. ihydro == 2 ) then
          !
          if ( LMXF ) w1(nmf ,:) = w1(nmfu,:)
          if ( LMXL ) w1(nmlu,:) = w1(nml ,:)
          !
       endif
       !
    endif
    !
    ! store layer interfaces in u-point
    !
    zkuo = zku
    !
    ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    ! solve global continuity equation by means of a second order predictor-corrector scheme
    ! (jhk = 1: RK3, jhk = 2: Hancock, jhk = 3: first order Euler)
    !
    do j = jhk, 3
       !
       if ( j <= 2 ) then
          !
          ! predictor phases
          !
          fh = 0.5
          !
       else
          !
          ! corrector phase
          !
          fh = 1.
          !
       endif
       !
       ! compute the mass flux per layer
       !
       do k = 1, kmax
          !
          do m = mf, ml
             !
             nm = kgrpnt(m)
             !
             qx(nm,k) = hku(nm,k)*u1(nm,k)
             !
          enddo
          !
       enddo
       !
       ! compute the depth-integrated mass flux
       !
       work(:,1) = 0.
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             do k = 1, kmax
                !
                work(nm,1) = work(nm,1) + qx(nm,k)
                !
             enddo
             !
          endif
          !
       enddo
       !
       ! compute the depth-averaged u-velocity
       !
       if ( j == jhk ) then
          !
          udep = 0.
          !
          do m = mf, ml
             !
             nm = kgrpnt(m)
             !
             if ( wetu(nm) == 1 ) then
                !
                udep(nm) = work(nm,1) / hu(nm)
                !
             endif
             !
          enddo
          !
       endif
       !
       ! compute the water level
       !
       if ( .not.linswe ) then
          !
          do m = mfu, mend
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             s1(nm) = s0(nm) - fh * dt * rdx * ( work(nm,1) - work(nmd,1) )
             !
          enddo
          !
       else
          !
          do m = mfu, mend
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             s1(nm) = s0(nm) - fh * dt * rdx * ( dpu(nm)*udep(nm) - dpu(nmd)*udep(nmd) )
             !
          enddo
          !
       endif
       !
       ! exchange water levels with neighbouring subdomains
       !
       call SWEXCHG ( s1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! copy to virtual cells at boundaries except for both water level and Riemann invariant openings
       !
       if ( ibl /= 2 .and. ibl /= 6 .and. LMXF ) s1(nmf ) = s1(nmfu)
       if ( ibr /= 2 .and. ibr /= 6 .and. LMXL ) s1(nmlu) = s1(nml )
       !
       ! in case of outflow condition, no change in water depth
       !
       if ( ibl == 10 ) s1(nmf ) = s1(nmfu) - dps(nmfu) + dps(nmfuu)
       if ( ibr == 10 ) s1(nmlu) = s1(nml ) - dps(nml ) + dps(nmld )
       !
       if ( j <= 2 ) then
          !
          ! correct the water depth for mass flux
          !
          call SwashUpdDepu ( udep )
          !
          ! recompute layer interfaces and layer thicknesses at u-point
          !
          zku(:,0) = zku(:,kmax) + hu(:)
          call sigmacoor( zku, hu, mcgrd )
          !
          do k = 1, kmax
             hku(:,k) = zku(:,k-1) - zku(:,k)
          enddo
          !
       endif
       !
    enddo
    !
    ! add mass source due to internal wave generation
    !
    if ( iwvgen /= 0 ) then
       !
       do m = mfu, mend
          !
          nm = kgrpnt(m)
          !
          s1(nm) = s1(nm) + dt * srcm(nm)
          !
       enddo
       !
       ! exchange water levels with neighbouring subdomains
       !
       call SWEXCHG ( s1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    endif
    !
    ! include effect of ambient current
    !
    if ( iamb == 1 ) then
       !
       so = s1
       !
       ! solve global continuity equation by means of the Heun method
       !
       do j = 1, 2
          !
          work(:,2) = 0.
          !
          if ( j == 1 ) then
             !
             do m = mfu, mend
                !
                md = m - 1
                mu = m + 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                nmu = kgrpnt(mu)
                !
                work(nm,2) = 0.25 * rdx * ( ( Ubck(nmu,1) + Ubck(nm,1) ) * ( so(nmu) + so(nm) ) - ( Ubck(nm,1) + Ubck(nmd,1) ) * ( so(nm) + so(nmd) ) )
                !
             enddo
             !
          else
             !
             do m = mfu, mend
                !
                md = m - 1
                mu = m + 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                nmu = kgrpnt(mu)
                !
                fac1 = thetamb * s1(nmd) + (1.-thetamb) * so(nmd)
                fac  = thetamb * s1(nm ) + (1.-thetamb) * so(nm )
                fac2 = thetamb * s1(nmu) + (1.-thetamb) * so(nmu)
                !
                work(nm,2) = 0.25 * rdx * ( ( Ubck(nmu,1) + Ubck(nm,1) ) * ( fac2 + fac ) - ( Ubck(nm,1) + Ubck(nmd,1) ) * ( fac + fac1 ) )
                !
             enddo
             !
          endif
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m)
             !
             s1(nm) = so(nm) - dt * work(nm,2)
             !
          enddo
          !
          ! exchange water levels with neighbouring subdomains
          !
          call SWEXCHG ( s1, kgrpnt, 1, 1 )
          if (STPNOW()) return
          !
          ! copy to virtual cells at boundaries except for both water level and Riemann invariant openings
          !
          if ( ibl /= 2 .and. ibl /= 6 .and. LMXF ) s1(nmf ) = s1(nmfu)
          if ( ibr /= 2 .and. ibr /= 6 .and. LMXL ) s1(nmlu) = s1(nml )
          !
          ! in case of outflow condition, no change in water depth
          !
          if ( ibl == 10 ) s1(nmf ) = s1(nmfu) - dps(nmfu) + dps(nmfuu)
          if ( ibr == 10 ) s1(nmlu) = s1(nml ) - dps(nml ) + dps(nmld )
          !
       enddo
       !
    else if ( iamb == 2 ) then
       !
       so = s1
       !
       ! solve global continuity equation by means of the Heun method
       !
       do j = 1, 2
          !
          work(:,2) = 0.
          !
          if ( j == 1 ) then
             !
             do m = mfu, mend
                !
                md = m - 1
                mu = m + 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                nmu = kgrpnt(mu)
                !
                work(nm,2) = 0.5 * rdx * ( Ubck(nm,1) * ( so(nmu) + so(nm) ) - Ubck(nmd,1) * ( so(nm) + so(nmd) ) )
                !
             enddo
             !
          else
             !
             do m = mfu, mend
                !
                md = m - 1
                mu = m + 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                nmu = kgrpnt(mu)
                !
                fac1 = thetamb * s1(nmd) + (1.-thetamb) * so(nmd)
                fac  = thetamb * s1(nm ) + (1.-thetamb) * so(nm )
                fac2 = thetamb * s1(nmu) + (1.-thetamb) * so(nmu)
                !
                work(nm,2) = 0.5 * rdx * ( Ubck(nm,1) * ( fac2 + fac ) - Ubck(nmd,1) * ( fac + fac1 ) )
                !
             enddo
             !
          endif
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m)
             !
             s1(nm) = so(nm) - dt * work(nm,2)
             !
          enddo
          !
          ! exchange water levels with neighbouring subdomains
          !
          call SWEXCHG ( s1, kgrpnt, 1, 1 )
          if (STPNOW()) return
          !
          ! copy to virtual cells at boundaries except for both water level and Riemann invariant openings
          !
          if ( ibl /= 2 .and. ibl /= 6 .and. LMXF ) s1(nmf ) = s1(nmfu)
          if ( ibr /= 2 .and. ibr /= 6 .and. LMXL ) s1(nmlu) = s1(nml )
          !
          ! in case of outflow condition, no change in water depth
          !
          if ( ibl == 10 ) s1(nmf ) = s1(nmfu) - dps(nmfu) + dps(nmfuu)
          if ( ibr == 10 ) s1(nmlu) = s1(nml ) - dps(nml ) + dps(nmld )
          !
       enddo
       !
    endif
    !
    ! calculate net mass outflow based on local continuity equation
    !
    if ( ITEST >= 30 ) then
       !
       moutf = 0.
       !
       if ( ihydro == 1 .or. ihydro == 2 ) then
          !
          do m = mfu, mend
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             do k = 1, kmax
                !
                kd = max(k-1,1   )
                ku = min(k+1,kmax)
                !
                kwu = 1.
                if ( k == kmax ) kwu = 0.
                !
                moutf = moutf + dmat(nm,k,2) * u1(nm,k ) + dmat(nm,k,4) * u1(nmd,k ) +  &
                                dmat(nm,k,1) * u1(nm,kd) + dmat(nm,k,3) * u1(nmd,kd) +  &
                                dmat(nm,k,5) * u1(nm,ku) + dmat(nm,k,6) * u1(nmd,ku) +  &
                                w1(nm,k-1) - kwu*w1(nm,k)
                !
             enddo
             !
          enddo
          !
       endif
       !
       ! accumulate net mass outflow over all subdomains
       !
       call SWREDUCE ( moutf, 1, SWREAL, SWSUM )
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
       do m = mfu, mend
          !
          md = m - 1
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          nmu = kgrpnt(mu)
          !
          ! compute displaced volume in the cell center
          vol = vol + s1(nm)
          !
          ! compute potential energy in the cell center
          ener = ener + 0.5 * grav * dx * s1(nm) * s1(nm)
          !
          ! compute kinetic energy in the cell center
          fac1 = ( s1(nm ) + dps(nm ) ) * dx
          fac2 = ( s1(nmd) + dps(nmd) ) * dx
          ener = ener + 0.125 * ( fac1 + fac2 ) * udep(nmd) * udep(nmd)
          fac2 = ( s1(nmu) + dps(nmu) ) * dx
          ener = ener + 0.125 * ( fac1 + fac2 ) * udep(nm ) * udep(nm )
          !
       enddo
       !
       ! accumulate displaced volume and energy over all subdomains
       !
       call SWREDUCE ( vol, 1, SWREAL, SWSUM )
       call SWREDUCE (ener, 1, SWREAL, SWSUM )
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
    call sigmacoor ( zksnew, work, mcgrd )
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
          do m = mfu, mend
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             if ( wets(nm) == 1 ) then
                !
                wom(nm,k) = wom(nm,k+1) - rdx * ( hku(nm,k+1)*u1(nm,k+1) - hku(nmd,k+1)*u1(nmd,k+1) ) - ( hksnew(nm,k+1) - hks(nm,k+1) ) / dt
                !
             else
                !
                wom(nm,k) = 0.
                !
             endif
             !
          enddo
          !
          if ( LMXF ) wom(nmf ,k) = wom(nmfu,k)
          if ( LMXL ) wom(nmlu,k) = wom(nml ,k)
          !
       enddo
       !
       ! check if relative vertical velocity at surface is zero
       !
       do m = mfu, mend
          !
          nm = kgrpnt(m)
          !
          if ( abs(wom(nm,0)) > epswom ) then
             !
             write (msgstr,'(a,i5,a,e9.3,a)') 'nonzero relative vertical velocity at surface in m=',m+MXF-2,'; omega = ',wom(nm,0),' m/s'
             call msgerr (2, trim(msgstr) )
             !
             wom(nm,0) = 0.
             !
          endif
          !
       enddo
       !
    else
       !
       do m = mfu, mend
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          if ( wets(nm) == 1 ) then
             !
             do k = 1, kmax-1
                !
                fac1 = 0.5 * ( u1(nmd,k  ) + u1(nm,k  ) )
                fac2 = 0.5 * ( u1(nmd,k+1) + u1(nm,k+1) )
                !
                fac  = ( fac2 * hks(nm,k) + fac1 * hks(nm,k+1) ) / ( hks(nm,k) + hks(nm,k+1) )
                !
                wom(nm,k) = w1(nm,k) - ( zksnew(nm,k) - zks(nm,k) ) / dt - fac * rdx *( zkuo(nm,k) - zkuo(nmd,k) )
                !
             enddo
             !
             wom(nm,0   ) = 0.
             wom(nm,kmax) = 0.
             !
          else if ( brks(nm) == 1 ) then
             !
             ! hydrostatic pressure is assumed at steep front of breaking wave, so relative vertical velocity is derived from local continuity equation
             !
             wom(nm,kmax) = 0.
             !
             do k = kmax-1, 0, -1
                !
                wom(nm,k) = wom(nm,k+1) - rdx * ( hku(nm,k+1)*u1(nm,k+1) - hku(nmd,k+1)*u1(nmd,k+1) ) - ( hksnew(nm,k+1) - hks(nm,k+1) ) / dt
                !
             enddo
             !
             if ( abs(wom(nm,0)) > epswom ) then
                !
                write (msgstr,'(a,i5,a,e9.3,a)') 'nonzero relative vertical velocity at surface in m=',m+MXF-2,'; omega = ',wom(nm,0),' m/s'
                if ( ITEST >= 50 ) call msgerr (1, trim(msgstr) )
                !
                wom(nm,0) = 0.
                !
             endif
             !
          else
             !
             wom(nm,:) = 0.
             !
          endif
          !
       enddo
       !
       if ( LMXF ) wom(nmf ,:) = wom(nmfu,:)
       if ( LMXL ) wom(nmlu,:) = wom(nml ,:)
       !
    endif
    !
    ! compute the maximum CFL number
    !
    cflmax = -999.
    !
    do m = mf, ml
       !
       mu = m + 1
       !
       nm  = kgrpnt(m )
       nmu = kgrpnt(mu)
       !
       dep = 0.5 * ( s1(nm) + dps(nm) + s1(nmu) + dps(nmu) )
       !
       cfl = rdx * dt * ( sqrt(grav*dep) + maxval(abs(u1(nm,:))) )
       if ( cfl > cflmax ) cflmax = cfl
       !
    enddo
    !
    ! find maximum of CFL number over all subdomains
    !
    call SWREDUCE ( cflmax, 1, SWREAL, SWMAX )
    !
    ! exchange velocities with neighbouring subdomains (if appropriate)
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       call SWEXCHG ( u1, kgrpnt, 1, kmax )
       if ( horwinc ) call SWEXCHG ( w1, kgrpnt, 0, kmax )
       !
    endif
    call SWEXCHG ( wom, kgrpnt, 0, kmax )
    if (STPNOW()) return
    !
    ! apply wave absorption by means of sponge layer, if appropriate
    !
    if ( spwidl > 0. ) then
       !
       tbndx = 0.
       if ( LMXF ) then
          tbndx(1,1:kmax) = u1(nmf ,1:kmax)
          tbndx(1,kmax+1) = s1(nmfu)
       endif
       call SWREDUCE ( tbndx, kmax+1, SWREAL, SWSUM )
       if (STPNOW()) return
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             u1(nm,1:kmax) = (1.-sponxl(nm))*u1(nm,1:kmax) + sponxl(nm)*tbndx(1,1:kmax)
             !
          endif
          !
       enddo
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          if ( wets(nm) == 1 ) then
             !
             s1(nm) = (1.-sponxl(nmd))*s1(nm) + sponxl(nmd)*tbndx(1,kmax+1)
             !
          endif
          !
       enddo
       !
    endif
    !
    if ( spwidr > 0. ) then
       !
       tbndx = 0.
       if ( LMXL ) then
          tbndx(1,1:kmax) = u1(nml,1:kmax)
          tbndx(1,kmax+1) = s1(nml)
       endif
       call SWREDUCE ( tbndx, kmax+1, SWREAL, SWSUM )
       if (STPNOW()) return
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             u1(nm,1:kmax) = (1.-sponxr(nm))*u1(nm,1:kmax) + sponxr(nm)*tbndx(1,1:kmax)
             !
          endif
          !
       enddo
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m)
          !
          if ( wets(nm) == 1 ) then
             !
             s1(nm) = (1.-sponxr(nm))*s1(nm) + sponxr(nm)*tbndx(1,kmax+1)
             !
          endif
          !
       enddo
       !
    endif
    !
 101 format (2x,'the net mass outflow is ',e14.8e2)
 102 format (2x,'the total displaced volume of water is ',e14.8e2)
 103 format (2x,'the total energy is ',e14.8e2)
    !
end subroutine SwashExpLay1DHflow
