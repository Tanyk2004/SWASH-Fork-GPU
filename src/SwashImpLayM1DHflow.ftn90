subroutine SwashImpLayM1DHflow ( ibl, ibr, kgrpnt )
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
!    1.00,    March 2015: New subroutine
!
!   Purpose
!
!   Performs the time integration for the non-hydrostatic, layer-averaged 1D shallow water equations
!
!   Note: this subroutine solves the layered-averaged 1D shallow water equations, which are strictly mass and
!         momentum conservative at the discrete level, and strictly energy conservative in discrete space,
!         but nearly energy conservative for discrete time. The discrete advective operator is skew-symmetric
!         (both time and space), and the negative transpose of the discrete surface gradient operator is equal
!         to the divergence operator. However, due to the time splitting between the surface gradient and other
!         terms of the momentum equation, the energy is not strictly conservative in time. Nevertheless, the
!         current approach has good conservation properties and is stable, even for high waves.
!
!   Method
!
!   The time integration with respect to the continuity equation and the water level gradient of the
!   u-momentum equation is based on a theta-scheme. Only a value of theta = 0.5 is taken.
!
!   The time integration with respect to the horizontal advective and viscosity terms is based on
!   the Crank-Nicolson scheme, while that for the bottom friction is based on Euler implicit and for
!   the non-hydrostatic pressure gradient a semi-implicit approach is employed (theta-scheme).
!
!   Both vertical advective and viscosity terms are treated using the Crank-Nicolson scheme as well.
!   This results in a tri-diagonal system.
!
!   The space discretization of the horizontal advective and viscosity terms is strictly momentum
!   conservative and are approximated with central differences, so that the advection term is
!   skew-symmetric and the viscosity term is symmetric.
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
!   on central differences in a finite volume fashion.
!
!   The non-hydrostatic pressure is obtained by means of the second order accurate pressure correction technique.
!
!   Note: in order to preserve skew-symmetry of the advective operator, only velocity is prescribed at the boundary.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
    use m_genarr, only: work, wrk
    use m_parall
    use SwashSolvedata, only: icond, iamout
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
    integer                         :: jj       ! iteration counter
    integer                         :: k        ! loop counter over vertical layers
    integer                         :: kd       ! index of layer k-1
    integer                         :: kinf     ! k-index of layer interface with largest error in solution
    integer                         :: ku       ! index of layer k+1
    integer                         :: l        ! loop counter
    integer                         :: m        ! loop counter over horizontal grid points
    integer                         :: maxit    ! maximum number of iterations
    integer                         :: md       ! index of point m-1
    integer                         :: mdd      ! index of point m-2
    integer                         :: mend     ! end index of loop over u-points
    integer                         :: minf     ! m-index of point with largest error in solution
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
    real                            :: ener     ! total energy of closed system
    real                            :: epslin   ! required accuracy in the linear solver
    real                            :: fac      ! a factor
    real                            :: fac1     ! another factor
    real                            :: fac2     ! some other factor
    real                            :: fac3     ! auxiliary factor
    real                            :: fac4     ! auxiliary factor
    real                            :: kwd      ! =1. if layer k-1 exists otherwise 0.
    real                            :: kwu      ! =1. if layer k+1 exists otherwise 0.
    real                            :: moutf    ! net mass outflow
    real                            :: qf       ! updated mass flux at present face
    real                            :: rdx      ! reciprocal of mesh size
    real                            :: reps     ! accuracy of the final approximation
    real                            :: res      ! residual
    real                            :: resm     ! maximum error
    real                            :: rhou     ! density of water in velocity point
    real                            :: rval     ! auxiliary real
    real                            :: s0mx     ! maximum value of water level
    real                            :: sumqf    ! sum of outgoing mass fluxes over present cell
    real                            :: theta3   ! implicitness factor for non-hydrostatic pressure gradient
    real                            :: thetau   ! implicitness factor for vertical terms in u-momentum equation
    real                            :: thetaw   ! implicitness factor for vertical terms in w-momentum equation
    real                            :: u        ! u-velocity at point different from its point of definition
    real                            :: u0mx     ! maximum value of horizontal velocity
    real                            :: utot     ! velocity magnitude
    real                            :: vol      ! total displaced volume of water
    real                            :: w        ! w-velocity at point different from its point of definition
    !
    logical                         :: EQREAL   ! compares two reals
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
    if (ltrace) call strace (ient,'SwashImpLayM1DHflow')
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
    teta   = 0.5
    teta2  = 0.5
    theta3 = pnums(5)
    thetau = 0.5
    thetaw = 0.5
    !
    ! adapt theta values underneath the floating object, if appropriate
    !
    if ( ifloat /= 0 ) then
       !
       do m = mf, ml
          nm = kgrpnt(m)
          if ( presu(nm) == 1 ) teta(nm) = pship(2)
       enddo
       !
       call SWEXCHG ( teta, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       do m = mfu, mend
          nm = kgrpnt(m)
          if ( presp(nm) == 1 ) teta2(nm) = pship(2)
       enddo
       !
       call SWEXCHG ( teta2, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    endif
    !
    ! build the u-momentum equation
    !
    msta = mf + 1
    mend = ml - 1
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
          bdx  (nm,:  ) = 0.
          bux  (nm,:  ) = 0.
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
    ! compute horizontal advection term (momentum conservative) at internal u-point (implicit)
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
          if ( hkumn(nm,k) > 0. ) then
             !
             fac1 = 0.5 * rdx * qm(nmu,k) / hkumn(nm,k)
             fac2 = 0.5 * rdx * qm(nm ,k) / hkumn(nm,k)
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
             bdx(nm,k) = -fac2
             bux(nm,k) =  fac1
             !
             amatu(nm,k,1) = amatu(nm,k,1) - fac1 + fac2
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
             if ( wetu(nm) * wetu(nmu) == 1 .and. (m /= mf .or. .not.LMXF) .and. (m /= ml .or. .not.LMXL) ) then
                !
                fac = 0.5 * rdx / hkum(nm,k)
                !
                fac1 = zks(nmu,k-1) - zks(nm,k-1)
                fac2 = zks(nmu,k  ) - zks(nm,k  )
                !
                if ( presu(nm) == 0 ) then
                   !
                   ! free surface flow
                   !
                   gmatu(nm,k,1) = (-hks(nm ,k) - fac1) * fac
                   gmatu(nm,k,2) = (-hks(nm ,k) + fac2) * fac
                   gmatu(nm,k,3) = ( hks(nmu,k) - fac1) * fac
                   gmatu(nm,k,4) = ( hks(nmu,k) + fac2) * fac
                   gmatu(nm,k,5) = 0.
                   gmatu(nm,k,6) = 0.
                   !
                else
                   !
                   ! pressurized flow
                   !
                   if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
                      !
                      ! changeover from free surface to pressurized
                      !
                      fac3 = 1. / ( hks(nmu,kd) + hks(nmu,k ) )
                      fac4 = 1. / ( hks(nmu,k ) + hks(nmu,ku) )
                      !
                      gmatu(nm,k,1) = (-hks(nm,k) - fac1) * fac
                      gmatu(nm,k,2) = (-hks(nm,k) + fac2) * fac
                      gmatu(nm,k,3) = -fac1 * hks(nmu,k)*fac3 * fac
                      gmatu(nm,k,4) = (2.*hks(nmu,k) - fac1 * hks(nmu,kd)*fac3 + fac2 * hks(nmu,ku)*fac4 ) * fac
                      gmatu(nm,k,5) = 0.
                      gmatu(nm,k,6) = fac2 * hks(nmu,k)*fac4 * fac
                      !
                   else if ( presp(nm) == 1 .and. presp(nmu) == 0 ) then
                      !
                      ! changeover from pressurized to free surface
                      !
                      fac3 = 1. / ( hks(nm,kd) + hks(nm,k ) )
                      fac4 = 1. / ( hks(nm,k ) + hks(nm,ku) )
                      !
                      gmatu(nm,k,1) = -fac1 * hks(nm,k)*fac3 * fac
                      gmatu(nm,k,2) = (-2.*hks(nm,k) - fac1 * hks(nm,kd)*fac3 + fac2 * hks(nm,ku)*fac4 ) * fac
                      gmatu(nm,k,3) = ( hks(nmu,k) - fac1) * fac
                      gmatu(nm,k,4) = ( hks(nmu,k) + fac2) * fac
                      gmatu(nm,k,5) = fac2 * hks(nm,k)*fac4 * fac
                      gmatu(nm,k,6) = 0.
                      !
                   else
                      !
                      fac1 = 0.5 * rdx * fac1 / ( hkum(nm,kd) + hkum(nm,k ) )
                      fac2 = 0.5 * rdx * fac2 / ( hkum(nm,k ) + hkum(nm,ku) )
                      !
                      gmatu(nm,k,1) = -fac1
                      gmatu(nm,k,2) = ( -rdx*hks(nm ,k) - hkum(nm,kd)*fac1 + hkum(nm,ku)*fac2 ) / hkum(nm,k)
                      gmatu(nm,k,3) = -fac1
                      gmatu(nm,k,4) = (  rdx*hks(nmu,k) - hkum(nm,kd)*fac1 + hkum(nm,ku)*fac2 ) / hkum(nm,k)
                      gmatu(nm,k,5) =  fac2
                      gmatu(nm,k,6) =  fac2
                      !
                   endif
                   !
                endif
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
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmu = kgrpnt(mu)
          !
          if ( presu(nm) == 1 ) then
             !
             ! pressurized flow
             !
             if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
                !
                ! changeover from free surface to pressurized
                !
                gmatu(nm,kmax,4) = gmatu(nm,kmax,4) + 2.*gmatu(nm,kmax,6)
                gmatu(nm,kmax,3) = gmatu(nm,kmax,3) -    gmatu(nm,kmax,6)
                gmatu(nm,kmax,6) = 0.
                !
                gmatu(nm,1,4) = gmatu(nm,1,4) + 2.*gmatu(nm,1,3)
                gmatu(nm,1,6) = gmatu(nm,1,6) -    gmatu(nm,1,3)
                !
             else if ( presp(nm) == 1 .and. presp(nmu) == 0 ) then
                !
                ! changeover from pressurized to free surface
                !
                gmatu(nm,kmax,2) = gmatu(nm,kmax,2) + 2.*gmatu(nm,kmax,5)
                gmatu(nm,kmax,1) = gmatu(nm,kmax,1) -    gmatu(nm,kmax,5)
                gmatu(nm,kmax,5) = 0.
                !
                gmatu(nm,1,2) = gmatu(nm,1,2) + 2.*gmatu(nm,1,1)
                gmatu(nm,1,5) = gmatu(nm,1,5) -    gmatu(nm,1,1)
                !
             else
                !
                gmatu(nm,kmax,2) = gmatu(nm,kmax,2) + 2.*gmatu(nm,kmax,5)
                gmatu(nm,kmax,1) = gmatu(nm,kmax,1) -    gmatu(nm,kmax,5)
                gmatu(nm,kmax,4) = gmatu(nm,kmax,4) + 2.*gmatu(nm,kmax,6)
                gmatu(nm,kmax,3) = gmatu(nm,kmax,3) -    gmatu(nm,kmax,6)
                gmatu(nm,kmax,5) = 0.
                gmatu(nm,kmax,6) = 0.
                !
                gmatu(nm,1,2) = gmatu(nm,1,2) + 2.*gmatu(nm,1,1)
                gmatu(nm,1,5) = gmatu(nm,1,5) -    gmatu(nm,1,1)
                gmatu(nm,1,4) = gmatu(nm,1,4) + 2.*gmatu(nm,1,3)
                gmatu(nm,1,6) = gmatu(nm,1,6) -    gmatu(nm,1,3)
                !
             endif
             !
          endif
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
             if ( wetu(nm) * wetu(nmu) == 1 .and. (m /= mf .or. .not.LMXF) .and. (m /= ml .or. .not.LMXL) ) then
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
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmu = kgrpnt(mu)
          !
          gmatu(nm,kmax,2) = gmatu(nm,kmax,2) + 2.*gmatu(nm,kmax,5)
          gmatu(nm,kmax,1) = gmatu(nm,kmax,1) -    gmatu(nm,kmax,5)
          gmatu(nm,kmax,4) = gmatu(nm,kmax,4) + 2.*gmatu(nm,kmax,6)
          gmatu(nm,kmax,3) = gmatu(nm,kmax,3) -    gmatu(nm,kmax,6)
          gmatu(nm,kmax,5) = 0.
          gmatu(nm,kmax,6) = 0.
          !
          if ( presu(nm) == 0 ) then
             !
             ! free surface flow
             !
             gmatu(nm,1,2) = gmatu(nm,1,2) - gmatu(nm,1,1)
             gmatu(nm,1,4) = gmatu(nm,1,4) - gmatu(nm,1,3)
             !
          else
             !
             ! pressurized flow
             !
             if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
                !
                ! changeover from free surface to pressurized
                !
                gmatu(nm,1,2) = gmatu(nm,1,2) -    gmatu(nm,1,1)
                gmatu(nm,1,4) = gmatu(nm,1,4) + 2.*gmatu(nm,1,3)
                gmatu(nm,1,6) = gmatu(nm,1,6) -    gmatu(nm,1,3)
                !
             else if ( presp(nm) == 1 .and. presp(nmu) == 0 ) then
                !
                ! changeover from pressurized to free surface
                !
                gmatu(nm,1,2) = gmatu(nm,1,2) + 2.*gmatu(nm,1,1)
                gmatu(nm,1,5) = gmatu(nm,1,5) -    gmatu(nm,1,1)
                gmatu(nm,1,4) = gmatu(nm,1,4) -    gmatu(nm,1,3)
                !
             else
                !
                gmatu(nm,1,2) = gmatu(nm,1,2) + 2.*gmatu(nm,1,1)
                gmatu(nm,1,5) = gmatu(nm,1,5) -    gmatu(nm,1,1)
                gmatu(nm,1,4) = gmatu(nm,1,4) + 2.*gmatu(nm,1,3)
                gmatu(nm,1,6) = gmatu(nm,1,6) -    gmatu(nm,1,3)
                !
             endif
             !
          endif
          !
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
       else if ( iproj == 2 .and. theta3 /= 1. ) then
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
                   rhsu(nm,k) = rhsu(nm,k) - (1.-theta3) * ( gmatu(nm,k,1)*q(nm ,kd) + gmatu(nm,k,2)*q(nm,k ) + gmatu(nm,k,3)*q(nmu,kd)  &
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
    ! compute horizontal viscosity term at internal u-point (implicit)
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
             if ( hkum(nm,k) > 0. ) then
                !
                fac1 = rdx * rdx * hvisc * hks(nmu,k) / hkum(nm,k)
                fac2 = rdx * rdx * hvisc * hks(nm ,k) / hkum(nm,k)
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
                bdx(nm,k) = bdx(nm,k) - fac2
                bux(nm,k) = bux(nm,k) - fac1
                !
                amatu(nm,k,1) = amatu(nm,k,1) + fac1 + fac2
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
             if ( hkum(nm,k) > 0. ) then
                !
                fac1 = rdx * rdx * vnu2d(nmu) * hks(nmu,k) / hkum(nm,k)
                fac2 = rdx * rdx * vnu2d(nm ) * hks(nm ,k) / hkum(nm,k)
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
                bdx(nm,k) = bdx(nm,k) - fac2
                bux(nm,k) = bux(nm,k) - fac1
                !
                amatu(nm,k,1) = amatu(nm,k,1) + fac1 + fac2
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
    ! compute vertical terms (implicit)
    !
    do k = 1, kmax
       !
       kd = max(k-1,1   )
       ku = min(k+1,kmax)
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
             ctrkt = 0.5 * kwd * ( wom(nm,k-1) + wom(nmu,k-1) ) / ( hkumn(nm,kd) + hkumn(nm,k ) )
             ctrkb = 0.5 * kwu * ( wom(nm,k  ) + wom(nmu,k  ) ) / ( hkumn(nm,k ) + hkumn(nm,ku) )
             !
             amatu(nm,k,1) = amatu(nm,k,1) - ctrkt + ctrkb
             amatu(nm,k,2) =  ctrkt
             amatu(nm,k,3) = -ctrkb
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
    ! incorporate boundary conditions
    !
    if ( LMXF ) then
       !
       rhsu(nmfu,:) = rhsu(nmfu,:) - bdx(nmfu,:) * u1(nmf,:)
       bdx (nmfu,:) = 0.
       !
    endif
    !
    if ( LMXL ) then
       !
       rhsu(nmld,:) = rhsu(nmld,:) - bux(nmld,:) * u1(nml,:)
       bux (nmld,:) = 0.
       !
    endif
    !
    ! first, a LU decomposition is carried out in vertical direction
    !
    do k = 2, kmax
       do m = msta, mend
          !
          nm = kgrpnt(m)
          !
          amatu(nm,k,2) = amatu(nm,k,2) / amatu(nm,k-1,1)
          amatu(nm,k,1) = amatu(nm,k,1) - amatu(nm,k-1,3) * amatu(nm,k,2)
          !
       enddo
    enddo
    !
    ! next, solve the u-momentum equation using an iterative red-black Jacobi method
    ! in horizontal direction and Gaussian elimination in the vertical direction
    !
    ! determine maximum of velocity and required accuracy
    !
    u0mx   = maxval(abs(u0))
    reps   = 1.e-4
    epslin = max(reps,reps*u0mx)
    call SWREDUCE ( epslin, 1, SWREAL, SWMAX )
    !
    ! initialize
    !
    maxit = 100
    jj    = 0
    resm  = epslin + 1.
    !
    do m = msta, mend
       !
       nm = kgrpnt(m)
       !
       u1(nm,:) = u0(nm,:)
       !
    enddo
    !
    ! start iteration process
    !
 10 if ( resm > epslin .and. jj < maxit ) then
       !
       jj = jj + 1
       !
       resm = 0.
       minf = 0
       kinf = 0
       !
       ! for all odd points:
       !
       ! compute right hand side
       !
       do k = 1, kmax
          !
          do m = msta+1, mend, 2
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             nmu = kgrpnt(mu)
             !
             wrk(nm,k) = rhsu(nm,k) - bux(nm,k)*u1(nmu,k) - bdx(nm,k)*u1(nmd,k)
             !
          enddo
          !
       enddo
       !
       ! solve the equation in vertical direction using Gaussian elimination
       !
       do m = msta+1, mend, 2
          !
          nm = kgrpnt(m)
          !
          ui(nm,1) = wrk(nm,1)
          !
       enddo
       !
       do k = 2, kmax
          !
          do m = msta+1, mend, 2
             !
             nm = kgrpnt(m)
             !
             ui(nm,k) = wrk(nm,k) - amatu(nm,k,2)*ui(nm,k-1)
             !
          enddo
          !
       enddo
       !
       do m = msta+1, mend, 2
          !
          nm = kgrpnt(m)
          !
          ui(nm,kmax) = ui(nm,kmax) / amatu(nm,kmax,1)
          !
       enddo
       !
       do k = kmax-1, 1, -1
          !
          do m = msta+1, mend, 2
             !
             nm = kgrpnt(m)
             !
             ui(nm,k) = ( ui(nm,k) - amatu(nm,k,3)*ui(nm,k+1) ) / amatu(nm,k,1)
             !
          enddo
          !
       enddo
       !
       ! determine maximum error to check convergence
       !
       do k = 1, kmax
          !
          do m = msta+1, mend, 2
             !
             nm = kgrpnt(m)
             !
             res = abs(u1(nm,k) - ui(nm,k))
             if ( res > resm ) then
                resm = res
                minf = m
                kinf = k
             endif
             u1(nm,k) = ui(nm,k)
             !
          enddo
          !
       enddo
       !
       ! exchange u-velocities with neighbouring subdomains
       !
       call SWEXCHG ( u1, kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
       ! for all even points:
       !
       ! compute right hand side
       !
       do k = 1, kmax
          !
          do m = msta, mend, 2
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             nmu = kgrpnt(mu)
             !
             wrk(nm,k) = rhsu(nm,k) - bux(nm,k)*u1(nmu,k) - bdx(nm,k)*u1(nmd,k)
             !
          enddo
          !
       enddo
       !
       ! solve the equation in vertical direction using Gaussian elimination
       !
       do m = msta, mend, 2
          !
          nm = kgrpnt(m)
          !
          ui(nm,1) = wrk(nm,1)
          !
       enddo
       !
       do k = 2, kmax
          !
          do m = msta, mend, 2
             !
             nm = kgrpnt(m)
             !
             ui(nm,k) = wrk(nm,k) - amatu(nm,k,2)*ui(nm,k-1)
             !
          enddo
          !
       enddo
       !
       do m = msta, mend, 2
          !
          nm = kgrpnt(m)
          !
          ui(nm,kmax) = ui(nm,kmax) / amatu(nm,kmax,1)
          !
       enddo
       !
       do k = kmax-1, 1, -1
          !
          do m = msta, mend, 2
             !
             nm = kgrpnt(m)
             !
             ui(nm,k) = ( ui(nm,k) - amatu(nm,k,3)*ui(nm,k+1) ) / amatu(nm,k,1)
             !
          enddo
          !
       enddo
       !
       ! determine maximum error to check convergence
       !
       do k = 1, kmax
          !
          do m = msta, mend, 2
             !
             nm = kgrpnt(m)
             !
             res = abs(u1(nm,k) - ui(nm,k))
             if ( res > resm ) then
                resm = res
                minf = m
                kinf = k
             endif
             u1(nm,k) = ui(nm,k)
             !
          enddo
          !
       enddo
       !
       ! exchange u-velocities with neighbouring subdomains
       !
       call SWEXCHG ( u1, kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
       rval = resm
       call SWREDUCE ( resm, 1, SWREAL, SWMAX )
       !
       if ( ITEST >= 30 .and. EQREAL(resm,rval) ) then
          !
          write (PRINTF,'(a,i3,a,e12.6,a,i5,a,i3)') ' ++ rb Jacobi: iter = ',jj,' res = ',resm,' in m=',minf+MXF-2,', layer=',kinf
          !
       endif
       !
       goto 10
       !
    endif
    !
    ! do a check for convergence
    !
    if ( jj >= maxit ) call msgerr (2, 'no convergence in red-black Jacobi iteration')
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
    ui = u1
    !
    ! determine maximum of water level and required accuracy
    !
    s0mx   = maxval(abs(s0))
    reps   = pnums(58)
    epslin = max(reps,reps*s0mx)
    call SWREDUCE ( epslin, 1, SWREAL, SWMAX )
    !
    ! initialize
    !
    maxit = nint(pnums(59))
    jj    = 0
    resm  = epslin + 1.
    !
    ! start iteration process to obtain 2nd order accuracy in pressure projection method
    !
 20 if ( resm > epslin .and. jj < maxit ) then
       !
       jj = jj + 1
       !
       resm = 0.
       minf = 0
       !
       ! compute the water level
       !
       fac = grav * dt * dt * rdx * rdx
       !
       ! first, build the equation for water level correction
       !
       if ( inewt == 0 ) then
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
             if ( presp(nm) == 1 ) then
                !
                fac1 = pship(1)
                !
             else
                !
                fac1 = 1.
                !
             endif
             !
             if ( m /= mfu .or. .not.LMXF ) then
                !
                a(nm) = -fac * teta(nmd) * hu(nmd)
                !
             else
                !
                a(nm) = 0.
                !
             endif
             !
             if ( m /= ml .or. .not.LMXL ) then
                !
                c(nm) = -fac * teta(nm) * hu(nm)
                !
             else
                !
                c(nm) = 0.
                !
             endif
             !
             b(nm) = fac1 - teta2(nm) * ( a(nm) + c(nm) )
             a(nm) = teta2(nmd) * a(nm)
             c(nm) = teta2(nmu) * c(nm)
             !
             fac1 = 0.
             fac2 = 0.
             !
             do k = 1, kmax
                !
                fac1 = fac1 + hku(nmd,k)* ( teta(nmd)*ui(nmd,k) + (1.-teta(nmd))*u0(nmd,k) )
                fac2 = fac2 + hku(nm ,k)* ( teta(nm )*ui(nm ,k) + (1.-teta(nm ))*u0(nm ,k) )
                !
             enddo
             !
             d(nm) = dt * rdx * ( fac1 - fac2 )
             !
             ! add contribution from rigid body motions
             !
             d(nm) = d(nm) - dt * skc(nm)
             !
          enddo
          !
          ! add mass source due to internal wave generation
          !
          if ( iwvgen /= 0 ) then
             !
             do m = mfu, ml
                !
                nm = kgrpnt(m)
                !
                d(nm) = d(nm) + dt * srcm(nm)
                !
             enddo
             !
          endif
          !
          ! add non-hydrostatic pressure part in the equation for water level correction (not piezometric head), if appropriate
          !
          if ( lpproj ) then
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                !
                if ( presp(nm) == 0 ) then
                   !
                   fac1 = 0.
                   fac2 = 0.
                   !
                   do k = 1, kmax
                      !
                      fac1 = fac1 + hku(nmd,k) * dqgrd(nmd,k)
                      fac2 = fac2 + hku(nm ,k) * dqgrd(nm ,k)
                      !
                   enddo
                   !
                   d(nm) = d(nm) + dt * dt * rdx * theta3 * ( teta(nm) * fac2 - teta(nmd) * fac1 )
                   !
                endif
                !
             enddo
             !
          endif
          !
          ! next, solve the equation for water level correction
          !
          call tridiag ( a, b, c, d, ds, kgrpnt )
          if (STPNOW()) return
          !
       else
          !
          lon = -dps  - s0
          upn = -flos - s0
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
             if ( presp(nm) == 1 ) then
                !
                fac1 = pship(1)
                !
             else
                !
                fac1 = 0.
                !
             endif
             !
             if ( m /= mfu .or. .not.LMXF ) then
                !
                an(nm) = -fac * teta(nmd) * hu(nmd)
                !
             else
                !
                an(nm) = 0.
                !
             endif
             !
             if ( m /= ml .or. .not.LMXL ) then
                !
                cn(nm) = -fac * teta(nm) * hu(nm)
                !
             else
                !
                cn(nm) = 0.
                !
             endif
             !
             bn(nm) = fac1 - teta2(nm) * ( an(nm) + cn(nm) )
             an(nm) = teta2(nmd) * an(nm)
             cn(nm) = teta2(nmu) * cn(nm)
             !
             fac1 = 0.
             fac2 = 0.
             !
             do k = 1, kmax
                !
                fac1 = fac1 + hku(nmd,k)* ( teta(nmd)*ui(nmd,k) + (1.-teta(nmd))*u0(nmd,k) )
                fac2 = fac2 + hku(nm ,k)* ( teta(nm )*ui(nm ,k) + (1.-teta(nm ))*u0(nm ,k) )
                !
             enddo
             !
             dn(nm) = dt * rdx * ( fac1 - fac2 )
             !
             ! add contribution from rigid body motions
             !
             dn(nm) = dn(nm) - dt * skc(nm)
             !
             dn(nm) = dn(nm) + max( lon(nm), min( upn(nm), 0.) )
             !
          enddo
          !
          ! add mass source due to internal wave generation
          !
          if ( iwvgen /= 0 ) then
             !
             do m = mfu, ml
                !
                nm = kgrpnt(m)
                !
                dn(nm) = dn(nm) + dt * srcm(nm)
                !
             enddo
             !
          endif
          !
          ! add non-hydrostatic pressure part in the equation for water level correction (not piezometric head), if appropriate
          !
          if ( lpproj ) then
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                !
                if ( presp(nm) == 0 ) then
                   !
                   fac1 = 0.
                   fac2 = 0.
                   !
                   do k = 1, kmax
                      !
                      fac1 = fac1 + hku(nmd,k) * dqgrd(nmd,k)
                      fac2 = fac2 + hku(nm ,k) * dqgrd(nm ,k)
                      !
                   enddo
                   !
                   dn(nm) = dn(nm) + dt * dt * rdx * theta3 * ( teta(nm) * fac2 - teta(nmd) * fac1 )
                   !
                endif
                !
             enddo
             !
          endif
          !
          ! next, solve the equation for water level correction
          !
          call newton1D ( an, bn, cn, dn, lon, upn, ds, kgrpnt )
          if (STPNOW()) return
          !
       endif
       !
       ! exchange water level corrections with neighbouring subdomains
       !
       call SWEXCHG ( ds, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! determine maximum error to check convergence projection pressure method
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m)
          !
          res = abs(s0(nm)+ds(nm) - s1(nm))
          if ( res > resm ) then
             resm = res
             minf = m
          endif
          !
       enddo
       !
       rval = resm
       call SWREDUCE ( resm, 1, SWREAL, SWMAX )
       !
       if ( iamout == 2 .and. EQREAL(resm,rval) .and. maxit > 1 ) then
          !
          write (PRINTF,'(a,i4,a,e12.6,a,i5)') ' ++ pressure projection: iter = ',jj,' res = ',resm,' in m=',minf+MXF-2
          !
       endif
       !
       ! correct horizontal velocity
       !
       msta = mf + 1
       mend = ml - 1
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
                u1(nm,k) = ui(nm,k) - grav * dt * rdx * ( teta2(nmu)*ds(nmu) - teta2(nm)*ds(nm) )
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! exchange u-velocities with neighbouring subdomains
       !
       call SWEXCHG ( u1, kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
       ! update water level
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m)
          !
          s1(nm) = s0(nm) + ds(nm)
          !
       enddo
       !
       ! check positivity of the water depth in each cell of the domain, if appropriate
       !
       if ( inewt == 0 ) then
          !
          cflmax = -999.
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             ! compute the sum of mass fluxes leaving the cell
             !
             sumqf = 0.
             !
             qf = 0.
             do k = 1, kmax
                qf = qf + hku(nm,k) * ( teta(nm)*u1(nm,k) + (1.-teta(nm))*u0(nm,k) )
             enddo
             !
             if ( qf > 0. ) sumqf = sumqf + qf
             !
             qf = 0.
             do k = 1, kmax
                qf = qf + hku(nmd,k) * ( teta(nmd)*u1(nmd,k) + (1.-teta(nmd))*u0(nmd,k) )
             enddo
             !
             if ( qf < 0. ) sumqf = sumqf - qf
             !
             ! compute the "flow" Courant number
             !
             if ( hs(nm) > epsdry ) then
                !
                cfl = sumqf * dt * rdx / hs(nm)
                if ( cfl > cflmax ) cflmax = cfl
                !
             endif
             !
          enddo
          !
          ! find maximum of CFL number over all subdomains
          !
          call SWREDUCE ( cflmax, 1, SWREAL, SWMAX )
          !
          ! give warning in case of CFL > 1
          !
          if ( .not. cflmax < 1. .and. INODE == MASTER ) then
             !
             call msgerr ( 1, 'positivity of the water depth cannot be guaranteed')
             call msgerr ( 0, 'It is advised to reduce the time step' )
             !
          endif
          !
       endif
       !
       ! copy to virtual cells at boundaries except for Riemann invariant openings
       !
       if ( ibl /= 6 .and. LMXF ) s1(nmf ) = s1(nmfu)
       if ( ibr /= 6 .and. LMXL ) s1(nmlu) = s1(nml )
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
                      rhsw(nm,k) = rhsw(nm,k) - 0.5 * ( fac1 * (w0(nm,k) - w0(nmd,k)) + fac2 * (w0(nmu,k) - w0(nm,k)) )
                      !
                   endif
                   !
                enddo
                !
             enddo
             !
             ! compute horizontal viscosity term (explicit)
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
          if ( verwinc ) then
             !
             ! compute vertical terms (implicit)
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
                      ! advection term
                      !
                      ctrkt = ( wom(nm,kd-1) + wom(nm,kd  ) ) / ( hks(nm,kd) + hks(nm,k+1) )
                      ctrkb = ( wom(nm,k   ) + wom(nm,k +1) ) / ( hks(nm,kd) + hks(nm,k+1) )
                      !
                      if ( k == 0 ) ctrkt = 0.
                      !
                      amatw(nm,k,1) = amatw(nm,k,1) - 0.5*ctrkt + 0.5*ctrkb
                      amatw(nm,k,2) =  0.5*ctrkt
                      amatw(nm,k,3) = -0.5*ctrkb
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
                   if ( presp(nm) == 0 ) then
                      !
                      ! free surface flow
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
                      ! pressurized flow
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
                   endif
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
             fac = -fac
             !
             do m = mfu, ml
                !
                nm = kgrpnt(m)
                !
                if ( wets(nm) == 1 ) then
                   !
                   if ( ihydro == 1 .and. presp(nm) == 0 ) then
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
                      !rhsw(nm,k) = rhsw(nm,k) - fac * ( w1(nm,kmax) - w0(nm,kmax) ) / (dt*thetaw)
                      !
                   else
                      !
                      kd = max(k,1)
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
          ! hull (pressurized flow):
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
             if ( wets(nm) == 1 .and. presp(nm) == 1 ) then
                !
                amatw(nm,0,1) = 1.
                amatw(nm,0,2) = 0.
                amatw(nm,0,3) = 0.
                w1   (nm,0  ) = 0.25 * rdx * ( 3.*(u1(nm,1)+u1(nmd,1)) - u1(nm,2) - u1(nmd,2) ) * ( zku(nm,0) - zku(nmd,0) )
                !
                ! add body motion
                w1  (nm,0) = w1(nm,0) + skc(nm)
                rhsw(nm,0) = w1(nm,0)
                !
             endif
             !
          enddo
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
             do m = mfu, ml
                !
                nm = kgrpnt(m)
                !
                if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
                   !
                   w1(nm,0) = ( w1(nm,0) - (1.-thetaw) * w0(nm,0) ) / thetaw
                   !
                endif
                !
             enddo
             !
             do k = 1, kmax-1
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
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             nmu = kgrpnt(mu)
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
                   if ( k == 1 .and. presp(nm) == 1 ) fac1 = 0.
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
                kwd = 1.
                if ( k == 1 .and. presp(nm) == 1 ) kwd = 0.
                !
                amatp(nm,k, 1) = dmat(nm,k,2) * gmatu(nm ,k ,2) + dmat(nm,k,4) * gmatu(nmd,k ,4) + dmat(nm,k,1) * gmatu(nm ,kd,5) +  &
                                 dmat(nm,k,3) * gmatu(nmd,kd,6) + dmat(nm,k,5) * gmatu(nm ,ku,1) + dmat(nm,k,6) * gmatu(nmd,ku,3)
                amatp(nm,k, 2) = dmat(nm,k,4) * gmatu(nmd,k ,2) + dmat(nm,k,3) * gmatu(nmd,kd,5) + dmat(nm,k,6) * gmatu(nmd,ku,1)
                amatp(nm,k, 3) = dmat(nm,k,2) * gmatu(nm ,k ,4) + dmat(nm,k,1) * gmatu(nm ,kd,6) + dmat(nm,k,5) * gmatu(nm ,ku,3)
                amatp(nm,k, 4) = dmat(nm,k,2) * gmatu(nm ,k ,1) + dmat(nm,k,4) * gmatu(nmd,k ,3) +  &
                                 dmat(nm,k,1) * gmatu(nm ,kd,2) + dmat(nm,k,3) * gmatu(nmd,kd,4)
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
                if ( ihydro == 1 .and. presp(nm) == 0 ) then
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
                   amatp(nm,k,4) = amatp(nm,k,4) + gmatw(nm,k,1)
                   !
                else
                   !
                   amatp(nm,k,1) = amatp(nm,k,1) + kwd * gmatw(nm,k,2) - kwu * gmatw(nm,ku,1)
                   amatp(nm,k,4) = amatp(nm,k,4) + kwd * gmatw(nm,k,1)
                   amatp(nm,k,5) = amatp(nm,k,5)                       - kwu * gmatw(nm,ku,2)
                   !
                endif
                !
                rhsp(nm,k) = ( dmat(nm,k,2) * u1(nm,k ) + dmat(nm,k,4) * u1(nmd,k ) +  &
                               dmat(nm,k,1) * u1(nm,kd) + dmat(nm,k,3) * u1(nmd,kd) +  &
                               dmat(nm,k,5) * u1(nm,ku) + dmat(nm,k,6) * u1(nmd,ku) +  &
                               kwd*w1(nm,k-1) - kwu*w1(nm,k) ) / (dt*theta3)
                !
                ! add contribution from rigid body motions (top layer only)
                !
                if ( k == 1 ) rhsp(nm,k) = rhsp(nm,k) + skc(nm) / (dt*theta3)
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
!TIMG                call SWTSTA(296)
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
!TIMG                call SWTSTO(296)
                !
             endif
             !
!TIMG             call SWTSTA(297)
             call bicgstab ( amatp(1:mcgrd,1:qmax,1:nconct), rhsp(1:mcgrd,1:qmax), dq(1:mcgrd,1:qmax) )
!TIMG             call SWTSTO(297)
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
       ! update pressure gradient in case of pressure projection
       !
       if ( lpproj ) then
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
                if ( wetu(nm) == 1 .and. presu(nm) == 0 ) then
                   !
                   dqgrd(nm,k) = gmatu(nm,k,1)*dq(nm ,kd) + gmatu(nm,k,2)*dq(nm,k ) + gmatu(nm,k,3)*dq(nmu,kd)  &
                               + gmatu(nm,k,4)*dq(nmu,k ) + gmatu(nm,k,5)*dq(nm,ku) + gmatu(nm,k,6)*dq(nmu,ku)
                   !
                else
                   ! no correction applied on piezometric head
                   !
                   dqgrd(nm,k) = 0.
                   !
                endif
                !
             enddo
             !
          enddo
          !
       endif
       !
       goto 20
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
                u1(nm,k) = u1(nm,k) - dt*theta3*( gmatu(nm,k,1)*dq(nm ,kd) + gmatu(nm,k,2)*dq(nm,k ) + gmatu(nm,k,3)*dq(nmu,kd) + &
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
          do m = mfu, ml
             !
             nm = kgrpnt(m)
             !
             if ( wets(nm) == 1 ) then
                !
                if ( ihydro == 1 .and. presp(nm) == 0 ) then
                   !
                   do j = 0, kmax-1
                      !
                      kd = max(k+j  ,1   )
                      kd = min(kd   ,kmax)
                      ku = min(k+j+1,kmax)
                      !
                      w1(nm,k) = w1(nm,k) - dt*theta3*( gmatw(nm,k+1,2*j+1)*dq(nm,kd) + gmatw(nm,k+1,2*j+2)*dq(nm,ku) )
                      !
                   enddo
                   !
                else
                   !
                   kd = max(k,1)
                   !
                   w1(nm,k) = w1(nm,k) - dt*theta3*( gmatw(nm,k+1,1)*dq(nm,kd) + gmatw(nm,k+1,2)*dq(nm,k+1) )
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
          if ( wets(nm) == 1 .and. presp(nm) == 1 ) then
             !
             w1(nm,0) = 0.25 * rdx * ( 3.*(u1(nm,1)+u1(nmd,1)) - u1(nm,2) - u1(nmd,2) ) * ( zku(nm,0) - zku(nmd,0) )
             !
             ! add body motion
             w1(nm,0) = w1(nm,0) + skc(nm)
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute the mass flux
    !
    do k = 1, kmax
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          qx(nm,k) = hku(nm,k) * ( teta(nm)*u1(nm,k) + (1.-teta(nm))*u0(nm,k) )
          !
       enddo
       !
    enddo
    !
    ! compute the depth-averaged u-velocity
    !
    udep = 0.
    !
    do m = mf, ml
       !
       nm = kgrpnt(m)
       !
       if ( wetu(nm) == 1 ) then
          !
          do k = 1, kmax
             !
             udep(nm) = udep(nm) + hku(nm,k)*u1(nm,k)
             !
          enddo
          !
          udep(nm) = udep(nm) / hu(nm)
          !
       endif
       !
    enddo
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
    ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
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
                kwd = 1.
                if ( k == 1 .and. presp(nm) == 1 ) kwd = 0.
                !
                moutf = moutf + dmat(nm,k,2) * u1(nm,k ) + dmat(nm,k,4) * u1(nmd,k ) +  &
                                dmat(nm,k,1) * u1(nm,kd) + dmat(nm,k,3) * u1(nmd,kd) +  &
                                dmat(nm,k,5) * u1(nm,ku) + dmat(nm,k,6) * u1(nmd,ku) +  &
                                kwd*w1(nm,k-1) - kwu*w1(nm,k)
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
          if ( presp(nm) == 0 ) vol = vol + s1(nm)
          !
          ! compute potential energy in the cell center
          ener = ener + 0.5 * grav * dx * s1(nm) * s1(nm)
          !
          ! compute kinetic energy in the cell center
          fac1 = ( min( -flos(nm ), s1(nm ) ) + dps(nm ) ) * dx
          fac2 = ( min( -flos(nmd), s1(nmd) ) + dps(nmd) ) * dx
          ener = ener + 0.125 * ( fac1 + fac2 ) * udep(nmd) * udep(nmd)
          fac2 = ( min( -flos(nmu), s1(nmu) ) + dps(nmu) ) * dx
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
    ! adapt water depth to include floating object
    !
    if ( ifloat /= 0 ) work(:,1) = min( dps(:)-flos(:), work(:,1) )
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
                wom(nm,k) = wom(nm,k+1) - rdx * (      teta(nm)*hku(nm,k+1)*u1(nm,k+1) -      teta(nmd)*hku(nmd,k+1)*u1(nmd,k+1) )  &
                                        - rdx * ( (1.-teta(nm))*hku(nm,k+1)*u0(nm,k+1) - (1.-teta(nmd))*hku(nmd,k+1)*u0(nmd,k+1) )  &
                                        - ( hksnew(nm,k+1) - hks(nm,k+1) ) / dt
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
                wom(nm,k) = w1(nm,k) - ( zksnew(nm,k) - zks(nm,k) ) / dt - fac * rdx *( zku(nm,k) - zku(nmd,k) )
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
                wom(nm,k) = wom(nm,k+1) - rdx * (      teta(nm)*hku(nm,k+1)*u1(nm,k+1) -      teta(nmd)*hku(nmd,k+1)*u1(nmd,k+1) )  &
                                        - rdx * ( (1.-teta(nm))*hku(nm,k+1)*u0(nm,k+1) - (1.-teta(nmd))*hku(nmd,k+1)*u0(nmd,k+1) )  &
                                        - ( hksnew(nm,k+1) - hks(nm,k+1) ) / dt
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
    ! compute the maximum CFL number
    !
    cflmax = -999.
    !
    do k = 1, kmax
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          cfl = rdx * dt * abs(u1(nm,k))
          if ( cfl > cflmax ) cflmax = cfl
          !
       enddo
       !
    enddo
    !
    ! find maximum of CFL number over all subdomains
    !
    call SWREDUCE ( cflmax, 1, SWREAL, SWMAX )
    !
    ! give warning in case of CFL > 1
    !
    if ( .not. cflmax < 1. .and. INODE == MASTER ) then
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
end subroutine SwashImpLayM1DHflow
