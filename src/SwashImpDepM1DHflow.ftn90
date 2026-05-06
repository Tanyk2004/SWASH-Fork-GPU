subroutine SwashImpDepM1DHflow ( u1, u0, qx, qm, q, dq, gmat, rho, ui, dqgrd, ibl, ibr, kgrpnt )
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
!   Performs the time integration for the non-hydrostatic, depth-averaged 1D shallow water equations
!
!   Note: this subroutine solves the depth-averaged 1D shallow water equations, which are strictly mass and
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
!   The time integration with respect to the advective and viscosity terms is based on the
!   Crank-Nicolson scheme, while that for the bottom friction is based on Euler implicit and for
!   the non-hydrostatic pressure gradient a semi-implicit approach is employed (theta-scheme).
!
!   The space discretization of the horizontal advective and viscosity terms is strictly momentum
!   conservative and are approximated with central differences, so that the advection term is
!   skew-symmetric and the viscosity term is symmetric.
!
!   The w-momentum equation only contains the z-gradient of the non-hydrostatic pressure and is
!   discretized by means of the Keller-box scheme.
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
    use m_parall
    use SwashSolvedata, only: iamout
    use SwashFlowdata, ibltmp   => ibl  , &
                       ibrtmp   => ibr  , &
                       rhotmp   => rho  , &
                       u1tmp    => u1   , &
                       u0tmp    => u0   , &
                       uitmp    => ui   , &
                       qxtmp    => qx   , &
                       qmtmp    => qm   , &
                       qtmp     => q    , &
                       dqtmp    => dq   , &
                       dqgrdtmp => dqgrd
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                        :: ibl    ! boundary condition type at left boundary
    integer, intent(in)                        :: ibr    ! boundary condition type at right boundary
    !
    integer, dimension(mxc)    , intent(in)    :: kgrpnt ! index table containing the address of each (active) grid point
                                                         ! =1: not active grid point
                                                         ! >1: active grid point
    !
    real   , dimension(mcgrd)  , intent(out)   :: dq     ! non-hydrostatic pressure correction
    real   , dimension(mcgrd)  , intent(inout) :: dqgrd  ! gradient of pressure correction to be employed in pressure projection method
    real   , dimension(mcgrd,2), intent(out)   :: gmat   ! gradient matrix for pressure in u-point
    real   , dimension(mcgrd)  , intent(inout) :: q      ! non-hydrostatic pressure
    real   , dimension(mcgrd)  , intent(out)   :: qm     ! mass flux in wl-point based on averaging
    real   , dimension(mcgrd)  , intent(out)   :: qx     ! mass flux in x-direction
    real   , dimension(mcgrd)  , intent(in)    :: rho    ! density of water
    real   , dimension(mcgrd)  , intent(in)    :: u0     ! u-velocity at previous time level
    real   , dimension(mcgrd)  , intent(inout) :: u1     ! u-velocity at current time level
    real   , dimension(mcgrd)  , intent(out)   :: ui     ! intermediate u-velocity in iterative process of pressure projection method
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! iteration counter
    integer       :: m        ! loop counter
    integer       :: maxit    ! maximum number of iterations
    integer       :: md       ! index of point m-1
    integer       :: mdd      ! index of point m-2
    integer       :: mend     ! end index of loop over u-points
    integer       :: minf     ! m-index of point with largest error in solution
    integer       :: msta     ! start index of loop over u-points
    integer       :: mu       ! index of point m+1
    integer       :: muu      ! index of point m+2
    integer       :: nm       ! pointer to m
    integer       :: nmd      ! pointer to m-1
    integer       :: nmdd     ! pointer to m-2
    integer       :: nmf      ! pointer to mf
    integer       :: nmfu     ! pointer to mfu
    integer       :: nmfuu    ! pointer to mfuu
    integer       :: nml      ! pointer to ml
    integer       :: nmld     ! pointer to mld
    integer       :: nmlu     ! pointer to mlu
    integer       :: nmu      ! pointer to m+1
    integer       :: nmuu     ! pointer to m+2
    !
    real          :: cfl      ! CFL number
    real          :: contrib  ! total contribution of momentum equation
    real          :: cporfr   ! flow resistance terms inside porous medium
    real          :: cveg     ! vegetation term
    real          :: denom    ! a denominator
    real          :: ener     ! total energy of closed system
    real          :: epslin   ! required accuracy in the linear solver
    real          :: fac      ! a factor
    real          :: fac1     ! another factor
    real          :: fac2     ! some other factor
    real          :: moutf    ! net mass outflow
    real          :: qf       ! updated mass flux at present face
    real          :: rdx      ! reciprocal of mesh size
    real          :: reps     ! accuracy of the final approximation
    real          :: res      ! residual
    real          :: resm     ! maximum error
    real          :: rval     ! auxiliary real
    real          :: s0mx     ! maximum value of water level
    real          :: sumqf    ! sum of outgoing mass fluxes over present cell
    real          :: theta3   ! implicitness factor for non-hydrostatic pressure gradient
    real          :: vol      ! total displaced volume of water
    real          :: zgrad    ! water level gradient term
    !
    logical       :: EQREAL   ! compares two reals
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashImpDepM1DHflow')
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
    ! compute the mass flux
    !
    do m = mf, ml
       !
       nm = kgrpnt(m)
       !
       qx(nm) = hu(nm)*u0(nm)
       !
    enddo
    !
    ! compute the mass flux in wl-point based on averaging
    !
    do m = mfu, ml
       !
       md = m - 1
       !
       nm  = kgrpnt(m )
       nmd = kgrpnt(md)
       !
       qm(nm) = 0.5 * ( qx(nm) + qx(nmd) )
       !
    enddo
    !
    ! build gradient matrix for non-hydrostatic pressure
    !
    if ( ihydro == 1 ) then
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
             fac = 0.5 * rdx / hum(nm)
             !
             if ( presu(nm) == 0 ) then
                !
                ! free surface flow
                !
                gmat(nm,1) = -fac * ( s0(nm ) + dps(nmu) )
                gmat(nm,2) =  fac * ( s0(nmu) + dps(nm ) )
                !
             else
                !
                ! pressurized flow
                !
                if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
                   !
                   ! changeover from free surface to pressurized
                   !
                   gmat(nm,1) = -0.5 * fac * ( hs(nm) + hs(nmu) )
                   gmat(nm,2) =        fac * ( hs(nm) + hs(nmu) )
                   !
                else if ( presp(nm) == 1 .and. presp(nmu) == 0 ) then
                   !
                   ! changeover from pressurized to free surface
                   !
                   gmat(nm,1) =       -fac * ( hs(nm) + hs(nmu) )
                   gmat(nm,2) =  0.5 * fac * ( hs(nm) + hs(nmu) )
                   !
                else
                   !
                   gmat(nm,1) = -fac * ( hs(nm) + hs(nmu) )
                   gmat(nm,2) =  fac * ( hs(nm) + hs(nmu) )
                   !
                endif
                !
             endif
             !
          else
             !
             gmat(nm,1) = 0.
             gmat(nm,2) = 0.
             !
          endif
          !
       enddo
       !
       if ( LMXF ) then
          !
          gmat(nmf,2) = gmat(nmf,2) - gmat(nmf,1)
          gmat(nmf,1) = 0.
          !
       endif
       !
       if ( LMXL ) then
          !
          gmat(nml,1) = gmat(nml,1) - gmat(nml,2)
          gmat(nml,2) = 0.
          !
       endif
       !
    endif
    !
    ! compute gradient of non-hydrostatic pressure in u-momentum equation
    !
    if ( ihydro == 1 ) then
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
             qgrad(nm) = gmat(nm,1) * q(nm) + gmat(nm,2) * q(nmu)
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
    ! compute implicit part of wind stress term, if appropriate
    !
    if ( relwnd ) then
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             wndimp(nm) = cwndu(nm) / max(1.e-3,hum(nm))
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
    ! compute bottom friction term, if appropriate
    !
    if ( irough == 11 ) then
       !
       ! linear bottom friction
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             cbot(nm) = cfricu(nm) / hum(nm)
             !
          endif
          !
       enddo
       !
    else if ( irough /= 0 ) then
       !
       ! quadratic bottom friction
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             cbot(nm) = cfricu(nm) * abs(u0(nm)) / hum(nm)
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
    ! compute baroclinic forcing
    !
    if ( idens /= 0 ) then
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
             pgrad(nm) = 0.5 * grav * hum(nm) * rdx * (rho(nmu) - rho(nm)) / rhow
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
       do m = mf, ml
          !
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmu = kgrpnt(mu)
          !
          if ( wetu(nm) == 1 ) then
             !
             pgrad(nm) = pgrad(nm) + rdx * (patm(nmu) - patm(nm)) / rhow
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute intermediate u-velocity
    !
    msta = mf + 1   ! first internal u-point
    mend = ml - 1   ! last  internal u-point
    !
    do m = msta, mend
       !
       md = m - 1
       mu = m + 1
       !
       nm  = kgrpnt(m )
       nmd = kgrpnt(md)
       nmu = kgrpnt(mu)
       !
       if ( humn(nm) > 0. ) then
          !
          fac1 = 0.25 * dt * rdx * qm(nmu) / humn(nm)
          fac2 = 0.25 * dt * rdx * qm(nm ) / humn(nm)
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
          ! compute water level gradient
          !
          zgrad = grav * rdx * (s0(nmu) - s0(nm))
          !
          ! compute flow resistance inside porous medium
          !
          cporfr = apomu(nm,1) + bpomu(nm,1) * abs(u0(nm))
          !
          ! compute friction induced by vegetation
          !
          cveg = cvegu(nm,1,1) * abs(u0(nm))
          !
          ! compute total contributions of the momentum equation (except advection and viscosity)
          !
          contrib = zgrad + pgrad(nm) + qgrad(nm) - windu(nm)/max(1.e-3,hum(nm))
          !
          denom = 1. + cvegu(nm,1,2) + cpomu(nm,1) + dt * ( cbot(nm) + cporfr + cveg + wndimp(nm) )
          !
          ! compute advection term using central differences (momentum conservative and skew-symmetric) and build system of equations
          !
          a(nm) = -fac2
          c(nm) =  fac1
          b(nm) =  denom - a(nm) - c(nm)
          d(nm) = (1.+cvegu(nm,1,2)+cpomu(nm,1)+a(nm)+c(nm)) * u0(nm) - a(nm) * u0(nmd) - c(nm) * u0(nmu) - dt * contrib
          !
       else
          !
          a(nm) = 0.
          b(nm) = 1.
          c(nm) = 0.
          d(nm) = 0.
          !
       endif
       !
    enddo
    !
    ! add viscosity term, if appropriate
    !
    if ( ihvisc == 1 .and. hvisc > 0. ) then
       !
       do m = msta, mend
          !
          md = m - 1
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          nmu = kgrpnt(mu)
          !
          if ( hum(nm) > 0. ) then
             !
             fac1 = 0.5 * dt * rdx * rdx * hvisc * hs(nmu) / hum(nm)
             fac2 = 0.5 * dt * rdx * rdx * hvisc * hs(nm ) / hum(nm)
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
             ! build system of equations
             !
             a(nm) = a(nm) - fac2
             c(nm) = c(nm) - fac1
             b(nm) = b(nm) + fac1 + fac2
             d(nm) = d(nm) - (fac1+fac2) * u0(nm) + fac2 * u0(nmd) + fac1 * u0(nmu)
             !
          endif
          !
       enddo
       !
    else if ( ihvisc > 1 ) then
       !
       do m = msta, mend
          !
          md = m - 1
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          nmu = kgrpnt(mu)
          !
          if ( hum(nm) > 0. ) then
             !
             fac1 = 0.5 * dt * rdx * rdx * vnu2d(nmu) * hs(nmu) / hum(nm)
             fac2 = 0.5 * dt * rdx * rdx * vnu2d(nm ) * hs(nm ) / hum(nm)
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
             ! build system of equations
             !
             a(nm) = a(nm) - fac2
             c(nm) = c(nm) - fac1
             b(nm) = b(nm) + fac1 + fac2
             d(nm) = d(nm) - (fac1+fac2) * u0(nm) + fac2 * u0(nmd) + fac1 * u0(nmu)
             !
          endif
          !
       enddo
       !
    endif
    !
    ! incorporate boundary conditions
    !
    if ( LMXF ) then
       !
       d(nmfu) = d(nmfu) - a(nmfu) * u1(nmf)
       a(nmfu) = 0.
       !
    endif
    !
    if ( LMXL ) then
       !
       a(nml) = 0.
       b(nml) = 1.
       c(nml) = 0.
       d(nml) = u1(nml)
       !
    endif
    !
    ! solve system of equations
    !
    call tridiag ( a, b, c, d, u1, kgrpnt )
    if (STPNOW()) return
    !
    ! exchange u-velocities with neighbouring subdomains
    !
    call SWEXCHG ( u1, kgrpnt, 1, 1 )
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
    j     = 0
    resm  = epslin + 1.
    !
    ! start iteration process to obtain 2nd order accuracy in pressure projection method
    !
 10 if ( resm > epslin .and. j < maxit ) then
       !
       j = j + 1
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
             d(nm) = dt * rdx * ( hu(nmd) * ( teta(nmd)*ui(nmd) + (1.-teta(nmd))*u0(nmd) ) - hu(nm) * ( teta(nm)*ui(nm) + (1.-teta(nm))*u0(nm) ) )
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
                   d(nm) = d(nm) + dt * dt * rdx * theta3 * ( teta(nm) * hu(nm) * dqgrd(nm) - teta(nmd) * hu(nmd) * dqgrd(nmd) )
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
             dn(nm) = dt * rdx * ( hu(nmd) * ( teta(nmd)*ui(nmd) + (1.-teta(nmd))*u0(nmd) ) - hu(nm) * ( teta(nm)*ui(nm) + (1.-teta(nm))*u0(nm) ) )
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
                   dn(nm) = dn(nm) + dt * dt * rdx * theta3 * ( teta(nm) * hu(nm) * dqgrd(nm) - teta(nmd) * hu(nmd) * dqgrd(nmd) )
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
          write (PRINTF,'(a,i4,a,e12.6,a,i5)') ' ++ pressure projection: iter = ',j,' res = ',resm,' in m=',minf+MXF-2
          !
       endif
       !
       ! correct horizontal velocity
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
             u1(nm) = ui(nm) - grav * dt * rdx * ( teta2(nmu)*ds(nmu) - teta2(nm)*ds(nm) )
             !
          endif
          !
       enddo
       !
       ! exchange u-velocities with neighbouring subdomains
       !
       call SWEXCHG ( u1, kgrpnt, 1, 1 )
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
             qf = hu(nm) * ( teta(nm)*u1(nm) + (1.-teta(nm))*u0(nm) )
             !
             if ( qf > 0. ) sumqf = sumqf + qf
             !
             qf = hu(nmd) * ( teta(nmd)*u1(nmd) + (1.-teta(nmd))*u0(nmd) )
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
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             if ( wets(nm) == 1 ) then
                !
                w1bot(nm) = -0.5 * rdx * ( u1(nm) + u1(nmd) ) * ( dpu(nm) - dpu(nmd) )
                !
                if ( presp(nm) == 0 ) then
                   !
                   ! free surface flow
                   !
                   w1top(nm) = w0top(nm) + w0bot(nm) - w1bot(nm) + 2.*dt*fac*q(nm)/hs(nm)
                   !
                else
                   !
                   ! pressurized flow
                   !
                   w1top(nm) = -0.5 * rdx * ( u1(nm) + u1(nmd) ) * ( flou(nm) - flou(nmd) )
                   !
                   ! add body motion
                   w1top(nm) = w1top(nm) + skc(nm)
                   !
                endif
                !
             else
                !
                w1bot(nm) = 0.
                w1top(nm) = 0.
                !
             endif
             !
          enddo
          !
          if ( horwinc ) then
             !
             ! compute horizontal advection term (momentum conservative and explicit)
             !
             fac = 2. * dt * rdx
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
                if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
                   !
                   fac1 = fac * qx(nmd) / hs(nm)
                   fac2 = fac * qx(nm ) / hs(nm)
                   !
                   w1top(nm) = w1top(nm) - 0.5 * ( fac1 * (w0top(nm) - w0top(nmd)) + fac2 * (w0top(nmu) - w0top(nm)) )
                   !
                endif
                !
             enddo
             !
          endif
          !
          ! add vegetation dissipation due to horizontal cylinders
          !
          if ( ivegw /= 0 ) then
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m )
                nmd = kgrpnt(md)
                !
                if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
                   !
                   ! compute vegetation factor
                   !
                   cveg = 0.25 * ( cvegu(nm,1,1) + cvegu(nmd,1,1) ) * dt * abs( w0top(nm) + w0bot(nm) )
                   !
                   w1top(nm) = w1top(nm) - cveg * w1bot(nm)
                   w1top(nm) = w1top(nm) / ( 1. + cveg )
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
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
                !
                fac1 = rdx * ( dpu(nm) - dpu(nmd) + 0.5 * ( hu(nm) + hu(nmd) ) )
                fac2 = rdx * ( dpu(nm) - dpu(nmd) - 0.5 * ( hu(nm) + hu(nmd) ) )
                !
                a(nm) = fac2 * gmat(nmd,1)
                b(nm) = fac1 * gmat(nm ,1) + fac2 * gmat(nmd,2) - 2./hs(nm)
                c(nm) = fac1 * gmat(nm ,2)
                d(nm) = ( fac1*u1(nm) + fac2*u1(nmd) + w1top(nm) + w1bot(nm) )/(dt*theta3)
                !
             else if ( wets(nm) == 1 .and. presp(nm) == 1 ) then
                !
                fac1 =  hu(nm )
                fac2 = -hu(nmd)
                !
                a(nm) = fac2 * gmat(nmd,1)
                b(nm) = fac1 * gmat(nm ,1) + fac2 * gmat(nmd,2)
                c(nm) = fac1 * gmat(nm ,2)
                d(nm) = ( fac1*u1(nm) + fac2*u1(nmd) )/(dt*theta3)
                !
                ! add contribution from rigid body motions
                !
                d(nm) = d(nm) + dx * skc(nm) / (dt*theta3)
                !
             else
                !
                a(nm) = 0.
                b(nm) = 1.
                c(nm) = 0.
                d(nm) = 0.
                !
             endif
             !
          enddo
          !
          ! solve the Poisson equation
          !
          call tridiag ( a, b, c, d, dq, kgrpnt )
          if (STPNOW()) return
          !
          ! exchange pressure corrections with neighbouring subdomains
          !
          call SWEXCHG ( dq, kgrpnt, 1, 1 )
          if (STPNOW()) return
          !
       endif
       !
       ! update pressure gradient in case of pressure projection
       !
       if ( lpproj ) then
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
                dqgrd(nm) = gmat(nm,1) * dq(nm) + gmat(nm,2) * dq(nmu)
                !
             else
                ! no correction applied on piezometric head
                !
                dqgrd(nm) = 0.
                !
             endif
             !
          enddo
          !
       endif
       !
       goto 10
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
       do m = mf, ml
          !
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmu = kgrpnt(mu)
          !
          if ( wetu(nm) == 1 ) then
             !
             u1(nm) = u1(nm) - dt*theta3*( gmat(nm,1) * dq(nm) + gmat(nm,2) * dq(nmu) )
             !
          endif
          !
       enddo
       !
       ! w-velocity
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
             if ( presp(nm) == 0 ) w1top(nm) = w1top(nm) + w1bot(nm) + 2.*dt*theta3*dq(nm)/hs(nm)
             w1bot(nm) = -0.5 * rdx * ( u1(nm) + u1(nmd) ) * ( dpu(nm) - dpu(nmd) )
             if ( presp(nm) == 0 ) then
                w1top(nm) = w1top(nm) - w1bot(nm)
             else
                w1top(nm) = -0.5 * rdx * ( u1(nm) + u1(nmd) ) * ( flou(nm) - flou(nmd) )
                ! add body motion
                w1top(nm) = w1top(nm) + skc(nm)
             endif
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute the mass flux
    !
    do m = mf, ml
       !
       nm = kgrpnt(m)
       !
       qx(nm) = hu(nm) * ( teta(nm)*u1(nm) + (1.-teta(nm))*u0(nm) )
       !
    enddo
    !
    ! impose Neumann condition for w-velocity at boundaries, if appropriate
    !
    if ( horwinc .and. ihydro == 1 ) then
       !
       if ( LMXF ) w1top(nmf ) = w1top(nmfu)
       if ( LMXL ) w1top(nmlu) = w1top(nml )
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
    if ( ITEST >= 30 .and. ihydro == 1 ) then
       !
       moutf = 0.
       !
       do m = mfu, mend
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          moutf = moutf + qx(nm) - qx(nmd) - 0.5 * ( u1(nm) + u1(nmd) ) * ( hu(nm) - hu(nmd) ) + dx * ( w1top(nm) - w1bot(nm) )
          !
       enddo
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
          ener = ener + 0.125 * ( fac1 + fac2 ) * u1(nmd) * u1(nmd)
          fac2 = ( min( -flos(nmu), s1(nmu) ) + dps(nmu) ) * dx
          ener = ener + 0.125 * ( fac1 + fac2 ) * u1(nm ) * u1(nm )
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
    ! exchange u-velocities with neighbouring subdomains (if appropriate)
    !
    if ( ihydro == 1 ) then
       !
       call SWEXCHG ( u1, kgrpnt, 1, 1 )
       if ( horwinc ) call SWEXCHG ( w1top, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    endif
    !
    ! apply wave absorption by means of sponge layer, if appropriate
    !
    if ( spwidl > 0. ) then
       !
       tbndx = 0.
       if ( LMXF ) then
          tbndx(1,1) = u1(nmf )
          tbndx(1,2) = s1(nmfu)
       endif
       call SWREDUCE ( tbndx, 2, SWREAL, SWSUM )
       if (STPNOW()) return
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             u1(nm) = (1.-sponxl(nm))*u1(nm) + sponxl(nm)*tbndx(1,1)
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
             s1(nm) = (1.-sponxl(nmd))*s1(nm) + sponxl(nmd)*tbndx(1,2)
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
          tbndx(1,1) = u1(nml)
          tbndx(1,2) = s1(nml)
       endif
       call SWREDUCE ( tbndx, 2, SWREAL, SWSUM )
       if (STPNOW()) return
       !
       do m = mf, ml
          !
          nm = kgrpnt(m)
          !
          if ( wetu(nm) == 1 ) then
             !
             u1(nm) = (1.-sponxr(nm))*u1(nm) + sponxr(nm)*tbndx(1,1)
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
             s1(nm) = (1.-sponxr(nm))*s1(nm) + sponxr(nm)*tbndx(1,2)
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
    do m = mf, ml
       !
       nm = kgrpnt(m)
       !
       cfl = rdx * dt * abs(u1(nm))
       if ( cfl > cflmax ) cflmax = cfl
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
end subroutine SwashImpDepM1DHflow
