subroutine SwashImpDepM2DHflow ( u1, u0, qx, qxm, v1, v0, qy, qym, q, dq, gmatu, gmatv, rho, ui, dqgrdu, vi, dqgrdv )
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
!   Performs the time integration for the non-hydrostatic, depth-averaged 2D shallow water equations
!
!   Note: this subroutine solves the depth-averaged 2D shallow water equations, which are strictly mass and
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
!   momentum equations is based on a theta-scheme. Only a value of theta = 0.5 is taken.
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
    use m_genarr, only: kgrpnt, guu, guv, gvu, gvv, gsqs, gsqsu, gsqsv
    use m_parall
    use SwashSolvedata, only: iamout
    use SwashFlowdata, rhotmp    => rho   , &
                       u1tmp     => u1    , &
                       u0tmp     => u0    , &
                       uitmp     => ui    , &
                       qxtmp     => qx    , &
                       qxmtmp    => qxm   , &
                       v1tmp     => v1    , &
                       v0tmp     => v0    , &
                       vitmp     => vi    , &
                       qytmp     => qy    , &
                       qymtmp    => qym   , &
                       qtmp      => q     , &
                       dqtmp     => dq    , &
                       dqgrdutmp => dqgrdu, &
                       dqgrdvtmp => dqgrdv, &
                       gmatutmp  => gmatu , &
                       gmatvtmp  => gmatv
!
    implicit none
!
!   Argument variables
!
    real   , dimension(mcgrd)  , intent(out)   :: dq       ! non-hydrostatic pressure correction
    real   , dimension(mcgrd)  , intent(inout) :: dqgrdu   ! gradient of pressure correction in x-direction to be employed in pressure projection method
    real   , dimension(mcgrd)  , intent(inout) :: dqgrdv   ! gradient of pressure correction in y-direction to be employed in pressure projection method
    real   , dimension(mcgrd,2), intent(out)   :: gmatu    ! gradient matrix for pressure in u-point
    real   , dimension(mcgrd,2), intent(out)   :: gmatv    ! gradient matrix for pressure in v-point
    real   , dimension(mcgrd)  , intent(inout) :: q        ! non-hydrostatic pressure
    real   , dimension(mcgrd)  , intent(out)   :: qx       ! mass flux per width in x-direction
    real   , dimension(mcgrd)  , intent(out)   :: qxm      ! mass flux per width in x-direction based on averaging
    real   , dimension(mcgrd)  , intent(out)   :: qy       ! mass flux per width in y-direction
    real   , dimension(mcgrd)  , intent(out)   :: qym      ! mass flux per width in y-direction based on averaging
    real   , dimension(mcgrd)  , intent(in)    :: rho      ! density of water
    real   , dimension(mcgrd)  , intent(in)    :: u0       ! u-velocity at previous time level
    real   , dimension(mcgrd)  , intent(inout) :: u1       ! u-velocity at current time level
    real   , dimension(mcgrd)  , intent(out)   :: ui       ! intermediate u-velocity in iterative process of pressure projection method
    real   , dimension(mcgrd)  , intent(in)    :: v0       ! v-velocity at previous time level
    real   , dimension(mcgrd)  , intent(inout) :: v1       ! v-velocity at current time level
    real   , dimension(mcgrd)  , intent(out)   :: vi       ! intermediate v-velocity in iterative process of pressure projection method
!
!   Local variables
!
    integer, save                              :: ient = 0 ! number of entries in this subroutine
    integer                                    :: j        ! iteration counter
    integer                                    :: m        ! loop counter
    integer                                    :: maxit    ! maximum number of iterations
    integer                                    :: md       ! index of point m-1
    integer                                    :: mdd      ! index of point m-2
    integer                                    :: mend     ! end index of loop over u-points
    integer                                    :: minf     ! m-index of point with largest error in solution
    integer                                    :: mlast    ! index of last internal u-point
    integer                                    :: msta     ! start index of loop over u-points
    integer                                    :: mu       ! index of point m+1
    integer                                    :: muu      ! index of point m+2
    integer                                    :: n        ! loop counter
    integer                                    :: nd       ! index of point n-1
    integer                                    :: ndd      ! index of point n-2
    integer                                    :: nddm     ! pointer to m,n-2
    integer                                    :: ndm      ! pointer to m,n-1
    integer                                    :: ndmd     ! pointer to m-1,n-1
    integer                                    :: ndmu     ! pointer to m+1,n-1
    integer                                    :: nend     ! end index of loop over v-points
    integer                                    :: nfm      ! pointer to m,nf
    integer                                    :: nfum     ! pointer to m,nfu
    integer                                    :: nfuum    ! pointer to m,nfuu
    integer                                    :: ninf     ! n-index of point with largest error in solution
    integer                                    :: nlast    ! index of last internal v-point
    integer                                    :: nlm      ! pointer to m,nl
    integer                                    :: nldm     ! pointer to m,nld
    integer                                    :: nlum     ! pointer to m,nlu
    integer                                    :: nm       ! pointer to m,n
    integer                                    :: nmd      ! pointer to m-1,n
    integer                                    :: nmdd     ! pointer to m-2,n
    integer                                    :: nmf      ! pointer to mf,n
    integer                                    :: nmfu     ! pointer to mfu,n
    integer                                    :: nmfuu    ! pointer to mfuu,n
    integer                                    :: nml      ! pointer to ml,n
    integer                                    :: nmld     ! pointer to mld,n
    integer                                    :: nmlu     ! pointer to mlu,n
    integer                                    :: nmu      ! pointer to m+1,n
    integer                                    :: nmuu     ! pointer to m+2,n
    integer                                    :: nsta     ! start index of loop over v-points
    integer                                    :: nu       ! index of point n+1
    integer                                    :: num      ! pointer to m,n+1
    integer                                    :: numd     ! pointer to m-1,n+1
    integer                                    :: nuu      ! index of point n+2
    integer                                    :: nuum     ! pointer to m,n+2
    !
    real                                       :: cfl      ! CFL number
    real                                       :: contrib  ! total contribution of momentum equation
    real                                       :: cporfr   ! flow resistance terms inside porous medium
    real                                       :: cveg     ! vegetation term
    real                                       :: denom    ! a denominator
    real                                       :: dxl      ! local mesh size in x-direction
    real                                       :: dyl      ! local mesh size in y-direction
    real                                       :: ener     ! total energy of closed system
    real                                       :: epslin   ! required accuracy in the linear solver
    real                                       :: fac      ! a factor
    real                                       :: fac1     ! another factor
    real                                       :: fac2     ! some other factor
    real                                       :: fac3     ! auxiliary factor
    real                                       :: fac4     ! auxiliary factor
    real                                       :: moutf    ! net mass outflow
    real                                       :: qf       ! updated mass flux at present face
    real                                       :: reps     ! accuracy of the final approximation
    real                                       :: res      ! residual
    real                                       :: resm     ! maximum error
    real                                       :: rval     ! auxiliary real
    real                                       :: s0mx     ! maximum value of water level
    real                                       :: sumqf    ! sum of outgoing mass fluxes over present cell
    real                                       :: theta3   ! implicitness factor for non-hydrostatic pressure gradient
    real                                       :: u        ! u-velocity at point different from its point of definition
    real                                       :: utot     ! velocity magnitude
    real                                       :: v        ! v-velocity at point different from its point of definition
    real                                       :: vol      ! total displaced volume of water
    real                                       :: zgrad    ! water level gradient term
    !
    logical                                    :: EQREAL   ! compares two reals
    logical                                    :: lmask    ! masks points for including curvature stress terms
    logical                                    :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashImpDepM2DHflow')
    !
    tetau  = 0.5
    tetav  = 0.5
    teta2  = 0.5
    theta3 = pnums(5)
    !
    ! adapt theta values underneath the floating object, if appropriate
    !
    if ( ifloat /= 0 ) then
       !
       do n = nfu, nl
          do m = mf, ml
             nm = kgrpnt(m,n)
             if ( presu(nm) == 1 ) tetau(nm) = pship(2)
          enddo
       enddo
       !
       call SWEXCHG ( tetau, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       do m = mfu, ml
          do n = nf, nl
             nm = kgrpnt(m,n)
             if ( presv(nm) == 1 ) tetav(nm) = pship(2)
          enddo
       enddo
       !
       call SWEXCHG ( tetav, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       nend = nl - 1
       if ( LMYL ) nend = nl
       !
       do m = mfu, mend
          do n = nfu, nend
             nm = kgrpnt(m,n)
             if ( presp(nm) == 1 ) teta2(nm) = pship(2)
          enddo
       enddo
       !
       call SWEXCHG ( teta2, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize theta values at appropriate boundaries in case of repeating grid
       !
       call periodic ( tetau, kgrpnt, 1, 1 )
       call periodic ( tetav, kgrpnt, 1, 1 )
       call periodic ( teta2, kgrpnt, 1, 1 )
       !
    endif
    !
    ! compute the mass flux per width in x-direction
    !
    do n = nfu, nl
       do m = mf, ml
          !
          nm = kgrpnt(m,n)
          !
          qx(nm) = hu(nm)*u0(nm)
          !
       enddo
    enddo
    !
    ! synchronize mass fluxes at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( qx, kgrpnt, 1, 1 )
    !
    ! compute the mass flux per width in x-direction in wl-point based on averaging
    !
    do n = nfu, nl
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m ,n)
          nmd = kgrpnt(md,n)
          !
          qxm(nm) = 0.5 * ( qx(nm) + qx(nmd) )
          !
       enddo
    enddo
    !
    ! synchronize mass fluxes at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( qxm, kgrpnt, 1, 1 )
    !
    ! compute the mass flux per width in y-direction
    !
    do n = nf, nl
       do m = mfu, ml
          !
          nm = kgrpnt(m,n)
          !
          qy(nm) = hv(nm)*v0(nm)
          !
       enddo
    enddo
    !
    ! synchronize mass fluxes at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( qy, kgrpnt, 1, 1 )
    !
    ! compute the mass flux per width in y-direction in depth point based on averaging
    !
    msta = mf + 1               ! first internal u-point
    if ( lreptx ) then
       mend = ml                ! last  internal u-point in case of repeating grid
    else
       mend = ml - 1            ! last  internal u-point
    endif
    !
    do n = nf, nl
       do m = msta, mend
          !
          mu = m + 1
          !
          nm  = kgrpnt(m ,n)
          nmu = kgrpnt(mu,n)
          !
          qym(nm) = 0.5 * ( qy(nm) + qy(nmu) )
          !
       enddo
    enddo
    !
    ! build gradient matrix for non-hydrostatic pressure in u-point
    !
    if ( ihydro == 1 ) then
       !
       do n = nfu, nl
          !
          nmf = kgrpnt(mf,n)
          nml = kgrpnt(ml,n)
          !
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( wetu(nm) * wetu(nmu) == 1 .and. (m /= mf .or. .not.LMXF) .and. (m /= ml .or. .not.LMXL) ) then
                !
                fac = 1. / ( 2. * gvu(nm) * hum(nm) )
                !
                if ( presu(nm) == 0 ) then
                   !
                   ! free surface flow
                   !
                   gmatu(nm,1) = -fac * ( s0(nm ) + dps(nmu) )
                   gmatu(nm,2) =  fac * ( s0(nmu) + dps(nm ) )
                   !
                else
                   !
                   ! pressurized flow
                   !
                   if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
                      !
                      ! changeover from free surface to pressurized
                      !
                      gmatu(nm,1) = -0.5 * fac * ( hs(nm) + hs(nmu) )
                      gmatu(nm,2) =        fac * ( hs(nm) + hs(nmu) )
                      !
                   else if ( presp(nm) == 1 .and. presp(nmu) == 0 ) then
                      !
                      ! changeover from pressurized to free surface
                      !
                      gmatu(nm,1) =       -fac * ( hs(nm) + hs(nmu) )
                      gmatu(nm,2) =  0.5 * fac * ( hs(nm) + hs(nmu) )
                      !
                   else
                      !
                      gmatu(nm,1) = -fac * ( hs(nm) + hs(nmu) )
                      gmatu(nm,2) =  fac * ( hs(nm) + hs(nmu) )
                      !
                   endif
                   !
                endif
                !
             else
                !
                gmatu(nm,1) = 0.
                gmatu(nm,2) = 0.
                !
             endif
             !
          enddo
          !
          if ( LMXF ) then
             !
             gmatu(nmf,2) = gmatu(nmf,2) - gmatu(nmf,1)
             gmatu(nmf,1) = 0.
             !
          endif
          !
          if ( LMXL ) then
             !
             gmatu(nml,1) = gmatu(nml,1) - gmatu(nml,2)
             gmatu(nml,2) = 0.
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute gradient of non-hydrostatic pressure in u-momentum equation
    !
    if ( ihydro == 1 ) then
       !
       do n = nfu, nl
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( wetu(nm) == 1 ) then
                !
                qgrad(nm) = gmatu(nm,1) * q(nm) + gmatu(nm,2) * q(nmu)
                !
             endif
             !
          enddo
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
       do n = nfu, nl
          do m = mf, ml
             !
             nm = kgrpnt(m,n)
             !
             if ( wetu(nm) == 1 ) then
                !
                wndimp(nm) = cwndu(nm) / max(1.e-3,hum(nm))
                !
             endif
             !
          enddo
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
       do n = nfu, nl
          do m = mf, ml
             !
             nm = kgrpnt(m,n)
             !
             if ( wetu(nm) == 1 ) then
                !
                cbot(nm) = cfricu(nm) / hum(nm)
                !
             endif
             !
          enddo
       enddo
       !
    else if ( irough /= 0 ) then
       !
       ! quadratic bottom friction
       !
       do n = nfu, nl
          do m = mf, ml
             !
             mu = m + 1
             nd = n - 1
             !
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             ndmu = kgrpnt(mu,nd)
             !
             if ( wetu(nm) == 1 ) then
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmu  == 1 ) nmu  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmu == 1 ) ndmu = ndm
                !
                ! compute v-velocity component and velocity magnitude
                !
                v = 0.25 * ( v0(nm) + v0(nmu) + v0(ndm) + v0(ndmu) )
                !
                utot = sqrt( u0(nm)*u0(nm) + v*v )
                !
                cbot(nm) = cfricu(nm) * utot / hum(nm)
                !
             endif
             !
          enddo
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
       do n = nfu, nl
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( wetu(nm) == 1 ) then
                !
                pgrad(nm) = 0.5 * grav * hum(nm) * (rho(nmu) - rho(nm)) / ( rhow*gvu(nm) )
                !
             endif
             !
          enddo
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
       do n = nfu, nl
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( wetu(nm) == 1 ) then
                !
                pgrad(nm) = pgrad(nm) + (patm(nmu) - patm(nm)) / ( rhow*gvu(nm) )
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute intermediate u-velocity
    !
    do n = nfu, nl
       !
       nmf  = kgrpnt(mf ,n)
       nmfu = kgrpnt(mfu,n)
       nml  = kgrpnt(ml ,n)
       !
       do m = msta, mend
          !
          md = m - 1
          mu = m + 1
          nd = n - 1
          nu = n + 1
          !
          nm   = kgrpnt(m ,n )
          nmd  = kgrpnt(md,n )
          nmu  = kgrpnt(mu,n )
          ndm  = kgrpnt(m ,nd)
          num  = kgrpnt(m ,nu)
          ndmu = kgrpnt(mu,nd)
          !
          if ( humn(nm) > 0. ) then
             !
             fac1 = 0.25 * dt * qxm(nmu) / ( gvu(nm) * humn(nm) )
             fac2 = 0.25 * dt * qxm(nm ) / ( gvu(nm) * humn(nm) )
             fac3 = 0.25 * dt * qym(nm ) / ( guu(nm) * humn(nm) )
             fac4 = 0.25 * dt * qym(ndm) / ( guu(nm) * humn(nm) )
             !
          else
             !
             fac1 = 0.
             fac2 = 0.
             fac3 = 0.
             fac4 = 0.
             !
          endif
          !
          if ( wetu(nm) == 1 .and. nmu /= 1 ) then
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmd  == 1 ) nmd  = nm
             if ( nmu  == 1 ) nmu  = nm
             if ( ndm  == 1 ) ndm  = nm
             if ( num  == 1 ) num  = nm
             if ( ndmu == 1 ) ndmu = ndm
             !
             ! compute v-velocity component and velocity magnitude
             !
             v = 0.25 * ( v0(nm) + v0(nmu) + v0(ndm) + v0(ndmu) )
             !
             utot = sqrt( u0(nm)*u0(nm) + v*v )
             !
             ! compute water level gradient
             !
             zgrad = grav * (s0(nmu) - s0(nm)) / gvu(nm)
             !
             ! compute flow resistance inside porous medium
             !
             cporfr = apomu(nm,1) + bpomu(nm,1) * utot
             !
             ! compute friction induced by vegetation
             !
             cveg = cvegu(nm,1,1) * utot
             !
             ! compute total contributions of the momentum equation (except advection)
             !
             contrib = zgrad + pgrad(nm) + qgrad(nm) - windu(nm)/max(1.e-3,hum(nm))
             !
             denom = 1. + cvegu(nm,1,2) + cpomu(nm,1) + dt * ( cbot(nm) + cporfr + cveg + wndimp(nm) )
             !
             ! compute advection term using central differences (momentum conservative and skew-symmetric) and build system of equations
             !
             amat(nm,2) = -fac2
             amat(nm,3) =  fac1
             amat(nm,4) = -fac4
             amat(nm,5) =  fac3
             amat(nm,1) =  denom - amat(nm,2) - amat(nm,3) - amat(nm,4) - amat(nm,5)
             rhs (nm  ) = (1.+cvegu(nm,1,2)+cpomu(nm,1)+amat(nm,2)+amat(nm,3)+amat(nm,4)+amat(nm,5))*u0(nm) - amat(nm,2)*u0(nmd) - amat(nm,3)*u0(nmu) - amat(nm,4)*u0(ndm) - amat(nm,5)*u0(num) - dt*contrib
             !
          else
             !
             amat(nm,1) = 1.
             amat(nm,2) = 0.
             amat(nm,3) = 0.
             amat(nm,4) = 0.
             amat(nm,5) = 0.
             rhs (nm  ) = 0.
             !
          endif
          !
       enddo
       !
       ! incorporate boundary conditions
       !
       if ( LMXF ) then
          !
          rhs (nmfu  ) = rhs(nmfu) - amat(nmfu,2)*u1(nmf)
          amat(nmfu,2) = 0.
          !
       endif
       !
       if ( .not.lreptx .and. LMXL ) then
          !
          amat(nml,1) = 1.
          amat(nml,2) = 0.
          amat(nml,3) = 0.
          amat(nml,4) = 0.
          amat(nml,5) = 0.
          rhs (nml  ) = u1(nml)
          !
       endif
       !
    enddo
    !
    if ( .not.lrepty ) then
       !
       do m = mfu, ml
          !
          nfm  = kgrpnt(m,nf )
          nfum = kgrpnt(m,nfu)
          nlm  = kgrpnt(m,nl )
          nlum = kgrpnt(m,nlu)
          !
          if ( ibb(m) == 1 .and. LMYF ) then
             !
             ! velocity is zero at closed boundary
             !
             amat(nfum,1) = amat(nfum,1) - amat(nfum,4)
             amat(nfum,4) = 0.
             !
          else if ( ibb(m) > 1 .and. LMYF ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             amat(nfum,1) = amat(nfum,1) + amat(nfum,4)
             amat(nfum,4) = 0.
             !
          endif
          !
          if ( ibt(m) == 1 .and. LMYL ) then
             !
             ! velocity is zero at closed boundary
             !
             amat(nlm,1) = amat(nlm,1) - amat(nlm,5)
             amat(nlm,5) = 0.
             !
          else if ( ibt(m) > 1 .and. LMYL ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             amat(nlm,1) = amat(nlm,1) + amat(nlm,5)
             amat(nlm,5) = 0.
             !
          endif
          !
       enddo
       !
    endif
    !
    ! solve system of equations
    !
    call sip( amat, rhs, u1 )
    if (STPNOW()) return
    !
    ! exchange u-velocities with neighbouring subdomains
    !
    call SWEXCHG ( u1, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize u-velocities at appropriate boundaries in case of repeating grid
    !
    call periodic ( u1, kgrpnt, 1, 1 )
    !
    ! compute the mass flux per width in x-direction in depth point based on averaging
    !
    nsta = nf + 1               ! first internal v-point
    if ( lrepty ) then
       nend = nl                ! last  internal v-point in case of repeating grid
    else
       nend = nl - 1            ! last  internal v-point
    endif
    !
    do m = mf, ml
       do n = nsta, nend
          !
          nu = n + 1
          !
          nm  = kgrpnt(m,n )
          num = kgrpnt(m,nu)
          !
          qxm(nm) = 0.5 * ( qx(nm) + qx(num) )
          !
       enddo
    enddo
    !
    ! compute the mass flux per width in y-direction in wl-point based on averaging
    !
    do m = mfu, ml
       do n = nfu, nl
          !
          nd = n - 1
          !
          nm  = kgrpnt(m,n )
          ndm = kgrpnt(m,nd)
          !
          qym(nm) = 0.5 * ( qy(nm) + qy(ndm) )
          !
       enddo
    enddo
    !
    ! synchronize mass fluxes at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( qym, kgrpnt, 1, 1 )
    !
    ! build gradient matrix for non-hydrostatic pressure in v-point
    !
    if ( ihydro == 1 ) then
       !
       do m = mfu, ml
          !
          nfm = kgrpnt(m,nf)
          nlm = kgrpnt(m,nl)
          !
          do n = nf, nl
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( wetv(nm) * wetv(num) == 1 .and. (n /= nf .or. .not.LMYF) .and. (n /= nl .or. .not.LMYL) ) then
                !
                fac = 1. / ( 2. * guv(nm) * hvm(nm) )
                !
                if ( presv(nm) == 0 ) then
                   !
                   ! free surface flow
                   !
                   gmatv(nm,1) = -fac * ( s0(nm ) + dps(num) )
                   gmatv(nm,2) =  fac * ( s0(num) + dps(nm ) )
                   !
                else
                   !
                   ! pressurized flow
                   !
                   if ( presp(nm) == 0 .and. presp(num) == 1 ) then
                      !
                      ! changeover from free surface to pressurized
                      !
                      gmatv(nm,1) = -0.5 * fac * ( hs(nm) + hs(num) )
                      gmatv(nm,2) =        fac * ( hs(nm) + hs(num) )
                      !
                   else if ( presp(nm) == 1 .and. presp(num) == 0 ) then
                      !
                      ! changeover from pressurized to free surface
                      !
                      gmatv(nm,1) =       -fac * ( hs(nm) + hs(num) )
                      gmatv(nm,2) =  0.5 * fac * ( hs(nm) + hs(num) )
                      !
                   else
                      !
                      gmatv(nm,1) = -fac * ( hs(nm) + hs(num) )
                      gmatv(nm,2) =  fac * ( hs(nm) + hs(num) )
                      !
                   endif
                   !
                endif
                !
             else
                !
                gmatv(nm,1) = 0.
                gmatv(nm,2) = 0.
                !
             endif
             !
          enddo
          !
          if ( LMYF ) then
             !
             gmatv(nfm,2) = gmatv(nfm,2) - gmatv(nfm,1)
             gmatv(nfm,1) = 0.
             !
          endif
          !
          if ( LMYL ) then
             !
             gmatv(nlm,1) = gmatv(nlm,1) - gmatv(nlm,2)
             gmatv(nlm,2) = 0.
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute gradient of non-hydrostatic pressure in v-momentum equation
    !
    if ( ihydro == 1 ) then
       !
       do m = mfu, ml
          do n = nf, nl
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( wetv(nm) == 1 ) then
                !
                qgrad(nm) = gmatv(nm,1) * q(nm) + gmatv(nm,2) * q(num)
                !
             endif
             !
          enddo
       enddo
       !
       if ( iproj == 2 ) qgrad = (1.-theta3) * qgrad
       !
    endif
    !
    ! compute implicit part of wind stress term, if appropriate
    !
    if ( relwnd ) then
       !
       do m = mfu, ml
          do n = nf, nl
             !
             nm = kgrpnt(m,n)
             !
             if ( wetv(nm) == 1 ) then
                !
                wndimp(nm) = cwndv(nm) / max(1.e-3,hvm(nm))
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute bottom friction term, if appropriate
    !
    if ( irough == 11 ) then
       !
       ! linear bottom friction
       !
       do m = mfu, ml
          do n = nf, nl
             !
             nm = kgrpnt(m,n)
             !
             if ( wetv(nm) == 1 ) then
                !
                cbot(nm) = cfricv(nm) / hvm(nm)
                !
             endif
             !
          enddo
       enddo
       !
    else if ( irough /= 0 ) then
       !
       ! quadratic bottom friction
       !
       do m = mfu, ml
          do n = nf, nl
             !
             nu = n + 1
             md = m - 1
             !
             nm   = kgrpnt(m ,n )
             num  = kgrpnt(m ,nu)
             nmd  = kgrpnt(md,n )
             numd = kgrpnt(md,nu)
             !
             if ( wetv(nm) == 1 ) then
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( num  == 1 ) num  = nm
                if ( nmd  == 1 ) nmd  = nm
                if ( numd == 1 ) numd = nmd
                !
                ! compute u-velocity component and velocity magnitude
                !
                u = 0.25 * ( u0(nm) + u0(num) + u0(nmd) + u0(numd) )
                !
                utot = sqrt( u*u + v0(nm)*v0(nm) )
                !
                cbot(nm) = cfricv(nm) * utot / hvm(nm)
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute baroclinic forcing
    !
    if ( idens /= 0 ) then
       !
       do m = mfu, ml
          do n = nf, nl
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( wetv(nm) == 1 ) then
                !
                pgrad(nm) = 0.5 * grav * hvm(nm) * (rho(num) - rho(nm)) / ( rhow*guv(nm) )
                !
             endif
             !
          enddo
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
       do m = mfu, ml
          do n = nf, nl
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( wetv(nm) == 1 ) then
                !
                pgrad(nm) = pgrad(nm) + (patm(num) - patm(nm)) / ( rhow*guv(nm) )
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute intermediate v-velocity
    !
    do m = mfu, ml
       !
       nfm  = kgrpnt(m,nf )
       nfum = kgrpnt(m,nfu)
       nlm  = kgrpnt(m,nl )
       !
       do n = nsta, nend
          !
          nd = n - 1
          nu = n + 1
          md = m - 1
          mu = m + 1
          !
          nm   = kgrpnt(m ,n )
          ndm  = kgrpnt(m ,nd)
          num  = kgrpnt(m ,nu)
          nmd  = kgrpnt(md,n )
          nmu  = kgrpnt(mu,n )
          numd = kgrpnt(md,nu)
          !
          if ( hvmn(nm) > 0. ) then
             !
             fac1 = 0.25 * dt * qym(num) / ( guv(nm) * hvmn(nm) )
             fac2 = 0.25 * dt * qym(nm ) / ( guv(nm) * hvmn(nm) )
             fac3 = 0.25 * dt * qxm(nm ) / ( gvv(nm) * hvmn(nm) )
             fac4 = 0.25 * dt * qxm(nmd) / ( gvv(nm) * hvmn(nm) )
             !
          else
             !
             fac1 = 0.
             fac2 = 0.
             fac3 = 0.
             fac4 = 0.
             !
          endif
          !
          if ( wetv(nm) == 1 .and. num /= 1 ) then
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( ndm  == 1 ) ndm  = nm
             if ( num  == 1 ) num  = nm
             if ( nmd  == 1 ) nmd  = nm
             if ( nmu  == 1 ) nmu  = nm
             if ( numd == 1 ) numd = nmd
             !
             ! compute u-velocity component and velocity magnitude
             !
             u = 0.25 * ( u0(nm) + u0(num) + u0(nmd) + u0(numd) )
             !
             utot = sqrt( u*u + v0(nm)*v0(nm) )
             !
             ! compute water level gradient
             !
             zgrad = grav * (s0(num) - s0(nm)) / guv(nm)
             !
             ! compute flow resistance inside porous medium
             !
             cporfr = apomv(nm,1) + bpomv(nm,1) * utot
             !
             ! compute friction induced by vegetation
             !
             cveg = cvegv(nm,1,1) * utot
             !
             ! compute total contributions of the momentum equation
             !
             contrib = zgrad + pgrad(nm) + qgrad(nm) - windv(nm)/max(1.e-3,hvm(nm))
             !
             denom = 1. + cvegv(nm,1,2) + cpomv(nm,1) + dt * ( cbot(nm) + cporfr + cveg + wndimp(nm) )
             !
             ! compute advection term using central differences (momentum conservative and skew-symmetric) and build system of equations
             !
             amat(nm,2) = -fac4
             amat(nm,3) =  fac3
             amat(nm,4) = -fac2
             amat(nm,5) =  fac1
             amat(nm,1) =  denom - amat(nm,2) - amat(nm,3) - amat(nm,4) - amat(nm,5)
             rhs (nm  ) = (1.+cvegv(nm,1,2)+cpomv(nm,1)+amat(nm,2)+amat(nm,3)+amat(nm,4)+amat(nm,5))*v0(nm) - amat(nm,2)*v0(nmd) - amat(nm,3)*v0(nmu) - amat(nm,4)*v0(ndm) - amat(nm,5)*v0(num) - dt*contrib
             !
          else
             !
             amat(nm,1) = 1.
             amat(nm,2) = 0.
             amat(nm,3) = 0.
             amat(nm,4) = 0.
             amat(nm,5) = 0.
             rhs (nm  ) = 0.
             !
          endif
          !
       enddo
       !
       ! incorporate boundary conditions
       !
       if ( LMYF ) then
          !
          rhs (nfum  ) = rhs(nfum) - amat(nfum,4)*v1(nfm)
          amat(nfum,4) = 0.
          !
       endif
       !
       if ( .not.lrepty .and. LMYL ) then
          !
          amat(nlm,1) = 1.
          amat(nlm,2) = 0.
          amat(nlm,3) = 0.
          amat(nlm,4) = 0.
          amat(nlm,5) = 0.
          rhs (nlm  ) = v1(nlm)
          !
       endif
       !
    enddo
    !
    if ( .not.lreptx ) then
       !
       do n = nfu, nl
          !
          nmf  = kgrpnt(mf, n)
          nmfu = kgrpnt(mfu,n)
          nml  = kgrpnt(ml ,n)
          nmlu = kgrpnt(mlu,n)
          !
          if ( ibl(n) == 1 .and. LMXF ) then
             !
             ! velocity is zero at closed boundary
             !
             amat(nmfu,1) = amat(nmfu,1) - amat(nmfu,2)
             amat(nmfu,2) = 0.
             !
          else if ( ibl(n) > 1 .and. LMXF ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             amat(nmfu,1) = amat(nmfu,1) + amat(nmfu,2)
             amat(nmfu,2) = 0.
             !
          endif
          !
          if ( ibr(n) == 1 .and. LMXL ) then
             !
             ! velocity is zero at closed boundary
             !
             amat(nml,1) = amat(nml,1) - amat(nml,3)
             amat(nml,3) = 0.
             !
          else if ( ibr(n) > 1 .and. LMXL ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             amat(nml,1) = amat(nml,1) + amat(nml,3)
             amat(nml,3) = 0.
             !
          endif
          !
       enddo
       !
    endif
    !
    ! solve system of equations
    !
    call sip( amat, rhs, v1 )
    if (STPNOW()) return
    !
    ! exchange v-velocities with neighbouring subdomains
    !
    call SWEXCHG ( v1, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize v-velocities at appropriate boundaries in case of repeating grid
    !
    call periodic ( v1, kgrpnt, 1, 1 )
    !
    ui = u1
    vi = v1
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
       ninf = 0
       !
       ! compute the water level
       !
       fac = grav * dt * dt
       !
       ! first, build the equation for water level correction
       !
       if ( inewt == 0 ) then
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
                if ( ( m /= mfu .or. .not.LMXF ) .and. nmd /= 1 ) then
                   !
                   amat(nm,2) = -fac * tetau(nmd) * guu(nmd) * hu(nmd) / gvu(nmd)
                   !
                else
                   !
                   amat(nm,2) = 0.
                   !
                endif
                !
                if ( ( m /= ml .or. .not.LMXL ) .and. nmu /= 1 ) then
                   !
                   amat(nm,3) = -fac * tetau(nm) * guu(nm) * hu(nm) / gvu(nm)
                   !
                else
                   !
                   amat(nm,3) = 0.
                   !
                endif
                !
                if ( ( n /= nfu .or. .not.LMYF ) .and. ndm /= 1 ) then
                   !
                   amat(nm,4) = -fac * tetav(ndm) * gvv(ndm) * hv(ndm) / guv(ndm)
                   !
                else
                   !
                   amat(nm,4) = 0.
                   !
                endif
                !
                if ( ( n /= nl .or. .not.LMYL ) .and. num /= 1 ) then
                   !
                   amat(nm,5) = -fac * tetav(nm) * gvv(nm) * hv(nm) / guv(nm)
                   !
                else
                   !
                   amat(nm,5) = 0.
                   !
                endif
                !
                amat(nm,1) = fac1 * gsqs(nm) - teta2(nm) * ( amat(nm,2) + amat(nm,3) + amat(nm,4) + amat(nm,5) )
                amat(nm,2) = teta2(nmd) * amat(nm,2)
                amat(nm,3) = teta2(nmu) * amat(nm,3)
                amat(nm,4) = teta2(ndm) * amat(nm,4)
                amat(nm,5) = teta2(num) * amat(nm,5)
                rhs (nm  ) = dt * ( guu(nmd)*hu(nmd) * ( tetau(nmd)*ui(nmd) + (1.-tetau(nmd))*u0(nmd) ) - guu(nm)*hu(nm) * ( tetau(nm)*ui(nm) + (1.-tetau(nm))*u0(nm) ) + &
                                    gvv(ndm)*hv(ndm) * ( tetav(ndm)*vi(ndm) + (1.-tetav(ndm))*v0(ndm) ) - gvv(nm)*hv(nm) * ( tetav(nm)*vi(nm) + (1.-tetav(nm))*v0(nm) ) )
                !
                ! add contribution from rigid body motions
                !
                rhs(nm) = rhs(nm) - dt * gsqs(nm) * skc(nm)
                !
             enddo
          enddo
          !
          ! add mass source due to internal wave generation
          !
          if ( iwvgen /= 0 ) then
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   rhs(nm) = rhs(nm) + dt * gsqs(nm) * srcm(nm)
                   !
                enddo
             enddo
             !
          endif
          !
          ! add non-hydrostatic pressure part in the equation for water level correction (not piezometric head), if appropriate
          !
          if ( lpproj ) then
             !
             ! not computed for point ml or nl at periodic boundary
             !
             if ( lreptx ) then
                mlast = ml - 1
             else
                mlast = ml
             endif
             !
             if ( lrepty ) then
                nlast = nl - 1
             else
                nlast = nl
             endif
             !
             do n = nfu, nlast
                do m = mfu, mlast
                   !
                   md = m - 1
                   nd = n - 1
                   !
                   nm  = kgrpnt(m ,n )
                   nmd = kgrpnt(md,n )
                   ndm = kgrpnt(m ,nd)
                   !
                   if ( presp(nm) == 0 ) then
                      rhs(nm) = rhs(nm) + dt * dt * theta3 * ( tetau(nm) * guu(nm) * hu(nm) * dqgrdu(nm) - tetau(nmd) * guu(nmd) * hu(nmd) * dqgrdu(nmd) + &
                                                               tetav(nm) * gvv(nm) * hv(nm) * dqgrdv(nm) - tetav(ndm) * gvv(ndm) * hv(ndm) * dqgrdv(ndm) )
                   endif
                   !
                enddo
             enddo
             !
          endif
          !
          ! next, solve the equation for water level correction
          !
          call pcg( amat, rhs, ds )
          if (STPNOW()) return
          !
       else
          !
          lon = -dps  - s0
          upn = -flos - s0
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
                if ( ( m /= mfu .or. .not.LMXF ) .and. nmd /= 1 ) then
                   !
                   amatn(nm,2) = -fac * tetau(nmd) * guu(nmd) * hu(nmd) / gvu(nmd)
                   !
                else
                   !
                   amatn(nm,2) = 0.
                   !
                endif
                !
                if ( ( m /= ml .or. .not.LMXL ) .and. nmu /= 1 ) then
                   !
                   amatn(nm,3) = -fac * tetau(nm) * guu(nm) * hu(nm) / gvu(nm)
                   !
                else
                   !
                   amatn(nm,3) = 0.
                   !
                endif
                !
                if ( ( n /= nfu .or. .not.LMYF ) .and. ndm /= 1 ) then
                   !
                   amatn(nm,4) = -fac * tetav(ndm) * gvv(ndm) * hv(ndm) / guv(ndm)
                   !
                else
                   !
                   amatn(nm,4) = 0.
                   !
                endif
                !
                if ( ( n /= nl .or. .not.LMYL ) .and. num /= 1 ) then
                   !
                   amatn(nm,5) = -fac * tetav(nm) * gvv(nm) * hv(nm) / guv(nm)
                   !
                else
                   !
                   amatn(nm,5) = 0.
                   !
                endif
                !
                amatn(nm,1) = fac1 * gsqs(nm) - teta2(nm) * ( amatn(nm,2) + amatn(nm,3) + amatn(nm,4) + amatn(nm,5) )
                amatn(nm,2) = teta2(nmd) * amatn(nm,2)
                amatn(nm,3) = teta2(nmu) * amatn(nm,3)
                amatn(nm,4) = teta2(ndm) * amatn(nm,4)
                amatn(nm,5) = teta2(num) * amatn(nm,5)
                rhsn (nm  ) = dt * ( guu(nmd)*hu(nmd) * ( tetau(nmd)*ui(nmd) + (1.-tetau(nmd))*u0(nmd) ) - guu(nm)*hu(nm) * ( tetau(nm)*ui(nm) + (1.-tetau(nm))*u0(nm) ) + &
                                     gvv(ndm)*hv(ndm) * ( tetav(ndm)*vi(ndm) + (1.-tetav(ndm))*v0(ndm) ) - gvv(nm)*hv(nm) * ( tetav(nm)*vi(nm) + (1.-tetav(nm))*v0(nm) ) )
                !
                ! add contribution from rigid body motions
                !
                rhsn(nm) = rhsn(nm) - dt * gsqs(nm) * skc(nm)
                !
                amatn(nm,:) = amatn(nm,:) / gsqs(nm)
                rhsn (nm  ) = rhsn (nm  ) / gsqs(nm)
                !
                rhsn(nm) = rhsn(nm) + max( lon(nm), min( upn(nm), 0.) )
                !
             enddo
          enddo
          !
          ! add mass source due to internal wave generation
          !
          if ( iwvgen /= 0 ) then
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   rhsn(nm) = rhsn(nm) + dt * srcm(nm)
                   !
                enddo
             enddo
             !
          endif
          !
          ! add non-hydrostatic pressure part in the equation for water level correction (not piezometric head), if appropriate
          !
          if ( lpproj ) then
             !
             ! not computed for point ml or nl at periodic boundary
             !
             if ( lreptx ) then
                mlast = ml - 1
             else
                mlast = ml
             endif
             !
             if ( lrepty ) then
                nlast = nl - 1
             else
                nlast = nl
             endif
             !
             do n = nfu, nlast
                do m = mfu, mlast
                   !
                   md = m - 1
                   nd = n - 1
                   !
                   nm  = kgrpnt(m ,n )
                   nmd = kgrpnt(md,n )
                   ndm = kgrpnt(m ,nd)
                   !
                   if ( presp(nm) == 0 ) then
                      rhsn(nm) = rhsn(nm) + dt * dt * theta3 * ( tetau(nm) * guu(nm) * hu(nm) * dqgrdu(nm) - tetau(nmd) * guu(nmd) * hu(nmd) * dqgrdu(nmd) + &
                                                                 tetav(nm) * gvv(nm) * hv(nm) * dqgrdv(nm) - tetav(ndm) * gvv(ndm) * hv(ndm) * dqgrdv(ndm) ) / gsqs(nm)
                   endif
                   !
                enddo
             enddo
             !
          endif
          !
          ! next, solve the equation for water level correction
          !
          call newton2D ( amatn, rhsn, lon, upn, ds )
          if (STPNOW()) return
          !
       endif
       !
       ! set to zero for permanently dry points
       !
       ds(1) = 0.
       !
       ! exchange water level corrections with neighbouring subdomains
       !
       call SWEXCHG ( ds, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize water level correction at appropriate boundaries in case of repeating grid
       !
       call periodic ( ds, kgrpnt, 1, 1 )
       !
       ! determine maximum error to check convergence projection pressure method
       !
       do m = mfu, ml
          do n = nfu, nl
             !
             nm = kgrpnt(m,n)
             !
             res = abs(s0(nm)+ds(nm) - s1(nm))
             if ( res > resm ) then
                resm = res
                minf = m
                ninf = n
             endif
             !
          enddo
       enddo
       !
       rval = resm
       call SWREDUCE ( resm, 1, SWREAL, SWMAX )
       !
       if ( iamout == 2 .and. EQREAL(resm,rval) .and. maxit > 1 ) then
          !
          write (PRINTF,'(a,i4,a,e12.6,a,i5,a,i5,a)') ' ++ pressure projection: iter = ',j,' res = ',resm,' in (m,n)=(',minf+MXF-2,',',ninf+MYF-2,')'
          !
       endif
       !
       ! correct horizontal velocities
       !
       do n = nfu, nl
          !
          do m = msta, mend
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( wetu(nm) == 1 .and. nmu /= 1 ) then
                !
                u1(nm) = ui(nm) - grav * dt * ( teta2(nmu) * ds(nmu) - teta2(nm) * ds(nm) ) / gvu(nm)
                !
             endif
             !
          enddo
       enddo
       !
       do m = mfu, ml
          !
          do n = nsta, nend
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( wetv(nm) == 1 .and. num /= 1 ) then
                !
                v1(nm) = vi(nm) - grav * dt * ( teta2(num) * ds(num) - teta2(nm) * ds(nm) ) / guv(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! exchange velocities with neighbouring subdomains
       !
       call SWEXCHG ( u1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       call SWEXCHG ( v1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize velocities at appropriate boundaries in case of repeating grid
       !
       call periodic ( u1, kgrpnt, 1, 1 )
       call periodic ( v1, kgrpnt, 1, 1 )
       !
       ! update water level
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             s1(nm) = s0(nm) + ds(nm)
             !
          enddo
       enddo
       !
       ! check positivity of the water depth in each cell of the domain, if appropriate
       !
       if ( inewt == 0 ) then
          !
          cflmax = -999.
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
                !
                ! compute the sum of mass fluxes leaving the cell
                !
                sumqf = 0.
                !
                qf = hu(nm) * guu(nm) * ( tetau(nm)*u1(nm) + (1.-tetau(nm))*u0(nm) )
                !
                if ( qf > 0. ) sumqf = sumqf + qf
                !
                qf = hu(nmd) * guu(nmd) * ( tetau(nmd)*u1(nmd) + (1.-tetau(nmd))*u0(nmd) )
                !
                if ( qf < 0. ) sumqf = sumqf - qf
                !
                qf = hv(nm) * gvv(nm) * ( tetav(nm)*v1(nm) + (1.-tetav(nm))*v0(nm) )
                !
                if ( qf > 0. ) sumqf = sumqf + qf
                !
                qf = hv(ndm) * gvv(ndm) * ( tetav(ndm)*v1(ndm) + (1.-tetav(ndm))*v0(ndm) )
                !
                if ( qf < 0. ) sumqf = sumqf - qf
                !
                ! compute the "flow" Courant number
                !
                if ( hs(nm) > epsdry ) then
                   !
                   cfl = sumqf * dt / hs(nm) / gsqs(nm)
                   if ( cfl > cflmax ) cflmax = cfl
                   !
                endif
                !
             enddo
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
       if ( .not.lreptx ) then
          !
          do n = nfu, nl
             !
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             nml  = kgrpnt(ml ,n)
             nmlu = kgrpnt(mlu,n)
             !
             if ( ibl(n) /= 6 .and. LMXF ) s1(nmf ) = s1(nmfu)
             if ( ibr(n) /= 6 .and. LMXL ) s1(nmlu) = s1(nml )
             !
          enddo
          !
       endif
       !
       if ( .not.lrepty ) then
          !
          do m = mfu, ml
             !
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             nlm  = kgrpnt(m,nl )
             nlum = kgrpnt(m,nlu)
             !
             if ( ibb(m) /= 6 .and. LMYF ) s1(nfm ) = s1(nfum)
             if ( ibt(m) /= 6 .and. LMYL ) s1(nlum) = s1(nlm )
             !
          enddo
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
          do n = nfu, nl
             do m = mfu, ml
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd == 1 ) nmd = nm
                if ( ndm == 1 ) ndm = nm
                !
                if ( wets(nm) == 1 ) then
                   !
                   dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
                   dyl = 0.5 * ( guu(nm) + guu(nmd) )
                   !
                   w1bot(nm) = -0.5 * ( ( u1(nm) + u1(nmd) ) * ( dpu(nm) - dpu(nmd) ) / dxl + ( v1(nm) + v1(ndm) ) * ( dpv(nm) - dpv(ndm) ) / dyl )
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
                      w1top(nm) = -0.5 * ( ( u1(nm) + u1(nmd) ) * ( flou(nm) - flou(nmd) ) / dxl + ( v1(nm) + v1(ndm) ) * ( flov(nm) - flov(ndm) ) / dyl )
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
          enddo
          !
          if ( horwinc ) then
             !
             ! compute horizontal advection terms (momentum conservative and explicit)
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   md  = m  - 1
                   mdd = md - 1
                   mu  = m  + 1
                   muu = mu + 1
                   !
                   if (                   LMXF .and. mdd < mf  ) mdd = mf
                   if ( .not.lreptx .and. LMXL .and. muu > ml  ) muu = ml
                   if (      lreptx .and. LMXL .and. muu > mlu ) muu = mlu
                   !
                   nd  = n  - 1
                   ndd = nd - 1
                   nu  = n  + 1
                   nuu = nu + 1
                   !
                   if (                   LMYF .and. ndd < nf  ) ndd = nf
                   if ( .not.lrepty .and. LMYL .and. nuu > nl  ) nuu = nl
                   if (      lrepty .and. LMYL .and. nuu > nlu ) nuu = nlu
                   !
                   nm   = kgrpnt(m  ,n  )
                   nmd  = kgrpnt(md ,n  )
                   nmdd = kgrpnt(mdd,n  )
                   nmu  = kgrpnt(mu ,n  )
                   nmuu = kgrpnt(muu,n  )
                   ndm  = kgrpnt(m  ,nd )
                   nddm = kgrpnt(m  ,ndd)
                   num  = kgrpnt(m  ,nu )
                   nuum = kgrpnt(m  ,nuu)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd  == 1 ) nmd  = nm
                   if ( nmdd == 1 ) nmdd = nmd
                   if ( nmu  == 1 ) nmu  = nm
                   if ( nmuu == 1 ) nmuu = nmu
                   if ( ndm  == 1 ) ndm  = nm
                   if ( nddm == 1 ) nddm = ndm
                   if ( num  == 1 ) num  = nm
                   if ( nuum == 1 ) nuum = num
                   !
                   if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
                      !
                      ! advection term in x-direction
                      !
                      fac = 2. * dt / gsqs(nm)
                      !
                      fac1 = fac * guu(nmd) * qx(nmd) / hs(nm)
                      fac2 = fac * guu(nm ) * qx(nm ) / hs(nm)
                      !
                      w1top(nm) = w1top(nm) - 0.5 * ( fac1 * (w0top(nm) - w0top(nmd)) + fac2 * (w0top(nmu) - w0top(nm)) )
                      !
                      ! advection term in y-direction
                      !
                      fac1 = fac * gvv(ndm) * qy(ndm) / hs(nm)
                      fac2 = fac * gvv(nm ) * qy(nm ) / hs(nm)
                      !
                      w1top(nm) = w1top(nm) - 0.5 * ( fac1 * (w0top(nm) - w0top(ndm)) + fac2 * (w0top(num) - w0top(nm)) )
                      !
                   endif
                   !
                enddo
             enddo
             !
          endif
          !
          ! add vegetation dissipation due to horizontal cylinders
          !
          if ( ivegw /= 0 ) then
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   md = m - 1
                   nd = n - 1
                   !
                   nm  = kgrpnt(m ,n )
                   nmd = kgrpnt(md,n )
                   ndm = kgrpnt(m ,nd)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd == 1 ) nmd = nm
                   if ( ndm == 1 ) ndm = nm
                   !
                   if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
                      !
                      ! compute vegetation factor
                      !
                      cveg = 0.125 * ( cvegu(nm,1,1) + cvegu(nmd,1,1) + cvegv(nm,1,1) + cvegv(ndm,1,1) ) * dt * abs( w0top(nm) + w0bot(nm) )
                      !
                      w1top(nm) = w1top(nm) - cveg * w1bot(nm)
                      w1top(nm) = w1top(nm) / ( 1. + cveg )
                      !
                   endif
                   !
                enddo
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
          do n = nfu, nl
             do m = mfu, ml
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd == 1 ) nmd = nm
                if ( ndm == 1 ) ndm = nm
                !
                if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
                   !
                   dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
                   dyl = 0.5 * ( guu(nm) + guu(nmd) )
                   !
                   ! u-contribution
                   !
                   fac1 =  hu(nm ) * guu(nm ) + gsqs(nm) * ( dpu(nm) - dpu(nmd) - 0.5 * ( hu(nm) - hu(nmd) ) ) / dxl
                   fac2 = -hu(nmd) * guu(nmd) + gsqs(nm) * ( dpu(nm) - dpu(nmd) - 0.5 * ( hu(nm) - hu(nmd) ) ) / dxl
                   !
                   amat(nm,2) = fac2 * gmatu(nmd,1)
                   amat(nm,3) = fac1 * gmatu(nm ,2)
                   amat(nm,1) = fac1 * gmatu(nm ,1) + fac2 * gmatu(nmd,2)
                   rhs (nm  ) = fac1*u1(nm) + fac2*u1(nmd)
                   !
                   ! v-contribution
                   !
                   fac1 =  hv(nm ) * gvv(nm ) + gsqs(nm) * ( dpv(nm) - dpv(ndm) - 0.5 * ( hv(nm) - hv(ndm) ) ) / dyl
                   fac2 = -hv(ndm) * gvv(ndm) + gsqs(nm) * ( dpv(nm) - dpv(ndm) - 0.5 * ( hv(nm) - hv(ndm) ) ) / dyl
                   !
                   amat(nm,4) = fac2 * gmatv(ndm,1)
                   amat(nm,5) = fac1 * gmatv(nm ,2)
                   amat(nm,1) = amat(nm,1) + fac1 * gmatv(nm,1) + fac2 * gmatv(ndm,2)
                   rhs (nm  ) = rhs (nm  ) + fac1*v1(nm) + fac2*v1(ndm)
                   !
                   ! w-contribution
                   !
                   amat(nm,1) = amat(nm,1) - 2.*gsqs(nm) / hs(nm)
                   rhs (nm  ) = rhs (nm  ) + gsqs(nm) * ( w1top(nm) + w1bot(nm) )
                   !
                   rhs(nm) = rhs(nm)/(dt*theta3)
                   !
                else if ( wets(nm) == 1 .and. presp(nm) == 1 ) then
                   !
                   ! u-contribution
                   !
                   fac1 =  hu(nm ) * guu(nm )
                   fac2 = -hu(nmd) * guu(nmd)
                   !
                   amat(nm,2) = fac2 * gmatu(nmd,1)
                   amat(nm,3) = fac1 * gmatu(nm ,2)
                   amat(nm,1) = fac1 * gmatu(nm ,1) + fac2 * gmatu(nmd,2)
                   rhs (nm  ) = fac1*u1(nm) + fac2*u1(nmd)
                   !
                   ! v-contribution
                   !
                   fac1 =  hv(nm ) * gvv(nm )
                   fac2 = -hv(ndm) * gvv(ndm)
                   !
                   amat(nm,4) = fac2 * gmatv(ndm,1)
                   amat(nm,5) = fac1 * gmatv(nm ,2)
                   amat(nm,1) = amat(nm,1) + fac1 * gmatv(nm,1) + fac2 * gmatv(ndm,2)
                   rhs (nm  ) = rhs (nm  ) + fac1*v1(nm) + fac2*v1(ndm)
                   !
                   ! add contribution from rigid body motions
                   !
                   rhs(nm) = rhs(nm) + gsqs(nm) * skc(nm)
                   !
                   rhs(nm) = rhs(nm)/(dt*theta3)
                   !
                else
                   !
                   amat(nm,1) = 1.
                   amat(nm,2) = 0.
                   amat(nm,3) = 0.
                   amat(nm,4) = 0.
                   amat(nm,5) = 0.
                   rhs (nm  ) = 0.
                   !
                endif
                !
             enddo
          enddo
          !
          ! solve the Poisson equation
          !
          call sip( amat, rhs, dq )
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
          do n = nfu, nl
             do m = mf, ml
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                if ( wetu(nm) == 1 .and. presu(nm) == 0 ) then
                   !
                   dqgrdu(nm) = gmatu(nm,1) * dq(nm) + gmatu(nm,2) * dq(nmu)
                   !
                else
                   ! no correction applied on piezometric head
                   !
                   dqgrdu(nm) = 0.
                   !
                endif
                !
             enddo
          enddo
          !
          do m = mfu, ml
             do n = nf, nl
                !
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                if ( wetv(nm) == 1 .and. presv(nm) == 0 ) then
                   !
                   dqgrdv(nm) = gmatv(nm,1) * dq(nm) + gmatv(nm,2) * dq(num)
                   !
                else
                   ! no correction applied on piezometric head
                   !
                   dqgrdv(nm) = 0.
                   !
                endif
                !
             enddo
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
       do n = nfu, nl
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( wetu(nm) == 1 .and. nmu /= 1 ) then
                !
                u1(nm) = u1(nm) - dt*theta3*( gmatu(nm,1) * dq(nm) + gmatu(nm,2) * dq(nmu) )
                !
             endif
             !
          enddo
       enddo
       !
       ! v-velocity
       !
       do m = mfu, ml
          do n = nf, nl
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( wetv(nm) == 1 .and. num /= 1 ) then
                !
                v1(nm) = v1(nm) - dt*theta3*( gmatv(nm,1) * dq(nm) + gmatv(nm,2) * dq(num) )
                !
             endif
             !
          enddo
       enddo
       !
       ! w-velocity
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             if ( nmd == 1 ) nmd = nm
             if ( ndm == 1 ) ndm = nm
             !
             if ( wets(nm) == 1 ) then
                !
                dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
                dyl = 0.5 * ( guu(nm) + guu(nmd) )
                !
                if ( presp(nm) == 0 ) w1top(nm) = w1top(nm) + w1bot(nm) + 2.*dt*theta3*dq(nm)/hs(nm)
                w1bot(nm) = -0.5 * ( ( u1(nm) + u1(nmd) ) * ( dpu(nm) - dpu(nmd) ) / dxl + ( v1(nm) + v1(ndm) ) * ( dpv(nm) - dpv(ndm) ) / dyl )
                if ( presp(nm) == 0 ) then
                   w1top(nm) = w1top(nm) - w1bot(nm)
                else
                   w1top(nm) = -0.5 * ( ( u1(nm) + u1(nmd) ) * ( flou(nm) - flou(nmd) ) / dxl + ( v1(nm) + v1(ndm) ) * ( flov(nm) - flov(ndm) ) / dyl )
                   ! add body motion
                   w1top(nm) = w1top(nm) + skc(nm)
                endif
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! set tangential velocity in virtual cells
    !
    if ( .not.lreptx ) then
       !
       do n = nf, nl
          !
          nmf  = kgrpnt(mf, n)
          nmfu = kgrpnt(mfu,n)
          nml  = kgrpnt(ml ,n)
          nmlu = kgrpnt(mlu,n)
          !
          if ( ibl(n) == 1 .and. LMXF ) then
             !
             ! velocity is zero at closed boundary
             !
             v1(nmf) = -v1(nmfu)
             !
          else if ( ibl(n) > 1 .and. LMXF ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             v1(nmf) = v1(nmfu)
             !
          endif
          !
          if ( ibr(n) == 1 .and. LMXL ) then
             !
             ! velocity is zero at closed boundary
             !
             v1(nmlu) = -v1(nml)
             !
          else if ( ibr(n) > 1 .and. LMXL ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             v1(nmlu) = v1(nml)
             !
          endif
          !
       enddo
       !
    else
       !
       call periodic ( v1, kgrpnt, 1, 1 )
       !
    endif
    !
    if ( .not.lrepty ) then
       !
       do m = mf, ml
          !
          nfm  = kgrpnt(m,nf )
          nfum = kgrpnt(m,nfu)
          nlm  = kgrpnt(m,nl )
          nlum = kgrpnt(m,nlu)
          !
          if ( ibb(m) == 1 .and. LMYF ) then
             !
             ! velocity is zero at closed boundary
             !
             u1(nfm) = -u1(nfum)
             !
          else if ( ibb(m) > 1 .and. LMYF ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             u1(nfm) = u1(nfum)
             !
          endif
          !
          if ( ibt(m) == 1 .and. LMYL ) then
             !
             ! velocity is zero at closed boundary
             !
             u1(nlum) = -u1(nlm)
             !
          else if ( ibt(m) > 1 .and. LMYL ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             u1(nlum) = u1(nlm)
             !
          endif
          !
       enddo
       !
    else
       !
       call periodic ( u1, kgrpnt, 1, 1 )
       !
    endif
    !
    ! compute the mass flux in x-direction
    !
    do n = nfu, nl
       do m = mf, ml
          !
          nm = kgrpnt(m,n)
          !
          qx(nm) = hu(nm) * guu(nm) * ( tetau(nm)*u1(nm) + (1.-tetau(nm))*u0(nm) )
          !
       enddo
    enddo
    !
    ! compute the mass flux in y-direction
    !
    do n = nf, nl
       do m = mfu, ml
          !
          nm = kgrpnt(m,n)
          !
          qy(nm) = hv(nm) * gvv(nm) * ( tetav(nm)*v1(nm) + (1.-tetav(nm))*v0(nm) )
          !
       enddo
    enddo
    !
    ! impose Neumann / periodic condition for w-velocity at boundaries, if appropriate
    !
    if ( horwinc .and. ihydro == 1 ) then
       !
       if ( .not.lreptx ) then
          !
          do n = nfu, nl
             !
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             nml  = kgrpnt(ml ,n)
             nmlu = kgrpnt(mlu,n)
             !
             if ( LMXF ) w1top(nmf ) = w1top(nmfu)
             if ( LMXL ) w1top(nmlu) = w1top(nml )
             !
          enddo
          !
       else
          !
          call periodic ( w1top, kgrpnt, 1, 1 )
          !
       endif
       !
       if ( .not.lrepty ) then
          !
          do m = mfu, ml
             !
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             nlm  = kgrpnt(m,nl )
             nlum = kgrpnt(m,nlu)
             !
             if ( LMYF ) w1top(nfm ) = w1top(nfum)
             if ( LMYL ) w1top(nlum) = w1top(nlm )
             !
          enddo
          !
       else
          !
          call periodic ( w1top, kgrpnt, 1, 1 )
          !
       endif
       !
    endif
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    ! calculate net mass outflow based on local continuity equation
    !
    if ( ITEST >= 30 .and. ihydro == 1 ) then
       !
       moutf = 0.
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             if ( nmd == 1 ) nmd = nm
             if ( ndm == 1 ) ndm = nm
             !
             dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
             dyl = 0.5 * ( guu(nm) + guu(nmd) )
             !
             if ( nm > 1 ) moutf = moutf + qx(nm) - qx(nmd) - 0.5 * gsqs(nm) * ( u1(nm) + u1(nmd) ) * ( hu(nm) - hu(nmd) ) / dxl +  &
                                           qy(nm) - qy(ndm) - 0.5 * gsqs(nm) * ( v1(nm) + v1(ndm) ) * ( hv(nm) - hv(ndm) ) / dyl +  &
                                           gsqs(nm) * ( w1top(nm) - w1bot(nm) )
             !
          enddo
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
       do n = nfu, nend
          do m = mfu, mend
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
             if ( nm > 1 ) then
                !
                ! compute displaced volume in the cell center
                if ( presp(nm) == 0 ) vol = vol + s1(nm)
                !
                ! compute potential energy in the cell center
                ener = ener + 0.5 * grav * gsqs(nm) * s1(nm) * s1(nm)
                !
                ! compute kinetic energy in the cell center
                fac1 = ( min( -flos(nm ), s1(nm ) ) + dps(nm ) ) * gsqs(nm )
                fac2 = ( min( -flos(nmd), s1(nmd) ) + dps(nmd) ) * gsqs(nmd)
                ener = ener + 0.125 * ( fac1 + fac2 ) * u1(nmd) * u1(nmd)
                fac2 = ( min( -flos(nmu), s1(nmu) ) + dps(nmu) ) * gsqs(nmu)
                ener = ener + 0.125 * ( fac1 + fac2 ) * u1(nm ) * u1(nm )
                !
                fac2 = ( min( -flos(ndm), s1(ndm) ) + dps(ndm) ) * gsqs(ndm)
                ener = ener + 0.125 * ( fac1 + fac2 ) * v1(ndm) * v1(ndm)
                fac2 = ( min( -flos(num), s1(num) ) + dps(num) ) * gsqs(num)
                ener = ener + 0.125 * ( fac1 + fac2 ) * v1(nm ) * v1(nm )
                !
             endif
             !
          enddo
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
    ! exchange velocities with neighbouring subdomains (if appropriate)
    !
    if ( ihydro == 1 ) then
       !
       call SWEXCHG ( u1, kgrpnt, 1, 1 )
       call SWEXCHG ( v1, kgrpnt, 1, 1 )
       if ( horwinc ) call SWEXCHG ( w1top, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    endif
    !
    ! synchronize velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    if ( ihydro == 1 ) then
       !
       call periodic ( u1, kgrpnt, 1, 1 )
       call periodic ( v1, kgrpnt, 1, 1 )
       if ( horwinc ) call periodic ( w1top, kgrpnt, 1, 1 )
       !
    endif
    !
    ! apply wave absorption by means of sponge layer, if appropriate
    !
    if ( spwidl > 0. ) then
       !
       tbndx = 0.
       if ( LMXF ) then
          do n = nfu, nl
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             tbndx(n,1) = u1(nmf )
             tbndx(n,2) = v1(nmfu)
             tbndx(n,3) = s1(nmfu)
          enddo
       endif
       !
       if ( kpart == 4 ) then
          call SWREDUCE ( tbndx, 3*myc, SWREAL, SWSUM )
          if (STPNOW()) return
       endif
       !
       do n = nfu, nl
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,n)
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm) = (1.-sponxl(nm))*u1(nm) + sponxl(nm)*tbndx(n,1)
                !
             endif
             !
          enddo
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmd == 1 ) nmd = nm
             !
             if ( wets(nm) == 1 ) then
                !
                v1(nm) = (1.-sponxl(nmd))*v1(nm) + sponxl(nmd)*tbndx(n,2)
                s1(nm) = (1.-sponxl(nmd))*s1(nm) + sponxl(nmd)*tbndx(n,3)
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    if ( spwidr > 0. ) then
       !
       tbndx = 0.
       if ( LMXL ) then
          do n = nfu, nl
             nml = kgrpnt(ml,n)
             tbndx(n,1) = u1(nml)
             tbndx(n,2) = v1(nml)
             tbndx(n,3) = s1(nml)
          enddo
       endif
       !
       if ( kpart == 4 ) then
          call SWREDUCE ( tbndx, 3*myc, SWREAL, SWSUM )
          if (STPNOW()) return
       endif
       !
       do n = nfu, nl
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,n)
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm) = (1.-sponxr(nm))*u1(nm) + sponxr(nm)*tbndx(n,1)
                !
             endif
             !
          enddo
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             if ( wets(nm) == 1 ) then
                !
                v1(nm) = (1.-sponxr(nm))*v1(nm) + sponxr(nm)*tbndx(n,2)
                s1(nm) = (1.-sponxr(nm))*s1(nm) + sponxr(nm)*tbndx(n,3)
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    if ( spwidb > 0. ) then
       !
       tbndy = 0.
       if ( LMYF ) then
          do m = mfu, ml
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             tbndy(m,1) = v1(nfm )
             tbndy(m,2) = u1(nfum)
             tbndy(m,3) = s1(nfum)
          enddo
       endif
       !
       if ( kpart == 3 ) then
          call SWREDUCE ( tbndy, 3*mxc, SWREAL, SWSUM )
          if (STPNOW()) return
       endif
       !
       do m = mfu, ml
          !
          do n = nf, nl
             !
             nm = kgrpnt(m,n)
             !
             if ( wetv(nm) == 1 ) then
                !
                v1(nm) = (1.-sponyb(nm))*v1(nm) + sponyb(nm)*tbndy(m,1)
                !
             endif
             !
          enddo
          !
          do n = nfu, nl
             !
             nd = n - 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( ndm == 1 ) ndm = nm
             !
             if ( wets(nm) == 1 ) then
                !
                u1(nm) = (1.-sponyb(ndm))*u1(nm) + sponyb(ndm)*tbndy(m,2)
                s1(nm) = (1.-sponyb(ndm))*s1(nm) + sponyb(ndm)*tbndy(m,3)
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    if ( spwidt > 0. ) then
       !
       tbndy = 0.
       if ( LMYL ) then
          do m = mfu, ml
             nlm = kgrpnt(m,nl)
             tbndy(m,1) = v1(nlm)
             tbndy(m,2) = u1(nlm)
             tbndy(m,3) = s1(nlm)
          enddo
       endif
       !
       if ( kpart == 3 ) then
          call SWREDUCE ( tbndy, 3*mxc, SWREAL, SWSUM )
          if (STPNOW()) return
       endif
       !
       do m = mfu, ml
          !
          do n = nf, nl
             !
             nm = kgrpnt(m,n)
             !
             if ( wetv(nm) == 1 ) then
                !
                v1(nm) = (1.-sponyt(nm))*v1(nm) + sponyt(nm)*tbndy(m,1)
                !
             endif
             !
          enddo
          !
          do n = nfu, nl
             !
             nm = kgrpnt(m,n)
             !
             if ( wets(nm) == 1 ) then
                !
                u1(nm) = (1.-sponyt(nm))*u1(nm) + sponyt(nm)*tbndy(m,2)
                s1(nm) = (1.-sponyt(nm))*s1(nm) + sponyt(nm)*tbndy(m,3)
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! compute the maximum CFL number
    !
    cflmax = -999.
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          md = m - 1
          nd = n - 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          ndm = kgrpnt(m ,nd)
          !
          ! for permanently dry neighbours, corresponding values will be mirrored
          !
          if ( nmd == 1 ) nmd = nm
          if ( ndm == 1 ) ndm = nm
          !
          if ( nm > 1 ) then
             !
             u    = 0.5 * ( u1(nm) + u1(nmd) )
             v    = 0.5 * ( v1(nm) + v1(ndm) )
             utot = sqrt( u*u + v*v )
             !
             dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
             dyl = 0.5 * ( guu(nm) + guu(nmd) )
             !
             cfl = dt * utot * sqrt(1./(dxl*dxl) + 1./(dyl*dyl))
             if ( cfl > cflmax ) cflmax = cfl
             !
          endif
          !
       enddo
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
end subroutine SwashImpDepM2DHflow
