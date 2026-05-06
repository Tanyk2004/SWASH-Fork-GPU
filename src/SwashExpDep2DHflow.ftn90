subroutine SwashExpDep2DHflow ( u1, u0, ua, ua2, up, qx, qxm, v1, v0, qy, qym, q, dq, gmatu, gmatv, rho )
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
!    8.01: Marcel Zijlema
!   10.05: Dirk Rijnsdorp
!
!   Updates
!
!    1.00,   April 2010: New subroutine
!    1.10,   March 2012: extension porous media
!    1.15, October 2012: extension repeating grid
!    6.01,    June 2019: extension internal wave generation
!    8.01, October 2021: Hancock scheme for global continuity equation
!   10.05,    July 2023: extension ambient current
!
!   Purpose
!
!   Performs the time integration for the non-hydrostatic, depth-averaged 2D shallow water equations
!
!   Method
!
!   The time integration with respect to the continuity equation and the water level gradient of the
!   momentum equations is based on the leap-frog technique (or the Hansen scheme).
!
!   The time integration with respect to the advective term is based on the predictor-corrector scheme
!   of MacCormack, while that for the bottom friction is based on Euler implicit and for the
!   non-hydrostatic pressure gradient a semi-implicit approach is employed (theta-scheme).
!
!   The space discretization of the advective term is momentum conservative (or energy head conservative
!   in case of flow contraction) and is approximated by either first order upwind or higher order
!   (flux-limited) scheme (CDS, Fromm, BDF, QUICK, MUSCL, Koren, etc.). The higher order scheme is
!   treated by defect correction consistent with the MacCormack scheme.
!
!   The w-momentum equation only contains the z-gradient of the non-hydrostatic pressure and is
!   discretized by means of the Keller-box scheme.
!
!   The non-hydrostatic pressure is obtained by means of the second order accurate pressure correction technique.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
    use m_genarr, only: kgrpnt, guu, guv, gvu, gvv, gsqs, gsqsu, gsqsv, sbimp, slimp, srimp, stimp, work
    use m_parall
    use SwashFlowdata, rhotmp   => rho  , &
                       u1tmp    => u1   , &
                       u0tmp    => u0   , &
                       uatmp    => ua   , &
                       ua2tmp   => ua2  , &
                       uptmp    => up   , &
                       qxtmp    => qx   , &
                       qxmtmp   => qxm  , &
                       v1tmp    => v1   , &
                       v0tmp    => v0   , &
                       qytmp    => qy   , &
                       qymtmp   => qym  , &
                       qtmp     => q    , &
                       dqtmp    => dq   , &
                       gmatutmp => gmatu, &
                       gmatvtmp => gmatv
!
    implicit none
!
!   Argument variables
!
    real   , dimension(mcgrd)  , intent(out)   :: dq       ! non-hydrostatic pressure correction
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
    real   , dimension(mcgrd)  , intent(out)   :: ua       ! advective velocity based on finite differencing
    real   , dimension(mcgrd)  , intent(out)   :: ua2      ! another advective velocity based on finite differencing
    real   , dimension(mcgrd)  , intent(out)   :: up       ! upwind velocity (for energy head computation)
    real   , dimension(mcgrd)  , intent(in)    :: v0       ! v-velocity at previous time level
    real   , dimension(mcgrd)  , intent(inout) :: v1       ! v-velocity at current time level
!
!   Local variables
!
    integer, save                              :: ient = 0 ! number of entries in this subroutine
    integer                                    :: j        ! loop counter
    integer                                    :: m        ! loop counter
    integer                                    :: md       ! index of point m-1
    integer                                    :: mdd      ! index of point m-2
    integer                                    :: mend     ! end index of loop over u-points
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
    integer                                    :: numu     ! pointer to m+1,n+1
    integer                                    :: nuu      ! index of point n+2
    integer                                    :: nuum     ! pointer to m,n+2
    !
    real                                       :: cfl      ! CFL number
    real                                       :: contrib  ! total contribution of momentum equation
    real                                       :: cporfr   ! flow resistance terms inside porous medium
    real                                       :: cveg     ! vegetation term
    real                                       :: denom    ! a denominator
    real                                       :: dep      ! local water depth
    real                                       :: dxl      ! local mesh size in x-direction
    real                                       :: dyl      ! local mesh size in y-direction
    real                                       :: ener     ! total energy of closed system
    real                                       :: epsab2   ! small amount to modify weights of AB2 scheme
    real                                       :: fac      ! a factor
    real                                       :: fac1     ! another factor
    real                                       :: fac2     ! some other factor
    real                                       :: fh       ! factor for Hancock/RK3 scheme
    real                                       :: fluxlim  ! flux limiter
    real                                       :: grad1    ! solution gradient
    real                                       :: grad2    ! another solution gradient
    real                                       :: gsqsd    ! Jacobian in depth point (m+1/2, n+1/2)
    real                                       :: gsqsdb   ! Jacobian in depth point (m+1/2, n-1/2)
    real                                       :: gsqsdl   ! Jacobian in depth point (m-1/2, n+1/2)
    real                                       :: hd       ! water depth in depth point (m+1/2, n+1/2)
    real                                       :: hdb      ! water depth in depth point (m+1/2, n-1/2)
    real                                       :: hdl      ! water depth in depth point (m-1/2, n+1/2)
    real                                       :: moutf    ! net mass outflow
    real                                       :: rxxs     ! component of Reynolds stress tensor Rxx in current wl-point
    real                                       :: rxxsu    ! component of Reynolds stress tensor Rxx in wl-point right from point of consideration
    real                                       :: rxxu     ! component of Reynolds stress tensor Rxx in current u/v-point
    real                                       :: rxyd     ! component of Reynolds stress tensor Rxy in current depth point
    real                                       :: rxydd    ! component of Reynolds stress tensor Rxy in depth point left/down from point of consideration
    real                                       :: rxyu     ! component of Reynolds stress tensor Rxy in current u/v-point
    real                                       :: ryys     ! component of Reynolds stress tensor Ryy in current wl-point
    real                                       :: ryysu    ! component of Reynolds stress tensor Ryy in wl-point up from point of consideration
    real                                       :: ryyu     ! component of Reynolds stress tensor Ryy in current u/v-point
    real                                       :: theta    ! implicitness factor for non-hydrostatic pressure gradient
    real                                       :: thetamb  ! implicitness factor for ambient current in global continuity eq
    real                                       :: u        ! u-velocity at point different from its point of definition
    real                                       :: utot     ! velocity magnitude
    real                                       :: v        ! v-velocity at point different from its point of definition
    real                                       :: vol      ! total displaced volume of water
    real                                       :: zgrad    ! water level gradient term
    !
    logical                                    :: lmask    ! masks points for including curvature stress terms
    logical                                    :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashExpDep2DHflow')
    !
    theta   = pnums( 5)
    thetamb = pnums(71)
    !
    sbimp = .false.
    slimp = .false.
    srimp = .false.
    stimp = .false.
    !
    ! check if water level is prescribed at the boundaries
    !
    do n = nfu, nl
       !
       if ( ibl(n) == 2 .or. ibl(n) == 10 ) then
          !
          slimp(n) = .true.
          !
       endif
       !
       if ( ibr(n) == 2 .or. ibr(n) == 10 ) then
          !
          srimp(n) = .true.
          !
       endif
       !
    enddo
    !
    do m = mfu, ml
       !
       if ( ibb(m) == 2 .or. ibb(m) == 10 ) then
          !
          sbimp(m) = .true.
          !
       endif
       !
       if ( ibt(m) == 2 .or. ibt(m) == 10 ) then
          !
          stimp(m) = .true.
          !
       endif
       !
    enddo
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
    ! compute the advective u-velocity in wl-point based on first order upwind
    !
    do n = nfu, nl
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m ,n)
          nmd = kgrpnt(md,n)
          !
          if ( ( u0(nm) > u0(nmd) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac = 0.5 * ( u0(nm) + u0(nmd) )
             !
          else
             !
             fac = qxm(nm)
             !
          endif
          !
          if ( .not. fac < 0. ) then
             !
             ua(nm) = u0(nmd)
             !
          else
             !
             ua(nm) = u0(nm )
             !
          endif
          !
       enddo
    enddo
    !
    ! synchronize advective u-velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( ua, kgrpnt, 1, 1 )
    !
    ! compute advection term in x-direction (momentum/energy head conservative) in u-momentum equation (prediction)
    !
    msta = mf + 1               ! first internal u-point
    if ( lreptx ) then
       mend = ml                ! last  internal u-point in case of repeating grid
    else
       mend = ml - 1            ! last  internal u-point
    endif
    !
    do n = nfu, nl
       !
       do m = msta, mend
          !
          md = m - 1
          mu = m + 1
          !
          nm  = kgrpnt(m ,n)
          nmd = kgrpnt(md,n)
          nmu = kgrpnt(mu,n)
          !
          if ( linswe ) then
             !
             fac1 = 0.
             fac2 = 0.
             !
          else if ( ( u0(nm) > u0(nmd) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac1 = 0.5 * ( ua(nmu) + ua(nm) )
             fac2 = fac1
             !
          else if ( humn(nm) > 0. ) then
             !
             fac1 = qxm(nmu) / humn(nm)
             fac2 = qxm(nm ) / humn(nm)
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
             advecx(nm) = ( fac1 * (ua(nmu) - u0(nm)) - fac2 * (ua(nm) - u0(nm)) ) / gvu(nm)
             !
          endif
          !
       enddo
       !
       if ( slimp(n) .and. .not.linswe ) then    ! described water level at left boundary point
          !
          nmf  = kgrpnt(mf ,n)
          nmfu = kgrpnt(mfu,n)
          !
          if ( .not. qxm(nmfu) < 0. ) then
             !
             advecx(nmf) = 0.
             !
          else
             !
             advecx(nmf) = qxm(nmfu) * (u0(nmfu) - u0(nmf)) / ( gvu(nmf) * hum(nmf) )
             !
          endif
          !
       endif
       !
       if ( srimp(n) .and. .not.linswe ) then    ! described water level at right boundary point
          !
          nml  = kgrpnt(ml  ,n)
          nmld = kgrpnt(ml-1,n)
          !
          if ( .not. qxm(nml) < 0. ) then
             !
             advecx(nml) = qxm(nml) * (u0(nml) - u0(nmld)) / ( gvu(nml) * hum(nml) )
             !
          else
             !
             advecx(nml) = 0.
             !
          endif
          !
       endif
       !
    enddo
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
    ! compute the advective u-velocity in depth point based on first order upwind
    !
    do n = nf, nl
       do m = msta, mend
          !
          md = m - 1
          mu = m + 1
          nu = n + 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          nmu = kgrpnt(mu,n )
          num = kgrpnt(m ,nu)
          !
          if ( ( u0(nm) > u0(nmd) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac = 0.5 * ( v0(nm) + v0(nmu) )
             !
          else
             !
             fac = qym(nm)
             !
          endif
          !
          if ( .not. fac < 0. ) then
             !
             ua(nm) = u0(nm )
             !
          else
             !
             ua(nm) = u0(num)
             !
          endif
          !
       enddo
    enddo
    !
    ! compute advection term in y-direction (momentum/energy head conservative) in u-momentum equation (prediction)
    !
    do n = nfu, nl
       !
       do m = msta, mend
          !
          md = m - 1
          mu = m + 1
          nd = n - 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          nmu = kgrpnt(mu,n )
          ndm = kgrpnt(m ,nd)
          !
          if ( linswe ) then
             !
             fac1 = 0.
             fac2 = 0.
             !
          else if ( ( u0(nm) > u0(nmd) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac1 = 0.5 * ( v0(nm) + v0(nmu) )
             fac2 = fac1
             !
          else if ( humn(nm) > 0. ) then
             !
             fac1 = qym(nm ) / humn(nm)
             fac2 = qym(ndm) / humn(nm)
             !
          else
             !
             fac1 = 0.
             fac2 = 0.
             !
          endif
          !
          ! to avoid boundary layer effect, set advection near lateral open boundary to zero
          !
          if ( ( .not.lrepty .and. LMYF .and. n == nfu ) .or. ( .not.lrepty .and. LMYL .and. n == nl ) ) then
             !
             fac1 = 0.
             fac2 = 0.
             !
          endif
          !
          if ( wetu(nm) == 1 ) then
             !
             advecy(nm) = ( fac1 * (ua(nm) - u0(nm)) - fac2 * (ua(ndm) - u0(nm)) ) / guu(nm)
             !
          endif
          !
       enddo
       !
       if ( slimp(n) ) then    ! described water level at left boundary point
          !
          nmf = kgrpnt(mf,n)
          !
          advecy(nmf) = 0.
          !
       endif
       !
       if ( srimp(n) ) then    ! described water level at right boundary point
          !
          nml = kgrpnt(ml,n)
          !
          advecy(nml) = 0.
          !
       endif
       !
    enddo
    !
    ! compute the curvature terms associated with advection terms (internal points only!)
    !
    if ( optg == 3 .and. .not.linswe ) then
       !
       do n = nfu, nl
          do m = msta, mend
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             !
             nm   = kgrpnt(m ,n )
             nmd  = kgrpnt(md,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             ndmu = kgrpnt(mu,nd)
             !
             if ( wetu(nm) == 1 ) then
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd  == 1 ) nmd  = nm
                if ( nmu  == 1 ) nmu  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmu == 1 ) ndmu = ndm
                !
                v = 0.25 * ( v0(nm) + v0(nmu) + v0(ndm) + v0(ndmu) )
                !
                if ( advecx(nm) /= 0. ) advecx(nm) = advecx(nm) + 0.5 * u0(nm) * v * ( gvv(nm) + gvv(nmu) - gvv(ndm) - gvv(ndmu) ) / gsqsu(nm)
                !
                if ( advecy(nm) /= 0. ) advecy(nm) = advecy(nm) - 0.5 * v * v * ( guu(nmu) - guu(nmd) ) / gsqsu(nm)
                !
             endif
             !
          enddo
       enddo
      !
    endif
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
             if ( wetu(nm) * wetu(nmu) == 1 .and. (m /= mf .or. slimp(n) .or. .not.LMXF) .and. (m /= ml .or. srimp(n) .or. .not.LMXL) ) then
                !
                fac = 1. / ( 2. * gvu(nm) * hum(nm) )
                !
                gmatu(nm,1) = -fac * ( s0(nm ) + dps(nmu) )
                gmatu(nm,2) =  fac * ( s0(nmu) + dps(nm ) )
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
       if ( iproj == 2 ) qgrad = (1.-theta) * qgrad
       !
    else
       !
       qgrad = 0.
       !
    endif
    !
    ! compute viscosity term
    !
    if ( ihvisc == 1 .and. hvisc > 0. ) then
       !
       do n = nfu, nl
          !
          nmf = kgrpnt(mf,n)
          nml = kgrpnt(ml,n)
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
             if ( wetu(nm) == 1 ) then
                !
                if ( wetu(nmd) * wetu(nmu) == 1 ) then
                   !
                   fac1 = hs(nmu) * ( guu(nmu)*u0(nmu) - guu(nm)*u0(nm) )/gsqs(nmu) - hs(nm) * ( guu(nm)*u0(nm) - guu(nmd)*u0(nmd) )/gsqs(nm)
                   !
                else
                   !
                   fac1 = 0.
                   !
                endif
                !
                if ( wetu(ndm) * wetu(num) == 1 ) then
                   !
                   gsqsd  = .25 * (guu(nm ) + guu(num)) * (gvv(nm ) + gvv(nmu ))
                   gsqsdb = .25 * (guu(ndm) + guu(nm )) * (gvv(ndm) + gvv(ndmu))
                   !
                   hd  = .5 * (hum(num) + hum(nm ))
                   hdb = .5 * (hum(nm ) + hum(ndm))
                   !
                   fac2 = hd * ( gvu(num)*u0(num) - gvu(nm)*u0(nm) )/gsqsd - hdb * ( gvu(nm)*u0(nm) - gvu(ndm)*u0(ndm) )/gsqsdb
                   !
                else
                   !
                   fac2 = 0.
                   !
                endif
                !
                visc(nm) = hvisc * ( fac1/gvu(nm) + fac2/guu(nm) ) / hum(nm)
                !
             endif
             !
          enddo
          !
          if ( slimp(n) ) visc(nmf) = 0.
          !
          if ( srimp(n) ) visc(nml) = 0.
          !
       enddo
       !
    else if ( ihvisc > 1 ) then
       !
       do n = nfu, nl
          !
          nmf = kgrpnt(mf,n)
          nml = kgrpnt(ml,n)
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
             ndmd = kgrpnt(md,nd)
             ndmu = kgrpnt(mu,nd)
             !
             if ( wetu(nm) == 1 ) then
                !
                lmask = .true.
                !
                hd  = .5 * (hum(num) + hum(nm ))
                hdb = .5 * (hum(nm ) + hum(ndm))
                !
                ! first, compute the components of the Reynolds stress tensor in different points of consideration
                !
                if ( wetu(nmd) * wetu(nmu) == 1 ) then
                   !
                   rxxs  = ( vnu2d(nm ) + vnu2d(ndm ) + vnu2d(nmd) + vnu2d(ndmd) ) * ( u0(nm ) - u0(nmd) ) / ( gvv(nm ) + gvv(ndm ) )
                   !
                   rxxsu = ( vnu2d(nmu) + vnu2d(ndmu) + vnu2d(nm ) + vnu2d(ndm ) ) * ( u0(nmu) - u0(nm ) ) / ( gvv(nmu) + gvv(ndmu) )
                   !
                else
                   !
                   rxxs  = 0.
                   rxxsu = 0.
                   !
                   lmask = .false.
                   !
                endif
                !
                if ( wetu(ndm) * wetu(num) * wetv(nm) * wetv(ndm) * wetv(nmu) * wetv(ndmu) == 1 ) then
                   !
                   rxyd  = 2. * vnu2d(nm ) * ( ( v0(nmu ) - v0(nm ) )/( gvv(nmu ) + gvv(nm ) ) + ( u0(num) - u0(nm ) )/( guu(num) + guu(nm ) ) )
                   !
                   rxydd = 2. * vnu2d(ndm) * ( ( v0(ndmu) - v0(ndm) )/( gvv(ndmu) + gvv(ndm) ) + ( u0(nm ) - u0(ndm) )/( guu(nm ) + guu(ndm) ) )
                   !
                else
                   !
                   rxyd  = 0.
                   rxydd = 0.
                   !
                   lmask = .false.
                   !
                endif
                !
                if ( lmask ) then
                   !
                   rxxu  = 0.50 * ( vnu2d(nm) + vnu2d(ndm) ) * ( u0(nmu) - u0(nmd) ) / gvu(nm)
                   !
                   rxyu  = 0.25 * ( vnu2d(nm) + vnu2d(ndm) ) * ( ( v0(nmu) + v0(ndmu) - v0(nm) - v0(ndm) )/gvu(nm) + ( u0(num) - u0(ndm) )/guu(nm) )
                   !
                   ryyu  = 0.50 * ( vnu2d(nm) + vnu2d(ndm) ) * ( v0(nmu) + v0(nm) - v0(ndm) -v0(ndmu) ) / guu(nm)
                   !
                else
                   !
                   rxxu = 0.
                   rxyu = 0.
                   ryyu = 0.
                   !
                endif
                !
                ! next, compute the divergence of the Reynolds stresses
                !
                visc(nm) = ( hs(nmu)*rxxsu - hs(nm)*rxxs ) / (hum(nm)*gvu(nm)) + ( hd*rxyd - hdb*rxydd ) / (hum(nm)*guu(nm))
                !
                ! finally, add the curvature terms associated with those stresses
                !
                visc(nm) = visc(nm) + ( rxyu * ( gvv(nm) + gvv(nmu) - gvv(ndm) - gvv(ndmu) ) + 0.5 * (rxxu - ryyu) * ( guu(nmu) - guu(nmd) ) ) / gsqsu(nm)
                !
             endif
             !
          enddo
          !
          if ( slimp(n) ) visc(nmf) = 0.
          !
          if ( srimp(n) ) visc(nml) = 0.
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
    ! compute Coriolis force
    !
    if ( coriolis ) then
       !
       epsab2 = pcor(2)
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
             if ( wetu(nm) == 1 .and. fcor(nm,1) /= 0. ) then
                !
                ! at time level n-1
                !
                ! --- modified quasi second order Adams-Bashforth scheme which is conditionally stable for the inviscid case (see Marshall et al., 1997)
                !
                corf(nm) = -(0.5+epsab2) * cfu(nm,1)
                !
                ! at time level n
                !
                ! --- interpolate v-velocity at u-point using energy-conserved weighting factors according to Espelid et al. (2000)
                !
                v = sqrt(abs(fcor(nm,2))*dpv(nm))*v0(nm) + sqrt(abs(fcor(nmu,2))*dpv(nmu))*v0(nmu) + sqrt(abs(fcor(ndm,2))*dpv(ndm))*v0(ndm) + sqrt(abs(fcor(ndmu,2))*dpv(ndmu))*v0(ndmu)
                v = v / (4.*sqrt(abs(fcor(nm,1))*dpu(nm)))
                !
                cfu(nm,1) = -fcor(nm,1) * v
                !
                corf(nm) = corf(nm) + (1.5+epsab2) * cfu(nm,1)
                !
             endif
             !
          enddo
       enddo
       !
    else
       !
       corf = 0.
       !
    endif
    !
    ! include effect of ambient current - first order approximation
    !
    if ( iamb == 1 ) then
       !
       ! compute the advective u-velocity in wl-point based on first order upwind
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             !
             if ( .not. Ubck(nm,1) < 0. ) then
                !
                ua(nm) = u0(nmd)
                !
             else
                !
                ua(nm) = u0(nm )
                !
             endif
             !
          enddo
       enddo
       !
       ! synchronize advective u-velocities at appropriate boundaries in case of repeating grid (if appropriate)
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       !
       ! compute advection term in x-direction (momentum conservative) at internal u-point (prediction)
       !
       do n = nfu, nl
          !
          do m = mf+1, mend
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( wetu(nm) == 1 ) then
                !
                advecx(nm) = advecx(nm) + ( Ubck(nmu,1) * ua(nmu) - Ubck(nm,1) * ua(nm) ) / gvu(nm)
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! compute the advective u-velocity in depth point based on first order upwind
       !
       do n = nf, nl
          do m = mf+1, mend
             !
             mu = m + 1
             nu = n + 1
             !
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             num  = kgrpnt(m ,nu)
             numu = kgrpnt(mu,nu)
             !
             fac = 0.25 * ( Vbck(nm,1) + Vbck(nmu,1) + Vbck(num,1) + Vbck(numu,1) )
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm) = u0(nm )
                !
             else
                !
                ua(nm) = u0(num)
                !
             endif
             !
          enddo
       enddo
       !
       ! compute advection term in y-direction (momentum conservative) at internal u-point (prediction)
       !
       do n = nfu, nl
          !
          do m = mf+1, mend
             !
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             num  = kgrpnt(m ,nu)
             ndmu = kgrpnt(mu,nd)
             numu = kgrpnt(mu,nu)
             !
             fac1 = 0.25 * ( Vbck(nm ,1) + Vbck(nmu ,1) + Vbck(num,1) + Vbck(numu,1) )
             fac2 = 0.25 * ( Vbck(ndm,1) + Vbck(ndmu,1) + Vbck(nm ,1) + Vbck(nmu ,1) )
             !
             ! to avoid boundary layer effect, set advection near lateral open boundary to zero
             !
             if ( ( .not.lrepty .and. LMYF .and. n == nfu ) .or. ( .not.lrepty .and. LMYL .and. n == nl ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetu(nm) == 1 ) then
                !
                advecy(nm) = advecy(nm) + ( fac1 * ( ua(nm) - u0(nm) ) - fac2 * ( ua(ndm) - u0(nm) ) ) / guu(nm)
                !
             endif
             !
             ! also add v dU/dy
             !
             v = 0.25 * ( v0(nm) + v0(nmu) + v0(ndm) + v0(ndmu) )
             if ( ( .not.lrepty .and. LMYF .and. n == nfu ) .or. ( .not.lrepty .and. LMYL .and. n == nl ) ) v = 0.
             if ( wetu(nm) == 1 ) advecy(nm) = advecy(nm) + 0.25 * v * ( Ubck(numu,1) + Ubck(num,1) - Ubck(ndmu,1) - Ubck(ndm,1) ) / guu(nm)
             !
          enddo
          !
       enddo
       !
    else if ( iamb == 2 ) then
       !
       ! compute the advective u-velocity in wl-point based on first order upwind
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             !
             fac = 0.5 * ( Ubck(nm,1) + Ubck(nmd,1) )
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm) = u0(nmd)
                !
             else
                !
                ua(nm) = u0(nm )
                !
             endif
             !
          enddo
       enddo
       !
       ! synchronize advective u-velocities at appropriate boundaries in case of repeating grid (if appropriate)
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       !
       ! compute advection term in x-direction (momentum conservative) at internal u-point (prediction)
       !
       do n = nfu, nl
          !
          do m = mf+1, mend
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             nmu = kgrpnt(mu,n)
             !
             fac1 = 0.5 * ( Ubck(nmu,1) + Ubck(nm ,1) )
             fac2 = 0.5 * ( Ubck(nm ,1) + Ubck(nmd,1) )
             !
             if ( wetu(nm) == 1 ) then
                !
                advecx(nm) = advecx(nm) + ( fac1 * ua(nmu) - fac2 * ua(nm) ) / gvu(nm)
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! compute the advective u-velocity in depth point based on first order upwind
       !
       do n = nf, nl
          do m = mf+1, mend
             !
             mu = m + 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmu = kgrpnt(mu,n )
             num = kgrpnt(m ,nu)
             !
             fac = 0.5 * ( Vbck(nm,1) + Vbck(nmu,1) )
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm) = u0(nm )
                !
             else
                !
                ua(nm) = u0(num)
                !
             endif
             !
          enddo
       enddo
       !
       !
       ! compute advection term in y-direction (momentum conservative) at internal u-point (prediction)
       !
       do n = nfu, nl
          !
          do m = mf+1, mend
             !
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             ndmu = kgrpnt(mu,nd)
             num  = kgrpnt(m ,nu)
             !
             fac1 = 0.5 * ( Vbck(nm ,1) + Vbck(nmu ,1) )
             fac2 = 0.5 * ( Vbck(ndm,1) + Vbck(ndmu,1) )
             !
             ! to avoid boundary layer effect, set advection near lateral open boundary to zero
             !
             if ( ( .not.lrepty .and. LMYF .and. n == nfu ) .or. ( .not.lrepty .and. LMYL .and. n == nl ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetu(nm) == 1 ) then
                !
                advecy(nm) = advecy(nm) + ( fac1 * ( ua(nm) - u0(nm) ) - fac2 * ( ua(ndm) - u0(nm) ) ) / guu(nm)
                !
             endif
             !
             ! also add v dU/dy
             !
             v = 0.25 * ( v0(nm) + v0(nmu) + v0(ndm) + v0(ndmu) )
             if ( ( .not.lrepty .and. LMYF .and. n == nfu ) .or. ( .not.lrepty .and. LMYL .and. n == nl ) ) v = 0.
             if ( wetu(nm) == 1 ) advecy(nm) = advecy(nm) + 0.5 * v * ( Ubck(num,1) - Ubck(ndm,1) ) / guu(nm)
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! compute intermediate u-velocity
    !
    do n = nfu, nl
       !
       nmf = kgrpnt(mf,n)
       nml = kgrpnt(ml,n)
       !
       if ( slimp(n) ) then       ! described water level at left boundary point
          !
          msta = mf
          !
          ! neglect wind, pressure gradients and Coriolis force at water level boundary
          !
          windu (nmf) = 0.
          wndimp(nmf) = 0.
          pgrad (nmf) = 0.
          corf  (nmf) = 0.
          !
       else
          !
          msta = mf + 1
          !
       endif
       !
       if ( srimp(n) ) then       ! described water level at right boundary point
          !
          mend = ml
          !
          ! neglect wind, pressure gradients and Coriolis force at water level boundary
          !
          windu (nml) = 0.
          wndimp(nml) = 0.
          pgrad (nml) = 0.
          corf  (nml) = 0.
          !
       else if ( lreptx ) then    ! periodicity at right boundary point
          !
          mend = ml
          !
       else
          !
          mend = ml - 1
          !
       endif
       !
       do m = msta, mend
          !
          mu = m + 1
          nd = n - 1
          !
          nm   = kgrpnt(m ,n )
          nmu  = kgrpnt(mu,n )
          ndm  = kgrpnt(m ,nd)
          ndmu = kgrpnt(mu,nd)
          !
          if ( wetu(nm) == 1 .and. nmu /= 1 ) then
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
             ! compute total contributions of the momentum equation
             !
             contrib = advecx(nm) + advecy(nm) + zgrad + pgrad(nm) + qgrad(nm) - windu(nm)/max(1.e-3,hum(nm)) - visc(nm) + corf(nm)
             !
             denom = 1. + cvegu(nm,1,2) + cpomu(nm,1) + dt * ( cbot(nm) + cporfr + cveg + wndimp(nm) )
             !
             ! compute flow velocity
             !
             u1(nm) = ( (1.+cvegu(nm,1,2)+cpomu(nm,1)) * u0(nm) - dt * contrib ) / denom
             !
          else
             !
             u1(nm) = 0.
             !
          endif
          !
       enddo
       !
    enddo
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
       do n = nfu, nl
          do m = mfu, ml
             !
             md  = m  - 1
             mdd = md - 1
             mu  = m  + 1
             !
             if ( LMXF .and. mdd < mf ) mdd = mf
             if ( LMXL .and. mu  > ml ) mu  = ml
             !
             nm   = kgrpnt(m  ,n)
             nmd  = kgrpnt(md ,n)
             nmdd = kgrpnt(mdd,n)
             nmu  = kgrpnt(mu ,n)
             !
             if ( ( u0(nm) > u0(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac = 0.5 * ( u0(nm) + u0(nmd) )
                !
             else
                !
                fac = qxm(nm)
                !
             endif
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = u1(nm ) - u0(nmd )
                grad2 = u1(nmd) - u0(nmdd)
                !
                up(nm) = u0(nmd)
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = u1(nmd) - u0(nm )
                grad2 = u1(nm ) - u0(nmu)
                !
                up(nm) = u0(nm)
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       call periodic ( up, kgrpnt, 1, 1 )
       !
    endif
    !
    ! compute higher order correction to the advective velocity in depth point (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       msta = mf + 1            ! first internal u-point
       if ( lreptx ) then
          mend = ml             ! last  internal u-point in case of repeating grid
       else
          mend = ml - 1         ! last  internal u-point
       endif
       !
       do n = nf, nl
          do m = msta, mend
             !
             md  = m  - 1
             mu  = m  + 1
             nd  = n  - 1
             nu  = n  + 1
             nuu = nu + 1
             !
             if ( LMYF .and. nd  < nfu ) nd  = nfu
             if ( LMYL .and. nu  > nl  ) nu  = nl
             if ( LMYL .and. nuu > nl  ) nuu = nl
             !
             nm   = kgrpnt(m ,n  )
             nmd  = kgrpnt(md,n  )
             nmu  = kgrpnt(mu,n  )
             ndm  = kgrpnt(m ,nd )
             num  = kgrpnt(m ,nu )
             nuum = kgrpnt(m ,nuu)
             !
             if ( ( u0(nm) > u0(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac = 0.5 * ( v0(nm) + v0(nmu) )
                !
             else
                !
                fac = qym(nm)
                !
             endif
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = u1(num) - u0(nm )
                grad2 = u1(nm ) - u0(ndm)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = u1(nm ) - u0(num )
                grad2 = u1(num) - u0(nuum)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! correct u-velocity associated with higher order approximation of advection term in x-direction (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do n = nfu, nl
          do m = msta, mend
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             nmu = kgrpnt(mu,n)
             !
             if ( ( u0(nm) > u0(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                !fac1 = up(nmu) + 0.5*ua(nmu)
                !fac2 = up(nm ) + 0.5*ua(nm )
                !jumps back to old implementation (020117)
                fac1 = 0.
                fac2 = 0.
                !
             else if ( humn(nm) > 0. ) then
                !
                fac1 = qxm(nmu) / humn(nm)
                fac2 = qxm(nm ) / humn(nm)
                !
             else
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetu(nm) == 1 .and. nmu /= 1 ) then
                !
                u1(nm) = u1(nm) - dt * ( fac1 * ua(nmu)  - fac2 * ua(nm) ) / gvu(nm)
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! correct u-velocity associated with higher order approximation of advection term in y-direction (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do n = nfu, nl
          do m = msta, mend
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             nmu = kgrpnt(mu,n )
             ndm = kgrpnt(m ,nd)
             !
             if ( ( u0(nm) > u0(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                !fac1 = 0.5 * ( v0(nm) + v0(nmu) )
                !fac2 = fac1
                !jumps back to old implementation (020117)
                fac1 = 0.
                fac2 = 0.
                !
             else if ( humn(nm) > 0. ) then
                !
                fac1 = qym(nm ) / humn(nm)
                fac2 = qym(ndm) / humn(nm)
                !
             else
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             ! set higher order approximation near lateral boundary to zero
             !
             if ( ( LMYF .and. n == nfu ) .or. ( LMYL .and. n == nl ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetu(nm) == 1 .and. nmu /= 1 ) then
                !
                u1(nm) = u1(nm) - dt * ( fac1 * ua2(nm) - fac2 * ua2(ndm) ) / guu(nm)
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! exchange u-velocities with neighbouring subdomains (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       call SWEXCHG ( u1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    endif
    !
    ! synchronize u-velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    if ( propsc /= 1 ) call periodic ( u1, kgrpnt, 1, 1 )
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
       do n = nfu, nl
          do m = mfu, ml
             !
             md  = m  - 1
             mdd = md - 1
             mu  = m  + 1
             !
             if ( LMXF .and. mdd < mf ) mdd = mf
             if ( LMXL .and. mu  > ml ) mu  = ml
             !
             nm   = kgrpnt(m  ,n)
             nmd  = kgrpnt(md ,n)
             nmdd = kgrpnt(mdd,n)
             nmu  = kgrpnt(mu ,n)
             !
             if ( .not. Ubck(nm,1) < 0. ) then
                !
                grad1 = u1(nm ) - u0(nmd )
                grad2 = u1(nmd) - u0(nmdd)
                !
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = u1(nmd) - u0(nm )
                grad2 = u1(nm ) - u0(nmu)
                !
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       !
       ! compute higher order correction to the advective velocity in depth point (if appropriate)
       !
       if ( lreptx ) then
          mend = ml             ! last internal u-point in case of repeating grid
       else
          mend = ml - 1         ! last internal u-point
       endif
       !
       do n = nf, nl
          do m = mf+1, mend
             !
             mu  = m  + 1
             nd  = n  - 1
             nu  = n  + 1
             nuu = nu + 1
             !
             if ( LMYF .and. nd  < nfu ) nd  = nfu
             if ( LMYL .and. nu  > nl  ) nu  = nl
             if ( LMYL .and. nuu > nl  ) nuu = nl
             !
             nm   = kgrpnt(m ,n  )
             nmu  = kgrpnt(mu,n  )
             ndm  = kgrpnt(m ,nd )
             num  = kgrpnt(m ,nu )
             numu = kgrpnt(mu,nu )
             nuum = kgrpnt(m ,nuu)
             !
             fac = 0.25 * ( Vbck(nm,1) + Vbck(nmu,1) + Vbck(num,1) + Vbck(numu,1) )
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = u1(num) - u0(nm )
                grad2 = u1(nm ) - u0(ndm)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = u1(nm ) - u0(num )
                grad2 = u1(num) - u0(nuum)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       ! correct u-velocity associated with higher order approximation of advection term in x-direction at internal u-point (if appropriate)
       !
       do n = nfu, nl
          do m = mf+1, mend
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm) = u1(nm) - dt * ( Ubck(nmu,1) * ua(nmu)  - Ubck(nm,1) * ua(nm) ) / gvu(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! correct u-velocity associated with higher order approximation of advection term in y-direction at internal u-point (if appropriate)
       !
       do n = nfu, nl
          do m = mf+1, mend
             !
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             num  = kgrpnt(m ,nu)
             ndmu = kgrpnt(mu,nd)
             numu = kgrpnt(mu,nu)
             !
             fac1 = 0.25 * ( Vbck(nm ,1) + Vbck(nmu ,1) + Vbck(num,1) + Vbck(numu,1) )
             fac2 = 0.25 * ( Vbck(ndm,1) + Vbck(ndmu,1) + Vbck(nm ,1) + Vbck(nmu ,1) )
             !
             ! set higher order approximation near lateral boundary to zero
             !
             if ( ( LMYF .and. n == nfu ) .or. ( LMYL .and. n == nl ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm) = u1(nm) - dt * ( fac1 * ua2(nm) - fac2 * ua2(ndm) ) / guu(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! exchange u-velocities with neighbouring subdomains (if appropriate)
       !
       call SWEXCHG ( u1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize u-velocities at appropriate boundaries in case of repeating grid (if appropriate)
       !
       call periodic ( u1, kgrpnt, 1, 1 )
       !
    else if ( iamb == 2 .and. propsc /= 1 ) then
       !
       ! compute higher order correction to the advective velocity in wl-point (if appropriate)
       !
       kappa  = pnums(62)
       mbound = pnums(63)
       phieby = pnums(64)
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md  = m  - 1
             mdd = md - 1
             mu  = m  + 1
             !
             if ( LMXF .and. mdd < mf ) mdd = mf
             if ( LMXL .and. mu  > ml ) mu  = ml
             !
             nm   = kgrpnt(m  ,n)
             nmd  = kgrpnt(md ,n)
             nmdd = kgrpnt(mdd,n)
             nmu  = kgrpnt(mu ,n)
             !
             fac = 0.5 * ( Ubck(nm,1) + Ubck(nmd,1) )
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = u1(nm ) - u0(nmd )
                grad2 = u1(nmd) - u0(nmdd)
                !
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = u1(nmd) - u0(nm )
                grad2 = u1(nm ) - u0(nmu)
                !
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       !
       ! compute higher order correction to the advective velocity in depth point (if appropriate)
       !
       if ( lreptx ) then
          mend = ml             ! last internal u-point in case of repeating grid
       else
          mend = ml - 1         ! last internal u-point
       endif
       !
       do n = nf, nl
          do m = mf+1, mend
             !
             mu  = m  + 1
             nd  = n  - 1
             nu  = n  + 1
             nuu = nu + 1
             !
             if ( LMYF .and. nd  < nfu ) nd  = nfu
             if ( LMYL .and. nu  > nl  ) nu  = nl
             if ( LMYL .and. nuu > nl  ) nuu = nl
             !
             nm   = kgrpnt(m ,n  )
             nmu  = kgrpnt(mu,n  )
             ndm  = kgrpnt(m ,nd )
             num  = kgrpnt(m ,nu )
             nuum = kgrpnt(m ,nuu)
             !
             fac = 0.5 * ( Vbck(nm,1) + Vbck(nmu,1) )
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = u1(num) - u0(nm )
                grad2 = u1(nm ) - u0(ndm)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = u1(nm ) - u0(num )
                grad2 = u1(num) - u0(nuum)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       ! correct u-velocity associated with higher order approximation of advection term in x-direction at internal u-point (if appropriate)
       !
       do n = nfu, nl
          do m = mf+1, mend
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             nmu = kgrpnt(mu,n)
             !
             fac1 = 0.5 * ( Ubck(nmu,1) + Ubck(nm ,1) )
             fac2 = 0.5 * ( Ubck(nm ,1) + Ubck(nmd,1) )
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm) = u1(nm) - dt * ( fac1 * ua(nmu)  - fac2 * ua(nm) ) / gvu(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! correct u-velocity associated with higher order approximation of advection term in y-direction at internal u-point (if appropriate)
       !
       do n = nfu, nl
          do m = mf+1, mend
             !
             mu = m + 1
             nd = n - 1
             !
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             ndmu = kgrpnt(mu,nd)
             !
             fac1 = 0.5 * ( Vbck(nm ,1) + Vbck(nmu ,1) )
             fac2 = 0.5 * ( Vbck(ndm,1) + Vbck(ndmu,1) )
             !
             ! set higher order approximation near lateral boundary to zero
             !
             if ( ( LMYF .and. n == nfu ) .or. ( LMYL .and. n == nl ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetu(nm) == 1 ) then
                !
                u1(nm) = u1(nm) - dt * ( fac1 * ua2(nm) - fac2 * ua2(ndm) ) / guu(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! exchange u-velocities with neighbouring subdomains (if appropriate)
       !
       call SWEXCHG ( u1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize u-velocities at appropriate boundaries in case of repeating grid (if appropriate)
       !
       call periodic ( u1, kgrpnt, 1, 1 )
       !
    endif
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
    ! compute the advective v-velocity in depth point based on first order upwind
    !
    do m = mf, ml
       do n = nsta, nend
          !
          nd = n - 1
          nu = n + 1
          mu = m + 1
          !
          nm  = kgrpnt(m ,n )
          ndm = kgrpnt(m ,nd)
          num = kgrpnt(m ,nu)
          nmu = kgrpnt(mu,n )
          !
          if ( ( v0(nm) > v0(ndm) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac = 0.5 * ( u0(nm) + u0(num) )
             !
          else
             !
             fac = qxm(nm)
             !
          endif
          !
          if ( .not. fac < 0. ) then
             !
             ua(nm) = v0(nm )
             !
          else
             !
             ua(nm) = v0(nmu)
             !
          endif
          !
       enddo
    enddo
    !
    ! compute advection term in x-direction (momentum/energy head conservative) in v-momentum equation (prediction)
    !
    do m = mfu, ml
       !
       do n = nsta, nend
          !
          nd = n - 1
          nu = n + 1
          md = m - 1
          !
          nm  = kgrpnt(m ,n )
          ndm = kgrpnt(m ,nd)
          num = kgrpnt(m ,nu)
          nmd = kgrpnt(md,n )
          !
          if ( linswe ) then
             !
             fac1 = 0.
             fac2 = 0.
             !
          else if ( ( v0(nm) > v0(ndm) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac1 = 0.5 * ( u0(nm) + u0(num) )
             fac2 = fac1
             !
          else if ( hvmn(nm) > 0. ) then
             !
             fac1 = qxm(nm ) / hvmn(nm)
             fac2 = qxm(nmd) / hvmn(nm)
             !
          else
             !
             fac1 = 0.
             fac2 = 0.
             !
          endif
          !
          ! to avoid boundary layer effect, set advection near lateral open boundary to zero
          !
          if ( ( .not.lreptx .and. LMXF .and. m == mfu ) .or. ( .not.lreptx .and. LMXL .and. m == ml ) ) then
             !
             fac1 = 0.
             fac2 = 0.
             !
          endif
          !
          if ( wetv(nm) == 1 ) then
             !
             advecx(nm) = ( fac1 * (ua(nm) - v0(nm)) - fac2 * (ua(nmd) - v0(nm)) ) / gvv(nm)
             !
          endif
          !
       enddo
       !
       if ( sbimp(m) ) then    ! described water level at lower boundary point
          !
          nfm = kgrpnt(m,nf)
          !
          advecx(nfm) = 0.
          !
       endif
       !
       if ( stimp(m) ) then    ! described water level at upper boundary point
          !
          nlm = kgrpnt(m,nl)
          !
          advecx(nlm) = 0.
          !
       endif
       !
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
    ! compute the advective v-velocity in wl-point based on first order upwind
    !
    do m = mfu, ml
       do n = nfu, nl
          !
          nd = n - 1
          !
          nm  = kgrpnt(m,n )
          ndm = kgrpnt(m,nd)
          !
          if ( ( v0(nm) > v0(ndm) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac = 0.5 * ( v0(nm) + v0(ndm) )
             !
          else
             !
             fac = qym(nm)
             !
          endif
          !
          if ( .not. fac < 0. ) then
             !
             ua(nm) = v0(ndm)
             !
          else
             !
             ua(nm) = v0(nm )
             !
          endif
          !
       enddo
    enddo
    !
    ! synchronize advective v-velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( ua, kgrpnt, 1, 1 )
    !
    ! compute advection term in y-direction (momentum/energy head conservative) in v-momentum equation (prediction)
    !
    do m = mfu, ml
       !
       do n = nsta, nend
          !
          nd = n - 1
          nu = n + 1
          !
          nm  = kgrpnt(m,n )
          ndm = kgrpnt(m,nd)
          num = kgrpnt(m,nu)
          !
          if ( linswe ) then
             !
             fac1 = 0.
             fac2 = 0.
             !
          else if ( ( v0(nm) > v0(ndm) .or. stricthead ) .and. .not.strictmom ) then
             !
             fac1 = 0.5 * ( ua(num) + ua(nm) )
             fac2 = fac1
             !
          else if ( hvmn(nm) > 0. ) then
             !
             fac1 = qym(num) / hvmn(nm)
             fac2 = qym(nm ) / hvmn(nm)
             !
          else
             !
             fac1 = 0.
             fac2 = 0.
             !
          endif
          !
          if ( wetv(nm) == 1 ) then
             !
             advecy(nm) = ( fac1 * (ua(num) - v0(nm)) - fac2 * (ua(nm) - v0(nm)) ) / guv(nm)
             !
          endif
          !
       enddo
       !
       if ( sbimp(m) .and. .not.linswe ) then    ! described water level at lower boundary point
          !
          nfm  = kgrpnt(m,nf )
          nfum = kgrpnt(m,nfu)
          !
          if ( .not. qym(nfum) < 0. ) then
             !
             advecy(nfm) = 0.
             !
          else
             !
             advecy(nfm) = qym(nfum) * (v0(nfum) - v0(nfm)) / ( guv(nfm) * hvm(nfm) )
             !
          endif
          !
       endif
       !
       if ( stimp(m) .and. .not.linswe ) then    ! described water level at upper boundary point
          !
          nlm  = kgrpnt(m,nl  )
          nldm = kgrpnt(m,nl-1)
          !
          if ( .not. qym(nlm) < 0. ) then
             !
             advecy(nlm) = qym(nlm) * (v0(nlm) - v0(nldm)) / ( guv(nlm) * hvm(nlm) )
             !
          else
             !
             advecy(nlm) = 0.
             !
          endif
          !
       endif
       !
    enddo
    !
    ! compute the curvature terms associated with advection terms (internal points only!)
    !
    if ( optg == 3 .and. .not.linswe ) then
       !
       do m = mfu, ml
          do n = nsta, nend
             !
             nd = n - 1
             nu = n + 1
             md = m - 1
             !
             nm   = kgrpnt(m ,n )
             ndm  = kgrpnt(m ,nd)
             num  = kgrpnt(m ,nu)
             nmd  = kgrpnt(md,n )
             numd = kgrpnt(md,nu)
             !
             if ( wetv(nm) == 1 ) then
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( ndm  == 1 ) ndm  = nm
                if ( num  == 1 ) num  = nm
                if ( nmd  == 1 ) nmd  = nm
                if ( numd == 1 ) numd = nmd
                !
                u = 0.25 * ( u0(nm) + u0(num) + u0(nmd) + u0(numd) )
                !
                if ( advecx(nm) /= 0. ) advecx(nm) = advecx(nm) - 0.5 * u * u * ( gvv(num) - gvv(ndm) ) / gsqsv(nm)
                !
                if ( advecy(nm) /= 0. ) advecy(nm) = advecy(nm) + 0.5 * v0(nm) * u * ( guu(nm) + guu(num) - guu(nmd) - guu(numd) ) / gsqsv(nm)
                !
             endif
             !
          enddo
       enddo
       !
    endif
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
             if ( wetv(nm) * wetv(num) == 1 .and. (n /= nf .or. sbimp(m) .or. .not.LMYF) .and. (n /= nl .or. stimp(m) .or. .not.LMYL) ) then
                !
                fac = 1. / ( 2. * guv(nm) * hvm(nm) )
                !
                gmatv(nm,1) = -fac * ( s0(nm ) + dps(num) )
                gmatv(nm,2) =  fac * ( s0(num) + dps(nm ) )
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
       if ( iproj == 2 ) qgrad = (1.-theta) * qgrad
       !
    endif
    !
    ! compute viscosity term
    !
    if ( ihvisc == 1 .and. hvisc > 0. ) then
       !
       do m = mfu, ml
          !
          nfm = kgrpnt(m,nf)
          nlm = kgrpnt(m,nl)
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
             if ( wetv(nm) == 1 ) then
                !
                if ( wetv(ndm) * wetv(num) == 1 ) then
                   !
                   fac1 = hs(num) * ( gvv(num)*v0(num) - gvv(nm)*v0(nm) )/gsqs(num) - hs(nm) * ( gvv(nm)*v0(nm) - gvv(ndm)*v0(ndm) )/gsqs(nm)
                   !
                else
                   !
                   fac1 = 0.
                   !
                endif
                !
                if ( wetv(nmd) * wetv(nmu) == 1 ) then
                   !
                   gsqsd  = .25 * (gvv(nm ) + gvv(nmu)) * (guu(nm ) + guu(num ))
                   gsqsdl = .25 * (gvv(nmd) + gvv(nm )) * (guu(nmd) + guu(numd))
                   !
                   hd  = .5 * (hvm(nmu) + hvm(nm ))
                   hdl = .5 * (hvm(nm ) + hvm(nmd))
                   !
                   fac2 = hd * ( guv(nmu)*v0(nmu) - guv(nm)*v0(nm) )/gsqsd - hdl * ( guv(nm)*v0(nm) - guv(nmd)*v0(nmd) )/gsqsdl
                   !
                else
                   !
                   fac2 = 0.
                   !
                endif
                !
                visc(nm) = hvisc * ( fac1/guv(nm) + fac2/gvv(nm) ) / hvm(nm)
                !
             endif
             !
          enddo
          !
          if ( sbimp(m) ) visc(nfm) = 0.
          !
          if ( stimp(m) ) visc(nlm) = 0.
          !
       enddo
       !
    else if ( ihvisc > 1 ) then
       !
       do m = mfu, ml
          !
          nfm = kgrpnt(m,nf)
          nlm = kgrpnt(m,nl)
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
             ndmd = kgrpnt(md,nd)
             numd = kgrpnt(md,nu)
             !
             if ( wetv(nm) == 1 ) then
                !
                lmask = .true.
                !
                hd  = .5 * (hvm(nmu) + hvm(nm ))
                hdl = .5 * (hvm(nm ) + hvm(nmd))
                !
                ! first, compute the components of the Reynolds stress tensor in different points of consideration
                !
                if ( wetv(ndm) * wetv(num) == 1 ) then
                   !
                   ryys  = ( vnu2d(nm ) + vnu2d(nmd ) + vnu2d(ndm) + vnu2d(ndmd) ) * ( v0(nm ) - v0(ndm) ) / ( guu(nm ) + guu(nmd ) )
                   !
                   ryysu = ( vnu2d(num) + vnu2d(numd) + vnu2d(nm ) + vnu2d(nmd ) ) * ( v0(num) - v0(nm ) ) / ( guu(num) + guu(numd) )
                   !
                else
                   !
                   ryys  = 0.
                   ryysu = 0.
                   !
                   lmask = .false.
                   !
                endif
                !
                if ( wetv(nmd) * wetv(nmu) * wetu(nm) * wetu(nmd) * wetu(num) * wetu(numd) == 1 ) then
                   !
                   rxyd  = 2. * vnu2d(nm ) * ( ( u0(num ) - u0(nm ) )/( guu(num ) + guu(nm ) ) + ( v0(nmu) - v0(nm ) )/( gvv(nmu) + gvv(nm ) ) )
                   !
                   rxydd = 2. * vnu2d(nmd) * ( ( u0(numd) - u0(nmd) )/( guu(numd) + guu(nmd) ) + ( v0(nm ) - v0(nmd) )/( gvv(nm ) + gvv(nmd) ) )
                   !
                else
                   !
                   rxyd  = 0.
                   rxydd = 0.
                   !
                   lmask = .false.
                   !
                endif
                !
                if ( lmask ) then
                   !
                   rxxu  = 0.50 * ( vnu2d(nm) + vnu2d(nmd) ) * ( u0(num) + u0(nm) - u0(nmd) -u0(numd) ) / gvv(nm)
                   !
                   rxyu  = 0.25 * ( vnu2d(nm) + vnu2d(nmd) ) * ( ( v0(nmu) - v0(nmd) )/gvv(nm) + ( u0(num) + u0(numd) - u0(nm) - u0(nmd) )/guv(nm) )
                   !
                   ryyu  = 0.50 * ( vnu2d(nm) + vnu2d(nmd) ) * ( v0(num) - v0(ndm) ) / guv(nm)
                   !
                else
                   !
                   rxxu = 0.
                   rxyu = 0.
                   ryyu = 0.
                   !
                endif
                !
                ! next, compute the divergence of the Reynolds stresses
                !
                visc(nm) = ( hd*rxyd - hdl*rxydd ) / (hvm(nm)*gvv(nm)) + ( hs(num)*ryysu - hs(nm)*ryys ) / (hvm(nm)*guv(nm))
                !
                ! finally, add the curvature terms associated with those stresses
                !
                visc(nm) = visc(nm) + ( rxyu * ( guu(nm) + guu(num) - guu(nmd) - guu(numd) ) + 0.5 * (ryyu - rxxu) * ( gvv(num) - gvv(ndm) ) ) / gsqsv(nm)
                !
             endif
             !
          enddo
          !
          if ( sbimp(m) ) visc(nfm) = 0.
          !
          if ( stimp(m) ) visc(nlm) = 0.
          !
       enddo
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
    ! compute Coriolis force
    !
    if ( coriolis ) then
       !
       epsab2 = pcor(2)
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
             if ( wetv(nm) == 1 .and. fcor(nm,2) /= 0. ) then
                !
                ! at time level n-1
                !
                ! --- modified quasi second order Adams-Bashforth scheme which is conditionally stable for the inviscid case (see Marshall et al., 1997)
                !
                corf(nm) = -(0.5+epsab2) * cfv(nm,1)
                !
                ! at time level n
                !
                ! --- interpolate u-velocity at v-point using energy-conserved weighting factors according to Espelid et al. (2000)
                !
                u = sqrt(abs(fcor(nm,1))*dpu(nm))*u0(nm) + sqrt(abs(fcor(num,1))*dpu(num))*u0(num) + sqrt(abs(fcor(nmd,1))*dpu(nmd))*u0(nmd) + sqrt(abs(fcor(numd,1))*dpu(numd))*u0(numd)
                u = u / (4.*sqrt(abs(fcor(nm,2))*dpv(nm)))
                !
                cfv(nm,1) = fcor(nm,2) * u
                !
                corf(nm) = corf(nm) + (1.5+epsab2) * cfv(nm,1)
                !
             endif
             !
          enddo
       enddo
       !
    else
       !
       corf = 0.
       !
    endif
    !
    ! include effect of ambient current - first order approximation
    !
    if ( iamb == 1 ) then
       !
       ! compute the advective v-velocity in depth point based on first order upwind
       !
       do m = mf, ml
          do n = nf+1, nend
             !
             mu = m + 1
             nu = n + 1
             !
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             num  = kgrpnt(m ,nu)
             numu = kgrpnt(mu,nu)
             !
             fac = 0.25 * ( Ubck(nm,1) + Ubck(num,1) + Ubck(nmu,1) + Ubck(numu,1) )
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm) = v0(nm )
                !
             else
                !
                ua(nm) = v0(nmu)
                !
             endif
             !
          enddo
       enddo
       !
       ! compute advection term in x-direction (momentum conservative) at internal v-point (prediction)
       !
       do m = mfu, ml
          !
          do n = nf+1, nend
             !
             md = m - 1
             mu = m + 1
             nu = n + 1
             !
             nm   = kgrpnt(m ,n )
             nmd  = kgrpnt(md,n )
             nmu  = kgrpnt(mu,n )
             num  = kgrpnt(m ,nu)
             numd = kgrpnt(md,nu)
             numu = kgrpnt(mu,nu)
             !
             fac1 = 0.25 * ( Ubck(nm ,1) + Ubck(num ,1) + Ubck(nmu,1) + Ubck(numu,1) )
             fac2 = 0.25 * ( Ubck(nmd,1) + Ubck(numd,1) + Ubck(nm ,1) + Ubck(num ,1) )
             !
             ! to avoid boundary layer effect, set advection near lateral open boundary to zero
             !
             if ( ( .not.lreptx .and. LMXF .and. m == mfu ) .or. ( .not.lreptx .and. LMXL .and. m == ml ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetv(nm) == 1 ) then
                !
                advecx(nm) = advecx(nm) + ( fac1 * ( ua(nm) - v0(nm) ) - fac2 * ( ua(nmd) -v0(nm) ) ) / gvv(nm)
                !
             endif
             !
             ! also add u dV/dx
             !
             u = 0.25 * ( u0(nm) + u0(num) + u0(nmd) + u0(numd) )
             if ( ( .not.lreptx .and. LMXF .and. m == mfu ) .or. ( .not.lreptx .and. LMXL .and. m == ml ) ) u = 0.
             if ( wetv(nm) == 1 ) advecx(nm) = advecx(nm) + 0.25 * u * ( Vbck(numu,1) + Vbck(nmu,1) - Vbck(numd,1) - Vbck(nmd,1) ) / gvv(nm)
             !
          enddo
          !
       enddo
       !
       ! compute the advective v-velocity in wl-point based on first order upwind
       !
       do m = mfu, ml
          do n = nfu, nl
             !
             nd = n - 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             !
             if ( .not. Vbck(nm,1) < 0. ) then
                !
                ua(nm) = v0(ndm)
                !
             else
                !
                ua(nm) = v0(nm )
                !
             endif
             !
          enddo
       enddo
       !
       ! synchronize advective v-velocities at appropriate boundaries in case of repeating grid (if appropriate)
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       !
       ! compute advection term in y-direction (momentum/energy head conservative) at internal v-point (prediction)
       !
       do m = mfu, ml
          !
          do n = nf+1, nend
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( wetv(nm) == 1 ) then
                !
                advecy(nm) = advecy(nm) + ( Vbck(num,1) * ua(num) - Vbck(nm,1) * ua(nm) ) / guv(nm)
                !
             endif
             !
          enddo
          !
       enddo
       !
    else if ( iamb == 2 ) then
       !
       ! compute the advective v-velocity in depth point based on first order upwind
       !
       do m = mf, ml
          do n = nf+1, nend
             !
             mu = m + 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmu = kgrpnt(mu,n )
             num = kgrpnt(m ,nu)
             !
             fac = 0.5 * ( Ubck(nm,1) + Ubck(num,1) )
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm) = v0(nm )
                !
             else
                !
                ua(nm) = v0(nmu)
                !
             endif
             !
          enddo
       enddo
       !
       ! compute advection term in x-direction (momentum conservative) at internal v-point (prediction)
       !
       do m = mfu, ml
          !
          do n = nf+1, nend
             !
             md = m - 1
             mu = m + 1
             nu = n + 1
             !
             nm   = kgrpnt(m ,n )
             nmd  = kgrpnt(md,n )
             nmu  = kgrpnt(mu,n )
             num  = kgrpnt(m ,nu)
             numd = kgrpnt(md,nu)
             !
             fac1 = 0.5 * ( Ubck(nm ,1) + Ubck(num ,1) )
             fac2 = 0.5 * ( Ubck(nmd,1) + Ubck(numd,1) )
             !
             ! to avoid boundary layer effect, set advection near lateral open boundary to zero
             !
             if ( ( .not.lreptx .and. LMXF .and. m == mfu ) .or. ( .not.lreptx .and. LMXL .and. m == ml ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetv(nm) == 1 ) then
                !
                advecx(nm) = advecx(nm) + ( fac1 * ( ua(nm) - v0(nm) ) - fac2 * ( ua(nmd) -v0(nm) ) ) / gvv(nm)
                !
             endif
             !
             ! also add u dV/dx
             !
             u = 0.25 * ( u0(nm) + u0(num) + u0(nmd) + u0(numd) )
             if ( ( .not.lreptx .and. LMXF .and. m == mfu ) .or. ( .not.lreptx .and. LMXL .and. m == ml ) ) u = 0.
             if ( wetv(nm) == 1 ) advecx(nm) = advecx(nm) + 0.5 * u * ( Vbck(nmu,1) - Vbck(nmd,1) ) / gvv(nm)
             !
          enddo
          !
       enddo
       !
       ! compute the advective v-velocity in wl-point based on first order upwind
       !
       do m = mfu, ml
          do n = nfu, nl
             !
             nd = n - 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             !
             fac = 0.5 * ( Vbck(nm,1) + Vbck(ndm,1) )
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm) = v0(ndm)
                !
             else
                !
                ua(nm) = v0(nm )
                !
             endif
             !
          enddo
       enddo
       !
       ! synchronize advective v-velocities at appropriate boundaries in case of repeating grid (if appropriate)
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       !
       ! compute advection term in y-direction (momentum/energy head conservative) at internal v-point (prediction)
       !
       do m = mfu, ml
          !
          do n = nf+1, nend
             !
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             num = kgrpnt(m,nu)
             !
             fac1 = 0.5 * ( Vbck(num,1) + Vbck(nm ,1) )
             fac2 = 0.5 * ( Vbck(nm ,1) + Vbck(ndm,1) )
             !
             if ( wetv(nm) == 1 ) then
                !
                advecy(nm) = advecy(nm) + ( fac1 * ua(num) - fac2 * ua(nm) ) / guv(nm)
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! compute intermediate v-velocity
    !
    do m = mfu, ml
       !
       nfm = kgrpnt(m,nf)
       nlm = kgrpnt(m,nl)
       !
       if ( sbimp(m) ) then       ! described water level at lower boundary point
          !
          nsta = nf
          !
          ! neglect wind, pressure gradients and Coriolis force at water level boundary
          !
          windv (nfm) = 0.
          wndimp(nfm) = 0.
          pgrad (nfm) = 0.
          corf  (nfm) = 0.
          !
       else
          !
          nsta = nf + 1
          !
       endif
       !
       if ( stimp(m) ) then       ! described water level at upper boundary point
          !
          nend = nl
          !
          ! neglect wind, pressure gradients and Coriolis force at water level boundary
          !
          windv (nlm) = 0.
          wndimp(nlm) = 0.
          pgrad (nlm) = 0.
          corf  (nlm) = 0.
          !
       else if ( lrepty ) then    ! periodicity at upper boundary point
          !
          nend = nl
          !
       else
          !
          nend = nl - 1
          !
       endif
       !
       do n = nsta, nend
          !
          nu = n + 1
          md = m - 1
          !
          nm   = kgrpnt(m ,n )
          num  = kgrpnt(m ,nu)
          nmd  = kgrpnt(md,n )
          numd = kgrpnt(md,nu)
          !
          if ( wetv(nm) == 1 .and. num /= 1 ) then
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
             contrib = advecx(nm) + advecy(nm) + zgrad + pgrad(nm) + qgrad(nm) - windv(nm)/max(1.e-3,hvm(nm)) - visc(nm) + corf(nm)
             !
             denom = 1. + cvegv(nm,1,2) + cpomv(nm,1) + dt * ( cbot(nm) + cporfr + cveg + wndimp(nm) )
             !
             ! compute flow velocity
             !
             v1(nm) = ( (1.+cvegv(nm,1,2)+cpomv(nm,1)) * v0(nm) - dt * contrib ) / denom
             !
          else
             !
             v1(nm) = 0.
             !
          endif
          !
       enddo
       !
    enddo
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
    ! compute higher order correction to the advective velocity in depth point (if appropriate)
    !
    propsc = nint(pnums(6))
    !
    if ( propsc /= 1 ) then
       !
       nsta = nf + 1            ! first internal v-point
       if ( lrepty ) then
          nend = nl             ! last  internal v-point in case of repeating grid
       else
          nend = nl - 1         ! last  internal v-point
       endif
       !
       do m = mf, ml
          do n = nsta, nend
             !
             md  = m  - 1
             mu  = m  + 1
             muu = mu + 1
             nd  = n  - 1
             nu  = n  + 1
             !
             if ( LMXF .and. md  < mfu ) md  = mfu
             if ( LMXL .and. mu  > ml  ) mu  = ml
             if ( LMXL .and. muu > ml  ) muu = ml
             !
             nm   = kgrpnt(m  ,n )
             nmd  = kgrpnt(md ,n )
             nmu  = kgrpnt(mu ,n )
             nmuu = kgrpnt(muu,n )
             ndm  = kgrpnt(m  ,nd)
             num  = kgrpnt(m  ,nu)
             !
             if ( ( v0(nm) > v0(ndm) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac = 0.5 * ( u0(nm) + u0(num) )
                !
             else
                !
                fac = qxm(nm)
                !
             endif
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = v1(nmu) - v0(nm )
                grad2 = v1(nm ) - v0(nmd)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = v1(nm ) - v0(nmu )
                grad2 = v1(nmu) - v0(nmuu)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute higher order correction to the advective velocity in wl-point (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do m = mfu, ml
          do n = nfu, nl
             !
             nd  = n  - 1
             ndd = nd - 1
             nu  = n  + 1
             !
             if ( LMYF .and. ndd < nf ) ndd = nf
             if ( LMYL .and. nu  > nl ) nu  = nl
             !
             nm   = kgrpnt(m,n  )
             ndm  = kgrpnt(m,nd )
             nddm = kgrpnt(m,ndd)
             num  = kgrpnt(m,nu )
             !
             if ( ( v0(nm) > v0(ndm) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac = 0.5 * ( v0(nm) + v0(ndm) )
                !
             else
                !
                fac = qym(nm)
                !
             endif
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = v1(nm ) - v0(ndm )
                grad2 = v1(ndm) - v0(nddm)
                !
                up(nm) = v0(ndm)
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = v1(ndm) - v0(nm )
                grad2 = v1(nm ) - v0(num)
                !
                up(nm) = v0(nm)
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       call periodic ( up, kgrpnt, 1, 1 )
       !
    endif
    !
    ! correct v-velocity associated with higher order approximation of advection term in x-direction (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do m = mfu, ml
          do n = nsta, nend
             !
             nd = n - 1
             nu = n + 1
             md = m - 1
             !
             nm  = kgrpnt(m ,n )
             ndm = kgrpnt(m ,nd)
             num = kgrpnt(m ,nu)
             nmd = kgrpnt(md,n )
             !
             if ( ( v0(nm) > v0(ndm) .or. stricthead ) .and. .not.strictmom ) then
                !
                !fac1 = 0.5 * ( u0(nm) + u0(num) )
                !fac2 = fac1
                !jumps back to old implementation (020117)
                fac1 = 0.
                fac2 = 0.
                !
             else if ( hvmn(nm) > 0. ) then
                !
                fac1 = qxm(nm ) / hvmn(nm)
                fac2 = qxm(nmd) / hvmn(nm)
                !
             else
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             ! set higher order approximation near lateral boundary to zero
             !
             if ( ( LMXF .and. m == mfu ) .or. ( LMXL .and. m == ml ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetv(nm) == 1 .and. num /= 1 ) then
                !
                v1(nm) = v1(nm) - dt * ( fac1 * ua2(nm) - fac2 * ua2(nmd) ) / gvv(nm)
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! correct v-velocity associated with higher order approximation of advection term in y-direction (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do m = mfu, ml
          do n = nsta, nend
             !
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             num = kgrpnt(m,nu)
             !
             if ( ( v0(nm) > v0(ndm) .or. stricthead ) .and. .not.strictmom ) then
                !
                !fac1 = up(num) + 0.5*ua(num)
                !fac2 = up(nm ) + 0.5*ua(nm )
                !jumps back to old implementation (020117)
                fac1 = 0.
                fac2 = 0.
                !
             else if ( hvmn(nm) > 0. ) then
                !
                fac1 = qym(num) / hvmn(nm)
                fac2 = qym(nm ) / hvmn(nm)
                !
             else
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetv(nm) == 1 .and. num /= 1 ) then
                !
                v1(nm) = v1(nm) - dt * ( fac1 * ua(num)  - fac2 * ua(nm) ) / guv(nm)
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! exchange v-velocities with neighbouring subdomains (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       call SWEXCHG ( v1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    endif
    !
    ! synchronize v-velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    if ( propsc /= 1 ) call periodic ( v1, kgrpnt, 1, 1 )
    !
    ! include effect of ambient current - high order correction
    !
    propsc = nint(pnums(61))
    !
    if ( iamb == 1 .and. propsc /= 1 ) then
       !
       ! compute higher order correction to the advective velocity in depth point (if appropriate)
       !
       kappa  = pnums(62)
       mbound = pnums(63)
       phieby = pnums(64)
       !
       if ( lrepty ) then
          nend = nl             ! last internal v-point in case of repeating grid
       else
          nend = nl - 1         ! last internal v-point
       endif
       !
       do m = mf, ml
          do n = nf+1, nend
             !
             md  = m  - 1
             mu  = m  + 1
             muu = mu + 1
             nu  = n  + 1
             !
             if ( LMXF .and. md  < mfu ) md  = mfu
             if ( LMXL .and. mu  > ml  ) mu  = ml
             if ( LMXL .and. muu > ml  ) muu = ml
             !
             nm   = kgrpnt(m  ,n )
             nmd  = kgrpnt(md ,n )
             nmu  = kgrpnt(mu ,n )
             nmuu = kgrpnt(muu,n )
             num  = kgrpnt(m  ,nu)
             numu = kgrpnt(mu ,nu)
             !
             fac = 0.25 * ( Ubck(nm,1) + Ubck(num,1) + Ubck(nmu,1) + Ubck(numu,1) )
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = v1(nmu) - v0(nm )
                grad2 = v1(nm ) - v0(nmd)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = v1(nm ) - v0(nmu )
                grad2 = v1(nmu) - v0(nmuu)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       ! compute higher order correction to the advective velocity in wl-point (if appropriate)
       !
       do m = mfu, ml
          do n = nfu, nl
             !
             nd  = n  - 1
             ndd = nd - 1
             nu  = n  + 1
             !
             if ( LMYF .and. ndd < nf ) ndd = nf
             if ( LMYL .and. nu  > nl ) nu  = nl
             !
             nm   = kgrpnt(m,n  )
             ndm  = kgrpnt(m,nd )
             nddm = kgrpnt(m,ndd)
             num  = kgrpnt(m,nu )
             !
             if ( .not. Vbck(nm,1) < 0. ) then
                !
                grad1 = v1(nm ) - v0(ndm )
                grad2 = v1(ndm) - v0(nddm)
                !
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = v1(ndm) - v0(nm )
                grad2 = v1(nm ) - v0(num)
                !
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       !
       ! correct v-velocity associated with higher order approximation of advection term in x-direction (if appropriate)
       !
       do m = mfu, ml
          do n = nf+1, nend
             !
             nu = n + 1
             md = m - 1
             mu = m + 1
             !
             nm   = kgrpnt(m ,n )
             num  = kgrpnt(m ,nu)
             nmd  = kgrpnt(md,n )
             nmu  = kgrpnt(mu,n )
             numd = kgrpnt(md,nu)
             numu = kgrpnt(mu,nu)
             !
             fac1 = 0.25 * ( Ubck(nm ,1) + Ubck(num ,1) + Ubck(nmu,1) + Ubck(numu,1) )
             fac2 = 0.25 * ( Ubck(nmd,1) + Ubck(numd,1) + Ubck(nm ,1) + Ubck(num ,1) )
             !
             ! set higher order approximation near lateral boundary to zero
             !
             if ( ( LMXF .and. m == mfu ) .or. ( LMXL .and. m == ml ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetv(nm) == 1 ) then
                !
                v1(nm) = v1(nm) - dt * ( fac1 * ua2(nm) - fac2 * ua2(nmd) ) / gvv(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! correct v-velocity associated with higher order approximation of advection term in y-direction (if appropriate)
       !
       do m = mfu, ml
          do n = nf+1, nend
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( wetv(nm) == 1 ) then
                !
                v1(nm) = v1(nm) - dt * ( Vbck(num,1) * ua(num)  - Vbck(nm,1) * ua(nm) ) / guv(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! exchange v-velocities with neighbouring subdomains (if appropriate)
       !
       call SWEXCHG ( v1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize v-velocities at appropriate boundaries in case of repeating grid (if appropriate)
       !
       call periodic ( v1, kgrpnt, 1, 1 )
       !
    else if ( iamb == 2 .and. propsc /= 1 ) then
       !
       ! compute higher order correction to the advective velocity in depth point (if appropriate)
       !
       kappa  = pnums(62)
       mbound = pnums(63)
       phieby = pnums(64)
       !
       if ( lrepty ) then
          nend = nl             ! last internal v-point in case of repeating grid
       else
          nend = nl - 1         ! last internal v-point
       endif
       !
       do m = mf, ml
          do n = nf+1, nend
             !
             md  = m  - 1
             mu  = m  + 1
             muu = mu + 1
             nu  = n  + 1
             !
             if ( LMXF .and. md  < mfu ) md  = mfu
             if ( LMXL .and. mu  > ml  ) mu  = ml
             if ( LMXL .and. muu > ml  ) muu = ml
             !
             nm   = kgrpnt(m  ,n )
             nmd  = kgrpnt(md ,n )
             nmu  = kgrpnt(mu ,n )
             nmuu = kgrpnt(muu,n )
             num  = kgrpnt(m  ,nu)
             !
             fac = 0.5 * ( Ubck(nm,1) + Ubck(num,1) )
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = v1(nmu) - v0(nm )
                grad2 = v1(nm ) - v0(nmd)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = v1(nm ) - v0(nmu )
                grad2 = v1(nmu) - v0(nmuu)
                !
                ua2(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       ! compute higher order correction to the advective velocity in wl-point (if appropriate)
       !
       do m = mfu, ml
          do n = nfu, nl
             !
             nd  = n  - 1
             ndd = nd - 1
             nu  = n  + 1
             !
             if ( LMYF .and. ndd < nf ) ndd = nf
             if ( LMYL .and. nu  > nl ) nu  = nl
             !
             nm   = kgrpnt(m,n  )
             ndm  = kgrpnt(m,nd )
             nddm = kgrpnt(m,ndd)
             num  = kgrpnt(m,nu )
             !
             fac = 0.5 * ( Vbck(nm,1) + Vbck(ndm,1) )
             !
             if ( .not. fac < 0. ) then
                !
                grad1 = v1(nm ) - v0(ndm )
                grad2 = v1(ndm) - v0(nddm)
                !
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             else
                !
                grad1 = v1(ndm) - v0(nm )
                grad2 = v1(nm ) - v0(num)
                !
                ua(nm) = 0.5 * fluxlim(grad1,grad2)
                !
             endif
             !
          enddo
       enddo
       !
       call periodic ( ua, kgrpnt, 1, 1 )
       !
       ! correct v-velocity associated with higher order approximation of advection term in x-direction (if appropriate)
       !
       do m = mfu, ml
          do n = nf+1, nend
             !
             nu = n + 1
             md = m - 1
             !
             nm   = kgrpnt(m ,n )
             num  = kgrpnt(m ,nu)
             nmd  = kgrpnt(md,n )
             numd = kgrpnt(md,nu)
             !
             fac1 = 0.5 * ( Ubck(nm ,1) + Ubck(num ,1) )
             fac2 = 0.5 * ( Ubck(nmd,1) + Ubck(numd,1) )
             !
             ! set higher order approximation near lateral boundary to zero
             !
             if ( ( LMXF .and. m == mfu ) .or. ( LMXL .and. m == ml ) ) then
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( wetv(nm) == 1 ) then
                !
                v1(nm) = v1(nm) - dt * ( fac1 * ua2(nm) - fac2 * ua2(nmd) ) / gvv(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! correct v-velocity associated with higher order approximation of advection term in y-direction (if appropriate)
       !
       do m = mfu, ml
          do n = nf+1, nend
             !
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             num = kgrpnt(m,nu)
             !
             fac1 = 0.5 * ( Vbck(num,1) + Vbck(nm ,1) )
             fac2 = 0.5 * ( Vbck(nm ,1) + Vbck(ndm,1) )
             !
             if ( wetv(nm) == 1 ) then
                !
                v1(nm) = v1(nm) - dt * ( fac1 * ua(num)  - fac2 * ua(nm) ) / guv(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! exchange v-velocities with neighbouring subdomains (if appropriate)
       !
       call SWEXCHG ( v1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize v-velocities at appropriate boundaries in case of repeating grid (if appropriate)
       !
       call periodic ( v1, kgrpnt, 1, 1 )
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
          fac = 1. - theta
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
                w1top(nm) = w0top(nm) + w0bot(nm) - w1bot(nm) + 2.*dt*fac*q(nm)/hs(nm)
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
          propsc = nint(pnums(16))
          kappa  = pnums(17)
          mbound = pnums(18)
          phieby = pnums(19)
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
                if ( wets(nm) == 1 ) then
                   !
                   ! advection term in x-direction
                   !
                   fac = 2. * dt / gsqs(nm)
                   !
                   fac1 = fac * guu(nmd) * qx(nmd) / hs(nm)
                   fac2 = fac * guu(nm ) * qx(nm ) / hs(nm)
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      w1top(nm) = w1top(nm) - 0.5 * ( fac1 * (w0top(nm) - w0top(nmd)) + fac2 * (w0top(nmu) - w0top(nm)) )
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      w1top(nm) = w1top(nm) - ( max(fac1,0.) * (w0top(nm) - w0top(nmd)) + min(fac2,0.) * (w0top(nmu) - w0top(nm)) )
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         if ( fac2 > 0. ) then
                            !
                            grad1 = w0top(nmu) - w0top(nm )
                            grad2 = w0top(nm ) - w0top(nmd)
                            !
                            w1top(nm) = w1top(nm) - 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         else if ( fac2 < 0. ) then
                            !
                            grad1 = w0top(nmu ) - w0top(nm )
                            grad2 = w0top(nmuu) - w0top(nmu)
                            !
                            w1top(nm) = w1top(nm) + 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                         if ( fac1 > 0. ) then
                            !
                            grad1 = w0top(nm ) - w0top(nmd )
                            grad2 = w0top(nmd) - w0top(nmdd)
                            !
                            w1top(nm) = w1top(nm) + 0.5 * fac1 * fluxlim(grad1,grad2)
                            !
                         else if ( fac1 < 0. ) then
                            !
                            grad1 = w0top(nm ) - w0top(nmd)
                            grad2 = w0top(nmu) - w0top(nm )
                            !
                            w1top(nm) = w1top(nm) - 0.5 * fac1 * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                      endif
                      !
                   endif
                   !
                   ! advection term in y-direction
                   !
                   fac1 = fac * gvv(ndm) * qy(ndm) / hs(nm)
                   fac2 = fac * gvv(nm ) * qy(nm ) / hs(nm)
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      w1top(nm) = w1top(nm) - 0.5 * ( fac1 * (w0top(nm) - w0top(ndm)) + fac2 * (w0top(num) - w0top(nm)) )
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      w1top(nm) = w1top(nm) - ( max(fac1,0.) * (w0top(nm) - w0top(ndm)) + min(fac2,0.) * (w0top(num) - w0top(nm)) )
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         if ( fac2 > 0. ) then
                            !
                            grad1 = w0top(num) - w0top(nm )
                            grad2 = w0top(nm ) - w0top(ndm)
                            !
                            w1top(nm) = w1top(nm) - 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         else if ( fac2 < 0. ) then
                            !
                            grad1 = w0top(num ) - w0top(nm )
                            grad2 = w0top(nuum) - w0top(num)
                            !
                            w1top(nm) = w1top(nm) + 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                         if ( fac1 > 0. ) then
                            !
                            grad1 = w0top(nm ) - w0top(ndm )
                            grad2 = w0top(ndm) - w0top(nddm)
                            !
                            w1top(nm) = w1top(nm) + 0.5 * fac1 * fluxlim(grad1,grad2)
                            !
                         else if ( fac1 < 0. ) then
                            !
                            grad1 = w0top(nm ) - w0top(ndm)
                            grad2 = w0top(num) - w0top(nm )
                            !
                            w1top(nm) = w1top(nm) - 0.5 * fac1 * fluxlim(grad1,grad2)
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
          enddo
          !
       endif
       !
       if ( iamb /= 0 ) then
          !
          ! compute horizontal advection terms (momentum conservative and explicit)
          !
          propsc = nint(pnums(66))
          kappa  = pnums(67)
          mbound = pnums(68)
          phieby = pnums(69)
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
                if ( wets(nm) == 1 ) then
                   !
                   ! advection term in x-direction
                   !
                   fac = 2. * dt / gsqs(nm)
                   !
                   if ( iamb == 1 ) then
                      fac1 = 0.5 * ( Ubck(nm ,1) + Ubck(nmd,1) )
                      fac2 = 0.5 * ( Ubck(nmu,1) + Ubck(nm ,1) )
                   else if ( iamb == 2 ) then
                      fac1 = Ubck(nmd,1)
                      fac2 = Ubck(nm ,1)
                   endif
                   fac1 = fac * guu(nmd) * hu(nmd) * fac1 / hs(nm)
                   fac2 = fac * guu(nm ) * hu(nm ) * fac2 / hs(nm)
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      w1top(nm) = w1top(nm) - 0.5 * ( fac1 * (w0top(nm) - w0top(nmd)) + fac2 * (w0top(nmu) - w0top(nm)) )
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      w1top(nm) = w1top(nm) - ( max(fac1,0.) * (w0top(nm) - w0top(nmd)) + min(fac2,0.) * (w0top(nmu) - w0top(nm)) )
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         if ( fac2 > 0. ) then
                            !
                            grad1 = w0top(nmu) - w0top(nm )
                            grad2 = w0top(nm ) - w0top(nmd)
                            !
                            w1top(nm) = w1top(nm) - 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         else if ( fac2 < 0. ) then
                            !
                            grad1 = w0top(nmu ) - w0top(nm )
                            grad2 = w0top(nmuu) - w0top(nmu)
                            !
                            w1top(nm) = w1top(nm) + 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                         if ( fac1 > 0. ) then
                            !
                            grad1 = w0top(nm ) - w0top(nmd )
                            grad2 = w0top(nmd) - w0top(nmdd)
                            !
                            w1top(nm) = w1top(nm) + 0.5 * fac1 * fluxlim(grad1,grad2)
                            !
                         else if ( fac1 < 0. ) then
                            !
                            grad1 = w0top(nm ) - w0top(nmd)
                            grad2 = w0top(nmu) - w0top(nm )
                            !
                            w1top(nm) = w1top(nm) - 0.5 * fac1 * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                      endif
                      !
                   endif
                   !
                   ! advection term in y-direction
                   !
                   if ( iamb == 1 ) then
                      fac1 = 0.5 * ( Vbck(nm ,1) + Vbck(ndm,1) )
                      fac2 = 0.5 * ( Vbck(num,1) + Vbck(nm ,1) )
                   else if ( iamb == 2 ) then
                      fac1 = Vbck(ndm,1)
                      fac2 = Vbck(nm ,1)
                   endif
                   fac1 = fac * gvv(ndm) * hv(ndm) * fac1 / hs(nm)
                   fac2 = fac * gvv(nm ) * hv(nm ) * fac2 / hs(nm)
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      w1top(nm) = w1top(nm) - 0.5 * ( fac1 * (w0top(nm) - w0top(ndm)) + fac2 * (w0top(num) - w0top(nm)) )
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      w1top(nm) = w1top(nm) - ( max(fac1,0.) * (w0top(nm) - w0top(ndm)) + min(fac2,0.) * (w0top(num) - w0top(nm)) )
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         if ( fac2 > 0. ) then
                            !
                            grad1 = w0top(num) - w0top(nm )
                            grad2 = w0top(nm ) - w0top(ndm)
                            !
                            w1top(nm) = w1top(nm) - 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         else if ( fac2 < 0. ) then
                            !
                            grad1 = w0top(num ) - w0top(nm )
                            grad2 = w0top(nuum) - w0top(num)
                            !
                            w1top(nm) = w1top(nm) + 0.5 * fac2 * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                         if ( fac1 > 0. ) then
                            !
                            grad1 = w0top(nm ) - w0top(ndm )
                            grad2 = w0top(ndm) - w0top(nddm)
                            !
                            w1top(nm) = w1top(nm) + 0.5 * fac1 * fluxlim(grad1,grad2)
                            !
                         else if ( fac1 < 0. ) then
                            !
                            grad1 = w0top(nm ) - w0top(ndm)
                            grad2 = w0top(num) - w0top(nm )
                            !
                            w1top(nm) = w1top(nm) - 0.5 * fac1 * fluxlim(grad1,grad2)
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
                if ( wets(nm) == 1 ) then
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
             if ( wets(nm) == 1 ) then
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
                rhs(nm) = rhs(nm)/(dt*theta)
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
                u1(nm) = u1(nm) - dt*theta*( gmatu(nm,1) * dq(nm) + gmatu(nm,2) * dq(nmu) )
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
                v1(nm) = v1(nm) - dt*theta*( gmatv(nm,1) * dq(nm) + gmatv(nm,2) * dq(num) )
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
                w1top(nm) = w1top(nm) + w1bot(nm) + 2.*dt*theta*dq(nm)/hs(nm)
                w1bot(nm) = -0.5 * ( ( u1(nm) + u1(nmd) ) * ( dpu(nm) - dpu(nmd) ) / dxl + ( v1(nm) + v1(ndm) ) * ( dpv(nm) - dpv(ndm) ) / dyl )
                w1top(nm) = w1top(nm) - w1bot(nm)
                !
             endif
             !
          enddo
       enddo
       !
    endif
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
             ! set to zero at closed boundary
             !
             v1(nmf) = 0.
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
             ! set to zero at closed boundary
             !
             v1(nmlu) = 0.
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
             ! set to zero at closed boundary
             !
             u1(nfm) = 0.
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
             ! set to zero at closed boundary
             !
             u1(nlum) = 0.
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
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
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
       ! compute the mass flux in x-direction
       !
       if ( .not.linswe ) then
          !
          do n = nfu, nl
             do m = mf, ml
                !
                nm = kgrpnt(m,n)
                !
                qx(nm) = hu(nm)*guu(nm)*u1(nm)
                !
             enddo
          enddo
          !
       else
          !
          do n = nfu, nl
             do m = mf, ml
                !
                nm = kgrpnt(m,n)
                !
                qx(nm) = dpu(nm)*guu(nm)*u1(nm)
                !
             enddo
          enddo
          !
       endif
       !
       ! compute the mass flux in y-direction
       !
       if ( .not.linswe ) then
          !
          do n = nf, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                qy(nm) = hv(nm)*gvv(nm)*v1(nm)
                !
             enddo
          enddo
          !
       else
          !
          do n = nf, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                qy(nm) = dpv(nm)*gvv(nm)*v1(nm)
                !
             enddo
          enddo
          !
       endif
       !
       ! compute the water level
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
             s1(nm) = s0(nm) - fh * dt * ( qx(nm) - qx(nmd) + qy(nm) - qy(ndm) ) / gsqs(nm)
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       s1(1) = 0.
       !
       ! exchange water levels with neighbouring subdomains
       !
       call SWEXCHG ( s1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize water level at appropriate boundaries in case of repeating grid
       !
       call periodic ( s1, kgrpnt, 1, 1 )
       !
       ! copy to virtual cells at boundaries except for both water level and Riemann invariant openings
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
             if ( ibl(n) /= 2 .and. ibl(n) /= 6 .and. LMXF ) s1(nmf ) = s1(nmfu)
             if ( ibr(n) /= 2 .and. ibr(n) /= 6 .and. LMXL ) s1(nmlu) = s1(nml )
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
             if ( ibb(m) /= 2 .and. ibb(m) /= 6 .and. LMYF ) s1(nfm ) = s1(nfum)
             if ( ibt(m) /= 2 .and. ibt(m) /= 6 .and. LMYL ) s1(nlum) = s1(nlm )
             !
          enddo
          !
       endif
       !
       ! in case of outflow condition, no change in water depth
       !
       do n = nfu, nl
          !
          nmf   = kgrpnt(mf   ,n)
          nmfu  = kgrpnt(mfu  ,n)
          nmfuu = kgrpnt(mfu+1,n)
          nml   = kgrpnt(ml   ,n)
          nmld  = kgrpnt(ml-1 ,n)
          nmlu  = kgrpnt(mlu  ,n)
          !
          if ( ibl(n) == 10 ) s1(nmf ) = s1(nmfu) - dps(nmfu) + dps(nmfuu)
          if ( ibr(n) == 10 ) s1(nmlu) = s1(nml ) - dps(nml ) + dps(nmld )
          !
       enddo
       !
       do m = mfu, ml
          !
          nfm   = kgrpnt(m,nf   )
          nfum  = kgrpnt(m,nfu  )
          nfuum = kgrpnt(m,nfu+1)
          nlm   = kgrpnt(m,nl   )
          nldm  = kgrpnt(m,nl-1 )
          nlum  = kgrpnt(m,nlu  )
          !
          if ( ibb(m) == 10 ) s1(nfm ) = s1(nfum) - dps(nfum) + dps(nfuum)
          if ( ibt(m) == 10 ) s1(nlum) = s1(nlm ) - dps(nlm ) + dps(nldm )
          !
       enddo
       !
       ! correct the water depths for mass fluxes
       !
       if ( j <= 2 ) call SwashUpdDepu ( u1, v1 )
       !
    enddo
    !
    ! add mass source due to internal wave generation
    !
    if ( iwvgen /= 0 ) then
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             s1(nm) = s1(nm) + dt * srcm(nm)
             !
          enddo
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
                   work(nm,2) = work(nm,2) + ( Ubck(nmu,1) + Ubck(nm,1) ) * guu(nm) * ( so(nmu) + so(nm) ) - ( Ubck(nm,1) + Ubck(nmd,1) ) * guu(nmd) * ( so(nm) + so(nmd) )
                   work(nm,2) = work(nm,2) + ( Vbck(num,1) + Vbck(nm,1) ) * gvv(nm) * ( so(num) + so(nm) ) - ( Vbck(nm,1) + Vbck(ndm,1) ) * gvv(ndm) * ( so(nm) + so(ndm) )
                   work(nm,2) = 0.25 * work(nm,2)
                   !
                enddo
             enddo
             !
          else
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
                   fac1 = thetamb * s1(nmd) + (1.-thetamb) * so(nmd)
                   fac  = thetamb * s1(nm ) + (1.-thetamb) * so(nm )
                   fac2 = thetamb * s1(nmu) + (1.-thetamb) * so(nmu)
                   !
                   work(nm,2) = work(nm,2) + ( Ubck(nmu,1) + Ubck(nm,1) ) * guu(nm) * ( fac2 + fac ) - ( Ubck(nm,1) + Ubck(nmd,1) ) * guu(nmd) * ( fac + fac1 )
                   !
                   fac1 = thetamb * s1(ndm) + (1.-thetamb) * so(ndm)
                   fac2 = thetamb * s1(num) + (1.-thetamb) * so(num)
                   !
                   work(nm,2) = work(nm,2) + ( Vbck(num,1) + Vbck(nm,1) ) * gvv(nm) * ( fac2 + fac ) - ( Vbck(nm,1) + Vbck(ndm,1) ) * gvv(ndm) * ( fac + fac1 )
                   !
                   work(nm,2) = 0.25 * work(nm,2)
                   !
                enddo
             enddo
             !
          endif
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                s1(nm) = so(nm) - dt * work(nm,2) / gsqs(nm)
                !
             enddo
          enddo
          !
          ! set to zero for permanently dry points
          !
          s1(1) = 0.
          !
          ! exchange water levels with neighbouring subdomains
          !
          call SWEXCHG ( s1, kgrpnt, 1, 1 )
          if (STPNOW()) return
          !
          ! synchronize water level at appropriate boundaries in case of repeating grid
          !
          call periodic ( s1, kgrpnt, 1, 1 )
          !
          ! copy to virtual cells at boundaries except for both water level and Riemann invariant openings
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
                if ( ibl(n) /= 2 .and. ibl(n) /= 6 .and. LMXF ) s1(nmf ) = s1(nmfu)
                if ( ibr(n) /= 2 .and. ibr(n) /= 6 .and. LMXL ) s1(nmlu) = s1(nml )
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
                if ( ibb(m) /= 2 .and. ibb(m) /= 6 .and. LMYF ) s1(nfm ) = s1(nfum)
                if ( ibt(m) /= 2 .and. ibt(m) /= 6 .and. LMYL ) s1(nlum) = s1(nlm )
                !
             enddo
             !
          endif
          !
          ! in case of outflow condition, no change in water depth
          !
          do n = nfu, nl
             !
             nmf   = kgrpnt(mf   ,n)
             nmfu  = kgrpnt(mfu  ,n)
             nmfuu = kgrpnt(mfu+1,n)
             nml   = kgrpnt(ml   ,n)
             nmld  = kgrpnt(ml-1 ,n)
             nmlu  = kgrpnt(mlu  ,n)
             !
             if ( ibl(n) == 10 ) s1(nmf ) = s1(nmfu) - dps(nmfu) + dps(nmfuu)
             if ( ibr(n) == 10 ) s1(nmlu) = s1(nml ) - dps(nml ) + dps(nmld )
             !
          enddo
          !
          do m = mfu, ml
             !
             nfm   = kgrpnt(m,nf   )
             nfum  = kgrpnt(m,nfu  )
             nfuum = kgrpnt(m,nfu+1)
             nlm   = kgrpnt(m,nl   )
             nldm  = kgrpnt(m,nl-1 )
             nlum  = kgrpnt(m,nlu  )
             !
             if ( ibb(m) == 10 ) s1(nfm ) = s1(nfum) - dps(nfum) + dps(nfuum)
             if ( ibt(m) == 10 ) s1(nlum) = s1(nlm ) - dps(nlm ) + dps(nldm )
             !
          enddo
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
                   work(nm,2) = work(nm,2) + Ubck(nm,1) * guu(nm) * ( so(nmu) + so(nm) ) - Ubck(nmd,1) * guu(nmd) * ( so(nm) + so(nmd) )
                   work(nm,2) = work(nm,2) + Vbck(nm,1) * gvv(nm) * ( so(num) + so(nm) ) - Vbck(ndm,1) * gvv(ndm) * ( so(nm) + so(ndm) )
                   work(nm,2) = 0.5 * work(nm,2)
                   !
                enddo
             enddo
             !
          else
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
                   fac1 = thetamb * s1(nmd) + (1.-thetamb) * so(nmd)
                   fac  = thetamb * s1(nm ) + (1.-thetamb) * so(nm )
                   fac2 = thetamb * s1(nmu) + (1.-thetamb) * so(nmu)
                   !
                   work(nm,2) = work(nm,2) + Ubck(nm,1) * guu(nm) * ( fac2 + fac ) - Ubck(nmd,1) * guu(nmd) * ( fac + fac1 )
                   !
                   fac1 = thetamb * s1(ndm) + (1.-thetamb) * so(ndm)
                   fac2 = thetamb * s1(num) + (1.-thetamb) * so(num)
                   !
                   work(nm,2) = work(nm,2) + Vbck(nm,1) * gvv(nm) * ( fac2 + fac ) - Vbck(ndm,1) * gvv(ndm) * ( fac + fac1 )
                   !
                   work(nm,2) = 0.5 * work(nm,2)
                   !
                enddo
             enddo
             !
          endif
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                s1(nm) = so(nm) - dt * work(nm,2) / gsqs(nm)
                !
             enddo
          enddo
          !
          ! set to zero for permanently dry points
          !
          s1(1) = 0.
          !
          ! exchange water levels with neighbouring subdomains
          !
          call SWEXCHG ( s1, kgrpnt, 1, 1 )
          if (STPNOW()) return
          !
          ! synchronize water level at appropriate boundaries in case of repeating grid
          !
          call periodic ( s1, kgrpnt, 1, 1 )
          !
          ! copy to virtual cells at boundaries except for both water level and Riemann invariant openings
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
                if ( ibl(n) /= 2 .and. ibl(n) /= 6 .and. LMXF ) s1(nmf ) = s1(nmfu)
                if ( ibr(n) /= 2 .and. ibr(n) /= 6 .and. LMXL ) s1(nmlu) = s1(nml )
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
                if ( ibb(m) /= 2 .and. ibb(m) /= 6 .and. LMYF ) s1(nfm ) = s1(nfum)
                if ( ibt(m) /= 2 .and. ibt(m) /= 6 .and. LMYL ) s1(nlum) = s1(nlm )
                !
             enddo
             !
          endif
          !
          ! in case of outflow condition, no change in water depth
          !
          do n = nfu, nl
             !
             nmf   = kgrpnt(mf   ,n)
             nmfu  = kgrpnt(mfu  ,n)
             nmfuu = kgrpnt(mfu+1,n)
             nml   = kgrpnt(ml   ,n)
             nmld  = kgrpnt(ml-1 ,n)
             nmlu  = kgrpnt(mlu  ,n)
             !
             if ( ibl(n) == 10 ) s1(nmf ) = s1(nmfu) - dps(nmfu) + dps(nmfuu)
             if ( ibr(n) == 10 ) s1(nmlu) = s1(nml ) - dps(nml ) + dps(nmld )
             !
          enddo
          !
          do m = mfu, ml
             !
             nfm   = kgrpnt(m,nf   )
             nfum  = kgrpnt(m,nfu  )
             nfuum = kgrpnt(m,nfu+1)
             nlm   = kgrpnt(m,nl   )
             nldm  = kgrpnt(m,nl-1 )
             nlum  = kgrpnt(m,nlu  )
             !
             if ( ibb(m) == 10 ) s1(nfm ) = s1(nfum) - dps(nfum) + dps(nfuum)
             if ( ibt(m) == 10 ) s1(nlum) = s1(nlm ) - dps(nlm ) + dps(nldm )
             !
          enddo
          !
       enddo
       !
    endif
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
                vol = vol + s1(nm)
                !
                ! compute potential energy in the cell center
                ener = ener + 0.5 * grav * gsqs(nm) * s1(nm) * s1(nm)
                !
                ! compute kinetic energy in the cell center
                fac1 = ( s1(nm ) + dps(nm ) ) * gsqs(nm )
                fac2 = ( s1(nmd) + dps(nmd) ) * gsqs(nmd)
                ener = ener + 0.125 * ( fac1 + fac2 ) * u1(nmd) * u1(nmd)
                fac2 = ( s1(nmu) + dps(nmu) ) * gsqs(nmu)
                ener = ener + 0.125 * ( fac1 + fac2 ) * u1(nm ) * u1(nm )
                !
                fac2 = ( s1(ndm) + dps(ndm) ) * gsqs(ndm)
                ener = ener + 0.125 * ( fac1 + fac2 ) * v1(ndm) * v1(ndm)
                fac2 = ( s1(num) + dps(num) ) * gsqs(num)
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
             dep = s1(nm) + dps(nm)
             !
             u    = 0.5 * ( u1(nm) + u1(nmd) )
             v    = 0.5 * ( v1(nm) + v1(ndm) )
             utot = sqrt( u*u + v*v )
             !
             dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
             dyl = 0.5 * ( guu(nm) + guu(nmd) )
             !
             cfl = dt * ( sqrt(grav*dep) + utot ) * sqrt(1./(dxl*dxl) + 1./(dyl*dyl))
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
 101 format (2x,'the net mass outflow is ',e14.8e2)
 102 format (2x,'the total displaced volume of water is ',e14.8e2)
 103 format (2x,'the total energy is ',e14.8e2)
    !
end subroutine SwashExpDep2DHflow
