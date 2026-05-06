subroutine SwashImpLay2DHflow
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
!    4.05: Dirk Rijnsdorp
!    6.02: Tom Bogaard
!
!   Updates
!
!    1.00, December 2013: New subroutine
!    4.05,  January 2017: extension floating objects
!    6.02,     July 2019: extension 3D (non)linear k-eps model
!
!   Purpose
!
!   Performs the time integration for the non-hydrostatic, layer-averaged 2D shallow water equations
!
!   Method
!
!   The time integration with respect to the continuity equation and the water level gradient of the
!   momentum equations is based on a theta-scheme. Only a value of 0.5 <= theta <= 1 will be taken.
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
    use m_genarr, only: kgrpnt, guu, guv, gvu, gvv, gsqs, gsqsu, gsqsv, msta, mend, nsta, nend, sbimp, slimp, srimp, stimp, work
    use m_parall
    use SwashSolvedata, only: icond, iamout
    use SwashFlowdata
!
    implicit none
!
!   Parameter variables
!
    real   , parameter :: epswom = 0.001 ! tolerance for relative vertical velocity at surface
!
!   Local variables
!
    integer, save                      :: ient = 0 ! number of entries in this subroutine
    integer                            :: j        ! loop counter
    integer                            :: jj       ! iteration counter
    integer                            :: k        ! loop counter over vertical layers
    integer                            :: kd       ! index of layer k-1
    integer                            :: kdd      ! index of layer k-2
    integer                            :: ku       ! index of layer k+1
    integer                            :: kuu      ! index of layer k+2
    integer                            :: l        ! loop counter
    integer                            :: m        ! loop counter
    integer                            :: maxit    ! maximum number of iterations
    integer                            :: md       ! index of point m-1
    integer                            :: mdd      ! index of point m-2
    integer                            :: minf     ! m-index of point with largest error in solution
    integer                            :: mlast    ! index of last internal u-point
    integer                            :: mu       ! index of point m+1
    integer                            :: muu      ! index of point m+2
    integer                            :: n        ! loop counter
    integer                            :: nd       ! index of point n-1
    integer                            :: ndd      ! index of point n-2
    integer                            :: nddm     ! pointer to m,n-2
    integer                            :: ndm      ! pointer to m,n-1
    integer                            :: ndmd     ! pointer to m-1,n-1
    integer                            :: ndmu     ! pointer to m+1,n-1
    integer                            :: nfm      ! pointer to m,nf
    integer                            :: nfum     ! pointer to m,nfu
    integer                            :: nfuum    ! pointer to m,nfuu
    integer                            :: ninf     ! n-index of point with largest error in solution
    integer                            :: nlast    ! index of last internal v-point
    integer                            :: nldm     ! pointer to m,nld
    integer                            :: nlm      ! pointer to m,nl
    integer                            :: nlum     ! pointer to m,nlu
    integer                            :: nm       ! pointer to m,n
    integer                            :: nmd      ! pointer to m-1,n
    integer                            :: nmdd     ! pointer to m-2,n
    integer                            :: nmf      ! pointer to mf,n
    integer                            :: nmfu     ! pointer to mfu,n
    integer                            :: nmfuu    ! pointer to mfuu,n
    integer                            :: nml      ! pointer to ml,n
    integer                            :: nmld     ! pointer to mld,n
    integer                            :: nmlu     ! pointer to mlu,n
    integer                            :: nmu      ! pointer to m+1,n
    integer                            :: nmuu     ! pointer to m+2,n
    integer                            :: nu       ! index of point n+1
    integer                            :: num      ! pointer to m,n+1
    integer                            :: numd     ! pointer to m-1,n+1
    integer                            :: nuu      ! index of point n+2
    integer                            :: nuum     ! pointer to m,n+2
    !
    real                               :: bi       ! inverse of main diagonal of the matrix
    real                               :: cfl      ! CFL number
    real                               :: ctrkb    ! contribution of vertical terms below considered point
    real                               :: ctrkt    ! contribution of vertical terms above considered point
    real                               :: dxl      ! local mesh size in x-direction
    real                               :: dyl      ! local mesh size in y-direction
    real                               :: ener     ! total energy of closed system
    real                               :: epsab2   ! small amount to modify weights of AB2 scheme
    real                               :: epslin   ! required accuracy in the linear solver
    real                               :: fac      ! a factor
    real                               :: fac1     ! another factor
    real                               :: fac2     ! some other factor
    real                               :: fac3     ! auxiliary factor
    real                               :: fac4     ! auxiliary factor
    real, dimension(5)                 :: fc       ! energy-conserved weighted interpolation factors for Coriolis force term
    real                               :: fluxlim  ! flux limiter
    real                               :: grad1    ! solution gradient
    real                               :: grad2    ! another solution gradient
    real                               :: gsqsd    ! Jacobian in depth point (m+1/2, n+1/2)
    real                               :: gsqsdb   ! Jacobian in depth point (m+1/2, n-1/2)
    real                               :: gsqsdl   ! Jacobian in depth point (m-1/2, n+1/2)
    real                               :: hd       ! water depth in depth point (m+1/2, n+1/2)
    real                               :: hdb      ! water depth in depth point (m+1/2, n-1/2)
    real                               :: hdl      ! water depth in depth point (m-1/2, n+1/2)
    real                               :: kwd      ! =1. if layer k-1 exists otherwise 0.
    real                               :: kwu      ! =1. if layer k+1 exists otherwise 0.
    real                               :: moutf    ! net mass outflow
    real                               :: qf       ! updated mass flux at present face
    real                               :: reps     ! accuracy of the final approximation
    real                               :: res      ! residual
    real                               :: resm     ! maximum error
    real                               :: rhou     ! density of water in velocity point
    real                               :: rval     ! auxiliary real
    real                               :: rxxs     ! component of Reynolds stress tensor Rxx in current wl-point
    real                               :: rxxsu    ! component of Reynolds stress tensor Rxx in wl-point right from point of consideration
    real                               :: rxxu     ! component of Reynolds stress tensor Rxx in current u/v-point
    real                               :: rxyd     ! component of Reynolds stress tensor Rxy in current depth point
    real                               :: rxydd    ! component of Reynolds stress tensor Rxy in depth point left/down from point of consideration
    real                               :: rxyu     ! component of Reynolds stress tensor Rxy in current u/v-point
    real                               :: ryys     ! component of Reynolds stress tensor Ryy in current wl-point
    real                               :: ryysu    ! component of Reynolds stress tensor Ryy in wl-point up from point of consideration
    real                               :: ryyu     ! component of Reynolds stress tensor Ryy in current u/v-point
    real                               :: s0mx     ! maximum value of water level
    real                               :: sumqf    ! sum of outgoing mass fluxes over present cell
    real                               :: theta3   ! implicitness factor for non-hydrostatic pressure gradient
    real                               :: thetau   ! implicitness factor for vertical terms in u- and v-momentum equations
    real                               :: thetaw   ! implicitness factor for vertical terms in w-momentum equation
    real                               :: u        ! u-velocity at point different from its point of definition
    real                               :: utot     ! velocity magnitude
    real                               :: v        ! v-velocity at point different from its point of definition
    real                               :: vol      ! total displaced volume of water
    real                               :: w        ! w-velocity at point different from its point of definition
    !
    logical                            :: EQREAL   ! compares two reals
    logical                            :: lmask    ! masks points for including curvature stress terms
    logical                            :: STPNOW   ! indicates that program must stop
    !
    character(90)                      :: msgstr   ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashImpLay2DHflow')
    !
    tetau  = pnums( 1)
    tetav  = pnums( 1)
    teta2  = pnums( 4)
    theta3 = pnums( 5)
    thetau = pnums(31)
    thetaw = pnums(32)
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
       mend(1) = ml - 1
       if ( LMXL ) mend(1) = ml
       !
       nend(1) = nl - 1
       if ( LMYL ) nend(1) = nl
       !
       do n = nfu, nend(1)
          do m = mfu, mend(1)
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
    sbimp = .false.
    slimp = .false.
    srimp = .false.
    stimp = .false.
    !
    uwetp = 1
    vwetp = 1
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
    ! build the u-momentum equation
    !
    do n = nfu, nl
       !
       if ( slimp(n) ) then    ! described water level at left boundary point
          !
          msta(n) = mf
          !
       else
          !
          msta(n) = mf + 1
          !
       endif
       !
       if ( srimp(n) .or. lreptx ) then    ! described water level or periodicity at right boundary point
          !
          mend(n) = ml
          !
       else
          !
          mend(n) = ml - 1
          !
       endif
       !
    enddo
    !
    ! initialize system of equations in (permanent) dry points
    !
    do n = nfu, nl
       do m = msta(n), mend(n)
          !
          mu = m + 1
          !
          nm  = kgrpnt(m ,n)
          nmu = kgrpnt(mu,n)
          !
          if ( wetu(nm) /= 1 .or. nmu == 1 ) then
             !
             amatu(nm,:,1) = 1.
             amatu(nm,:,2) = 0.
             amatu(nm,:,3) = 0.
             rhsu (nm,:  ) = 0.
             !
             uwetp(nm) = 0
             !
          endif
          !
       enddo
    enddo
    !
    ! set to zero for permanently dry points
    !
    uwetp(1) = 0
    !
    ! exchange mask array uwetp with neighbouring subdomains
    !
    call SWEXCHGI ( uwetp, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize uwetp at appropriate boundaries in case of repeating grid
    !
    call periodici ( uwetp, kgrpnt, 1, 1 )
    !
    ! compute the time derivative
    !
    do k = 1, kmax
       !
       do n = nfu, nl
          do m = msta(n), mend(n)
             !
             nm = kgrpnt(m,n)
             !
             if ( uwetp(nm) == 1 ) then
                !
                fac = dt * thetau
                !
                amatu(nm,k,1) = 1. / fac
                rhsu (nm,k  ) = u0(nm,k) / fac
                !
             endif
             !
          enddo
       enddo
       !
    enddo
    !
    ! compute the mass flux per width in x-direction
    !
    do k = 1, kmax
       !
       do n = nfu, nl
          do m = mf, ml
             !
             nm = kgrpnt(m,n)
             !
             qx(nm,k) = hku(nm,k)*u0(nm,k)
             !
          enddo
       enddo
       !
    enddo
    !
    ! synchronize mass fluxes at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( qx, kgrpnt, 1, kmax )
    !
    ! compute the mass flux per width in x-direction in wl-point based on averaging
    !
    do k = 1, kmax
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             !
             qxm(nm,k) = 0.5 * ( qx(nm,k) + qx(nmd,k) )
             !
          enddo
       enddo
       !
    enddo
    !
    ! synchronize mass fluxes at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( qxm, kgrpnt, 1, kmax )
    !
    ! compute the advective u-velocity in wl-point based on first order upwind
    !
    do k = 1, kmax
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             !
             if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac = 0.5 * ( u0(nm,k) + u0(nmd,k) )
                !
             else
                !
                fac = qxm(nm,k)
                !
             endif
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm,k) = u0(nmd,k)
                !
             else
                !
                ua(nm,k) = u0(nm ,k)
                !
             endif
             !
          enddo
       enddo
       !
    enddo
    !
    ! synchronize advective u-velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( ua, kgrpnt, 1, kmax )
    !
    ! compute advection term in x-direction (momentum/energy head conservative) at internal u-point (prediction)
    !
    if ( lreptx ) then
       mlast = ml               ! last  internal u-point in case of repeating grid
    else
       mlast = ml - 1           ! last  internal u-point
    endif
    !
    do k = 1, kmax
       !
       do n = nfu, nl
          !
          do m = mf+1, mlast
             !
             md = m - 1
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             nmu = kgrpnt(mu,n)
             !
             if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac1 = 0.5 * ( ua(nmu,k) + ua(nm,k) )
                fac2 = fac1
                !
             else if ( hkumn(nm,k) > 0. ) then
                !
                fac1 = qxm(nmu,k) / hkumn(nm,k)
                fac2 = qxm(nm ,k) / hkumn(nm,k)
                !
             else
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( uwetp(nm) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) - ( fac1 * (ua(nmu,k) - u0(nm,k)) - fac2 * (ua(nm,k) - u0(nm,k)) ) / gvu(nm)
                !
             endif
             !
          enddo
          !
          if ( slimp(n) ) then    ! described water level at left boundary point
             !
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             !
             if ( qxm(nmfu,k) < 0. ) then
                !
                rhsu(nmf,k) = rhsu(nmf,k) - qxm(nmfu,k) * (u0(nmfu,k) - u0(nmf,k)) /  ( gvu(nmf) * hkum(nmf,k) )
                !
             endif
             !
          endif
          !
          if ( srimp(n) ) then    ! described water level at right boundary point
             !
             nml  = kgrpnt(ml  ,n)
             nmld = kgrpnt(ml-1,n)
             !
             if ( .not. qxm(nml,k) < 0. ) then
                !
                rhsu(nml,k) = rhsu(nml,k) - qxm(nml,k) * (u0(nml,k) - u0(nmld,k)) / ( gvu(nml) * hkum(nml,k) )
                !
             endif
             !
          endif
          !
       enddo
       !
    enddo
    !
    ! compute the mass flux per width in y-direction
    !
    do k = 1, kmax
       !
       do n = nf, nl
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             qy(nm,k) = hkv(nm,k)*v0(nm,k)
             !
          enddo
       enddo
       !
    enddo
    !
    ! synchronize mass fluxes at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( qy, kgrpnt, 1, kmax )
    !
    ! compute the mass flux per width in y-direction in depth point based on averaging
    !
    do k = 1, kmax
       !
       do n = nf, nl
          do m = mf+1, mlast
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             qym(nm,k) = 0.5 * ( qy(nm,k) + qy(nmu,k) )
             !
          enddo
       enddo
       !
    enddo
    !
    ! compute the advective u-velocity in depth point based on first order upwind
    !
    do k = 1, kmax
       !
       do n = nf, nl
          do m = mf+1, mlast
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
             if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac = 0.5 * ( v0(nm,k) + v0(nmu,k) )
                !
             else
                !
                fac = qym(nm,k)
                !
             endif
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm,k) = u0(nm ,k)
                !
             else
                !
                ua(nm,k) = u0(num,k)
                !
             endif
             !
          enddo
       enddo
       !
    enddo
    !
    ! compute advection term in y-direction (momentum/energy head conservative) at internal u-point (prediction)
    !
    do k = 1, kmax
       !
       do n = nfu, nl
          !
          do m = mf+1, mlast
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
             if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac1 = 0.5 * ( v0(nm,k) + v0(nmu,k) )
                fac2 = fac1
                !
             else if ( hkumn(nm,k) > 0. ) then
                !
                fac1 = qym(nm ,k) / hkumn(nm,k)
                fac2 = qym(ndm,k) / hkumn(nm,k)
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
             if ( uwetp(nm) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) - ( fac1 * (ua(nm,k) - u0(nm,k)) - fac2 * (ua(ndm,k) - u0(nm,k)) ) / guu(nm)
                !
             endif
             !
          enddo
          !
       enddo
       !
    enddo
    !
    ! compute the curvature terms associated with advection terms (internal points only!)
    !
    if ( optg == 3 ) then
       !
       do k = 1, kmax
          !
          do n = nfu, nl
             do m = mf+1, mlast
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
                if ( uwetp(nm) == 1 ) then
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd  == 1 ) nmd  = nm
                   if ( nmu  == 1 ) nmu  = nm
                   if ( ndm  == 1 ) ndm  = nm
                   if ( ndmu == 1 ) ndmu = ndm
                   !
                   v = 0.25 * ( v0(nm,k) + v0(nmu,k) + v0(ndm,k) + v0(ndmu,k) )
                   !
                   rhsu(nm,k) = rhsu(nm,k) - 0.5 * u0(nm,k) * v * ( gvv(nm) + gvv(nmu) - gvv(ndm) - gvv(ndmu) ) / gsqsu(nm)
                   !
                   if ( (ibb(m) == 1 .or. n /= nfu .or. .not.LMYF) .and. (ibt(m) == 1 .or. n /= nl .or. .not.LMYL) ) rhsu(nm,k) = rhsu(nm,k) + 0.5 * v * v * ( guu(nmu) - guu(nmd) ) / gsqsu(nm)
                   !
                endif
                !
             enddo
          enddo
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
          do n = nfu, nl
             do m = msta(n), mend(n)
                !
                mu = m + 1
                nd = n - 1
                !
                nm   = kgrpnt(m ,n )
                nmu  = kgrpnt(mu,n )
                ndm  = kgrpnt(m ,nd)
                ndmu = kgrpnt(mu,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmu  == 1 ) nmu  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmu == 1 ) ndmu = ndm
                !
                if ( uwetp(nm) == 1 ) then
                   !
                   v = 0.25 * ( v0(nm,k  ) + v0(nmu,k  ) + v0(ndm,k) + v0(ndmu,k) )
                   w = 0.25 * ( w0(nm,k-1) + w0(nmu,k-1) + w0(nm ,k) + w0(nmu ,k) )
                   !
                   utot = sqrt( u0(nm,k)*u0(nm,k) + v*v + w*w )
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
          do n = nfu, nl
             do m = msta(n), mend(n)
                !
                mu = m + 1
                nd = n - 1
                !
                nm   = kgrpnt(m ,n )
                nmu  = kgrpnt(mu,n )
                ndm  = kgrpnt(m ,nd)
                ndmu = kgrpnt(mu,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmu  == 1 ) nmu  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmu == 1 ) ndmu = ndm
                !
                if ( uwetp(nm) == 1 ) then
                   !
                   v = 0.25 * ( v0(nm,k) + v0(nmu,k) + v0(ndm,k) + v0(ndmu,k) )
                   !
                   utot = sqrt( u0(nm,k)*u0(nm,k) + v*v )
                   !
                   amatu(nm,k,1) = amatu(nm,k,1) + cvegu(nm,k,1) * utot
                   !
                endif
                !
             enddo
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
             do n = nfu, nl
                do m = msta(n), mend(n)
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( uwetp(nm) == 1 ) then
                      !
                      fac = dt * thetau
                      !
                      amatu(nm,k,1) = amatu(nm,k,1) + cvegu(nm,k,2) / fac
                      rhsu (nm,k  ) = rhsu (nm,k  ) + cvegu(nm,k,2) * u0(nm,k) / fac
                      !
                   endif
                   !
                enddo
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
    do n = nfu, nl
       do m = msta(n), mend(n)
          !
          mu = m + 1
          !
          nm  = kgrpnt(m ,n)
          nmu = kgrpnt(mu,n)
          !
          if ( uwetp(nm) == 1 ) then
             !
             rhsu(nm,:) = rhsu(nm,:) - grav * (s0(nmu) - s0(nm)) / gvu(nm)
             !
          endif
          !
       enddo
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
                   fac = 1. / ( 2. * gvu(nm) * hkum(nm,k) )
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
                         fac1 = 0.5 * fac1 / ( gvu(nm) * ( hkum(nm,kd) + hkum(nm,k ) ) )
                         fac2 = 0.5 * fac2 / ( gvu(nm) * ( hkum(nm,k ) + hkum(nm,ku) ) )
                         !
                         gmatu(nm,k,1) = -fac1
                         gmatu(nm,k,2) = ( -hks(nm ,k)/gvu(nm) - hkum(nm,kd)*fac1 + hkum(nm,ku)*fac2 ) / hkum(nm,k)
                         gmatu(nm,k,3) = -fac1
                         gmatu(nm,k,4) = (  hks(nmu,k)/gvu(nm) - hkum(nm,kd)*fac1 + hkum(nm,ku)*fac2 ) / hkum(nm,k)
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
       enddo
       !
       do n = nfu, nl
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
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
       enddo
       !
       ! to reduce the pressure Poisson equation set pressure of bottom face to that of top face for a number of layers
       !
       do l = 1, qlay
          !
          do n = nfu, nl
             do m = mf, ml
                !
                nm = kgrpnt(m,n)
                !
                gmatu(nm,kmax-l+1,1) = gmatu(nm,kmax-l+1,1) + gmatu(nm,kmax-l+1,2)
                gmatu(nm,kmax-l+1,2) = 0.
                gmatu(nm,kmax-l+1,3) = gmatu(nm,kmax-l+1,3) + gmatu(nm,kmax-l+1,4)
                gmatu(nm,kmax-l+1,4) = 0.
                !
             enddo
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
                   fac1 = 0.5 * ( zks(nmu,k-1) - zks(nm,k-1) ) / ( gvu(nm) * ( hkum(nm,kd) + hkum(nm,k ) ) )
                   fac2 = 0.5 * ( zks(nmu,k  ) - zks(nm,k  ) ) / ( gvu(nm) * ( hkum(nm,k ) + hkum(nm,ku) ) )
                   !
                   gmatu(nm,k,1) = -fac1
                   gmatu(nm,k,2) = ( -hks(nm ,k)/gvu(nm) - hkum(nm,kd)*fac1 + hkum(nm,ku)*fac2 ) / hkum(nm,k)
                   gmatu(nm,k,3) = -fac1
                   gmatu(nm,k,4) = (  hks(nmu,k)/gvu(nm) - hkum(nm,kd)*fac1 + hkum(nm,ku)*fac2 ) / hkum(nm,k)
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
       enddo
       !
       do n = nfu, nl
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
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
             do n = nfu, nl
                do m = msta(n), mend(n)
                   !
                   mu = m + 1
                   !
                   nm  = kgrpnt(m ,n)
                   nmu = kgrpnt(mu,n)
                   !
                   if ( uwetp(nm) == 1 ) then
                      !
                      rhsu(nm,k) = rhsu(nm,k) - gmatu(nm,k,1)*q(nm ,kd) - gmatu(nm,k,2)*q(nm,k ) - gmatu(nm,k,3)*q(nmu,kd)  &
                                              - gmatu(nm,k,4)*q(nmu,k ) - gmatu(nm,k,5)*q(nm,ku) - gmatu(nm,k,6)*q(nmu,ku)
                      !
                   endif
                   !
                enddo
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
             do n = nfu, nl
                do m = msta(n), mend(n)
                   !
                   mu = m + 1
                   !
                   nm  = kgrpnt(m ,n)
                   nmu = kgrpnt(mu,n)
                   !
                   if ( uwetp(nm) == 1 ) then
                      !
                      rhsu(nm,k) = rhsu(nm,k) - (1.-theta3) * ( gmatu(nm,k,1)*q(nm ,kd) + gmatu(nm,k,2)*q(nm,k ) + gmatu(nm,k,3)*q(nmu,kd)  &
                                                              + gmatu(nm,k,4)*q(nmu,k ) + gmatu(nm,k,5)*q(nm,ku) + gmatu(nm,k,6)*q(nmu,ku) )
                      !
                   endif
                   !
                enddo
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
          do n = nfu, nl
             do m = mf+1, mlast
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                if ( uwetp(nm) == 1 ) then
                   !
                   rhsu(nm,k) = rhsu(nm,k) - 0.5 * grav * hkum(nm,k) * (rho(nmu,k) - rho(nm,k)) / ( rhow*gvu(nm) )
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       do k = 2, kmax
          !
          do n = nfu, nl
             do m = mf+1, mlast
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                if ( uwetp(nm) == 1 ) then
                   !
                   rhou = rhow + 0.5 * (rho(nmu,k-1) + rho(nm,k-1))
                   !
                   work(nm,1) = work(nm,1) + rhou * (hks(nmu,k-1) - hks(nm,k-1)) + hkum(nm,k-1) * (rho(nmu,k-1) - rho(nm,k-1))
                   !
                   work(nm,2) = work(nm,2) + hks(nmu,k-1) - hks(nm,k-1)
                   !
                   rhou = rhow + 0.5 * (rho(nmu,k) + rho(nm,k))
                   !
                   rhsu(nm,k) = rhsu(nm,k) - grav * ( work(nm,1) - rhou*work(nm,2) ) / ( rhow*gvu(nm) )
                   !
                endif
                !
             enddo
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
       do n = nfu, nl
          do m = mf+1, mlast
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( uwetp(nm) == 1 ) then
                !
                rhsu(nm,:) = rhsu(nm,:) - (patm(nmu) - patm(nm)) / ( rhow*gvu(nm) )
                !
             endif
             !
          enddo
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
          do n = nfu, nl
             !
             do m = mf+1, mlast
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
                if ( uwetp(nm) == 1 ) then
                   !
                   if ( wetu(nmd) * wetu(nmu) == 1 ) then
                      !
                      fac1 = hks(nmu,k) * ( guu(nmu)*u0(nmu,k) - guu(nm)*u0(nm,k) )/gsqs(nmu) - hks(nm,k) * ( guu(nm)*u0(nm,k) - guu(nmd)*u0(nmd,k) )/gsqs(nm)
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
                      hd  = .5 * (hkum(num,k) + hkum(nm ,k))
                      hdb = .5 * (hkum(nm ,k) + hkum(ndm,k))
                      !
                      fac2 = hd * ( gvu(num)*u0(num,k) - gvu(nm)*u0(nm,k) )/gsqsd - hdb * ( gvu(nm)*u0(nm,k) - gvu(ndm)*u0(ndm,k) )/gsqsdb
                      !
                   else
                      !
                      fac2 = 0.
                      !
                   endif
                   !
                   rhsu(nm,k) = rhsu(nm,k) + hvisc * ( fac1/gvu(nm) + fac2/guu(nm) ) / hkum(nm,k)
                   !
                endif
                !
             enddo
             !
          enddo
          !
       enddo
       !
    else if ( ihvisc > 1 ) then
       !
       do k = 1, kmax
          !
          do n = nfu, nl
             !
             do m = mf+1, mlast
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
                if ( uwetp(nm) == 1 ) then
                   !
                   lmask = .true.
                   !
                   hd  = .5 * (hkum(num,k) + hkum(nm ,k))
                   hdb = .5 * (hkum(nm ,k) + hkum(ndm,k))
                   !
                   ! first, compute the components of the Reynolds stress tensor in different points of consideration
                   !
                   if ( wetu(nmd) * wetu(nmu) == 1 ) then
                      !
                      rxxs  = ( vnu2d(nm ) + vnu2d(ndm ) + vnu2d(nmd) + vnu2d(ndmd) ) * ( u0(nm ,k) - u0(nmd,k) ) / ( gvv(nm ) + gvv(ndm ) )
                      !
                      rxxsu = ( vnu2d(nmu) + vnu2d(ndmu) + vnu2d(nm ) + vnu2d(ndm ) ) * ( u0(nmu,k) - u0(nm ,k) ) / ( gvv(nmu) + gvv(ndmu) )
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
                      rxyd  = 2. * vnu2d(nm ) * ( ( v0(nmu ,k) - v0(nm ,k) )/( gvv(nmu ) + gvv(nm ) ) + ( u0(num,k) - u0(nm ,k) )/( guu(num) + guu(nm ) ) )
                      !
                      rxydd = 2. * vnu2d(ndm) * ( ( v0(ndmu,k) - v0(ndm,k) )/( gvv(ndmu) + gvv(ndm) ) + ( u0(nm ,k) - u0(ndm,k) )/( guu(nm ) + guu(ndm) ) )
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
                      rxxu  = 0.50 * ( vnu2d(nm) + vnu2d(ndm) ) * ( u0(nmu,k) - u0(nmd,k) ) / gvu(nm)
                      !
                      rxyu  = 0.25 * ( vnu2d(nm) + vnu2d(ndm) ) * ( ( v0(nmu,k) + v0(ndmu,k) - v0(nm,k) - v0(ndm,k) )/gvu(nm) + ( u0(num,k) - u0(ndm,k) )/guu(nm) )
                      !
                      ryyu  = 0.50 * ( vnu2d(nm) + vnu2d(ndm) ) * ( v0(nmu,k) + v0(nm,k) - v0(ndm,k) -v0(ndmu,k) ) / guu(nm)
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
                   rhsu(nm,k) = rhsu(nm,k) + ( hks(nmu,k)*rxxsu - hks(nm,k)*rxxs ) / (hkum(nm,k)*gvu(nm)) + ( hd*rxyd - hdb*rxydd ) / (hkum(nm,k)*guu(nm))
                   !
                   ! finally, add the curvature terms associated with those stresses
                   !
                   rhsu(nm,k) = rhsu(nm,k) + ( rxyu * ( gvv(nm) + gvv(nmu) - gvv(ndm) - gvv(ndmu) ) + 0.5 * (rxxu - ryyu) * ( guu(nmu) - guu(nmd) ) ) / gsqsu(nm)
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
    ! add contribution due to full 3D turbulence
    !
    if ( iturb > 1 ) then
       !
       ! compute divergence of Reynolds stress tensor at internal u-point
       !
       do k = 1, kmax
          !
          do n = nfu, nl
             do m = mf+1, mlast
                !
                mu = m + 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmu = kgrpnt(mu,n )
                ndm = kgrpnt(m ,nd)
                !
                if ( uwetp(nm) == 1 ) then
                   !
                   rhsu(nm,k) = rhsu(nm,k) + ( rsuu(nmu,k  ) - rsuu(nm ,k) ) / gvu(nm)     &
                                           + ( rsuv(nm ,k  ) - rsuv(ndm,k) ) / guu(nm)     &
                                           + ( rsuw(nm ,k-1) - rsuw(nm ,k) ) / hkum(nm,k)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! impose wall shear stress at closed boundaries
       !
       if ( .not.lrepty ) then
          !
          do m = mf+1, mlast
             !
             nfum = kgrpnt(m,nfu)
             nlm  = kgrpnt(m,nl )
             !
             if ( ibb(m) == 1 .and. LMYF ) then
                !
                amatu(nfum,:,1) = amatu(nfum,:,1) + logwlb(m,:) / guu(nfum)
                !
             endif
             !
             if ( ibt(m) == 1 .and. LMYL ) then
                !
                amatu(nlm,:,1) = amatu(nlm,:,1) + logwlt(m,:) / guu(nlm)
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
    ! compute explicit part of wind stress term at internal u-point (top layer only), if appropriate
    !
    if ( iwind /= 0 ) then
       !
       do n = nfu, nl
          do m = mf+1, mlast
             !
             nm = kgrpnt(m,n)
             !
             if ( uwetp(nm) == 1 ) then
                !
                rhsu(nm,1) = rhsu(nm,1) + windu(nm)/max(1.e-3,hkum(nm,1))
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute implicit part of wind stress term at internal u-point (top layer only), if appropriate
    !
    if ( relwnd ) then
       !
       do n = nfu, nl
          do m = mf+1, mlast
             !
             nm = kgrpnt(m,n)
             !
             if ( uwetp(nm) == 1 ) then
                !
                amatu(nm,1,1) = amatu(nm,1,1) + cwndu(nm) / max(1.e-3,hkum(nm,1))
                !
             endif
             !
          enddo
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
       do n = nfu, nl
          do m = msta(n), mend(n)
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( uwetp(nm) == 1 ) then
                !
                amatu(nm,kmax,1) = amatu(nm,kmax,1) + 0.5 * ( logfrc(nm,1) + logfrc(nmu,1) ) / hkum(nm,kmax)
                !
             endif
             !
          enddo
       enddo
       !
    else if ( irough == 11 ) then
       !
       ! linear bottom friction
       !
       do n = nfu, nl
          do m = msta(n), mend(n)
             !
             nm = kgrpnt(m,n)
             !
             if ( uwetp(nm) == 1 ) then
                !
                amatu(nm,kmax,1) = amatu(nm,kmax,1) + cfricu(nm) / hkum(nm,kmax)
                !
             endif
             !
          enddo
       enddo
       !
    else if ( irough /= 0 ) then
       !
       do n = nfu, nl
          do m = msta(n), mend(n)
             !
             mu = m + 1
             nd = n - 1
             !
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             ndmu = kgrpnt(mu,nd)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmu  == 1 ) nmu  = nm
             if ( ndm  == 1 ) ndm  = nm
             if ( ndmu == 1 ) ndmu = ndm
             !
             v = 0.25 * ( v0(nm,kmax) + v0(nmu,kmax) + v0(ndm,kmax) + v0(ndmu,kmax) )
             !
             utot = sqrt( u0(nm,kmax)*u0(nm,kmax) + v*v )
             !
             if ( uwetp(nm) == 1 .and. utot > 1.e-8 ) then
                !
                v = 0.25 * ( vdep(nm) + vdep(nmu) + vdep(ndm) + vdep(ndmu) )
                !
                fac = udep(nm)*udep(nm) + v*v
                !
                amatu(nm,kmax,1) = amatu(nm,kmax,1) + cfricu(nm) * fac / ( utot * hkum(nm,kmax) )
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute Coriolis force, if appropriate
    !
    if ( coriolis ) then
       !
       epsab2 = pcor(2)
       !
       do n = nfu, nl
          do m = msta(n), mend(n)
             !
             mu = m + 1
             nd = n - 1
             !
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             ndmu = kgrpnt(mu,nd)
             !
             if ( uwetp(nm) == 1 .and. fcor(nm,1) /= 0. ) then
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmu  == 1 ) nmu  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmu == 1 ) ndmu = ndm
                !
                ! energy-conserved weighting factors according to Espelid et al. (2000)
                !
                fc(1) = sqrt(abs(fcor(nm  ,2))*dpv(nm  ))
                fc(2) = sqrt(abs(fcor(nmu ,2))*dpv(nmu ))
                fc(3) = sqrt(abs(fcor(ndm ,2))*dpv(ndm ))
                fc(4) = sqrt(abs(fcor(ndmu,2))*dpv(ndmu))
                fc(5) = 1. / (4.*sqrt(abs(fcor(nm,1))*dpu(nm)))
                !
                do k = 1, kmax
                   !
                   ! at time level n-1
                   !
                   ! --- modified quasi second order Adams-Bashforth scheme which is conditionally stable for the inviscid case (see Marshall et al., 1997)
                   !
                   rhsu(nm,k) = rhsu(nm,k) + (0.5+epsab2) * cfu(nm,k)
                   !
                   ! at time level n
                   !
                   ! --- interpolate v-velocity at u-point using energy-conserved weighting factors
                   !
                   v = fc(1) * v0(nm,k) + fc(2) * v0(nmu,k) + fc(3) * v0(ndm,k) + fc(4) * v0(ndmu,k)
                   v = fc(5) * v
                   !
                   cfu(nm,k) = -fcor(nm,1) * v
                   !
                   rhsu(nm,k) = rhsu(nm,k) - (1.5+epsab2) * cfu(nm,k)
                   !
                enddo
                !
             endif
             !
          enddo
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
       do n = nfu, nl
          do m = msta(n), mend(n)
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( uwetp(nm) == 1 ) then
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
       enddo
       !
    enddo
    !
    ! solve the u-momentum equation
    !
    do n = nfu, nl
       do m = msta(n), mend(n)
          !
          nm = kgrpnt(m,n)
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
    enddo
    !
    ! re-update the solution in case of thetau <> 1
    !
    if ( thetau /= 1. ) then
       !
       do k = 1, kmax
          !
          do n = nfu, nl
             do m = msta(n), mend(n)
                !
                nm = kgrpnt(m,n)
                !
                if ( uwetp(nm) == 1 ) then
                   !
                   u1(nm,k) = ( u1(nm,k) - (1.-thetau) * u0(nm,k) ) / thetau
                   !
                endif
                !
             enddo
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
    ! synchronize u-velocities at appropriate boundaries in case of repeating grid
    !
    call periodic ( u1, kgrpnt, 1, kmax )
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
                if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
                   !
                   fac = 0.5 * ( u0(nm,k) + u0(nmd,k) )
                   !
                else
                   !
                   fac = qxm(nm,k)
                   !
                endif
                !
                if ( .not. fac < 0. ) then
                   !
                   grad1 = u1(nm ,k) - u0(nmd ,k)
                   grad2 = u1(nmd,k) - u0(nmdd,k)
                   !
                   up(nm,k) = u0(nmd,k)
                   ua(nm,k) = 0.5 * fluxlim(grad1,grad2)
                   !
                else
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
          enddo
          !
       enddo
       !
       call periodic ( ua, kgrpnt, 1, kmax )
       call periodic ( up, kgrpnt, 1, kmax )
       !
    endif
    !
    ! compute higher order correction to the advective velocity in depth point (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       if ( lreptx ) then
          mlast = ml            ! last  internal u-point in case of repeating grid
       else
          mlast = ml - 1        ! last  internal u-point
       endif
       !
       do k = 1, kmax
          !
          do n = nf, nl
             do m = mf+1, mlast
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
                if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
                   !
                   fac = 0.5 * ( v0(nm,k) + v0(nmu,k) )
                   !
                else
                   !
                   fac = qym(nm,k)
                   !
                endif
                !
                if ( .not. fac < 0. ) then
                   !
                   grad1 = u1(num,k) - u0(nm ,k)
                   grad2 = u1(nm ,k) - u0(ndm,k)
                   !
                   ua2(nm,k) = 0.5 * fluxlim(grad1,grad2)
                   !
                else
                   !
                   grad1 = u1(nm ,k) - u0(num ,k)
                   grad2 = u1(num,k) - u0(nuum,k)
                   !
                   ua2(nm,k) = 0.5 * fluxlim(grad1,grad2)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
    endif
    !
    ! correct u-velocity associated with higher order approximation of advection term in x-direction at internal u-point (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do k = 1, kmax
          !
          do n = nfu, nl
             do m = mf+1, mlast
                !
                md = m - 1
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmd = kgrpnt(md,n)
                nmu = kgrpnt(mu,n)
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
                   fac1 = qxm(nmu,k) / hkumn(nm,k)
                   fac2 = qxm(nm ,k) / hkumn(nm,k)
                   !
                else
                   !
                   fac1 = 0.
                   fac2 = 0.
                   !
                endif
                !
                if ( uwetp(nm) == 1 .and. presu(nm) /= 1 ) then
                   !
                   u1(nm,k) = u1(nm,k) - dt * ( fac1 * ua(nmu,k)  - fac2 * ua(nm,k) ) / gvu(nm)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
    endif
    !
    ! correct u-velocity associated with higher order approximation of advection term in y-direction at internal u-point (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do k = 1, kmax
          !
          do n = nfu, nl
             do m = mf+1, mlast
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
                if ( ( udep(nm) > udep(nmd) .or. stricthead ) .and. .not.strictmom ) then
                   !
                   !fac1 = 0.5 * ( v0(nm,k) + v0(nmu,k) )
                   !fac2 = fac1
                   !jumps back to old implementation (020117)
                   fac1 = 0.
                   fac2 = 0.
                   !
                else if ( hkumn(nm,k) > 0. ) then
                   !
                   fac1 = qym(nm ,k) / hkumn(nm,k)
                   fac2 = qym(ndm,k) / hkumn(nm,k)
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
                if ( uwetp(nm) == 1 .and. presu(nm) /= 1 ) then
                   !
                   u1(nm,k) = u1(nm,k) - dt * ( fac1 * ua2(nm,k) - fac2 * ua2(ndm,k) ) / guu(nm)
                   !
                endif
                !
             enddo
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
    ! synchronize u-velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    if ( propsc /= 1 ) call periodic ( u1, kgrpnt, 1, kmax )
    !
    ! build the v-momentum equation
    !
    do m = mfu, ml
       !
       if ( sbimp(m) ) then    ! described water level at lower boundary point
          !
          nsta(m) = nf
          !
       else
          !
          nsta(m) = nf + 1
          !
       endif
       !
       if ( stimp(m) .or. lrepty ) then    ! described water level or periodicity at upper boundary point
          !
          nend(m) = nl
          !
       else
          !
          nend(m) = nl - 1
          !
       endif
       !
    enddo
    !
    ! initialize system of equations in dry points
    !
    do m = mfu, ml
       do n = nsta(m), nend(m)
          !
          nu = n + 1
          !
          nm  = kgrpnt(m,n )
          num = kgrpnt(m,nu)
          !
          if ( wetv(nm) /= 1 .or. num == 1 ) then
             !
             amatu(nm,:,1) = 1.
             amatu(nm,:,2) = 0.
             amatu(nm,:,3) = 0.
             rhsu (nm,:  ) = 0.
             !
             vwetp(nm) = 0
             !
          endif
          !
       enddo
    enddo
    !
    ! set to zero for permanently dry points
    !
    vwetp(1) = 0
    !
    ! exchange mask array vwetp with neighbouring subdomains
    !
    call SWEXCHGI ( vwetp, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize vwetp at appropriate boundaries in case of repeating grid
    !
    call periodici ( vwetp, kgrpnt, 1, 1 )
    !
    ! compute the time derivative
    !
    do k = 1, kmax
       !
       do m = mfu, ml
          do n = nsta(m), nend(m)
             !
             nm = kgrpnt(m,n)
             !
             if ( vwetp(nm) == 1 ) then
                !
                fac = dt * thetau
                !
                amatu(nm,k,1) = 1. / fac
                rhsu (nm,k  ) = v0(nm,k) / fac
                !
             endif
             !
          enddo
       enddo
       !
    enddo
    !
    ! compute the mass flux per width in x-direction in depth point based on averaging
    !
    if ( lrepty ) then
       nlast = nl               ! last  internal v-point in case of repeating grid
    else
       nlast = nl - 1           ! last  internal v-point
    endif
    !
    do k = 1, kmax
       !
       do m = mf, ml
          do n = nf+1, nlast
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             qxm(nm,k) = 0.5 * ( qx(nm,k) + qx(num,k) )
             !
          enddo
       enddo
       !
    enddo
    !
    ! compute the advective v-velocity in depth point based on first order upwind
    !
    do k = 1, kmax
       !
       do m = mf, ml
          do n = nf+1, nlast
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
             if ( ( vdep(nm) > vdep(ndm) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac = 0.5 * ( u0(nm,k) + u0(num,k) )
                !
             else
                !
                fac = qxm(nm,k)
                !
             endif
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm,k) = v0(nm ,k)
                !
             else
                !
                ua(nm,k) = v0(nmu,k)
                !
             endif
             !
          enddo
       enddo
       !
    enddo
    !
    ! compute advection term in x-direction (momentum/energy head conservative) at internal v-point (prediction)
    !
    do k = 1, kmax
       !
       do m = mfu, ml
          !
          do n = nf+1, nlast
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
             if ( ( vdep(nm) > vdep(ndm) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac1 = 0.5 * ( u0(nm,k) + u0(num,k) )
                fac2 = fac1
                !
             else if ( hkvmn(nm,k) > 0. ) then
                !
                fac1 = qxm(nm ,k) / hkvmn(nm,k)
                fac2 = qxm(nmd,k) / hkvmn(nm,k)
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
             if ( vwetp(nm) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) - ( fac1 * (ua(nm,k) - v0(nm,k)) - fac2 * (ua(nmd,k) - v0(nm,k)) ) / gvv(nm)
                !
             endif
             !
          enddo
          !
       enddo
       !
    enddo
    !
    ! compute the mass flux per width in y-direction in wl-point based on averaging
    !
    do k = 1, kmax
       !
       do m = mfu, ml
          do n = nfu, nl
             !
             nd = n - 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             !
             qym(nm,k) = 0.5 * ( qy(nm,k) + qy(ndm,k) )
             !
          enddo
       enddo
       !
    enddo
    !
    ! synchronize mass fluxes at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( qym, kgrpnt, 1, kmax )
    !
    ! compute the advective v-velocity in wl-point based on first order upwind
    !
    do k = 1, kmax
       !
       do m = mfu, ml
          do n = nfu, nl
             !
             nd = n - 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             !
             if ( ( vdep(nm) > vdep(ndm) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac = 0.5 * ( v0(nm,k) + v0(ndm,k) )
                !
             else
                !
                fac = qym(nm,k)
                !
             endif
             !
             if ( .not. fac < 0. ) then
                !
                ua(nm,k) = v0(ndm,k)
                !
             else
                !
                ua(nm,k) = v0(nm ,k)
                !
             endif
             !
          enddo
       enddo
       !
    enddo
    !
    ! synchronize advective v-velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( ua, kgrpnt, 1, kmax )
    !
    ! compute advection term in y-direction (momentum/energy head conservative) at internal v-point (prediction)
    !
    do k = 1, kmax
       !
       do m = mfu, ml
          !
          do n = nf+1, nlast
             !
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             num = kgrpnt(m,nu)
             !
             if ( ( vdep(nm) > vdep(ndm) .or. stricthead ) .and. .not.strictmom ) then
                !
                fac1 = 0.5 * ( ua(num,k) + ua(nm,k) )
                fac2 = fac1
                !
             else if ( hkvmn(nm,k) > 0. ) then
                !
                fac1 = qym(num,k) / hkvmn(nm,k)
                fac2 = qym(nm ,k) / hkvmn(nm,k)
                !
             else
                !
                fac1 = 0.
                fac2 = 0.
                !
             endif
             !
             if ( vwetp(nm) == 1 ) then
                !
                rhsu(nm,k) = rhsu(nm,k) - ( fac1 * (ua(num,k) - v0(nm,k)) - fac2 * (ua(nm,k) - v0(nm,k)) ) / guv(nm)
                !
             endif
             !
          enddo
          !
          if ( sbimp(m) ) then    ! described water level at lower boundary point
             !
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             !
             if ( qym(nfum,k) < 0. ) then
                !
                rhsu(nfm,k) = rhsu(nfm,k) - qym(nfum,k) * (v0(nfum,k) - v0(nfm,k)) / ( guv(nfm) * hkvm(nfm,k) )
                !
             endif
             !
          endif
          !
          if ( stimp(m) ) then    ! described water level at upper boundary point
             !
             nlm  = kgrpnt(m,nl  )
             nldm = kgrpnt(m,nl-1)
             !
             if ( .not. qym(nlm,k) < 0. ) then
                !
                rhsu(nlm,k) = rhsu(nlm,k) - qym(nlm,k) * (v0(nlm,k) - v0(nldm,k)) / ( guv(nlm) * hkvm(nlm,k) )
                !
             endif
             !
          endif
          !
       enddo
       !
    enddo
    !
    ! compute the curvature terms associated with advection terms (internal points only!)
    !
    if ( optg == 3 ) then
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             do n = nf+1, nlast
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
                if ( vwetp(nm) == 1 ) then
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( ndm  == 1 ) ndm  = nm
                   if ( num  == 1 ) num  = nm
                   if ( nmd  == 1 ) nmd  = nm
                   if ( numd == 1 ) numd = nmd
                   !
                   u = 0.25 * ( u0(nm,k) + u0(num,k) + u0(nmd,k) + u0(numd,k) )
                   !
                   if ( (ibl(n) == 1 .or. m /= mfu .or. .not.LMXF) .and. (ibr(n) == 1 .or. m /= ml .or. .not.LMXL) ) rhsu(nm,k) = rhsu(nm,k) + 0.5 * u * u * ( gvv(num) - gvv(ndm) ) / gsqsv(nm)
                   !
                   rhsu(nm,k) = rhsu(nm,k) - 0.5 * v0(nm,k) * u * ( guu(nm) + guu(num) - guu(nmd) - guu(numd) ) / gsqsv(nm)
                   !
                endif
                !
             enddo
          enddo
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
          do m = mfu, ml
             do n = nsta(m), nend(m)
                !
                nu = n + 1
                md = m - 1
                !
                nm   = kgrpnt(m ,n )
                num  = kgrpnt(m ,nu)
                nmd  = kgrpnt(md,n )
                numd = kgrpnt(md,nu)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( num  == 1 ) num  = nm
                if ( nmd  == 1 ) nmd  = nm
                if ( numd == 1 ) numd = nmd
                !
                if ( vwetp(nm) == 1 ) then
                   !
                   u = 0.25 * ( u0(nm,k  ) + u0(num,k  ) + u0(nmd,k) + u0(numd,k) )
                   w = 0.25 * ( w0(nm,k-1) + w0(num,k-1) + w0(nm ,k) + w0(num ,k) )
                   !
                   utot = sqrt( u*u + v0(nm,k)*v0(nm,k) +w*w )
                   !
                   amatu(nm,k,1) = amatu(nm,k,1) + apomv(nm,k) + bpomv(nm,k) * utot
                   !
                   fac = dt * thetau
                   !
                   amatu(nm,k,1) = amatu(nm,k,1) + cpomv(nm,k) / fac
                   rhsu (nm,k  ) = rhsu (nm,k  ) + cpomv(nm,k) * v0(nm,k) / fac
                   !
                endif
                !
             enddo
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
          do m = mfu, ml
             do n = nsta(m), nend(m)
                !
                nu = n + 1
                md = m - 1
                !
                nm   = kgrpnt(m ,n )
                num  = kgrpnt(m ,nu)
                nmd  = kgrpnt(md,n )
                numd = kgrpnt(md,nu)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( num  == 1 ) num  = nm
                if ( nmd  == 1 ) nmd  = nm
                if ( numd == 1 ) numd = nmd
                !
                if ( vwetp(nm) == 1 ) then
                   !
                   u = 0.25 * ( u0(nm,k) + u0(num,k) + u0(nmd,k) + u0(numd,k) )
                   !
                   utot = sqrt( u*u + v0(nm,k)*v0(nm,k) )
                   !
                   amatu(nm,k,1) = amatu(nm,k,1) + cvegv(nm,k,1) * utot
                   !
                endif
                !
             enddo
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
             do m = mfu, ml
                do n = nsta(m), nend(m)
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( vwetp(nm) == 1 ) then
                      !
                      fac = dt * thetau
                      !
                      amatu(nm,k,1) = amatu(nm,k,1) + cvegv(nm,k,2) / fac
                      rhsu (nm,k  ) = rhsu (nm,k  ) + cvegv(nm,k,2) * v0(nm,k) / fac
                      !
                   endif
                   !
                enddo
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
    do m = mfu, ml
       do n = nsta(m), nend(m)
          !
          nu = n + 1
          !
          nm  = kgrpnt(m,n )
          num = kgrpnt(m,nu)
          !
          if ( vwetp(nm) == 1 ) then
             !
             rhsu(nm,:) = rhsu(nm,:) - grav * (s0(num) - s0(nm)) / guv(nm)
             !
          endif
          !
       enddo
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
                   fac = 1. / ( 2. * guv(nm) * hkvm(nm,k) )
                   !
                   fac1 = zks(num,k-1) - zks(nm,k-1)
                   fac2 = zks(num,k  ) - zks(nm,k  )
                   !
                   if ( presv(nm) == 0 ) then
                      !
                      ! free surface flow
                      !
                      gmatv(nm,k,1) = (-hks(nm ,k) - fac1) * fac
                      gmatv(nm,k,2) = (-hks(nm ,k) + fac2) * fac
                      gmatv(nm,k,3) = ( hks(num,k) - fac1) * fac
                      gmatv(nm,k,4) = ( hks(num,k) + fac2) * fac
                      gmatv(nm,k,5) = 0.
                      gmatv(nm,k,6) = 0.
                      !
                   else
                      !
                      ! pressurized flow
                      !
                      if ( presp(nm) == 0 .and. presp(num) == 1 ) then
                         !
                         ! changeover from free surface to pressurized
                         !
                         fac3 = 1. / ( hks(num,kd) + hks(num,k ) )
                         fac4 = 1. / ( hks(num,k ) + hks(num,ku) )
                         !
                         gmatv(nm,k,1) = (-hks(nm,k) - fac1) * fac
                         gmatv(nm,k,2) = (-hks(nm,k) + fac2) * fac
                         gmatv(nm,k,3) = -fac1 * hks(num,k)*fac3 * fac
                         gmatv(nm,k,4) = (2.*hks(num,k) - fac1 * hks(num,kd)*fac3 + fac2 * hks(num,ku)*fac4 ) * fac
                         gmatv(nm,k,5) = 0.
                         gmatv(nm,k,6) = fac2 * hks(num,k)*fac4 * fac
                         !
                      else if ( presp(nm) == 1 .and. presp(num) == 0 ) then
                         !
                         ! changeover from pressurized to free surface
                         !
                         fac3 = 1. / ( hks(nm,kd) + hks(nm,k ) )
                         fac4 = 1. / ( hks(nm,k ) + hks(nm,ku) )
                         !
                         gmatv(nm,k,1) = -fac1 * hks(nm,k)*fac3 * fac
                         gmatv(nm,k,2) = (-2.*hks(nm,k) - fac1 * hks(nm,kd)*fac3 + fac2 * hks(nm,ku)*fac4 ) * fac
                         gmatv(nm,k,3) = ( hks(num,k) - fac1) * fac
                         gmatv(nm,k,4) = ( hks(num,k) + fac2) * fac
                         gmatv(nm,k,5) = fac2 * hks(nm,k)*fac4 * fac
                         gmatv(nm,k,6) = 0.
                         !
                      else
                         !
                         fac1 = 0.5 * fac1 / ( guv(nm) * ( hkvm(nm,kd) + hkvm(nm,k ) ) )
                         fac2 = 0.5 * fac2 / ( guv(nm) * ( hkvm(nm,k ) + hkvm(nm,ku) ) )
                         !
                         gmatv(nm,k,1) = -fac1
                         gmatv(nm,k,2) = ( -hks(nm ,k)/guv(nm) - hkvm(nm,kd)*fac1 + hkvm(nm,ku)*fac2 ) / hkvm(nm,k)
                         gmatv(nm,k,3) = -fac1
                         gmatv(nm,k,4) = (  hks(num,k)/guv(nm) - hkvm(nm,kd)*fac1 + hkvm(nm,ku)*fac2 ) / hkvm(nm,k)
                         gmatv(nm,k,5) =  fac2
                         gmatv(nm,k,6) =  fac2
                         !
                      endif
                      !
                   endif
                   !
                else
                   !
                   gmatv(nm,k,1) = 0.
                   gmatv(nm,k,2) = 0.
                   gmatv(nm,k,3) = 0.
                   gmatv(nm,k,4) = 0.
                   gmatv(nm,k,5) = 0.
                   gmatv(nm,k,6) = 0.
                   !
                endif
                !
             enddo
             !
             if ( LMYF ) then
                !
                gmatv(nfm,k,3) = gmatv(nfm,k,3) - gmatv(nfm,k,1)
                gmatv(nfm,k,1) = 0.
                gmatv(nfm,k,4) = gmatv(nfm,k,4) - gmatv(nfm,k,2)
                gmatv(nfm,k,2) = 0.
                gmatv(nfm,k,6) = gmatv(nfm,k,6) - gmatv(nfm,k,5)
                gmatv(nfm,k,5) = 0.
                !
             endif
             !
             if ( LMYL ) then
                !
                gmatv(nlm,k,1) = gmatv(nlm,k,1) - gmatv(nlm,k,3)
                gmatv(nlm,k,3) = 0.
                gmatv(nlm,k,2) = gmatv(nlm,k,2) - gmatv(nlm,k,4)
                gmatv(nlm,k,4) = 0.
                gmatv(nlm,k,5) = gmatv(nlm,k,5) - gmatv(nlm,k,6)
                gmatv(nlm,k,6) = 0.
                !
             endif
             !
          enddo
          !
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
             if ( presv(nm) == 1 ) then
                !
                ! pressurized flow
                !
                if ( presp(nm) == 0 .and. presp(num) == 1 ) then
                   !
                   ! changeover from free surface to pressurized
                   !
                   gmatv(nm,kmax,4) = gmatv(nm,kmax,4) + 2.*gmatv(nm,kmax,6)
                   gmatv(nm,kmax,3) = gmatv(nm,kmax,3) -    gmatv(nm,kmax,6)
                   gmatv(nm,kmax,6) = 0.
                   !
                   gmatv(nm,1,4) = gmatv(nm,1,4) + 2.*gmatv(nm,1,3)
                   gmatv(nm,1,6) = gmatv(nm,1,6) -    gmatv(nm,1,3)
                   !
                else if ( presp(nm) == 1 .and. presp(num) == 0 ) then
                   !
                   ! changeover from pressurized to free surface
                   !
                   gmatv(nm,kmax,2) = gmatv(nm,kmax,2) + 2.*gmatv(nm,kmax,5)
                   gmatv(nm,kmax,1) = gmatv(nm,kmax,1) -    gmatv(nm,kmax,5)
                   gmatv(nm,kmax,5) = 0.
                   !
                   gmatv(nm,1,2) = gmatv(nm,1,2) + 2.*gmatv(nm,1,1)
                   gmatv(nm,1,5) = gmatv(nm,1,5) -    gmatv(nm,1,1)
                   !
                else
                   !
                   gmatv(nm,kmax,2) = gmatv(nm,kmax,2) + 2.*gmatv(nm,kmax,5)
                   gmatv(nm,kmax,1) = gmatv(nm,kmax,1) -    gmatv(nm,kmax,5)
                   gmatv(nm,kmax,4) = gmatv(nm,kmax,4) + 2.*gmatv(nm,kmax,6)
                   gmatv(nm,kmax,3) = gmatv(nm,kmax,3) -    gmatv(nm,kmax,6)
                   gmatv(nm,kmax,5) = 0.
                   gmatv(nm,kmax,6) = 0.
                   !
                   gmatv(nm,1,2) = gmatv(nm,1,2) + 2.*gmatv(nm,1,1)
                   gmatv(nm,1,5) = gmatv(nm,1,5) -    gmatv(nm,1,1)
                   gmatv(nm,1,4) = gmatv(nm,1,4) + 2.*gmatv(nm,1,3)
                   gmatv(nm,1,6) = gmatv(nm,1,6) -    gmatv(nm,1,3)
                   !
                endif
                !
             endif
             !
             gmatv(nm,1,1) = 0.
             gmatv(nm,1,3) = 0.
             !
          enddo
       enddo
       !
       ! to reduce the pressure Poisson equation set pressure of bottom face to that of top face for a number of layers
       !
       do l = 1, qlay
          !
          do m = mfu, ml
             do n = nf, nl
                !
                nm = kgrpnt(m,n)
                !
                gmatv(nm,kmax-l+1,1) = gmatv(nm,kmax-l+1,1) + gmatv(nm,kmax-l+1,2)
                gmatv(nm,kmax-l+1,2) = 0.
                gmatv(nm,kmax-l+1,3) = gmatv(nm,kmax-l+1,3) + gmatv(nm,kmax-l+1,4)
                gmatv(nm,kmax-l+1,4) = 0.
                !
             enddo
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
                   fac1 = 0.5 * ( zks(num,k-1) - zks(nm,k-1) ) / ( guv(nm) * ( hkvm(nm,kd) + hkvm(nm,k ) ) )
                   fac2 = 0.5 * ( zks(num,k  ) - zks(nm,k  ) ) / ( guv(nm) * ( hkvm(nm,k ) + hkvm(nm,ku) ) )
                   !
                   gmatv(nm,k,1) = -fac1
                   gmatv(nm,k,2) = ( -hks(nm ,k)/guv(nm) - hkvm(nm,kd)*fac1 + hkvm(nm,ku)*fac2 ) / hkvm(nm,k)
                   gmatv(nm,k,3) = -fac1
                   gmatv(nm,k,4) = (  hks(num,k)/guv(nm) - hkvm(nm,kd)*fac1 + hkvm(nm,ku)*fac2 ) / hkvm(nm,k)
                   gmatv(nm,k,5) =  fac2
                   gmatv(nm,k,6) =  fac2
                   !
                else
                   !
                   gmatv(nm,k,1) = 0.
                   gmatv(nm,k,2) = 0.
                   gmatv(nm,k,3) = 0.
                   gmatv(nm,k,4) = 0.
                   gmatv(nm,k,5) = 0.
                   gmatv(nm,k,6) = 0.
                   !
                endif
                !
             enddo
             !
             if ( LMYF ) then
                !
                gmatv(nfm,k,3) = gmatv(nfm,k,3) - gmatv(nfm,k,1)
                gmatv(nfm,k,1) = 0.
                gmatv(nfm,k,4) = gmatv(nfm,k,4) - gmatv(nfm,k,2)
                gmatv(nfm,k,2) = 0.
                gmatv(nfm,k,6) = gmatv(nfm,k,6) - gmatv(nfm,k,5)
                gmatv(nfm,k,5) = 0.
                !
             endif
             !
             if ( LMYL ) then
                !
                gmatv(nlm,k,1) = gmatv(nlm,k,1) - gmatv(nlm,k,3)
                gmatv(nlm,k,3) = 0.
                gmatv(nlm,k,2) = gmatv(nlm,k,2) - gmatv(nlm,k,4)
                gmatv(nlm,k,4) = 0.
                gmatv(nlm,k,5) = gmatv(nlm,k,5) - gmatv(nlm,k,6)
                gmatv(nlm,k,6) = 0.
                !
             endif
             !
          enddo
          !
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
             gmatv(nm,kmax,2) = gmatv(nm,kmax,2) + 2.*gmatv(nm,kmax,5)
             gmatv(nm,kmax,1) = gmatv(nm,kmax,1) -    gmatv(nm,kmax,5)
             gmatv(nm,kmax,4) = gmatv(nm,kmax,4) + 2.*gmatv(nm,kmax,6)
             gmatv(nm,kmax,3) = gmatv(nm,kmax,3) -    gmatv(nm,kmax,6)
             gmatv(nm,kmax,5) = 0.
             gmatv(nm,kmax,6) = 0.
             !
             if ( presv(nm) == 0 ) then
                !
                ! free surface flow
                !
                gmatv(nm,1,2) = gmatv(nm,1,2) - gmatv(nm,1,1)
                gmatv(nm,1,4) = gmatv(nm,1,4) - gmatv(nm,1,3)
                !
             else
                !
                ! pressurized flow
                !
                if ( presp(nm) == 0 .and. presp(num) == 1 ) then
                   !
                   ! changeover from free surface to pressurized
                   !
                   gmatv(nm,1,2) = gmatv(nm,1,2) -    gmatv(nm,1,1)
                   gmatv(nm,1,4) = gmatv(nm,1,4) + 2.*gmatv(nm,1,3)
                   gmatv(nm,1,6) = gmatv(nm,1,6) -    gmatv(nm,1,3)
                   !
                else if ( presp(nm) == 1 .and. presp(num) == 0 ) then
                   !
                   ! changeover from pressurized to free surface
                   !
                   gmatv(nm,1,2) = gmatv(nm,1,2) + 2.*gmatv(nm,1,1)
                   gmatv(nm,1,5) = gmatv(nm,1,5) -    gmatv(nm,1,1)
                   gmatv(nm,1,4) = gmatv(nm,1,4) -    gmatv(nm,1,3)
                   !
                else
                   !
                   gmatv(nm,1,2) = gmatv(nm,1,2) + 2.*gmatv(nm,1,1)
                   gmatv(nm,1,5) = gmatv(nm,1,5) -    gmatv(nm,1,1)
                   gmatv(nm,1,4) = gmatv(nm,1,4) + 2.*gmatv(nm,1,3)
                   gmatv(nm,1,6) = gmatv(nm,1,6) -    gmatv(nm,1,3)
                   !
                endif
                !
             endif
             !
             gmatv(nm,1,1) = 0.
             gmatv(nm,1,3) = 0.
             !
          enddo
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
             do m = mfu, ml
                do n = nsta(m), nend(m)
                   !
                   nu = n + 1
                   !
                   nm  = kgrpnt(m,n )
                   num = kgrpnt(m,nu)
                   !
                   if ( vwetp(nm) == 1 ) then
                      !
                      rhsu(nm,k) = rhsu(nm,k) - gmatv(nm,k,1)*q(nm ,kd) - gmatv(nm,k,2)*q(nm,k ) - gmatv(nm,k,3)*q(num,kd)  &
                                              - gmatv(nm,k,4)*q(num,k ) - gmatv(nm,k,5)*q(nm,ku) - gmatv(nm,k,6)*q(num,ku)
                      !
                   endif
                   !
                enddo
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
             do m = mfu, ml
                do n = nsta(m), nend(m)
                   !
                   nu = n + 1
                   !
                   nm  = kgrpnt(m,n )
                   num = kgrpnt(m,nu)
                   !
                   if ( vwetp(nm) == 1 ) then
                      !
                      rhsu(nm,k) = rhsu(nm,k) - (1.-theta3) * ( gmatv(nm,k,1)*q(nm ,kd) + gmatv(nm,k,2)*q(nm,k ) + gmatv(nm,k,3)*q(num,kd)  &
                                                              + gmatv(nm,k,4)*q(num,k ) + gmatv(nm,k,5)*q(nm,ku) + gmatv(nm,k,6)*q(num,ku) )
                      !
                   endif
                   !
                enddo
             enddo
             !
          enddo
          !
       endif
       !
    endif
    !
    ! compute baroclinic forcing at internal v-point
    !
    if ( idens /= 0 ) then
       !
       work(:,1) = 0.
       work(:,2) = 0.
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             do n = nf+1, nlast
                !
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                if ( vwetp(nm) == 1 ) then
                   !
                   rhsu(nm,k) = rhsu(nm,k) - 0.5 * grav * hkvm(nm,k) * (rho(num,k) - rho(nm,k)) / ( rhow*guv(nm) )
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       do k = 2, kmax
          !
          do m = mfu, ml
             do n = nf+1, nlast
                !
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                if ( vwetp(nm) == 1 ) then
                   !
                   rhou = rhow + 0.5 * (rho(num,k-1) + rho(nm,k-1))
                   !
                   work(nm,1) = work(nm,1) + rhou * (hks(num,k-1) - hks(nm,k-1)) + hkvm(nm,k-1) * (rho(num,k-1) - rho(nm,k-1))
                   !
                   work(nm,2) = work(nm,2) + hks(num,k-1) - hks(nm,k-1)
                   !
                   rhou = rhow + 0.5 * (rho(num,k) + rho(nm,k))
                   !
                   rhsu(nm,k) = rhsu(nm,k) - grav * ( work(nm,1) - rhou*work(nm,2) ) / ( rhow*guv(nm) )
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
    endif
    !
    ! compute atmospheric pressure gradient at internal v-point
    !
    if ( svwp ) then
       !
       do m = mfu, ml
          do n = nf+1, nlast
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( vwetp(nm) == 1 ) then
                !
                rhsu(nm,:) = rhsu(nm,:) - (patm(num) - patm(nm)) / ( rhow*guv(nm) )
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute horizontal viscosity term at internal v-point
    !
    if ( ihvisc == 1 .and. hvisc > 0. ) then
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             !
             do n = nf+1, nlast
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
                if ( vwetp(nm) == 1 ) then
                   !
                   if ( wetv(ndm) * wetv(num) == 1 ) then
                      !
                      fac1 = hks(num,k) * ( gvv(num)*v0(num,k) - gvv(nm)*v0(nm,k) )/gsqs(num) - hks(nm,k) * ( gvv(nm)*v0(nm,k) - gvv(ndm)*v0(ndm,k) )/gsqs(nm)
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
                      hd  = .5 * (hkvm(nmu,k) + hkvm(nm ,k))
                      hdl = .5 * (hkvm(nm ,k) + hkvm(nmd,k))
                      !
                      fac2 = hd * ( guv(nmu)*v0(nmu,k) - guv(nm)*v0(nm,k) )/gsqsd - hdl * ( guv(nm)*v0(nm,k) - guv(nmd)*v0(nmd,k) )/gsqsdl
                      !
                   else
                      !
                      fac2 = 0.
                      !
                   endif
                   !
                   rhsu(nm,k) = rhsu(nm,k) + hvisc * ( fac1/guv(nm) + fac2/gvv(nm) ) / hkvm(nm,k)
                   !
                endif
                !
             enddo
             !
          enddo
          !
       enddo
       !
    else if ( ihvisc > 1 ) then
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             !
             do n = nf+1, nlast
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
                if ( vwetp(nm) == 1 ) then
                   !
                   lmask = .true.
                   !
                   hd  = .5 * (hkvm(nmu,k) + hkvm(nm ,k))
                   hdl = .5 * (hkvm(nm ,k) + hkvm(nmd,k))
                   !
                   ! first, compute the components of the Reynolds stress tensor in different points of consideration
                   !
                   if ( wetv(ndm) * wetv(num) == 1 ) then
                      !
                      ryys  = ( vnu2d(nm ) + vnu2d(nmd ) + vnu2d(ndm) + vnu2d(ndmd) ) * ( v0(nm ,k) - v0(ndm,k) ) / ( guu(nm ) + guu(nmd ) )
                      !
                      ryysu = ( vnu2d(num) + vnu2d(numd) + vnu2d(nm ) + vnu2d(nmd ) ) * ( v0(num,k) - v0(nm ,k) ) / ( guu(num) + guu(numd) )
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
                      rxyd  = 2. * vnu2d(nm ) * ( ( u0(num ,k) - u0(nm ,k) )/( guu(num ) + guu(nm ) ) + ( v0(nmu,k) - v0(nm ,k) )/( gvv(nmu) + gvv(nm ) ) )
                      !
                      rxydd = 2. * vnu2d(nmd) * ( ( u0(numd,k) - u0(nmd,k) )/( guu(numd) + guu(nmd) ) + ( v0(nm ,k) - v0(nmd,k) )/( gvv(nm ) + gvv(nmd) ) )
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
                      rxxu  = 0.50 * ( vnu2d(nm) + vnu2d(nmd) ) * ( u0(num,k) + u0(nm,k) - u0(nmd,k) -u0(numd,k) ) / gvv(nm)
                      !
                      rxyu  = 0.25 * ( vnu2d(nm) + vnu2d(nmd) ) * ( ( v0(nmu,k) - v0(nmd,k) )/gvv(nm) + ( u0(num,k) + u0(numd,k) - u0(nm,k) - u0(nmd,k) )/guv(nm) )
                      !
                      ryyu  = 0.50 * ( vnu2d(nm) + vnu2d(nmd) ) * ( v0(num,k) - v0(ndm,k) ) / guv(nm)
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
                   rhsu(nm,k) = rhsu(nm,k) + ( hd*rxyd - hdl*rxydd ) / (hkvm(nm,k)*gvv(nm)) + ( hks(num,k)*ryysu - hks(nm,k)*ryys ) / (hkvm(nm,k)*guv(nm))
                   !
                   ! finally, add the curvature terms associated with those stresses
                   !
                   rhsu(nm,k) = rhsu(nm,k) + ( rxyu * ( guu(nm) + guu(num) - guu(nmd) - guu(numd) ) + 0.5 * (ryyu - rxxu) * ( gvv(num) - gvv(ndm) ) ) / gsqsv(nm)
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
    ! add contribution due to full 3D turbulence
    !
    if ( iturb > 1 ) then
       !
       ! compute divergence of Reynolds stress tensor at internal v-point
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             do n = nf+1, nlast
                !
                nu = n + 1
                md = m - 1
                !
                nm  = kgrpnt(m ,n )
                num = kgrpnt(m ,nu)
                nmd = kgrpnt(md,n )
                !
                if ( vwetp(nm) == 1 ) then
                   !
                   rhsu(nm,k) = rhsu(nm,k) + ( rsvu(nm ,k  ) - rsvu(nmd,k) ) / gvv(nm)     &
                                           + ( rsvv(num,k  ) - rsvv(nm ,k) ) / guv(nm)     &
                                           + ( rsvw(nm ,k-1) - rsvw(nm ,k) ) / hkvm(nm,k)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! impose wall shear stress at closed boundaries
       !
       if ( .not.lreptx ) then
          !
          do n = nf+1, nlast
             !
             nmfu = kgrpnt(mfu,n)
             nml  = kgrpnt(ml ,n)
             !
             if ( ibl(n) == 1 .and. LMXF ) then
                !
                amatu(nmfu,:,1) = amatu(nmfu,:,1) + logwll(n,:) / gvv(nmfu)
                !
             endif
             !
             if ( ibr(n) == 1 .and. LMXL ) then
                !
                amatu(nml,:,1) = amatu(nml,:,1) + logwlr(n,:) / gvv(nml)
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
    ! compute explicit part of wind stress term at internal v-point (top layer only), if appropriate
    !
    if ( iwind /= 0 ) then
       !
       do m = mfu, ml
          do n = nf+1, nlast
             !
             nm = kgrpnt(m,n)
             !
             if ( vwetp(nm) == 1 ) then
                !
                rhsu(nm,1) = rhsu(nm,1) + windv(nm)/max(1.e-3,hkvm(nm,1))
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute implicit part of wind stress term at internal v-point (top layer only), if appropriate
    !
    if ( relwnd ) then
       !
       do m = mfu, ml
          do n = nf+1, nlast
             !
             nm = kgrpnt(m,n)
             !
             if ( vwetp(nm) == 1 ) then
                !
                amatu(nm,1,1) = amatu(nm,1,1) + cwndv(nm) / max(1.e-3,hkvm(nm,1))
                !
             endif
             !
          enddo
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
       do m = mfu, ml
          do n = nsta(m), nend(m)
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( vwetp(nm) == 1 ) then
                !
                amatu(nm,kmax,1) = amatu(nm,kmax,1) + 0.5 * ( logfrc(nm,1) + logfrc(num,1) ) / hkvm(nm,kmax)
                !
             endif
             !
          enddo
       enddo
       !
    else if ( irough == 11 ) then
       !
       ! linear bottom friction
       !
       do m = mfu, ml
          do n = nsta(m), nend(m)
             !
             nm = kgrpnt(m,n)
             !
             if ( vwetp(nm) == 1 ) then
                !
                amatu(nm,kmax,1) = amatu(nm,kmax,1) + cfricv(nm) / hkvm(nm,kmax)
                !
             endif
             !
          enddo
       enddo
       !
    else if ( irough /= 0 ) then
       !
       do m = mfu, ml
          do n = nsta(m), nend(m)
             !
             nu = n + 1
             md = m - 1
             !
             nm   = kgrpnt(m ,n )
             num  = kgrpnt(m ,nu)
             nmd  = kgrpnt(md,n )
             numd = kgrpnt(md,nu)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( num  == 1 ) num  = nm
             if ( nmd  == 1 ) nmd  = nm
             if ( numd == 1 ) numd = nmd
             !
             u = 0.25 * ( u0(nm,kmax) + u0(num,kmax) + u0(nmd,kmax) + u0(numd,kmax) )
             !
             utot = sqrt( u*u + v0(nm,kmax)*v0(nm,kmax) )
             !
             if ( vwetp(nm) == 1 .and. utot > 1.e-8 ) then
                !
                u = 0.25 * ( udep(nm) + udep(num) + udep(nmd) + udep(numd) )
                !
                fac = u*u + vdep(nm)*vdep(nm)
                !
                amatu(nm,kmax,1) = amatu(nm,kmax,1) + cfricv(nm) * fac / ( utot * hkvm(nm,kmax) )
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
    ! compute Coriolis force, if appropriate
    !
    if ( coriolis ) then
       !
       epsab2 = pcor(2)
       !
       do m = mfu, ml
          do n = nsta(m), nend(m)
             !
             nu = n + 1
             md = m - 1
             !
             nm   = kgrpnt(m ,n )
             num  = kgrpnt(m ,nu)
             nmd  = kgrpnt(md,n )
             numd = kgrpnt(md,nu)
             !
             if ( vwetp(nm) == 1 .and. fcor(nm,2) /= 0. ) then
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( num  == 1 ) num  = nm
                if ( nmd  == 1 ) nmd  = nm
                if ( numd == 1 ) numd = nmd
                !
                ! energy-conserved weighting factors according to Espelid et al. (2000)
                !
                fc(1) = sqrt(abs(fcor(nm  ,1))*dpu(nm  ))
                fc(2) = sqrt(abs(fcor(num ,1))*dpu(num ))
                fc(3) = sqrt(abs(fcor(nmd ,1))*dpu(nmd ))
                fc(4) = sqrt(abs(fcor(numd,1))*dpu(numd))
                fc(5) = 1. / (4.*sqrt(abs(fcor(nm,2))*dpv(nm)))
                !
                do k = 1, kmax
                   !
                   ! at time level n-1
                   !
                   ! --- modified quasi second order Adams-Bashforth scheme which is conditionally stable for the inviscid case (see Marshall et al., 1997)
                   !
                   rhsu(nm,k) = rhsu(nm,k) + (0.5+epsab2) * cfv(nm,k)
                   !
                   ! at time level n
                   !
                   ! --- interpolate u-velocity at v-point using energy-conserved weighting factors
                   !
                   u = fc(1) * u0(nm,k) + fc(2) * u0(num,k) + fc(3) * u0(nmd,k) + fc(4) * u0(numd,k)
                   u = fc(5) * u
                   !
                   cfv(nm,k) = fcor(nm,2) * u
                   !
                   rhsu(nm,k) = rhsu(nm,k) - (1.5+epsab2) * cfv(nm,k)
                   !
                enddo
                !
             endif
             !
          enddo
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
       do m = mfu, ml
          do n = nsta(m), nend(m)
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( vwetp(nm) == 1 ) then
                !
                ! advection term
                !
                if ( propsc == 3 .and. kappa == 1. ) then
                   !
                   ! central differences
                   !
                   ctrkt = 0.5 * kwd * ( wom(nm,k-1) + wom(num,k-1) ) / ( hkvm(nm,kd) + hkvm(nm,k ) )
                   ctrkb = 0.5 * kwu * ( wom(nm,k  ) + wom(num,k  ) ) / ( hkvm(nm,k ) + hkvm(nm,ku) )
                   !
                   amatu(nm,k,1) = amatu(nm,k,1) - ctrkt + ctrkb
                   amatu(nm,k,2) =  ctrkt
                   amatu(nm,k,3) = -ctrkb
                   !
                else
                   !
                   ! first order upwind scheme
                   !
                   ctrkt = 0.5 * kwd * ( wom(nm,k-1) + wom(num,k-1) ) / hkvm(nm,k)
                   ctrkb = 0.5 * kwu * ( wom(nm,k  ) + wom(num,k  ) ) / hkvm(nm,k)
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
                         grad1 = v0(nm,kd) - v0(nm,k )
                         grad2 = v0(nm,k ) - v0(nm,ku)
                         !
                         rhsu(nm,k) = rhsu(nm,k) - 0.5 * ctrkt * fluxlim(grad1,grad2)
                         !
                      else if ( ctrkt < 0. ) then
                         !
                         grad1 = v0(nm,kd ) - v0(nm,k )
                         grad2 = v0(nm,kdd) - v0(nm,kd)
                         !
                         rhsu(nm,k) = rhsu(nm,k) + 0.5 * ctrkt * fluxlim(grad1,grad2)
                         !
                      endif
                      !
                      if ( ctrkb > 0. ) then
                         !
                         grad1 = v0(nm,k ) - v0(nm,ku )
                         grad2 = v0(nm,ku) - v0(nm,kuu)
                         !
                         rhsu(nm,k) = rhsu(nm,k) + 0.5 * ctrkb * fluxlim(grad1,grad2)
                         !
                      else if ( ctrkb < 0. ) then
                         !
                         grad1 = v0(nm,k ) - v0(nm,ku)
                         grad2 = v0(nm,kd) - v0(nm,k )
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
                   ctrkt = ( vnu3d(nm,k-1) + vnu3d(num,k-1) ) / ( hkvm(nm,k)*( hkvm(nm,kd) + hkvm(nm,k ) ) )
                   ctrkb = ( vnu3d(nm,k  ) + vnu3d(num,k  ) ) / ( hkvm(nm,k)*( hkvm(nm,k ) + hkvm(nm,ku) ) )
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
       enddo
       !
    enddo
    !
    ! solve the v-momentum equation
    !
    do m = mfu, ml
       do n = nsta(m), nend(m)
          !
          nm = kgrpnt(m,n)
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
          v1(nm,kmax) = rhsu(nm,kmax)
          do k = kmax-1, 1, -1
             v1(nm,k) = rhsu(nm,k) - amatu(nm,k,3)*v1(nm,k+1)
          enddo
          !
       enddo
    enddo
    !
    ! re-update the solution in case of thetau <> 1
    !
    if ( thetau /= 1. ) then
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             do n = nsta(m), nend(m)
                !
                nm = kgrpnt(m,n)
                !
                if ( vwetp(nm) == 1 ) then
                   !
                   v1(nm,k) = ( v1(nm,k) - (1.-thetau) * v0(nm,k) ) / thetau
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
    endif
    !
    ! exchange v-velocities with neighbouring subdomains
    !
    call SWEXCHG ( v1, kgrpnt, 1, kmax )
    if (STPNOW()) return
    !
    ! synchronize v-velocities at appropriate boundaries in case of repeating grid
    !
    call periodic ( v1, kgrpnt, 1, kmax )
    !
    ! compute higher order correction to the advective velocity in depth point (if appropriate)
    !
    propsc = nint(pnums(6))
    !
    if ( propsc /= 1 ) then
       !
       kappa  = pnums(7)
       mbound = pnums(8)
       phieby = pnums(9)
       !
       if ( lrepty ) then
          nlast = nl            ! last  internal v-point in case of repeating grid
       else
          nlast = nl - 1        ! last  internal v-point
       endif
       !
       do k = 1, kmax
          !
          do m = mf, ml
             do n = nf+1, nlast
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
                if ( ( vdep(nm) > vdep(ndm) .or. stricthead ) .and. .not.strictmom ) then
                   !
                   fac = 0.5 * ( u0(nm,k) + u0(num,k) )
                   !
                else
                   !
                   fac = qxm(nm,k)
                   !
                endif
                !
                if ( .not. fac < 0. ) then
                   !
                   grad1 = v1(nmu,k) - v0(nm ,k)
                   grad2 = v1(nm ,k) - v0(nmd,k)
                   !
                   ua2(nm,k) = 0.5 * fluxlim(grad1,grad2)
                   !
                else
                   !
                   grad1 = v1(nm ,k) - v0(nmu ,k)
                   grad2 = v1(nmu,k) - v0(nmuu,k)
                   !
                   ua2(nm,k) = 0.5 * fluxlim(grad1,grad2)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
    endif
    !
    ! compute higher order correction to the advective velocity in wl-point (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do k = 1, kmax
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
                if ( ( vdep(nm) > vdep(ndm) .or. stricthead ) .and. .not.strictmom ) then
                   !
                   fac = 0.5 * ( v0(nm,k) + v0(ndm,k) )
                   !
                else
                   !
                   fac = qym(nm,k)
                   !
                endif
                !
                if ( .not. fac < 0. ) then
                   !
                   grad1 = v1(nm ,k) - v0(ndm ,k)
                   grad2 = v1(ndm,k) - v0(nddm,k)
                   !
                   up(nm,k) = v0(ndm,k)
                   ua(nm,k) = 0.5 * fluxlim(grad1,grad2)
                   !
                else
                   !
                   grad1 = v1(ndm,k) - v0(nm ,k)
                   grad2 = v1(nm ,k) - v0(num,k)
                   !
                   up(nm,k) = v0(nm,k)
                   ua(nm,k) = 0.5 * fluxlim(grad1,grad2)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       call periodic ( ua, kgrpnt, 1, kmax )
       call periodic ( up, kgrpnt, 1, kmax )
       !
    endif
    !
    ! correct v-velocity associated with higher order approximation of advection term in x-direction (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             do n = nf+1, nlast
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
                if ( ( vdep(nm) > vdep(ndm) .or. stricthead ) .and. .not.strictmom ) then
                   !
                   !fac1 = 0.5 * ( u0(nm,k) + u0(num,k) )
                   !fac2 = fac1
                   !jumps back to old implementation (020117)
                   fac1 = 0.
                   fac2 = 0.
                   !
                else if ( hkvmn(nm,k) > 0. ) then
                   !
                   fac1 = qxm(nm ,k) / hkvmn(nm,k)
                   fac2 = qxm(nmd,k) / hkvmn(nm,k)
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
                if ( vwetp(nm) == 1 .and. presv(nm) /= 1 ) then
                   !
                   v1(nm,k) = v1(nm,k) - dt * ( fac1 * ua2(nm,k) - fac2 * ua2(nmd,k) ) / gvv(nm)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
    endif
    !
    ! correct v-velocity associated with higher order approximation of advection term in y-direction (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             do n = nf+1, nlast
                !
                nd = n - 1
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                ndm = kgrpnt(m,nd)
                num = kgrpnt(m,nu)
                !
                if ( ( vdep(nm) > vdep(ndm) .or. stricthead ) .and. .not.strictmom ) then
                   !
                   !fac1 = up(num,k) + 0.5*ua(num,k)
                   !fac2 = up(nm ,k) + 0.5*ua(nm ,k)
                   !jumps back to old implementation (020117)
                   fac1 = 0.
                   fac2 = 0.
                   !
                else if ( hkvmn(nm,k) > 0. ) then
                   !
                   fac1 = qym(num,k) / hkvmn(nm,k)
                   fac2 = qym(nm ,k) / hkvmn(nm,k)
                   !
                else
                   !
                   fac1 = 0.
                   fac2 = 0.
                   !
                endif
                !
                if ( vwetp(nm) == 1 .and. presv(nm) /= 1 ) then
                   !
                   v1(nm,k) = v1(nm,k) - dt * ( fac1 * ua(num,k)  - fac2 * ua(nm,k) ) / guv(nm)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
    endif
    !
    ! exchange v-velocities with neighbouring subdomains (if appropriate)
    !
    if ( propsc /= 1 ) then
       !
       call SWEXCHG ( v1, kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
    endif
    !
    ! synchronize v-velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    if ( propsc /= 1 ) call periodic ( v1, kgrpnt, 1, kmax )
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
    jj    = 0
    resm  = epslin + 1.
    !
    ! start iteration process to obtain 2nd order accuracy in pressure projection method
    !
 10 if ( resm > epslin .and. jj < maxit ) then
       !
       jj = jj + 1
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
                if ( ( m /= mfu .or. slimp(n) .or. .not.LMXF ) .and. nmd /= 1 ) then
                   !
                   amat(nm,2) = -fac * tetau(nmd) * guu(nmd) * hu(nmd) / gvu(nmd)
                   !
                else
                   !
                   amat(nm,2) = 0.
                   !
                endif
                !
                if ( ( m /= ml .or. srimp(n) .or. .not.LMXL ) .and. nmu /= 1 ) then
                   !
                   amat(nm,3) = -fac * tetau(nm) * guu(nm) * hu(nm) / gvu(nm)
                   !
                else
                   !
                   amat(nm,3) = 0.
                   !
                endif
                !
                if ( ( n /= nfu .or. sbimp(m) .or. .not.LMYF ) .and. ndm /= 1 ) then
                   !
                   amat(nm,4) = -fac * tetav(ndm) * gvv(ndm) * hv(ndm) / guv(ndm)
                   !
                else
                   !
                   amat(nm,4) = 0.
                   !
                endif
                !
                if ( ( n /= nl .or. stimp(m) .or. .not.LMYL ) .and. num /= 1 ) then
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
                !
                fac1 = 0.
                fac2 = 0.
                !
                do k = 1, kmax
                   !
                   fac1 = fac1 + hku(nmd,k)* ( tetau(nmd)*ui(nmd,k) + (1.-tetau(nmd))*u0(nmd,k) )
                   fac2 = fac2 + hku(nm ,k)* ( tetau(nm )*ui(nm ,k) + (1.-tetau(nm ))*u0(nm ,k) )
                   !
                enddo
                !
                rhs(nm) = dt * ( guu(nmd) * fac1 - guu(nm) * fac2 )
                !
                fac1 = 0.
                fac2 = 0.
                !
                do k = 1, kmax
                   !
                   fac1 = fac1 + hkv(ndm,k)* ( tetav(ndm)*vi(ndm,k) + (1.-tetav(ndm))*v0(ndm,k) )
                   fac2 = fac2 + hkv(nm ,k)* ( tetav(nm )*vi(nm ,k) + (1.-tetav(nm ))*v0(nm ,k) )
                   !
                enddo
                !
                rhs(nm) = rhs(nm) + dt * ( gvv(ndm) * fac1 - gvv(nm) * fac2 )
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
                      !
                      fac1 = 0.
                      fac2 = 0.
                      !
                      do k = 1, kmax
                         !
                         fac1 = fac1 + hku(nmd,k) * dqgrdu(nmd,k)
                         fac2 = fac2 + hku(nm ,k) * dqgrdu(nm ,k)
                         !
                      enddo
                      !
                      rhs(nm) = rhs(nm) + dt * dt * theta3 * ( tetau(nm) * guu(nm) * fac2 - tetau(nmd) * guu(nmd) * fac1 )
                      !
                      fac1 = 0.
                      fac2 = 0.
                      !
                      do k = 1, kmax
                         !
                         fac1 = fac1 + hkv(ndm,k) * dqgrdv(ndm,k)
                         fac2 = fac2 + hkv(nm ,k) * dqgrdv(nm ,k)
                         !
                      enddo
                      !
                      rhs(nm) = rhs(nm) + dt * dt * theta3 * ( tetav(nm) * gvv(nm) * fac2 - tetav(ndm) * gvv(ndm) * fac1 )
                      !
                   endif
                   !
                enddo
             enddo
             !
          endif
          !
          do n = nfu, nl
             !
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             nml  = kgrpnt(ml ,n)
             nmlu = kgrpnt(mlu,n)
             !
             if ( slimp(n) ) then
                !
                ds  (nmf   ) = s1(nmf) - s0(nmf)
                rhs (nmfu  ) = rhs(nmfu) - amat(nmfu,2)*ds(nmf)
                amat(nmfu,2) = 0.
                !
             endif
             !
             if ( srimp(n) ) then
                !
                ds  (nmlu ) = s1(nmlu) - s0(nmlu)
                rhs (nml  ) = rhs(nml) - amat(nml,3)*ds(nmlu)
                amat(nml,3) = 0.
                !
             endif
             !
          enddo
          !
          do m = mfu, ml
             !
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             nlm  = kgrpnt(m,nl )
             nlum = kgrpnt(m,nlu)
             !
             if ( sbimp(m) ) then
                !
                ds  (nfm   ) = s1(nfm) - s0(nfm)
                rhs (nfum  ) = rhs(nfum) - amat(nfum,4)*ds(nfm)
                amat(nfum,4) = 0.
                !
             endif
             !
             if ( stimp(m) ) then
                !
                ds  (nlum ) = s1(nlum) - s0(nlum)
                rhs (nlm  ) = rhs(nlm) - amat(nlm,5)*ds(nlum)
                amat(nlm,5) = 0.
                !
             endif
             !
          enddo
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
                if ( ( m /= mfu .or. slimp(n) .or. .not.LMXF ) .and. nmd /= 1 ) then
                   !
                   amatn(nm,2) = -fac * tetau(nmd) * guu(nmd) * hu(nmd) / gvu(nmd)
                   !
                else
                   !
                   amatn(nm,2) = 0.
                   !
                endif
                !
                if ( ( m /= ml .or. srimp(n) .or. .not.LMXL ) .and. nmu /= 1 ) then
                   !
                   amatn(nm,3) = -fac * tetau(nm) * guu(nm) * hu(nm) / gvu(nm)
                   !
                else
                   !
                   amatn(nm,3) = 0.
                   !
                endif
                !
                if ( ( n /= nfu .or. sbimp(m) .or. .not.LMYF ) .and. ndm /= 1 ) then
                   !
                   amatn(nm,4) = -fac * tetav(ndm) * gvv(ndm) * hv(ndm) / guv(ndm)
                   !
                else
                   !
                   amatn(nm,4) = 0.
                   !
                endif
                !
                if ( ( n /= nl .or. stimp(m) .or. .not.LMYL ) .and. num /= 1 ) then
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
                !
                fac1 = 0.
                fac2 = 0.
                !
                do k = 1, kmax
                   !
                   fac1 = fac1 + hku(nmd,k)* ( tetau(nmd)*ui(nmd,k) + (1.-tetau(nmd))*u0(nmd,k) )
                   fac2 = fac2 + hku(nm ,k)* ( tetau(nm )*ui(nm ,k) + (1.-tetau(nm ))*u0(nm ,k) )
                   !
                enddo
                !
                rhsn(nm) = dt * ( guu(nmd) * fac1 - guu(nm) * fac2 )
                !
                fac1 = 0.
                fac2 = 0.
                !
                do k = 1, kmax
                   !
                   fac1 = fac1 + hkv(ndm,k)* ( tetav(ndm)*vi(ndm,k) + (1.-tetav(ndm))*v0(ndm,k) )
                   fac2 = fac2 + hkv(nm ,k)* ( tetav(nm )*vi(nm ,k) + (1.-tetav(nm ))*v0(nm ,k) )
                   !
                enddo
                !
                rhsn(nm) = rhsn(nm) + dt * ( gvv(ndm) * fac1 - gvv(nm) * fac2 )
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
                      !
                      fac1 = 0.
                      fac2 = 0.
                      !
                      do k = 1, kmax
                         !
                         fac1 = fac1 + hku(nmd,k) * dqgrdu(nmd,k)
                         fac2 = fac2 + hku(nm ,k) * dqgrdu(nm ,k)
                         !
                      enddo
                      !
                      rhsn(nm) = rhsn(nm) + dt * dt * theta3 * ( tetau(nm) * guu(nm) * fac2 - tetau(nmd) * guu(nmd) * fac1 )
                      !
                      fac1 = 0.
                      fac2 = 0.
                      !
                      do k = 1, kmax
                         !
                         fac1 = fac1 + hkv(ndm,k) * dqgrdv(ndm,k)
                         fac2 = fac2 + hkv(nm ,k) * dqgrdv(nm ,k)
                         !
                      enddo
                      !
                      rhsn(nm) = rhsn(nm) + dt * dt * theta3 * ( tetav(nm) * gvv(nm) * fac2 - tetav(ndm) * gvv(ndm) * fac1 )
                      rhsn(nm) = rhsn(nm) / gsqs(nm)
                      !
                   endif
                   !
                enddo
             enddo
             !
          endif
          !
          do n = nfu, nl
             !
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             nml  = kgrpnt(ml ,n)
             nmlu = kgrpnt(mlu,n)
             !
             if ( slimp(n) ) then
                !
                ds   (nmf   ) = s1(nmf) - s0(nmf)
                rhsn (nmfu  ) = rhsn(nmfu) - amatn(nmfu,2)*ds(nmf)
                amatn(nmfu,2) = 0.
                !
             endif
             !
             if ( srimp(n) ) then
                !
                ds   (nmlu ) = s1(nmlu) - s0(nmlu)
                rhsn (nml  ) = rhsn(nml) - amatn(nml,3)*ds(nmlu)
                amatn(nml,3) = 0.
                !
             endif
             !
          enddo
          !
          do m = mfu, ml
             !
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             nlm  = kgrpnt(m,nl )
             nlum = kgrpnt(m,nlu)
             !
             if ( sbimp(m) ) then
                !
                ds   (nfm   ) = s1(nfm) - s0(nfm)
                rhsn (nfum  ) = rhsn(nfum) - amatn(nfum,4)*ds(nfm)
                amatn(nfum,4) = 0.
                !
             endif
             !
             if ( stimp(m) ) then
                !
                ds   (nlum ) = s1(nlum) - s0(nlum)
                rhsn (nlm  ) = rhsn(nlm) - amatn(nlm,5)*ds(nlum)
                amatn(nlm,5) = 0.
                !
             endif
             !
          enddo
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
          write (PRINTF,'(a,i4,a,e12.6,a,i5,a,i5,a)') ' ++ pressure projection: iter = ',jj,' res = ',resm,' in (m,n)=(',minf+MXF-2,',',ninf+MYF-2,')'
          !
       endif
       !
       ! correct horizontal velocities
       !
       do k = 1, kmax
          !
          do n = nfu, nl
             do m = msta(n), mend(n)
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                if ( wetu(nm) == 1 .and. nmu /= 1 ) then
                   !
                   u1(nm,k) = ui(nm,k) - grav * dt * ( teta2(nmu) * ds(nmu) - teta2(nm) * ds(nm) ) / gvu(nm)
                   !
                endif
                !
             enddo
          enddo
          !
          do m = mfu, ml
             do n = nsta(m), nend(m)
                !
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                if ( wetv(nm) == 1 .and. num /= 1 ) then
                   !
                   v1(nm,k) = vi(nm,k) - grav * dt * ( teta2(num) * ds(num) - teta2(nm) * ds(nm) ) / guv(nm)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! exchange velocities with neighbouring subdomains
       !
       call SWEXCHG ( u1, kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
       call SWEXCHG ( v1, kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
       ! synchronize velocities at appropriate boundaries in case of repeating grid
       !
       call periodic ( u1, kgrpnt, 1, kmax )
       call periodic ( v1, kgrpnt, 1, kmax )
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
                qf = 0.
                do k = 1, kmax
                   qf = qf + hku(nm,k) * guu(nm) * ( tetau(nm)*u1(nm,k) + (1.-tetau(nm))*u0(nm,k) )
                enddo
                !
                if ( qf > 0. ) sumqf = sumqf + qf
                !
                qf = 0.
                do k = 1, kmax
                   qf = qf + hku(nmd,k) * guu(nmd) * ( tetau(nmd)*u1(nmd,k) + (1.-tetau(nmd))*u0(nmd,k) )
                enddo
                !
                if ( qf < 0. ) sumqf = sumqf - qf
                !
                qf = 0.
                do k = 1, kmax
                   qf = qf + hkv(nm,k) * gvv(nm) * ( tetav(nm)*v1(nm,k) + (1.-tetav(nm))*v0(nm,k) )
                enddo
                !
                if ( qf > 0. ) sumqf = sumqf + qf
                !
                qf = 0.
                do k = 1, kmax
                   qf = qf + hkv(ndm,k) * gvv(ndm) * ( tetav(ndm)*v1(ndm,k) + (1.-tetav(ndm))*v0(ndm,k) )
                enddo
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
       ! build the w-momentum equation
       !
       if ( ihydro == 1 .or. ihydro == 2 ) then
          !
          ! initialize system of equations in dry points
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
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
          enddo
          !
          ! bottom:
          !
          ! the kinematic condition is imposed
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
                   amatw(nm,kmax,1) = 1.
                   amatw(nm,kmax,2) = 0.
                   amatw(nm,kmax,3) = 0.
                   w1   (nm,kmax  ) = 0.25 * ( 3.*(u1(nm,kmax)+u1(nmd,kmax)) - u1(nm,kmax-1) - u1(nmd,kmax-1) ) * ( zku(nm,kmax) - zku(nmd,kmax) ) / dxl +  &
                                      0.25 * ( 3.*(v1(nm,kmax)+v1(ndm,kmax)) - v1(nm,kmax-1) - v1(ndm,kmax-1) ) * ( zkv(nm,kmax) - zkv(ndm,kmax) ) / dyl
                   rhsw (nm,kmax  ) = w1(nm,kmax)
                   !
                endif
                !
             enddo
          enddo
          !
          ! free surface and interior part:
          !
          ! compute the time derivative
          !
          do k = 0, kmax-1
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wets(nm) == 1 ) then
                      !
                      amatw(nm,k,1) = 1. / (dt*thetaw)
                      rhsw (nm,k  ) = w0(nm,k) / (dt*thetaw)
                      !
                   endif
                   !
                enddo
             enddo
             !
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
             do k = 0, kmax-1
                !
                kd = max(k,1)
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
                         fac = 0.5 * gsqs(nm) * ( hks(nm,kd) + hks(nm,k+1) )
                         !
                         if ( uwetp(nmd) == 1 ) then
                            fac1 = guu(nmd) * ( qx(nmd,k+1)*hku(nmd,kd) + qx(nmd,kd)*hku(nmd,k+1) ) / ( fac*( hku(nmd,kd) + hku(nmd,k+1) ) )
                         else
                            fac1 = 0.
                         endif
                         if ( uwetp(nm) == 1 ) then
                            fac2 = guu(nm ) * ( qx(nm ,k+1)*hku(nm ,kd) + qx(nm ,kd)*hku(nm ,k+1) ) / ( fac*( hku(nm ,kd) + hku(nm ,k+1) ) )
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
                         ! advection term in y-direction
                         !
                         if ( vwetp(ndm) == 1 ) then
                            fac1 = gvv(ndm) * ( qy(ndm,k+1)*hkv(ndm,kd) + qy(ndm,kd)*hkv(ndm,k+1) ) / ( fac*( hkv(ndm,kd) + hkv(ndm,k+1) ) )
                         else
                            fac1 = 0.
                         endif
                         if ( vwetp(nm) == 1 ) then
                            fac2 = gvv(nm ) * ( qy(nm ,k+1)*hkv(nm ,kd) + qy(nm ,kd)*hkv(nm ,k+1) ) / ( fac*( hkv(nm ,kd) + hkv(nm ,k+1) ) )
                         else
                            fac2 = 0.
                         endif
                         !
                         if ( propsc == 3 .and. kappa == 1. ) then
                            !
                            ! central differences
                            !
                            rhsw(nm,k) = rhsw(nm,k) - 0.5 * ( fac1 * (w0(nm,k) - w0(ndm,k)) + fac2 * (w0(num,k) - w0(nm,k)) )
                            !
                         else
                            !
                            ! first order upwind scheme
                            !
                            rhsw(nm,k) = rhsw(nm,k) - ( max(fac1,0.) * (w0(nm,k) - w0(ndm,k)) + min(fac2,0.) * (w0(num,k) - w0(nm,k)) )
                            !
                            ! add higher order (flux-limited) correction, if appropriate
                            !
                            if ( propsc /= 1 ) then
                               !
                               if ( fac2 > 0. ) then
                                  !
                                  grad1 = w0(num,k) - w0(nm ,k)
                                  grad2 = w0(nm ,k) - w0(ndm,k)
                                  !
                                  rhsw(nm,k) = rhsw(nm,k) - 0.5 * fac2 * fluxlim(grad1,grad2)
                                  !
                               else if ( fac2 < 0. ) then
                                  !
                                  grad1 = w0(num ,k) - w0(nm ,k)
                                  grad2 = w0(nuum,k) - w0(num,k)
                                  !
                                  rhsw(nm,k) = rhsw(nm,k) + 0.5 * fac2 * fluxlim(grad1,grad2)
                                  !
                               endif
                               !
                               if ( fac1 > 0. ) then
                                  !
                                  grad1 = w0(nm ,k) - w0(ndm ,k)
                                  grad2 = w0(ndm,k) - w0(nddm,k)
                                  !
                                  rhsw(nm,k) = rhsw(nm,k) + 0.5 * fac1 * fluxlim(grad1,grad2)
                                  !
                               else if ( fac1 < 0. ) then
                                  !
                                  grad1 = w0(nm ,k) - w0(ndm,k)
                                  grad2 = w0(num,k) - w0(nm ,k)
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
                         ! viscosity term in x-direction
                         !
                         if ( wets(nm) * wets(nmd) * wets(nmu) == 1 ) then
                            !
                            fac = 0.25 * ( gvv(nm) + gvv(ndm) ) * ( hks(nm,kd) + hks(nm,k+1) )
                            !
                            fac1 = 0.5 * hvisc * guu(nm ) * ( hkum(nm ,kd) + hkum(nm ,k+1) ) / gvu(nm )
                            !
                            fac2 = 0.5 * hvisc * guu(nmd) * ( hkum(nmd,kd) + hkum(nmd,k+1) ) / gvu(nmd)
                            !
                            rhsw(nm,k) = rhsw(nm,k) + ( fac1 * (w0(nmu,k) - w0(nm,k)) - fac2 * (w0(nm,k) - w0(nmd,k)) ) / fac
                            !
                         endif
                         !
                         ! viscosity term in y-direction
                         !
                         if ( wets(nm) * wets(ndm) * wets(num) == 1 ) then
                            !
                            fac = 0.25 * ( guu(nm) + guu(nmd) ) * ( hks(nm,kd) + hks(nm,k+1) )
                            !
                            fac1 = 0.5 * hvisc * gvv(nm ) * ( hkvm(nm ,kd) + hkvm(nm ,k+1) ) / guv(nm )
                            !
                            fac2 = 0.5 * hvisc * gvv(ndm) * ( hkvm(ndm,kd) + hkvm(ndm,k+1) ) / guv(ndm)
                            !
                            rhsw(nm,k) = rhsw(nm,k) + ( fac1 * (w0(num,k) - w0(nm,k)) - fac2 * (w0(nm,k) - w0(ndm,k)) ) / fac
                            !
                         endif
                         !
                      enddo
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
                   do n = nfu, nl
                      do m = mfu, ml
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
                         !
                         ! for permanently dry neighbours, corresponding values will be mirrored
                         !
                         if ( nmd  == 1 ) nmd  = nm
                         if ( nmu  == 1 ) nmu  = nm
                         if ( ndm  == 1 ) ndm  = nm
                         if ( num  == 1 ) num  = nm
                         if ( ndmd == 1 ) ndmd = ndm
                         !
                         ! viscosity term in x-direction
                         !
                         if ( wets(nm) * wets(nmd) * wets(nmu) == 1 ) then
                            !
                            fac = 0.25 * ( gvv(nm) + gvv(ndm) ) * ( hks(nm,kd) + hks(nm,k+1) )
                            !
                            fac1 = 0.25 * ( vnu2d(nm ) + vnu2d(ndm ) ) * guu(nm ) * ( hkum(nm ,kd) + hkum(nm ,k+1) ) / gvu(nm )
                            !
                            fac2 = 0.25 * ( vnu2d(nmd) + vnu2d(ndmd) ) * guu(nmd) * ( hkum(nmd,kd) + hkum(nmd,k+1) ) / gvu(nmd)
                            !
                            rhsw(nm,k) = rhsw(nm,k) + ( fac1 * (w0(nmu,k) - w0(nm,k)) - fac2 * (w0(nm,k) - w0(nmd,k)) ) / fac
                            !
                         endif
                         !
                         ! viscosity term in y-direction
                         !
                         if ( wets(nm) * wets(ndm) * wets(num) == 1 ) then
                            !
                            fac = 0.25 * ( guu(nm) + guu(nmd) ) * ( hks(nm,kd) + hks(nm,k+1) )
                            !
                            fac1 = 0.25 * ( vnu2d(nm ) + vnu2d(nmd ) ) * gvv(nm ) * ( hkvm(nm ,kd) + hkvm(nm ,k+1) ) / guv(nm )
                            !
                            fac2 = 0.25 * ( vnu2d(ndm) + vnu2d(ndmd) ) * gvv(ndm) * ( hkvm(ndm,kd) + hkvm(ndm,k+1) ) / guv(ndm)
                            !
                            rhsw(nm,k) = rhsw(nm,k) + ( fac1 * (w0(num,k) - w0(nm,k)) - fac2 * (w0(nm,k) - w0(ndm,k)) ) / fac
                            !
                         endif
                         !
                      enddo
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
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      nm = kgrpnt(m,n)
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
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
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
             enddo
             !
             ! to reduce the pressure Poisson equation set pressure of bottom face to that of top face for a number of layers
             !
             do l = 1, qlay
                !
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      nm = kgrpnt(m,n)
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
                enddo
                !
             enddo
             !
          else if ( ihydro == 2 ) then
             !
             ! central differences, so non-hydrostatic pressure is located at the cell centers
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
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
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
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
                         u = 0.25 * ( u0(nm,kd) + u0(nmd,kd) + u0(nm,k+1) + u0(nmd,k+1) )
                         v = 0.25 * ( v0(nm,kd) + v0(ndm,kd) + v0(nm,k+1) + v0(ndm,k+1) )
                         !
                         utot = sqrt( u*u + v*v + w0(nm,k)*w0(nm,k) )
                         !
                         amatw(nm,k,1) = amatw(nm,k,1) + apoks(nm,k) + bpoks(nm,k) * utot
                         !
                         amatw(nm,k,1) = amatw(nm,k,1) + cpoks(nm,k) / (dt*thetaw)
                         rhsw (nm,k  ) = rhsw (nm,k  ) + cpoks(nm,k) * w0(nm,k) / (dt*thetaw)
                         !
                      endif
                      !
                   enddo
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
                         amatw(nm,k,1) = amatw(nm,k,1) + 0.125 * ( cvegu(nm,kd,1) + cvegu(nmd,kd,1) + cvegu(nm,k+1,1) + cvegu(nmd,k+1,1) +  &
                                                                   cvegv(nm,kd,1) + cvegv(ndm,kd,1) + cvegv(nm,k+1,1) + cvegv(ndm,k+1,1) ) * abs(w0(nm,k))
                         !
                      endif
                      !
                   enddo
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
                      if ( wets(nm) == 1 ) then
                         !
                         rhsw(nm,k) = rhsw(nm,k) + 2. * ( rswu(nm,k) - rswu(nmd,k  ) ) / ( gvv(nm) + gvv(ndm) )    &
                                                 + 2. * ( rswv(nm,k) - rswv(ndm,k  ) ) / ( guu(nm) + guu(nmd) )    &
                                                 + 2. * ( rsww(nm,k) - rsww(nm ,k+1) ) / ( hks(nm,k) + hks(nm,k+1) )
                         !
                      endif
                      !
                   enddo
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
                if ( wets(nm) == 1 .and. presp(nm) == 1 ) then
                   !
                   dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
                   dyl = 0.5 * ( guu(nm) + guu(nmd) )
                   !
                   amatw(nm,0,1) = 1.
                   amatw(nm,0,2) = 0.
                   amatw(nm,0,3) = 0.
                   w1   (nm,0  ) = 0.25 * ( 3.*(u1(nm,1)+u1(nmd,1)) - u1(nm,2) - u1(nmd,2) ) * ( zku(nm,0) - zku(nmd,0) ) / dxl +  &
                                   0.25 * ( 3.*(v1(nm,1)+v1(ndm,1)) - v1(nm,2) - v1(ndm,2) ) * ( zkv(nm,0) - zkv(ndm,0) ) / dyl
                   !
                   ! add body motion
                   w1  (nm,0) = w1(nm,0) + skc(nm)
                   rhsw(nm,0) = w1(nm,0)
                   !
                endif
                !
             enddo
          enddo
          !
       endif
       !
       ! solve the w-momentum equation
       !
       if ( ihydro == 1 .or. ihydro == 2 ) then
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
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
          enddo
          !
          ! re-update the solution in case of thetaw <> 1
          !
          if ( thetaw /= 1. ) then
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
                      !
                      w1(nm,0) = ( w1(nm,0) - (1.-thetaw) * w0(nm,0) ) / thetaw
                      !
                   endif
                   !
                enddo
             enddo
             !
             do k = 1, kmax-1
                !
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      nm = kgrpnt(m,n)
                      !
                      if ( wets(nm) == 1 ) then
                         !
                         w1(nm,k) = ( w1(nm,k) - (1.-thetaw) * w0(nm,k) ) / thetaw
                         !
                      endif
                      !
                   enddo
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
                if ( wets(nm) == 1 ) then
                   !
                   dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
                   dyl = 0.5 * ( guu(nm) + guu(nmd) )
                   !
                   do k = 1, kmax
                      !
                      kd = max(k-1,1   )
                      ku = min(k+1,kmax)
                      !
                      fac1 = ( zku(nm,k-1) - zku(nmd,k-1) ) / ( 2.* dxl * ( hks(nm,k ) + hks(nm,kd) ) )
                      fac2 = ( zku(nm,k  ) - zku(nmd,k  ) ) / ( 2.* dxl * ( hks(nm,ku) + hks(nm,k ) ) )
                      !
                      if ( k == kmax ) fac2 = 0.
                      if ( k == 1 .and. presp(nm) == 1 ) fac1 = 0.
                      !
                      dmat(nm,k,1) =     -fac1 * hks(nm ,k)
                      dmat(nm,k,2) =   guu(nm) * hku(nm ,k) / gsqs(nm) - fac1 * hks(nm,kd) + fac2 * hks(nm,ku)
                      dmat(nm,k,3) =     -fac1 * hks(nm ,k)
                      dmat(nm,k,4) = -guu(nmd) * hku(nmd,k) / gsqs(nm) - fac1 * hks(nm,kd) + fac2 * hks(nm,ku)
                      dmat(nm,k,5) =      fac2 * hks(nm ,k)
                      dmat(nm,k,6) =      fac2 * hks(nm ,k)
                      !
                      fac1 = ( zkv(nm,k-1) - zkv(ndm,k-1) ) / ( 2.* dyl * ( hks(nm,k ) + hks(nm,kd) ) )
                      fac2 = ( zkv(nm,k  ) - zkv(ndm,k  ) ) / ( 2.* dyl * ( hks(nm,ku) + hks(nm,k ) ) )
                      !
                      if ( k == kmax ) fac2 = 0.
                      if ( k == 1 .and. presp(nm) == 1 ) fac1 = 0.
                      !
                      dmat(nm,k, 7) =     -fac1 * hks(nm ,k)
                      dmat(nm,k, 8) =   gvv(nm) * hkv(nm ,k) / gsqs(nm) - fac1 * hks(nm,kd) + fac2 * hks(nm,ku)
                      dmat(nm,k, 9) =     -fac1 * hks(nm ,k)
                      dmat(nm,k,10) = -gvv(ndm) * hkv(ndm,k) / gsqs(nm) - fac1 * hks(nm,kd) + fac2 * hks(nm,ku)
                      dmat(nm,k,11) =      fac2 * hks(nm ,k)
                      dmat(nm,k,12) =      fac2 * hks(nm ,k)
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
          enddo
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                dmat(nm,1, 2) = dmat(nm,1, 2) + 2.*dmat(nm,1,1)
                dmat(nm,1, 5) = dmat(nm,1, 5) -    dmat(nm,1,1)
                dmat(nm,1, 4) = dmat(nm,1, 4) + 2.*dmat(nm,1,3)
                dmat(nm,1, 6) = dmat(nm,1, 6) -    dmat(nm,1,3)
                dmat(nm,1, 1) = 0.
                dmat(nm,1, 3) = 0.
                !
                dmat(nm,1, 8) = dmat(nm,1, 8) + 2.*dmat(nm,1,7)
                dmat(nm,1,11) = dmat(nm,1,11) -    dmat(nm,1,7)
                dmat(nm,1,10) = dmat(nm,1,10) + 2.*dmat(nm,1,9)
                dmat(nm,1,12) = dmat(nm,1,12) -    dmat(nm,1,9)
                dmat(nm,1, 7) = 0.
                dmat(nm,1, 9) = 0.
                !
             enddo
          enddo
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
                   amatp(nm,k, 1) = dmat(nm,k, 2) * gmatu(nm ,k ,2) + dmat(nm,k, 4) * gmatu(nmd,k ,4) + dmat(nm,k, 1) * gmatu(nm ,kd,5) +  &
                                    dmat(nm,k, 3) * gmatu(nmd,kd,6) + dmat(nm,k, 5) * gmatu(nm ,ku,1) + dmat(nm,k, 6) * gmatu(nmd,ku,3) +  &
                                    dmat(nm,k, 8) * gmatv(nm ,k ,2) + dmat(nm,k,10) * gmatv(ndm,k ,4) + dmat(nm,k, 7) * gmatv(nm ,kd,5) +  &
                                    dmat(nm,k, 9) * gmatv(ndm,kd,6) + dmat(nm,k,11) * gmatv(nm ,ku,1) + dmat(nm,k,12) * gmatv(ndm,ku,3)
                   amatp(nm,k, 2) = dmat(nm,k, 4) * gmatu(nmd,k ,2) + dmat(nm,k, 3) * gmatu(nmd,kd,5) + dmat(nm,k, 6) * gmatu(nmd,ku,1)
                   amatp(nm,k, 3) = dmat(nm,k, 2) * gmatu(nm ,k ,4) + dmat(nm,k, 1) * gmatu(nm ,kd,6) + dmat(nm,k, 5) * gmatu(nm ,ku,3)
                   amatp(nm,k, 4) = dmat(nm,k, 2) * gmatu(nm ,k ,1) + dmat(nm,k, 4) * gmatu(nmd,k ,3) + dmat(nm,k, 1) * gmatu(nm ,kd,2) +  &
                                    dmat(nm,k, 3) * gmatu(nmd,kd,4) + dmat(nm,k, 8) * gmatv(nm ,k ,1) + dmat(nm,k,10) * gmatv(ndm,k ,3) +  &
                                    dmat(nm,k, 7) * gmatv(nm ,kd,2) + dmat(nm,k, 9) * gmatv(ndm,kd,4)
                   amatp(nm,k, 5) = dmat(nm,k, 2) * gmatu(nm ,k ,5) + dmat(nm,k, 4) * gmatu(nmd,k ,6) + dmat(nm,k, 5) * gmatu(nm ,ku,2) +  &
                                    dmat(nm,k, 6) * gmatu(nmd,ku,4) + dmat(nm,k, 8) * gmatv(nm ,k ,5) + dmat(nm,k,10) * gmatv(ndm,k ,6) +  &
                                    dmat(nm,k,11) * gmatv(nm ,ku,2) + dmat(nm,k,12) * gmatv(ndm,ku,4)
                   amatp(nm,k, 6) = dmat(nm,k,10) * gmatv(ndm,k ,2) + dmat(nm,k, 9) * gmatv(ndm,kd,5) + dmat(nm,k,12) * gmatv(ndm,ku,1)
                   amatp(nm,k, 7) = dmat(nm,k, 8) * gmatv(nm ,k ,4) + dmat(nm,k, 7) * gmatv(nm ,kd,6) + dmat(nm,k,11) * gmatv(nm ,ku,3)
                   amatp(nm,k, 8) = dmat(nm,k, 2) * gmatu(nm ,k ,3) + dmat(nm,k, 1) * gmatu(nm ,kd,4)
                   amatp(nm,k, 9) = dmat(nm,k, 8) * gmatv(nm ,k ,3) + dmat(nm,k, 7) * gmatv(nm ,kd,4)
                   amatp(nm,k,10) = dmat(nm,k, 4) * gmatu(nmd,k ,1) + dmat(nm,k, 3) * gmatu(nmd,kd,2)
                   amatp(nm,k,11) = dmat(nm,k,10) * gmatv(ndm,k ,1) + dmat(nm,k, 9) * gmatv(ndm,kd,2)
                   amatp(nm,k,12) = dmat(nm,k, 2) * gmatu(nm ,k ,6) + dmat(nm,k, 5) * gmatu(nm ,ku,4)
                   amatp(nm,k,13) = dmat(nm,k, 8) * gmatv(nm ,k ,6) + dmat(nm,k,11) * gmatv(nm ,ku,4)
                   amatp(nm,k,14) = dmat(nm,k, 4) * gmatu(nmd,k ,5) + dmat(nm,k, 6) * gmatu(nmd,ku,2)
                   amatp(nm,k,15) = dmat(nm,k,10) * gmatv(ndm,k ,5) + dmat(nm,k,12) * gmatv(ndm,ku,2)
                   amatp(nm,k,16) = dmat(nm,k, 1) * gmatu(nm ,kd,1) + dmat(nm,k, 3) * gmatu(nmd,kd,3) + dmat(nm,k, 7) * gmatv(nm,kd,1) + dmat(nm,k, 9) * gmatv(ndm,kd,3)
                   amatp(nm,k,17) = dmat(nm,k, 5) * gmatu(nm ,ku,5) + dmat(nm,k, 6) * gmatu(nmd,ku,6) + dmat(nm,k,11) * gmatv(nm,ku,5) + dmat(nm,k,12) * gmatv(ndm,ku,6)
                   amatp(nm,k,18) = dmat(nm,k, 1) * gmatu(nm ,kd,3)
                   amatp(nm,k,19) = dmat(nm,k, 7) * gmatv(nm ,kd,3)
                   amatp(nm,k,20) = dmat(nm,k, 3) * gmatu(nmd,kd,1)
                   amatp(nm,k,21) = dmat(nm,k, 9) * gmatv(ndm,kd,1)
                   amatp(nm,k,22) = dmat(nm,k, 5) * gmatu(nm ,ku,6)
                   amatp(nm,k,23) = dmat(nm,k,11) * gmatv(nm ,ku,6)
                   amatp(nm,k,24) = dmat(nm,k, 6) * gmatu(nmd,ku,5)
                   amatp(nm,k,25) = dmat(nm,k,12) * gmatv(ndm,ku,5)
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
                   rhsp(nm,k) = ( dmat(nm,k,2)*u1(nm ,k ) + dmat(nm,k, 4)*u1(nmd,k ) + dmat(nm,k, 1)*u1(nm ,kd) +  &
                                  dmat(nm,k,3)*u1(nmd,kd) + dmat(nm,k, 5)*u1(nm ,ku) + dmat(nm,k, 6)*u1(nmd,ku) +  &
                                  dmat(nm,k,8)*v1(nm ,k ) + dmat(nm,k,10)*v1(ndm,k ) + dmat(nm,k, 7)*v1(nm ,kd) +  &
                                  dmat(nm,k,9)*v1(ndm,kd) + dmat(nm,k,11)*v1(nm ,ku) + dmat(nm,k,12)*v1(ndm,ku) +  &
                                  kwd*w1(nm,k-1) - kwu*w1(nm,k) ) / (dt*theta3)
                   !
                   ! add contribution from rigid body motions (top layer only)
                   !
                   if ( k == 1 ) rhsp(nm,k) = rhsp(nm,k) + skc(nm) / (dt*theta3)
                   !
                enddo
                !
             enddo
          enddo
          !
          ! reduce the pressure Poisson equation
          !
          do l = 1, qlay
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   amatp(nm,qmax, 1) = amatp(nm,qmax, 1) + real(qlay+2-l)*amatp(nm,kmax-l+1, 4)
                   amatp(nm,qmax, 2) = amatp(nm,qmax, 2) + real(qlay+2-l)*amatp(nm,kmax-l+1,10)
                   amatp(nm,qmax, 3) = amatp(nm,qmax, 3) + real(qlay+2-l)*amatp(nm,kmax-l+1, 8)
                   amatp(nm,qmax, 6) = amatp(nm,qmax, 6) + real(qlay+2-l)*amatp(nm,kmax-l+1,11)
                   amatp(nm,qmax, 7) = amatp(nm,qmax, 7) + real(qlay+2-l)*amatp(nm,kmax-l+1, 9)
                   amatp(nm,qmax, 4) = amatp(nm,qmax, 4) + real(qlay+2-l)*amatp(nm,kmax-l+1,16)
                   amatp(nm,qmax, 8) = amatp(nm,qmax, 8) + real(qlay+2-l)*amatp(nm,kmax-l+1,18)
                   amatp(nm,qmax,10) = amatp(nm,qmax,10) + real(qlay+2-l)*amatp(nm,kmax-l+1,20)
                   amatp(nm,qmax, 9) = amatp(nm,qmax, 9) + real(qlay+2-l)*amatp(nm,kmax-l+1,19)
                   amatp(nm,qmax,11) = amatp(nm,qmax,11) + real(qlay+2-l)*amatp(nm,kmax-l+1,21)
                   rhsp (nm,qmax   ) = rhsp (nm,qmax   ) + real(qlay+2-l)*rhsp (nm,kmax-l+1   )
                   !
                enddo
             enddo
             !
          enddo
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
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
          enddo
          !
          ! solve the Poisson equation
          !
          if ( qmax == 1 ) then
             !
             amat(:,1) = amatp(:,1,1)
             amat(:,2) = amatp(:,1,2)
             amat(:,3) = amatp(:,1,3)
             amat(:,4) = amatp(:,1,6)
             amat(:,5) = amatp(:,1,7)
             call sip( amat, rhsp(1,1), dq(1,1) )
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
             do n = nfu, nl
                do m = mf, ml
                   !
                   mu = m + 1
                   !
                   nm  = kgrpnt(m ,n)
                   nmu = kgrpnt(mu,n)
                   !
                   if ( uwetp(nm) == 1 .and. presu(nm) == 0 ) then
                      !
                      dqgrdu(nm,k) = gmatu(nm,k,1)*dq(nm ,kd) + gmatu(nm,k,2)*dq(nm,k ) + gmatu(nm,k,3)*dq(nmu,kd)  &
                                   + gmatu(nm,k,4)*dq(nmu,k ) + gmatu(nm,k,5)*dq(nm,ku) + gmatu(nm,k,6)*dq(nmu,ku)
                      !
                   else
                      ! no correction applied on piezometric head
                      !
                      dqgrdu(nm,k) = 0.
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
                   if ( vwetp(nm) == 1 .and. presv(nm) == 0 ) then
                      !
                      dqgrdv(nm,k) = gmatv(nm,k,1)*dq(nm ,kd) + gmatv(nm,k,2)*dq(nm,k ) + gmatv(nm,k,3)*dq(num,kd)  &
                                   + gmatv(nm,k,4)*dq(num,k ) + gmatv(nm,k,5)*dq(nm,ku) + gmatv(nm,k,6)*dq(num,ku)
                      !
                   else
                      ! no correction applied on piezometric head
                      !
                      dqgrdv(nm,k) = 0.
                      !
                   endif
                   !
                enddo
             enddo
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
          do n = nfu, nl
             do m = mf, ml
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                if ( uwetp(nm) == 1 ) then
                   !
                   u1(nm,k) = u1(nm,k) - dt*theta3*( gmatu(nm,k,1)*dq(nm ,kd) + gmatu(nm,k,2)*dq(nm,k ) + gmatu(nm,k,3)*dq(nmu,kd) + &
                                                     gmatu(nm,k,4)*dq(nmu,k ) + gmatu(nm,k,5)*dq(nm,ku) + gmatu(nm,k,6)*dq(nmu,ku) )
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! v-velocity
       !
       do k = 1, kmax
          !
          kd = max(k-1,1   )
          ku = min(k+1,kmax)
          !
          do m = mfu, ml
             do n = nf, nl
                !
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                if ( vwetp(nm) == 1 ) then
                   !
                   v1(nm,k) = v1(nm,k) - dt*theta3*( gmatv(nm,k,1)*dq(nm ,kd) + gmatv(nm,k,2)*dq(nm,k ) + gmatv(nm,k,3)*dq(num,kd) + &
                                                     gmatv(nm,k,4)*dq(num,k ) + gmatv(nm,k,5)*dq(nm,ku) + gmatv(nm,k,6)*dq(num,ku) )
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! w-velocity
       !
       do k = 0, kmax-1
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
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
          enddo
          !
       enddo
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
                w1(nm,kmax) = 0.25 * ( 3.*(u1(nm,kmax)+u1(nmd,kmax)) - u1(nm,kmax-1) - u1(nmd,kmax-1) ) * ( zku(nm,kmax) - zku(nmd,kmax) ) / dxl +  &
                              0.25 * ( 3.*(v1(nm,kmax)+v1(ndm,kmax)) - v1(nm,kmax-1) - v1(ndm,kmax-1) ) * ( zkv(nm,kmax) - zkv(ndm,kmax) ) / dyl
                !
             endif
             !
             if ( wets(nm) == 1 .and. presp(nm) == 1 ) then
                !
                dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
                dyl = 0.5 * ( guu(nm) + guu(nmd) )
                !
                w1(nm,0) = 0.25 * ( 3.*(u1(nm,1)+u1(nmd,1)) - u1(nm,2) - u1(nmd,2) ) * ( zku(nm,0) - zku(nmd,0) ) / dxl +  &
                           0.25 * ( 3.*(v1(nm,1)+v1(ndm,1)) - v1(nm,2) - v1(ndm,2) ) * ( zkv(nm,0) - zkv(ndm,0) ) / dyl
                !
                ! add body motion
                w1(nm,0) = w1(nm,0) + skc(nm)
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
          if ( ibl(n) == 1 .and. LMXF .and. iturb < 2 ) then
             !
             ! set to zero at closed boundary
             !
             v1(nmf,:) = 0.
             !
          else if ( LMXF ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             v1(nmf,:) = v1(nmfu,:)
             !
          endif
          !
          if ( ibr(n) == 1 .and. LMXL .and. iturb < 2 ) then
             !
             ! set to zero at closed boundary
             !
             v1(nmlu,:) = 0.
             !
          else if ( LMXL ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             v1(nmlu,:) = v1(nml,:)
             !
          endif
          !
       enddo
       !
    else
       !
       call periodic ( v1, kgrpnt, 1, kmax )
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
          if ( ibb(m) == 1 .and. LMYF .and. iturb < 2 ) then
             !
             ! set to zero at closed boundary
             !
             u1(nfm,:) = 0.
             !
          else if ( LMYF ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             u1(nfm,:) = u1(nfum,:)
             !
          endif
          !
          if ( ibt(m) == 1 .and. LMYL .and. iturb < 2 ) then
             !
             ! set to zero at closed boundary
             !
             u1(nlum,:) = 0.
             !
          else if ( LMYL ) then
             !
             ! zero normal gradient condition at open boundaries
             !
             u1(nlum,:) = u1(nlm,:)
             !
          endif
          !
       enddo
       !
    else
       !
       call periodic ( u1, kgrpnt, 1, kmax )
       !
    endif
    !
    ! compute the mass flux in x-direction
    !
    do k = 1, kmax
       !
       do n = nfu, nl
          do m = mf, ml
             !
             nm = kgrpnt(m,n)
             !
             qx(nm,k) = hku(nm,k) * guu(nm) * ( tetau(nm)*u1(nm,k) + (1.-tetau(nm))*u0(nm,k) )
             !
          enddo
       enddo
       !
    enddo
    !
    ! compute the mass flux in y-direction
    !
    do k = 1, kmax
       !
       do n = nf, nl
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             qy(nm,k) = hkv(nm,k) * gvv(nm) * ( tetav(nm)*v1(nm,k) + (1.-tetav(nm))*v0(nm,k) )
             !
          enddo
       enddo
       !
    enddo
    !
    ! compute the depth-averaged velocities
    !
    udep = 0.
    !
    do n = nfu, nl
       do m = mf, ml
          !
          nm = kgrpnt(m,n)
          !
          if ( uwetp(nm) == 1 ) then
             !
             do k = 1, kmax
                !
                udep(nm) = udep(nm) + hku(nm,k)*guu(nm)*u1(nm,k)
                !
             enddo
             !
             udep(nm) = udep(nm) / ( hu(nm)*guu(nm) )
             !
          endif
          !
       enddo
    enddo
    !
    vdep = 0.
    !
    do m = mfu, ml
       do n = nf, nl
          !
          nm = kgrpnt(m,n)
          !
          if ( vwetp(nm) == 1 ) then
             !
             do k = 1, kmax
                !
                vdep(nm) = vdep(nm) + hkv(nm,k)*gvv(nm)*v1(nm,k)
                !
             enddo
             !
             vdep(nm) = vdep(nm) / ( hv(nm)*gvv(nm) )
             !
          endif
          !
       enddo
    enddo
    !
    ! impose Neumann condition for w-velocity at boundaries, if appropriate
    !
    if ( horwinc ) then
       !
       if ( ihydro == 1 .or. ihydro == 2 ) then
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
                if ( LMXF ) w1(nmf ,:) = w1(nmfu,:)
                if ( LMXL ) w1(nmlu,:) = w1(nml ,:)
                !
             enddo
             !
          else
             !
             call periodic ( w1, kgrpnt, 0, kmax )
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
                if ( LMYF ) w1(nfm ,:) = w1(nfum,:)
                if ( LMYL ) w1(nlum,:) = w1(nlm ,:)
                !
             enddo
             !
          else
             !
             call periodic ( w1, kgrpnt, 0, kmax )
             !
          endif
          !
       endif
       !
    endif
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend(1) = ml - 1
    if ( LMXL ) mend(1) = ml
    !
    nend(1) = nl - 1
    if ( LMYL ) nend(1) = nl
    !
    ! calculate net mass outflow based on local continuity equation
    !
    if ( ITEST >= 30 ) then
       !
       moutf = 0.
       !
       if ( ihydro == 1 .or. ihydro == 2 ) then
          !
          do n = nfu, nend(1)
             do m = mfu, mend(1)
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
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
                   if ( nm > 1 ) moutf = moutf + dmat(nm,k,2)*u1(nm ,k ) + dmat(nm,k, 4)*u1(nmd,k ) + dmat(nm,k, 1)*u1(nm ,kd) +  &
                                                 dmat(nm,k,3)*u1(nmd,kd) + dmat(nm,k, 5)*u1(nm ,ku) + dmat(nm,k, 6)*u1(nmd,ku) +  &
                                                 dmat(nm,k,8)*v1(nm ,k ) + dmat(nm,k,10)*v1(ndm,k ) + dmat(nm,k, 7)*v1(nm ,kd) +  &
                                                 dmat(nm,k,9)*v1(ndm,kd) + dmat(nm,k,11)*v1(nm ,ku) + dmat(nm,k,12)*v1(ndm,ku) +  &
                                                 kwd*w1(nm,k-1) - kwu*w1(nm,k)
                   !
                enddo
                !
             enddo
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
       do n = nfu, nend(1)
          do m = mfu, mend(1)
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
                ener = ener + 0.125 * ( fac1 + fac2 ) * udep(nmd) * udep(nmd)
                fac2 = ( min( -flos(nmu), s1(nmu) ) + dps(nmu) ) * gsqs(nmu)
                ener = ener + 0.125 * ( fac1 + fac2 ) * udep(nm ) * udep(nm )
                !
                fac2 = ( min( -flos(ndm), s1(ndm) ) + dps(ndm) ) * gsqs(ndm)
                ener = ener + 0.125 * ( fac1 + fac2 ) * vdep(ndm) * vdep(ndm)
                fac2 = ( min( -flos(num), s1(num) ) + dps(num) ) * gsqs(num)
                ener = ener + 0.125 * ( fac1 + fac2 ) * vdep(nm ) * vdep(nm )
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
          do n = nfu, nend(1)
             !
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             nml  = kgrpnt(ml ,n)
             nmlu = kgrpnt(mlu,n)
             !
             do m = mfu, mend(1)
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
                !
                if ( wets(nm) == 1 ) then
                   !
                   wom(nm,k) = wom(nm,k+1) - (      tetau(nm)*hku(nm,k+1)*guu(nm)*u1(nm,k+1) -      tetau(nmd)*hku(nmd,k+1)*guu(nmd)*u1(nmd,k+1) +  &
                                                    tetav(nm)*hkv(nm,k+1)*gvv(nm)*v1(nm,k+1) -      tetav(ndm)*hkv(ndm,k+1)*gvv(ndm)*v1(ndm,k+1) ) / gsqs(nm)  &
                                           - ( (1.-tetau(nm))*hku(nm,k+1)*guu(nm)*u0(nm,k+1) - (1.-tetau(nmd))*hku(nmd,k+1)*guu(nmd)*u0(nmd,k+1) +  &
                                               (1.-tetav(nm))*hkv(nm,k+1)*gvv(nm)*v0(nm,k+1) - (1.-tetav(ndm))*hkv(ndm,k+1)*gvv(ndm)*v0(ndm,k+1) ) / gsqs(nm)  &
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
             if ( .not.lreptx ) then
                !
                if ( LMXF ) wom(nmf ,k) = wom(nmfu,k)
                if ( LMXL ) wom(nmlu,k) = wom(nml ,k)
                !
             endif
             !
          enddo
          !
          if ( .not.lrepty ) then
             !
             do m = mf, mlu
                !
                nfm  = kgrpnt(m,nf )
                nfum = kgrpnt(m,nfu)
                nlm  = kgrpnt(m,nl )
                nlum = kgrpnt(m,nlu)
                !
                if ( LMYF ) wom(nfm ,k) = wom(nfum,k)
                if ( LMYL ) wom(nlum,k) = wom(nlm ,k)
                !
             enddo
             !
          endif
          !
       enddo
       !
       ! check if relative vertical velocity at surface is zero
       !
       do n = nfu, nend(1)
          do m = mfu, mend(1)
             !
             nm = kgrpnt(m,n)
             !
             if ( abs(wom(nm,0)) > epswom ) then
                !
                write (msgstr,'(a,i4,a,i4,a,e9.3,a)') 'nonzero relative vertical velocity at surface in (m,n)=(',m+MXF-2,',',n+MYF-2,'); omega = ',wom(nm,0),' m/s'
                call msgerr (2, trim(msgstr) )
                !
                wom(nm,0) = 0.
                !
             endif
             !
          enddo
       enddo
       !
    else
       !
       do n = nfu, nend(1)
          !
          nmf  = kgrpnt(mf ,n)
          nmfu = kgrpnt(mfu,n)
          nml  = kgrpnt(ml ,n)
          nmlu = kgrpnt(mlu,n)
          !
          do m = mfu, mend(1)
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
                do k = 1, kmax-1
                   !
                   fac1 = 0.5 * ( u1(nmd,k  ) + u1(nm,k  ) )
                   fac2 = 0.5 * ( u1(nmd,k+1) + u1(nm,k+1) )
                   !
                   u    = ( fac2 * hks(nm,k) + fac1 * hks(nm,k+1) ) / ( hks(nm,k) + hks(nm,k+1) )
                   !
                   fac1 = 0.5 * ( v1(ndm,k  ) + v1(nm,k  ) )
                   fac2 = 0.5 * ( v1(ndm,k+1) + v1(nm,k+1) )
                   !
                   v    = ( fac2 * hks(nm,k) + fac1 * hks(nm,k+1) ) / ( hks(nm,k) + hks(nm,k+1) )
                   !
                   wom(nm,k) = w1(nm,k) - ( zksnew(nm,k) - zks(nm,k) ) / dt - u * ( zku(nm,k) - zku(nmd,k) ) / dxl - v * ( zkv(nm,k) - zkv(ndm,k) ) / dyl
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
                   wom(nm,k) = wom(nm,k+1) - (      tetau(nm)*hku(nm,k+1)*guu(nm)*u1(nm,k+1) -      tetau(nmd)*hku(nmd,k+1)*guu(nmd)*u1(nmd,k+1) +  &
                                                    tetav(nm)*hkv(nm,k+1)*gvv(nm)*v1(nm,k+1) -      tetav(ndm)*hkv(ndm,k+1)*gvv(ndm)*v1(ndm,k+1) ) / gsqs(nm)  &
                                           - ( (1.-tetau(nm))*hku(nm,k+1)*guu(nm)*u0(nm,k+1) - (1.-tetau(nmd))*hku(nmd,k+1)*guu(nmd)*u0(nmd,k+1) +  &
                                               (1.-tetav(nm))*hkv(nm,k+1)*gvv(nm)*v0(nm,k+1) - (1.-tetav(ndm))*hkv(ndm,k+1)*gvv(ndm)*v0(ndm,k+1) ) / gsqs(nm)  &
                                           - ( hksnew(nm,k+1) - hks(nm,k+1) ) / dt
                   !
                enddo
                !
                if ( abs(wom(nm,0)) > epswom ) then
                   !
                   write (msgstr,'(a,i4,a,i4,a,e9.3,a)') 'nonzero relative vertical velocity at surface in (m,n)=(',m+MXF-2,',',n+MYF-2,'); omega = ',wom(nm,0),' m/s'
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
          if ( .not.lreptx ) then
             !
             if ( LMXF ) wom(nmf ,:) = wom(nmfu,:)
             if ( LMXL ) wom(nmlu,:) = wom(nml ,:)
             !
          endif
          !
       enddo
       !
       if ( .not.lrepty ) then
          !
          do m = mf, mlu
             !
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             nlm  = kgrpnt(m,nl )
             nlum = kgrpnt(m,nlu)
             !
             if ( LMYF ) wom(nfm ,:) = wom(nfum,:)
             if ( LMYL ) wom(nlum,:) = wom(nlm ,:)
             !
          enddo
          !
       endif
       !
    endif
    !
    ! exchange velocities with neighbouring subdomains (if appropriate)
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       call SWEXCHG ( u1, kgrpnt, 1, kmax )
       call SWEXCHG ( v1, kgrpnt, 1, kmax )
       if ( horwinc ) call SWEXCHG ( w1, kgrpnt, 0, kmax )
       !
    endif
    call SWEXCHG ( wom, kgrpnt, 0, kmax )
    if (STPNOW()) return
    !
    ! synchronize velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       call periodic ( u1, kgrpnt, 1, kmax )
       call periodic ( v1, kgrpnt, 1, kmax )
       if ( horwinc ) call periodic ( w1, kgrpnt, 0, kmax )
       !
    endif
    call periodic ( wom, kgrpnt, 0, kmax )
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
             tbndx(n,       1:  kmax) = u1(nmf ,1:kmax)
             tbndx(n,  kmax+1:2*kmax) = v1(nmfu,1:kmax)
             tbndx(n,2*kmax+1       ) = s1(nmfu)
          enddo
       endif
       !
       if ( kpart == 4 ) then
          call SWREDUCE ( tbndx, (2*kmax+1)*myc, SWREAL, SWSUM )
          if (STPNOW()) return
       endif
       !
       do n = nfu, nl
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,n)
             !
             if ( uwetp(nm) == 1 ) then
                !
                u1(nm,1:kmax) = (1.-sponxl(nm))*u1(nm,1:kmax) + sponxl(nm)*tbndx(n,1:kmax)
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
                v1(nm,1:kmax) = (1.-sponxl(nmd))*v1(nm,1:kmax) + sponxl(nmd)*tbndx(n,  kmax+1:2*kmax)
                s1(nm       ) = (1.-sponxl(nmd))*s1(nm       ) + sponxl(nmd)*tbndx(n,2*kmax+1       )
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
             tbndx(n,       1:  kmax) = u1(nml,1:kmax)
             tbndx(n,  kmax+1:2*kmax) = v1(nml,1:kmax)
             tbndx(n,2*kmax+1       ) = s1(nml)
          enddo
       endif
       !
       if ( kpart == 4 ) then
          call SWREDUCE ( tbndx, (2*kmax+1)*myc, SWREAL, SWSUM )
          if (STPNOW()) return
       endif
       !
       do n = nfu, nl
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,n)
             !
             if ( uwetp(nm) == 1 ) then
                !
                u1(nm,1:kmax) = (1.-sponxr(nm))*u1(nm,1:kmax) + sponxr(nm)*tbndx(n,1:kmax)
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
                v1(nm,1:kmax) = (1.-sponxr(nm))*v1(nm,1:kmax) + sponxr(nm)*tbndx(n,  kmax+1:2*kmax)
                s1(nm       ) = (1.-sponxr(nm))*s1(nm       ) + sponxr(nm)*tbndx(n,2*kmax+1       )
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
             tbndy(m,       1:  kmax) = v1(nfm ,1:kmax)
             tbndy(m,  kmax+1:2*kmax) = u1(nfum,1:kmax)
             tbndy(m,2*kmax+1       ) = s1(nfum)
          enddo
       endif
       !
       if ( kpart == 3 ) then
          call SWREDUCE ( tbndy, (2*kmax+1)*mxc, SWREAL, SWSUM )
          if (STPNOW()) return
       endif
       !
       do m = mfu, ml
          !
          do n = nf, nl
             !
             nm = kgrpnt(m,n)
             !
             if ( vwetp(nm) == 1 ) then
                !
                v1(nm,1:kmax) = (1.-sponyb(nm))*v1(nm,1:kmax) + sponyb(nm)*tbndy(m,1:kmax)
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
                u1(nm,1:kmax) = (1.-sponyb(ndm))*u1(nm,1:kmax) + sponyb(ndm)*tbndy(m,  kmax+1:2*kmax)
                s1(nm       ) = (1.-sponyb(ndm))*s1(nm       ) + sponyb(ndm)*tbndy(m,2*kmax+1       )
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
             tbndy(m,       1:  kmax) = v1(nlm,1:kmax)
             tbndy(m,  kmax+1:2*kmax) = u1(nlm,1:kmax)
             tbndy(m,2*kmax+1       ) = s1(nlm)
          enddo
       endif
       !
       if ( kpart == 3 ) then
          call SWREDUCE ( tbndy, (2*kmax+1)*mxc, SWREAL, SWSUM )
          if (STPNOW()) return
       endif
       !
       do m = mfu, ml
          !
          do n = nf, nl
             !
             nm = kgrpnt(m,n)
             !
             if ( vwetp(nm) == 1 ) then
                !
                v1(nm,1:kmax) = (1.-sponyt(nm))*v1(nm,1:kmax) + sponyt(nm)*tbndy(m,1:kmax)
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
                u1(nm,1:kmax) = (1.-sponyt(nm))*u1(nm,1:kmax) + sponyt(nm)*tbndy(m,  kmax+1:2*kmax)
                s1(nm       ) = (1.-sponyt(nm))*s1(nm       ) + sponyt(nm)*tbndy(m,2*kmax+1       )
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
    do k = 1, kmax
       !
       do n = nfu, nend(1)
          do m = mfu, mend(1)
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
                u    = 0.5 * ( u1(nm,k) + u1(nmd,k) )
                v    = 0.5 * ( v1(nm,k) + v1(ndm,k) )
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
end subroutine SwashImpLay2DHflow
