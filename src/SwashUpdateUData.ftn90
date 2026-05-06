subroutine SwashUpdateUData ( it )
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
!   11.01: Panagiotis Vasarmidis
!   11.02: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, February 2020: New subroutine
!   11.01,  October 2023: extension first and second order Stokes wave at boundaries
!   11.02, November 2023: extension second order transfer functions at boundaries
!
!   Purpose
!
!   Updates flow data, boundary conditions and input fields in case of flexible mesh
!
!   Method
!
!   The following boundary conditions for flow can be imposed:
!
!   1) water level
!   2) velocity or discharge
!   3) Riemann invariant or Sommerfeld radiation or weakly reflective
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use SwashTimecomm
    use m_bndspec
    use m_genarr
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)   :: it       ! integration step counter
!
!   Local variables
!
    integer               :: btype    ! boundary type
                                      ! = 1; closed
                                      ! = 2; water level opening
                                      ! = 3; velocity opening
                                      ! = 5; discharge opening
                                      ! = 6; Riemann invariant opening
                                      ! = 7; weakly reflective opening
                                      ! = 8; Sommerfeld radiation condition
                                      ! =10; outflow condition
                                      ! < 0; =-btype, where both normal and tangential components are described
    integer               :: i        ! loop counter
    integer               :: ibloc    ! counter for boundary point
    integer               :: ibnd     ! bound long wave is added (=1) or not added (=0)
    integer               :: icell    ! loop counter over cells
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer               :: iface    ! face index
    integer               :: ifpco    ! frequency with which preconditioner is applied
    integer               :: indx     ! boundary point index
    integer               :: indxf    ! index of boundary face
    integer               :: indxn    ! index of next boundary point
    integer               :: indxs    ! index of boundary cell
    integer               :: istok    ! indicates the order of Stokes wave theory
                                      ! = 0; use hyperbolic cosine distribution for velocity (Airy wave theory)
                                      ! = 1; first order Stokes wave (Airy wave theory)
                                      ! = 2; second order Stokes wave
                                      ! = 3; second order sub- and super-harmonic transfer functions
    integer               :: itmp     ! temporary stored integer
    integer               :: ix       ! index of point
    integer               :: j        ! loop counter
    integer               :: jb       ! loop counter over boundary faces/cells
    integer               :: jj       ! loop counter
    integer               :: k        ! loop counter over vertical layers
    integer               :: k1       ! user-defined location of boundary point
    integer               :: k2       ! user-defined location of subsequent boundary point
    integer               :: klay     ! layer number
    integer, dimension(1) :: kx       ! location of minimum value in array of x-coordinates of boundary vertices
    integer, dimension(1) :: ky       ! location of minimum value in array of y-coordinates of boundary vertices
    integer               :: nfreq    ! number of components in Fourier series
    integer               :: shape    ! spectral shape
                                      ! = 1; Pierson Moskowitz
                                      ! = 2; Jonswap
                                      ! = 3; TMA
    integer, dimension(2) :: v        ! vertices of boundary face
    integer               :: vm       ! boundary marker
    !
    real                  :: ampl     ! amplitude of a Fourier component
    real                  :: azero    ! amplitude at zero frequency in Fourier series
    real                  :: beta     ! shape factor beta of source area
    real                  :: bval     ! actual boundary value
    real                  :: cgsrc    ! center of gravity of source area
    real                  :: d1       ! distance of a point to origin
    real                  :: d2       ! distance of another point to origin
    real                  :: dcor     ! correction in depth integration to get zero mass outflow
    real                  :: dep      ! local water depth
    real                  :: dist     ! distance to boundary
    real, save            :: epsab2   ! parameter to modify the weights of Adams-Bashforth scheme
    real                  :: fac      ! a factor
    real                  :: fac1     ! another factor
    real                  :: fac2     ! some other factor
    real                  :: fsmo     ! ramp function for smoothing incident waves
    real                  :: fval     ! value from a Fourier series
    real                  :: fvalu    ! value from a Fourier series for velocity component
    real                  :: fvalb    ! value from bound long wave
    real                  :: fvalbu   ! value from bound long wave for velocity component
    real                  :: hwidt    ! half width of source area
    real                  :: kwav     ! wave number
    real                  :: nx       ! x-component of normal to face
    real                  :: ny       ! y-component of normal to face
    real                  :: omega    ! angular frequency of a Fourier component
    real                  :: omegsb   ! angular frequency of bound sub-harmonic component
    real                  :: omegsp   ! angular frequency of bound super-harmonic component
    real                  :: phase    ! phase of a Fourier component
    real                  :: rdx      ! reciprocal of distance between cell circumcenter and boundary face
    real                  :: rsgn     ! sign for indicating in- and outflowing depending on boundary
                                      ! =+1; refers to inflowing at left and lower boundaries
                                      ! =-1; refers to outflowing at right and upper boundaries
    real                  :: sfval    ! value from an internal-generated series of wave components
    real                  :: swd      ! still water depth
    real                  :: theta    ! randomly chosen directions assigned to each frequency
    real                  :: tsmo     ! period for smoothing boundary values during cold start
    real                  :: wdir     ! incident or peak wave direction with respect to problem coordinates
    real                  :: wf       ! weighting factor
    real                  :: wf1      ! weighting factor for given boundary value to be interpolated in space
    real                  :: wf2      ! weighting factor for subsequent boundary value to be interpolated in space
    real                  :: x0       ! x-coordinate of reference point
    real                  :: xb       ! x-coordinate of boundary cell
    real, save            :: xp       ! x-coordinate of grid point
    real                  :: y0       ! y-coordinate of reference point
    real                  :: yb       ! y-coordinate of boundary cell
    real, save            :: yp       ! y-coordinate of grid point
    real                  :: z0       ! roughness height for logarithmic velocity profile
    real                  :: z1       ! interface at bottom of layer
    real                  :: z2       ! interface at top of layer
    !
    logical               :: lpb      ! indicates whether test point is a boundary point
    logical               :: lriem    ! indicates whether linearized Riemann invariant is imposed or not
    logical               :: STPNOW   ! indicates that program must stop
    logical               :: vdir     ! indicates direction of in- or outcoming velocity on boundary
                                      ! =.true.; eastward
                                      ! =.false.; northward
    !
    character(80)         :: msgstr   ! string to pass message
    !
    type(bfldat), pointer :: curbfl   ! current item in list of boundary condition file
    type(bfsdat), pointer :: curbfs   ! current item in list of Fourier series parameters
    !
    type(celltype), dimension(:), pointer :: cell ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdateUData')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       ! indicate whether preconditioner should be applied or not
       !
       ifpco = nint(pnums(28))
       !
       if ( ifpco < 0 ) then
          lprecon = .false.
       else if ( mod(it-1,ifpco) == 0 ) then
          lprecon = .true.
       else
          lprecon = .false.
       endif
       !
    endif
    !
    ! store flow variables to previous time level
    !
    if ( nstatc == 1 ) then
       !
       if ( relwav ) so = s0
       s0 = s1
       u0 = u1
       !
       if ( ihydro /= 0 ) then
          !
          if ( kmax == 1 .or. ihydro == 3 ) then
             !
             w0bot = w1bot
             w0top = w1top
             !
          else
             !
             w0 = w1
             !
          endif
          !
       endif
       !
       hso   = hs
       bcso  = bcs
       if ( kmax   > 1 ) hkso = hks
       if ( ltrans > 0 ) bcrpo = bcrp
       !
       if ( it == 1 ) then
          epsab2 = pcor(2)
          ! at first time step use explicit Euler scheme to Coriolis forcing term
          pcor(2) = -0.5
       else if ( it > 1 ) then
          ! from second time step on use modified AB2 scheme to Coriolis forcing term
          pcor(2) = epsab2
       endif
       !
    endif
    !
    ! update boundary conditions
    !
    ! compute boundary values based on Fourier series
    !
    curbfs => fbfs
    do
       ibloc = curbfs%nbfs
       if ( ibloc == -999 ) exit       ! no Fourier series given
       !
       nfreq = curbfs%nfreq
       azero = curbfs%azero
       if ( .not. azero < -0.9e10 ) then
          !
          fval = azero
          do i = 1, nfreq
             ampl  = curbfs%ampl (i)
             omega = curbfs%omega(i)
             phase = curbfs%phase(i)
             fval = fval + ampl * cos( omega*timco )*cos( phase ) + ampl * sin( omega*timco )*sin( phase )
          enddo
          bndval(ibloc,1) = fval
          !
       endif
       !
       if ( .not.associated(curbfs%nextbfs) ) exit
       curbfs => curbfs%nextbfs
       !
    enddo
    !
    ! compute boundary values based on time series
    !
    if ( ITEST >= 80 ) write (PRTEST,*) ' number of boundary files with time series ', nbfils
    !
    curbfl => fbndfil
    !
    do i = 1, nbfils
       !
       ! read boundary values from file(s) and interpolate in time
       !
       call SwashReadBndval ( curbfl%bfiled, curbfl%bctime, curbfl%bcloc, bndval )
       if (STPNOW()) return
       !
       if (.not.associated(curbfl%nextbfl)) exit
       curbfl => curbfl%nextbfl
       !
    enddo
    !
    ! determine boundary values at boundary faces of computational grid
    !
    if ( nbgrpt > 1 ) then
       !
       if ( ITEST >= 80 ) write (PRTEST,*) ' number of boundary faces ', nbgrpt-1
       !
       do i = 1, nbgrpt-1
          !
          indx  = int(bgridp(9*i-8),4)
          indxn = int(bgridp(9*i+1),4)
          !
          ! search for corresponding boundary face and cell
          !
          faceloop: do jb = 1, nfacesb
             j = jbface(jb)
             v(1) = face(j)%atti(FACEV1)
             v(2) = face(j)%atti(FACEV2)
             if ( (v(1) == indx .and. v(2) == indxn) .or. (v(2) == indx .and. v(1) == indxn) ) then
                indxf = j
                indxs = face(j)%atti(FACEC1)
                vm    = vmark(v(1))
                exit faceloop
             endif
          enddo faceloop
          !
          btype = int(bgridp(9*i-7),4)
          btype = abs(btype)
          !
          face(indxf)%atti(FBTYPE) = btype
          !
          fac  = face(indxf)%attr(FACEDISTC)
          rsgn = sign(1.,-fac)
          !
          nx = face(indxf)%attr(FACENORMX)
          ny = face(indxf)%attr(FACENORMY)
          !
          if ( .not. ny /= 0. ) then
             vdir = .true.
          else if ( .not. nx /= 0. ) then
             vdir = .false.
          endif
          !
          ! obtain boundary value from interpolation in space
          !
          wf1 = 0.001 * real(bgridp(9*i-6))
          k1  = int(bgridp(9*i-5),4)
          wf2 = 1. - wf1
          k2  = int(bgridp(9*i-3),4)
          !
          bval = wf1 * bndval(k1,1) + wf2 * bndval(k2,1)
          !
          ! if appropriate, ramp function is applied to prevent initially short large waves
          !
          tsmo = 0.00001 * real(bgridp(9*i-2))
          !
          if ( lrampf .and. tsmo /= 0. ) then
             !
             fsmo = .5 * ( 1. + tanh( timco/tsmo - 3. ) )
             !
          else
             !
             fsmo = 1.
             !
          endif
          !
          bval = fsmo * bval
          !
          ! get layer number and order of Stokes wave
          !
          itmp  = int(bgridp(9*i-1),4)
          if ( itmp < 0 ) then
             klay  = mod(itmp,10)
             istok = -itmp/10
          else
             klay  = itmp
             istok = 0
          endif
          !
          ! get coordinates of the midface with respect to first vertex of boundary
          !
          if ( it == 0 .and. klay == -1 ) then
             !
             kx = minloc(xcugrd(:), vmark(:)==vm)
             ky = minloc(ycugrd(:), vmark(:)==vm)
             !
             if ( kx(1) == ky(1) ) then
                x0 = xcugrd(kx(1))
                y0 = ycugrd(ky(1))
             else
                !
                d1 = sqrt((xcugrd(kx(1)))**2+(ycugrd(kx(1)))**2)
                d2 = sqrt((xcugrd(ky(1)))**2+(ycugrd(ky(1)))**2)
                !
                if ( d1 < d2 ) then
                   x0 = xcugrd(kx(1))
                   y0 = ycugrd(kx(1))
                else
                   x0 = xcugrd(ky(1))
                   y0 = ycugrd(ky(1))
                endif
                !
             endif
             !
             xp = face(indxf)%attr(FACEMX) - x0
             yp = face(indxf)%attr(FACEMY) - x0
             !
          endif
          !
          ! indicator for addition of bound long wave or linearized Riemann invariant
          !
          itmp = int(bgridp(9*i),4)
          if ( itmp /= 2 ) then
             ibnd  = itmp
             lriem = .false.
          else
             ibnd  = 0
             lriem = .true.
          endif
          !
          ! impose boundary condition for flow
          !
          if ( btype == 1 ) then
             !
             ! closed boundary
             !
             u1(indxf,:) = 0.
             !
          else if ( btype == 2 ) then
             !
             ! water level imposed
             !
             bcs(indxf) = bval
             !
             if ( wetu(indxf) == 0 ) then
                !
                ! dry boundary face
                !
                bcs(indxf) = max ( bval, -dps(indxs)+0.99*epsdry )
                !
             endif
             !
          else if ( btype == 3 ) then
             !
             ! flow velocity imposed
             !
             if ( klay == -2 ) then
                !
                ! logarithmic distribution in vertical
                !
                if ( hu(indxf) > 0. ) then
                   if ( irough == 4 ) then
                      z1 = 0.5 * ( zks(indxs,kmax-1) - zks(indxs,kmax) )
                      z0 = z1 / exp( vonkar*logfrc(indxs,2) )
                   else if ( cfricu(indxs) > 0. ) then
                      z0 = hu(indxf) / exp( vonkar/sqrt(cfricu(indxs)) + 1. )
                   endif
                   fac = bval * hu(indxf) / ( hu(indxf)*(log(hu(indxf)/z0)-1.) + z0 )
                   do k = 1, kmax
                      z2 = zku(indxf,k-1) - zku(indxf,kmax)
                      z1 = max(z0, zku(indxf,k) - zku(indxf,kmax) )
                      u1(indxf,k) = fac * ( z2*(log(z2/z0)-1.) - z1*(log(z1/z0)-1.) ) / ( z2 - z1 )
                   enddo
                endif
                !
             else if ( klay == -1 ) then
                !
                ! hyperbolic cosine distribution / Stokes wave solution in vertical
                !
                u1(indxf,:) = 0.
                !
                dep = hu(indxf)
                !
                if ( dep > 0. ) then
                   !
                   swd  = dps(indxs)
                   dcor = swd / dep
                   !
                   curbfs => fbfs
                   do
                      ibloc = curbfs%nbfs
                      if ( ibloc /= k1 .and. ibloc /= k2 ) then
                         if ( .not.associated(curbfs%nextbfs) ) exit
                         curbfs => curbfs%nextbfs
                         cycle
                      endif
                      !
                      if ( ibloc == k2 ) wf = wf2
                      if ( ibloc == k1 ) wf = wf1
                      if (    k1 == k2 ) wf = 1.
                      !
                      nfreq = curbfs%nfreq
                      wdir  = curbfs%spparm(3)
                      !
                      shape = abs(spshape(2))
                      !
                      ! determine free short wave components based on first order Fourier series
                      !
                      if ( it == 0 ) call SwashBCshortwave ( curbfs, nfreq, xp, yp, i, swd, wdir, rsgn, vdir, shape )
                      !
                      ! determine bound long wave using Hasselmann (1962) transfer functions (obsolete)
                      !
                      if ( ibnd == 1 .and. it == 0 ) call SwashBCboundwave ( curbfs, nfreq, xp, yp, i, swd, wdir, rsgn, vdir, shape )
                      !
                      ! determine first and second order Stokes waves
                      !
                      if ( istok /= 0 .and. it == 0 ) call SwashBCStokeswave ( curbfs, nfreq, i, swd, istok )
                      !
                      ! for irregular waves, compute sub- and super-harmonic transfer functions
                      !
                      if ( istok == 3 .and. it == 0 ) call SwashBCtransferfnc ( curbfs, nfreq, xp, yp, i, swd, wdir, rsgn, vdir, shape )
                      !
                      if ( kmax == 1 .or. ihydro == 3 ) then
                         !
                         if ( istok == 0 ) then
                            !
                            ! hyperbolic cosine distribution
                            !
                            fvalu = 0.
                            do j = 1, nfreq
                               !
                               omega = curbfs%omega(j)
                               kwav  = kwave(i,j)
                               !
                               if ( kwav /= 0. ) then
                                  fac = omega / ( kwav * dep )
                               else
                                  fac = sqrt( grav / dep )
                               endif
                               !
                               ! synthesize time series
                               !
                               fvalu = fvalu + fac * comp1(i,j) * cos( omega*timco ) + fac * comp2(i,j) * sin( omega*timco )
                               !
                            enddo
                            !
                         else if ( istok == 1 ) then
                            !
                            ! first order Stokes wave solution
                            !
                            fvalu = 0.
                            do j = 1, nfreq
                               !
                               omega = curbfs%omega(j)
                               !
                               ! synthesize time series
                               !
                               fvalu = fvalu + stku1c(i,j,1) * cos( omega*timco ) + stku1s(i,j,1) * sin( omega*timco )
                               !
                            enddo
                            !
                         else if ( istok > 1 ) then
                            !
                            ! second order Stokes wave solution
                            !
                            fvalu = 0.
                            do j = 1, nfreq
                               !
                               omega = curbfs%omega(j)
                               !
                               ! synthesize time series
                               !
                               fvalu = fvalu + stku1c(i,j,1) * cos(    omega*timco ) + stku1s(i,j,1) * sin(    omega*timco ) &
                                             + stku2c(i,j,1) * cos( 2.*omega*timco ) + stku2s(i,j,1) * sin( 2.*omega*timco )
                               !
                            enddo
                            !
                            ! add cross-interacting sub- and super-harmonics
                            ! (note: self-interacting super-harmonics are already included in the foregoing, e.g. stku2c, stkz2c)
                            !
                            if ( istok == 3 ) then
                               !
                               do j = 1, nfreq-1
                                  !
                                  omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                                  omegsp = curbfs%omega(j+1) + curbfs%omega(1)
                                  !
                                  fvalu = fvalu + subuc(i,j,1) * cos( omegsb*timco ) + subus(i,j,1) * sin( omegsb*timco ) &
                                                + supuc(i,j,1) * cos( omegsp*timco ) + supus(i,j,1) * sin( omegsp*timco )
                                  !
                               enddo
                               !
                            endif
                            !
                         endif
                         !
                         ! add bound long wave based on transfer functions of Hasselmann (obsolete)
                         ! note: mass flux should be divided by the water depth, but for reasons of robustness the still water depth is chosen
                         !
                         if ( ibnd == 1 ) then
                            !
                            do j = 1, nfreq-1
                               !
                               omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                               !
                               fvalu = fvalu + real ( fluxbu(i,j) / swd * exp( (0.,1.)*omegsb*timco ) )
                               !
                            enddo
                            !
                         endif
                         !
                         u1(indxf,:) = u1(indxf,:) + rsgn * wf * fsmo * fvalu
                         !
                      else
                         !
                         ! compute bound long wave (uniform distribution in vertical) based on transfer functions of Hasselmann (obsolete)
                         ! note: mass flux should be divided by the water depth, but for reasons of robustness the still water depth is chosen
                         !
                         fvalbu = 0.
                         !
                         if ( ibnd == 1 ) then
                            !
                            do j = 1, nfreq-1
                               !
                               omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                               !
                               fvalbu = fvalbu + real ( fluxbu(i,j) / swd * exp( (0.,1.)*omegsb*timco ) )
                               !
                            enddo
                            !
                         endif
                         !
                         if ( istok == 0 ) then
                            !
                            ! hyperbolic cosine distribution
                            !
                            do k = 1, kmax
                               !
                               z2 = zku(indxf,k-1) - zku(indxf,kmax)
                               z1 = zku(indxf,k  ) - zku(indxf,kmax)
                               !
                               fvalu = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  kwav  = kwave(i,j)
                                  !
                                  ! compute hyperbolic cosine distribution
                                  !
                                  fac1 = sinh(min(30.,kwav*dcor*z2)) - sinh(min(30.,kwav*dcor*z1))
                                  fac2 = kwav * (z2-z1) * sinh(min(30.,kwav*swd))
                                  !
                                  if ( fac2 /= 0. ) then
                                     fac = omega * fac1 / fac2
                                  else
                                     fac = sqrt( grav / dep )
                                  endif
                                  !
                                  ! synthesize time series
                                  !
                                  fvalu = fvalu + fac * comp1(i,j) * cos( omega*timco ) + fac * comp2(i,j) * sin( omega*timco )
                                  !
                               enddo
                               !
                               fvalu = fsmo * ( fvalu + fvalbu )
                               !
                               u1(indxf,k) = u1(indxf,k) + rsgn * wf * fvalu
                               !
                            enddo
                            !
                         else if ( istok == 1 ) then
                            !
                            ! first order Stokes wave solution
                            !
                            do k = 1, kmax
                               !
                               fvalu = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  !
                                  ! synthesize time series
                                  !
                                  fvalu = fvalu + stku1c(i,j,k) * cos( omega*timco ) + stku1s(i,j,k) * sin( omega*timco )
                                  !
                               enddo
                               !
                               fvalu = fsmo * ( fvalu + fvalbu )
                               !
                               u1(indxf,k) = u1(indxf,k) + rsgn * wf * fvalu
                               !
                            enddo
                            !
                         else if ( istok > 1 ) then
                            !
                            ! second order Stokes wave solution
                            !
                            do k = 1, kmax
                               !
                               fvalu = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  !
                                  ! synthesize time series
                                  !
                                  fvalu = fvalu + stku1c(i,j,k) * cos(    omega*timco ) + stku1s(i,j,k) * sin(    omega*timco ) &
                                                + stku2c(i,j,k) * cos( 2.*omega*timco ) + stku2s(i,j,k) * sin( 2.*omega*timco )
                                  !
                               enddo
                               !
                               ! add cross-interacting sub- and super-harmonics
                               ! (note: self-interacting super-harmonics are already included in the foregoing, e.g. stku2c, stkz2c)
                               !
                               if ( istok == 3 ) then
                                  !
                                  do j = 1, nfreq-1
                                     !
                                     omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                                     omegsp = curbfs%omega(j+1) + curbfs%omega(1)
                                     !
                                     fvalu = fvalu + subuc(i,j,k) * cos( omegsb*timco ) + subus(i,j,k) * sin( omegsb*timco ) &
                                                   + supuc(i,j,k) * cos( omegsp*timco ) + supus(i,j,k) * sin( omegsp*timco )
                                     !
                                  enddo
                                  !
                               endif
                               !
                               fvalu = fsmo * ( fvalu + fvalbu )
                               !
                               u1(indxf,k) = u1(indxf,k) + rsgn * wf * fvalu
                               !
                            enddo
                            !
                         endif
                         !
                      endif
                      !
                      if ( .not.associated(curbfs%nextbfs) ) exit
                      curbfs => curbfs%nextbfs
                      !
                   enddo
                   !
                endif
                !
             else if ( klay == 0 ) then
                !
                ! uniform distribution in vertical
                !
                u1(indxf,:) = bval
                !
             else
                !
                ! for each layer
                !
                u1(indxf,klay) = bval
                !
             endif
             !
          else if ( btype == 5 ) then
             !
             ! (specific) discharge imposed
             !
             if ( klay == -2 ) then
                !
                ! logarithmic distribution in vertical
                !
                if ( hu(indxf) > 0. ) then
                   if ( irough == 4 ) then
                      z1 = 0.5 * ( zks(indxs,kmax-1) - zks(indxs,kmax) )
                      z0 = z1 / exp( vonkar*logfrc(indxs,2) )
                   else if ( cfricu(indxs) > 0. ) then
                      z0 = hu(indxf) / exp( vonkar/sqrt(cfricu(indxs)) + 1. )
                   endif
                   fac = bval / ( hu(indxf)*(log(hu(indxf)/z0)-1.) + z0 )
                   do k = 1, kmax
                      z2 = zku(indxf,k-1) - zku(indxf,kmax)
                      z1 = max(z0, zku(indxf,k) - zku(indxf,kmax) )
                      u1(indxf,k) = fac * ( z2*(log(z2/z0)-1.) - z1*(log(z1/z0)-1.) ) / ( z2 - z1 )
                   enddo
                endif
                !
             else if ( klay == 0 ) then
                !
                ! uniform distribution in vertical
                !
                if ( hu(indxf) > 0. ) u1(indxf,:) = bval / hu(indxf)
                !
             else
                !
                ! for each layer
                !
                if ( hku(indxf,klay) > 0. ) u1(indxf,klay) = bval / hku(indxf,klay)
                !
             endif
             !
          else if ( btype == 6 ) then
             !
             ! Riemann invariant imposed (explicit approach)
             !
             if ( lriem ) then
                !
                ! linearized Riemann invariant
                !
                u1(indxf,:) = bval - rsgn * sqrt( grav/dps(indxs) ) * s0(indxs)
                !
             else
                !
                u1(indxf,:) = bval - 2.*rsgn * sqrt( grav * hu(indxf) )
                !
             endif
             !
          else if ( btype == 7 ) then
             !
             ! weakly reflective boundary condition imposed (explicit approach)
             !
             if ( klay == -2 ) then
                !
                ! logarithmic distribution in vertical
                !
                if ( hu(indxf) > 0. ) then
                   if ( irough == 4 ) then
                      z1 = 0.5 * ( zks(indxs,kmax-1) - zks(indxs,kmax) )
                      z0 = z1 / exp( vonkar*logfrc(indxs,2) )
                   else if ( cfricu(indxs) > 0. ) then
                      z0 = hu(indxf) / exp( vonkar/sqrt(cfricu(indxs)) + 1. )
                   endif
                   fac = rsgn * sqrt( grav / hu(indxf) ) * ( 2.*bval - s0(indxs) )
                   fac = fac * hu(indxf) / ( hu(indxf)*(log(hu(indxf)/z0)-1.) + z0 )
                   do k = 1, kmax
                      z2 = zku(indxf,k-1) - zku(indxf,kmax)
                      z1 = max(z0, zku(indxf,k) - zku(indxf,kmax) )
                      u1(indxf,k) = fac * ( z2*(log(z2/z0)-1.) - z1*(log(z1/z0)-1.) ) / ( z2 - z1 )
                   enddo
                endif
                !
             else if ( klay == -1 ) then
                !
                ! hyperbolic cosine distribution / Stokes wave solution in vertical
                !
                u1(indxf,:) = 0.
                !
                dep = hu(indxf)
                !
                if ( dep > 0. ) then
                   !
                   swd  = dps(indxs)
                   dcor = swd / dep
                   !
                   curbfs => fbfs
                   do
                      ibloc = curbfs%nbfs
                      if ( ibloc /= k1 .and. ibloc /= k2 ) then
                         if ( .not.associated(curbfs%nextbfs) ) exit
                         curbfs => curbfs%nextbfs
                         cycle
                      endif
                      !
                      if ( ibloc == k2 ) wf = wf2
                      if ( ibloc == k1 ) wf = wf1
                      if (    k1 == k2 ) wf = 1.
                      !
                      nfreq = curbfs%nfreq
                      wdir  = curbfs%spparm(3)
                      !
                      shape = abs(spshape(2))
                      !
                      ! determine free short wave components based on first order Fourier series
                      !
                      if ( it == 0 ) call SwashBCshortwave ( curbfs, nfreq, xp, yp, i, swd, wdir, rsgn, vdir, shape )
                      !
                      ! determine bound long wave using Hasselmann (1962) transfer functions (obsolete)
                      !
                      if ( ibnd == 1 .and. it == 0 ) call SwashBCboundwave ( curbfs, nfreq, xp, yp, i, swd, wdir, rsgn, vdir, shape )
                      !
                      ! determine first and second order Stokes waves
                      !
                      if ( istok /= 0 .and. it == 0 ) call SwashBCStokeswave ( curbfs, nfreq, i, swd, istok )
                      !
                      ! for irregular waves, compute sub- and super-harmonic transfer functions
                      !
                      if ( istok == 3 .and. it == 0 ) call SwashBCtransferfnc ( curbfs, nfreq, xp, yp, i, swd, wdir, rsgn, vdir, shape )
                      !
                      if ( kmax == 1 .or. ihydro == 3 ) then
                         !
                         if ( istok == 0 ) then
                            !
                            ! hyperbolic cosine distribution
                            !
                            fvalu = 0.
                            do j = 1, nfreq
                               !
                               omega = curbfs%omega(j)
                               kwav  = kwave(i,j)
                               !
                               if ( kwav /= 0. ) then
                                  fac = omega / ( kwav * dep ) + sqrt( grav / dep )
                               else
                                  fac = 2. * sqrt( grav / dep )
                               endif
                               !
                               ! synthesize time series
                               !
                               fvalu = fvalu + fac * comp1(i,j) * cos( omega*timco ) + fac * comp2(i,j) * sin( omega*timco )
                               !
                            enddo
                            !
                         else if ( istok == 1 ) then
                            !
                            ! first order Stokes wave solution
                            !
                            fac = sqrt( grav / dep )
                            !
                            fvalu = 0.
                            do j = 1, nfreq
                               !
                               omega = curbfs%omega(j)
                               !
                               ! synthesize time series
                               !
                               fvalu = fvalu + ( stku1c(i,j,1) + fac*stkz1c(i,j) ) * cos( omega*timco ) + ( stku1s(i,j,1) + fac*stkz1s(i,j) ) * sin( omega*timco )
                               !
                            enddo
                            !
                         else if ( istok > 1 ) then
                            !
                            ! second order Stokes wave solution
                            !
                            fac = sqrt( grav / dep )
                            !
                            fvalu = 0.
                            do j = 1, nfreq
                               !
                               omega = curbfs%omega(j)
                               !
                               ! synthesize time series
                               !
                               fvalu = fvalu + ( stku1c(i,j,1) + fac*stkz1c(i,j) ) * cos(    omega*timco ) + ( stku1s(i,j,1) + fac*stkz1s(i,j) ) * sin(    omega*timco ) &
                                             + ( stku2c(i,j,1) + fac*stkz2c(i,j) ) * cos( 2.*omega*timco ) + ( stku2s(i,j,1) + fac*stkz2s(i,j) ) * sin( 2.*omega*timco )
                               !
                            enddo
                            !
                            ! add cross-interacting sub- and super-harmonics
                            ! (note: self-interacting super-harmonics are already included in the foregoing, e.g. stku2c, stkz2c)
                            !
                            if ( istok == 3 ) then
                               !
                               do j = 1, nfreq-1
                                  !
                                  omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                                  omegsp = curbfs%omega(j+1) + curbfs%omega(1)
                                  !
                                  fvalu = fvalu + ( subuc(i,j,1) + fac*subzc(i,j) ) * cos( omegsb*timco ) + ( subus(i,j,1) + fac*subzs(i,j) ) * sin( omegsb*timco ) &
                                                + ( supuc(i,j,1) + fac*supzc(i,j) ) * cos( omegsp*timco ) + ( supus(i,j,1) + fac*supzs(i,j) ) * sin( omegsp*timco )
                                  !
                               enddo
                               !
                            endif
                            !
                         endif
                         !
                         ! add bound long wave based on transfer functions of Hasselmann (obsolete)
                         ! note: mass flux should be divided by the water depth, but for reasons of robustness the still water depth is chosen
                         !
                         if ( ibnd == 1 ) then
                            !
                            do j = 1, nfreq-1
                               !
                               omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                               !
                               fvalu = fvalu + real ( ( fluxbu(i,j) / swd + sqrt( grav / dep ) * zetab(i,j) ) * exp( (0.,1.)*omegsb*timco ) )
                               !
                            enddo
                            !
                         endif
                         !
                         u1(indxf,:) = u1(indxf,:) + rsgn * wf * ( fsmo * fvalu - sqrt( grav / hu(indxf) ) * s0(indxs) )
                         !
                      else
                         !
                         ! compute bound long wave (uniform distribution in vertical) based on transfer functions of Hasselmann (obsolete)
                         ! note: mass flux should be divided by the water depth, but for reasons of robustness the still water depth is chosen
                         !
                         fvalbu = 0.
                         !
                         if ( ibnd == 1 ) then
                            !
                            do j = 1, nfreq-1
                               !
                               omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                               !
                               fvalbu = fvalbu + real ( ( fluxbu(i,j) / swd + sqrt( grav / dep ) * zetab(i,j) ) * exp( (0.,1.)*omegsb*timco ) )
                               !
                            enddo
                            !
                         endif
                         !
                         if ( istok == 0 ) then
                            !
                            ! hyperbolic cosine distribution
                            !
                            do k = 1, kmax
                               !
                               z2 = zku(indxf,k-1) - zku(indxf,kmax)
                               z1 = zku(indxf,k  ) - zku(indxf,kmax)
                               !
                               fvalu = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  kwav  = kwave(i,j)
                                  !
                                  ! compute hyperbolic cosine distribution
                                  !
                                  fac1 = sinh(min(30.,kwav*dcor*z2)) - sinh(min(30.,kwav*dcor*z1))
                                  fac2 = kwav * (z2-z1) * sinh(min(30.,kwav*swd))
                                  !
                                  if ( fac2 /= 0. ) then
                                     fac = omega * fac1 / fac2 + sqrt( grav / dep )
                                  else
                                     fac = 2. * sqrt( grav / dep )
                                  endif
                                  !
                                  ! synthesize time series
                                  !
                                  fvalu = fvalu + fac * comp1(i,j) * cos( omega*timco ) + fac * comp2(i,j) * sin( omega*timco )
                                  !
                               enddo
                               !
                               fvalu = fsmo * ( fvalu + fvalbu )
                               !
                               u1(indxf,k) = u1(indxf,k) + rsgn * wf * ( fvalu - sqrt( grav / hu(indxf) ) * s0(indxs) )
                               !
                            enddo
                            !
                         else if ( istok == 1 ) then
                            !
                            ! first order Stokes wave solution
                            !
                            fac = sqrt( grav / dep )
                            !
                            do k = 1, kmax
                               !
                               fvalu = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  !
                                  ! synthesize time series
                                  !
                                  fvalu = fvalu + ( stku1c(i,j,k) + fac*stkz1c(i,j) ) * cos( omega*timco ) + ( stku1s(i,j,k) + fac*stkz1s(i,j) ) * sin( omega*timco )
                                  !
                               enddo
                               !
                               fvalu = fsmo * ( fvalu + fvalbu )
                               !
                               u1(indxf,k) = u1(indxf,k) + rsgn * wf * ( fvalu - sqrt( grav / hu(indxf) ) * s0(indxs) )
                               !
                            enddo
                            !
                         else if ( istok > 1 ) then
                            !
                            ! second order Stokes wave solution
                            !
                            fac = sqrt( grav / dep )
                            !
                            do k = 1, kmax
                               !
                               fvalu = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  !
                                  ! synthesize time series
                                  !
                                  fvalu = fvalu + ( stku1c(i,j,k) + fac*stkz1c(i,j) ) * cos(    omega*timco ) + ( stku1s(i,j,k) + fac*stkz1s(i,j) ) * sin(    omega*timco ) &
                                                + ( stku2c(i,j,k) + fac*stkz2c(i,j) ) * cos( 2.*omega*timco ) + ( stku2s(i,j,k) + fac*stkz2s(i,j) ) * sin( 2.*omega*timco )
                                  !
                               enddo
                               !
                               ! add cross-interacting sub- and super-harmonics
                               ! (note: self-interacting super-harmonics are already included in the foregoing, e.g. stku2c, stkz2c)
                               !
                               if ( istok == 3 ) then
                                  !
                                  do j = 1, nfreq-1
                                     !
                                     omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                                     omegsp = curbfs%omega(j+1) + curbfs%omega(1)
                                     !
                                     fvalu = fvalu + ( subuc(i,j,k) + fac*subzc(i,j) ) * cos( omegsb*timco ) + ( subus(i,j,k) + fac*subzs(i,j) ) * sin( omegsb*timco ) &
                                                   + ( supuc(i,j,k) + fac*supzc(i,j) ) * cos( omegsp*timco ) + ( supus(i,j,k) + fac*supzs(i,j) ) * sin( omegsp*timco )
                                     !
                                  enddo
                                  !
                               endif
                               !
                               fvalu = fsmo * ( fvalu + fvalbu )
                               !
                               u1(indxf,k) = u1(indxf,k) + rsgn * wf * ( fvalu - sqrt( grav / hu(indxf) ) * s0(indxs) )
                               !
                            enddo
                            !
                         endif
                         !
                      endif
                      !
                      if ( .not.associated(curbfs%nextbfs) ) exit
                      curbfs => curbfs%nextbfs
                      !
                   enddo
                   !
                endif
                !
             else if ( klay == 0 ) then
                !
                ! uniform distribution in vertical
                !
                if ( hu(indxf) > 0. ) u1(indxf,:) = rsgn * sqrt( grav / hu(indxf) ) * ( 2.*bval - s0(indxs) )
                !
             endif
             !
          else if ( btype == 8 ) then
             !
             ! Sommerfeld radiation condition applied to water level
             !
             rdx = face(indxf)%attr(FACEDISTC)
             !
             fac = abs(rdx) * dt * sqrt( grav * hu(indxf) )
             !
             bcs(indxf) = ( 1. - fac ) * bcso(indxf) + fac * s0(indxs)
             !
          else if ( btype == 10 ) then
             !
             ! outflow condition
             !
             ! implemented in routines SwashExpDepUflow, SwashImpDepUflow, SwashExpLayUflow, SwashImpLayUflow
             !
          else
             !
             write (msgstr,'(a,i2,a,i8)') 'unknown boundary type ',btype,' at face with index ',indxf
             call msgerr (2, trim(msgstr) )
             !
          endif
          !
          ! test output: flow variables in test points on boundary
          !
          if ( nptst > 0 ) then
             do j = 1, nptst
                ix  = xytst(j)
                lpb = indx == ix
                if ( lpb ) then
                   write (PRTEST,202) i, ix, wf1, k1, wf2, k2
                   write (PRTEST,203) s1(indxs), u1(indxf,1)
                endif
             enddo
          endif
          !
       enddo
       !
    endif
    !
    ! update source function by means of internal wave generation
    !
    if ( iwvgen /= 0 ) then
       !
       curbfs => fbfs
       do
          nfreq = curbfs%nfreq
          wdir  = curbfs%spparm(3)
          !
          ! if appropriate, ramp function is applied to prevent initially short large waves
          !
          tsmo = piwg(5)
          if ( lrampf .and. tsmo /= 0. ) then
             !
             fsmo = .5 * ( 1. + tanh( timco/tsmo - 3. ) )
             !
          else
             !
             fsmo = 1.
             !
          endif
          !
          shape = abs(spshape(2))
          !
          ! compute the source function amplitude and shape factor of source area
          !
          if ( it == 0 ) call SwashIntWavgen ( curbfs, nfreq, wdir, shape )
          !
          ! compute shortest distance to boundary
          !
          if ( it == 0 ) then
             !
             vm = nint(piwg(6))
             !
             do icell = 1, ncells
                !
                x0 = cell(icell)%attr(CELLCCX)
                y0 = cell(icell)%attr(CELLCCY)
                !
                dsrc(icell) = 100.*max(xclen,yclen)
                !
                do jb = 1, ncellsb
                   !
                   j = jbcell(jb)
                   !
                   do jj = 1, cell(j)%nof
                      !
                      ! face identifier
                      !
                      iface = cell(j)%face(jj)%atti(FACEID)
                      !
                      ! get vertices of present face
                      !
                      v(1) = face(iface)%atti(FACEV1)
                      v(2) = face(iface)%atti(FACEV2)
                      !
                      ! check boundary marker
                      !
                      if ( vmark(v(1)) == vm .and. vmark(v(2)) == vm ) then
                         !
                         xb = face(iface)%attr(FACEMX)
                         yb = face(iface)%attr(FACEMY)
                         !
                         dist = sqrt ( ( x0 - xb )**2 + ( y0 - yb )**2 )
                         !
                         if ( dist < dsrc(icell) ) dsrc(icell) = dist
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
          ! compute source function
          !
          cgsrc =       piwg(1)
          hwidt = 0.5 * piwg(2)
          !
          do icell = 1, ncells
             !
             if ( .not. dsrc(icell) < cgsrc - hwidt .and. .not. dsrc(icell) > cgsrc + hwidt ) then
                !
                sfval = 0.
                !
                do jj = 1, nfreq
                   !
                   ampl  = curbfs%sfamp(jj)
                   beta  = curbfs%bshap(jj)
                   !
                   omega = curbfs%omega(jj)
                   theta = curbfs%theta(jj)
                   kwav  = kwave(jj,1)
                   phase = curbfs%phase(jj)
                   !
                   sfval = sfval + ampl * exp( -beta * (cgsrc-dsrc(icell))**2 ) * ( cos( phase ) * cos( omega*timco ) + sin( phase ) * sin( omega*timco ) )
                   !
                enddo
                !
                srcm(icell) = fsmo * sfval
                !
             endif
             !
          enddo
          !
          if ( .not.associated(curbfs%nextbfs) ) exit
          curbfs => curbfs%nextbfs
          !
       enddo
       !
    endif
    !
    ! check if water level is prescribed at boundary faces
    !
    do j = 1, nfaces
       !
       if ( face(j)%atti(FBTYPE) == 2 ) then
          !
          wlimp(j) = .true.
          ! note: rdx in pressure gradient is reciprocal of the distance between boundary face and cell center
          !       here we take twice the distance in case of water level boundary
          pfac (j) = 0.5
          !
       else if ( face(j)%atti(FBTYPE) == 8 ) then
          !
          wlimp(j) = .true.
          pfac (j) = 1.
          !
       else
          !
          wlimp(j) = .false.
          pfac (j) = 0.
          !
       endif
       !
    enddo
    !
    ! update user-defined input fields, map onto computational grid and interpolate in time
    !
    !
    ! flow velocities
    !
    if ( ifldyn(2) == 1 ) then
       call SwashUpdateFld ( 2, 3, uxb   , uyb   , cosvc, sinvc, uxf   , uyf   , xcgrid, ycgrid, kgrpnt )
       if (STPNOW()) return
    endif
    !
    ! bottom friction coefficient
    !
    if ( ifldyn(4) == 1 ) then
       call SwashUpdateFld ( 4, 0, fric  , (/0./),    1.,    0., fricf , (/0./), xcgrid, ycgrid, kgrpnt )
       if (STPNOW()) return
    endif
    !
    ! wind velocities
    !
    if ( ifldyn(5) == 1 ) then
       call SwashUpdateFld ( 5, 6, wxi   , wyi   , coswc, sinwc, wxf   , wyf   , xcgrid, ycgrid, kgrpnt )
       if (STPNOW()) return
    endif
    !
    ! water level
    !
    if ( ifldyn(7) == 1 ) then
       call SwashUpdateFld ( 7, 0, wlevl , (/0./),    1.,    0., wlevf , (/0./), xcgrid, ycgrid, kgrpnt )
       if (STPNOW()) return
    endif
    !
    ! atmospheric pressure
    !
    if ( ifldyn(13) == 1 ) then
       call SwashUpdateFld (13, 0, pres  , (/0./),    1.,    0., presf , (/0./), xcgrid, ycgrid, kgrpnt )
       if (STPNOW()) return
    endif
    !
    ! update flow variables based on space varying input fields
    !
    call SwashUpdUFlowFlds
    if (STPNOW()) return
    !
    ! update bottom friction coefficient
    !
!TIMG    call SWTSTA(55)
    if ( irough == 4 .and. kmax > 1 ) then
       call SwashULogLaw
    else if ( it > 0 .and. (ifldyn(4) == 1 .or. irough == 3 .or. irough == 4) ) then
       call SwashUBotFrict ( u0 )
    endif
!TIMG    call SWTSTO(55)
    if (STPNOW()) return
    !
    ! calculate wind stresses, if appropriate
    !
!TIMG    call SWTSTA(56)
    if ( iwind /= 0 .and. (it == 0 .or. ifldyn(5) == 1 .or. relwnd .or. relwav) ) call SwashUWindStress
!TIMG    call SWTSTO(56)
    if (STPNOW()) return
    !
    ! update atmospheric pressure based on space varying input field and correct water level on open boundaries, if appropriate
    !
    if ( it == 0 ) then
       itmp       = ifldyn(13)
       ifldyn(13) = 1
    endif
    if ( svwp .and. (ifldyn(13) == 1 .or. prmean > 0.) ) call SwashUpdUPress
    if (STPNOW()) return
    if ( it == 0 ) ifldyn(13) = itmp
    !
    ! calculate horizontal eddy viscosity coefficient, if appropriate
    !
!TIMG    call SWTSTA(57)
    if ( ihvisc /= 0 ) then
       if ( kmax == 1 ) then
          call SwashUHorzVisc (   u0 )
       else
          call SwashUHorzVisc ( udep )
       endif
    endif
!TIMG    call SWTSTO(57)
    if (STPNOW()) return
    !
    ! calculate vertical eddy viscosity coefficient, if appropriate
    !
!TIMG    call SWTSTA(58)
    if ( iturb /= 0 ) call SwashVertVisc
!TIMG    call SWTSTO(58)
    !
    ! compute density, if appropriate
    !
    if ( idens /= 0 ) call SwashDensity
    !
 202 format (' boundary face', 2i8, 2(f8.3, i4))
 203 format (' wl, u: ', 2e12.4)
    !
end subroutine SwashUpdateUData
