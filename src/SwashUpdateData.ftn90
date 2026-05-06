subroutine SwashUpdateData ( it )
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
!    1.16: Dirk Rijnsdorp
!    3.02: Floris de Wit
!    6.01: Panagiotis Vasarmidis
!    6.02: Tom Bogaard
!    9.01: Dirk Rijnsdorp
!   11.01: Panagiotis Vasarmidis
!   11.02: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, February 2010: New subroutine
!    1.16,  October 2012: extension second order bound long wave correction
!    3.02,    April 2016: take into account heading per frequency in case of btype = 3 and 7
!    6.01,     June 2019: extension internal wave generation
!    6.02,     July 2019: extension 3D (non)linear k-eps model
!    9.01,  October 2022: extension moving rigid bodies
!   11.01,  October 2023: extension first and second order Stokes wave at boundaries
!   11.02, November 2023: extension second order transfer functions at boundaries
!
!   Purpose
!
!   Updates flow data, boundary conditions and input fields
!
!   Method
!
!   The following boundary conditions for flow can be imposed:
!
!   1) water level
!   2) velocity or discharge
!   3) Riemann invariant or Sommerfeld radiation or weakly reflective
!   4) outflow in case of supercritical flow
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
    use SwashRigBoddata
    use m_parall
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
                                      ! = 1; closed (not used here)
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
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer               :: ifpco    ! frequency with which preconditioner is applied
    integer               :: indx     ! point index
    integer               :: indxd    ! index of downward point
    integer               :: indxs    ! index of water level point
    integer               :: indxu    ! index of upward point
    integer               :: isrcb    ! start point index of source area (left/lower border)
    integer               :: isrce    ! end point index of source area (right/upper border)
    integer               :: istok    ! indicates the order of Stokes wave theory
                                      ! = 0; use hyperbolic cosine distribution for velocity (Airy wave theory)
                                      ! = 1; first order Stokes wave (Airy wave theory)
                                      ! = 2; second order Stokes wave
                                      ! = 3; second order sub- and super-harmonic transfer functions
    integer               :: itmp     ! temporary stored integer
    integer               :: ix       ! index of point in x-direction
    integer               :: iy       ! index of point in y-direction
    integer               :: j        ! loop counter
    integer               :: jj       ! loop counter
    integer               :: k        ! loop counter over vertical layers
    integer               :: k1       ! user-defined location of boundary point
    integer               :: k2       ! user-defined location of subsequent boundary point
    integer               :: klay     ! layer number
    integer               :: nfreq    ! number of components in Fourier series
    integer               :: shape    ! spectral shape
                                      ! = 1; Pierson Moskowitz
                                      ! = 2; Jonswap
                                      ! = 3; TMA
    !
    real                  :: ampl     ! amplitude of a Fourier component
    real                  :: azero    ! amplitude at zero frequency in Fourier series
    real                  :: beta     ! shape factor beta of source area
    real                  :: bval     ! actual boundary value
    real                  :: cgsrc    ! center of gravity of source area
    real                  :: cosbnd   ! cosine of inwarding normal to boundary side
    real                  :: coslat   ! cosine of latitude
    real                  :: dcor     ! correction in depth integration to get zero mass outflow
    real                  :: dep      ! local water depth
    real                  :: dx1      ! first component of covariant base vector a_(1)
    real                  :: dx2      ! first component of covariant base vector a_(2)
    real                  :: dy1      ! second component of covariant base vector a_(1)
    real                  :: dy2      ! second component of covariant base vector a_(2)
    real, save            :: epsab2   ! parameter to modify the weights of Adams-Bashforth scheme
    real                  :: fac      ! a factor
    real                  :: fac1     ! another factor
    real                  :: fac2     ! some other factor
    real                  :: fsmo     ! ramp function for smoothing incident waves
    real                  :: fval     ! value from a Fourier series
    real                  :: fval2    ! another value from a Fourier series
    real                  :: fvalu    ! value from a Fourier series for U-velocity component
    real                  :: fvalv    ! value from a Fourier series for V-velocity component
    real                  :: fvalb    ! value from bound long wave
    real                  :: fvalbu   ! value from bound long wave for U-velocity component
    real                  :: fvalbv   ! value from bound long wave for V-velocity component
    real                  :: hwidt    ! half width of source area
    real                  :: kwav     ! wave number
    real                  :: omega    ! angular frequency of a Fourier component
    real                  :: omegsb   ! angular frequency of bound sub-harmonic component
    real                  :: omegsp   ! angular frequency of bound super-harmonic component
    real                  :: phase    ! phase of a Fourier component
    real                  :: rsgn     ! sign for indicating in- and outflowing depending on boundary
                                      ! =+1; refers to inflowing at left and lower boundaries
                                      ! =-1; refers to outflowing at right and upper boundaries
    real                  :: s        ! sign to preserve symmetry at boundary
    real                  :: sinbnd   ! sine of inwarding normal to boundary side
    real                  :: sfval    ! value from an internal-generated series of wave components
    real                  :: srcb     ! start point of source area (left/lower border)
    real                  :: srce     ! end point of source area (right/upper border)
    real                  :: swd      ! still water depth
    real                  :: theta    ! randomly chosen directions assigned to each frequency
    real                  :: tsmo     ! period for smoothing boundary values during cold start
    real, dimension(kmax) :: utmp     ! temporary array to store Cartesian U-velocity component
    real, dimension(kmax) :: vtmp     ! temporary array to store Cartesian V-velocity component
    real                  :: wdir     ! incident or peak wave direction with respect to problem coordinates
    real                  :: wf       ! weighting factor
    real                  :: wf1      ! weighting factor for given boundary value to be interpolated in space
    real                  :: wf2      ! weighting factor for subsequent boundary value to be interpolated in space
    real                  :: xp       ! x-coordinate of grid point
    real                  :: yp       ! y-coordinate of grid point
    real                  :: z0       ! roughness height for logarithmic velocity profile
    real                  :: z1       ! interface at bottom of layer
    real                  :: z2       ! interface at top of layer
    !
    logical               :: isright  ! indicates right boundary is treated
    logical               :: lpb      ! indicates whether test point is a boundary point
    logical               :: lriem    ! indicates whether linearized Riemann invariant is imposed or not
    logical               :: lspon    ! indicates whether boundary point is in sponge layer
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
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdateData')
    !
    if ( (ihydro == 1 .or. ihydro == 2) .and. kmax > 1 ) then
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
       if ( .not.oned ) v0 = v1
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
             if ( lsubg ) then
                u0p = u1p
                if ( .not.oned ) v0p = v1p
                w0p = w1p
             endif
             !
          endif
          !
       endif
       !
       hso = hs
       if ( kmax > 1    ) hkso = hks
       if ( icreep /= 0 ) zkso = zks
       !
       ! moving rigid bodies
       if ( ifloat == 2 ) then
          !
          ! displacements, velocities and accelerations
          xfot0 = xfot1
          xfor0 = xfor1
          vfot0 = vfot1
          vfor0 = vfor1
          afot0 = afot1
          afor0 = afor1
          !
          ! hydrodynamic forces and torques
          fhyd0 = fhyd1
          thyd0 = thyd1
          !
          ! body-related forces and torques
          fbod0 = fbod1
          tbod0 = tbod1
          !
          ! extension mooring lines
          extml0 = extml1
          !
       endif
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
    if ( INODE == MASTER ) then
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
    endif
    !
    ! scatter array bndval to all nodes
    !
    call SWBROADC ( bndval, 2*nbval, SWREAL )
    if (STPNOW()) return
    !
    ! determine boundary values at boundary points of computational grid
    !
    if ( nbgrpt > 0 ) then
       !
       if ( ITEST >= 80 ) write (PRTEST,*) ' number of boundary points ', nbgrpt
       isright = .false.
       lspon   = .false.
       !
       do i = 1, nbgrpt
          !
          indx  = int(bgridp(9*i-8),4)
          btype = int(bgridp(9*i-7),4)
          !
          ! obtain indices of down- and upwards points next to boundary point
          !
          if ( oned ) then
             if ( indx == kgrpnt(1,1) ) then
                ! left boundary
                indxd  = indx
                indxu  = kgrpnt(2,1)
                indxs  = indxu
                rsgn   = +1.
                ibl(1) = btype
             else if ( indx == kgrpnt(mxc-1,1) ) then
                ! right boundary
                indxd  = kgrpnt(mxc  ,1)
                indxu  = kgrpnt(mxc-2,1)
                indxs  = indx
                rsgn   = -1.
                ibr(1) = btype
             endif
             !
             vdir = .true.
             xp   = 0.
             yp   = 0.
          else
             do j = 2, mxc-1
                if ( indx == kgrpnt(j,1) ) then
                   ! lower boundary
                   indxd  = indx
                   indxu  = kgrpnt(j,2)
                   indxs  = indxu
                   rsgn   = +1.
                   vdir   = .false.
                   cosbnd = -sinpc
                   sinbnd =  cospc
                   ibb(j) = btype
                   xp     = XGRDGL(j+MXF-1,1) - XGRDGL(2,1)
                   yp     = YGRDGL(j+MXF-1,1) - YGRDGL(2,1)
                   dx1    = xcgrid(j,1) - xcgrid(j-1,1)
                   dy1    = ycgrid(j,1) - ycgrid(j-1,1)
                   dx2    = xcgrid(j,2) - xcgrid(j  ,1)
                   dy2    = ycgrid(j,2) - ycgrid(j  ,1)
                   if ( kspher > 0 ) then
                      coslat = cos(degrad*(0.5*(ycgrid(j,1)+ycgrid(j-1,1)) + yoffs))
                      xp  = lendeg * coslat * xp
                      yp  = lendeg * yp
                      dx1 = lendeg * coslat * dx1
                      dy1 = lendeg * dy1
                      dx2 = lendeg * coslat * dx2
                      dy2 = lendeg * dy2
                   endif
                   ! check sponge layer at left/right
                   if ( sponxl(indxu) > 0. .or. sponxr(indxu) > 0. ) then
                      lspon = .true.
                   else
                      lspon = .false.
                   endif
                   goto 100
                else if ( j == mxc-1 .and. isright ) then
                   ! corner point (mxc-1,myc-1) belongs to the right boundary
                else if ( indx == kgrpnt(j,myc-1) ) then
                   ! upper boundary
                   indxd  = kgrpnt(j,myc  )
                   indxu  = kgrpnt(j,myc-2)
                   indxs  = indx
                   rsgn   = -1.
                   vdir   = .false.
                   cosbnd =  sinpc
                   sinbnd = -cospc
                   ibt(j) = btype
                   xp     = XGRDGL(j+MXF-1,myc-1) - XGRDGL(2,myc-1)
                   yp     = YGRDGL(j+MXF-1,myc-1) - YGRDGL(2,myc-1)
                   dx1    = xcgrid(j,myc-1) - xcgrid(j-1,myc-1)
                   dy1    = ycgrid(j,myc-1) - ycgrid(j-1,myc-1)
                   dx2    = xcgrid(j,myc  ) - xcgrid(j  ,myc-1)
                   dy2    = ycgrid(j,myc  ) - ycgrid(j  ,myc-1)
                   if ( kspher > 0 ) then
                      coslat = cos(degrad*(0.5*(ycgrid(j,myc-1)+ycgrid(j-1,myc-1)) + yoffs))
                      xp  = lendeg * coslat * xp
                      yp  = lendeg * yp
                      dx1 = lendeg * coslat * dx1
                      dy1 = lendeg * dy1
                      dx2 = lendeg * coslat * dx2
                      dy2 = lendeg * dy2
                   endif
                   ! check sponge layer at left/right
                   if ( sponxl(indxu) > 0. .or. sponxr(indxu) > 0. ) then
                      lspon = .true.
                   else
                      lspon = .false.
                   endif
                   goto 100
                endif
             enddo
             do j = 2, myc-1
                if ( indx == kgrpnt(1,j) ) then
                   ! left boundary
                   indxd  = indx
                   indxu  = kgrpnt(2,j)
                   indxs  = indxu
                   rsgn   = +1.
                   vdir   = .true.
                   cosbnd =  cospc
                   sinbnd =  sinpc
                   ibl(j) = btype
                   xp     = XGRDGL(1,j+MYF-1) - XGRDGL(1,2)
                   yp     = YGRDGL(1,j+MYF-1) - YGRDGL(1,2)
                   dx1    = xcgrid(2,j) - xcgrid(1,j  )
                   dy1    = ycgrid(2,j) - ycgrid(1,j  )
                   dx2    = xcgrid(1,j) - xcgrid(1,j-1)
                   dy2    = ycgrid(1,j) - ycgrid(1,j-1)
                   if ( kspher > 0 ) then
                      coslat = cos(degrad*(0.5*(ycgrid(1,j)+ycgrid(1,j-1)) + yoffs))
                      xp  = lendeg * coslat * xp
                      yp  = lendeg * yp
                      dx1 = lendeg * coslat * dx1
                      dy1 = lendeg * dy1
                      dx2 = lendeg * coslat * dx2
                      dy2 = lendeg * dy2
                   endif
                   ! check sponge layer at bottom/top
                   if ( sponyb(indxu) > 0. .or. sponyt(indxu) > 0. ) then
                      lspon = .true.
                   else
                      lspon = .false.
                   endif
                   goto 100
                else if ( indx == kgrpnt(mxc-1,j) ) then
                   ! right boundary
                   indxd  = kgrpnt(mxc  ,j)
                   indxu  = kgrpnt(mxc-2,j)
                   indxs  = indx
                   rsgn   = -1.
                   vdir   = .true.
                   cosbnd = -cospc
                   sinbnd = -sinpc
                   ibr(j) = btype
                   xp     = XGRDGL(mxc-1,j+MYF-1) - XGRDGL(mxc-1,2)
                   yp     = YGRDGL(mxc-1,j+MYF-1) - YGRDGL(mxc-1,2)
                   dx1    = xcgrid(mxc  ,j) - xcgrid(mxc-1,j  )
                   dy1    = ycgrid(mxc  ,j) - ycgrid(mxc-1,j  )
                   dx2    = xcgrid(mxc-1,j) - xcgrid(mxc-1,j-1)
                   dy2    = ycgrid(mxc-1,j) - ycgrid(mxc-1,j-1)
                   if ( kspher > 0 ) then
                      coslat = cos(degrad*(0.5*(ycgrid(mxc-1,j)+ycgrid(mxc-1,j-1)) + yoffs))
                      xp  = lendeg * coslat * xp
                      yp  = lendeg * yp
                      dx1 = lendeg * coslat * dx1
                      dy1 = lendeg * dy1
                      dx2 = lendeg * coslat * dx2
                      dy2 = lendeg * dy2
                   endif
                   isright= .true.
                   ! check sponge layer at bottom/top
                   if ( sponyb(indxu) > 0. .or. sponyt(indxu) > 0. ) then
                      lspon = .true.
                   else
                      lspon = .false.
                   endif
                   goto 100
                endif
             enddo
             cycle
          endif
          !
 100      continue
          !
          btype = abs(btype)
          !
          indxd = max(1,indxd)
          indxs = max(1,indxs)
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
          if ( btype == 2 ) then
             !
             ! water level imposed
             !
             s1(indxd) = bval
             !
             if ( wetu(indx) == 0 ) then
                !
                ! dry u-point
                !
                s1(indxd) = max ( bval, -dps(indxd)+0.99*epsdry )
                !
             endif
             !
             if ( .not.oned ) then
                if ( wetv(indx) == 0 ) then
                   !
                   ! dry v-point
                   !
                   s1(indxd) = max ( bval, -dps(indxd)+0.99*epsdry )
                   !
                endif
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
                if ( vdir ) then
                   if ( hu(indx) > 0. ) then
                      if ( irough == 4 ) then
                         z1 = 0.5 * ( zks(indx,kmax-1) - zks(indx,kmax) )
                         z0 = z1 / exp( vonkar*logfrc(indx,2) )
                      else if ( cfricu(indx) > 0. ) then
                         z0 = hu(indx) / exp( vonkar/sqrt(cfricu(indx)) + 1. )
                      endif
                      fac = bval * hu(indx) / ( hu(indx)*(log(hu(indx)/z0)-1.) + z0 )
                      do k = 1, kmax
                         z2 = zku(indx,k-1) - zku(indx,kmax)
                         z1 = max(z0, zku(indx,k) - zku(indx,kmax) )
                         u1(indx,k) = fac * ( z2*(log(z2/z0)-1.) - z1*(log(z1/z0)-1.) ) / ( z2 - z1 )
                      enddo
                   endif
                else
                   if ( hv(indx) > 0. ) then
                      if ( irough == 4 ) then
                         z1 = 0.5 * ( zks(indx,kmax-1) - zks(indx,kmax) )
                         z0 = z1 / exp( vonkar*logfrc(indx,2) )
                      else if ( cfricv(indx) > 0. ) then
                         z0 = hv(indx) / exp( vonkar/sqrt(cfricv(indx)) + 1. )
                      endif
                      fac = bval * hv(indx) / ( hv(indx)*(log(hv(indx)/z0)-1.) + z0 )
                      do k = 1, kmax
                         z2 = zkv(indx,k-1) - zkv(indx,kmax)
                         z1 = max(z0, zkv(indx,k) - zkv(indx,kmax) )
                         v1(indx,k) = fac * ( z2*(log(z2/z0)-1.) - z1*(log(z1/z0)-1.) ) / ( z2 - z1 )
                      enddo
                   endif
                endif
                !
             else if ( klay == -1 ) then
                !
                ! hyperbolic cosine distribution / Stokes wave solution in vertical
                !
                if ( oned ) then
                   u1(indx,:) = 0.
                else
                   utmp(:) = 0.
                   vtmp(:) = 0.
                endif
                !
                if ( vdir ) then
                   dep = hu(indx)
                else
                   dep = hv(indx)
                endif
                !
                if ( dep > 0. .and. .not.lspon ) then
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
                      ! if bichromatic wave, then determine incident angles with respect to normal on the boundary
                      !
                      if ( wdir == -888. .and. nfreq == 2 ) then
                         if ( it == 0 ) then
                            if ( vdir ) then
                               theta = curbfs%theta(1)
                               if ( theta /= -999. ) then
                                  curbfs%theta(1) = theta - 0.5*(rsgn-1.)*pi - alpc
                               else
                                  curbfs%theta(1) = 0.
                               endif
                               theta = curbfs%theta(2)
                               if ( theta /= -999. ) then
                                  curbfs%theta(2) = theta - 0.5*(rsgn-1.)*pi - alpc
                               else
                                  curbfs%theta(2) = 0.
                               endif
                            else
                               theta = curbfs%theta(1)
                               if ( theta /= -999. ) then
                                  curbfs%theta(1) = theta - 0.5*rsgn*pi - alpc
                               else
                                  curbfs%theta(1) = 0.
                               endif
                               theta = curbfs%theta(2)
                               if ( theta /= -999. ) then
                                  curbfs%theta(2) = theta - 0.5*rsgn*pi - alpc
                               else
                                  curbfs%theta(2) = 0.
                               endif
                            endif
                         endif
                         wdir = -999.
                      endif
                      !
                      if ( wdir == -999. ) then
                         ! incident or peak wave direction normal to the boundary
                         if ( vdir ) then
                            wdir = 0.5*(rsgn-1.)*pi + alpc
                         else
                            wdir = 0.5*rsgn*pi + alpc
                         endif
                      endif
                      !
                      ! direction of each wave component in line with incident direction
                      ! in order to preserve symmetry at boundaries
                      !
                      if ( vdir ) then
                         if ( .not. sin(wdir-alpc) < -0.0001 ) then
                            s = +1.
                         else
                            s = -1.
                         endif
                      else
                         if ( .not. cos(wdir-alpc) < -0.0001 ) then
                            s = +1.
                         else
                            s = -1.
                         endif
                      endif
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
                         if ( oned ) then
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
                            u1(indx,:) = u1(indx,:) + rsgn * wf * fsmo * fvalu
                            !
                         else
                            !
                            if ( istok == 0 ) then
                               !
                               ! hyperbolic cosine distribution
                               !
                               fvalu = 0.
                               fvalv = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  theta = curbfs%theta(j)
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
                                  fval = fac * comp1(i,j) * cos( omega*timco ) + fac * comp2(i,j) * sin( omega*timco )
                                  !
                                  fvalu = fvalu + fval * cos( wdir + s*theta )
                                  fvalv = fvalv + fval * sin( wdir + s*theta )
                                  !
                               enddo
                               !
                            else if ( istok == 1 ) then
                               !
                               ! first order Stokes wave solution
                               !
                               fvalu = 0.
                               fvalv = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  theta = curbfs%theta(j)
                                  !
                                  ! synthesize time series
                                  !
                                  fval = stku1c(i,j,1) * cos( omega*timco ) + stku1s(i,j,1) * sin( omega*timco )
                                  !
                                  fvalu = fvalu + fval * cos( wdir + s*theta )
                                  fvalv = fvalv + fval * sin( wdir + s*theta )
                                  !
                               enddo
                               !
                            else if ( istok > 1 ) then
                               !
                               ! second order Stokes wave solution
                               !
                               fvalu = 0.
                               fvalv = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  theta = curbfs%theta(j)
                                  !
                                  ! synthesize time series
                                  !
                                  fval = stku1c(i,j,1) * cos(    omega*timco ) + stku1s(i,j,1) * sin(    omega*timco ) &
                                       + stku2c(i,j,1) * cos( 2.*omega*timco ) + stku2s(i,j,1) * sin( 2.*omega*timco )
                                  !
                                  fvalu = fvalu + fval * cos( wdir + s*theta )
                                  fvalv = fvalv + fval * sin( wdir + s*theta )
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
                                     fvalv = fvalv + subvc(i,j,1) * cos( omegsb*timco ) + subvs(i,j,1) * sin( omegsb*timco ) &
                                                   + supvc(i,j,1) * cos( omegsp*timco ) + supvs(i,j,1) * sin( omegsp*timco )
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
                                  fvalv = fvalv + real ( fluxbv(i,j) / swd * exp( (0.,1.)*omegsb*timco ) )
                                  !
                               enddo
                               !
                            endif
                            !
                            utmp(:) = utmp(:) + wf * fsmo * fvalu
                            vtmp(:) = vtmp(:) + wf * fsmo * fvalv
                            !
                         endif
                         !
                      else
                         !
                         if ( oned ) then
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
                                  if ( vdir ) then
                                     z2 = zku(indx,k-1) - zku(indx,kmax)
                                     z1 = zku(indx,k  ) - zku(indx,kmax)
                                  else
                                     z2 = zkv(indx,k-1) - zkv(indx,kmax)
                                     z1 = zkv(indx,k  ) - zkv(indx,kmax)
                                  endif
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
                                  u1(indx,k) = u1(indx,k) + rsgn * wf * fvalu
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
                                  u1(indx,k) = u1(indx,k) + rsgn * wf * fvalu
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
                                  u1(indx,k) = u1(indx,k) + rsgn * wf * fvalu
                                  !
                               enddo
                               !
                            endif
                            !
                         else
                            !
                            ! compute bound long wave (uniform distribution in vertical) based on transfer functions of Hasselmann (obsolete)
                            ! note: mass flux should be divided by the water depth, but for reasons of robustness the still water depth is chosen
                            !
                            fvalbu = 0.
                            fvalbv = 0.
                            !
                            if ( ibnd == 1 ) then
                               !
                               do j = 1, nfreq-1
                                  !
                                  omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                                  !
                                  fvalbu = fvalbu + real ( fluxbu(i,j) / swd * exp( (0.,1.)*omegsb*timco ) )
                                  fvalbv = fvalbv + real ( fluxbv(i,j) / swd * exp( (0.,1.)*omegsb*timco ) )
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
                                  if ( vdir ) then
                                     z2 = zku(indx,k-1) - zku(indx,kmax)
                                     z1 = zku(indx,k  ) - zku(indx,kmax)
                                  else
                                     z2 = zkv(indx,k-1) - zkv(indx,kmax)
                                     z1 = zkv(indx,k  ) - zkv(indx,kmax)
                                  endif
                                  !
                                  fvalu = 0.
                                  fvalv = 0.
                                  do j = 1, nfreq
                                     !
                                     omega = curbfs%omega(j)
                                     theta = curbfs%theta(j)
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
                                     fval = fac * comp1(i,j) * cos( omega*timco ) + fac * comp2(i,j) * sin( omega*timco )
                                     !
                                     fvalu = fvalu + fval * cos( wdir + s*theta )
                                     fvalv = fvalv + fval * sin( wdir + s*theta )
                                     !
                                  enddo
                                  !
                                  fvalu = fsmo * ( fvalu + fvalbu )
                                  fvalv = fsmo * ( fvalv + fvalbv )
                                  !
                                  utmp(k) = utmp(k) + wf * fvalu
                                  vtmp(k) = vtmp(k) + wf * fvalv
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
                                  fvalv = 0.
                                  do j = 1, nfreq
                                     !
                                     omega = curbfs%omega(j)
                                     theta = curbfs%theta(j)
                                     !
                                     ! synthesize time series
                                     !
                                     fval = stku1c(i,j,k) * cos( omega*timco ) + stku1s(i,j,k) * sin( omega*timco )
                                     !
                                     fvalu = fvalu + fval * cos( wdir + s*theta )
                                     fvalv = fvalv + fval * sin( wdir + s*theta )
                                     !
                                  enddo
                                  !
                                  fvalu = fsmo * ( fvalu + fvalbu )
                                  fvalv = fsmo * ( fvalv + fvalbv )
                                  !
                                  utmp(k) = utmp(k) + wf * fvalu
                                  vtmp(k) = vtmp(k) + wf * fvalv
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
                                  fvalv = 0.
                                  do j = 1, nfreq
                                     !
                                     omega = curbfs%omega(j)
                                     theta = curbfs%theta(j)
                                     !
                                     ! synthesize time series
                                     !
                                     fval = stku1c(i,j,k) * cos(    omega*timco ) + stku1s(i,j,k) * sin(    omega*timco ) &
                                          + stku2c(i,j,k) * cos( 2.*omega*timco ) + stku2s(i,j,k) * sin( 2.*omega*timco )
                                     !
                                     fvalu = fvalu + fval * cos( wdir + s*theta )
                                     fvalv = fvalv + fval * sin( wdir + s*theta )
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
                                        fvalv = fvalv + subvc(i,j,k) * cos( omegsb*timco ) + subvs(i,j,k) * sin( omegsb*timco ) &
                                                      + supvc(i,j,k) * cos( omegsp*timco ) + supvs(i,j,k) * sin( omegsp*timco )
                                        !
                                     enddo
                                     !
                                  endif
                                  !
                                  fvalu = fsmo * ( fvalu + fvalbu )
                                  fvalv = fsmo * ( fvalv + fvalbv )
                                  !
                                  utmp(k) = utmp(k) + wf * fvalu
                                  vtmp(k) = vtmp(k) + wf * fvalv
                                  !
                               enddo
                               !
                            endif
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
                ! transform to grid-oriented velocity components
                !
                if ( .not.oned ) then
                   !
                   if ( cvleft ) then
                      if ( vdir ) then
                         u1(indx ,:) = ( dy2 * utmp(:) - dx2 * vtmp(:) ) / guu(indx )
                         v1(indxd,:) = ( dx1 * vtmp(:) - dy1 * utmp(:) ) / gvv(indxd)
                      else
                         u1(indxd,:) = ( dy2 * utmp(:) - dx2 * vtmp(:) ) / guu(indxd)
                         v1(indx ,:) = ( dx1 * vtmp(:) - dy1 * utmp(:) ) / gvv(indx )
                      endif
                   else
                      if ( vdir ) then
                         u1(indx ,:) = ( dx2 * vtmp(:) - dy2 * utmp(:) ) / guu(indx )
                         v1(indxd,:) = ( dy1 * utmp(:) - dx1 * vtmp(:) ) / gvv(indxd)
                      else
                         u1(indxd,:) = ( dx2 * vtmp(:) - dy2 * utmp(:) ) / guu(indxd)
                         v1(indx ,:) = ( dy1 * utmp(:) - dx1 * vtmp(:) ) / gvv(indx )
                      endif
                   endif
                   !
                endif
                !
             else if ( klay == 0 ) then
                !
                ! uniform distribution in vertical
                !
                if ( vdir ) then
                   u1(indx,:) = bval
                else
                   v1(indx,:) = bval
                endif
                !
             else
                !
                ! for each layer
                !
                if ( vdir ) then
                   u1(indx,klay) = bval
                else
                   v1(indx,klay) = bval
                endif
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
                if ( vdir ) then
                   if ( hu(indx) > 0. ) then
                      if ( irough == 4 ) then
                         z1 = 0.5 * ( zks(indx,kmax-1) - zks(indx,kmax) )
                         z0 = z1 / exp( vonkar*logfrc(indx,2) )
                      else if ( cfricu(indx) > 0. ) then
                         z0 = hu(indx) / exp( vonkar/sqrt(cfricu(indx)) + 1. )
                      endif
                      fac = bval / ( hu(indx)*(log(hu(indx)/z0)-1.) + z0 )
                      do k = 1, kmax
                         z2 = zku(indx,k-1) - zku(indx,kmax)
                         z1 = max(z0, zku(indx,k) - zku(indx,kmax) )
                         u1(indx,k) = fac * ( z2*(log(z2/z0)-1.) - z1*(log(z1/z0)-1.) ) / ( z2 - z1 )
                      enddo
                   endif
                else
                   if ( hv(indx) > 0. ) then
                      if ( irough == 4 ) then
                         z1 = 0.5 * ( zks(indx,kmax-1) - zks(indx,kmax) )
                         z0 = z1 / exp( vonkar*logfrc(indx,2) )
                      else if ( cfricv(indx) > 0. ) then
                         z0 = hv(indx) / exp( vonkar/sqrt(cfricv(indx)) + 1. )
                      endif
                      fac = bval / ( hv(indx)*(log(hv(indx)/z0)-1.) + z0 )
                      do k = 1, kmax
                         z2 = zkv(indx,k-1) - zkv(indx,kmax)
                         z1 = max(z0, zkv(indx,k) - zkv(indx,kmax) )
                         v1(indx,k) = fac * ( z2*(log(z2/z0)-1.) - z1*(log(z1/z0)-1.) ) / ( z2 - z1 )
                      enddo
                   endif
                endif
                !
             else if ( klay == 0 ) then
                !
                ! uniform distribution in vertical
                !
                if ( vdir ) then
                   if ( hu(indx) > 0. ) u1(indx,:) = bval / hu(indx)
                else
                   if ( hv(indx) > 0. ) v1(indx,:) = bval / hv(indx)
                endif
                !
             else
                !
                ! for each layer
                !
                if ( vdir ) then
                   if ( hku(indx,klay) > 0. ) u1(indx,klay) = bval / hku(indx,klay)
                else
                   if ( hkv(indx,klay) > 0. ) v1(indx,klay) = bval / hkv(indx,klay)
                endif
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
                if ( vdir ) then
                   u1(indx,:) = bval - rsgn * sqrt( grav/dps(indxd) ) * s0(indxs)
                else
                   v1(indx,:) = bval - rsgn * sqrt( grav/dps(indxd) ) * s0(indxs)
                endif
                !
             else
                !
                if ( vdir ) then
                   u1(indx,:) = bval - 2.*rsgn * sqrt( grav * hu(indx) )
                else
                   v1(indx,:) = bval - 2.*rsgn * sqrt( grav * hv(indx) )
                endif
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
                if ( vdir ) then
                   if ( hu(indx) > 0. ) then
                      if ( irough == 4 ) then
                         z1 = 0.5 * ( zks(indx,kmax-1) - zks(indx,kmax) )
                         z0 = z1 / exp( vonkar*logfrc(indx,2) )
                      else if ( cfricu(indx) > 0. ) then
                         z0 = hu(indx) / exp( vonkar/sqrt(cfricu(indx)) + 1. )
                      endif
                      fac = rsgn * sqrt( grav / hu(indx) ) * ( 2.*bval - s0(indxs) )
                      fac = fac * hu(indx) / ( hu(indx)*(log(hu(indx)/z0)-1.) + z0 )
                      do k = 1, kmax
                         z2 = zku(indx,k-1) - zku(indx,kmax)
                         z1 = max(z0, zku(indx,k) - zku(indx,kmax) )
                         u1(indx,k) = fac * ( z2*(log(z2/z0)-1.) - z1*(log(z1/z0)-1.) ) / ( z2 - z1 )
                      enddo
                   endif
                else
                   if ( hv(indx) > 0. ) then
                      if ( irough == 4 ) then
                         z1 = 0.5 * ( zks(indx,kmax-1) - zks(indx,kmax) )
                         z0 = z1 / exp( vonkar*logfrc(indx,2) )
                      else if ( cfricv(indx) > 0. ) then
                         z0 = hv(indx) / exp( vonkar/sqrt(cfricv(indx)) + 1. )
                      endif
                      fac = rsgn * sqrt( grav / hv(indx) ) * ( 2.*bval - s0(indxs) )
                      fac = fac * hv(indx) / ( hv(indx)*(log(hv(indx)/z0)-1.) + z0 )
                      do k = 1, kmax
                         z2 = zkv(indx,k-1) - zkv(indx,kmax)
                         z1 = max(z0, zkv(indx,k) - zkv(indx,kmax) )
                         v1(indx,k) = fac * ( z2*(log(z2/z0)-1.) - z1*(log(z1/z0)-1.) ) / ( z2 - z1 )
                      enddo
                   endif
                endif
                !
             else if ( klay == -1 ) then
                !
                ! hyperbolic cosine distribution / Stokes wave solution in vertical
                !
                if ( oned ) then
                   u1(indx,:) = 0.
                else
                   utmp(:) = 0.
                   vtmp(:) = 0.
                endif
                !
                if ( vdir ) then
                   dep = hu(indx)
                else
                   dep = hv(indx)
                endif
                !
                if ( dep > 0. .and. .not.lspon ) then
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
                      ! if bichromatic wave, then determine incident angles with respect to normal on the boundary
                      !
                      if ( wdir == -888. .and. nfreq == 2 ) then
                         if ( it == 0 ) then
                            if ( vdir ) then
                               theta = curbfs%theta(1)
                               if ( theta /= -999. ) then
                                  curbfs%theta(1) = theta - 0.5*(rsgn-1.)*pi - alpc
                               else
                                  curbfs%theta(1) = 0.
                               endif
                               theta = curbfs%theta(2)
                               if ( theta /= -999. ) then
                                  curbfs%theta(2) = theta - 0.5*(rsgn-1.)*pi - alpc
                               else
                                  curbfs%theta(2) = 0.
                               endif
                            else
                               theta = curbfs%theta(1)
                               if ( theta /= -999. ) then
                                  curbfs%theta(1) = theta - 0.5*rsgn*pi - alpc
                               else
                                  curbfs%theta(1) = 0.
                               endif
                               theta = curbfs%theta(2)
                               if ( theta /= -999. ) then
                                  curbfs%theta(2) = theta - 0.5*rsgn*pi - alpc
                               else
                                  curbfs%theta(2) = 0.
                               endif
                            endif
                         endif
                         wdir = -999.
                      endif
                      !
                      if ( wdir == -999. ) then
                         ! incident or peak wave direction normal to the boundary
                         if ( vdir ) then
                            wdir = 0.5*(rsgn-1.)*pi + alpc
                         else
                            wdir = 0.5*rsgn*pi + alpc
                         endif
                      endif
                      !
                      ! direction of each wave component in line with incident direction
                      ! in order to preserve symmetry at boundaries
                      !
                      if ( vdir ) then
                         if ( .not. sin(wdir-alpc) < -0.0001 ) then
                            s = +1.
                         else
                            s = -1.
                         endif
                      else
                         if ( .not. cos(wdir-alpc) < -0.0001 ) then
                            s = +1.
                         else
                            s = -1.
                         endif
                      endif
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
                         if ( oned ) then
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
                            u1(indx,:) = u1(indx,:) + rsgn * wf * ( fsmo * fvalu - sqrt( grav / hu(indx) ) * s0(indxs) )
                            !
                         else
                            !
                            if ( istok == 0 ) then
                               !
                               ! hyperbolic cosine distribution
                               !
                               fvalu = 0.
                               fvalv = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  theta = curbfs%theta(j)
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
                                  fval = comp1(i,j) * cos( omega*timco ) + comp2(i,j) * sin( omega*timco )                       ! surface elevation
                                  !
                                  fvalu = fvalu + fac * cos( wdir + s*theta ) * fval + sqrt( grav / dep ) * fval * cosbnd
                                  fvalv = fvalv + fac * sin( wdir + s*theta ) * fval + sqrt( grav / dep ) * fval * sinbnd
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
                               fvalv = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  theta = curbfs%theta(j)
                                  !
                                  ! synthesize time series
                                  !
                                  fval  = stku1c(i,j,1) * cos( omega*timco ) + stku1s(i,j,1) * sin( omega*timco )
                                  fval2 = stkz1c(i,j  ) * cos( omega*timco ) + stkz1s(i,j  ) * sin( omega*timco )
                                  !
                                  fvalu = fvalu + fval * cos( wdir + s*theta ) + fac * fval2 * cosbnd
                                  fvalv = fvalv + fval * sin( wdir + s*theta ) + fac * fval2 * sinbnd
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
                               fvalv = 0.
                               do j = 1, nfreq
                                  !
                                  omega = curbfs%omega(j)
                                  theta = curbfs%theta(j)
                                  !
                                  ! synthesize time series
                                  !
                                  fval  = stku1c(i,j,1) * cos(    omega*timco ) + stku1s(i,j,1) * sin(    omega*timco ) &
                                        + stku2c(i,j,1) * cos( 2.*omega*timco ) + stku2s(i,j,1) * sin( 2.*omega*timco )
                                  fval2 = stkz1c(i,j  ) * cos(    omega*timco ) + stkz1s(i,j  ) * sin(    omega*timco ) &
                                        + stkz2c(i,j  ) * cos( 2.*omega*timco ) + stkz2s(i,j  ) * sin( 2.*omega*timco )
                                  !
                                  fvalu = fvalu + fval * cos( wdir + s*theta ) + fac * fval2 * cosbnd
                                  fvalv = fvalv + fval * sin( wdir + s*theta ) + fac * fval2 * sinbnd
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
                                     fvalu = fvalu + ( subuc(i,j,1) + fac*subzc(i,j) * cosbnd ) * cos( omegsb*timco ) + ( subus(i,j,1) + fac*subzs(i,j) * cosbnd ) * sin( omegsb*timco ) &
                                                   + ( supuc(i,j,1) + fac*supzc(i,j) * cosbnd ) * cos( omegsp*timco ) + ( supus(i,j,1) + fac*supzs(i,j) * cosbnd ) * sin( omegsp*timco )
                                     !
                                     fvalv = fvalv + ( subvc(i,j,1) + fac*subzc(i,j) * sinbnd ) * cos( omegsb*timco ) + ( subvs(i,j,1) + fac*subzs(i,j) * sinbnd ) * sin( omegsb*timco ) &
                                                   + ( supvc(i,j,1) + fac*supzc(i,j) * sinbnd ) * cos( omegsp*timco ) + ( supvs(i,j,1) + fac*supzs(i,j) * sinbnd ) * sin( omegsp*timco )
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
                                  fvalu = fvalu + real ( ( fluxbu(i,j) / swd + sqrt( grav / dep ) * zetab(i,j) * cosbnd ) * exp( (0.,1.)*omegsb*timco ) )
                                  fvalv = fvalv + real ( ( fluxbv(i,j) / swd + sqrt( grav / dep ) * zetab(i,j) * sinbnd ) * exp( (0.,1.)*omegsb*timco ) )
                                  !
                               enddo
                               !
                            endif
                            !
                            utmp(:) = utmp(:) + wf * ( fsmo * fvalu - sqrt( grav / dep ) * s0(indxs) * cosbnd )
                            vtmp(:) = vtmp(:) + wf * ( fsmo * fvalv - sqrt( grav / dep ) * s0(indxs) * sinbnd )
                            !
                         endif
                         !
                      else
                         !
                         if ( oned ) then
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
                                  if ( vdir ) then
                                     z2 = zku(indx,k-1) - zku(indx,kmax)
                                     z1 = zku(indx,k  ) - zku(indx,kmax)
                                  else
                                     z2 = zkv(indx,k-1) - zkv(indx,kmax)
                                     z1 = zkv(indx,k  ) - zkv(indx,kmax)
                                  endif
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
                                  u1(indx,k) = u1(indx,k) + rsgn * wf * ( fvalu - sqrt( grav / hu(indx) ) * s0(indxs) )
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
                                  u1(indx,k) = u1(indx,k) + rsgn * wf * ( fvalu - sqrt( grav / hu(indx) ) * s0(indxs) )
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
                                  u1(indx,k) = u1(indx,k) + rsgn * wf * ( fvalu - sqrt( grav / hu(indx) ) * s0(indxs) )
                                  !
                               enddo
                               !
                            endif
                            !
                         else
                            !
                            ! compute bound long wave (uniform distribution in vertical) based on transfer functions of Hasselmann (obsolete)
                            ! note: mass flux should be divided by the water depth, but for reasons of robustness the still water depth is chosen
                            !
                            fvalbu = 0.
                            fvalbv = 0.
                            !
                            if ( ibnd == 1 ) then
                               !
                               do j = 1, nfreq-1
                                  !
                                  omegsb = curbfs%omega(j+1) - curbfs%omega(1)
                                  !
                                  fvalbu = fvalbu + real ( ( fluxbu(i,j) / swd + sqrt( grav / dep ) * zetab(i,j) * cosbnd ) * exp( (0.,1.)*omegsb*timco ) )
                                  fvalbv = fvalbv + real ( ( fluxbv(i,j) / swd + sqrt( grav / dep ) * zetab(i,j) * sinbnd ) * exp( (0.,1.)*omegsb*timco ) )
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
                                  if ( vdir ) then
                                     z2 = zku(indx,k-1) - zku(indx,kmax)
                                     z1 = zku(indx,k  ) - zku(indx,kmax)
                                  else
                                     z2 = zkv(indx,k-1) - zkv(indx,kmax)
                                     z1 = zkv(indx,k  ) - zkv(indx,kmax)
                                  endif
                                  !
                                  fvalu = 0.
                                  fvalv = 0.
                                  do j = 1, nfreq
                                     !
                                     omega = curbfs%omega(j)
                                     theta = curbfs%theta(j)
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
                                     fval = comp1(i,j) * cos( omega*timco ) + comp2(i,j) * sin( omega*timco )                       ! surface elevation
                                     !
                                     fvalu = fvalu + fac * cos( wdir + s*theta ) * fval + sqrt( grav / dep ) * fval * cosbnd
                                     fvalv = fvalv + fac * sin( wdir + s*theta ) * fval + sqrt( grav / dep ) * fval * sinbnd
                                     !
                                  enddo
                                  !
                                  fvalu = fsmo * ( fvalu + fvalbu )
                                  fvalv = fsmo * ( fvalv + fvalbv )
                                  !
                                  utmp(k) = utmp(k) + wf * ( fvalu - sqrt( grav / dep ) * s0(indxs) * cosbnd )
                                  vtmp(k) = vtmp(k) + wf * ( fvalv - sqrt( grav / dep ) * s0(indxs) * sinbnd )
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
                                  fvalv = 0.
                                  do j = 1, nfreq
                                     !
                                     omega = curbfs%omega(j)
                                     theta = curbfs%theta(j)
                                     !
                                     ! synthesize time series
                                     !
                                     fval  = stku1c(i,j,k) * cos( omega*timco ) + stku1s(i,j,k) * sin( omega*timco )
                                     fval2 = stkz1c(i,j  ) * cos( omega*timco ) + stkz1s(i,j  ) * sin( omega*timco )
                                     !
                                     fvalu = fvalu + fval * cos( wdir + s*theta ) + fac * fval2 * cosbnd
                                     fvalv = fvalv + fval * sin( wdir + s*theta ) + fac * fval2 * sinbnd
                                     !
                                  enddo
                                  !
                                  fvalu = fsmo * ( fvalu + fvalbu )
                                  fvalv = fsmo * ( fvalv + fvalbv )
                                  !
                                  utmp(k) = utmp(k) + wf * ( fvalu - sqrt( grav / dep ) * s0(indxs) * cosbnd )
                                  vtmp(k) = vtmp(k) + wf * ( fvalv - sqrt( grav / dep ) * s0(indxs) * sinbnd )
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
                                  fvalv = 0.
                                  do j = 1, nfreq
                                     !
                                     omega = curbfs%omega(j)
                                     theta = curbfs%theta(j)
                                     !
                                     ! synthesize time series
                                     !
                                     fval  = stku1c(i,j,k) * cos(    omega*timco ) + stku1s(i,j,k) * sin(    omega*timco ) &
                                           + stku2c(i,j,k) * cos( 2.*omega*timco ) + stku2s(i,j,k) * sin( 2.*omega*timco )
                                     fval2 = stkz1c(i,j  ) * cos(    omega*timco ) + stkz1s(i,j  ) * sin(    omega*timco ) &
                                           + stkz2c(i,j  ) * cos( 2.*omega*timco ) + stkz2s(i,j  ) * sin( 2.*omega*timco )
                                     !
                                     fvalu = fvalu + fval * cos( wdir + s*theta ) + fac * fval2 * cosbnd
                                     fvalv = fvalv + fval * sin( wdir + s*theta ) + fac * fval2 * sinbnd
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
                                        fvalu = fvalu + ( subuc(i,j,k) + fac*subzc(i,j) * cosbnd ) * cos( omegsb*timco ) + ( subus(i,j,k) + fac*subzs(i,j) * cosbnd ) * sin( omegsb*timco ) &
                                                      + ( supuc(i,j,k) + fac*supzc(i,j) * cosbnd ) * cos( omegsp*timco ) + ( supus(i,j,k) + fac*supzs(i,j) * cosbnd ) * sin( omegsp*timco )
                                        !
                                        fvalv = fvalv + ( subvc(i,j,k) + fac*subzc(i,j) * sinbnd ) * cos( omegsb*timco ) + ( subvs(i,j,k) + fac*subzs(i,j) * sinbnd ) * sin( omegsb*timco ) &
                                                      + ( supvc(i,j,k) + fac*supzc(i,j) * sinbnd ) * cos( omegsp*timco ) + ( supvs(i,j,k) + fac*supzs(i,j) * sinbnd ) * sin( omegsp*timco )
                                        !
                                     enddo
                                     !
                                  endif
                                  !
                                  fvalu = fsmo * ( fvalu + fvalbu )
                                  fvalv = fsmo * ( fvalv + fvalbv )
                                  !
                                  utmp(k) = utmp(k) + wf * ( fvalu - sqrt( grav / dep ) * s0(indxs) * cosbnd )
                                  vtmp(k) = vtmp(k) + wf * ( fvalv - sqrt( grav / dep ) * s0(indxs) * sinbnd )
                                  !
                               enddo
                               !
                            endif
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
                ! transform to grid-oriented velocity components
                !
                if ( .not.oned ) then
                   !
                   if ( cvleft ) then
                      if ( vdir ) then
                         u1(indx ,:) = ( dy2 * utmp(:) - dx2 * vtmp(:) ) / guu(indx )
                         v1(indxd,:) = ( dx1 * vtmp(:) - dy1 * utmp(:) ) / gvv(indxd)
                      else
                         u1(indxd,:) = ( dy2 * utmp(:) - dx2 * vtmp(:) ) / guu(indxd)
                         v1(indx ,:) = ( dx1 * vtmp(:) - dy1 * utmp(:) ) / gvv(indx )
                      endif
                   else
                      if ( vdir ) then
                         u1(indx ,:) = ( dx2 * vtmp(:) - dy2 * utmp(:) ) / guu(indx )
                         v1(indxd,:) = ( dy1 * utmp(:) - dx1 * vtmp(:) ) / gvv(indxd)
                      else
                         u1(indxd,:) = ( dx2 * vtmp(:) - dy2 * utmp(:) ) / guu(indxd)
                         v1(indx ,:) = ( dy1 * utmp(:) - dx1 * vtmp(:) ) / gvv(indx )
                      endif
                   endif
                   !
                endif
                !
             else if ( klay == 0 ) then
                !
                ! uniform distribution in vertical
                !
                if ( vdir ) then
                   if ( hu(indx) > 0. ) u1(indx,:) = rsgn * sqrt( grav / hu(indx) ) * ( 2.*bval - s0(indxs) )
                else
                   if ( hv(indx) > 0. ) v1(indx,:) = rsgn * sqrt( grav / hv(indx) ) * ( 2.*bval - s0(indxs) )
                endif
                !
             endif
             !
          else if ( btype == 8 ) then
             !
             ! Sommerfeld radiation condition applied to flow velocity
             !
             if ( oned ) then
                !
                fac = dt * sqrt( grav * hu(indx) ) / dx
                !
                u1(indx,:) = ( 1. - fac ) * u0(indx,:) + fac * u0(indxu,:)
                !
             else
                !
                if ( vdir ) then
                   !
                   fac = dt * sqrt( grav * hu(indx) ) / gvu(indx)
                   !
                   u1(indx,:) = ( 1. - fac ) * u0(indx,:) + fac * u0(indxu,:)
                   !
                else
                   !
                   fac = dt * sqrt( grav * hv(indx) ) / guv(indx)
                   !
                   v1(indx,:) = ( 1. - fac ) * v0(indx,:) + fac * v0(indxu,:)
                   !
                endif
                !
             endif
             !
          else if ( btype == 10 ) then
             !
             ! outflow condition
             !
             ! implemented in routines SwashExpDep1DHflow, SwashImpDep1DHflow, SwashExpDep2DHflow and SwashImpDep2DHflow
             !
          else
             !
             write (msgstr,'(a,i2,a,i4)') 'unknown boundary type ',btype,' at point with index ',indx
             call msgerr (2, trim(msgstr) )
             !
          endif
          !
          ! test output: flow variables in test points on boundary
          !
          if ( nptst > 0 ) then
             do j = 1, nptst
                ix  = xytst(2*j-1)
                iy  = xytst(2*j  )
                lpb = indx == kgrpnt(ix,iy)
                if ( lpb ) then
                   write (PRTEST,201) i, ix, iy, wf1, k1, wf2, k2
                   if ( oned ) then
                      write (PRTEST, 203) s1(indx), u1(indx,1)
                   else
                      write (PRTEST, 204) s1(indx), u1(indx,1), v1(indx,1)
                   endif
                endif
             enddo
          endif
          !
       enddo
       !
    endif
    !
    ! at closed boundaries, set velocities to zero
    !
    if ( oned ) then
       !
       if ( ibl(1) == 1 .and. LMXF ) u1(kgrpnt(1    ,1),:) = 0.
       if ( ibr(1) == 1 .and. LMXL ) u1(kgrpnt(mxc-1,1),:) = 0.
       !
    else
       !
       if ( .not.lreptx ) then
          !
          do i = 2, myc-1
             !
             if ( ibl(i) == 1 .and. LMXF ) u1(kgrpnt(1    ,i),:) = 0.
             if ( ibr(i) == 1 .and. LMXL ) u1(kgrpnt(mxc-1,i),:) = 0.
             !
          enddo
          !
       endif
       !
       if ( .not.lrepty ) then
          !
          do i = 2, mxc-1
             !
             if ( ibb(i) == 1 .and. LMYF ) v1(kgrpnt(i,1    ),:) = 0.
             if ( ibt(i) == 1 .and. LMYL ) v1(kgrpnt(i,myc-1),:) = 0.
             !
          enddo
          !
       endif
       !
    endif
    !
    ! after update boundary conditions, exchange water levels and velocities with neighbouring subdomains (if appropriate)
    !
    call SWEXCHG ( s1, kgrpnt, 1, 1    )
    call SWEXCHG ( u1, kgrpnt, 1, kmax )
    if ( .not.oned ) call SWEXCHG ( v1, kgrpnt, 1, kmax )
    if (STPNOW()) return
    !
    ! after update boundary conditions, synchronize water levels and velocities at appropriate boundaries in case of repeating grid (if appropriate)
    !
    call periodic ( s1, kgrpnt, 1, 1    )
    call periodic ( u1, kgrpnt, 1, kmax )
    if ( .not.oned ) call periodic ( v1, kgrpnt, 1, kmax )
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
          if ( wdir == -999. ) then
             ! incident or peak wave direction normal to the boundary
             if ( lsrcfy ) then
                wdir = alpc
             else if ( lsrcfx ) then
                wdir = 0.5*pi + alpc
             endif
          endif
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
          ! compute source function
          !
          if ( oned ) then
             !
             cgsrc = piwg(1) - dble(xoffs)
             hwidt = 0.5 * piwg(2)
             !
             isrcb = nint(cgsrc/dx - hwidt/dx - real(MXF) + 101.) - 99
             isrce = nint(cgsrc/dx + hwidt/dx - real(MXF) + 101.) - 99
             !
             do i = 2, mxc-1
                !
                indx = kgrpnt(i,1)
                !
                if ( .not. i < isrcb .and. .not. i > isrce ) then
                   !
                   xp = real(i-isrcb) * dx - hwidt
                   !
                   sfval = 0.
                   !
                   do jj = 1, nfreq
                      !
                      ampl  = curbfs%sfamp(jj)
                      beta  = curbfs%bshap(jj)
                      !
                      omega = curbfs%omega(jj)
                      phase = curbfs%phase(jj)
                      !
                      sfval = sfval + ampl * exp( -beta * xp**2 ) * cos( omega * timco + phase )
                      !
                   enddo
                   !
                   srcm(indx) = fsmo * sfval
                   !
                endif
                !
             enddo
             !
          else
             !
             if ( optg == 1 ) then
                !
                ! grid domain is rectilinear and uniform
                !
                if ( lsrcfy ) then
                   !
                   cgsrc = piwg(1) - dble(xoffs)
                   hwidt = 0.5 * piwg(2)
                   !
                   isrcb = nint(cgsrc/dx - hwidt/dx - real(MXF) + 101.) - 99
                   isrce = nint(cgsrc/dx + hwidt/dx - real(MXF) + 101.) - 99
                   !
                   do j = 2, myc-1
                      do i = 2, mxc-1
                         !
                         indx = kgrpnt(i,j)
                         !
                         if ( .not. i < isrcb .and. .not. i > isrce ) then
                            !
                            xp = real(i-isrcb) * dx - hwidt
                            yp = real(j+MYF-2) * dy
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
                               sfval = sfval + ampl * exp( -beta * xp**2 ) * cos( kwav * yp * sin(theta) + omega * timco + phase )
                               !
                            enddo
                            !
                            srcm(indx) = fsmo * sfval
                            !
                         endif
                         !
                      enddo
                   enddo
                   !
                else if ( lsrcfx ) then
                   !
                   cgsrc = piwg(1) - dble(yoffs)
                   hwidt = 0.5 * piwg(2)
                   !
                   isrcb = nint(cgsrc/dy - hwidt/dy - real(MYF) + 101.) - 99
                   isrce = nint(cgsrc/dy + hwidt/dy - real(MYF) + 101.) - 99
                   !
                   do i = 2, mxc-1
                      do j = 2, myc-1
                         !
                         indx = kgrpnt(i,j)
                         !
                         if ( .not. j < isrcb .and. .not. j > isrce ) then
                            !
                            yp = real(j-isrcb) * dy - hwidt
                            xp = real(i+MXF-2) * dx
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
                               sfval = sfval + ampl * exp( -beta * yp**2 ) * cos( kwav * xp * cos(theta) + omega * timco + phase )
                               !
                            enddo
                            !
                            srcm(indx) = fsmo * sfval
                            !
                         endif
                         !
                      enddo
                   enddo
                   !
                endif
                !
             else if ( optg == 3 ) then
                !
                ! domain is rectilinear but grid may not be uniform
                !
                if ( lsrcfy ) then
                   !
                   cgsrc = piwg(1) - dble(xoffs)
                   hwidt = 0.5 * piwg(2)
                   !
                   srcb = cgsrc - hwidt
                   srce = cgsrc + hwidt
                   !
                   do j = 2, myc-1
                      do i = 2, mxc-1
                         !
                         indx = kgrpnt(i,j)
                         !
                         xp = 0.5 * ( xcgrid(i,j) + xcgrid(i-1,j) )
                         yp = 0.5 * ( ycgrid(i,j) + ycgrid(i,j-1) )
                         !
                         if ( .not. xp < srcb .and. .not. xp > srce ) then
                            !
                            xp = xp - cgsrc
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
                               sfval = sfval + ampl * exp( -beta * xp**2 ) * cos( kwav * yp * sin(theta) + omega * timco + phase )
                               !
                            enddo
                            !
                            srcm(indx) = fsmo * sfval
                            !
                         endif
                         !
                      enddo
                   enddo
                   !
                else if ( lsrcfx ) then
                   !
                   cgsrc = piwg(1) - dble(yoffs)
                   hwidt = 0.5 * piwg(2)
                   !
                   srcb = cgsrc - hwidt
                   srce = cgsrc + hwidt
                   !
                   do i = 2, mxc-1
                      do j = 2, myc-1
                         !
                         indx = kgrpnt(i,j)
                         !
                         xp = 0.5 * ( xcgrid(i,j) + xcgrid(i-1,j) )
                         yp = 0.5 * ( ycgrid(i,j) + ycgrid(i,j-1) )
                         !
                         if ( .not. yp < srcb .and. .not. yp > srce ) then
                            !
                            yp = yp - cgsrc
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
                               sfval = sfval + ampl * exp( -beta * yp**2 ) * cos( kwav * xp * cos(theta) + omega * timco + phase )
                               !
                            enddo
                            !
                            srcm(indx) = fsmo * sfval
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
          endif
          !
          if ( .not.associated(curbfs%nextbfs) ) exit
          curbfs => curbfs%nextbfs
          !
       enddo
       !
    endif
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
    call SwashUpdFlowFlds
    if (STPNOW()) return
    !
    ! update bottom friction coefficient
    !
!TIMG    call SWTSTA(55)
    if ( irough == 4 .and. kmax > 1 ) then
       call SwashLogLaw
    else if ( it > 0 .and. (ifldyn(4) == 1 .or. irough == 3 .or. irough == 4) ) then
       call SwashBotFrict ( u0, v0 )
    endif
!TIMG    call SWTSTO(55)
    if (STPNOW()) return
    !
    ! calculate wind stresses, if appropriate
    !
!TIMG    call SWTSTA(56)
    if ( iwind /= 0 .and. (it == 0 .or. ifldyn(5) == 1 .or. relwnd .or. relwav) ) call SwashWindStress
!TIMG    call SWTSTO(56)
    if (STPNOW()) return
    !
    ! update atmospheric pressure based on space varying input field and correct water level on open boundaries, if appropriate
    !
    if ( it == 0 ) then
       itmp       = ifldyn(13)
       ifldyn(13) = 1
    endif
    if ( svwp .and. (ifldyn(13) == 1 .or. prmean > 0.) ) call SwashUpdPress
    if (STPNOW()) return
    if ( it == 0 ) ifldyn(13) = itmp
    !
    ! calculate horizontal eddy viscosity coefficient, if appropriate
    !
!TIMG    call SWTSTA(57)
    if ( ihvisc > 1 ) then
       if ( kmax == 1 ) then
          call SwashHorzVisc (   u0,   v0 )
       else
          call SwashHorzVisc ( udep, vdep )
       endif
       if (STPNOW()) return
    endif
!TIMG    call SWTSTO(57)
    !
    ! calculate vertical eddy viscosity coefficient, if appropriate
    !
!TIMG    call SWTSTA(58)
    if ( iturb /= 0 ) call SwashVertVisc
!TIMG    call SWTSTO(58)
    !
    ! calculate the (non)linear Reynolds stress tensor, if appropriate
    !
    if ( iturb > 1 ) call SwashReynoldsStress
    !
    ! adapt frictional forces inside porous structures, if appropriate
    !
    if ( iporos /= 0 ) then
       if ( kmax == 1 ) then
          call SwashPorFricDep
       else
          call SwashPorFricLay
       endif
    endif
    !
    ! compute vegetation coefficients, if appropriate
    !
    if ( iveg /= 0 ) call SwashVeget
    !
    ! compute density, if appropriate
    !
    if ( idens /= 0 ) call SwashDensity
    !
 201 format (' boundary point', 3i8, 2(f8.3, i4))
 202 format (' boundary vertex', 2i8, 2(f8.3, i4))
 203 format (' wl, u: ', 2e12.4)
 204 format (' wl, u, v: ', 3e12.4)
    !
end subroutine SwashUpdateData
