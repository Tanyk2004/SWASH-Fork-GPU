subroutine SwashCheckPrep
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
!    6.01: Panagiotis Vasarmidis
!    6.02: Tom Bogaard
!    7.01: Marcel Zijlema
!    9.01: Marcel Zijlema
!   10.05: Dirk Rijnsdorp
!
!   Updates
!
!    1.00,  February 2010: New subroutine
!    1.05,   January 2012: extension logarithmic wall-law
!    1.10,     March 2012: extension porous media
!    1.21,      July 2013: extension transport constituents
!    1.22,    August 2013: extension vegetation
!    1.31,      July 2014: extension sediment transport
!    4.01, September 2016: subgrid approach
!    4.05,   January 2017: extension floating objects
!    6.01,      June 2019: extension internal wave generation
!    6.02,      July 2019: extension 3D (non)linear k-eps model
!    7.01,   January 2020: extension to unstructured grids
!    9.01,   October 2022: extension moving rigid bodies
!   10.05,      July 2023: extension ambient current
!
!   Purpose
!
!   Checks inconsistencies in settings and changes if necessary
!   Does some preparations before computation is started
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use SwashFlowdata
    use SwashSolvedata
    use SwashRigBoddata
    use m_genarr
    use outp_data
    use m_bndspec
    use m_parall
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Parameter variables
!
    real, parameter                             :: eps = 1.e-3 ! a small criterion number
!
!   Local variables
!
    integer                                     :: i           ! loop counter
    integer                                     :: icell       ! loop counter over cells
    integer                                     :: icelll      ! left cell of present face
    integer                                     :: icellr      ! right cell of present face
    integer, save                               :: ient = 0    ! number of entries in this subroutine
    integer                                     :: iface       ! loop counter over faces
    integer                                     :: indx        ! grid point index
    integer                                     :: indxo       ! grid point index in own subdomain
    integer                                     :: indxr       ! index of point right from point of consideration
    integer                                     :: indxu       ! index of point up from point of consideration
    integer                                     :: istat       ! indicate status of allocation
    integer                                     :: istok       ! indicates the order of Stokes wave theory
                                                               ! = 0; use hyperbolic cosine distribution for velocity (Airy wave theory)
                                                               ! = 1; first order Stokes wave (Airy wave theory)
                                                               ! = 2; second order Stokes wave
                                                               ! = 3; second order sub- and super-harmonic transfer functions
    integer                                     :: itmp        ! temporary stored integer
    integer                                     :: ix          ! loop counter in x-direction of subdomain
    integer                                     :: ixgl        ! counter in x-direction of global domain
    integer                                     :: iy          ! loop counter in y-direction of subdomain
    integer                                     :: iygl        ! counter in y-direction of global domain
    integer                                     :: j           ! (loop) counter
    integer                                     :: jlay        ! layer numbers counted so far
    integer                                     :: jlow        ! lower bound of input array
    integer                                     :: jupp        ! upper bound of input array
    integer                                     :: k           ! counter
    integer                                     :: l           ! counter
    integer                                     :: m           ! body index
    integer                                     :: nfen        ! number of fenders
    integer                                     :: nfreq       ! number of frequencies
    integer                                     :: nmli        ! number of mooring lines
    integer                                     :: sumlay      ! sum of layer numbers
    integer, dimension(3)                       :: vc          ! vertices of present cell
    !
    real                                        :: altmp       ! an auxiliary real to temporarily store a direction
    real                                        :: cf          ! local bottom friction coefficient / Coriolis parameter
    real                                        :: cosva       ! =cos(-alpg(20)),
                                                               ! cosine of angle of ambient velocity input grid w.r.t. computational grid
    real                                        :: dep         ! local bottom level
    real                                        :: draft       ! local draft of floating object
    real                                        :: dsiz        ! characteristic grain size
    real                                        :: fac         ! a factor
    real                                        :: hlayp       ! layer thickness of pressure layer
    real                                        :: hps         ! porous structure height
    real                                        :: labl        ! label of rigid body
    real*8                                      :: lambda      ! polynomial function of temperature and salinity (for sea water) in Eckart's formula
    real                                        :: n           ! volumetric porosity
    real                                        :: npl         ! local number of plants per square meter
    real                                        :: omega       ! Earth's angular velocity
    real                                        :: pp          ! local atmospheric pressure
    real*8                                      :: pf0         ! another polynomial function of temperature and salinity (for sea water) in Eckart's formula
    real                                        :: rhoinf      ! spectral radius for infinite frequency to control generalized-alpha scheme
    real                                        :: sinva       ! =sine(-alpg(20)),
                                                               ! sine of angle of ambient velocity input grid w.r.t. computational grid
    real                                        :: ss          ! salinity/sediment
    real                                        :: sumh        ! sum of the layer thicknesses
    real                                        :: SVALQI      ! interpolated value of an input array for point given in problem coordinates
    real                                        :: tt          ! temperature
    real                                        :: uu          ! local u-velocity component
    real                                        :: vv          ! local v-velocity component
    real                                        :: wlvl        ! local water level
    real                                        :: xcl         ! length of computational domain in x-direction
    real                                        :: xp          ! x-coordinate of grid point
    real                                        :: yc          ! y-coordinate of cell center
    real                                        :: ycl         ! length of computational domain in y-direction
    real                                        :: yp          ! y-coordinate of grid point
    !
    real, dimension(:), allocatable             :: tempxl      ! temporary array to store global left sponge layer in x-direction
    real, dimension(:), allocatable             :: tempxr      ! temporary array to store global right sponge layer in x-direction
    real, dimension(:), allocatable             :: tempyb      ! temporary array to store global lower sponge layer in y-direction
    real, dimension(:), allocatable             :: tempyt      ! temporary array to store global upper sponge layer in y-direction
    !
    logical                                     :: bperc       ! indicates whether at least one layer thickness is variable
    logical                                     :: EQREAL      ! compares two reals
    logical                                     :: logbnd      ! true if logarithmic velocity profile is specified at boundary
    logical                                     :: nonequi     ! true if the layers are non-equidistant distributed
    logical                                     :: STPNOW      ! indicates that program must stop
    !
    character(80)                               :: msgstr      ! string to pass message
    !
    type(bfsdat), pointer                       :: curbfs      ! current item in list of Fourier series parameters
    type(bgpdat), pointer                       :: curbgp      ! current item in list of boundary grid points
    type(ribdat), pointer                       :: currib      ! current item in list of rigid body parameters
    !
    type(celltype), dimension(:), pointer       :: cell        ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer       :: face        ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashCheckPrep')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    if ( mxc <= 1   .and. optg /= 5   ) call msgerr (3, 'no valid computational grid; check command CGRID')
    if ( mcgrd <= 1 .and. nverts <= 0 ) call msgerr (3, 'no valid computational grid; check command READ BOT or READ UNSTRU')
    !
    ! set periodicity for incident wave boundary condition
    !
    if ( lreptx ) then
       bcperx = .true.
    else if ( lrepty ) then
       bcpery = .true.
    endif
    !
    ! check repeating grid for parallel runs
    !
    if ( LORB .and. ( lreptx .or. lrepty ) ) call msgerr (3,'repeating grid not allowed in case of orthogonal recursive bisection')
    if ( .not.LMXF .and. .not.LMXL ) lreptx = .false.
    if ( .not.LMYF .and. .not.LMYL ) lrepty = .false.
    !
    ! check combination of REPeating option and grid type and dimension
    !
    if ( lreptx ) then
       if ( oned ) call msgerr (3, 'repeating grid not allowed in 1D mode')
       if ( optg == 5 ) call msgerr (3, 'unstructured grid cannot be REPeating')
       if ( ( pnums(6) == 1. .and. pnums(16) == 1. .and. pnums(46) == 1. ) .and. mxc < 3 ) call msgerr (3, 'MXC must be larger or equal 1 for REPeating option')
       if ( ( pnums(6) /= 1. .or.  pnums(16) /= 1. .or.  pnums(46) /= 1. ) .and. mxc < 4 ) call msgerr (3, 'MXC must be larger or equal 2 for REPeating option')
    else if ( lrepty ) then
       if ( optg == 5 ) call msgerr (3, 'unstructured grid cannot be REPeating')
       if ( ( pnums(6) == 1. .and. pnums(16) == 1. .and. pnums(46) == 1. ) .and. myc < 3 ) call msgerr (3, 'MYC must be larger or equal 1 for REPeating option')
       if ( ( pnums(6) /= 1. .or.  pnums(16) /= 1. .or.  pnums(46) /= 1. ) .and. myc < 4 ) call msgerr (3, 'MYC must be larger or equal 2 for REPeating option')
    endif
    !
    ! check grid partitioning in case of 2DH depth-averaged, non-hydrostatic mode
    !
    if ( LORB .and. .not.oned .and. kmax == 1 .and. ihydro /= 0 ) call msgerr (3, 'stripwise grid partitioning is required in case of 2DH depth-averaged, non-hydrostatic mode')
    !
    ! determine bounds for loop counting in computational grid
    !
    !   Basic rules for loop counting:
    !
    !    wl-points running from mfu to ml
    !     u-points running from mf  to ml
    !
    mf  = 1
    mfu = mf  + 1
    mlu = mxc
    ml  = mlu - 1
    !
    if ( .not.oned ) then
       !
       nf  = 1
       nfu = nf  + 1
       nlu = myc
       nl  = nlu - 1
       !
    endif
    !
    ! for parallel runs, however, at a subdomain interface at the start of the loop not coinciding with a real boundary, the
    ! u-point mf is owned by the neighbouring subdomain and hence, should not be (re)calculated in this subdomain. In effect,
    ! this can be realised by setting mf = 3 assuming that the subdomain is extended with a halo area of 3 cells width.
    ! At a subdomain interface at the end of the loop not coinciding with a real boundary, the u-point ml is owned by the
    ! current subdomain and should be treated as an internal point. This can be effected by setting ml = mxc - 2, assuming
    ! a halosize of 3. As a consequence, the calculation of the water level and non-hydrostatic pressure at point ml is
    ! redundant. The same holds for the v- and wl-points in y-direction.
    !
    if ( .not.LMXF ) then
       mf  = IHALOX
       mfu = mf + 1
    endif
    if ( .not.LMYF ) then
       nf  = IHALOY
       nfu = nf + 1
    endif
    if ( .not.LMXL ) then
       ml  = mxc - IHALOX + 1
       mlu = ml + 1
    endif
    if ( .not.LMYL ) then
       nl  = myc - IHALOY + 1
       nlu = nl + 1
    endif
    !
    ! if no vertical layer schematisation specified for depth-averaged case, specify it
    !
    if ( .not.allocated(hlay) ) then
       allocate(hlay(kmax))
       hlay = 1.
    endif
    if ( .not.allocated(indlay) ) then
       allocate(indlay(kmax))
       indlay = 1
    endif
    !
    if ( leds(1) == 0 ) call msgerr (3, 'bottom grid not defined')
    if ( leds(1) == 1 ) call msgerr (3, 'no bottom levels read')
    !
    if ( leds(2) == 2 ) then
       if ( leds(3) /= 2 ) call msgerr (3, 'VY not read, while VX has been read')
       altmp = - alpg(2)
       cosvc = cos(altmp)
       sinvc = sin(altmp)
    endif
    !
    if ( leds(4) == 2 ) varfr = .true.
    !
    if ( leds(5) == 2 ) then
       if ( leds(6) /= 2 ) call msgerr (3, 'WY not read, while WX has been read')
       varwi = .true.
       if ( iwind == 0 ) iwind = 1
       altmp = - alpg(5)
       coswc = cos(altmp)
       sinwc = sin(altmp)
    endif
    !
    if ( iporos == 1 ) then
       if ( leds(10) == 0 ) call msgerr (3, 'porosity grid not defined')
       if ( leds(10) == 1 ) call msgerr (3, 'no porosity read')
       if ( leds(11) == 2 ) vargs = .true.
       if ( leds(12) == 2 ) varsh = .true.
    endif
    !
    if ( ifloat == 2 ) then
       if ( leds(19) == 2 ) then
          if ( leds(18) /= 2 ) call msgerr (3, 'DRAFT not read, while LABEL has been read')
       endif
    endif
    !
    if ( leds(20) == 2 ) then
       if ( leds(21) /= 2 ) call msgerr (3, 'AVY not read, while AVX has been read')
       if ( ifllay(20) /= kpmax .and. ifllay(20) /= 1 ) call msgerr (3, 'input field of ambient current has incorrect number of layers')
       varva = .true.
       if ( iamb == 0 ) iamb = 2
       altmp = - alpg(20)
       cosva = cos(altmp)
       sinva = sin(altmp)
    endif
    !
    if ( leds(22) == 2 ) then
       if ( leds(20) /= 2 ) call msgerr (3, 'ACURRENT not read, while MWL has been read')
       varea = .true.
    endif
    !
    if ( leds(13) == 2 ) svwp = .true.
    !
    if ( leds(14) == 2 ) varnpl = .true.
    !
    ! check consistency with respect to space varying wind and pressure
    !
    if ( svwp .and. .not.varwi ) call msgerr (3, 'space varying wind not read, while atmospheric pressure has been read')
    !
    ! check consistency with respect to space varying vegetation density
    !
    if ( iveg == 0 .and. varnpl ) then
       call msgerr (1, 'space varying vegetation density has been read but will be ignored')
       varnpl = .false.
    endif
    !
    ! check if there are instationary flow input fields while shallow water equations need to be solved
    !
    if ( .not.momskip ) then
       if ( ifldyn(2) /= 0 ) then
          call msgerr (1, 'instationary VX field has been read but will be ignored')
          ifldyn(2) = 0
       endif
       if ( ifldyn(3) /= 0 ) then
          call msgerr (1, 'instationary VY field has been read but will be ignored')
          ifldyn(3) = 0
       endif
       if ( ifldyn(7) /= 0 ) then
          call msgerr (1, 'instationary WLEV field has been read but will be ignored')
          ifldyn(7) = 0
       endif
    endif
    !
    ! check different parameters when momentum equations are skipped
    !
    if ( momskip ) then
       if ( itrans == 0 .or. ltrans /= lsed ) call msgerr (3, 'transport not specified, while momentum equations have been skipped')
       if ( kmax > 1 ) call msgerr (3, 'layer-averaged mode not allowed when momentum equations are skipped')
       lmixt  = .false.
       rhos   = 1.
       ihydro = 0
       irough = 0
       iwind  = 0
       isurf  = 0
       iturb  = 0
       ltur   = 0
    endif
    !
    ! check parameters with respect to the Coriolis force
    !
    if ( oned .and. coriolis ) then
       call msgerr (1, 'no Coriolis force in case of 1D computation')
       coriolis = .false.
    endif
    !
    if ( ( .not.pcor(1) /= 0. .or. .not.ylatc /= 0. ) .and. coriolis ) then
       call msgerr (0, 'Coriolis force will not be included')
       coriolis = .false.
    endif
    !
    if ( coriolis ) then
       if ( kspher > 0 .and. leds(23) == 2 ) then
          call msgerr (1, 'space varying Coriolis parameter has been read but will be ignored')
          call msgerr (1, 'since the grid is spherical, this parameter will be calculated from the local latitude coordinates')
          leds(23) = 0
       endif
    endif
    !
    ! check smoothing of boundary conditions/internal wave generation in case of restarting
    !
    if ( lrampf .and. restrt ) then
       if ( iwvgen == 0 ) then
          call msgerr (1, 'ramp function should not be applied to boundary conditions in case of hot-starting')
       else
          call msgerr (1, 'ramp function should not be applied to internal-generated waves in case of hot-starting')
       endif
       lrampf = .false.
    endif
    !
    ! check initial conditions
    !
    if ( instead .and. ( inituf .or. initvf ) ) then
       call msgerr (1, 'flow velocities have been initialized, no steady flow condition')
       instead = .false.
    else if ( instead .and. .not.initsf ) then
       call msgerr (1, 'water level should be initialized in case of steady flow condition')
       instead = .false.
    endif
    !
    if ( instead .and. irough == 4 .and. kmax > 1 ) then
       call msgerr (1, 'no steady flow condition in case of logarithmic wall-law')
       instead = .false.
    endif
    !
    ! check vertical grid schematisation
    !
    sumh    = 0.
    bperc   = .false.
    nonequi = .false.
    !
    do i = 1, kmax
       !
       if ( indlay(i) == 1 ) then
          sumh = sumh + hlay(i)
          if ( hlay(i)*real(kmax) /= 1. ) nonequi = .true.
          bperc = .true.
       endif
       !
    enddo
    !
    if ( abs(sumh - 1.0) > eps ) then
       !
       write (msgstr, '(a,i4,a)') 'the sum of layer thicknesses in percentage is ',nint(sumh * 100.0),' (should be 100)'
       call msgerr ( 3, trim(msgstr) )
       !
    endif
    !
    if ( .not.bperc ) call msgerr (3, 'at least one layer with variable thickness is required')
    !
    ! construct shadow sigma transformation, if appropriate
    !
    if (.not.allocated(hlaysh)) allocate(hlaysh(kmax))
    !
    if ( maxval(indlay) > 1 ) then
       fixlay = .true.
    else
       fixlay = .false.
    endif
    !
    if ( fixlay ) then
       !
       depfix = 0.
       !
       do i = 1, kmax
          !
          if ( indlay(i) /= 1 ) depfix = depfix + hlay(i)
          !
       enddo
       !
       ! threshold value equals twice the total thickness of fixed layers
       !
       epsfix = 2. * depfix
       !
       ! recalculate layer thicknesses for all layers
       !
       do i = 1, kmax
          !
          if ( indlay(i) == 1 ) then
             hlaysh(i) = hlay(i)/2.
          else
             hlaysh(i) = hlay(i)/epsfix
          endif
          !
       enddo
       !
    else
       !
       depfix = 0.
       epsfix = 0.
       hlaysh = hlay
       !
    endif
    !
    ! determine vertical pressure layer distribution for subgrid approach
    !
    if ( lsubg ) then
       !
       if ( fixlay ) call msgerr (3,'only layers with variable thickness must be applied in case of subgrid approach')
       !
       if ( ifloat /= 0 ) call msgerr (3,'subgrid approach not allowed in case of floating objects')
       !
       if ( optg == 5 ) call msgerr (3,'subgrid approach not supported in case of unstructured grids')
       !
       hlayp = 100./real(kpmax)
       !
       ! array npu contains number of velocity layers for each pressure layer
       ! note: pressure layers are equally distributed over the water depth as much as possible
       !
       if (.not.allocated(npu)) allocate(npu(kpmax))
       !
       sumh = 0.
       k    = 1
       l    = 0
       !
       do i = 1, kmax
          sumh = sumh + 100.*hlay(i)
          if ( .not. sumh < hlayp .or. i == kmax ) then
             if ( k > kpmax ) call msgerr (3, 'Internal error: number of elements of array npu exceeds a maximum')
             npu(k) = i  - l
             l      = l  + npu(k)
             k      = k  + 1
             sumh   = 0.
          endif
       enddo
       !
       ! array kup contains layer interfaces of coarse (pressure) grid as counted in fine (velocity) grid
       !
       if (.not.allocated(kup)) allocate(kup(0:kpmax))
       !
       kup(0) = 0
       do i = 1, kpmax
          kup(i) = sum(npu(1:i))
       enddo
       !
       ! pressure layers equally distributed?
       nonequi = .false.
       hlayp = sum(hlay(kup(0)+1:kup(1)))
       do i = 2, kpmax
          sumh = sum(hlay(kup(i-1)+1:kup(i)))
          if ( .not.EQREAL(sumh,hlayp) ) nonequi = .true.
          hlayp = sumh
       enddo
       !
       if ( nonequi .and. kpmax > 1 ) call msgerr (1,'pressure layers are not equally distributed over the water depth')
       !
       ! check total number of velocity layers in coarse (pressure) grid
       !
       if ( sum(npu) /= kmax ) then
          call msgerr (4,'inconsistency found in number of velocity layers due to subgrid approach')
          return
       endif
       !
       ! array kpv contains layer indices of fine (velocity) grid mapped to coarse (pressure) grid
       !
       if (.not.allocated(kpv)) allocate(kpv(1:kmax))
       !
       k = 1
       do i = 1, kpmax
          do j = 1, npu(i)
             kpv(k) = i
             k      = k + 1
          enddo
       enddo
       !
    endif
    !
    ! check non-hydrostatic schemes in case of depth-averaged mode or 1 pressure layer
    !
    if ( ihydro == 2 .and. kmax == 1 ) then
       call msgerr (1, 'only the box scheme for vertical pressure gradient should be applied in depth-averaged mode')
       ihydro = 1
    endif
    !
    if ( ihydro == 3 .and. kmax == 1 ) ihydro = 1
    !
    if ( qlay /= 0 .and. kpmax == 1 ) qlay = 0
    !
    if ( qlay /= 0 .and. ihydro /= 1 ) then
       call msgerr (1, 'reduced pressure Poisson equation method not supported in case of standard layout')
       qlay = 0
    endif
    !
    if ( qlay /= 0 .and. ifloat /= 0 ) then
       call msgerr (1, 'reduced pressure Poisson equation method not supported in case of floating objects')
       qlay = 0
    endif
    !
    if ( qlay /= 0 .and. qlay > kpmax - 1 ) then
       call msgerr (3, 'number of layers to reduce the pressure Poisson equation is too large')
    endif
    !
    qmax = kpmax - qlay
    !
    if ( LORB .and. .not.oned .and. qlay /= 0 .and. qmax == 1 ) call msgerr (3, 'stripwise grid partitioning is required in case of reduced method with one pressure layer')
    !
    ! check dispersion relation
    !
    if ( ( ihydro == 1 .and. kpmax < 5 .and. qlay == 0 ) .or. ihydro == 3 ) then
       numdisp = .true.
    else
       numdisp = .false.
    endif
    if ( nonequi .or. fixlay ) numdisp = .false.
    if ( .not. numdisp ) then
       ! use hyperbolic cosine profile in case of exact dispersion relation
       istok = 0
    else
       istok = -999
    endif
    if ( numdisp .and. kpmax > 3 ) numdisp = .false.
    !
    ! check default time integration with reference to computational mesh
    !
    if ( mtimei < 0 ) then
       if ( optg /= 5 ) then
          mtimei = 1
          if ( mimetic ) mtimei = 2
       else
          mtimei = 2
       endif
    endif
    !
    ! check Newton method with reference to computational mesh
    !
    if ( inewt < 0 ) then
       if ( optg /= 5 ) then
          inewt = 0
       else
          inewt = 1
       endif
    endif
    !
    if ( inewt /= 0 .and. mtimei /= 2 ) inewt = 0
    !
    ! check default projection method with reference to floating objects
    !
    if ( iproj < 0 ) then
       if ( ifloat == 0 ) then
          iproj = 1
       else
          iproj = 2
       endif
    endif
    !
    ! check default preconditioner with reference to floating objects
    !
    if ( icond < 0 ) then
       if ( ifloat == 0 ) then
          icond = 2
          if ( optg == 5 ) icond = 3
       else
          icond = 3
       endif
    endif
    !
    ! check inclusion of non-hydrostatic pressure part
    !
    if ( lpproj .and. mtimei /= 2 ) then
       call msgerr (1, 'inclusion of non-hydrostatic pressure part in continuity equation is superfluous')
       lpproj = .false.
    endif
    !
    if ( .not.lpproj ) pnums(59) = 1.
    !
    ! check sponge layer widths
    !
    if ( oned .and. ( spwidb > 0. .or. spwidt > 0. ) ) then
       call msgerr (1, '1D simulation: width of sponge layer in y-direction set to zero !')
       spwidb = 0.
       spwidt = 0.
    endif
    !
    if ( lreptx .and. ( spwidl > 0. .or. spwidr > 0. ) ) then
       call msgerr (1, 'no sponge layer in x-direction in case of repeating grid')
       spwidl = 0.
       spwidr = 0.
    else if ( lrepty .and. ( spwidb > 0. .or. spwidt > 0. ) ) then
       call msgerr (1, 'no sponge layer in y-direction in case of repeating grid')
       spwidb = 0.
       spwidt = 0.
    endif
    !
    if ( optg == 1 ) then
       !
       if ( kspher == 0 ) then
          !
          xcl = xclen
          ycl = yclen
          !
       else
          !
          xcl = lendeg * cos(degrad*yoffs) * xclen
          ycl = lendeg * yclen
          !
       endif
       !
       if ( spwidl + spwidr > 0.9*xcl ) then
          call msgerr (2, 'total width of sponge layers more than 90% of length of the computational domain in x-direction')
       endif
       !
       if ( spwidb + spwidt > 0.9*ycl ) then
          call msgerr (2, 'total width of sponge layers more than 90% of length of the computational domain in y-direction')
       endif
       !
    endif
    !
    if ( lwavoutp .and. ( spwidl > 0. .or. spwidr > 0. .or. spwidb > 0. .or. spwidt > 0. ) ) then
       call msgerr (1, 'computation of wave-induced setup will not be correct due to the sponge layers')
    endif
    !
    ! check internal wave generation
    !
    if ( iwvgen /= 0 ) then
       if ( nbgrpt > 0 ) call msgerr (3, 'concurrent specification of boundary conditions and internal wave generation is not allowed')
    endif
    !
    ! check mimetic discretizations
    !
    if ( mimetic .and. mtimei /= 2 ) call msgerr (3, 'semi-implicit time integration is required in case of mimetic discretizations')
    !
    ! check horizontal eddy viscosity coefficient
    !
    if ( oned .and. ihvisc == 2 ) then
       call msgerr (1, 'no Smagorinsky model for computing eddy viscosity in case of 1D computation')
       ihvisc = 0
    else if ( iturb > 1 .and. ihvisc /= 0 ) then
       call msgerr (1, 'no horizontal mixing in case of 3D turbulence')
       ihvisc = 0
    endif
    !
    ! check vertical eddy viscosity coefficient
    !
    if ( iturb /= 0 .and. kmax == 1 ) then
       call msgerr (1, 'no vertical mixing in depth-averaged mode')
       iturb = 0
       ltur  = 0
    endif
    !
    ! check turbulence model in conjunction with roughness method
    !
    if ( iturb /= 0 .and. irough == 0 ) then
       call msgerr (3, 'roughness method must be given for the applied 3D turbulence model')
    else if ( iturb /= 0 .and. irough == 11 ) then
       call msgerr (3, 'linear bottom friction not allowed in case of 3D turbulence model')
    else if ( iturb > 1 .and. irough /= 4 ) then
       call msgerr (3, 'logarithmic wall-law must be applied in case of 3D turbulence')
    else if ( irough == 4 .and. kmax > 1 .and. iturb == 0 ) then
       call msgerr (3, 'k-eps turbulence model must be employed in case of logarithmic wall-law')
    endif
    !
    ! check 3D turbulence in conjunction with non-hydrostatic mode
    !
    if ( iturb > 1 .and. ihydro /= 2 ) then
       call msgerr (3, 'non-hydrostatic mode (option STANDARD) must be switch on in case of 3D turbulence')
    endif
    !
    ! check availability of 3D turbulence model in case of unstructured grids
    !
    if ( iturb > 1 .and. optg == 5 ) then
       call msgerr (3,'3D turbulence model not supported in case of unstructured grids')
    endif
    !
    ! check outputting mean turbulence quantities
    !
    if ( iturb == 0 ) lturoutp = .false.
    !
    ! check necessity of advection terms of the w-momentum equation
    !
    if ( ihydro == 2 .and. ( iturb /= 0 .or. itrans == 1 ) ) then
       horwinc = .true.
       verwinc = .true.
    endif
    !
    ! check vegetation
    !
    if ( iveg /= 0 .and. lmax == 0 ) call msgerr (3,'input for vegetation is required')
    !
    if ( iveg /= 0 ) then
       if ( cvm < 0. ) then
          call msgerr (1, 'incorrect added mass coefficient')
          cvm = 0.
       endif
       if ( iveg == 2 .and. iporos /= 0 ) then
          call msgerr (1, 'porous model should not be applied in case of vegetation porosity')
          iporos = 0
       endif
    endif
    !
    ! check transport of constituent
    !
    if ( itrans /= 0 .and. ltrans == 0 ) call msgerr (3, 'at least one input grid for constituent must be given in case of transport')
    !
    ! check return time for unsteady salt intrusion
    !
    if ( lsal > 0 ) then
       if ( tcret < 0. ) then
          call msgerr (1, 'incorrect return time for salt intrusion')
          tcret = 0.
       endif
    endif
    !
    ! check anti-creep
    !
    if ( itrans == 0 .or. kmax == 1 ) icreep = 0
    !
    if ( icreep == 2 .and. ( lreptx .or. lrepty ) ) call msgerr (3,'repeating grid not allowed in case of Stelling and Van Kester method')
    !
    ! check if mass exchange is applicable for sediment transport
    !
    if ( lsed > 0 ) then
       if ( kmax == 1 ) then
          call msgerr (1, 'no mass exchange of sediment in depth-averaged mode')
          psed( 1) = 0.
          psed( 2) = 0.
          psed(10) = 0.
          psed(11) = 0.
          psed(12) = 0.
       endif
       if ( .not.lmixt ) then
          call msgerr (1, 'no mass exchange of sediment for tracer')
          psed( 1) = 0.
          psed( 2) = 0.
          psed(10) = 0.
          psed(11) = 0.
          psed(12) = 0.
       endif
    endif
    !
    ! check fall velocity and sediment size for suspended sediment
    !
    if ( lsed > 0 ) then
       if ( psed(1) < 0. ) then
          call msgerr (1, 'incorrect sediment settling velocity')
          psed(1) = 0.
       else
          ! fall velocity is given by user in millimeter per second
          psed(1) = psed(1) * 1.e-3
       endif
       if ( psed(2) < 0. ) then
          call msgerr (1, 'incorrect median sediment diameter')
          psed(2) = 0.
       else
          ! grain size is given by user in micrometer
          psed(2) = psed(2) * 1.e-6
       endif
    endif
    !
    ! check turbulence model and bottom friction related to sediment transport in 3D buoyancy flow
    !
    if ( lsed > 0 .and. (psed(2) /= 0. .or. psed(10) /= 0. .or. lsal > 0 .or. ltemp > 0) .and. kmax > 1 .and. lmixt ) then
       if ( iturb /= 1 ) then
          call msgerr (3, 'vertical mixing (k-epsilon model) must be employed in case of sediment transport')
       endif
       if ( irough /= 4 ) then
          call msgerr (3, 'logarithmic wall-law must be applied in case of sediment transport')
       endif
       if ( instead ) then
          call msgerr (1, 'no steady flow condition in case of sediment transport')
          instead = .false.
       endif
       ! calculate the Nikuradse roughness height in case of sand transport, if appropriate
       if ( .not.varfr .and. .not. pbot(2) /= 0. .and. psed(2) /= 0. ) then
          pbot(2) = psed(4) * psed(2)
       endif
    endif
    !
    ! check VARANS equations
    !
    if ( iporos == 1 .and. kmax > 1 ) then
       horwinc = .true.
       verwinc = .true.
       if ( ihydro /= 1 .and. ihydro /= 2 ) then
          call msgerr (3, 'non-hydrostatic mode must be switch on in case of VARANS equations')
       endif
    endif
    !
    ! check ambient current
    !
    if ( iamb /= 0 ) then
       if ( mtimei == 2 ) call msgerr (3, 'ambient current not supported in case of semi-implicit time integration')
       if (   optg == 5 ) call msgerr (3, 'ambient current not supported in case of unstructured grids')
       if ( pnums(71) < 0.5 .or. pnums(71) > 1. ) then
          call msgerr (1, 'incorrect theta value for time integration of global continuity equation due to ambient current')
          pnums(71) = 1.
       endif
       horwinc = .true.
       verwinc = .true.
    endif
    !
    ! check outputting mean transported constituents
    !
    if ( itrans == 0 .or. ltrans == 0 ) ltraoutp = .false.
    !
    ! check baroclinic forcing
    !
    if ( lsal > 0 .or. ltemp > 0 .or. (lsed > 0 .and. lmixt) ) idens = 1
    !
    ! if no baroclinic forcing, recompute water density based on temperature and salinity according to the Eckart's formula (1958)
    ! Note: the computed density is the actual density and not a reference density!
    !
    if ( idens == 0 .and. .not. rhow /= -999. ) then
       !
       pf0    =  5890d0  + (  38d0  - 375d-3*dble(tempw))*dble(tempw) +                        3d0*dble(salw)
       lambda = 17795d-1 + (1125d-2 - 745d-4*dble(tempw))*dble(tempw) - (38d-1 + 1d-2*dble(tempw))*dble(salw)
       ! density based on Eckart is g/ml, hence multiply by 1000 to get g/l
       rhow   = 1000.*pf0/(lambda+0.698*pf0)
       !
    else if ( .not. rhow /= -999. ) then
       !
       rhow = 1000.
       !
    endif
    !
    ! check background viscosity
    !
    if ( bvisc < 0. ) then
       bvisc = 0.
       ! note: including baroclinic forcing may amplify horizontal divergence checkerboard modes on unstructured mesh
       !       therefore, add background viscosity to filter these spurious modes
       if ( optg == 5 .and. idens /= 0 ) bvisc = 5.e-5
    endif
    !
    ! kinematic viscosity of water
    !
    kinvis = dynvis / rhow
    !
    ! compute some factors meant for noncohesive sediment sand
    !
    if ( lsed > 0 ) then
       if ( psed(2) > 0. ) then
          psed(8) = (-1.+rhos/rhow) * grav * psed(2)
          psed(9) = psed(6) * psed(8)**0.6 * psed(2)**0.2 / ( kinvis**0.2 )
          psed(5) = psed(5) * psed(8)
       else
          psed(8) = 0.
          psed(9) = 0.
       endif
       ! compute fall velocity based on sediment size, if appropriate
       if ( .not. psed(1) /= 0. .and. psed(2) > 0. ) then
          fac = 36.*kinvis*kinvis / ( psed(8) * psed(2)*psed(2) )
          psed(1) = sqrt(psed(8)) * ( sqrt(fac+2./3.) - sqrt(fac) )
       endif
    endif
    !
    ! coefficients for transformation from user coordinates to computational coordinates
    !
    cospc = cos(alpc)
    sinpc = sin(alpc)
    xcp   = -xpc * cospc - ypc * sinpc
    ycp   =  xpc * sinpc - ypc * cospc
    alcp  = -alpc
    !
    ! check location of output areas
    !
    call SPRCON ( xcgrid, ycgrid, kgrpnt )
    !
    ! check duration for wave, mean current, mean transported constituent and mean turbulence output
    !
    if ( lwavoutp .and. twavoutp < 0. ) call msgerr (2, 'parameter [dur] for setup or wave height is not specified')
    if ( lcuroutp .and. tcuroutp < 0. ) call msgerr (2, 'parameter [dur] for time-averaged current is not specified')
    if ( ltraoutp .and. ttraoutp < 0. ) call msgerr (2, 'parameter [dur] for time-averaged constituents is not specified')
    if ( lturoutp .and. tturoutp < 0. ) call msgerr (2, 'parameter [dur] for time-averaged turbulence quantities is not specified')
    !
    ! compute geometric quantities based on given grid (rotated 2D rectilinear or curvilinear grid)
    !
    if ( .not.oned .and. optg /= 5 ) then
       !
       istat = 0
       if ( .not.allocated(guu) ) allocate (guu(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(guv  ) ) allocate (guv  (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(gvu  ) ) allocate (gvu  (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(gvv  ) ) allocate (gvv  (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(gsqs ) ) allocate (gsqs (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(gsqsu) ) allocate (gsqsu(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(gsqsv) ) allocate (gsqsv(mcgrd), stat = istat)
       !
       if ( istat /= 0 ) then
          write (msgstr, '(a,i6)') 'allocation problem: geometric quantities and return code is ',istat
          call msgerr ( 4, trim(msgstr) )
          return
       endif
       !
       guu   = 1.
       guv   = 1.
       gvu   = 1.
       gvv   = 1.
       gsqs  = 1.
       gsqsu = 1.
       gsqsv = 1.
       !
       call SwashGeometrics
       !
       call SWEXCHG (   guu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG (   guv, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG (   gvu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG (   gvv, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG (  gsqs, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( gsqsu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( gsqsv, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize geometric quantities at appropriate boundaries in case of repeating grid
       !
       call periodic (   guu, kgrpnt, 1, 1 )
       call periodic (   guv, kgrpnt, 1, 1 )
       call periodic (   gvu, kgrpnt, 1, 1 )
       call periodic (   gvv, kgrpnt, 1, 1 )
       call periodic (  gsqs, kgrpnt, 1, 1 )
       call periodic ( gsqsu, kgrpnt, 1, 1 )
       call periodic ( gsqsv, kgrpnt, 1, 1 )
       !
    else
       !
       if ( .not.allocated(guu  ) ) allocate (guu  (0))
       if ( .not.allocated(guv  ) ) allocate (guv  (0))
       if ( .not.allocated(gvu  ) ) allocate (gvu  (0))
       if ( .not.allocated(gvv  ) ) allocate (gvv  (0))
       if ( .not.allocated(gsqs ) ) allocate (gsqs (0))
       if ( .not.allocated(gsqsu) ) allocate (gsqsu(0))
       if ( .not.allocated(gsqsv) ) allocate (gsqsv(0))
       !
    endif
    !
    ! wind direction w.r.t. computational grid
    !
    wdic = pi2 * (wdip/pi2 - nint(wdip/pi2))
    !
    ! check relativeness of wind velocity
    !
    if ( relwav .and. relwnd ) then
       call msgerr (1, 'wind velocity is taken relative to wave celerity')
       relwnd  = .false.
    endif
    !
    ! check rescaling factor for wind velocity
    !
    if ( (.not.relwnd .or. relwav) .and. pwnd(4) /= 1. ) then
       call msgerr (1, 'rescaling factor for wind velocity not applied')
       pwnd(4) = 1.
    else if ( .not. ( pwnd(4) > 0. .and. pwnd(4) <= 1. ) ) then
       call msgerr (1, 'incorrect rescaling factor for wind velocity')
       pwnd(4) = 1.
    endif
    !
    ! check ratio of forced crest height to maximum surface elevation
    !
    if ( pwnd(10) < 0. .or. pwnd(10) > 1. ) then
       call msgerr (1, 'incorrect ratio of forced crest height to maximum surface elevation')
       pwnd(10) = 0.4
    endif
    !
    ! check control wave breaking; no control in hydrostatic mode
    !
    if ( ihydro == 0 ) isurf = 0
    !
    ! check associated wave breaking parameters
    !
    if ( .not. psurf(1) > 0. ) then
       call msgerr (1, 'incorrect threshold value for initiation of wave breaking')
       if ( isurf == 1 ) then
          psurf(1) = 0.6
       else if ( isurf == 2 ) then
          psurf(1) = 0.4
       endif
    endif
    if ( isurf == 1 .and. .not. psurf(2) > 0. .and. psurf(2) /= -1. ) then
       call msgerr (1, 'incorrect threshold value for re-initiation of wave breaking in post-breaking area')
       psurf(2) = 0.3
    endif
    !
    ! link breaking parameters to horizontal momentum advection scheme
    !
    if ( isurf == 1 ) then
       !
       ! choose beta = 0.15 when BDF scheme
       !
       if ( pnums(6) == 3. .and. pnums(7) == -1. .and. psurf(2) == -1. ) then
          !
          psurf(2) = 0.15
          !
       else if ( psurf(2) == -1. ) then
          !
          psurf(2) = 0.3
          !
       endif
       !
       ! choose central differences when beta = 0.3
       !
       if ( .not. psurf(2) /= 0.3 .and. pnums(6) == -1. ) then
          !
          pnums(6) = 3.
          pnums(7) = 1.
          !
       else if ( pnums(6) == -1. ) then
          !
          pnums(6) =  3.
          pnums(7) = -1.
          !
       endif
       !
    else
       !
       if ( pnums(6) == -1. ) then
          !
          pnums(6) =  3.
          pnums(7) = -1.
          !
       endif
       !
    endif
    !
    ! check parameters of floating objects
    !
    if ( ifloat /= 0 ) then
       !
       horwinc = .false.
       !
       if ( mtimei /= 2 ) call msgerr (3, 'semi-implicit time integration is required in case of floating object or pressurized flow')
       !
       if ( pship(1) < 0. .or. pship(1) > 1. ) then
          call msgerr (1, 'incorrect compressibility factor for flow beneath fixed floating object')
          pship(1) = 0.
       endif
       !
       if ( pship(2) < 0.5 .or. pship(2) > 1. ) then
          call msgerr (1, 'incorrect theta value for time integration beneath floating object')
          pship(2) = 1.
       endif
       !
       if ( .not. pship(2) /= 0.5 ) pship(2) = 0.51
       !
       ! check motion rigid bodies
       !
       if ( ifloat /= 2 .and. mbod > 1 ) then
          call msgerr (3, 'specification of labels for rigid bodies as input grid is required')
       else if ( mbod == 1 ) then
          ifloat = 2
       else if ( mbod == 0 ) then
          ifloat = min( 1, ifloat )
       endif
       !
       ! in case of moving object, set compressibility factor to zero
       !
       if ( ifldyn(18) /= 0 .or. ifloat == 2 ) pship(1) = 0.
       !
       ! for computing hydrodynamic loads acting on (moving) objects accurately, only user-defined draft levels and an upwind depth interpolation is allowed
       !
       if ( ifloat /= 0 ) dpsopt = 4
       if ( ifloat == 2 ) then
          if ( pnums(11) == 3. .and. pnums(12) == 1. ) then
             pnums(11) = 6.
             pnums(12) = 0.
          endif
       endif
       !
    endif
    !
    ! close input files containing stationary input fields
    !
    do i = 1, numgrd
       if ( ifldyn(i) == 0 .and. iflnds(i) /= 0 ) then
          close (iflnds(i))
          iflnds(i) = 0
       endif
    enddo
    !
    ! check parameters for internal wave generation
    !
    if ( iwvgen /= 0 ) then
       !
       if ( .not. piwg(2) > 0. ) then
          call msgerr (2, 'incorrect width of source area')
       endif
       if ( .not. piwg(3) > 0. ) then
          call msgerr (2, 'incorrect still water depth at source area')
       endif
       if ( .not. piwg(4) > 0. ) then
          call msgerr (1, 'incorrect delta shape factor for internal wave generation')
          piwg(4) = 0.5
       endif
       !
    endif
    !
    ! check minimal depth for drying and flooding
    !
    if ( .not. epsdry > 0. ) then
       call msgerr (1, 'minimal depth for drying/wetting should not be negative or zero')
       epsdry = 5.e-5
    else if ( epsdry > 0.01 ) then
       call msgerr (1, 'minimal depth for drying/wetting should not be larger than 1 cm')
       epsdry = 0.01
    endif
    !
    ! check inundation depth
    !
    if ( hrunp /= -999. ) then
       if ( .not. hrunp > 0. ) then
          call msgerr (1, 'depth for horizontal runup should be larger than zero')
          hrunp = epsdry
       else if ( hrunp < epsdry ) then
          hrunp = epsdry
       endif
    else
       hrunp = epsdry
    endif
    !
    ! check threshold depth for runup level
    !
    if ( delrp /= -999. ) then
       if ( .not. delrp > 0. ) then
          call msgerr (1, 'threshold depth for wave runup should be larger than zero')
          delrp = epsdry
       else if ( delrp < epsdry ) then
          delrp = epsdry
       endif
    else
       delrp = epsdry
    endif
    !
    if ( inewt == 0 ) epshu = epsdry
    !
    ! check porosity parameters
    !
    if ( iporos == 1 ) then
       !
       if ( .not.vargs .and. .not. ppor(1) > 0. ) then
          call msgerr (2, 'incorrect grain size')
       endif
       if ( varsh ) ppor(2) = 99999.
       if ( .not.varsh .and. .not. ppor(2) > 0. ) then
          call msgerr (2, 'incorrect structure height')
       endif
       if ( ppor(3) < 0. ) then
          call msgerr (1, 'incorrect constant for laminar friction factor')
          ppor(3) = 200.
       endif
       if ( ppor(4) < 0. ) then
          call msgerr (1, 'incorrect constant for turbulent friction factor')
          ppor(4) = 1.1
       endif
       !
    endif
    !
    ! check CFL numbers
    !
    if ( .not. ( pnums(2) > 0. .and. pnums(2) < 1. ) ) then
       call msgerr (1, 'incorrect lowest CFL number')
       pnums(2) = 0.4
    endif
    if ( .not. ( pnums(3) > 0. .and. pnums(3) < 1. ) ) then
       call msgerr (1, 'incorrect highest CFL number')
       pnums(3) = 0.8
    endif
    if ( pnums(3) < pnums(2) ) then
       call msgerr (1, 'highest CFL number is always larger than lowest CFL number')
       pnums(2) = 0.4
       pnums(3) = 0.8
    endif
    !
    ! initialize maximum CFL number
    !
    cflmax = 0.5 * (pnums(2) + pnums(3))
    !
    ! check implicitness factors
    !
    if ( pnums(1) < 0.5 .or. pnums(1) > 1. ) then
       call msgerr (1, 'incorrect theta value for time integration of continuity equation')
       pnums(1) = 0.5
    endif
    if ( pnums(4) < 0.5 .or. pnums(4) > 1. ) then
       call msgerr (1, 'incorrect theta value for time integration of water level gradient')
       pnums(4) = 0.5
    endif
    if ( pnums(5) > 0. .and. pnums(5) < 0.5 .or. pnums(5) > 1. ) then
       call msgerr (1, 'incorrect theta value for time integration of non-hydrostatic pressure gradient')
       pnums(5) = -1.
    endif
    if ( pnums(5) < 0. ) then
       !if ( ihydro == 2 ) then
       !   pnums(5) = 0.50
       !else if ( kmax == 1 .or. ihydro == 3 ) then
       !   pnums(5) = 0.51
       !else
       !   pnums(5) = 0.60
       !endif
       pnums(5) = 1.
    endif
    if ( pnums(31) < 0.5 .or. pnums(31) > 1. ) then
       call msgerr (1, 'incorrect theta value for time integration of vertical terms in u-momentum equation')
       pnums(31) = 0.5
    endif
    if ( pnums(32) < 0.5 .or. pnums(32) > 1. ) then
       call msgerr (1, 'incorrect theta value for time integration of vertical terms in w-momentum equation')
       pnums(32) = 0.5
    endif
    if ( pnums(33) < 0.5 .or. pnums(33) > 1. ) then
       call msgerr (1, 'incorrect theta value for time integration of vertical terms in transport equation')
       pnums(33) = 0.5
    endif
    !
    if ( (irough == 0 .or. .not. pbot(1) /= 0.) .and. irough /= 4 .and. .not.mimetic .and. iwvgen == 0 ) then
       if ( .not. pnums(1) /= 0.5 ) pnums(1) = 0.51
       if ( .not. pnums(4) /= 0.5 ) pnums(4) = 0.51
    endif
    !
    ! check termination criterion for linear solvers
    !
    if ( pnums(21) < 0. ) then
       if ( inewt == 0 ) then
          pnums(21) = 0.01
       else
          pnums(21) = 0.001
       endif
    endif
    if ( pnums(22) < 0. ) then
       call msgerr (1, 'incorrect termination criterion for the SIP or BiCGSTAB solver')
       pnums(22) = 0.01
    endif
    if ( pnums(23) < 0. ) then
       call msgerr (1, 'incorrect termination criterion for the BiCGSTAB method')
       pnums(23) = 0.
    endif
    !
    ! check and set weight factors / relaxation parameter for linear solvers
    !
    if ( pnums(26) > 1. ) then
       call msgerr (1, 'incorrect weight factor for the CG method')
       pnums(26) = -1.
    endif
    if ( pnums(26) < 0. ) then
       if ( optg /= 5 ) then
          pnums(26) = 0.9
       else
          pnums(26) = 0.
       endif
    endif
    if ( pnums(27) > 1. ) then
       call msgerr (1, 'incorrect weight factor for the SIP or BiCGSTAB solver')
       pnums(27) = -1.
    endif
    if ( pnums(27) < 0. ) then
       if ( optg == 5 ) then
          pnums(27) = 0.
       else if ( kmax == 1 .or. ihydro == 3 ) then
          pnums(27) = 0.91
       else if ( qmax == 1 ) then
          pnums(27) = 0.93
       else
          if ( ihydro == 1 ) then
             if ( icond == 1 ) then
                pnums(27) = 0.
             else if ( icond == 2 ) then
                pnums(27) = 0.99
             else if ( icond == 3 ) then
                pnums(27) = 0.55
             endif
          else if ( ihydro == 2 ) then
             if ( icond == 3 .and. ifloat /= 0 ) then
                pnums(27) = 0.
             else
                pnums(27) = 0.9
             endif
          endif
       endif
    endif
    !
    ! check water depth correction
    !
    if ( corrdep ) then
       if ( nint(pnums(11)) == 1 ) corrdep = .false.
       if ( ( nint(pnums(11)) == 3 .and. pnums(12) == 1. ) .or. mimetic ) then
          corrdep = .false.
          depcds  = .true.
       endif
    endif
    !
    ! check generalized-alpha method
    !
    if ( mfoti == 1 ) then
       if ( pship(3) < 0. .or. pship(3) > 0.5 ) then
          call msgerr (1, 'incorrect beta for the Newmark scheme')
          pship(3) = 0.25
       endif
       if ( pship(4) < 0. .or. pship(4) > 1. ) then
          call msgerr (1, 'incorrect gamma for the Newmark scheme')
          pship(4) = 0.5
       endif
    else if ( mfoti > 1 ) then
       rhoinf = pship(5)
       if ( rhoinf < 0. .or. rhoinf > 1. ) then
          call msgerr (1, 'incorrect spectral radius for the generalized-alpha method')
          rhoinf = 1.
       endif
    endif
    !
    if ( mfoti == 1 ) then
       !
       ! Newmark scheme
       !
       alfam = 0.
       alfaf = 0.
       !
    else if ( mfoti == 2 ) then
       !
       ! Chung-Hulbert scheme
       !
       alfam = (2.*rhoinf-1.) / (1.+rhoinf)
       alfaf = rhoinf / (1.+rhoinf)
       !
    else if ( mfoti == 3 ) then
       !
       ! Hilber-Hughes-Taylor scheme
       !
       alfam = 0.
       alfaf = (1.-rhoinf) / (1.+rhoinf)
       !
    else if ( mfoti == 4 ) then
       !
       ! Wood-Bossak-Zienkiewicz scheme
       !
       alfam = (rhoinf-1.) / (1.+rhoinf)
       alfaf = 0.
       !
    endif
    !
    if ( mfoti > 1 ) then
       !
       fac      = 1. - alfam + alfaf
       pship(3) = 0.25 * fac * fac
       pship(4) = fac - 0.5
       !
    endif
    !
    ! check FSI coupling
    !
    if ( ifloat == 2 ) then
       if ( pship(10) < 0. .or. pship(10) > 1. ) then
          call msgerr (1, 'incorrect under-relaxation value for FSI coupling')
          pship(10) = 0.3
       endif
    else
       ! no iteration process required!
       pship(8) = 0.
       pship(9) = 1.
    endif
    !
    ! check implicitness factor for kinematic boundary conditions
    !
    if ( ifloat == 2 ) then
       if ( pship(11) < 0. .or. pship(11) > 1. ) then
          call msgerr (1, 'incorrect theta value for update kinematic boundary conditions')
          pship(11) = 0.5
       endif
    endif
    !
    ! check whether linearized shallow water equations are employed
    !
    if ( linswe ) then
       !
       jhk      = 3
       pnums(6) = 1.
       horwinc  = .false.
       verwinc  = .false.
       !
    endif
    !
    ! allocate arrays for storing flow variables at previous time level in case of nonstationary run
    !
    if ( nstatc == 1 ) then
       !
       istat = 0
       if ( .not.allocated(s0) ) then
          if ( optg /= 5 ) then
             allocate (s0(mcgrd), stat = istat)
          else
             allocate (s0(ncells), stat = istat)
          endif
       else if ( size(s0) == 0 ) then
          deallocate (s0)
          if ( optg /= 5 ) then
             allocate (s0(mcgrd), stat = istat)
          else
             allocate (s0(ncells), stat = istat)
          endif
       endif
       !
       if ( istat == 0 ) then
          if ( .not.allocated(u0) ) then
             if ( optg /= 5 ) then
                allocate (u0(mcgrd,kmax), stat = istat)
             else
                allocate (u0(nfaces,kmax), stat = istat)
             endif
          else if ( size(u0) == 0 ) then
             deallocate (u0)
             if ( optg /= 5 ) then
                allocate (u0(mcgrd,kmax), stat = istat)
             else
                allocate (u0(nfaces,kmax), stat = istat)
             endif
          endif
       endif
       !
       if ( optg /= 5 ) then
          if ( istat == 0 ) then
             if ( .not.allocated(v0) ) then
                if ( .not.oned ) then
                   allocate (v0(mcgrd,kmax), stat = istat)
                else
                   allocate (v0(0,0))
                endif
             else if ( size(v0) == 0 ) then
                deallocate (v0)
                if ( .not.oned ) then
                   allocate (v0(mcgrd,kmax), stat = istat)
                else
                   allocate (v0(0,0))
                endif
             endif
          endif
       endif
       !
       if ( istat /= 0 ) then
          write (msgstr, '(a,i6)') 'allocation problem: array s0, u0 or v0 and return code is ',istat
          call msgerr ( 4, trim(msgstr) )
          return
       endif
       !
       s0 = 0.
       u0 = 0.
       if ( .not.oned .and. optg /= 5 ) v0 = 0.
    else
       if ( .not.allocated(s0) ) allocate (s0(0))
       if ( .not.allocated(u0) ) allocate (u0(0,0))
       if ( .not.allocated(v0) ) allocate (v0(0,0))
    endif
    !
    if ( alobnd .and. allocated(bndval) ) deallocate(bndval)
    if ( .not.allocated(bndval) ) allocate(bndval(nbval,2))
    if ( alobnd ) bndval = 0.
    !
    ! compute incident spectrum and store wave components, if appropriate
    !
    nfreq = -999
    if ( alobnd ) then
       curbfs => fbfs
       do
          if ( curbfs%nfreq == -1 ) then
             call SwashBCspectrum ( curbfs, curbfs%spparm )
             if (STPNOW()) return
          endif
          nfreq = max( nfreq, curbfs%nfreq )
          if ( .not.associated(curbfs%nextbfs) ) exit
          curbfs => curbfs%nextbfs
       enddo
    endif
    !
    ! store boundary grid points and their data in own subdomain
    !
    if ( alobnd .and. allocated(bgridp) ) deallocate( bgridp )
    if ( .not.allocated(bgridp) ) allocate( bgridp(9*nbgrpt) )
    if ( alobnd ) then
       bgridp = 0
       if ( optg /= 5 ) then
          if ( NBGGL == 0 ) NBGGL = nbgrpt
          nbgrpt = 0
          do ix = MXF, MXL
             do iy = MYF, MYL
                indx = KGRPGL(ix,iy)
                curbgp => fbgp
                do i = 1, NBGGL
                   j = curbgp%bgp(1)
                   if ( j == indx ) then
                      nbgrpt = nbgrpt + 1
                      bgridp(9*nbgrpt-8) = int(kgrpnt(ix-MXF+1,iy-MYF+1),8)
                      bgridp(9*nbgrpt-7) = curbgp%bgp(2)
                      bgridp(9*nbgrpt-6) = curbgp%bgp(3)
                      bgridp(9*nbgrpt-5) = curbgp%bgp(4)
                      bgridp(9*nbgrpt-4) = curbgp%bgp(5)
                      bgridp(9*nbgrpt-3) = curbgp%bgp(6)
                      bgridp(9*nbgrpt-2) = curbgp%bgp(7)
                      bgridp(9*nbgrpt-1) = curbgp%bgp(8)
                      bgridp(9*nbgrpt  ) = curbgp%bgp(9)
                   endif
                   if (.not.associated(curbgp%nextbgp)) exit
                   curbgp => curbgp%nextbgp
                enddo
             enddo
          enddo
       else
          curbgp => fbgp
          do i = 1, nbgrpt
             bgridp(9*i-8) = curbgp%bgp(1)
             bgridp(9*i-7) = curbgp%bgp(2)
             bgridp(9*i-6) = curbgp%bgp(3)
             bgridp(9*i-5) = curbgp%bgp(4)
             bgridp(9*i-4) = curbgp%bgp(5)
             bgridp(9*i-3) = curbgp%bgp(6)
             bgridp(9*i-2) = curbgp%bgp(7)
             bgridp(9*i-1) = curbgp%bgp(8)
             bgridp(9*i  ) = curbgp%bgp(9)
             if (.not.associated(curbgp%nextbgp)) exit
             curbgp => curbgp%nextbgp
          enddo
       endif
    endif
    !
    ! check whether Stokes waves are imposed at open boundaries
    !
    do i = 1, nbgrpt
       itmp = int(bgridp(9*i-1),4)
       if ( itmp == -1 ) then
          if ( istok  < 0 ) istok = 0
       elseif ( itmp == -11 ) then
          if ( istok  < 0 ) istok = 1
          if ( istok == 0 ) bgridp(9*i-1) = -1
       elseif ( itmp == -21 ) then
          if ( istok  < 0 ) istok = 2
          if ( istok == 0 ) bgridp(9*i-1) = -1
       elseif ( itmp == -31 ) then
          if ( istok  < 0 ) istok = 3
          if ( istok == 0 ) bgridp(9*i-1) = -1
       endif
    enddo
    !
    if ( istok == 3 .and. lsubg ) call msgerr (3,'second order transfer functions not supported in case of subgrid approach')
    !
    ! allocate arrays for incident short and bound wave components, if appropriate
    !
    if ( nfreq > 0 ) then
       !
       if ( nbgrpt > 0 ) then
          !
          if ( alobnd .and. allocated(kwave) ) deallocate(kwave,comp1,comp2)
          !
          if ( .not.allocated(kwave) ) allocate(kwave(nbgrpt,nfreq))
          if ( .not.allocated(comp1) ) allocate(comp1(nbgrpt,nfreq))
          if ( .not.allocated(comp2) ) allocate(comp2(nbgrpt,nfreq))
          !
          kwave = 0.
          comp1 = 0.
          comp2 = 0.
          !
          if ( alobnd .and. allocated(zetab) ) deallocate(zetab,fluxbu,fluxbv)
          !
          if ( .not.allocated(zetab ) ) allocate(zetab (nbgrpt,nfreq))
          if ( .not.allocated(fluxbu) ) allocate(fluxbu(nbgrpt,nfreq))
          if ( .not.allocated(fluxbv) ) allocate(fluxbv(nbgrpt,nfreq))
          !
          zetab  = (0.,0.)
          fluxbu = (0.,0.)
          fluxbv = (0.,0.)
          !
          if ( istok /= 0 ) then
             !
             if ( alobnd .and. allocated(stku1c) ) deallocate(stku1c,stku1s,stkz1c,stkz1s)
             !
             if ( .not.allocated(stku1c) ) allocate(stku1c(nbgrpt,nfreq,kmax))
             if ( .not.allocated(stku1s) ) allocate(stku1s(nbgrpt,nfreq,kmax))
             !
             if ( .not.allocated(stkz1c) ) allocate(stkz1c(nbgrpt,nfreq))
             if ( .not.allocated(stkz1s) ) allocate(stkz1s(nbgrpt,nfreq))
             !
             stku1c = 0.
             stku1s = 0.
             !
             stkz1c = 0.
             stkz1s = 0.
             !
             if ( istok > 1 ) then
                !
                if ( alobnd .and. allocated(stku2c) ) deallocate(stku2c,stku2s,stkz2c,stkz2s)
                !
                if ( .not.allocated(stku2c) ) allocate(stku2c(nbgrpt,nfreq,kmax))
                if ( .not.allocated(stku2s) ) allocate(stku2s(nbgrpt,nfreq,kmax))
                !
                if ( .not.allocated(stkz2c) ) allocate(stkz2c(nbgrpt,nfreq))
                if ( .not.allocated(stkz2s) ) allocate(stkz2s(nbgrpt,nfreq))
                !
                stku2c = 0.
                stku2s = 0.
                !
                stkz2c = 0.
                stkz2s = 0.
                !
                if ( istok == 3 ) then
                   !
                   if ( alobnd .and. allocated(subuc) ) deallocate(subuc,supuc,subus,supus,subvc,supvc,subvs,supvs,subzc,supzc,subzs,supzs)
                   !
                   if ( .not.allocated(subuc) ) allocate(subuc(nbgrpt,nfreq,kmax))
                   if ( .not.allocated(supuc) ) allocate(supuc(nbgrpt,nfreq,kmax))
                   if ( .not.allocated(subus) ) allocate(subus(nbgrpt,nfreq,kmax))
                   if ( .not.allocated(supus) ) allocate(supus(nbgrpt,nfreq,kmax))
                   !
                   if ( .not.allocated(subvc) ) allocate(subvc(nbgrpt,nfreq,kmax))
                   if ( .not.allocated(supvc) ) allocate(supvc(nbgrpt,nfreq,kmax))
                   if ( .not.allocated(subvs) ) allocate(subvs(nbgrpt,nfreq,kmax))
                   if ( .not.allocated(supvs) ) allocate(supvs(nbgrpt,nfreq,kmax))
                   !
                   if ( .not.allocated(subzc) ) allocate(subzc(nbgrpt,nfreq))
                   if ( .not.allocated(supzc) ) allocate(supzc(nbgrpt,nfreq))
                   if ( .not.allocated(subzs) ) allocate(subzs(nbgrpt,nfreq))
                   if ( .not.allocated(supzs) ) allocate(supzs(nbgrpt,nfreq))
                   !
                   subuc = 0.
                   supuc = 0.
                   subus = 0.
                   supus = 0.
                   !
                   subvc = 0.
                   supvc = 0.
                   subvs = 0.
                   supvs = 0.
                   !
                   subzc = 0.
                   supzc = 0.
                   subzs = 0.
                   supzs = 0.
                   !
                endif
                !
             endif
             !
          endif
          !
       else if ( iwvgen /= 0 ) then
          !
          if ( alobnd .and. allocated(kwave) ) deallocate(kwave)
          if ( .not.allocated(kwave) ) allocate( kwave(nfreq,1) )
          !
          kwave = 0.
          !
       endif
       !
    endif
    !
    alobnd = .false.
    !
    ! check whether boundary conditions have been given for each layer, if appropriate
    !
    jlay = 0
    do i = 1, nbgrpt
       if ( bgridp(9*i-1) > 0 ) jlay = jlay + int(bgridp(9*i-1),4)
    enddo
    sumlay = kmax*(kmax+1)/2
    if ( mod(jlay,sumlay) /= 0 ) then
       call msgerr (4, 'boundary conditions must be given for each layer')
       return
    endif
    !
    ! check roughness method in case of logarithmic velocity profile at open boundaries
    !
    logbnd = .false.
    do i = 1, nbgrpt
       if ( bgridp(9*i-1) == -2 ) logbnd = .true.
    enddo
    !
    if ( logbnd ) then
       if ( irough == 0 ) then
          call msgerr (3, 'roughness method must be given for logarithmic velocity profile at open boundaries')
          return
       else if ( irough == 11 ) then
          call msgerr (3, 'linear bottom friction not allowed in case of logarithmic velocity profile at open boundaries')
          return
       endif
    endif
    !
    ! check boundary conditions in case of mimetic discretizations
    !
    if ( mimetic ) then
       !
       do i = 1, nbgrpt
          if ( bgridp(9*i-7) ==  2 ) bgridp(9*i-7) = 7
          if ( bgridp(9*i-7) == 10 ) then
             call msgerr (3, 'outflow condition not supported in case of mimetic discretizations')
             exit
          endif
       enddo
       !
    endif
    !
    ! allocate arrays and store COG and DOF for a fixed body case (only used for output purposes)
    !
    if ( ifloat == 1 ) then
       !
       if ( .not.allocated(bcog) ) allocate(bcog(1,3)  )
       if ( .not.allocated(bdof) ) allocate(bdof(1,3,2))
       !
       bcog(1,1) = cogx
       bcog(1,2) = cogy
       bcog(1,3) = cogz
       !
       bdof = .false.
       !
    else if ( ifloat == 0 ) then
       !
       if ( .not.allocated(bcog) ) allocate(bcog(0,0)  )
       if ( .not.allocated(bdof) ) allocate(bdof(0,0,0))
       !
    endif
    !
    ! allocate arrays and store data for rigid bodies, if appropriate
    !
    if ( ifloat == 2 ) then
       !
       ! physical properties
       !
       if ( .not.allocated(bmass) ) allocate(bmass(mbod)       )
       if ( .not.allocated(bmoi ) ) allocate(bmoi (mbod,ndim)  )
       if ( .not.allocated(bcog ) ) allocate(bcog (mbod,ndim)  )
       if ( .not.allocated(bdof ) ) allocate(bdof (mbod,ndim,2))
       if ( .not.allocated(barea) ) allocate(barea(mbod)       )
       if ( .not.allocated(bvol ) ) allocate(bvol (mbod)       )
       if ( .not.allocated(bmoax) ) allocate(bmoax(mbod)       )
       if ( .not.allocated(bmoay) ) allocate(bmoay(mbod)       )
       if ( .not.allocated(bcob ) ) allocate(bcob (mbod)       )
       !
       ! power of PTO system
       !
       if ( .not.allocated(ptop) ) allocate(ptop(mbod))
       !
       ! mooring lines
       !
       if ( .not.allocated(nmlb) ) allocate(nmlb(mbod))
       if ( mxmli > 0 ) then
          !
          if ( .not.allocated(bmli  ) ) allocate(bmli  (mbod,mxmli,10))
          if ( .not.allocated(cpto  ) ) allocate(cpto  (mbod)         )
          !
          if ( .not.allocated(extml0) ) allocate(extml0(mbod,mxmli)   )
          if ( .not.allocated(extml1) ) allocate(extml1(mbod,mxmli)   )
          !
       endif
       !
       ! fenders
       !
       if ( .not.allocated(nfeb) ) allocate(nfeb(mbod))
       if ( mxfen > 0 ) then
          !
          if ( .not.allocated(bfen) ) allocate(bfen(mbod,mxfen,4))
          !
       endif
       !
       bmass = 0.
       bmoi  = 0.
       bcog  = 0.
       bdof  = .false.
       barea = 0.
       bvol  = 0.
       bmoax = 0.
       bmoay = 0.
       bcob  = 0.
       ptop  = 0.
       nmlb  = 0
       nfeb  = 0
       if ( mxmli > 0 ) then
          bmli   = 0.
          cpto   = .false.
          extml0 = 0.
          extml1 = 0.
       endif
       if ( mxfen > 0 ) then
          bfen = 0.
       endif
       !
       currib => frigbod
       do i = 1, mbod
          m = currib%label
          bmass(m) = currib%parm(1)
          do j = 1, ndim
             bmoi(m,j) = currib%parm(j+1)
             bcog(m,j) = currib%parm(j+4)
             do k = 1, 2
                l = j + ndim*(k-1)
                bdof(m,j,k) = currib%dof(l)
             enddo
          enddo
          if ( mxmli > 0 ) then
             nmli = currib%nmli
             do j = 1, nmli
                do k = 1, 9
                   bmli(m,j,k) = currib%mli(j,k)
                enddo
             enddo
             nmlb(m) = nmli
             cpto(m) = currib%pret
          endif
          if ( mxfen > 0 ) then
             nfen = currib%nfen
             do j = 1, nfen
                do k = 1, 4
                   bfen(m,j,k) = currib%fen(j,k)
                enddo
             enddo
             nfeb(m) = nfen
          endif
          if (.not.associated(currib%nextrib)) exit
          currib => currib%nextrib
       enddo
       !
       if ( oned ) then
          !
          ! 1D mode - three degrees of freedom: surge, heave and pitch
          !
          bdof(:,2,1) = .false.
          bdof(:,1,2) = .false.
          bdof(:,3,2) = .false.
          !
       endif
       !
    else
       !
       if ( .not.allocated(bmass ) ) allocate(bmass (0)    )
       if ( .not.allocated(bmoi  ) ) allocate(bmoi  (0,0)  )
       if ( .not.allocated(barea ) ) allocate(barea (0)    )
       if ( .not.allocated(bvol  ) ) allocate(bvol  (0)    )
       if ( .not.allocated(bmoax ) ) allocate(bmoax (0)    )
       if ( .not.allocated(bmoay ) ) allocate(bmoay (0)    )
       if ( .not.allocated(bcob  ) ) allocate(bcob  (0)    )
       if ( .not.allocated(ptop  ) ) allocate(ptop  (0)    )
       if ( .not.allocated(nmlb  ) ) allocate(nmlb  (0)    )
       if ( .not.allocated(bmli  ) ) allocate(bmli  (0,0,0))
       if ( .not.allocated(cpto  ) ) allocate(cpto  (0)    )
       if ( .not.allocated(extml0) ) allocate(extml0(0,0)  )
       if ( .not.allocated(extml1) ) allocate(extml1(0,0)  )
       if ( .not.allocated(nfeb  ) ) allocate(nfeb  (0)    )
       if ( .not.allocated(bfen  ) ) allocate(bfen  (0,0,0))
       !
    endif
    !
    if ( ifloat == 1 .or. ifloat == 2 ) then
       !
       bcog(:,1) = bcog(:,1) - xoffs
       bcog(:,2) = bcog(:,2) - yoffs
       !
    endif
    !
    ! allocate work arrays to temporarily store data
    !
    if ( lwork /= 0 ) then
       if ( optg /= 5 ) then
          if ( .not.allocated(iwrk) ) allocate(iwrk(mcgrd,lwork))
          if ( .not.allocated(work) ) allocate(work(mcgrd,lwork))
       else
          if ( .not.allocated(iwrk) ) allocate(iwrk(ncells,lwork))
          if ( .not.allocated(work) ) allocate(work(ncells,lwork))
          if ( .not.allocated(uwrk) ) allocate(uwrk(nfaces))
       endif
    else
       if ( .not.allocated(iwrk) ) allocate(iwrk(0,0))
       if ( .not.allocated(work) ) allocate(work(0,0))
       if ( .not.allocated(uwrk) ) allocate(uwrk(0  ))
    endif
    !
    ! allocate arrays oparr and oparrk for outputting, if appropriate
    !
    if ( nreoq /= 0 ) then
       if ( optg /= 5 ) then
          if ( .not.allocated(oparr) ) allocate(oparr(mcgrd,mopa))
          if ( .not.allocated(oparrk) .and. kmax > 1 ) allocate(oparrk(mcgrd,0:kmax,mopak))
       else
          if ( .not.allocated(oparr) ) allocate(oparr(nverts,mopa))
          if ( .not.allocated(oparrk) .and. kmax > 1 ) allocate(oparrk(nverts,0:kmax,mopak))
       endif
    else
       if ( .not.allocated(oparr ) ) allocate(oparr (0,0  ))
       if ( .not.allocated(oparrk) ) allocate(oparrk(0,0,0))
    endif
    !
    if ( optg == 5 ) then
       !
       ! deallocate arrays kvertc and kvertf (we don't use them anymore!)
       !
       if (allocated(kvertc)) deallocate(kvertc)
       if (allocated(kvertf)) deallocate(kvertf)
       !
    endif
    !
    ! allocate and fill input fields by mapping onto computational grid
    !
    istat = 0
    if ( .not.allocated(depf) ) allocate ( depf(mcgrd), stat = istat )
    if ( istat == 0 ) depf = 0.
    !
    if ( istat == 0 ) then
       !
       if ( igtype(7) /= 0 ) then
          !
          if ( .not.allocated(wlevf) ) allocate ( wlevf(mcgrd,1:3), stat = istat )
          if ( istat == 0 ) wlevf = 0.
          !
       else
          !
          if ( .not.allocated(wlevf) ) allocate ( wlevf(0,0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( igtype(2) /= 0 ) then
          !
          if ( .not.allocated(uxf) ) allocate ( uxf(mcgrd,1:3), stat = istat )
          if ( istat == 0 ) then
             uxf = 0.
             if ( .not.allocated(uyf) ) allocate ( uyf(mcgrd,1:3), stat = istat )
             if ( istat == 0 ) uyf = 0.
          endif
          !
       else
          !
          if ( .not.allocated(uxf) ) allocate ( uxf(0,0) )
          if ( .not.allocated(uyf) ) allocate ( uyf(0,0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( varfr ) then
          !
          if ( .not.allocated(fricf) ) allocate ( fricf(mcgrd,1:3), stat = istat )
          if ( istat == 0 ) fricf = 0.
          !
       else
          !
          if ( .not.allocated(fricf) ) allocate ( fricf(0,0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( varwi ) then
          !
          if ( .not.allocated(wxf) ) allocate ( wxf(mcgrd,1:3), stat = istat )
          if ( istat == 0 ) then
             wxf = 0.
             if ( .not.allocated(wyf) ) allocate ( wyf(mcgrd,1:3), stat = istat )
             if ( istat == 0 ) wyf = 0.
          endif
          !
       else
          !
          if ( .not.allocated(wxf) ) allocate ( wxf(0,0) )
          if ( .not.allocated(wyf) ) allocate ( wyf(0,0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( svwp ) then
          !
          if ( .not.allocated(presf) ) allocate ( presf(mcgrd,1:3), stat = istat )
          if ( istat == 0 ) presf = 0.
          !
       else
          !
          if ( .not.allocated(presf) ) allocate ( presf(0,0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( iporos == 1 ) then
          !
          if ( .not.allocated(nporf) ) allocate ( nporf(mcgrd), stat = istat )
          if ( istat == 0 ) nporf = 1.
          !
       else
          !
          if ( .not.allocated(nporf) ) allocate ( nporf(0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( vargs ) then
          !
          if ( .not.allocated(gsizf) ) allocate ( gsizf(mcgrd), stat = istat )
          if ( istat == 0 ) gsizf = 0.
          !
       else
          !
          if ( .not.allocated(gsizf) ) allocate ( gsizf(0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( varsh ) then
          !
          if ( .not.allocated(hporf) ) allocate ( hporf(mcgrd), stat = istat )
          if ( istat == 0 ) hporf = ppor(2)
          !
       else
          !
          if ( .not.allocated(hporf) ) allocate ( hporf(0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( varnpl ) then
          !
          if ( .not.allocated(nplaf) ) allocate ( nplaf(mcgrd), stat = istat )
          if ( istat == 0 ) nplaf = 0.
          !
       else
          !
          if ( .not.allocated(nplaf) ) allocate ( nplaf(0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( lsal > 0 ) then
          !
          if ( .not.allocated(salf) ) allocate ( salf(mcgrd,kmax), stat = istat )
          if ( istat == 0 ) salf = 0.
          !
       else
          !
          if ( .not.allocated(salf) ) allocate ( salf(0,0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( ltemp > 0 ) then
          !
          if ( .not.allocated(tempf) ) allocate ( tempf(mcgrd,kmax), stat = istat )
          if ( istat == 0 ) tempf = 0.
          !
       else
          !
          if ( .not.allocated(tempf) ) allocate ( tempf(0,0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( lsed > 0 ) then
          !
          if ( .not.allocated(sedf) ) allocate ( sedf(mcgrd,kmax), stat = istat )
          if ( istat == 0 ) sedf = 0.
          !
       else
          !
          if ( .not.allocated(sedf) ) allocate ( sedf(0,0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( ifloat /= 0 ) then
          !
          if ( .not.allocated(flobjf) ) allocate ( flobjf(mcgrd), stat = istat )
          if ( istat == 0 ) flobjf = 0.
          !
       else
          !
          if ( .not.allocated(flobjf) ) allocate ( flobjf(0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( ifloat == 2 ) then
          !
          if ( .not.allocated(lbodf) ) allocate ( lbodf(mcgrd), stat = istat )
          if ( istat == 0 ) lbodf = 0.
          !
       else
          !
          if ( .not.allocated(lbodf) ) allocate ( lbodf(0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( varva ) then
          !
          if ( .not.allocated(avxf) ) allocate ( avxf(mcgrd,kpmax), stat = istat )
          if ( istat == 0 ) then
             avxf = 0.
             if ( .not.allocated(avyf) ) allocate ( avyf(mcgrd,kpmax), stat = istat )
             if ( istat == 0 ) avyf = 0.
          endif
          !
       else
          !
          if ( .not.allocated(avxf) ) allocate ( avxf(0,0) )
          if ( .not.allocated(avyf) ) allocate ( avyf(0,0) )
          !
       endif
       !
    endif
    !
    if ( istat == 0 ) then
       !
       if ( leds(23) == 2 ) then
          !
          if ( .not.allocated(corpf) ) allocate ( corpf(mcgrd), stat = istat )
          if ( istat == 0 ) corpf = 0.
          !
       else
          !
          if ( .not.allocated(corpf) ) allocate ( corpf(0) )
          !
       endif
       !
    endif
    !
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: input fields at computational grid and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    if ( optg /= 5 ) then
       !
       ! structured grid
       !
       do iy = 1, myc - 1
          do ix = 1, mxc - 1
             !
             indx = kgrpnt(ix,iy)
             !
             if ( indx > 1 ) then
                !
                ! take global indices instead of local ones
                !
                ixgl = ix + MXF - 1
                iygl = iy + MYF - 1
                !
                xp = xcgrid(ix,iy)
                yp = ycgrid(ix,iy)
                !
                ! retrieve bottom and water level and map to computational grid
                !
                if ( lflgrd(1) ) then
                   dep = depth((iygl-1)*mxg(1)+ixgl)
                else
                   dep = SVALQI (xp, yp, 1, depth, NEAREST)
                endif
                depf(indx) = dep + swl
                !
                if ( igtype(7) /= 0 ) then
                   if ( lflgrd(7) ) then
                      wlvl = wlevl((iygl-1)*mxg(7)+ixgl)
                   else
                      wlvl = SVALQI (xp, yp, 7, wlevl, NEAREST)
                   endif
                   wlevf(indx,2) = wlvl
                endif
                !
                ! retrieve flow velocities and map to computational grid
                !
                if ( igtype(2) /= 0 ) then
                   if ( dep + wlvl > epsdry ) then
                      if ( lflgrd(2) ) then
                         uu = uxb((iygl-1)*mxg(2)+ixgl)
                      else
                         uu = SVALQI (xp, yp, 2, uxb, 0.)
                      endif
                      if ( lflgrd(3) ) then
                         vv = uyb((iygl-1)*mxg(3)+ixgl)
                      else
                         vv = SVALQI (xp, yp, 3, uyb, 0.)
                      endif
                      uxf(indx,2) =  uu*cosvc + vv*sinvc
                      uyf(indx,2) = -uu*sinvc + vv*cosvc
                   else
                      uxf(indx,2) = 0.
                      uyf(indx,2) = 0.
                   endif
                endif
                !
                ! retrieve space varying bottom friction coefficient and map to computational grid
                !
                if ( varfr ) then
                   if ( lflgrd(4) ) then
                      cf = fric((iygl-1)*mxg(4)+ixgl)
                   else
                      cf = SVALQI (xp, yp, 4, fric, NEAREST)
                   endif
                   fricf(indx,2) = cf
                endif
                !
                ! retrieve space varying wind velocities and map to computational grid
                !
                if ( varwi ) then
                   if ( lflgrd(5) ) then
                      uu = wxi((iygl-1)*mxg(5)+ixgl)
                   else
                      uu = SVALQI (xp, yp, 5, wxi, 0.)
                   endif
                   if ( lflgrd(6) ) then
                      vv = wyi((iygl-1)*mxg(6)+ixgl)
                   else
                      vv = SVALQI (xp, yp, 6, wyi, 0.)
                   endif
                   wxf(indx,2) =  uu*coswc + vv*sinwc
                   wyf(indx,2) = -uu*sinwc + vv*coswc
                endif
                !
                ! retrieve space varying atmospheric pressure and map to computational grid
                !
                if ( svwp ) then
                   if ( lflgrd(13) ) then
                      pp = pres((iygl-1)*mxg(13)+ixgl)
                   else
                      pp = SVALQI (xp, yp, 13, pres, NEAREST)
                   endif
                   presf(indx,2) = pp
                endif
                !
                ! retrieve porosity and map to computational grid
                !
                if ( iporos == 1 ) then
                   if ( lflgrd(10) ) then
                      n = npor((iygl-1)*mxg(10)+ixgl)
                   else
                      n = SVALQI (xp, yp, 10, npor, 1.)
                   endif
                   nporf(indx) = max( n, 0.1 )
                endif
                !
                ! retrieve space varying grain sizes and map to computational grid
                !
                if ( vargs ) then
                   if ( lflgrd(11) ) then
                      dsiz = gsiz((iygl-1)*mxg(11)+ixgl)
                   else
                      dsiz = SVALQI (xp, yp, 11, gsiz, NEAREST)
                   endif
                   gsizf(indx) = dsiz
                endif
                !
                ! retrieve space varying structure heights and map to computational grid
                !
                if ( varsh ) then
                   if ( lflgrd(12) ) then
                      hps = hpor((iygl-1)*mxg(12)+ixgl)
                   else
                      hps = SVALQI (xp, yp, 12, hpor, ppor(2))
                   endif
                   hporf(indx) = hps
                endif
                !
                ! retrieve space varying number of plants and map to computational grid
                !
                if ( varnpl ) then
                   if ( lflgrd(14) ) then
                      npl = npla((iygl-1)*mxg(14)+ixgl)
                   else
                      npl = SVALQI (xp, yp, 14, npla, NEAREST)
                   endif
                   nplaf(indx) = npl
                endif
                !
                ! retrieve salinity and map to computational grid
                !
                if ( lsal > 0 ) then
                   if ( ifllay(15) == 1 ) then
                      if ( lflgrd(15) ) then
                         ss = sal((iygl-1)*mxg(15)+ixgl)
                      else
                         ss = SVALQI (xp, yp, 15, sal, NEAREST)
                      endif
                      salf(indx,1) = ss
                      if ( kmax > 1 ) then
                         do k = 2, kmax
                            salf(indx,k) = ss
                         enddo
                      endif
                   else
                      do k = 1, kmax
                         if ( lflgrd(15) ) then
                            ss = sal((iygl-1+(k-1)*myg(15))*mxg(15)+ixgl)
                         else
                            jlow = 1+(k-1)*mxg(15)*myg(15)
                            jupp = k*mxg(15)*myg(15)
                            ss = SVALQI (xp, yp, 15, sal(jlow:jupp), NEAREST)
                         endif
                         salf(indx,k) = ss
                      enddo
                   endif
                endif
                !
                ! retrieve temperature and map to computational grid
                !
                if ( ltemp > 0 ) then
                   if ( ifllay(16) == 1 ) then
                      if ( lflgrd(16) ) then
                         tt = temp((iygl-1)*mxg(16)+ixgl)
                      else
                         tt = SVALQI (xp, yp, 16, temp, NEAREST)
                      endif
                      tempf(indx,1) = tt
                      if ( kmax > 1 ) then
                         do k = 2, kmax
                            tempf(indx,k) = tt
                         enddo
                      endif
                   else
                      do k = 1, kmax
                         if ( lflgrd(16) ) then
                            tt = temp((iygl-1+(k-1)*myg(16))*mxg(16)+ixgl)
                         else
                            jlow = 1+(k-1)*mxg(16)*myg(16)
                            jupp = k*mxg(16)*myg(16)
                            tt = SVALQI (xp, yp, 16, temp(jlow:jupp), NEAREST)
                         endif
                         tempf(indx,k) = tt
                      enddo
                   endif
                endif
                !
                ! retrieve suspended sediment and map to computational grid
                !
                if ( lsed > 0 ) then
                   if ( ifllay(17) == 1 ) then
                      if ( lflgrd(17) ) then
                         ss = sed((iygl-1)*mxg(17)+ixgl)
                      else
                         ss = SVALQI (xp, yp, 17, sed, NEAREST)
                      endif
                      sedf(indx,1) = ss
                      if ( kmax > 1 ) then
                         do k = 2, kmax
                            sedf(indx,k) = ss
                         enddo
                      endif
                   else
                      do k = 1, kmax
                         if ( lflgrd(17) ) then
                            ss = sed((iygl-1+(k-1)*myg(17))*mxg(17)+ixgl)
                         else
                            jlow = 1+(k-1)*mxg(17)*myg(17)
                            jupp = k*mxg(17)*myg(17)
                            ss = SVALQI (xp, yp, 17, sed(jlow:jupp), NEAREST)
                         endif
                         sedf(indx,k) = ss
                      enddo
                   endif
                endif
                !
                ! retrieve draft of floating object and map to computational grid
                !
                if ( ifloat /= 0 ) then
                   if ( lflgrd(18) ) then
                      draft = flobj((iygl-1)*mxg(18)+ixgl)
                   else
                      draft = SVALQI (xp, yp, 18, flobj, NEAREST)
                   endif
                   flobjf(indx) = draft
                endif
                !
                ! retrieve label of rigid body and map to computational grid
                !
                if ( ifloat == 2 ) then
                   if ( allocated(lbod) ) then
                      if ( lflgrd(19) ) then
                         labl = lbod((iygl-1)*mxg(19)+ixgl)
                      else
                         labl = SVALQI (xp, yp, 19, lbod, 0.)
                      endif
                      lbodf(indx) = labl
                   else
                      if ( flobjf(indx) > 0. ) then
                         lbodf(indx) = 1.
                      else
                         lbodf(indx) = 0.
                      endif
                   endif
                endif
                !
                ! retrieve space varying ambient velocities and map to computational grid
                !
                if ( varva ) then
                   if ( ifllay(20) == 1 ) then
                      if ( lflgrd(20) ) then
                         uu = avxi((iygl-1)*mxg(20)+ixgl)
                      else
                         uu = SVALQI (xp, yp, 20, avxi, 0.)
                      endif
                      if ( lflgrd(21) ) then
                         vv = avyi((iygl-1)*mxg(21)+ixgl)
                      else
                         vv = SVALQI (xp, yp, 21, avyi, 0.)
                      endif
                      avxf(indx,1) =  uu*cosva + vv*sinva
                      avyf(indx,1) = -uu*sinva + vv*cosva
                      if ( kpmax > 1 ) then
                         do k = 2, kpmax
                            avxf(indx,k) = avxf(indx,1)
                            avyf(indx,k) = avyf(indx,1)
                         enddo
                      endif
                   else
                      do k = 1, kpmax
                         if ( lflgrd(20) ) then
                            uu = avxi((iygl-1+(k-1)*myg(20))*mxg(20)+ixgl)
                         else
                            jlow = 1+(k-1)*mxg(20)*myg(20)
                            jupp = k*mxg(20)*myg(20)
                            uu = SVALQI (xp, yp, 20, avxi(jlow:jupp), 0.)
                         endif
                         if ( lflgrd(21) ) then
                            vv = avyi((iygl-1+(k-1)*myg(21))*mxg(21)+ixgl)
                         else
                            jlow = 1+(k-1)*mxg(21)*myg(21)
                            jupp = k*mxg(21)*myg(21)
                            vv = SVALQI (xp, yp, 21, avyi(jlow:jupp), 0.)
                         endif
                         avxf(indx,k) =  uu*cosva + vv*sinva
                         avyf(indx,k) = -uu*sinva + vv*cosva
                      enddo
                   endif
                endif
                !
                ! retrieve space varying mean water level (tide/wind surge), map to computational grid and finally add to bottom elevation
                !
                if ( varea ) then
                   if ( lflgrd(22) ) then
                      wlvl = mwl((iygl-1)*mxg(22)+ixgl)
                   else
                      wlvl = SVALQI (xp, yp, 22, mwl, NEAREST)
                   endif
                   depf(indx) = depf(indx) + wlvl
                else
                   depf(indx) = depf(indx) + eamb
                endif
                !
                ! retrieve space varying Coriolis factor and map to computational grid
                !
                if ( leds(23) == 2 ) then
                   if ( lflgrd(23) ) then
                      cf = corp((iygl-1)*mxg(23)+ixgl)
                   else
                      cf = SVALQI (xp, yp, 23, corp, NEAREST)
                   endif
                   corpf(indx) = cf
                endif
                !
             endif
             !
          enddo
       enddo
       !
    else
       !
       ! unstructured grid
       !
       do indx = 1, nverts
          !
          xp = xcugrd(indx)
          yp = ycugrd(indx)
          !
          ! retrieve bottom and water level and map to computational grid
          !
          if ( lflgrd(1) ) then
             dep = depth(indx)
          else
             dep = SVALQI (xp, yp, 1, depth, NEAREST)
          endif
          depf(indx) = dep + swl
          !
          if ( igtype(7) /= 0 ) then
             if ( lflgrd(7) ) then
                wlvl = wlevl(indx)
             else
                wlvl = SVALQI (xp, yp, 7, wlevl, NEAREST)
             endif
             wlevf(indx,2) = wlvl
          endif
          !
          ! retrieve flow velocities and map to computational grid
          !
          if ( igtype(2) /= 0 ) then
             if ( dep + wlvl > epsdry ) then
                if ( lflgrd(2) ) then
                   uu = uxb(indx)
                else
                   uu = SVALQI (xp, yp, 2, uxb, 0.)
                endif
                if ( lflgrd(3) ) then
                   vv = uyb(indx)
                else
                   vv = SVALQI (xp, yp, 3, uyb, 0.)
                endif
                uxf(indx,2) =  uu*cosvc + vv*sinvc
                uyf(indx,2) = -uu*sinvc + vv*cosvc
             else
                uxf(indx,2) = 0.
                uyf(indx,2) = 0.
             endif
          endif
          !
          ! retrieve space varying bottom friction coefficient and map to computational grid
          !
          if ( varfr ) then
             if ( lflgrd(4) ) then
                cf = fric(indx)
             else
                cf = SVALQI (xp, yp, 4, fric, NEAREST)
             endif
             fricf(indx,2) = cf
          endif
          !
          ! retrieve space varying wind velocities and map to computational grid
          !
          if ( varwi ) then
             if ( lflgrd(5) ) then
                uu = wxi(indx)
             else
                uu = SVALQI (xp, yp, 5, wxi, 0.)
             endif
             if ( lflgrd(6) ) then
                vv = wyi(indx)
             else
                vv = SVALQI (xp, yp, 6, wyi, 0.)
             endif
             wxf(indx,2) =  uu*coswc + vv*sinwc
             wyf(indx,2) = -uu*sinwc + vv*coswc
          endif
          !
          ! retrieve space varying atmospheric pressure and map to computational grid
          !
          if ( svwp ) then
             if ( lflgrd(13) ) then
                pp = pres(indx)
             else
                pp = SVALQI (xp, yp, 13, pres, NEAREST)
             endif
             presf(indx,2) = pp
          endif
          !
          ! retrieve porosity and map to computational grid
          !
          if ( iporos == 1 ) then
             if ( lflgrd(10) ) then
                n = npor(indx)
             else
                n = SVALQI (xp, yp, 10, npor, 1.)
             endif
             nporf(indx) = max( n, 0.1 )
          endif
          !
          ! retrieve space varying grain sizes and map to computational grid
          !
          if ( vargs ) then
             if ( lflgrd(11) ) then
                dsiz = gsiz(indx)
             else
                dsiz = SVALQI (xp, yp, 11, gsiz, NEAREST)
             endif
             gsizf(indx) = dsiz
          endif
          !
          ! retrieve space varying structure heights and map to computational grid
          !
          if ( varsh ) then
             if ( lflgrd(12) ) then
                hps = hpor(indx)
             else
                hps = SVALQI (xp, yp, 12, hpor, ppor(2))
             endif
             hporf(indx) = hps
          endif
          !
          ! retrieve space varying number of plants and map to computational grid
          !
          if ( varnpl ) then
             if ( lflgrd(14) ) then
                npl = npla(indx)
             else
                npl = SVALQI (xp, yp, 14, npla, NEAREST)
             endif
             nplaf(indx) = npl
          endif
          !
          ! retrieve salinity and map to computational grid
          !
          if ( lsal > 0 ) then
             if ( ifllay(15) == 1 ) then
                if ( lflgrd(15) ) then
                   ss = sal(indx)
                else
                   ss = SVALQI (xp, yp, 15, sal, NEAREST)
                endif
                salf(indx,1) = ss
                if ( kmax > 1 ) then
                   do k = 2, kmax
                      salf(indx,k) = ss
                   enddo
                endif
             else
                do k = 1, kmax
                   if ( lflgrd(15) ) then
                      ss = sal(indx+(k-1)*nverts)
                   else
                      jlow = 1+(k-1)*mxg(15)*myg(15)
                      jupp = k*mxg(15)*myg(15)
                      ss = SVALQI (xp, yp, 15, sal(jlow:jupp), NEAREST)
                   endif
                   salf(indx,k) = ss
                enddo
             endif
          endif
          !
          ! retrieve temperature and map to computational grid
          !
          if ( ltemp > 0 ) then
             if ( ifllay(16) == 1 ) then
                if ( lflgrd(16) ) then
                   tt = temp(indx)
                else
                   tt = SVALQI (xp, yp, 16, temp, NEAREST)
                endif
                tempf(indx,1) = tt
                if ( kmax > 1 ) then
                   do k = 2, kmax
                      tempf(indx,k) = tt
                   enddo
                endif
             else
                do k = 1, kmax
                   if ( lflgrd(16) ) then
                      tt = temp(indx+(k-1)*nverts)
                   else
                      jlow = 1+(k-1)*mxg(16)*myg(16)
                      jupp = k*mxg(16)*myg(16)
                      tt = SVALQI (xp, yp, 16, temp(jlow:jupp), NEAREST)
                   endif
                   tempf(indx,k) = tt
                enddo
             endif
          endif
          !
          ! retrieve suspended sediment and map to computational grid
          !
          if ( lsed > 0 ) then
             if ( ifllay(17) == 1 ) then
                if ( lflgrd(17) ) then
                   ss = sed(indx)
                else
                   ss = SVALQI (xp, yp, 17, sed, NEAREST)
                endif
                sedf(indx,1) = ss
                if ( kmax > 1 ) then
                   do k = 2, kmax
                      sedf(indx,k) = ss
                   enddo
                endif
             else
                do k = 1, kmax
                   if ( lflgrd(17) ) then
                      ss = sed(indx+(k-1)*nverts)
                   else
                      jlow = 1+(k-1)*mxg(17)*myg(17)
                      jupp = k*mxg(17)*myg(17)
                      ss = SVALQI (xp, yp, 17, sed(jlow:jupp), NEAREST)
                   endif
                   sedf(indx,k) = ss
                enddo
             endif
          endif
          !
          ! retrieve draft of floating object and map to computational grid
          !
          if ( ifloat /= 0 ) then
             if ( lflgrd(18) ) then
                draft = flobj(indx)
             else
                draft = SVALQI (xp, yp, 18, flobj, NEAREST)
             endif
             flobjf(indx) = draft
          endif
          !
          ! retrieve label of rigid body and map to computational grid
          !
          if ( ifloat == 2 ) then
             if ( allocated(lbod) ) then
                if ( lflgrd(19) ) then
                   labl = lbod(indx)
                else
                   labl = SVALQI (xp, yp, 19, lbod, 0.)
                endif
                lbodf(indx) = labl
             else
                if ( flobjf(indx) > 0. ) then
                   lbodf(indx) = 1.
                else
                   lbodf(indx) = 0.
                endif
             endif
          endif
          !
          ! retrieve space varying ambient velocities and map to computational grid
          !
          if ( varva ) then
             if ( ifllay(20) == 1 ) then
                if ( lflgrd(20) ) then
                   uu = avxi(indx)
                else
                   uu = SVALQI (xp, yp, 20, avxi, 0.)
                endif
                if ( lflgrd(21) ) then
                   vv = avyi(indx)
                else
                   vv = SVALQI (xp, yp, 21, avyi, 0.)
                endif
                avxf(indx,1) =  uu*cosva + vv*sinva
                avyf(indx,1) = -uu*sinva + vv*cosva
                if ( kpmax > 1 ) then
                   do k = 2, kpmax
                      avxf(indx,k) = avxf(indx,1)
                      avyf(indx,k) = avyf(indx,1)
                   enddo
                endif
             else
                do k = 1, kpmax
                   if ( lflgrd(20) ) then
                      uu = avxi(indx+(k-1)*nverts)
                   else
                      jlow = 1+(k-1)*mxg(20)*myg(20)
                      jupp = k*mxg(20)*myg(20)
                      uu = SVALQI (xp, yp, 20, avxi(jlow:jupp), 0.)
                   endif
                   if ( lflgrd(21) ) then
                      vv = avyi(indx+(k-1)*nverts)
                   else
                      jlow = 1+(k-1)*mxg(21)*myg(21)
                      jupp = k*mxg(21)*myg(21)
                      vv = SVALQI (xp, yp, 21, avyi(jlow:jupp), 0.)
                   endif
                   avxf(indx,k) =  uu*cosva + vv*sinva
                   avyf(indx,k) = -uu*sinva + vv*cosva
                enddo
             endif
          endif
          !
          ! retrieve space varying mean water level (tide/wind surge), map to computational grid and finally add to bottom elevation
          !
          if ( varea ) then
             if ( lflgrd(22) ) then
                wlvl = mwl(indx)
             else
                wlvl = SVALQI (xp, yp, 22, mwl, NEAREST)
             endif
             depf(indx) = depf(indx) + wlvl
          else
             depf(indx) = depf(indx) + eamb
          endif
          !
          ! retrieve space varying Coriolis factor and map to computational grid
          !
          if ( leds(23) == 2 ) then
             if ( lflgrd(23) ) then
                cf = corp(indx)
             else
                cf = SVALQI (xp, yp, 23, corp, NEAREST)
             endif
             corpf(indx) = cf
          endif
          !
       enddo
       !
    endif
    !
    ! initialize flow variables based on space varying input fields
    !
    if ( optg /= 5 ) then
       call SwashUpdFlowFlds
    else
       call SwashUpdUFlowFlds
    endif
    if (STPNOW()) return
    !
    ! allocate and initialize flow data (friction, wind, pressure, depths and mask arrays for wetting and drying)
    !
    if ( optg /= 5 ) then
       if ( .not.allocated(dps) ) allocate(dps(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(dpu) ) allocate (dpu(mcgrd), stat = istat)
    else
       if ( .not.allocated(dps) ) allocate(dps(ncells), stat = istat)
       if ( istat == 0 .and. .not.allocated(dpu) ) allocate (dpu(nfaces), stat = istat)
    endif
    !
    if ( optg /= 5 ) then
       if ( istat == 0 .and. .not.allocated(hs ) ) allocate (hs (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(hso) ) allocate (hso(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(hu ) ) allocate (hu (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(hum) ) allocate (hum(mcgrd), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(hs ) ) allocate (hs (ncells), stat = istat)
       if ( istat == 0 .and. .not.allocated(hso) ) allocate (hso(ncells), stat = istat)
       if ( istat == 0 .and. .not.allocated(hu ) ) allocate (hu (nfaces), stat = istat)
       if ( istat == 0 .and. .not.allocated(hum) ) allocate (hum(nfaces), stat = istat)
    endif
    if ( optg /= 5 ) then
       if ( istat == 0 .and. .not.allocated(humn) ) allocate (humn(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(humo) ) allocate (humo(mcgrd), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(humn) ) allocate (humn(nfaces), stat = istat)
       if ( istat == 0 .and. .not.allocated(humo) ) allocate (humo(nfaces), stat = istat)
    endif
    !
    if ( optg /= 5 ) then
       if ( istat == 0 .and. .not.allocated(cfricu) ) allocate (cfricu(mcgrd), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(cfricu) ) allocate (cfricu(ncells), stat = istat)
    endif
    if ( istat == 0 .and. .not.allocated(apors ) ) allocate (apors (mcgrd), stat = istat)
    if ( istat == 0 .and. .not.allocated(bpors ) ) allocate (bpors (mcgrd), stat = istat)
    if ( istat == 0 .and. .not.allocated(cpors ) ) allocate (cpors (mcgrd), stat = istat)
    if ( istat == 0 .and. .not.allocated(aporu ) ) allocate (aporu (mcgrd), stat = istat)
    if ( istat == 0 .and. .not.allocated(bporu ) ) allocate (bporu (mcgrd), stat = istat)
    if ( istat == 0 .and. .not.allocated(cporu ) ) allocate (cporu (mcgrd), stat = istat)
    !
    if ( optg /= 5 ) then
       if ( istat == 0 .and. .not.allocated(windu) ) allocate (windu(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(cwndu) ) allocate (cwndu(mcgrd), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(windu) ) allocate (windu(nfaces), stat = istat)
       if ( istat == 0 .and. .not.allocated(cwndu) ) allocate (cwndu(nfaces), stat = istat)
    endif
    !
    if ( optg /= 5 ) then
       if ( istat == 0 .and. .not.allocated(flos) ) allocate (flos(mcgrd), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(flos) ) allocate (flos(ncells), stat = istat)
    endif
    if ( istat == 0 .and. .not.allocated(flou) ) allocate (flou(mcgrd), stat = istat)
    !
    if ( istat == 0 .and. .not.allocated(lfbs) ) allocate (lfbs(mcgrd), stat = istat)
    !
    if ( svwp ) then
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(patm) ) allocate (patm(mcgrd), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(patm) ) allocate (patm(ncells), stat = istat)
       endif
    endif
    !
    if ( optg /= 5 ) then
       if ( istat == 0 .and. .not.allocated(wets) ) allocate (wets(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(wetu) ) allocate (wetu(mcgrd), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(wets) ) allocate (wets(ncells), stat = istat)
       if ( istat == 0 .and. .not.allocated(wetu) ) allocate (wetu(nfaces), stat = istat)
    endif
    !
    if ( istat == 0 .and. .not.allocated(presp) ) allocate (presp(mcgrd), stat = istat)
    if ( istat == 0 .and. .not.allocated(presu) ) allocate (presu(mcgrd), stat = istat)
    !
    if ( istat == 0 .and. .not.allocated(skc) ) allocate (skc(mcgrd), stat = istat)
    !
    if ( optg /= 5 ) then
       if ( istat == 0 .and. .not.allocated(brks) ) allocate (brks(mcgrd), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(brks) ) allocate (brks(ncells), stat = istat)
    endif
    !
    if ( istat == 0 .and. .not.allocated(npors) ) allocate (npors(mcgrd), stat = istat)
    if ( istat == 0 .and. .not.allocated(nporu) ) allocate (nporu(mcgrd), stat = istat)
    if ( istat == 0 .and. .not.allocated(hpors) ) allocate (hpors(mcgrd), stat = istat)
    if ( istat == 0 .and. .not.allocated(hporu) ) allocate (hporu(mcgrd), stat = istat)
    !
    if ( istat == 0 .and. .not.allocated(cvegu) ) allocate (cvegu(mcgrd,kmax,2), stat = istat)
    !
    if ( optg /= 5 ) then
       if ( istat == 0 .and. .not.allocated(hindun) ) allocate (hindun(mcgrd), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(hindun) ) allocate (hindun(ncells), stat = istat)
    endif
    !
    if ( relwav ) then
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(smax) ) allocate (smax(mcgrd), stat = istat)
          if ( istat == 0 .and. .not.allocated(so  ) ) allocate (so  (mcgrd), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(smax) ) allocate (smax(nfaces), stat = istat)
          if ( istat == 0 .and. .not.allocated(so  ) ) allocate (so  (ncells), stat = istat)
       endif
    endif
    !
    if ( .not.oned ) then
       !
       if ( istat == 0 .and. .not.allocated(dpv   ) ) allocate (dpv   (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(hv    ) ) allocate (hv    (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(hvm   ) ) allocate (hvm   (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(hvmn  ) ) allocate (hvmn  (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(hvmo  ) ) allocate (hvmo  (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(cfricv) ) allocate (cfricv(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(aporv ) ) allocate (aporv (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(bporv ) ) allocate (bporv (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(cporv ) ) allocate (cporv (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(windv ) ) allocate (windv (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(cwndv ) ) allocate (cwndv (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(flov  ) ) allocate (flov  (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(wetv  ) ) allocate (wetv  (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(presv ) ) allocate (presv (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(nporv ) ) allocate (nporv (mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(hporv ) ) allocate (hporv (mcgrd), stat = istat)
       !
       if ( istat == 0 .and. .not.allocated(cvegv ) ) allocate (cvegv (mcgrd,kmax,2), stat = istat)
       !
    else
       !
       if ( .not.allocated(dpv   ) ) allocate (dpv   (0))
       if ( .not.allocated(hv    ) ) allocate (hv    (0))
       if ( .not.allocated(hvm   ) ) allocate (hvm   (0))
       if ( .not.allocated(hvmn  ) ) allocate (hvmn  (0))
       if ( .not.allocated(hvmo  ) ) allocate (hvmo  (0))
       if ( .not.allocated(cfricv) ) allocate (cfricv(0))
       if ( .not.allocated(aporv ) ) allocate (aporv (0))
       if ( .not.allocated(bporv ) ) allocate (bporv (0))
       if ( .not.allocated(cporv ) ) allocate (cporv (0))
       if ( .not.allocated(windv ) ) allocate (windv (0))
       if ( .not.allocated(cwndv ) ) allocate (cwndv (0))
       if ( .not.allocated(flov  ) ) allocate (flov  (0))
       if ( .not.allocated(wetv  ) ) allocate (wetv  (0))
       if ( .not.allocated(presv ) ) allocate (presv (0))
       if ( .not.allocated(nporv ) ) allocate (nporv (0))
       if ( .not.allocated(hporv ) ) allocate (hporv (0))
       if ( .not.allocated(cvegv ) ) allocate (cvegv (0,0,0))
       !
    endif
    !
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: flow data (friction, wind, depths and wet/dry data) and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    dps    = 0.
    dpu    = 0.
    !
    hs     = 0.
    hso    = 0.
    hu     = 0.
    hum    = 0.
    humn   = 0.
    humo   = 0.
    !
    cfricu = 0.
    apors  = 0.
    bpors  = 0.
    cpors  = 0.
    aporu  = 0.
    bporu  = 0.
    cporu  = 0.
    windu  = 0.
    cwndu  = 0.
    !
    flos   = -9999.
    flou   = -9999.
    !
    lfbs   = 0
    !
    if ( svwp ) patm = 0.
    !
    wets   = 0
    wetu   = 0
    !
    presp  = 0
    presu  = 0
    !
    skc    = 0.
    !
    brks   = 0
    !
    npors  = 1.
    nporu  = 1.
    hpors  = ppor(2)
    hporu  = ppor(2)
    !
    cvegu  = 0.
    !
    hindun = 0.
    !
    if ( relwav ) smax = 0.
    !
    if ( .not.oned ) then
       dpv    = 0.
       hv     = 0.
       hvm    = 0.
       hvmn   = 0.
       hvmo   = 0.
       cfricv = 0.
       aporv  = 0.
       bporv  = 0.
       cporv  = 0.
       windv  = 0.
       cwndv  = 0.
       flov   = -9999.
       wetv   = 0
       presv  = 0
       nporv  = 1.
       hporv  = ppor(2)
       cvegv  = 0.
    endif
    !
    ! fill arrays dps (bottom values in wl-points) and dpu/dpv (tilded bottom values in velocity points)
    !
    if ( optg /= 5 ) then
       call SwashFlowDP
    else
       call SwashFlowUDP
    endif
    if (STPNOW()) return
    !
    ! determine some necessary parameters for porous structures
    !
    if ( optg /= 5 .and. iporos == 1 ) call SwashPorousStruc
    if (STPNOW()) return
    !
    ! determine some parameters for floating objects
    !
    if ( optg /= 5 .and. ifloat /= 0 ) call SwashFloatObjects
    if (STPNOW()) return
    !
    ! compute the Coriolis parameter, if appropriate
    !
    if ( coriolis ) then
       !
       if ( optg /= 5 ) then
          if ( .not.allocated(fcor) ) allocate(fcor(mcgrd,2   ))
          if ( .not.allocated(cfu ) ) allocate(cfu (mcgrd,kmax))
          if ( .not.allocated(cfv ) ) allocate(cfv (mcgrd,kmax))
       else
          if ( .not.allocated(fcor) ) allocate(fcor(ncells,1))
          if ( .not.allocated(cfu ) ) allocate(cfu (nfaces,kmax))
       endif
       fcor = 0.
       cfu  = 0.
       if ( optg /= 5 ) cfv  = 0.
       !
       if ( leds(23) /= 2 ) then
          !
          if ( pcor(1) /= -999. ) then
             !
             fcor = pcor(1)
             !
          else
             !
             if ( kspher == 0 ) then
                !
                ! Cartesian grid
                !
                fcor = sin(degrad*ylatc)
                !
             else
                !
                ! spherical grid
                !
                if ( optg /= 5 ) then
                   !
                   ! compute sine of latitude in u-points
                   !
                   do iy = 2, myc-1
                      do ix = 1, mxc-1
                         !
                         indx = kgrpnt(ix,iy)
                         !
                         if ( indx > 1 ) then
                            !
                            fcor(indx,1) = sin(degrad*(0.5*(ycgrid(ix,iy)+ycgrid(ix,iy-1)) + yoffs))
                            !
                         endif
                         !
                      enddo
                   enddo
                   !
                   call periodic ( fcor(1,1), kgrpnt, 1, 1 )
                   !
                   ! compute sine of latitude in v-points
                   !
                   do ix = 2, mxc-1
                      do iy = 1, myc-1
                         !
                         indx = kgrpnt(ix,iy)
                         !
                         if ( indx > 1 ) then
                            !
                            fcor(indx,2) = sin(degrad*(0.5*(ycgrid(ix,iy)+ycgrid(ix-1,iy)) + yoffs))
                            !
                         endif
                         !
                      enddo
                   enddo
                   !
                   call periodic ( fcor(1,2), kgrpnt, 1, 1 )
                   !
                else
                   !
                   ! compute sine of latitude in cells
                   !
                   do icell = 1, ncells
                      !
                      ! get y-coordinate of the cell center
                      !
                      yc = cell(icell)%attr(CELLCY)
                      !
                      fcor(icell,1) = sin(degrad*(yc + yoffs))
                      !
                   enddo
                   !
                endif
                !
             endif
             !
             ! compute Coriolis parameter
             !
             omega = pi2 / ( 24. * 3600. )
             !
             fcor  = 2. * omega * fcor
             !
          endif
          !
       else
          !
          if ( optg /= 5 ) then
             !
             ! interpolate Coriolis parameter from vertices to u-point
             !
             do iy = 2, myc-1
                do ix = 1, mxc-1
                   !
                   indx  = kgrpnt(ix,iy  )
                   indxu = kgrpnt(ix,iy-1)
                   if ( indxu == 1 ) indxu = indx
                   !
                   if ( indx > 1 ) fcor(indx,1) = 0.5 * ( corpf(indx) + corpf(indxu) )
                   !
                enddo
             enddo
             !
             call periodic ( fcor(1,1), kgrpnt, 1, 1 )
             !
             ! interpolate Coriolis parameter from vertices to v-point
             !
             do ix = 2, mxc-1
                do iy = 1, myc-1
                   !
                   indx  = kgrpnt(ix  ,iy)
                   indxr = kgrpnt(ix-1,iy)
                   if ( indxr == 1 ) indxr = indx
                   !
                   if ( indx > 1 ) fcor(indx,2) = 0.5 * ( corpf(indx) + corpf(indxr) )
                   !
                enddo
             enddo
             !
             call periodic ( fcor(1,2), kgrpnt, 1, 1 )
             !
          else
             !
             ! interpolate Coriolis parameter from vertices to cells
             !
             do icell = 1, ncells
                !
                ! get vertices of present cell
                !
                vc(1) = cell(icell)%atti(CELLV1)
                vc(2) = cell(icell)%atti(CELLV2)
                vc(3) = cell(icell)%atti(CELLV3)
                !
                fcor(icell,1) = ( corpf(vc(1)) + corpf(vc(2)) + corpf(vc(3)) )/ 3.
                !
             enddo
             !
          endif
          !
       endif
       !
    else
       if ( .not.allocated(fcor) ) allocate (fcor(0,0))
       if ( .not.allocated(cfu ) ) allocate (cfu (0,0))
       if ( .not.allocated(cfv ) ) allocate (cfv (0,0))
    endif
    !
    ! allocate and initialize other flow data (discharges and some help arrays)
    !
    if ( .not.allocated(qx) ) allocate(qx(mcgrd,kmax), stat = istat)
    if ( istat == 0 .and. .not.allocated(apomu ) ) allocate (apomu (mcgrd,kmax), stat = istat)
    if ( istat == 0 .and. .not.allocated(bpomu ) ) allocate (bpomu (mcgrd,kmax), stat = istat)
    if ( istat == 0 .and. .not.allocated(cpomu ) ) allocate (cpomu (mcgrd,kmax), stat = istat)
    if ( optg /= 5 ) then
       if ( istat == 0 .and. .not.allocated(wndimp) ) allocate (wndimp(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(sponxl) ) allocate (sponxl(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(sponxr) ) allocate (sponxr(mcgrd), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(wndimp) ) allocate (wndimp(nfaces     ), stat = istat)
       if ( istat == 0 .and. .not.allocated(sponu ) ) allocate (sponu (nfaces,nspl), stat = istat)
       if ( istat == 0 .and. .not.allocated(spons ) ) allocate (spons (ncells,nspl), stat = istat)
    endif
    !
    if ( .not.oned .and. optg /= 5 ) then
       !
       if ( istat == 0 .and. .not.allocated(qy    ) ) allocate (qy    (mcgrd,kmax), stat = istat)
       if ( istat == 0 .and. .not.allocated(apomv ) ) allocate (apomv (mcgrd,kmax), stat = istat)
       if ( istat == 0 .and. .not.allocated(bpomv ) ) allocate (bpomv (mcgrd,kmax), stat = istat)
       if ( istat == 0 .and. .not.allocated(cpomv ) ) allocate (cpomv (mcgrd,kmax), stat = istat)
       if ( istat == 0 .and. .not.allocated(sponyb) ) allocate (sponyb(mcgrd)     , stat = istat)
       if ( istat == 0 .and. .not.allocated(sponyt) ) allocate (sponyt(mcgrd)     , stat = istat)
       !
    else
       !
       if ( .not.allocated(qy    ) ) allocate (qy    (0,0))
       if ( .not.allocated(apomv ) ) allocate (apomv (0,0))
       if ( .not.allocated(bpomv ) ) allocate (bpomv (0,0))
       if ( .not.allocated(cpomv ) ) allocate (cpomv (0,0))
       if ( .not.allocated(sponyb) ) allocate (sponyb(0)  )
       if ( .not.allocated(sponyt) ) allocate (sponyt(0)  )
       !
    endif
    !
    if ( oned ) then
       if ( istat == 0 .and. .not.allocated(qm ) ) allocate (qm (mcgrd,kmax), stat = istat)
    else
       if ( istat == 0 .and. .not.allocated(qxm) ) allocate (qxm(mcgrd,kmax), stat = istat)
       if ( istat == 0 .and. .not.allocated(qym) ) allocate (qym(mcgrd,kmax), stat = istat)
    endif
    if ( istat == 0 .and. .not.allocated(ua) ) allocate (ua(mcgrd,kmax), stat = istat)
    if ( .not.oned ) then
       if ( istat == 0 .and. .not.allocated(ua2) ) allocate (ua2(mcgrd,kmax), stat = istat)
    endif
    if ( istat == 0 .and. .not.allocated(up) ) allocate (up(mcgrd,kmax), stat = istat)
    !
    if ( optg == 5 ) then
       if ( istat == 0 .and. .not.allocated(qn  ) ) allocate(qn  (nfaces,kmax  ), stat = istat)
       if ( istat == 0 .and. .not.allocated(quf ) ) allocate(quf (ncells,kmax,3), stat = istat)
       if ( istat == 0 .and. .not.allocated(uvc ) ) allocate(uvc (ncells,kmax,2), stat = istat)
       if ( istat == 0 .and. .not.allocated(divu) ) allocate(divu(ncells,kmax  ), stat = istat)
       if ( istat == 0 .and. .not.allocated(divq) ) allocate(divq(ncells,kmax  ), stat = istat)
    endif
    !
    if ( kmax == 1 ) then
       !
       if ( .not.mimetic ) then
          if ( oned ) then
             if ( istat == 0 .and. .not.allocated(advec ) ) allocate (advec (mcgrd), stat = istat)
          else if ( optg == 5 ) then
             if ( istat == 0 .and. .not.allocated(advec ) ) allocate (advec(nfaces), stat = istat)
          else
             if ( istat == 0 .and. .not.allocated(advecx) ) allocate (advecx(mcgrd), stat = istat)
             if ( istat == 0 .and. .not.allocated(advecy) ) allocate (advecy(mcgrd), stat = istat)
          endif
          if ( optg /= 5 ) then
             if ( istat == 0 .and. .not.allocated(visc) ) allocate (visc(mcgrd), stat = istat)
          else
             if ( istat == 0 .and. .not.allocated(visc) ) allocate (visc(nfaces), stat = istat)
          endif
       endif
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(cbot ) ) allocate (cbot (mcgrd), stat = istat)
          if ( istat == 0 .and. .not.allocated(pgrad) ) allocate (pgrad(mcgrd), stat = istat)
          if ( istat == 0 .and. .not.allocated(qgrad) ) allocate (qgrad(mcgrd), stat = istat)
          if ( istat == 0 .and. .not.allocated(corf ) ) allocate (corf (mcgrd), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(cbot ) ) allocate (cbot (nfaces), stat = istat)
          if ( istat == 0 .and. .not.allocated(pgrad) ) allocate (pgrad(nfaces), stat = istat)
          if ( istat == 0 .and. .not.allocated(qgrad) ) allocate (qgrad(nfaces), stat = istat)
          if ( istat == 0 .and. .not.allocated(corf ) ) allocate (corf (nfaces), stat = istat)
       endif
       !
    else
       !
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(wom) ) allocate (wom(mcgrd,0:kmax), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(wom) ) allocate (wom(ncells,0:kmax), stat = istat)
       endif
       !
    endif
    !
    if ( mtimei == 2 ) then
       !
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(ds) ) allocate (ds(mcgrd), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(ds) ) allocate (ds(ncells), stat = istat)
       endif
       !
       if ( istat == 0 .and. .not.allocated(ui) ) allocate (ui(mcgrd,kmax), stat = istat)
       if ( istat == 0 .and. .not.allocated(vi) ) allocate (vi(mcgrd,kmax), stat = istat)
       !
       if ( oned ) then
          if ( istat == 0 .and. .not.allocated(teta) ) allocate (teta(mcgrd), stat = istat)
          if ( lpproj ) then
             if ( istat == 0 .and. .not.allocated(dqgrd) ) allocate (dqgrd(mcgrd,kmax), stat = istat)
          endif
       else
          if ( istat == 0 .and. .not.allocated(tetau) ) allocate (tetau(mcgrd), stat = istat)
          if ( istat == 0 .and. .not.allocated(tetav) ) allocate (tetav(mcgrd), stat = istat)
          if ( lpproj ) then
             if ( istat == 0 .and. .not.allocated(dqgrdu) ) allocate (dqgrdu(mcgrd,kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(dqgrdv) ) allocate (dqgrdv(mcgrd,kmax), stat = istat)
          endif
       endif
       if ( istat == 0 .and. .not.allocated(teta2) ) allocate (teta2(mcgrd), stat = istat)
       !
    endif
    !
    if ( .not.oned .and. kmax > 1 ) then
       !
       if ( istat == 0 .and. .not.allocated(uwetp) ) allocate (uwetp(mcgrd), stat = istat)
       if ( istat == 0 .and. .not.allocated(vwetp) ) allocate (vwetp(mcgrd), stat = istat)
       !
    endif
    !
    if ( .not.oned ) then
       !
       if ( istat == 0 .and. .not.allocated(slimp) ) allocate (slimp(myc), stat = istat)
       if ( istat == 0 .and. .not.allocated(srimp) ) allocate (srimp(myc), stat = istat)
       if ( istat == 0 .and. .not.allocated(sbimp) ) allocate (sbimp(mxc), stat = istat)
       if ( istat == 0 .and. .not.allocated(stimp) ) allocate (stimp(mxc), stat = istat)
       !
    endif
    !
    if ( optg == 5 ) then
       !
       if ( istat == 0 .and. .not.allocated(wlimp) ) allocate (wlimp(nfaces), stat = istat)
       !
       wlimp = .false.
       !
    endif
    !
    if ( .not.oned .and. kmax > 1 ) then
       !
       if ( istat == 0 .and. .not.allocated(msta) ) allocate (msta(myc), stat = istat)
       if ( istat == 0 .and. .not.allocated(mend) ) allocate (mend(myc), stat = istat)
       if ( istat == 0 .and. .not.allocated(nsta) ) allocate (nsta(mxc), stat = istat)
       if ( istat == 0 .and. .not.allocated(nend) ) allocate (nend(mxc), stat = istat)
       !
    endif
    !
    if ( spwidl > 0. .or. spwidr > 0. ) then
       !
       if ( oned ) then
          !
          if ( istat == 0 .and. .not.allocated(tbndx) ) allocate (tbndx(  1,  kmax+1), stat = istat)
          !
       else
          !
          if ( istat == 0 .and. .not.allocated(tbndx) ) allocate (tbndx(myc,2*kmax+1), stat = istat)
          !
       endif
       !
    endif
    !
    if ( spwidb > 0. .or. spwidt > 0. ) then
       !
       if ( istat == 0 .and. .not.allocated(tbndy) ) allocate (tbndy(mxc,2*kmax+1), stat = istat)
       !
    endif
    !
    if ( iwvgen /= 0 ) then
       !
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(srcm) ) allocate (srcm(mcgrd), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(srcm) ) allocate (srcm(ncells), stat = istat)
          if ( istat == 0 .and. .not.allocated(dsrc) ) allocate (dsrc(ncells), stat = istat)
       endif
       !
    endif
    !
    if ( ihvisc /= 0 ) then
       !
       if ( optg /= 5 ) then
          if ( ihvisc > 1 ) then
             if ( istat == 0 .and. .not.allocated(vnu2d) ) allocate (vnu2d(mcgrd), stat = istat)
          endif
       else
          if ( istat == 0 .and. .not.allocated(vnu2d) ) allocate (vnu2d(nfaces), stat = istat)
          if ( ihvisc > 1 ) then
             if ( istat == 0 .and. .not.allocated(xshear) ) allocate (xshear(nverts), stat = istat)
          endif
       endif
       !
    endif
    !
    if ( idens /= 0 ) then
       !
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(rho) ) allocate (rho(mcgrd,kmax), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(rho) ) allocate (rho(ncells,kmax), stat = istat)
       endif
       !
    else
       !
       if ( .not.allocated(rho) ) allocate (rho(1,1))
       !
    endif
    !
    if ( ltrans > 0 ) then
       !
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(rp) ) allocate (rp(mcgrd,kmax,ltrans), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(rp) ) allocate (rp(ncells,kmax,ltrans), stat = istat)
       endif
       !
       if ( oned ) then
          !
          if ( istat == 0 .and. .not.allocated(flux) ) allocate (flux(mcgrd,kmax,1), stat = istat)
          !
          if ( .not.allocated(cbndl) ) allocate (cbndl(1,kmax,ltrans))
          if ( .not.allocated(cbndr) ) allocate (cbndr(1,kmax,ltrans))
          if ( .not.allocated(cbndb) ) allocate (cbndb(0,0,0))
          if ( .not.allocated(cbndt) ) allocate (cbndt(0,0,0))
          !
          if ( .not.allocated(coutl) ) allocate (coutl(1,kmax))
          if ( .not.allocated(coutr) ) allocate (coutr(1,kmax))
          if ( .not.allocated(coutb) ) allocate (coutb(0,0))
          if ( .not.allocated(coutt) ) allocate (coutt(0,0))
          !
          if ( .not.allocated(icretl) ) allocate (icretl(1,kmax))
          if ( .not.allocated(icretr) ) allocate (icretr(1,kmax))
          if ( .not.allocated(icretb) ) allocate (icretb(0,0))
          if ( .not.allocated(icrett) ) allocate (icrett(0,0))
          !
       else if ( optg == 5 ) then
          !
          if ( istat == 0 .and. .not.allocated(fluxt ) ) allocate (fluxt (nfaces ,kmax       ), stat = istat)
          if ( istat == 0 .and. .not.allocated(cbndu ) ) allocate (cbndu (nfacesb,kmax,ltrans), stat = istat)
          if ( istat == 0 .and. .not.allocated(coutu ) ) allocate (coutu (nfacesb,kmax       ), stat = istat)
          if ( istat == 0 .and. .not.allocated(icretu) ) allocate (icretu(nfacesb,kmax       ), stat = istat)
          !
       else
          !
          if ( istat == 0 .and. .not.allocated(flux) ) allocate (flux(mcgrd,kmax,2), stat = istat)
          !
          if ( istat == 0 .and. .not.allocated(cbndl) ) allocate (cbndl(myc,kmax,ltrans), stat = istat)
          if ( istat == 0 .and. .not.allocated(cbndr) ) allocate (cbndr(myc,kmax,ltrans), stat = istat)
          if ( istat == 0 .and. .not.allocated(cbndb) ) allocate (cbndb(mxc,kmax,ltrans), stat = istat)
          if ( istat == 0 .and. .not.allocated(cbndt) ) allocate (cbndt(mxc,kmax,ltrans), stat = istat)
          !
          if ( istat == 0 .and. .not.allocated(coutl) ) allocate (coutl(myc,kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(coutr) ) allocate (coutr(myc,kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(coutb) ) allocate (coutb(mxc,kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(coutt) ) allocate (coutt(mxc,kmax), stat = istat)
          !
          if ( istat == 0 .and. .not.allocated(icretl) ) allocate (icretl(myc,kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(icretr) ) allocate (icretr(myc,kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(icretb) ) allocate (icretb(mxc,kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(icrett) ) allocate (icrett(mxc,kmax), stat = istat)
          !
       endif
       !
    endif
    !
    if ( itrans /= 0 ) then
       !
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(rpo) ) allocate (rpo(mcgrd,kmax), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(rpo) ) allocate (rpo(ncells,kmax), stat = istat)
       endif
       !
    endif
    !
    ! local time stepping (unstructured-mesh transport module)
    !
    if ( itrans /= 0 .and. optg == 5 ) then
       !
       if (                  .not.allocated(mcell) ) allocate (mcell(ncells                  ), stat = istat)
       if ( istat == 0 .and. .not.allocated(mlev ) ) allocate (mlev (ncells                  ), stat = istat)
       if ( istat == 0 .and. .not.allocated(dtc  ) ) allocate (dtc  (ncells                  ), stat = istat)
       if ( istat == 0 .and. .not.allocated(rp0  ) ) allocate (rp0  (ncells,kmax             ), stat = istat)
       if ( istat == 0 .and. .not.allocated(rp1  ) ) allocate (rp1  (ncells,kmax             ), stat = istat)
       if ( istat == 0 .and. .not.allocated(rpi  ) ) allocate (rpi  (ncells,kmax,0:2**mtamx-1), stat = istat)
       !
    endif
    !
    if ( icreep /= 0 ) then
       !
       if ( istat == 0 .and. .not.allocated(wrk ) ) allocate (wrk (mcgrd,  kmax), stat = istat)
       if ( istat == 0 .and. .not.allocated(zkso) ) allocate (zkso(mcgrd,0:kmax), stat = istat)
       !
    endif
    !
    if ( kmax > 1 ) then
       !
       if ( optg /= 5 ) then
          if ( istat == 0 .and. .not.allocated(vnu3d) ) allocate (vnu3d(mcgrd,0:kmax), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(vnu3d) ) allocate (vnu3d(ncells,0:kmax), stat = istat)
       endif
       if ( ltur /= 0 ) then
          if ( optg /= 5 ) then
             if ( istat == 0 .and. .not.allocated(rtur) ) allocate (rtur(mcgrd,0:kmax,ltur), stat = istat)
          else
             if ( istat == 0 .and. .not.allocated(rtur) ) allocate (rtur(ncells,0:kmax,ltur), stat = istat)
          endif
       endif
       if ( iturb /= 0 ) then
          if ( optg /= 5 ) then
             if ( istat == 0 .and. .not.allocated(zshear) ) allocate (zshear(mcgrd,kmax-1), stat = istat)
          else
             if ( istat == 0 .and. .not.allocated(zshear) ) allocate (zshear(ncells,kmax-1), stat = istat)
          endif
       endif
       if ( irough == 4 ) then
          if ( optg /= 5 ) then
             if ( istat == 0 .and. .not.allocated(logfrc) ) allocate (logfrc(mcgrd,2), stat = istat)
          else
             if ( istat == 0 .and. .not.allocated(logfrc) ) allocate (logfrc(ncells,2), stat = istat)
          endif
       endif
       if ( iturb > 1 ) then
          if ( istat == 0 .and. .not.allocated(rsuu) ) allocate (rsuu(mcgrd,1:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(rsuw) ) allocate (rsuw(mcgrd,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(rswu) ) allocate (rswu(mcgrd,1:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(rsww) ) allocate (rsww(mcgrd,1:kmax), stat = istat)
          !
          if ( .not.oned ) then
             if ( istat == 0 .and. .not.allocated(rsuv) ) allocate (rsuv(mcgrd,1:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(rswv) ) allocate (rswv(mcgrd,1:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(rsvu) ) allocate (rsvu(mcgrd,1:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(rsvv) ) allocate (rsvv(mcgrd,1:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(rsvw) ) allocate (rsvw(mcgrd,0:kmax), stat = istat)
          else
             if ( .not.allocated(rsuv) ) allocate (rsuv(0,0))
             if ( .not.allocated(rswv) ) allocate (rswv(0,0))
             if ( .not.allocated(rsvu) ) allocate (rsvu(0,0))
             if ( .not.allocated(rsvv) ) allocate (rsvv(0,0))
             if ( .not.allocated(rsvw) ) allocate (rsvw(0,0))
          endif
          !
          if ( .not.oned ) then
             if ( istat == 0 .and. .not.allocated(logwll) ) allocate (logwll(myc,kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(logwlr) ) allocate (logwlr(myc,kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(logwlb) ) allocate (logwlb(mxc,kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(logwlt) ) allocate (logwlt(mxc,kmax), stat = istat)
          endif
       endif
       !
    endif
    !
    if ( iamb /= 0 ) then
       if ( oned ) then
          if ( istat == 0 .and. .not.allocated(Ubck) ) allocate (Ubck(mcgrd,kpmax), stat = istat)
       else if ( optg == 5 ) then
          if ( istat == 0 .and. .not.allocated(Ubck) ) allocate (Ubck(nfaces,kpmax), stat = istat)
       else
          if ( istat == 0 .and. .not.allocated(Ubck) ) allocate (Ubck(mcgrd,kpmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(Vbck) ) allocate (Vbck(mcgrd,kpmax), stat = istat)
       endif
       if ( istat == 0 .and. .not.allocated(so ) ) allocate (so (mcgrd      ), stat = istat)
       if ( istat == 0 .and. .not.allocated(wrk) ) allocate (wrk(mcgrd,kpmax), stat = istat)
    endif
    !
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: other flow data (discharges and some help arrays) and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    qx     = 0.
    wndimp = 0.
    apomu  = 0.
    bpomu  = 0.
    cpomu  = 0.
    if ( optg /= 5 ) then
       sponxl = 0.
       sponxr = 0.
    else
       sponu(:,:)%gamma = 0.
       sponu(:,:)%bface = 1
       spons(:,:)%gamma = 0.
       spons(:,:)%bcell = 1
    endif
    !
    if ( .not.oned ) then
       qy     = 0.
       apomv  = 0.
       bpomv  = 0.
       cpomv  = 0.
       sponyb = 0.
       sponyt = 0.
    endif
    !
    if ( oned ) then
       qm  = 0.
    else
       qxm = 0.
       qym = 0.
    endif
    ua = 0.
    if ( .not.oned ) ua2 = 0.
    up = 0.
    !
    if ( optg == 5 ) then
       qn   = 0.
       quf  = 0.
       uvc  = 0.
       divu = 0.
       divq = 0.
    endif
    !
    if ( kmax == 1 ) then
       !
       if ( .not.mimetic ) then
          if ( oned .or. optg == 5 ) then
             advec  = 0.
          else
             advecx = 0.
             advecy = 0.
          endif
          visc  = 0.
       endif
       cbot  = 0.
       pgrad = 0.
       qgrad = 0.
       corf  = 0.
       !
    else
       !
       wom = 0.
       !
    endif
    !
    if ( mtimei == 2 ) then
       ds = 0.
       if ( lpproj ) then
          if ( oned ) then
             dqgrd  = 0.
          else
             dqgrdu = 0.
             dqgrdv = 0.
          endif
       endif
    endif
    !
    if ( iwvgen /= 0 ) then
       srcm = 0.
       if ( optg == 5 ) dsrc = 0.
    endif
    !
    if ( ihvisc /= 0 ) then
       if ( optg /= 5 ) then
          if ( ihvisc > 1 ) vnu2d = 0.
       else
          vnu2d = 0.
          if ( ihvisc > 1 ) xshear = 0.
       endif
    endif
    !
    if ( idens /= 0 ) rho = 0.
    !
    if ( ltrans > 0 ) then
       rp = 0.
       if ( oned ) then
          flux   = 0.
          cbndl  = 0.
          cbndr  = 0.
          coutl  = 0.
          coutr  = 0.
          icretl = 0.
          icretr = 0.
       else if ( optg == 5 ) then
          fluxt  = 0.
          cbndu  = 0.
          coutu  = 0.
          icretu = 0.
       else
          flux   = 0.
          cbndl  = 0.
          cbndr  = 0.
          cbndb  = 0.
          cbndt  = 0.
          coutl  = 0.
          coutr  = 0.
          coutb  = 0.
          coutt  = 0.
          icretl = 0.
          icretr = 0.
          icretb = 0.
          icrett = 0.
       endif
    endif
    !
    if ( itrans /= 0 ) rpo = 0.
    !
    if ( itrans /= 0 .and. optg == 5 ) then
       !
       mcell = 0
       mlev  = 0
       dtc   = 0.
       rp0   = 0.
       rp1   = 0.
       rpi   = 0.
       !
    endif
    !
    if ( icreep /= 0 ) then
       wrk  = 0.
       zkso = 0.
    endif
    !
    if ( kmax > 1 ) then
       !if ( irough /= 0 .and. pbot(1) /= 0. ) then
       !   vnu3d = bvisc + kinvis
       !else
       !   vnu3d = 0.
       !endif
       if ( iturb > 1 ) then
          vnu3d = kinvis
       else
          vnu3d = bvisc + kinvis
       endif
       if ( ltur /= 0 ) then
          if ( iturb /= 0 ) then
             rtur(:,:,1) = tkeini
             rtur(:,:,2) = epsini
          endif
       endif
       if ( iturb /= 0 ) zshear = 0.
       if ( irough == 4 ) logfrc = 0.
       if ( iturb > 1 ) then
          rsuu = 0.
          rsuw = 0.
          rswu = 0.
          rsww = 0.
          if ( .not.oned ) then
             rsuv   = 0.
             rswv   = 0.
             rsvu   = 0.
             rsvv   = 0.
             rsvw   = 0.
             logwll = 0.
             logwlr = 0.
             logwlb = 0.
             logwlt = 0.
          endif
       endif
    endif
    !
    if ( iamb /= 0 ) then
       if ( oned .or. optg == 5 ) then
          Ubck = 0.
       else
          Ubck = 0.
          Vbck = 0.
       endif
       so  = 0.
       wrk = 0.
    endif
    !
    ! determine sponge layers
    !
    if ( optg /= 5 ) then
       !
       allocate(tempxl(MCGRDGL))
       allocate(tempxr(MCGRDGL))
       allocate(tempyb(MCGRDGL))
       allocate(tempyt(MCGRDGL))
       tempxl = 0.
       tempxr = 0.
       tempyb = 0.
       tempyt = 0.
       !
       call SwashSpongeLayer ( tempxl, tempxr, tempyb, tempyt )
       !
       ! create copies of parts of sponxl/r and sponyb/t for each subdomain
       !
       if ( oned ) then
          do ix = MXF, MXL
             indx  = KGRPGL(ix,1)
             indxo = kgrpnt(ix-MXF+1,1)
             sponxl(indxo) = tempxl(indx)
             sponxr(indxo) = tempxr(indx)
          enddo
       else
          do ix = MXF, MXL
             do iy = MYF, MYL
                indx  = KGRPGL(ix,iy)
                indxo = kgrpnt(ix-MXF+1,iy-MYF+1)
                sponxl(indxo) = tempxl(indx)
                sponxr(indxo) = tempxr(indx)
                sponyb(indxo) = tempyb(indx)
                sponyt(indxo) = tempyt(indx)
             enddo
          enddo
       endif
       !
       deallocate(tempxl,tempxr,tempyb,tempyt)
       !
    else
       call SwashUSpongLayer
    endif
    !
    ! determine ambient current (tide/wind surge)
    !
    if ( optg /= 5 .and. iamb /= 0 ) call SwashAmbCurrent
    if (STPNOW()) return
    !
    ! allocate and initialize layer interfaces and layer thicknesses
    !
    if ( kmax > 1 ) then
       !
       if ( optg /= 5 ) then
          if ( .not.allocated(zks) ) allocate (zks(mcgrd,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(zksnew) ) allocate (zksnew(mcgrd,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(zku   ) ) allocate (zku   (mcgrd,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(zkum  ) ) allocate (zkum  (mcgrd,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(zkuo  ) ) allocate (zkuo  (mcgrd,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hks   ) ) allocate (hks   (mcgrd,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hksnew) ) allocate (hksnew(mcgrd,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hkso  ) ) allocate (hkso  (mcgrd,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hku   ) ) allocate (hku   (mcgrd,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hkum  ) ) allocate (hkum  (mcgrd,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hkumn ) ) allocate (hkumn (mcgrd,  kmax), stat = istat)
          if ( .not.oned ) then
             if ( istat == 0 .and. .not.allocated(zkv  ) ) allocate (zkv  (mcgrd,0:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(zkvm ) ) allocate (zkvm (mcgrd,0:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(zkvo ) ) allocate (zkvo (mcgrd,0:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(hkv  ) ) allocate (hkv  (mcgrd,  kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(hkvm ) ) allocate (hkvm (mcgrd,  kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(hkvmn) ) allocate (hkvmn(mcgrd,  kmax), stat = istat)
          endif
       else
          if ( .not.allocated(zks) ) allocate (zks(ncells,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(zksnew) ) allocate (zksnew(ncells,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(zku   ) ) allocate (zku   (nfaces,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(zkum  ) ) allocate (zkum  (nfaces,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(zkuo  ) ) allocate (zkuo  (nfaces,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hks   ) ) allocate (hks   (ncells,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hksnew) ) allocate (hksnew(ncells,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hkso  ) ) allocate (hkso  (ncells,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hku   ) ) allocate (hku   (nfaces,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hkum  ) ) allocate (hkum  (nfaces,  kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hkumn ) ) allocate (hkumn (nfaces,  kmax), stat = istat)
       endif
       !
       if ( lsubg ) then
          if ( istat == 0 .and. .not.allocated(hksc ) ) allocate (hksc (mcgrd,kpmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hkuc ) ) allocate (hkuc (mcgrd,kpmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(hkumc) ) allocate (hkumc(mcgrd,kpmax), stat = istat)
          if ( .not.oned ) then
             if ( istat == 0 .and. .not.allocated(hkvc ) ) allocate (hkvc (mcgrd,kpmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(hkvmc) ) allocate (hkvmc(mcgrd,kpmax), stat = istat)
          endif
       endif
       !
    else
       !
       if ( .not.allocated(zks   ) ) allocate (zks   (0,0))
       if ( .not.allocated(zksnew) ) allocate (zksnew(0,0))
       if ( .not.allocated(zku   ) ) allocate (zku   (0,0))
       if ( .not.allocated(zkum  ) ) allocate (zkum  (0,0))
       if ( .not.allocated(zkuo  ) ) allocate (zkuo  (0,0))
       if ( .not.allocated(zkv   ) ) allocate (zkv   (0,0))
       if ( .not.allocated(zkvm  ) ) allocate (zkvm  (0,0))
       if ( .not.allocated(zkvo  ) ) allocate (zkvo  (0,0))
       if ( .not.allocated(hks   ) ) allocate (hks   (0,0))
       if ( .not.allocated(hksnew) ) allocate (hksnew(0,0))
       if ( .not.allocated(hkso  ) ) allocate (hkso  (0,0))
       if ( .not.allocated(hku   ) ) allocate (hku   (0,0))
       if ( .not.allocated(hkum  ) ) allocate (hkum  (0,0))
       if ( .not.allocated(hkumn ) ) allocate (hkumn (0,0))
       if ( .not.allocated(hkv   ) ) allocate (hkv   (0,0))
       if ( .not.allocated(hkvm  ) ) allocate (hkvm  (0,0))
       if ( .not.allocated(hkvmn ) ) allocate (hkvmn (0,0))
       if ( .not.allocated(hksc  ) ) allocate (hksc  (0,0))
       if ( .not.allocated(hkuc  ) ) allocate (hkuc  (0,0))
       if ( .not.allocated(hkumc ) ) allocate (hkumc (0,0))
       if ( .not.allocated(hkvc  ) ) allocate (hkvc  (0,0))
       if ( .not.allocated(hkvmc ) ) allocate (hkvmc (0,0))
       !
    endif
    !
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: layer interfaces and thicknesses and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    if ( kmax > 1 ) then
       !
       zks    = 0.
       zksnew = 0.
       zku    = 0.
       zkum   = 0.
       zkuo   = 0.
       hks    = 0.
       hksnew = 0.
       hkso   = 0.
       hku    = 0.
       hkum   = 0.
       hkumn  = 0.
       if ( .not.oned .and. optg /= 5 ) then
          zkv   = 0.
          zkvm  = 0.
          zkvo  = 0.
          hkv   = 0.
          hkvm  = 0.
          hkvmn = 0.
       endif
       !
       if ( lsubg ) then
          hksc  = 0.
          hkuc  = 0.
          hkumc = 0.
          if ( .not.oned ) then
             hkvc  = 0.
             hkvmc = 0.
          endif
       endif
       !
       ! fill bottom elevation
       !
       zks   (:,kmax) = -dps(:)
       zksnew(:,kmax) = -dps(:)
       !
       if ( optg /= 5 ) then
          !
          if ( oned ) then
             !
             do ix = 1, mxc-1
                !
                indx  = kgrpnt(ix  ,1)
                indxr = kgrpnt(ix+1,1)
                !
                zkum(indx,kmax) = -0.5 * ( dps(indx) + dps(indxr) )
                !
             enddo
             !
             if ( .not.depcds ) then
                !
                do ix = 1, mxc-1
                   !
                   indx  = kgrpnt(ix  ,1)
                   indxr = kgrpnt(ix+1,1)
                   !
                   zku(indx,kmax) = -min( dps(indx), dps(indxr) )
                   !
                enddo
                !
             else
                !
                zku(:,kmax) = zkum(:,kmax)
                !
             endif
             !
          else
             !
             do ix = 1, mxc-1
                do iy = 2, myc-1
                   !
                   indx  = kgrpnt(ix  ,iy)
                   indxr = kgrpnt(ix+1,iy)
                   if ( indxr == 1 ) indxr = indx
                   !
                   if ( indx > 1 ) zkum(indx,kmax) = -0.5 * ( dps(indx) + dps(indxr) )
                   !
                enddo
             enddo
             !
             ! synchronize layer interfaces at u-point at appropriate boundaries in case of repeating grid
             !
             call periodic ( zkum, kgrpnt, 0, kmax )
             !
             do iy = 1, myc-1
                do ix = 2, mxc-1
                   !
                   indx  = kgrpnt(ix,iy  )
                   indxu = kgrpnt(ix,iy+1)
                   if ( indxu == 1 ) indxu = indx
                   !
                   if ( indx > 1 ) zkvm(indx,kmax) = -0.5 * ( dps(indx) + dps(indxu) )
                   !
                enddo
             enddo
             !
             ! synchronize layer interfaces at v-point at appropriate boundaries in case of repeating grid
             !
             call periodic ( zkvm, kgrpnt, 0, kmax )
             !
             if ( .not.depcds ) then
                !
                do ix = 1, mxc-1
                   do iy = 2, myc-1
                      !
                      indx  = kgrpnt(ix  ,iy)
                      indxr = kgrpnt(ix+1,iy)
                      if ( indxr == 1 ) indxr = indx
                      !
                      if ( indx > 1 ) zku(indx,kmax) = -min( dps(indx), dps(indxr) )
                      !
                   enddo
                enddo
                !
                ! synchronize layer interfaces at u-point at appropriate boundaries in case of repeating grid
                !
                call periodic ( zku, kgrpnt, 0, kmax )
                !
                do iy = 1, myc-1
                   do ix = 2, mxc-1
                      !
                      indx  = kgrpnt(ix,iy  )
                      indxu = kgrpnt(ix,iy+1)
                      if ( indxu == 1 ) indxu = indx
                      !
                      if ( indx > 1 ) zkv(indx,kmax) = -min( dps(indx), dps(indxu) )
                      !
                   enddo
                enddo
                !
                ! synchronize layer interfaces at v-point at appropriate boundaries in case of repeating grid
                !
                call periodic ( zkv, kgrpnt, 0, kmax )
                !
             else
                !
                zku(:,kmax) = zkum(:,kmax)
                zkv(:,kmax) = zkvm(:,kmax)
                !
             endif
             !
          endif
          !
       else
          !
          do iface = 1, nfaces
             !
             fac = face(iface)%attr(FACELINPF)
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             if ( icelll /= 0 .and. icellr /= 0 ) then
                !
                zkum(iface,kmax) = - fac * dps(icelll) - (1.-fac) * dps(icellr)
                !
             else if ( icelll == 0 ) then
                !
                zkum(iface,kmax) = -dps(icellr)
                !
             else
                !
                zkum(iface,kmax) = -dps(icelll)
                !
             endif
             !
          enddo
          !
          if ( .not.depcds ) then
             !
             zku(:,kmax) = -dpu(:)
             !
          else
             !
             zku(:,kmax) = zkum(:,kmax)
             !
          endif
          !
       endif
       !
    endif
    !
    ! allocate and initialize arrays for boundary condition type
    !
    if ( oned ) then
       !
       if ( .not.allocated(ibl) ) allocate(ibl(1))
       if ( .not.allocated(ibr) ) allocate(ibr(1))
       if ( .not.allocated(ibb) ) allocate(ibb(0))
       if ( .not.allocated(ibt) ) allocate(ibt(0))
       !
       ibl = 1
       ibr = 1
       !
    else
       !
       if ( .not.allocated(ibl) ) allocate(ibl(myc))
       if ( .not.allocated(ibr) ) allocate(ibr(myc))
       if ( .not.allocated(ibb) ) allocate(ibb(mxc))
       if ( .not.allocated(ibt) ) allocate(ibt(mxc))
       !
       ibl = 1
       ibr = 1
       ibb = 1
       ibt = 1
       !
    endif
    !
    ! allocate and initialize arrays for storing water levels and constituents prescribed at boundary faces of triangular mesh
    !
    if ( optg == 5 ) then
       !
       if ( .not.allocated(bcs) ) allocate(bcs(nfaces), stat = istat)
       if ( istat == 0 .and. .not.allocated(bcso) ) allocate(bcso(nfaces), stat = istat)
       if ( istat == 0 .and. .not.allocated(pfac) ) allocate(pfac(nfaces), stat = istat)
       !
       if ( ltrans > 0 ) then
          !
          if ( istat == 0 .and. .not.allocated(bcrp ) ) allocate(bcrp (nfacesb,kmax,ltrans), stat = istat)
          if ( istat == 0 .and. .not.allocated(bcrpo) ) allocate(bcrpo(nfacesb,kmax,ltrans), stat = istat)
          !
       endif
       !
       if ( istat /= 0 ) then
          write (msgstr, '(a,i6)') 'allocation problem: variables for flexible mesh and return code is ',istat
          call msgerr ( 4, trim(msgstr) )
          return
       endif
       !
       bcs  = 0.
       bcso = 0.
       pfac = 0.
       !
       if ( ltrans > 0 ) then
          !
          bcrp  = 0.
          bcrpo = 0.
          !
       endif
       !
    endif
    !
    ! initialize transport constituents and store boundary values
    !
    if ( itrans /= 0 ) then
       if ( optg /= 5 ) then
          call SwashInitBCtrans
       else
          call SwashInitBCUtrans
       endif
    endif
    !
    ! allocate and initialize variables for non-hydrostatic computation
    !
    if ( ihydro /= 0 ) then
       !
       if ( kmax == 1 ) then
          !
          if ( optg /= 5 ) then
             !
             if ( .not.allocated(q) ) allocate(q(mcgrd,1), stat = istat)
             if ( istat == 0 .and. .not.allocated(dq   ) ) allocate (dq   (mcgrd,1  ), stat = istat)
             if ( istat == 0 .and. .not.allocated(gmatu) ) allocate (gmatu(mcgrd,1,2), stat = istat)
             if ( istat == 0 .and. .not.allocated(w1bot) ) allocate (w1bot(mcgrd    ), stat = istat)
             if ( istat == 0 .and. .not.allocated(w1top) ) allocate (w1top(mcgrd    ), stat = istat)
             !
             if ( .not.oned ) then
                !
                if ( istat == 0 .and. .not.allocated(gmatv) ) allocate (gmatv(mcgrd,1,2), stat = istat)
                !
             else
                !
                if ( .not.allocated(gmatv) ) allocate (gmatv(0,0,0))
                !
             endif
             !
          else
             !
             if ( .not.allocated(q) ) allocate(q(ncells,1), stat = istat)
             if ( istat == 0 .and. .not.allocated(dq   ) ) allocate (dq   (ncells,1  ), stat = istat)
             if ( istat == 0 .and. .not.allocated(gmatu) ) allocate (gmatu(nfaces,1,2), stat = istat)
             if ( istat == 0 .and. .not.allocated(w1bot) ) allocate (w1bot(ncells    ), stat = istat)
             if ( istat == 0 .and. .not.allocated(w1top) ) allocate (w1top(ncells    ), stat = istat)
             !
          endif
          !
       else
          !
          if ( optg /= 5 ) then
             if (.not.allocated(q)) allocate(q(mcgrd,kpmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(gmatu) ) allocate (gmatu(mcgrd,kpmax,6), stat = istat)
          else
             if (.not.allocated(q)) allocate(q(ncells,kpmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(gmatu) ) allocate (gmatu(nfaces,kpmax,6), stat = istat)
          endif
          !
          if ( optg /= 5 ) then
             if ( .not.oned ) then
                !
                if ( istat == 0 .and. .not.allocated(gmatv) ) allocate (gmatv(mcgrd,kpmax,6), stat = istat)
                !
             else
                !
                if ( .not.allocated(gmatv) ) allocate (gmatv(0,0,0))
                !
             endif
          endif
          !
          if ( lsubg ) then
             if ( istat == 0 .and. .not.allocated(qv    ) ) allocate (qv    (mcgrd,kmax  ), stat = istat)
             if ( istat == 0 .and. .not.allocated(gmatuv) ) allocate (gmatuv(mcgrd,kmax,6), stat = istat)
             if ( .not.oned ) then
                if ( istat == 0 .and. .not.allocated(gmatvv) ) allocate (gmatvv(mcgrd,kmax,6), stat = istat)
             else
                if ( .not.allocated(gmatvv) ) allocate (gmatvv(0,0,0))
             endif
          endif
          !
          if ( optg /= 5 ) then
             if ( ihydro == 1 ) then
                if ( istat == 0 .and. .not.allocated(gmatw) ) allocate (gmatw(mcgrd,kpmax,2*kpmax), stat = istat)
             else if ( ihydro == 2 ) then
                if ( istat == 0 .and. .not.allocated(gmatw) ) allocate (gmatw(mcgrd,kpmax,      2), stat = istat)
             endif
          else
             if ( ihydro == 1 ) then
                if ( istat == 0 .and. .not.allocated(gmatw) ) allocate (gmatw(ncells,kpmax,2*kpmax), stat = istat)
             else if ( ihydro == 2 ) then
                if ( istat == 0 .and. .not.allocated(gmatw) ) allocate (gmatw(ncells,kpmax,      2), stat = istat)
             endif
          endif
          !
          if ( ihydro == 1 .or. ihydro == 2 ) then
             !
             if ( optg /= 5 ) then
                if ( istat == 0 .and. .not.allocated(dq) ) allocate (dq(mcgrd,kpmax ), stat = istat)
                if ( istat == 0 .and. .not.allocated(w1) ) allocate (w1(mcgrd,0:kmax), stat = istat)
                if ( .not.oned ) then
                   if ( istat == 0 .and. .not.allocated(dmat) ) allocate (dmat(mcgrd,kpmax,12), stat = istat)
                else
                   if ( istat == 0 .and. .not.allocated(dmat) ) allocate (dmat(mcgrd,kpmax, 6), stat = istat)
                endif
             else
                if ( istat == 0 .and. .not.allocated(dq   ) ) allocate (dq   (ncells,kpmax    ), stat = istat)
                if ( istat == 0 .and. .not.allocated(w1   ) ) allocate (w1   (ncells,0:kmax   ), stat = istat)
                if ( istat == 0 .and. .not.allocated(dmatu) ) allocate (dmatu(ncells,kpmax,3,3), stat = istat)
             endif
             !
             if ( lsubg ) then
                if ( istat == 0 .and. .not.allocated(dqv ) ) allocate (dqv (mcgrd,kmax   ), stat = istat)
                if ( istat == 0 .and. .not.allocated(u1p ) ) allocate (u1p (mcgrd,1:kpmax), stat = istat)
                if ( istat == 0 .and. .not.allocated(w1p ) ) allocate (w1p (mcgrd,0:kpmax), stat = istat)
                if ( istat == 0 .and. .not.allocated(womp) ) allocate (womp(mcgrd,0:kpmax), stat = istat)
                if ( .not.oned ) then
                   if ( istat == 0 .and. .not.allocated(v1p) ) allocate (v1p(mcgrd,1:kpmax), stat = istat)
                else
                   if ( .not.allocated(v1p) ) allocate (v1p(0,0))
                endif
             endif
             !
             if ( istat == 0 .and. .not.allocated(apoks) ) allocate (apoks(mcgrd,0:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(bpoks) ) allocate (bpoks(mcgrd,0:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(cpoks) ) allocate (cpoks(mcgrd,0:kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(npoks) ) allocate (npoks(mcgrd,0:kmax), stat = istat)
             !
          else if ( ihydro == 3 ) then
             !
             if ( istat == 0 .and. .not.allocated(qbot  ) ) allocate (qbot  (mcgrd  ), stat = istat)
             if ( istat == 0 .and. .not.allocated(dq0   ) ) allocate (dq0   (mcgrd  ), stat = istat)
             if ( istat == 0 .and. .not.allocated(gmatu0) ) allocate (gmatu0(mcgrd,2), stat = istat)
             if ( istat == 0 .and. .not.allocated(w1bot ) ) allocate (w1bot (mcgrd  ), stat = istat)
             if ( istat == 0 .and. .not.allocated(w1top ) ) allocate (w1top (mcgrd  ), stat = istat)
             if ( .not.oned ) then
                if ( istat == 0 .and. .not.allocated(gmatv0) ) allocate (gmatv0(mcgrd,2), stat = istat)
             else
                if ( .not.allocated(gmatv0) ) allocate (gmatv0(0,0))
             endif
             !
          endif
          !
       endif
       !
       ! allocate arrays for storing variables at previous time level in case of nonstationary run
       !
       if ( nstatc == 1 ) then
          !
          if ( kmax == 1 .or. ihydro == 3 ) then
             !
             if ( optg /= 5 ) then
                if ( istat == 0 .and. .not.allocated(w0bot) ) allocate (w0bot(mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(w0top) ) allocate (w0top(mcgrd), stat = istat)
             else
                if ( istat == 0 .and. .not.allocated(w0bot) ) allocate (w0bot(ncells), stat = istat)
                if ( istat == 0 .and. .not.allocated(w0top) ) allocate (w0top(ncells), stat = istat)
             endif
             !
          else
             !
             if ( optg /= 5 ) then
                if ( istat == 0 .and. .not.allocated(w0) ) allocate (w0(mcgrd,0:kmax), stat = istat)
             else
                if ( istat == 0 .and. .not.allocated(w0) ) allocate (w0(ncells,0:kmax), stat = istat)
             endif
             if ( lsubg ) then
                if ( istat == 0 .and. .not.allocated(u0p) ) allocate (u0p(mcgrd,1:kpmax), stat = istat)
                if ( .not.oned ) then
                   if ( istat == 0 .and. .not.allocated(v0p) ) allocate (v0p(mcgrd,1:kpmax), stat = istat)
                endif
                if ( istat == 0 .and. .not.allocated(w0p) ) allocate (w0p(mcgrd,0:kpmax), stat = istat)
             endif
             !
          endif
          !
       else
          !
          if ( kmax == 1 .or. ihydro == 3 ) then
             if ( .not.allocated(w0bot) ) allocate (w0bot(0))
             if ( .not.allocated(w0top) ) allocate (w0top(0))
          else
             if ( .not.allocated(w0   ) ) allocate (w0 (0,0))
             if ( .not.allocated(u0p  ) ) allocate (u0p(0,0))
             if ( .not.allocated(v0p  ) ) allocate (v0p(0,0))
             if ( .not.allocated(w0p  ) ) allocate (w0p(0,0))
          endif
          !
       endif
       !
       if ( istat /= 0 ) then
          write (msgstr, '(a,i6)') 'allocation problem: variables for non-hydrostatic computation and return code is ',istat
          call msgerr ( 4, trim(msgstr) )
          return
       endif
       !
       q     = 0.
       gmatu = 0.
       if ( optg /= 5 ) then
          if ( .not.oned ) gmatv = 0.
       endif
       !
       if ( lsubg ) then
          qv     = 0.
          gmatuv = 0.
          if ( .not.oned ) gmatvv = 0.
       endif
       !
       if ( ihydro == 1 .or. ihydro == 2 ) then
          !
          dq = 0.
          !
       else if ( ihydro == 3 ) then
          !
          qbot   = 0.
          dq0    = 0.
          gmatu0 = 0.
          if ( .not.oned ) gmatv0 = 0.
          !
       endif
       !
       if ( kmax == 1 .or. ihydro == 3 ) then
          !
          w1bot = 0.
          w1top = 0.
          !
          if ( nstatc == 1 ) then
             !
             w0bot = 0.
             w0top = 0.
             !
          endif
          !
       else
          !
          gmatw = 0.
          if ( optg /= 5 ) then
             dmat  = 0.
          else
             dmatu = 0.
          endif
          w1    = 0.
          if ( nstatc == 1 ) w0 = 0.
          !
          if ( lsubg ) then
             dqv    = 0.
             w1p    = 0.
             womp   = 0.
             if ( nstatc == 1 ) then
                u0p = 0.
                if ( .not.oned ) v0p = 0.
                w0p = 0.
             endif
          endif
          !
          apoks = 0.
          bpoks = 0.
          cpoks = 0.
          npoks = 1.
          !
       endif
       !
    else
       !
       if ( .not.allocated(q     ) ) allocate (q     (0,0    ))
       if ( .not.allocated(dq    ) ) allocate (dq    (0,0    ))
       if ( .not.allocated(gmatu ) ) allocate (gmatu (0,0,0  ))
       if ( .not.allocated(gmatv ) ) allocate (gmatv (0,0,0  ))
       if ( .not.allocated(gmatw ) ) allocate (gmatw (0,0,0  ))
       if ( .not.allocated(dmat  ) ) allocate (dmat  (0,0,0  ))
       if ( .not.allocated(dmatu ) ) allocate (dmatu (0,0,0,0))
       if ( .not.allocated(u1p   ) ) allocate (u1p   (0,0    ))
       if ( .not.allocated(u0p   ) ) allocate (u0p   (0,0    ))
       if ( .not.allocated(v1p   ) ) allocate (v1p   (0,0    ))
       if ( .not.allocated(v0p   ) ) allocate (v0p   (0,0    ))
       if ( .not.allocated(w1bot ) ) allocate (w1bot (0      ))
       if ( .not.allocated(w1top ) ) allocate (w1top (0      ))
       if ( .not.allocated(w0bot ) ) allocate (w0bot (0      ))
       if ( .not.allocated(w0top ) ) allocate (w0top (0      ))
       if ( .not.allocated(w1    ) ) allocate (w1    (0,0    ))
       if ( .not.allocated(w1p   ) ) allocate (w1p   (0,0    ))
       if ( .not.allocated(w0    ) ) allocate (w0    (0,0    ))
       if ( .not.allocated(w0p   ) ) allocate (w0p   (0,0    ))
       if ( .not.allocated(womp  ) ) allocate (womp  (0,0    ))
       if ( .not.allocated(qbot  ) ) allocate (qbot  (0      ))
       if ( .not.allocated(dq0   ) ) allocate (dq0   (0      ))
       if ( .not.allocated(gmatu0) ) allocate (gmatu0(0,0    ))
       if ( .not.allocated(gmatv0) ) allocate (gmatv0(0,0    ))
       if ( .not.allocated(qv    ) ) allocate (qv    (0,0    ))
       if ( .not.allocated(dqv   ) ) allocate (dqv   (0,0    ))
       if ( .not.allocated(gmatuv) ) allocate (gmatuv(0,0,0  ))
       if ( .not.allocated(gmatvv) ) allocate (gmatvv(0,0,0  ))
       if ( .not.allocated(apoks ) ) allocate (apoks (0,0    ))
       if ( .not.allocated(bpoks ) ) allocate (bpoks (0,0    ))
       if ( .not.allocated(cpoks ) ) allocate (cpoks (0,0    ))
       if ( .not.allocated(npoks ) ) allocate (npoks (0,0    ))
       !
    endif
    !
    ! allocate and initialize variables for system of equations
    !
    if ( kmax == 1 ) then
       !
       if ( mtimei == 2 .or. ihydro == 1 .or. mimetic ) then
          !
          if ( oned ) then
             !
             if ( .not.allocated(a) ) allocate (a(mcgrd), stat = istat)
             if ( istat == 0 .and. .not.allocated(b) ) allocate (b(mcgrd), stat = istat)
             if ( istat == 0 .and. .not.allocated(c) ) allocate (c(mcgrd), stat = istat)
             if ( istat == 0 .and. .not.allocated(d) ) allocate (d(mcgrd), stat = istat)
             !
             if ( istat /= 0 ) then
                write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
                call msgerr ( 4, trim(msgstr) )
                return
             endif
             !
             a = 0.
             b = 0.
             c = 0.
             d = 0.
             !
             if ( NPROC > 2 ) then
                !
                if ( .not.allocated(z1) ) allocate (z1(mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(z2) ) allocate (z2(mcgrd), stat = istat)
                !
                if ( istat /= 0 ) then
                   write (msgstr, '(a,i6)') 'allocation problem: variables for DAC solver and return code is ',istat
                   call msgerr ( 4, trim(msgstr) )
                   return
                endif
                !
                z1 = 0.
                z2 = 0.
                !
                if ( inewt /= 0 ) then
                   !
                   if ( .not.allocated(z12) ) allocate (z12(mcgrd), stat = istat)
                   if ( istat == 0 .and. .not.allocated(z22) ) allocate (z22(mcgrd), stat = istat)
                   !
                   if ( istat /= 0 ) then
                      write (msgstr, '(a,i6)') 'allocation problem: variables for DAC solver and return code is ',istat
                      call msgerr ( 4, trim(msgstr) )
                      return
                   endif
                   !
                   z12 = 0d0
                   z22 = 0d0
                   !
                endif
                !
             endif
             !
             if ( inewt /= 0 ) then
                !
                if ( .not.allocated(an) ) allocate (an(mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(bn ) ) allocate (bn (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(cn ) ) allocate (cn (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(dn ) ) allocate (dn (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(lon) ) allocate (lon(mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(upn) ) allocate (upn(mcgrd), stat = istat)
                !
                if ( istat == 0 .and. .not.allocated(p0  ) ) allocate (p0  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(p1  ) ) allocate (p1  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(q0  ) ) allocate (q0  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(q1  ) ) allocate (q1  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(ba  ) ) allocate (ba  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(da  ) ) allocate (da  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(sola) ) allocate (sola(mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(vt  ) ) allocate (vt  (mcgrd), stat = istat)
                !
                if ( istat /= 0 ) then
                   write (msgstr, '(a,i6)') 'allocation problem: variables for Newton iteration and return code is ',istat
                   call msgerr ( 4, trim(msgstr) )
                   return
                endif
                !
                an  = 0d0
                bn  = 0d0
                cn  = 0d0
                dn  = 0d0
                lon = 0d0
                upn = 0d0
                !
                p0  = 1d0
                p1  = 1d0
                q0  = 0d0
                q1  = 0d0
                ba  = 1d0
                da  = 0d0
                sola= 0d0
                vt  = 0d0
                !
             endif
             !
          else
             !
             if ( optg /= 5 ) then
                if ( .not.allocated(amat) ) allocate (amat(mcgrd,5), stat = istat)
                if ( istat == 0 .and. .not.allocated(rhs ) ) allocate (rhs (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(resd) ) allocate (resd(mcgrd), stat = istat)
             else
                if ( .not.allocated(amat) ) allocate (amat(ncells,0:3), stat = istat)
                if ( istat == 0 .and. .not.allocated(rhs ) ) allocate (rhs (ncells), stat = istat)
                if ( istat == 0 .and. .not.allocated(resd) ) allocate (resd(ncells), stat = istat)
             endif
             !
             if ( istat /= 0 ) then
                write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
                call msgerr ( 4, trim(msgstr) )
                return
             endif
             !
             amat = 0.
             rhs  = 0.
             resd = 0.
             !
             if ( mtimei == 2 ) then
                !
                if ( optg /= 5 ) then
                   !
                   if ( inewt == 0 ) then
                      !
                      if ( .not.allocated(diag) ) allocate (diag(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(s) ) allocate (s(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(r) ) allocate (r(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(z) ) allocate (z(mcgrd), stat = istat)
                      !
                      if ( istat /= 0 ) then
                         write (msgstr, '(a,i6)') 'allocation problem: variables for CG solver and return code is ',istat
                         call msgerr ( 4, trim(msgstr) )
                         return
                      endif
                      !
                      diag = 1.
                      s    = 0.
                      r    = 0.
                      z    = 0.
                      !
                   else
                      !
                      if ( .not.allocated(diag2) ) allocate (diag2(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(resd2) ) allocate (resd2(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(s2   ) ) allocate (s2   (mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(r2   ) ) allocate (r2   (mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(z3   ) ) allocate (z3   (mcgrd), stat = istat)
                      !
                      if ( istat /= 0 ) then
                         write (msgstr, '(a,i6)') 'allocation problem: variables for CG solver and return code is ',istat
                         call msgerr ( 4, trim(msgstr) )
                         return
                      endif
                      !
                      diag2 = 1d0
                      resd2 = 0d0
                      s2    = 0d0
                      r2    = 0d0
                      z3    = 0d0
                      !
                   endif
                   !
                else
                   !
                   ! integer arrays to be used for the CSR format
                   !
                   if ( .not.allocated(ias) ) allocate (ias(ncells+1), stat = istat)
                   if ( istat == 0 .and. .not.allocated(irws) ) allocate (irws(  ncells,1:4), stat = istat)
                   if ( istat == 0 .and. .not.allocated(jas ) ) allocate (jas (4*ncells    ), stat = istat)
                   !
                   if ( istat == 0 .and. .not.allocated(dis ) ) allocate (dis (ncells), stat = istat)
                   if ( istat == 0 .and. .not.allocated(iplu) ) allocate (iplu(ncells), stat = istat)
                   !
                   if ( inewt == 0 ) then
                      !
                      if ( istat == 0 .and. .not.allocated(prec2) ) allocate (prec2(0:4*ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(s    ) ) allocate (s    (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(r    ) ) allocate (r    (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(z    ) ) allocate (z    (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(axs2 ) ) allocate (axs2 (  4*ncells), stat = istat)
                      !
                      if ( istat /= 0 ) then
                         write (msgstr, '(a,i6)') 'allocation problem: variables for CG solver and return code is ',istat
                         call msgerr ( 4, trim(msgstr) )
                         return
                      endif
                      !
                      prec2 = 0.
                      s     = 0.
                      r     = 0.
                      z     = 0.
                      axs2  = 0.
                      !
                   else
                      !
                      if ( istat == 0 .and. .not.allocated(prec3) ) allocate (prec3(0:4*ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(resd2) ) allocate (resd2(    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(s2   ) ) allocate (s2   (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(r2   ) ) allocate (r2   (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(z3   ) ) allocate (z3   (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(axs  ) ) allocate (axs  (  4*ncells), stat = istat)
                      !
                      if ( istat /= 0 ) then
                         write (msgstr, '(a,i6)') 'allocation problem: variables for CG solver and return code is ',istat
                         call msgerr ( 4, trim(msgstr) )
                         return
                      endif
                      !
                      prec3 = 0d0
                      resd2 = 0d0
                      s2    = 0d0
                      r2    = 0d0
                      z3    = 0d0
                      axs   = 0d0
                      !
                   endif
                   !
                   ! fill CSR arrays
                   !
                   call csrf ( ias, jas, irws, dis, 1, 3 )
                   !
                endif
                !
             endif
             !
             if ( inewt /= 0 ) then
                !
                if ( optg /= 5 ) then
                   if ( .not.allocated(amatn) ) allocate (amatn(mcgrd,5), stat = istat)
                   if ( istat == 0 .and. .not.allocated(rhsn ) ) allocate (rhsn (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(lon  ) ) allocate (lon  (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(upn  ) ) allocate (upn  (mcgrd  ), stat = istat)
                   !
                   if ( istat == 0 .and. .not.allocated(p0   ) ) allocate (p0   (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(p1   ) ) allocate (p1   (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(q0   ) ) allocate (q0   (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(q1   ) ) allocate (q1   (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(sola ) ) allocate (sola (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(amata) ) allocate (amata(mcgrd,5), stat = istat)
                   if ( istat == 0 .and. .not.allocated(rhsa ) ) allocate (rhsa (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(vt   ) ) allocate (vt   (mcgrd  ), stat = istat)
                else
                   if ( .not.allocated(amatn) ) allocate (amatn(ncells,0:3), stat = istat)
                   if ( istat == 0 .and. .not.allocated(rhsn ) ) allocate (rhsn (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(lon  ) ) allocate (lon  (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(upn  ) ) allocate (upn  (ncells    ), stat = istat)
                   !
                   if ( istat == 0 .and. .not.allocated(sol0 ) ) allocate (sol0 (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(sol1 ) ) allocate (sol1 (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(amata) ) allocate (amata(ncells,0:3), stat = istat)
                   if ( istat == 0 .and. .not.allocated(rhsa ) ) allocate (rhsa (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(vt   ) ) allocate (vt   (ncells    ), stat = istat)
                endif
                !
                if ( istat /= 0 ) then
                   write (msgstr, '(a,i6)') 'allocation problem: variables for Newton iteration and return code is ',istat
                   call msgerr ( 4, trim(msgstr) )
                   return
                endif
                !
                amatn = 0d0
                rhsn  = 0d0
                lon   = 0d0
                upn   = 0d0
                !
                if ( optg /= 5 ) then
                   p0   = 1d0
                   p1   = 1d0
                   q0   = 0d0
                   q1   = 0d0
                   sola = 0d0
                else
                   sol0 = 0d0
                   sol1 = 0d0
                endif
                !
                amata = 0d0
                rhsa  = 0d0
                vt    = 0d0
                !
             endif
             !
             if ( optg /= 5 ) then
                !
                if ( ihydro == 1 .or. mimetic ) then
                   !
                   if ( .not.allocated(cmat) ) allocate (cmat(mcgrd,5), stat = istat)
                   !
                   if ( istat /= 0 ) then
                      write (msgstr, '(a,i6)') 'allocation problem: variables for SIP solver and return code is ',istat
                      call msgerr ( 4, trim(msgstr) )
                      return
                   endif
                   !
                   cmat = 0.
                   !
                endif
                !
             else
                !
                if ( ihydro == 1 ) then
                   !
                   ! integer arrays to be used for the CSR format
                   !
                   if ( .not.allocated(ia) ) allocate (ia(ncells+1), stat = istat)
                   if ( istat == 0 .and. .not.allocated(irw) ) allocate (irw(  ncells,1:4), stat = istat)
                   if ( istat == 0 .and. .not.allocated(ja ) ) allocate (ja (4*ncells    ), stat = istat)
                   !
                   if ( istat == 0 .and. .not.allocated(di  ) ) allocate (di  (ncells), stat = istat)
                   if ( istat == 0 .and. .not.allocated(iplu) ) allocate (iplu(ncells), stat = istat)
                   !
                   if ( icond == 2 ) then
                      if ( istat == 0 .and. .not.allocated(diag) ) allocate (diag(ncells), stat = istat)
                   else if ( icond == 3 ) then
                      if ( istat == 0 .and. .not.allocated(prec2) ) allocate (prec2(0:4*ncells), stat = istat)
                   endif
                   if ( istat == 0 .and. .not.allocated(res ) ) allocate (res (ncells,1), stat = istat)
                   if ( istat == 0 .and. .not.allocated(res0) ) allocate (res0(ncells,1), stat = istat)
                   if ( istat == 0 .and. .not.allocated(p   ) ) allocate (p   (ncells,1), stat = istat)
                   if ( istat == 0 .and. .not.allocated(v   ) ) allocate (v   (ncells,1), stat = istat)
                   if ( istat == 0 .and. .not.allocated(u   ) ) allocate (u   (ncells,1), stat = istat)
                   if ( istat == 0 .and. .not.allocated(t   ) ) allocate (t   (ncells,1), stat = istat)
                   if ( istat == 0 .and. .not.allocated(ax2 ) ) allocate (ax2 (4*ncells), stat = istat)
                   if ( istat == 0 .and. .not.allocated(w2  ) ) allocate (w2  (ncells  ), stat = istat)
                   !
                   if ( icond == 2 ) then
                      diag = 1.
                   else if ( icond == 3 ) then
                      prec2 = 0.
                   endif
                   res  = 0.
                   res0 = 0.
                   p    = 0.
                   v    = 0.
                   u    = 0.
                   t    = 0.
                   ax2  = 0.
                   w2   = 0.
                   !
                   if ( istat /= 0 ) then
                      write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
                      call msgerr ( 4, trim(msgstr) )
                      return
                   endif
                   !
                   ! fill CSR arrays
                   !
                   call csrf ( ia, ja, irw, di, 1, 3 )
                   !
                endif
                !
             endif
             !
          endif
          !
       else
          !
          if ( .not.allocated(a) ) allocate (a(0))
          if ( .not.allocated(b) ) allocate (b(0))
          if ( .not.allocated(c) ) allocate (c(0))
          if ( .not.allocated(d) ) allocate (d(0))
          !
          if ( .not.allocated(z1) ) allocate (z1(0))
          if ( .not.allocated(z2) ) allocate (z2(0))
          !
          if ( .not.allocated(z12) ) allocate (z12(0))
          if ( .not.allocated(z22) ) allocate (z22(0))
          !
          if ( .not.allocated(an ) ) allocate (an (0))
          if ( .not.allocated(bn ) ) allocate (bn (0))
          if ( .not.allocated(cn ) ) allocate (cn (0))
          if ( .not.allocated(dn ) ) allocate (dn (0))
          if ( .not.allocated(lon) ) allocate (lon(0))
          if ( .not.allocated(upn) ) allocate (upn(0))
          !
          if ( .not.allocated(p0  ) ) allocate (p0  (0))
          if ( .not.allocated(p1  ) ) allocate (p1  (0))
          if ( .not.allocated(q0  ) ) allocate (q0  (0))
          if ( .not.allocated(q1  ) ) allocate (q1  (0))
          if ( .not.allocated(ba  ) ) allocate (ba  (0))
          if ( .not.allocated(da  ) ) allocate (da  (0))
          if ( .not.allocated(sola) ) allocate (sola(0))
          if ( .not.allocated(vt  ) ) allocate (vt  (0))
          !
          if ( .not.allocated(amat) ) allocate (amat(0,0))
          if ( .not.allocated(rhs ) ) allocate (rhs (0  ))
          if ( .not.allocated(resd) ) allocate (resd(0  ))
          !
          if ( .not.allocated(diag) ) allocate (diag(0))
          if ( .not.allocated(s   ) ) allocate (s   (0))
          if ( .not.allocated(r   ) ) allocate (r   (0))
          if ( .not.allocated(z   ) ) allocate (z   (0))
          !
          if ( .not.allocated(diag2) ) allocate (diag2(0))
          if ( .not.allocated(resd2) ) allocate (resd2(0))
          if ( .not.allocated(s2   ) ) allocate (s2   (0))
          if ( .not.allocated(r2   ) ) allocate (r2   (0))
          if ( .not.allocated(z3   ) ) allocate (z3   (0))
          !
          if ( .not.allocated(prec2) ) allocate (prec2(0))
          if ( .not.allocated(prec3) ) allocate (prec3(0))
          !
          if ( .not.allocated(ias ) ) allocate (ias (0))
          if ( .not.allocated(irws) ) allocate (irws(0,0))
          if ( .not.allocated(dis ) ) allocate (dis (0))
          if ( .not.allocated(jas ) ) allocate (jas (0))
          if ( .not.allocated(axs ) ) allocate (axs (0))
          if ( .not.allocated(axs2) ) allocate (axs2(0))
          !
          if ( .not.allocated(amatn) ) allocate (amatn(0,0))
          if ( .not.allocated(rhsn ) ) allocate (rhsn (0  ))
          !
          if ( .not.allocated(amata) ) allocate (amata(0,0))
          if ( .not.allocated(rhsa ) ) allocate (rhsa (0  ))
          if ( .not.allocated(sol0 ) ) allocate (sol0 (0  ))
          if ( .not.allocated(sol1 ) ) allocate (sol1 (0  ))
          !
          if ( .not.allocated(cmat) ) allocate (cmat(0,0))
          !
          if ( .not.allocated(res ) ) allocate (res (0,0))
          if ( .not.allocated(res0) ) allocate (res0(0,0))
          if ( .not.allocated(p   ) ) allocate (p   (0,0))
          if ( .not.allocated(v   ) ) allocate (v   (0,0))
          if ( .not.allocated(w   ) ) allocate (w   (0,0))
          if ( .not.allocated(u   ) ) allocate (u   (0,0))
          if ( .not.allocated(t   ) ) allocate (t   (0,0))
          !
          if ( .not.allocated(iplu) ) allocate (iplu(0))
          !
          if ( .not.allocated(ia ) ) allocate (ia (0))
          if ( .not.allocated(irw) ) allocate (irw(0,0))
          if ( .not.allocated(di ) ) allocate (di (0))
          if ( .not.allocated(ja ) ) allocate (ja (0))
          if ( .not.allocated(ax2) ) allocate (ax2(0))
          if ( .not.allocated(w2 ) ) allocate (w2 (0))
          !
       endif
       !
    else
       !
       if ( optg /= 5 ) then
          if ( .not.allocated(amatu) ) allocate (amatu(mcgrd,kmax,3), stat = istat)
          if ( istat == 0 .and. .not.allocated(rhsu) ) allocate (rhsu(mcgrd,kmax), stat = istat)
       else
          if ( .not.allocated(amatu) ) allocate (amatu(nfaces,kmax,3), stat = istat)
          if ( istat == 0 .and. .not.allocated(rhsu) ) allocate (rhsu(nfaces,kmax), stat = istat)
       endif
       !
       if ( istat /= 0 ) then
          write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
          call msgerr ( 4, trim(msgstr) )
          return
       endif
       !
       amatu = 0.
       rhsu  = 0.
       !
       if ( mtimei == 2 .or. ihydro == 3 .or. qmax == 1 ) then
          !
          if ( oned ) then
             !
             if ( .not.allocated(a) ) allocate (a(mcgrd), stat = istat)
             if ( istat == 0 .and. .not.allocated(b) ) allocate (b(mcgrd), stat = istat)
             if ( istat == 0 .and. .not.allocated(c) ) allocate (c(mcgrd), stat = istat)
             if ( istat == 0 .and. .not.allocated(d) ) allocate (d(mcgrd), stat = istat)
             !
             if ( istat /= 0 ) then
                write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
                call msgerr ( 4, trim(msgstr) )
                return
             endif
             !
             a = 0.
             b = 0.
             c = 0.
             d = 0.
             !
             if ( NPROC > 2 ) then
                !
                if ( .not.allocated(z1) ) allocate (z1(mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(z2) ) allocate (z2(mcgrd), stat = istat)
                !
                if ( istat /= 0 ) then
                   write (msgstr, '(a,i6)') 'allocation problem: variables for DAC solver and return code is ',istat
                   call msgerr ( 4, trim(msgstr) )
                   return
                endif
                !
                z1 = 0.
                z2 = 0.
                !
                if ( inewt /= 0 ) then
                   !
                   if ( .not.allocated(z12) ) allocate (z12(mcgrd), stat = istat)
                   if ( istat == 0 .and. .not.allocated(z22) ) allocate (z22(mcgrd), stat = istat)
                   !
                   if ( istat /= 0 ) then
                      write (msgstr, '(a,i6)') 'allocation problem: variables for DAC solver and return code is ',istat
                      call msgerr ( 4, trim(msgstr) )
                      return
                   endif
                   !
                   z12 = 0d0
                   z22 = 0d0
                   !
                endif
                !
             endif
             !
             if ( inewt /= 0 ) then
                !
                if ( .not.allocated(an) ) allocate (an(mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(bn ) ) allocate (bn (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(cn ) ) allocate (cn (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(dn ) ) allocate (dn (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(lon) ) allocate (lon(mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(upn) ) allocate (upn(mcgrd), stat = istat)
                !
                if ( istat == 0 .and. .not.allocated(p0  ) ) allocate (p0  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(p1  ) ) allocate (p1  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(q0  ) ) allocate (q0  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(q1  ) ) allocate (q1  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(ba  ) ) allocate (ba  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(da  ) ) allocate (da  (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(sola) ) allocate (sola(mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(vt  ) ) allocate (vt  (mcgrd), stat = istat)
                !
                if ( istat /= 0 ) then
                   write (msgstr, '(a,i6)') 'allocation problem: variables for Newton iteration and return code is ',istat
                   call msgerr ( 4, trim(msgstr) )
                   return
                endif
                !
                an  = 0d0
                bn  = 0d0
                cn  = 0d0
                dn  = 0d0
                lon = 0d0
                upn = 0d0
                !
                p0  = 1d0
                p1  = 1d0
                q0  = 0d0
                q1  = 0d0
                ba  = 1d0
                da  = 0d0
                sola= 0d0
                vt  = 0d0
                !
             endif
             !
          else
             !
             if ( optg /= 5 ) then
                if ( .not.allocated(amat) ) allocate (amat(mcgrd,5), stat = istat)
                if ( istat == 0 .and. .not.allocated(rhs ) ) allocate (rhs (mcgrd), stat = istat)
                if ( istat == 0 .and. .not.allocated(resd) ) allocate (resd(mcgrd), stat = istat)
             else
                if ( .not.allocated(amat) ) allocate (amat(ncells,0:3), stat = istat)
                if ( istat == 0 .and. .not.allocated(rhs ) ) allocate (rhs (ncells), stat = istat)
                if ( istat == 0 .and. .not.allocated(resd) ) allocate (resd(ncells), stat = istat)
             endif
             !
             if ( istat /= 0 ) then
                write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
                call msgerr ( 4, trim(msgstr) )
                return
             endif
             !
             amat = 0.
             rhs  = 0.
             resd = 0.
             !
             if ( mtimei == 2 ) then
                !
                if ( optg /= 5 ) then
                   !
                   if ( inewt == 0 ) then
                      !
                      if ( .not.allocated(diag) ) allocate (diag(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(s) ) allocate (s(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(r) ) allocate (r(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(z) ) allocate (z(mcgrd), stat = istat)
                      !
                      if ( istat /= 0 ) then
                         write (msgstr, '(a,i6)') 'allocation problem: variables for CG solver and return code is ',istat
                         call msgerr ( 4, trim(msgstr) )
                         return
                      endif
                      !
                      diag = 1.
                      s    = 0.
                      r    = 0.
                      z    = 0.
                      !
                   else
                      !
                      if ( .not.allocated(diag2) ) allocate (diag2(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(resd2) ) allocate (resd2(mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(s2   ) ) allocate (s2   (mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(r2   ) ) allocate (r2   (mcgrd), stat = istat)
                      if ( istat == 0 .and. .not.allocated(z3   ) ) allocate (z3   (mcgrd), stat = istat)
                      !
                      if ( istat /= 0 ) then
                         write (msgstr, '(a,i6)') 'allocation problem: variables for CG solver and return code is ',istat
                         call msgerr ( 4, trim(msgstr) )
                         return
                      endif
                      !
                      diag2 = 1d0
                      resd2 = 0d0
                      s2    = 0d0
                      r2    = 0d0
                      z3    = 0d0
                      !
                   endif
                   !
                else
                   !
                   ! integer arrays to be used for the CSR format
                   !
                   if ( .not.allocated(ias) ) allocate (ias(ncells+1), stat = istat)
                   if ( istat == 0 .and. .not.allocated(irws) ) allocate (irws(  ncells,1:4), stat = istat)
                   if ( istat == 0 .and. .not.allocated(jas ) ) allocate (jas (4*ncells    ), stat = istat)
                   !
                   if ( istat == 0 .and. .not.allocated(dis ) ) allocate (dis (ncells), stat = istat)
                   if ( istat == 0 .and. .not.allocated(iplu) ) allocate (iplu(ncells), stat = istat)
                   !
                   if ( inewt == 0 ) then
                      !
                      if ( istat == 0 .and. .not.allocated(prec2) ) allocate (prec2(0:4*ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(s    ) ) allocate (s    (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(r    ) ) allocate (r    (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(z    ) ) allocate (z    (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(axs2 ) ) allocate (axs2 (  4*ncells), stat = istat)
                      !
                      if ( istat /= 0 ) then
                         write (msgstr, '(a,i6)') 'allocation problem: variables for CG solver and return code is ',istat
                         call msgerr ( 4, trim(msgstr) )
                         return
                      endif
                      !
                      prec2 = 0.
                      s     = 0.
                      r     = 0.
                      z     = 0.
                      axs2  = 0.
                      !
                   else
                      !
                      if ( istat == 0 .and. .not.allocated(prec3) ) allocate (prec3(0:4*ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(resd2) ) allocate (resd2(    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(s2   ) ) allocate (s2   (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(r2   ) ) allocate (r2   (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(z3   ) ) allocate (z3   (    ncells), stat = istat)
                      if ( istat == 0 .and. .not.allocated(axs  ) ) allocate (axs  (  4*ncells), stat = istat)
                      !
                      if ( istat /= 0 ) then
                         write (msgstr, '(a,i6)') 'allocation problem: variables for CG solver and return code is ',istat
                         call msgerr ( 4, trim(msgstr) )
                         return
                      endif
                      !
                      prec3 = 0d0
                      resd2 = 0d0
                      s2    = 0d0
                      r2    = 0d0
                      z3    = 0d0
                      axs   = 0d0
                      !
                   endif
                   !
                   ! fill CSR arrays
                   !
                   call csrf ( ias, jas, irws, dis, 1, 3 )
                   !
                endif
                !
             endif
             !
             if ( inewt /= 0 ) then
                !
                if ( optg /= 5 ) then
                   if ( .not.allocated(amatn) ) allocate (amatn(mcgrd,5), stat = istat)
                   if ( istat == 0 .and. .not.allocated(rhsn ) ) allocate (rhsn (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(lon  ) ) allocate (lon  (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(upn  ) ) allocate (upn  (mcgrd  ), stat = istat)
                   !
                   if ( istat == 0 .and. .not.allocated(p0   ) ) allocate (p0   (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(p1   ) ) allocate (p1   (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(q0   ) ) allocate (q0   (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(q1   ) ) allocate (q1   (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(sola ) ) allocate (sola (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(amata) ) allocate (amata(mcgrd,5), stat = istat)
                   if ( istat == 0 .and. .not.allocated(rhsa ) ) allocate (rhsa (mcgrd  ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(vt   ) ) allocate (vt   (mcgrd  ), stat = istat)
                else
                   if ( .not.allocated(amatn) ) allocate (amatn(ncells,0:3), stat = istat)
                   if ( istat == 0 .and. .not.allocated(rhsn ) ) allocate (rhsn (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(lon  ) ) allocate (lon  (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(upn  ) ) allocate (upn  (ncells    ), stat = istat)
                   !
                   if ( istat == 0 .and. .not.allocated(sol0 ) ) allocate (sol0 (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(sol1 ) ) allocate (sol1 (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(amata) ) allocate (amata(ncells,0:3), stat = istat)
                   if ( istat == 0 .and. .not.allocated(rhsa ) ) allocate (rhsa (ncells    ), stat = istat)
                   if ( istat == 0 .and. .not.allocated(vt   ) ) allocate (vt   (ncells    ), stat = istat)
                endif
                !
                if ( istat /= 0 ) then
                   write (msgstr, '(a,i6)') 'allocation problem: variables for Newton iteration and return code is ',istat
                   call msgerr ( 4, trim(msgstr) )
                   return
                endif
                !
                amatn = 0d0
                rhsn  = 0d0
                lon   = 0d0
                upn   = 0d0
                !
                if ( optg /= 5 ) then
                   p0   = 1d0
                   p1   = 1d0
                   q0   = 0d0
                   q1   = 0d0
                   sola = 0d0
                else
                   sol0 = 0d0
                   sol1 = 0d0
                endif
                !
                amata = 0d0
                rhsa  = 0d0
                vt    = 0d0
                !
             endif
             !
             if ( ihydro == 3 .or. qmax == 1 ) then
                !
                if ( .not.allocated(cmat) ) allocate (cmat(mcgrd,5), stat = istat)
                !
                if ( istat /= 0 ) then
                   write (msgstr, '(a,i6)') 'allocation problem: variables for SIP solver and return code is ',istat
                   call msgerr ( 4, trim(msgstr) )
                   return
                endif
                !
                cmat = 0.
                !
             endif
             !
          endif
          !
       else
          !
          if ( .not.allocated(a) ) allocate (a(0))
          if ( .not.allocated(b) ) allocate (b(0))
          if ( .not.allocated(c) ) allocate (c(0))
          if ( .not.allocated(d) ) allocate (d(0))
          !
          if ( .not.allocated(z1) ) allocate (z1(0))
          if ( .not.allocated(z2) ) allocate (z2(0))
          !
          if ( .not.allocated(z12) ) allocate (z12(0))
          if ( .not.allocated(z22) ) allocate (z22(0))
          !
          if ( .not.allocated(an ) ) allocate (an (0))
          if ( .not.allocated(bn ) ) allocate (bn (0))
          if ( .not.allocated(cn ) ) allocate (cn (0))
          if ( .not.allocated(dn ) ) allocate (dn (0))
          if ( .not.allocated(lon) ) allocate (lon(0))
          if ( .not.allocated(upn) ) allocate (upn(0))
          !
          if ( .not.allocated(p0  ) ) allocate (p0  (0))
          if ( .not.allocated(p1  ) ) allocate (p1  (0))
          if ( .not.allocated(q0  ) ) allocate (q0  (0))
          if ( .not.allocated(q1  ) ) allocate (q1  (0))
          if ( .not.allocated(ba  ) ) allocate (ba  (0))
          if ( .not.allocated(da  ) ) allocate (da  (0))
          if ( .not.allocated(sola) ) allocate (sola(0))
          if ( .not.allocated(vt  ) ) allocate (vt  (0))
          !
          if ( .not.allocated(amat) ) allocate (amat(0,0))
          if ( .not.allocated(rhs ) ) allocate (rhs (0  ))
          if ( .not.allocated(resd) ) allocate (resd(0  ))
          !
          if ( .not.allocated(diag) ) allocate (diag(0))
          if ( .not.allocated(s   ) ) allocate (s   (0))
          if ( .not.allocated(r   ) ) allocate (r   (0))
          if ( .not.allocated(z   ) ) allocate (z   (0))
          !
          if ( .not.allocated(diag2) ) allocate (diag2(0))
          if ( .not.allocated(resd2) ) allocate (resd2(0))
          if ( .not.allocated(s2   ) ) allocate (s2   (0))
          if ( .not.allocated(r2   ) ) allocate (r2   (0))
          if ( .not.allocated(z3   ) ) allocate (z3   (0))
          !
          if ( .not.allocated(prec2) ) allocate (prec2(0))
          if ( .not.allocated(prec3) ) allocate (prec3(0))
          !
          if ( .not.allocated(ias ) ) allocate (ias (0))
          if ( .not.allocated(irws) ) allocate (irws(0,0))
          if ( .not.allocated(dis ) ) allocate (dis (0))
          if ( .not.allocated(jas ) ) allocate (jas (0))
          if ( .not.allocated(axs ) ) allocate (axs (0))
          if ( .not.allocated(axs2) ) allocate (axs2(0))
          !
          if ( .not.allocated(amatn) ) allocate (amatn(0,0))
          if ( .not.allocated(rhsn ) ) allocate (rhsn (0  ))
          !
          if ( .not.allocated(amata) ) allocate (amata(0,0))
          if ( .not.allocated(rhsa ) ) allocate (rhsa (0  ))
          if ( .not.allocated(sol0 ) ) allocate (sol0 (0  ))
          if ( .not.allocated(sol1 ) ) allocate (sol1 (0  ))
          !
          if ( .not.allocated(cmat) ) allocate (cmat(0,0))
          !
          if ( .not.allocated(iplu) ) allocate (iplu(0))
          !
       endif
       !
       if ( ihydro == 1 .or. ihydro == 2 ) then
          !
          if ( optg /= 5 ) then
             if ( ihydro == 1 ) then
                nconct = max(25,22+kpmax)
             else
                nconct = 25
             endif
          else
             if ( ihydro == 1 ) then
                nconct = max(19,16+kpmax)
             else
                nconct = 19
             endif
          endif
          !
          if ( ihydro == 1 .and. kpmax > 2 ) then
             if ( .not.allocated(ishif) ) allocate (ishif(0:kpmax-1), stat = istat)
          else
             if ( .not.allocated(ishif) ) allocate (ishif(0:2      ), stat = istat)
          endif
          !
          if ( optg /= 5 ) then
             if ( istat == 0 .and. .not.allocated(amatw) ) allocate (amatw(mcgrd,0:kpmax,3), stat = istat)
             if ( istat == 0 .and. .not.allocated(rhsw ) ) allocate (rhsw (mcgrd,0:kpmax  ), stat = istat)
          else
             if ( istat == 0 .and. .not.allocated(amatw) ) allocate (amatw(ncells,0:kpmax,3), stat = istat)
             if ( istat == 0 .and. .not.allocated(rhsw ) ) allocate (rhsw (ncells,0:kpmax  ), stat = istat)
          endif
          if ( optg /= 5 ) then
             if ( istat == 0 .and. .not.allocated(amatp) ) allocate (amatp(mcgrd,kpmax,nconct), stat = istat)
             if ( istat == 0 .and. .not.allocated(rhsp ) ) allocate (rhsp (mcgrd,kpmax       ), stat = istat)
          else
             if ( istat == 0 .and. .not.allocated(amatp) ) allocate (amatp(ncells,kpmax,0:nconct), stat = istat)
             if ( istat == 0 .and. .not.allocated(rhsp ) ) allocate (rhsp (ncells,kpmax         ), stat = istat)
          endif
          !
          if ( optg /= 5 ) then
             if ( icond == 1 .or. icond == 2 ) then
                if ( istat == 0 .and. .not.allocated(diagb) ) allocate (diagb(mcgrd,-1:qmax+nconct-23), stat = istat)
             endif
             if ( icond == 2 .or. icond == 3 ) then
                if ( istat == 0 .and. .not.allocated(prec ) ) allocate (prec (mcgrd,-1:qmax+2, nconct), stat = istat)
             endif
             if ( istat == 0 .and. .not.allocated(res  ) ) allocate (res  (mcgrd,   qmax          ), stat = istat)
             if ( istat == 0 .and. .not.allocated(res0 ) ) allocate (res0 (mcgrd,   qmax          ), stat = istat)
             if ( istat == 0 .and. .not.allocated(p    ) ) allocate (p    (mcgrd,   qmax          ), stat = istat)
             if ( istat == 0 .and. .not.allocated(v    ) ) allocate (v    (mcgrd,   qmax          ), stat = istat)
             if ( istat == 0 .and. .not.allocated(w    ) ) allocate (w    (mcgrd,-1:qmax+nconct-23), stat = istat)
             if ( istat == 0 .and. .not.allocated(u    ) ) allocate (u    (mcgrd,   qmax          ), stat = istat)
             if ( istat == 0 .and. .not.allocated(t    ) ) allocate (t    (mcgrd,   qmax          ), stat = istat)
             if ( istat == 0 .and. .not.allocated(vt2  ) ) allocate (vt2  (mcgrd,   qmax          ), stat = istat)
          else
             !
             ! integer arrays to be used for the CSR format
             !
             if ( .not.allocated(ia) ) allocate (ia(ncells*qmax+1), stat = istat)
             if ( istat == 0 .and. .not.allocated(irw) ) allocate (irw(           ncells,1:4 ), stat = istat)
             if ( istat == 0 .and. .not.allocated(ja ) ) allocate (ja ((nconct+1)*ncells*qmax), stat = istat)
             !
             if ( istat == 0 .and. .not.allocated(di   ) ) allocate (di   (ncells*qmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(iplu3) ) allocate (iplu3(ncells*qmax), stat = istat)
             !
             ! arrays to be used in bicgstab3 routine
             !
             if ( istat == 0 .and. .not.allocated(res ) ) allocate (res (ncells,qmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(res0) ) allocate (res0(ncells,qmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(p   ) ) allocate (p   (ncells,qmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(v   ) ) allocate (v   (ncells,qmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(u   ) ) allocate (u   (ncells,qmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(t   ) ) allocate (t   (ncells,qmax), stat = istat)
             !
             if ( istat == 0 .and. .not.allocated(prec4) ) allocate (prec4(0:(nconct+1)*ncells*qmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(w2   ) ) allocate (w2  (              ncells*qmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(ax2  ) ) allocate (ax2 (   (nconct+1)*ncells*qmax), stat = istat)
             !
             if ( istat == 0 .and. .not.allocated(sol3) ) then
                allocate (sol3(ncells*qmax), stat = istat)
             else if ( size(sol3) == 0 ) then
                deallocate(sol3)
                allocate (sol3(ncells*qmax), stat = istat)
             endif
             !
          endif
          !
          if ( istat /= 0 ) then
             write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
             call msgerr ( 4, trim(msgstr) )
             return
          endif
          !
          if ( optg /= 5 ) then
             ishif(0) = 1
             ishif(1) = 5
             ishif(2) = 17
             !
             if ( ihydro == 1 ) then
                !
                do i = 3, kpmax-1
                   ishif(i) = 23+i
                enddo
                !
             endif
          else
             ishif(0) = 0
             ishif(1) = 12
             ishif(2) = 16
             !
             if ( ihydro == 1 ) then
                !
                do i = 3, kpmax-1
                   ishif(i) = 17+i
                enddo
                !
             endif
          endif
          !
          amatw = 0.
          rhsw  = 0.
          amatp = 0.
          rhsp  = 0.
          if ( optg /= 5 ) then
             if ( icond == 1 .or. icond == 2 ) diagb = 1.
             if ( icond == 2 .or. icond == 3 ) prec  = 0.
             if ( icond == 3 .and. pnums(28) < 0. ) prec(:,:,1) = 1.
             res   = 0.
             res0  = 0.
             p     = 0.
             v     = 0.
             w     = 0.
             u     = 0.
             t     = 0.
             vt2   = 0.
          else
             iplu3 = 0
             prec4 = 0.
             res   = 0.
             res0  = 0.
             p     = 0.
             v     = 0.
             u     = 0.
             t     = 0.
             w2    = 0.
             ax2   = 0.
             sol3  = 0.
             !
             ! fill CSR arrays
             !
             call csrf ( ia, ja, irw, di, qmax, nconct )
          endif
          !
       else
          !
          if ( .not.allocated(ishif) ) allocate (ishif(0    ))
          if ( .not.allocated(amatw) ) allocate (amatw(0,0,0))
          if ( .not.allocated(rhsw ) ) allocate (rhsw (0,0  ))
          if ( .not.allocated(amatp) ) allocate (amatp(0,0,0))
          if ( .not.allocated(rhsp ) ) allocate (rhsp (0,0  ))
          !
          if ( .not.allocated(prec ) ) allocate (prec (0,0,0))
          if ( .not.allocated(diagb) ) allocate (diagb(0,0  ))
          if ( .not.allocated(res  ) ) allocate (res  (0,0  ))
          if ( .not.allocated(res0 ) ) allocate (res0 (0,0  ))
          if ( .not.allocated(p    ) ) allocate (p    (0,0  ))
          if ( .not.allocated(v    ) ) allocate (v    (0,0  ))
          if ( .not.allocated(w    ) ) allocate (w    (0,0  ))
          if ( .not.allocated(u    ) ) allocate (u    (0,0  ))
          if ( .not.allocated(t    ) ) allocate (t    (0,0  ))
          if ( .not.allocated(vt2  ) ) allocate (vt2  (0,0  ))
          !
          if ( .not.allocated(ia ) ) allocate (ia (0))
          if ( .not.allocated(ja ) ) allocate (ja (0))
          !
          if ( .not.allocated(iplu3) ) allocate (iplu3(0))
          if ( .not.allocated(prec4) ) allocate (prec4(0))
          if ( .not.allocated(ax2  ) ) allocate (ax2  (0))
          if ( .not.allocated(w2   ) ) allocate (w2   (0))
          if ( .not.allocated(sol3 ) ) allocate (sol3 (0))
          !
       endif
       !
       if ( itrans /= 0 ) then
          !
          if ( optg /= 5 ) then
             if ( .not.allocated(amatc) ) allocate (amatc(mcgrd,kmax,3), stat = istat)
             if ( istat == 0 .and. .not.allocated(rhsc) ) allocate (rhsc(mcgrd,kmax), stat = istat)
          else
             if ( .not.allocated(amatc) ) allocate (amatc(ncells,kmax,3), stat = istat)
             if ( istat == 0 .and. .not.allocated(rhsc) ) allocate (rhsc(ncells,kmax), stat = istat)
          endif
          !
          if ( istat /= 0 ) then
             write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
             call msgerr ( 4, trim(msgstr) )
             return
          endif
          !
          amatc = 0.
          rhsc  = 0.
          !
       else
          !
          if ( .not.allocated(amatc) ) allocate (amatc(0,0,0))
          if ( .not.allocated(rhsc ) ) allocate (rhsc (0,0  ))
          !
       endif
       !
       if ( iturb /= 0 ) then
          !
          if ( optg /= 5 ) then
             if ( oned ) then
                if ( .not.allocated(amatt) ) allocate (amatt(mcgrd,0:kmax,3), stat = istat)
             else
                if ( .not.allocated(amatt) ) allocate (amatt(mcgrd,0:kmax,5), stat = istat)
             endif
             if ( istat == 0 .and. .not.allocated(rhst) ) allocate (rhst(mcgrd,0:kmax), stat = istat)
          else
             if ( .not.allocated(amatt) ) allocate (amatt(ncells,0:kmax,3), stat = istat)
             if ( istat == 0 .and. .not.allocated(rhst) ) allocate (rhst(ncells,0:kmax), stat = istat)
          endif
          !
          if ( istat /= 0 ) then
             write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
             call msgerr ( 4, trim(msgstr) )
             return
          endif
          !
          amatt        = 0.
          amatt(:,:,1) = 1.
          rhst         = 1.e-7
          !
       else
          !
          if ( .not.allocated(amatt) ) allocate (amatt(0,0,0))
          if ( .not.allocated(rhst ) ) allocate (rhst (0,0  ))
          !
       endif
       !
       if ( mimetic ) then
          !
          if ( .not.allocated(bdx) ) allocate (bdx(mcgrd,kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(bux) ) allocate (bux(mcgrd,kmax), stat = istat)
          if ( .not.oned ) then
             if ( istat == 0 .and. .not.allocated(bdy) ) allocate (bdy(mcgrd,kmax), stat = istat)
             if ( istat == 0 .and. .not.allocated(buy) ) allocate (buy(mcgrd,kmax), stat = istat)
          endif
          if ( istat == 0 .and. .not.allocated(wrk) ) allocate (wrk(mcgrd,kmax), stat = istat)
          !
          if ( istat /= 0 ) then
             write (msgstr, '(a,i6)') 'allocation problem: variables for system of equations and return code is ',istat
             call msgerr ( 4, trim(msgstr) )
             return
          endif
          !
          bdx = 0.
          bux = 0.
          if ( .not.oned ) then
             bdy = 0.
             buy = 0.
          endif
          wrk = 0.
          !
       else
          !
          if ( .not.allocated(bdx) ) allocate (bdx(0,0))
          if ( .not.allocated(bux) ) allocate (bux(0,0))
          if ( .not.allocated(bdy) ) allocate (bdy(0,0))
          if ( .not.allocated(buy) ) allocate (buy(0,0))
          !
       endif
       !
    endif
    !
    ! allocate and initialize data for rigid body solver
    !
    if ( ifloat == 2 ) then
       !
       if ( istat == 0 .and. .not.allocated(fbod0) ) allocate (fbod0(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(tbod0) ) allocate (tbod0(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(fbod1) ) allocate (fbod1(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(tbod1) ) allocate (tbod1(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(fhyd0) ) allocate (fhyd0(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(thyd0) ) allocate (thyd0(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(fhyd1) ) allocate (fhyd1(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(thyd1) ) allocate (thyd1(mbod,ndim), stat = istat)
       !
       if ( istat == 0 .and. .not.allocated(afot0) ) allocate (afot0(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(afor0) ) allocate (afor0(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(vfot0) ) allocate (vfot0(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(vfor0) ) allocate (vfor0(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(xfot0) ) allocate (xfot0(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(xfor0) ) allocate (xfor0(mbod,ndim), stat = istat)
       !
       if ( istat == 0 .and. .not.allocated(afot1) ) allocate (afot1(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(afor1) ) allocate (afor1(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(vfot1) ) allocate (vfot1(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(vfor1) ) allocate (vfor1(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(xfot1) ) allocate (xfot1(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(xfor1) ) allocate (xfor1(mbod,ndim), stat = istat)
       !
       if ( istat == 0 .and. .not.allocated(xfoti) ) allocate (xfoti(mbod,ndim), stat = istat)
       if ( istat == 0 .and. .not.allocated(xfori) ) allocate (xfori(mbod,ndim), stat = istat)
       !
    else
       !
       if ( .not.allocated(fbod0) ) allocate (fbod0(0,0))
       if ( .not.allocated(tbod0) ) allocate (tbod0(0,0))
       if ( .not.allocated(fbod1) ) allocate (fbod1(0,0))
       if ( .not.allocated(tbod1) ) allocate (tbod1(0,0))
       if ( .not.allocated(fhyd0) ) allocate (fhyd0(0,0))
       if ( .not.allocated(thyd0) ) allocate (thyd0(0,0))
       if ( .not.allocated(fhyd1) ) allocate (fhyd1(0,0))
       if ( .not.allocated(thyd1) ) allocate (thyd1(0,0))
       !
       if ( .not.allocated(afot0) ) allocate (afot0(0,0))
       if ( .not.allocated(afor0) ) allocate (afor0(0,0))
       if ( .not.allocated(vfot0) ) allocate (vfot0(0,0))
       if ( .not.allocated(vfor0) ) allocate (vfor0(0,0))
       if ( .not.allocated(xfot0) ) allocate (xfot0(1,1))
       if ( .not.allocated(xfor0) ) allocate (xfor0(1,1))
       !
       if ( .not.allocated(afot1) ) allocate (afot1(0,0))
       if ( .not.allocated(afor1) ) allocate (afor1(0,0))
       if ( .not.allocated(vfot1) ) allocate (vfot1(0,0))
       if ( .not.allocated(vfor1) ) allocate (vfor1(0,0))
       if ( .not.allocated(xfot1) ) allocate (xfot1(1,1))
       if ( .not.allocated(xfor1) ) allocate (xfor1(1,1))
       !
       if ( .not.allocated(xfoti) ) allocate (xfoti(1,1))
       if ( .not.allocated(xfori) ) allocate (xfori(1,1))
       !
    endif
    !
    if ( ifloat == 2 ) then
       !
       fbod0 = 0.
       tbod0 = 0.
       fbod1 = 0.
       tbod1 = 0.
       fhyd0 = 0.
       thyd0 = 0.
       fhyd1 = 0.
       thyd1 = 0.
       !
       afot0 = 0.
       afor0 = 0.
       vfot0 = 0.
       vfor0 = 0.
       xfot0 = 0.
       xfor0 = 0.
       !
       afot1 = 0.
       afor1 = 0.
       vfot1 = 0.
       vfor1 = 0.
       xfot1 = 0.
       xfor1 = 0.
       !
       xfoti = 0.
       xfori = 0.
       !
    else
       !
       xfot0 = 0.
       xfor0 = 0.
       xfot1 = 0.
       xfor1 = 0.
       xfoti = 0.
       xfori = 0.
       !
    endif
    !
    ! allocate and initialize wave output
    !
    if ( lwavoutp ) then
       !
       if ( optg /= 5 ) then
          if ( .not.allocated(setup) ) allocate (setup(mcgrd), stat = istat)
          if ( istat == 0 .and. .not.allocated(etavar) ) allocate (etavar(mcgrd), stat = istat)
          if ( istat == 0 .and. .not.allocated(hsig  ) ) allocate (hsig  (mcgrd), stat = istat)
       else
          if ( .not.allocated(setup) ) allocate (setup(ncells), stat = istat)
          if ( istat == 0 .and. .not.allocated(etavar) ) allocate (etavar(ncells), stat = istat)
          if ( istat == 0 .and. .not.allocated(hsig  ) ) allocate (hsig  (ncells), stat = istat)
       endif
       !
    else
       !
       if ( .not.allocated(setup ) ) allocate (setup (0))
       if ( .not.allocated(etavar) ) allocate (etavar(0))
       if ( .not.allocated(hsig  ) ) allocate (hsig  (0))
       !
    endif
    !
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: wave output and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    if ( lwavoutp ) then
       !
       setup  = 0.
       etavar = 0.
       hsig   = 0.
       !
    endif
    !
    ! allocate and initialize mean current output
    !
    if ( lcuroutp ) then
       !
       if ( optg /= 5 ) then
          if ( .not.allocated(mvelu) ) allocate (mvelu(mcgrd,kmax), stat = istat)
          if ( .not.oned ) then
             if ( istat == 0 .and. .not.allocated(mvelv) ) allocate (mvelv(mcgrd,kmax), stat = istat)
          endif
       else
          if ( .not.allocated(mvelu) ) allocate (mvelu(nfaces,kmax), stat = istat)
       endif
       !
    else
       !
       if ( .not.allocated(mvelu) ) allocate (mvelu(0,0))
       if ( .not.allocated(mvelv) ) allocate (mvelv(0,0))
       !
    endif
    !
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: mean current output and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    if ( lcuroutp ) then
       !
       mvelu = 0.
       if ( .not.oned .and. optg /= 5 ) mvelv = 0.
       !
    endif
    !
    ! allocate and initialize mean transport constituent output
    !
    if ( ltraoutp ) then
       !
       if ( optg /= 5 ) then
          if ( .not.allocated(mcons) ) allocate (mcons(mcgrd,kmax,ltrans), stat = istat)
       else
          if ( .not.allocated(mcons) ) allocate (mcons(ncells,kmax,ltrans), stat = istat)
       endif
       !
    else
       !
       if ( .not.allocated(mcons) ) allocate (mcons(0,0,0))
       !
    endif
    !
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: mean transport constituent output and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    if ( ltraoutp ) then
       !
       mcons = 0.
       !
    endif
    !
    ! allocate and initialize mean turbulence output
    !
    if ( lturoutp ) then
       !
       if ( optg /= 5 ) then
          if ( .not.allocated(mtke) ) allocate (mtke(mcgrd,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(meps) ) allocate (meps(mcgrd,0:kmax), stat = istat)
       else
          if ( .not.allocated(mtke) ) allocate (mtke(ncells,0:kmax), stat = istat)
          if ( istat == 0 .and. .not.allocated(meps) ) allocate (meps(ncells,0:kmax), stat = istat)
       endif
       !
    else
       !
       if ( .not.allocated(mtke) ) allocate (mtke(0,0))
       if ( .not.allocated(meps) ) allocate (meps(0,0))
       !
    endif
    !
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: mean turbulence output and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    if ( lturoutp ) then
       !
       mtke = 0.
       meps = 0.
       !
    endif
    !
    ! allocate and initialize physical vertical velocity for output purposes (only meant for hydrostatic computation)
    !
    if ( ihydro == 0 .or. ihydro == 3 ) then
       !
       if ( optg /= 5 ) then
          if ( .not.allocated(wphy) ) allocate (wphy(mcgrd,0:kmax), stat = istat)
       else
          if ( .not.allocated(wphy) ) allocate (wphy(ncells,0:kmax), stat = istat)
       endif
       !
    else
       !
       if ( .not.allocated(wphy) ) allocate (wphy(0,0))
       !
    endif
    !
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: physical vertical velocity and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    if ( ihydro == 0 .or. ihydro == 3 ) then
       !
       wphy = 0.
       !
    endif
    !
    ! initialize water depths
    !
!TIMG    call SWTSTA(52)
    if ( optg /= 5 ) then
       if ( kmax == 1 ) then
          call SwashUpdateDepths (   u1,   v1 )
       else
          call SwashUpdateDepths ( udep, vdep )
       endif
       humn = hum
       if ( .not.oned ) hvmn = hvm
    else
       if ( kmax == 1 ) then
          call SwashUpdateUDepths ( u1(:,1) )
       else
          call SwashUpdateUDepths ( udep )
       endif
       humn = hum
    endif
!TIMG    call SWTSTO(52)
    if (STPNOW()) return
    !
    ! initialize layer interfaces
    !
!TIMG    call SWTSTA(53)
    if ( kmax > 1 ) then
       if ( optg /= 5 ) then
          call SwashLayerIntfaces
       else
          call SwashLayUIntfaces
       endif
    endif
!TIMG    call SWTSTO(53)
    !
    ! initialize mask arrays for wetting and drying
    !
!TIMG    call SWTSTA(54)
    if ( optg /= 5 ) then
       call SwashDryWet
    else
       call SwashUDryWet
    endif
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! update mask arrays for free surface / pressurized flow
    !
!TIMG    call SWTSTA(54)
    if ( optg /= 5 .and. ifloat /= 0 ) call SwashPresFlow
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! initialize bottom friction coefficient
    !
!TIMG    call SWTSTA(55)
    if ( irough /= 0 .and. ( irough /= 4 .or. kmax == 1 ) ) then
       if ( optg /= 5 ) then
          call SwashBotFrict ( u1, v1 )
       else
          call SwashUBotFrict ( u1 )
       endif
    endif
!TIMG    call SWTSTO(55)
    if (STPNOW()) return
    !
    ! in case of moving bodies initialize forces and moments
    !
!TIMG    call SWTSTA(62)
    call SwashForcesRigidBod
!TIMG    call SWTSTO(62)
    if (STPNOW()) return
    !
    ! in case of (quasi-)steady flow initialize flow velocities based on Chezy formula
    !
    if ( optg /= 5 .and. instead ) call SwashInitSteady
    if (STPNOW()) return
    !
    ! in case of subgrid approach initialize flow velocities on coarse (pressure) grid
    !
    if ( lsubg .and. ihydro /= 3 ) then
       !
       u1p = 0.
       if ( .not.oned ) v1p = 0.
       !
       k = 0
       !
       do i = 1, kpmax
          !
          do j = 1, npu(i)
             !
             k = k + 1
             !
             u1p(:,i) = u1p(:,i) + hku(:,k)*u1(:,k)
             !
             if ( .not.oned ) v1p(:,i) = v1p(:,i) + hkv(:,k)*v1(:,k)
             !
          enddo
          !
          u1p(:,i) = u1p(:,i) / hkuc(:,i)
          !
          if ( .not.oned ) v1p(:,i) = v1p(:,i) / hkvc(:,i)
          !
       enddo
       !
    endif
    !
end subroutine SwashCheckPrep
