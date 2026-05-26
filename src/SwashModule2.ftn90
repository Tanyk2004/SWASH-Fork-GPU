!                 ALLOCATABLE DATA RELATED MODULES, file 2 of 2
!
!     Contents of this file
!
!     outp_data          information for output data
!     m_bndspec          information for boundary condition specifications
!     m_genarr           contains a number of general arrays

module outp_data
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
!   41.95: Marcel Zijlema
!
!   Updates
!
!    1.00, February 2010: New Module
!   41.95,     July 2022: introduction VTK and PVD attributes
!
!   Purpose
!
!   Contains data needed during generation of output
!
!   Method
!
!   MODULE construct
!
!   Modules used
!
    use ocpcomm2
!
    implicit none
!
!   Module parameters
!
    integer, parameter               :: max_outp_req = 250         ! maximum number of output requests
    integer, parameter               :: mopa         = 31          ! maximum number of help arrays for outputting
    integer, parameter               :: mopak        = 30          ! maximum number of layer-dependent help arrays for outputting
!
!   Module variables
!
    integer                          :: dec_block =  6             ! number of decimals for fixed-point block
    integer                          :: fld_table = 15             ! field length for fixed-point table
    integer, save                    :: nreoq     =  0             ! actual number of requests saved
    integer, save                    :: ntvtk(1:max_outp_req) = -1 ! timestep counter for time-varying VTK files
    integer                          :: upvdf(1:max_outp_req)      ! unit number PVD files
    !
    real                             :: alpobj                     ! (counter clockwise) rotation angle of axis of floating body relative to computational grid
    real                             :: cogx                       ! x-coordinate of center of gravity of floating body relative to computational grid
    real                             :: cogy                       ! y-coordinate of center of gravity of floating body relative to computational grid
    real                             :: cogz                       ! z-coordinate of center of gravity of floating body relative to computational grid
    real                             :: delrp                      ! threshold depth for runup height calculation (for output quantity RUNUP)
    real                             :: hrunp                      ! depth where inundation takes place (for output quantity HRUN)
    real                             :: ncuroutp                   ! number of steps during integration of current for output
    real                             :: ntraoutp                   ! number of steps during integration of concentrations for output
    real                             :: nturoutp                   ! number of steps during integration of turbulence quantities for output
    real                             :: nwavoutp                   ! number of steps during integration of surface elevation for output
    real                             :: tcuroutp                   ! duration over which the mean current is computed
    real                             :: ttraoutp                   ! duration over which the mean constituent is computed
    real                             :: tturoutp                   ! duration over which the mean turbulence parameters are computed
    real                             :: twavoutp                   ! duration over which the wave output parameters are computed
    !
    logical                          :: lcompgrd                   ! indicates whether computational grid is used as output frame
    logical                          :: lcuroutp                   ! indicates whether mean current output is requested or not
    logical                          :: ltraoutp                   ! indicates whether mean constituent output is requested or not
    logical                          :: lturoutp                   ! indicates whether mean turbulence output is requested or not
    logical                          :: lwavoutp                   ! indicates whether wave output is requested or not
    !
    character (len=1)                :: out_comment = '%'          ! comment sign for heading lines
    !
    ! formats for output
    character (len=40)               :: flt_block = '(6e14.6)'     ! floating point block
    character (len=40)               :: flt_table = '(e15.8)'      ! floating point table
    character (len=40)               :: flt_swash = '(10e17.9)'    ! floating point for SWASH/STAB formatted table
    !
    ! format for block output per process in case of collecting data
    character (len=40)               :: flt_blkp = '(6e17.9)'      !
    !
    ! longer filenames for output requests
    character (len=lenfnm)           :: outp_files(1:max_outp_req) !
    !
    ! output directory containing a series of time-varying VTK files
    character (len=lenfnm)           :: vtkdir(1:max_outp_req)
    !
    ! textline containing VTK line in XML format
    character(len=1024)              :: vtkline
    !
    ! XML header lines
    character(len=25)                :: xmllin1 = '<?xml version="1.0"?>'
    character(len= 5)                :: xmllin2 = '<!--'
    character(len= 5)                :: xmllin3 = '-->'
    !
    ! default lines of PVD file
    character(len=80)                :: pvdlin1 = '<VTKFile type="Collection" '//              &
                                                  'version="0.1" byte_order="LittleEndian">'
    character(len=15)                :: pvdlin2 = '  <Collection>'
    character(len=15)                :: pvdlin3 = '  </Collection>'
    character(len=10)                :: pvdlin4 = '</VTKFile>'
    !
    real*8, dimension(2)             :: etorq =  9999.             ! indicates end time and final time step for output requests
    real*8, dimension(max_outp_req)  :: ourqt = -9999.             ! array indicating at what time requested output is processed
    !
    real, dimension(:,:), save, allocatable :: oparr               ! help array to save some quantities for outputting
                                                                   ! = 1; water depth
                                                                   ! = 2; water level
                                                                   ! = 3; contravariant U-velocity component
                                                                   ! = 4; contravariant V-velocity component
                                                                   ! = 5; Cartesian U-velocity component
                                                                   ! = 6; Cartesian V-velocity component
                                                                   ! = 7; contravariant U-discharge component
                                                                   ! = 8; contravariant V-discharge component
                                                                   ! = 9; Cartesian U-discharge component
                                                                   ! =10; Cartesian V-discharge component
                                                                   ! =11; U-component of friction velocity
                                                                   ! =12; V-component of friction velocity
                                                                   ! =13; non-hydrostatic pressure at bottom
                                                                   ! =14; pressure at bottom
                                                                   ! =15; auxiliary array
                                                                   ! =16; inundation depth
                                                                   ! =17; wave-induced setup
                                                                   ! =18; significant wave height
                                                                   ! =19; RMS wave height
                                                                   ! =20; time-averaged contravariant U-velocity component
                                                                   ! =21; time-averaged contravariant V-velocity component
                                                                   ! =22; time-averaged Cartesian U-velocity component
                                                                   ! =23; time-averaged Cartesian V-velocity component
                                                                   ! =24; wave breaking point
                                                                   ! =25; salinity
                                                                   ! =26; temperature
                                                                   ! =27; suspended sediment
                                                                   ! =28; time-averaged salinity
                                                                   ! =29; time-averaged temperature
                                                                   ! =30; time-averaged suspended sediment
                                                                   ! =31; vorticity
    real, dimension(:,:,:), save, allocatable :: oparrk            ! help array to save some layer-dependent quantities for outputting
                                                                   ! = 1; layer thickness
                                                                   ! = 2; layer interface
                                                                   ! = 3; contravariant U-velocity component
                                                                   ! = 4; contravariant V-velocity component
                                                                   ! = 5; Cartesian U-velocity component
                                                                   ! = 6; Cartesian V-velocity component
                                                                   ! = 7; contravariant U-discharge component
                                                                   ! = 8; contravariant V-discharge component
                                                                   ! = 9; Cartesian U-discharge component
                                                                   ! =10; Cartesian V-discharge component
                                                                   ! =11; relative vertical velocity
                                                                   ! =12; physical vertical velocity
                                                                   ! =13; non-hydrostatic pressure
                                                                   ! =14; pressure
                                                                   ! =15; turbulent kinetic energy
                                                                   ! =16; dissipation rate of turbulent energy
                                                                   ! =17; eddy viscosity
                                                                   ! =18; time-averaged contravariant U-velocity component
                                                                   ! =19; time-averaged contravariant V-velocity component
                                                                   ! =20; time-averaged Cartesian U-velocity component
                                                                   ! =21; time-averaged Cartesian V-velocity component
                                                                   ! =22; time-averaged turbulent kinetic energy
                                                                   ! =23; time-averaged dissipation rate of turbulent energy
                                                                   ! =24; time-averaged eddy viscosity
                                                                   ! =25; salinity
                                                                   ! =26; temperature
                                                                   ! =27; suspended sediment
                                                                   ! =28; time-averaged salinity
                                                                   ! =29; time-averaged temperature
                                                                   ! =30; time-averaged suspended sediment
    !
    real, dimension(:)    , save, allocatable :: etavar            ! variance of the surface elevation
    real, dimension(:)    , save, allocatable :: hsig              ! significant wave height
    real, dimension(:,:,:), save, allocatable :: mcons             ! time-averaged or mean constituents
    real, dimension(:,:)  , save, allocatable :: meps              ! time-averaged or mean dissipation rate of turbulent energy
    real, dimension(:,:)  , save, allocatable :: mtke              ! time-averaged or mean turbulent kinetic energy
    real, dimension(:,:)  , save, allocatable :: mvelu             ! time-averaged or mean U-velocity component
    real, dimension(:,:)  , save, allocatable :: mvelv             ! time-averaged or mean V-velocity component
    real, dimension(:)    , save, allocatable :: setup             ! wave-induced setup
    !
    type opsdat
       character (len=1)     :: pstype                             ! type (F(rame), C(urve), P(oints), ...)
       character (len=8)     :: psname                             ! name of point set
       integer               :: opi(2)                             ! integer coefficients
       real                  :: opr(6)                             ! real coefficients
       integer               :: mip                                ! number of points
       real, pointer         :: xp(:), yp(:), xq(:), yq(:)         ! point coordinates
       type(opsdat), pointer :: nextops                            ! pointer to next point set in list
    end type opsdat
    !
    type(opsdat), save, target  :: fops                            ! first item in list of point sets
    type(opsdat), save, pointer :: cops                            ! current item in list of point sets
    logical, save               :: lops = .false.                  ! indicates whether there is a list of point sets
    !
    type orqdat
       character (len=4)      :: rqtype                            ! type (BLK or TAB)
       character (len=8)      :: psname                            ! name of point set
       integer                :: oqi(4)                            ! integer coefficients
       real*8                 :: oqr(2)                            ! real coefficients
       integer, pointer       :: ivtyp(:)                          ! type of output variable
       real, pointer          :: fac(:)                            ! multiplication factor of block output
       type(orqdat), pointer  :: nextorq                           ! pointer to next output request in list
    end type orqdat
    !
    type(orqdat), save, target  :: forq                            ! first item in list of output requests
    logical, save               :: lorq = .false.                  ! indicates whether there is a list of output requests
!
!   Source text
!
    interface delete
      module procedure deleteops
      module procedure deleteorq
    end interface
    !
    contains

    recursive subroutine deleteops ( ops )

    type(opsdat) :: ops

    if ( lops ) then
       if ( associated( ops%xp ) ) then
          deallocate( ops%xp )
          nullify( ops%xp )
       endif
       if ( associated( ops%yp ) ) then
          deallocate( ops%yp )
          nullify( ops%yp )
       endif
       if ( ops%pstype.eq.'R' .and. associated( ops%xq ) ) then
          deallocate( ops%xq )
          nullify( ops%xq )
       endif
       if ( ops%pstype.eq.'R' .and. associated( ops%yq ) ) then
          deallocate( ops%yq )
          nullify( ops%yq )
       endif
       if ( associated( ops%nextops ) ) then
          call deleteops ( ops%nextops )
          deallocate( ops%nextops )
          nullify( ops%nextops )
       endif
    endif

    end subroutine
    !
    recursive subroutine deleteorq ( orq )

    type(orqdat) :: orq

    if ( lorq ) then
       if ( associated( orq%ivtyp ) ) then
          deallocate( orq%ivtyp )
          nullify( orq%ivtyp )
       endif
       if ( associated( orq%fac ) ) then
          deallocate( orq%fac )
          nullify( orq%fac )
       endif
       if ( associated( orq%nextorq ) ) then
          call deleteorq ( orq%nextorq )
          deallocate( orq%nextorq )
          nullify( orq%nextorq )
       endif
    endif

    end subroutine
    !
end module outp_data

module m_bndspec
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
!    1.00, February 2010: New Module
!
!   Purpose
!
!   Contains data with respect to specification of boundary conditions
!
!   Method
!
!   MODULE construct
!
!   Modules used
!
    use SwashCommdata3, only: iwvgen
!
    implicit none
!
!   Module variables
!
    integer                  :: ibound                        ! bound long wave is added (=1) or not added (=0)
    integer                  :: nseed                         ! number of seed integers
    !
    integer, dimension(3)    :: spshape = 2                   ! =1; indicates option for value of given wave height concerning the incident spectrum
                                                              !     = 1; root-mean-square (rms) value of wave height
                                                              !     = 2; significant wave height
                                                              ! =2; indicates option for computation of frequency distribution in the incident spectrum
                                                              !     = 1; Pierson Moskowitz
                                                              !     = 2; Jonswap
                                                              !     = 3; TMA
                                                              ! =3; indicates option for computation of directional distribution in the incident spectrum
                                                              !     = 1; directional spreading in degrees is given
                                                              !     = 2; power of cosine is given
    !
    real                     :: gamma = 3.3                   ! peak-enhancement factor
    !
    logical                  :: alobnd                        ! if true, user has specified boundary conditions
    logical                  :: bcperx = .false.              ! indicates whether incident wave at boundary along x-direction must make periodic or not
    logical                  :: bcpery = .false.              ! indicates whether incident wave at boundary along y-direction must make periodic or not
    logical, save            :: lseed  = .true.               ! indicates to initialize the seed or not
    !
    integer   , dimension(:)    , save, allocatable :: seedf  ! seed used to generate pseudorandom phases
    integer   , dimension(:)    , save, allocatable :: seedt  ! seed used to generate pseudorandom spectral directions
    !
    real      , dimension(:,:)  , save, allocatable :: kwave  ! wave number of a Fourier component
    real      , dimension(:,:)  , save, allocatable :: comp1  ! first free wave component for synthesizing time series
    real      , dimension(:,:)  , save, allocatable :: comp2  ! second free wave component for synthesizing time series
    complex(8), dimension(:,:)  , save, allocatable :: fluxbu ! mass flux of bound long wave for U-velocity component
    complex(8), dimension(:,:)  , save, allocatable :: fluxbv ! mass flux of bound long wave for V-velocity component
    real      , dimension(:,:,:), save, allocatable :: stku1c ! first order Stokes component of horizontal velocity associated with cosine of phase
    real      , dimension(:,:,:), save, allocatable :: stku2c ! second order Stokes component of horizontal velocity associated with cosine of phase
    real      , dimension(:,:,:), save, allocatable :: stku1s ! first order Stokes component of horizontal velocity associated with sine of phase
    real      , dimension(:,:,:), save, allocatable :: stku2s ! second order Stokes component of horizontal velocity associated with sine of phase
    real      , dimension(:,:)  , save, allocatable :: stkz1c ! first order Stokes component of surface elevation associated with cosine of phase
    real      , dimension(:,:)  , save, allocatable :: stkz2c ! second order Stokes component of surface elevation associated with cosine of phase
    real      , dimension(:,:)  , save, allocatable :: stkz1s ! first order Stokes component of surface elevation associated with sine of phase
    real      , dimension(:,:)  , save, allocatable :: stkz2s ! second order Stokes component of surface elevation associated with sine of phase
    real      , dimension(:,:,:), save, allocatable :: subuc  ! sub-harmonic transfer function of U-velocity component associated with cosine of phase
    real      , dimension(:,:,:), save, allocatable :: supuc  ! super-harmonic transfer function of U-velocity component associated with cosine of phase
    real      , dimension(:,:,:), save, allocatable :: subus  ! sub-harmonic transfer function of U-velocity component associated with sine of phase
    real      , dimension(:,:,:), save, allocatable :: supus  ! super-harmonic transfer function of U-velocity component associated with sine of phase
    real      , dimension(:,:,:), save, allocatable :: subvc  ! sub-harmonic transfer function of V-velocity component associated with cosine of phase
    real      , dimension(:,:,:), save, allocatable :: supvc  ! super-harmonic transfer function of V-velocity component associated with cosine of phase
    real      , dimension(:,:,:), save, allocatable :: subvs  ! sub-harmonic transfer function of V-velocity component associated with sine of phase
    real      , dimension(:,:,:), save, allocatable :: supvs  ! super-harmonic transfer function of V-velocity component associated with sine of phase
    real      , dimension(:,:)  , save, allocatable :: subzc  ! sub-harmonic transfer function of surface elevation associated with cosine of phase
    real      , dimension(:,:)  , save, allocatable :: supzc  ! super-harmonic transfer function of surface elevation associated with cosine of phase
    real      , dimension(:,:)  , save, allocatable :: subzs  ! sub-harmonic transfer function of surface elevation associated with sine of phase
    real      , dimension(:,:)  , save, allocatable :: supzs  ! super-harmonic transfer function of surface elevation associated with sine of phase
    complex(8), dimension(:,:)  , save, allocatable :: zetab  ! surface elevation of bound long wave
    !
    type bfldat
       integer               :: bfiled(20)                    ! data concerning boundary condition files
                                                              ! 1.  status; 0: stationary, 1: nonstationary, -1: exhausted
                                                              ! 2.  not used
                                                              ! 3.  not used
                                                              ! 4.  unit reference number of file containing filenames
                                                              ! 5.  unit reference number of file containing data
                                                              ! 6.  time coding option for reading time from data file
                                                              ! 7.  =0: time series, =1: nesting
                                                              ! 8.  number of locations for which boundary values are in the file
                                                              ! 9.  not used
                                                              ! 10. not used
                                                              ! 11. not used
                                                              ! 12. not used
                                                              ! 13. ordering of data on file
                                                              ! 14. number of heading lines per file
                                                              ! 15. number of heading lines per time step
                                                              ! 16. number of heading lines per flow variable
                                                              ! 17. not used
                                                              ! 18. =1: Cartesian direction, =2: Nautical direction
                                                              ! 19. not used
                                                              ! 20. not used
       real*8                :: bctime(2)                     ! data containing times of boundary values to be read
                                                              ! 1. time of boundary values read one before last
                                                              ! 2. time of boundary values read last
       integer, pointer      :: bcloc(:)                      ! place in array bndval where to store interpolated boundary value
       type(bfldat), pointer :: nextbfl                       ! pointer to next boundary condition file in list
    end type bfldat

    type(bfldat), save, target :: fbndfil                     ! first boundary condition file in list of files
    logical, save              :: lbfils = .false.            ! indicates whether linked list of boundary condition file is initialized or not

    type bfsdat
       integer                :: nbfs                         ! sequence number
       integer                :: nfreq                        ! number of Fourier components
       real                   :: azero                        ! amplitude at zero frequency in Fourier series
       real, pointer          :: ampl (:)                     ! amplitude of a Fourier component
       real, pointer          :: omega(:)                     ! angular frequency of a Fourier component
       real, pointer          :: phase(:)                     ! phase of a Fourier component
       real, pointer          :: theta(:)                     ! randomly chosen directions assigned to each frequency
       real, pointer          :: sfamp(:)                     ! source function amplitude
       real, pointer          :: bshap(:)                     ! shape factor beta of the source area
       real, dimension(5)     :: spparm                       ! parameters used for computation of incident spectrum and subsequent Fourier series
                                                              ! =1; wave height (significant or rms)
                                                              ! =2; wave period (peak or mean)
                                                              ! =3; incident or peak wave direction w.r.t. normal on the boundary
                                                              ! =4; directional distribution coefficient
                                                              ! =5; cyclic period of time series
       type(bfsdat), pointer  :: nextbfs                      ! pointer to next item in list of Fourier series parameters
    end type bfsdat

    type(bfsdat), save, target  :: fbfs                       ! first item in list of Fourier series parameters
    type(bfsdat), save, pointer :: cubfs                      ! current item in list of Fourier series parameters
    logical, save               :: lbfs   = .false.           ! indicates whether linked list of Fourier series parameters is initialized or not

    type bgpdat
       integer*8              :: bgp(9)                       ! array containing data w.r.t. boundary grid points
                                                              ! =1; grid point index
                                                              ! =2; boundary type
                                                              ! =3; weight factor belonging to start boundary value for interpolation
                                                              ! =4; sequence number of start boundary value
                                                              ! =5; weight factor belonging to end boundary value for interpolation
                                                              ! =6; sequence number of end boundary value
                                                              ! =7; smoothing period
                                                              ! =8; layer number
                                                              ! =9; indicator for addition of bound long wave or linearized Riemann invariant
       type(bgpdat), pointer  :: nextbgp                      ! pointer to next item in list of boundary grid points
    end type bgpdat

    type(bgpdat), save, target  :: fbgp                       ! first item in list of boundary grid points
    type(bgpdat), save, pointer :: cubgp                      ! current item in list of boundary grid points
    logical, save               :: lbgp   = .false.           ! indicates whether linked list of boundary grid points is initialized or not
!
!   Source text
!
    interface delete
      module procedure deletebfl
      module procedure deletebfs
      module procedure deletebgp
    end interface
    !
    contains

    recursive subroutine deletebfl ( bfl )

    type(bfldat) :: bfl

    if ( lbfils ) then
       if ( associated( bfl%bcloc ) ) then
          deallocate( bfl%bcloc )
          nullify( bfl%bcloc )
       endif
       if ( associated( bfl%nextbfl ) ) then
          call deletebfl ( bfl%nextbfl )
          deallocate( bfl%nextbfl )
          nullify( bfl%nextbfl )
       endif
    endif

    end subroutine
    !
    recursive subroutine deletebfs ( bfs )

    type(bfsdat) :: bfs
    integer :: nf

    nf = bfs%nfreq

    if ( lbfs ) then
       if ( nf /= 0 ) then
          if ( associated( bfs%ampl ) ) then
             deallocate( bfs%ampl )
             nullify( bfs%ampl )
          endif
          if ( associated( bfs%omega ) ) then
             deallocate( bfs%omega )
             nullify( bfs%omega )
          endif
          if ( associated( bfs%phase ) ) then
             deallocate( bfs%phase )
             nullify( bfs%phase )
          endif
          if ( associated( bfs%theta ) ) then
             deallocate( bfs%theta )
             nullify( bfs%theta )
          endif
          if ( iwvgen /= 0 .and. associated( bfs%sfamp ) ) then
             deallocate( bfs%sfamp )
             nullify( bfs%sfamp )
          endif
          if ( iwvgen /= 0 .and. associated( bfs%bshap ) ) then
             deallocate( bfs%bshap )
             nullify( bfs%bshap )
          endif
       endif
       if ( associated( bfs%nextbfs ) ) then
          call deletebfs ( bfs%nextbfs )
          deallocate( bfs%nextbfs )
          nullify( bfs%nextbfs )
       endif
    endif

    end subroutine
    !
    recursive subroutine deletebgp ( bgp )

    type(bgpdat) :: bgp

    if ( lbgp ) then
       if ( associated( bgp%nextbgp ) ) then
          call deletebgp ( bgp%nextbgp )
          deallocate( bgp%nextbgp )
          nullify( bgp%nextbgp )
       endif
    endif

    end subroutine
!
end module m_bndspec

module m_genarr
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
!    1.00, February 2010: New Module
!
!   Purpose
!
!   Creates several allocatable arrays for Swash computation
!
!   Method
!
!   MODULE construct
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
    integer, parameter                         :: lwork = 2   ! maximum number of work arrays
!
!   Module variables
!
    integer*8, dimension(:), save, allocatable :: bgridp      ! data concerning boundary grid points
    integer, dimension(:)  , save, allocatable :: indlay      ! layer thickness indicator
                                                              ! =1; thickness in percentage
                                                              ! =2; thickness in meter
    integer, dimension(:,:), save, allocatable :: iwrk        ! work array to temporarily store data
    integer, dimension(:,:), save, allocatable :: kgrpnt      ! index table containing the address of each (active) grid point
                                                              ! =1; not active grid point
                                                              ! >1; active grid point
    integer, dimension(:)  , save, allocatable :: kpv         ! layer index of pressure grid on finer velocity grid
    integer, dimension(:)  , save, allocatable :: kup         ! layer interfaces of pressure grid as counted in velocity grid
    integer, dimension(:)  , save, allocatable :: mend        ! end index of loop over u-points (for 2DH or 3D computation)
    integer, dimension(:)  , save, allocatable :: msta        ! start index of loop over u-points (for 2DH or 3D computation)
    integer, dimension(:)  , save, allocatable :: nend        ! end index of loop over v-points (for 2DH or 3D computation)
    integer, dimension(:)  , save, allocatable :: nsta        ! start index of loop over v-points (for 2DH or 3D computation)
    integer, dimension(:)  , save, allocatable :: npu         ! number of velocity layers per pressure layer
    integer, dimension(:)  , save, allocatable :: vmspon      ! boundary markers where sponge layers are to be applied
    integer, dimension(:)  , save, allocatable :: xytst       ! grid point indices of test points
    !
    real   , dimension(:,:), save, allocatable :: avxf        ! input field of ambient u-velocity mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: avxi        ! user-defined input field of ambient u-velocity
    real   , dimension(:,:), save, allocatable :: avyf        ! input field of ambient v-velocity mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: avyi        ! user-defined input field of ambient v-velocity
    real   , dimension(:,:), save, allocatable :: bndval      ! array containing boundary values at two time levels
    real   , dimension(:)  , save, allocatable :: bveg        ! stem diameter for each layer
    real   , dimension(:)  , save, allocatable :: cdveg       ! drag coefficient due to vegetation for each layer
    real   , dimension(:)  , save, allocatable :: corp        ! user-defined input field of Coriolis parameter
    real   , dimension(:)  , save, allocatable :: corpf       ! input field of Coriolis parameter mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: depf        ! input field of depth w.r.t. still water level mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: depth       ! user-defined input field of bottom level
    real   , dimension(:)  , save, allocatable :: flobj       ! user-defined input field of draft of floating object
    real   , dimension(:)  , save, allocatable :: flobjf      ! input field of draft of floating object mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: fric        ! user-defined input field of friction
    real   , dimension(:,:), save, allocatable :: fricf       ! input field of friction at three time levels mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: gsiz        ! user-defined input field of characterstic grain size of armour rock
    real   , dimension(:)  , save, allocatable :: gsizf       ! input field of characterstic grain size of armour rock mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: guu         ! square root of g_22 in u-point
    real   , dimension(:)  , save, allocatable :: guv         ! square root of g_22 in v-point
    real   , dimension(:)  , save, allocatable :: gvu         ! square root of g_11 in u-point
    real   , dimension(:)  , save, allocatable :: gvv         ! square root of g_11 in v-point
    real   , dimension(:)  , save, allocatable :: gsqs        ! Jacobian / volume cell in wl-point
    real   , dimension(:)  , save, allocatable :: gsqsu       ! Jacobian in u-point
    real   , dimension(:)  , save, allocatable :: gsqsv       ! Jacobian in v-point
    real   , dimension(:)  , save, allocatable :: hlay        ! layer thicknesses
    real   , dimension(:)  , save, allocatable :: hlaysh      ! layer thicknesses when a fixed layer exist (shadow sigma transformation)
    real   , dimension(:)  , save, allocatable :: hlayv       ! layer thickness for vegetation model
    real   , dimension(:)  , save, allocatable :: hpor        ! user-defined input field of porous structure height
    real   , dimension(:)  , save, allocatable :: hporf       ! input field of porous structure height mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: lbod        ! user-defined input field of label of rigid body
    real   , dimension(:)  , save, allocatable :: lbodf       ! input field of label of rigid body mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: mwl         ! user-defined input field of mean water level (tide/wind surge)
    real   , dimension(:)  , save, allocatable :: npla        ! user-defined input field of number of plants per square meter
    real   , dimension(:)  , save, allocatable :: nplaf       ! input field of number of plants per square meter mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: npor        ! user-defined input field of porosity
    real   , dimension(:)  , save, allocatable :: nporf       ! input field of porosity mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: nveg        ! number of plants / m2 for each layer
    real   , dimension(:)  , save, allocatable :: pres        ! user-defined input field of atmospheric pressure
    real   , dimension(:,:), save, allocatable :: presf       ! input field of atmospheric pressure at three time levels mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: sal         ! user-defined input field of salinity
    real   , dimension(:,:), save, allocatable :: salf        ! input field of salinity mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: sed         ! user-defined input field of suspended sediment
    real   , dimension(:,:), save, allocatable :: sedf        ! input field of suspended sediment mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: temp        ! user-defined input field of temperature
    real   , dimension(:,:), save, allocatable :: tempf       ! input field of temperature mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: uwrk        ! work array to temporarily store data (unstructured grids)
    real   , dimension(:)  , save, allocatable :: uxb         ! user-defined input field of u-velocity
    real   , dimension(:,:), save, allocatable :: uxf         ! input field of u-velocity at three time levels mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: uyb         ! user-defined input field of v-velocity
    real   , dimension(:,:), save, allocatable :: uyf         ! input field of v-velocity at three time levels mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: wdspon      ! width of sponge layers
    real   , dimension(:,:), save, allocatable :: wlevf       ! input field of water level at three time levels mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: wlevl       ! user-defined input field of water level
    real   , dimension(:,:), save, allocatable :: work        ! work array to temporarily store data
    real   , dimension(:,:), save, allocatable :: wrk         ! work array to temporarily store 3D data
    real   , dimension(:,:), save, allocatable :: wxf         ! input field of wind u-velocity at three time levels mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: wxi         ! user-defined input field of wind u-velocity
    real   , dimension(:,:), save, allocatable :: wyf         ! input field of wind v-velocity at three time levels mapped onto computational grid
    real   , dimension(:)  , save, allocatable :: wyi         ! user-defined input field of wind v-velocity
    real   , dimension(:,:), save, allocatable :: xcgrid      ! coordinates of computational grid in x-direction
    real   , dimension(:,:), save, allocatable :: ycgrid      ! coordinates of computational grid in y-direction
    !
    logical, dimension(:)  , save, allocatable :: sbimp       ! water level prescribed at lower boundary (for 2DH or 3D computation)
    logical, dimension(:)  , save, allocatable :: slimp       ! water level prescribed at left boundary (for 2DH or 3D computation)
    logical, dimension(:)  , save, allocatable :: srimp       ! water level prescribed at right boundary (for 2DH or 3D computation)
    logical, dimension(:)  , save, allocatable :: stimp       ! water level prescribed at upper boundary (for 2DH or 3D computation)
    logical, dimension(:)  , save, allocatable :: wlimp       ! water level prescribed at boundary face (for unstructured mesh)
!
!   Source text
!
end module m_genarr
