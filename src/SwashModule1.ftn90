!                 COMMON VARIABLES RELATED MODULES, file 1 of 2
!
!     Contents of this file
!
!     SwashCommdata1     contains common variables for Swash
!     SwashCommdata2     contains common variables for Swash
!     SwashCommdata3     contains common variables for Swash
!     SwashCommdata4     contains common variables for Swash
!     SwashTimecomm      contains common variables for Swash
!
module SwashCommdata1
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
!    1.00, February 2010: New Module based on module SWCOMM1
!
!   Purpose
!
!   Module containing data of units and output
!
!   Method
!
!   Data with respect to units and output
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
    integer, parameter                         :: nmovar = 120 ! maximum number of output variables
    integer, parameter                         :: moutpa =   5 ! number of output parameters
!
!   Module variables
!
    ! units
    !
    character(6)                               :: uc          ! unit of salinity
    character(6)                               :: ud          ! unit of density/concentration
    character(6)                               :: udi         ! unit of direction
    character(6)                               :: ue          ! unit of (eddy) viscosity
    character(6)                               :: uf          ! unit of stress
    character(6)                               :: uh          ! unit of vertical length
    character(6)                               :: uk          ! unit of temperature
    character(6)                               :: ul          ! unit of horizontal length
    character(6)                               :: up          ! unit of pressure
    character(6)                               :: uq          ! unit of discharge
    character(6)                               :: ut          ! unit of time
    character(6)                               :: uv          ! unit of velocity
    !
    ! information for output
    !
    integer                                    :: errpts      ! unit reference number of file containing coordinates of problem points
                                                              ! =16, unit reference number of the file
    integer, dimension(nmovar)                 :: ovsvty      ! type of the output variable
                                                              ! =1; scalar
                                                              ! =2; angle
                                                              ! =3; vector
                                                              ! =4; multiple
    !
    real                                       :: alcq        ! angle between x-axis of computational grid and output frame
    real                                       :: alpq        ! angle between x-axis of user coordinate system and output frame
    real                                       :: coscq       ! cos of alcq
    real                                       :: cospq       ! cos of alpq
    real                                       :: dxk         ! mesh size of output frame
    real                                       :: dyk         ! mesh size of output frame
    real, dimension(moutpa)                    :: outpar      ! contains parameters meant for outputting
    real, dimension(nmovar)                    :: ovexcv      ! exception value for output quantity
    real, dimension(nmovar)                    :: ovlexp      ! lower expected limit of output quantity
    real, dimension(nmovar)                    :: ovllim      ! lower limit of validity of output quantity
    real, dimension(nmovar)                    :: ovhexp      ! upper expected limit of output quantity
    real, dimension(nmovar)                    :: ovulim      ! upper limit of validity
    real                                       :: sincq       ! sin of alcq
    real                                       :: sinpq       ! sin of alpq
    real                                       :: xpq         ! x-coordinate of origin of output frame w.r.t. computational coordinates
    real                                       :: xqlen       ! length of x-side of output frame
    real                                       :: xqp         ! x-coordinate (user coordinates) of origin of output frame
    real                                       :: ypq         ! y-coordinate of origin of output frame w.r.t. computational coordinates
    real                                       :: yqlen       ! length of y-side of output frame
    real                                       :: yqp         ! y-coordinate (user coordinates) of origin of output frame
    !
    character(8),  dimension(nmovar)           :: ovkeyw      ! keyword identifying output quantity in a Swash command
    character(60), dimension(nmovar)           :: ovlnam      ! long name of output quantity
    character(6),  dimension(nmovar)           :: ovsnam      ! short name of output quantity
    character(16), dimension(nmovar)           :: ovunit      ! unit of of output quantity
    character(8)                               :: sname       ! name of output point set
!
!   Source text
!
end module SwashCommdata1
!
module SwashCommdata2
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
!    1.00, February 2010: New Module based on module SWCOMM2
!
!   Purpose
!
!   Module containing data of input grids
!
!   Method
!
!   Data with respect to input grids
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
    integer, parameter                         :: numgrd = 23 ! maximum number of input grids
    !
    real   , parameter                         :: NEAREST = -12345678. ! quantity outside input grid equals value nearest boundary of input grid
!
!   Module variables
!
    ! locations and dimensions
    !
    integer, dimension(numgrd)                 :: igtype      ! =0; not used
                                                              ! =1; rectilinear grid
                                                              ! =2; curvilinear grid
                                                              ! =3; unstructured grid
    integer, dimension(numgrd)                 :: leds        ! =0; values have not been read
                                                              ! =1; indicates INP command has been given
                                                              ! =2; indicates field is read
                                                              ! =3; indicates bottom and current are read
    integer, dimension(numgrd)                 :: mxg         ! number of grid points of input grid in x-direction
    integer, dimension(numgrd)                 :: myg         ! number of grid points of input grid in y-direction
    !
    real, dimension(numgrd)                    :: alpg        ! direction of the x-axis w.r.t. user coordinates
    real, dimension(numgrd)                    :: cospg       ! cos of alpg
    real                                       :: cosvc       ! =cos(-alpg(2)),
                                                              ! cos of angle of flow velocity input grid w.r.t. computational grid
    real                                       :: coswc       ! =cos(-alpg(5)),
                                                              ! cos of angle of wind velocity input grid w.r.t. computational grid
    real, dimension(numgrd)                    :: dxg         ! mesh size of input grid in x-direction
    real, dimension(numgrd)                    :: dyg         ! mesh size of input grid in y-direction
    real, dimension(numgrd)                    :: excfld      ! exception values for input grids
    real, dimension(numgrd)                    :: sinpg       ! sin of alpg
    real                                       :: sinvc       ! =sin(-alpg(2)),
                                                              ! sin of angle of flow velocity input grid w.r.t. computational grid
    real                                       :: sinwc       ! =sin(-alpg(5)),
                                                              ! sin of angle of wind velocity input grid w.r.t. computational grid
    real, dimension(numgrd)                    :: xpg         ! x-coordinate of origin of input grid
    real, dimension(numgrd)                    :: ypg         ! y-coordinate of origin of input grid
    !
    logical                                    :: initsf      ! indicates water level to be initialized based on input field
    logical                                    :: inituf      ! indicates u-velocity to be initialized based on input field
    logical                                    :: initvf      ! indicates v-velocity to be initialized based on input field
    logical                                    :: varea       ! mean water level associated with ambient flow is/is not variable over space
    logical                                    :: varfr       ! friction coefficient is/is not variable over space
    logical                                    :: vargs       ! grain size of porous structure is/is not variable over space
    logical                                    :: varnpl      ! number of plants / m2 is/is not variable over space
    logical                                    :: varsh       ! structure height is/is not variable over space
    logical                                    :: varva       ! ambient current is/is not variable over space
    logical                                    :: varwi       ! wind velocity is/is not variable over space
    logical, dimension(numgrd)                 :: lflgrd      ! indicates whether input grid equals computational grid
    logical, dimension(numgrd)                 :: lstag       ! staggering of curvilinear input grid w.r.t. computational grid
    !
    ! variables for input files
    !
    integer, dimension(numgrd)                 :: ifldyn      ! if =0: data is stationary, if =1: nonstationary
    integer, dimension(numgrd)                 :: iflidl      ! lay-out in input file
    integer, dimension(numgrd)                 :: iflifm      ! format identifier
    integer, dimension(numgrd)                 :: ifllay      ! number of layers for each input field
    integer, dimension(numgrd)                 :: iflndf      ! unit reference number of namelist file
    integer, dimension(numgrd)                 :: iflnds      ! unit reference number of data file
    integer, dimension(numgrd)                 :: iflnhd      ! number of heading lines per input field
    integer, dimension(numgrd)                 :: iflnhf      ! number of heading lines per file
    !
    real*8, dimension(numgrd)                  :: iflbeg      ! begin time of data on file
    real*8, dimension(numgrd)                  :: iflend      ! end time of data on file
    real  , dimension(numgrd)                  :: iflfac      ! multiplication factor
    real*8, dimension(numgrd)                  :: iflint      ! time interval of data on file
    real*8, dimension(numgrd)                  :: ifltim      ! time of last reading
    !
    character(40), dimension(numgrd)           :: iflfrm      ! format string
!
!   Source text
!
end module SwashCommdata2
!
module SwashCommdata3
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
!    1.00, February 2010: New Module based on module SWCOMM3
!
!   Purpose
!
!   Module containing data of physical and numerical parameters
!
!   Method
!
!   Data with respect to physical and numerical parameters
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
    integer, parameter                         :: mbot  =  5  ! dimension of array pbot
    integer, parameter                         :: mcor  =  5  ! dimension of array pcor
    integer, parameter                         :: miwg  = 10  ! dimension of array piwg
    integer, parameter                         :: mnums = 75  ! dimension of array pnums
    integer, parameter                         :: mpor  =  5  ! dimension of array ppor
    integer, parameter                         :: msed  = 15  ! dimension of array psed
    integer, parameter                         :: mship = 15  ! dimension of array pship
    integer, parameter                         :: msurf =  5  ! dimension of array psurf
    integer, parameter                         :: mturb =  5  ! dimension of array pturb
    integer, parameter                         :: mwnd  = 10  ! dimension of array pwnd
    !
    integer, parameter                         :: mtamx = 3   ! assumed maximum number of time-step level in local time stepping
                                                              ! (to be used in unstructured-mesh transport module)
!
!   Module variables
!
    ! location and dimensions of computational grid
    !
    integer                                    :: mcgrd       ! number of active grid points of the computational grid
    integer                                    :: mtc         ! number of time steps
    integer                                    :: mxc         ! number of cells in x-direction (including virtual ones)
    integer                                    :: myc         ! number of cells in y-direction (including virtual ones)
    integer                                    :: optg        ! type of computational grid
                                                              ! =1; rectilinear
                                                              ! =3; curvilinear
                                                              ! =5; unstructured grid
    !
    real                                       :: alcp        ! =-alpc;
                                                              ! direction of user coordinates w.r.t. computational coordinates
    real                                       :: alpc        ! direction of x-axis of computational grid w.r.t. user coordinates
    real                                       :: cospc       ! =cos(alpc)
    real                                       :: dx          ! =xclen/mxs; mesh size in x-direction of computational grid
    real                                       :: dy          ! =yclen/mys; mesh size in y-direction of computational grid
    real                                       :: sinpc       ! =sin(alpc)
    real                                       :: xcgmax      ! maximum x-coordinate of computational grid points
    real                                       :: xcgmin      ! minimum x-coordinate of computational grid points
    real                                       :: xclen       ! length of computational grid in x-direction
    real                                       :: xcp         ! =-xpc*cospc-ypc*sinpc;
                                                              ! origin of user coordinates w.r.t. computational coordinates
    real                                       :: xoffs       ! offset value in x-coordinate
                                                              ! (from user coordinate system to internal coordinate system)
    real                                       :: xpc         ! x-coordinate of origin of computational grid
    real                                       :: ycgmax      ! maximum y-coordinate of computational grid points
    real                                       :: ycgmin      ! minimum y-coordinate of computational grid points
    real                                       :: yclen       ! length of computational grid in y-direction
    real                                       :: ycp         ! =xpc*sinpc-ypc*cospc;
                                                              ! origin of user coordinates w.r.t. computational coordinates
    real                                       :: ylatc       ! geographic latitude coordinate of model area
    real                                       :: yoffs       ! offset value in y-coordinate
                                                              ! (from user coordinate system to internal coordinate system)
    real                                       :: ypc         ! y-coordinate of origin of computational grid
    !
    logical                                    :: cvleft      ! the computational grid is left/right-oriented
                                                              ! =.true. ; left-handed
                                                              ! =.false.; right-handed
    logical                                    :: lxoffs      ! offset values were/were not initialized already
    !
    ! vertical layer schematization
    !
    integer                                    :: kmax        ! number of vertical (velocity) layers
    integer                                    :: kpmax       ! number of vertical pressure layers
    !
    real                                       :: depfix      ! sum of the thicknesses of fixed layers (if any)
    real                                       :: epsfix      ! threshold depth below which the shadow sigma transformation is employed
    !
    logical                                    :: fixlay      ! true if there is at least one layer with fixed thickness
    logical                                    :: lsubg       ! subgrid approach is applied
    !
    ! information for initial conditions
    !
    real                                       :: tkeini      ! initial constant value for turbulent kinetic energy
    real                                       :: epsini      ! initial constant value for dissipation rate of turbulent kinetic energy
    !
    logical                                    :: instead     ! indicates whether initial flow velocities are based on steady flow condition
    logical                                    :: restrt      ! indicates whether the run starts cold or hot
    !
    ! information for boundary conditions
    !
    integer                                    :: nbfils      ! number of boundary condition files
    integer                                    :: nbval       ! number of boundary values
    integer                                    :: nbgrpt      ! number of boundary grid points
    integer                                    :: nspl        ! number of sponge layers (unstructured grids)
    !
    real                                       :: spwidb      ! width of sponge layer from lower boundary in y-direction
    real                                       :: spwidl      ! width of sponge layer from left boundary in x-direction
    real                                       :: spwidr      ! width of sponge layer from right boundary in x-direction
    real                                       :: spwidt      ! width of sponge layer from upper boundary in y-direction
    real                                       :: tcret       ! constituent return time to be used in Thatcher-Harleman's model of unsteady salt intrusion
    !
    logical                                    :: lrampf      ! ramp function is/is not applied
    !
    ! internal wave generation
    !
    integer                                    :: iwvgen      ! indicator internal wave generation
                                                              ! =0; no addition mass source to continuity equation
                                                              ! =1; add mass source to continuity equation
    logical                                    :: lsrcfx      ! indicates whether the source function is parallel to x-axis or not
    logical                                    :: lsrcfy      ! indicates whether the source function is parallel to y-axis or not
    !
    ! constants
    !
    integer                                    :: myseed      ! =12345678; constant to seed the random number generator
                                                              ! set by command SET ... [seed] ...
    !
    real                                       :: degrad      ! =pi/180; constant to transform degrees to radians
    real                                       :: dnorth      ! direction of the North w.r.t. x-axis of user coordinates
                                                              ! =nor; set by command SET ... [nor] ...
    real                                       :: pi          ! circular constant
    real                                       :: pi2         ! 2*pi
    real                                       :: swl         ! still water level
                                                              ! =swl; set by command SET ... [level] ...
    !
    ! physical parameters
    !
    real                                       :: dynvis      ! dynamic viscosity of water
                                                              ! =dynvis; set by command SET ... [dynvis] ...
    real                                       :: grav        ! acceleration due to gravity
                                                              ! =grav; set by command SET ... [grav] ...
    real                                       :: kinvis      ! kinematic viscosity of water
    real                                       :: rhoa        ! density of air
                                                              ! =rhoa; set by command SET ... [rhoair] ...
    real                                       :: rhow        ! (reference) density of water
                                                              ! =rhow; set by command SET ... [rhowat] ...
    real                                       :: vonkar      ! the Von Karman constant
                                                              ! =vonkar; set by command SET ... [kappa] ...
    !
    ! density parameters
    !
    integer                                    :: idens       ! indicator baroclinic forcing
                                                              ! =0; no baroclinic forcing
                                                              ! =1; include baroclinic forcing
    real                                       :: rhos        ! density of sediment
                                                              ! =rhos; set by command SET ... [rhosed] ...
    real                                       :: salw        ! (ambient) salinity of water
                                                              ! =salw; set by command SET ... [salinity] ...
    real                                       :: tempw       ! (ambient) temperature of water
                                                              ! =tempw; set by command SET ... [temp] ...
    logical                                    :: lmixt       ! mixture density due to sediment is/is not included
    !
    ! viscosity and diffusivity parameters
    !
    real                                       :: bvisc       ! background viscosity
                                                              ! =bvisc; set by command SET ... [backvisc] ...
    real                                       :: csmag       ! the Smagorinsky constant
    real                                       :: hdiff       ! horizontal diffusivity
    real                                       :: hvisc       ! horizontal viscosity
    real                                       :: lmix        ! mixing length
    !
    ! parameters for ambient current (riverine, tidal and wind-driven flows)
    !
    integer                                    :: iamb        ! indicator ambient current
                                                              ! =0; no ambient current
                                                              ! =1; cell-centered approach
                                                              ! =2; staggered approach
    real                                       :: eamb        ! mean water level associated with ambient flow
    real                                       :: uamb        ! ambient U-velocity
    real                                       :: vamb        ! ambient V-velocity
    !
    ! parameters for floating objects
    !
    integer                                    :: ifloat      ! indicator floating objects
                                                              ! =0; no floating objects
                                                              ! =1; include floating objects
    !
    ! parameters for physics
    !
    integer                                    :: ihvisc      ! indicator horizontal eddy viscosity
                                                              ! =0; no horizontal eddy viscosity
                                                              ! =1; constant eddy viscosity
                                                              ! =2; Smagorinsky model
                                                              ! =3; mixing length model
                                                              ! =4; water depth based mixing length model for wave breaking
                                                              !     (not used if isurf=2)
    integer                                    :: iporos      ! indicator porosity
                                                              ! =0; no porosity
                                                              ! =1; porosity included
    integer                                    :: irough      ! indicator roughness method
                                                              ! =0; no bottom friction
                                                              ! =1; constant friction parameter
                                                              ! =2; Chezy roughness method
                                                              ! =3; Manning roughness method
                                                              ! =4; logarithmic wall-law
                                                              ! =5; Colebrook-White roughness method
                                                              ! =11; linear bottom friction
    integer                                    :: isurf       ! indicator control of wave breaking
                                                              ! =0; no control of wave breaking
                                                              ! =1; controls wave breaking using bore formation concept (2nd version)
                                                              ! =2; controls wave breaking using bore formation concept (1st version)
    integer                                    :: itrans      ! indicator transport of constituent
                                                              ! =0; no constituent
                                                              ! =1; transport of constituent
    integer                                    :: iturb       ! indicator turbulence modelling
                                                              ! =0; no turbulence
                                                              ! =1; standard k-eps model (vertical)
                                                              ! =2: linear k-eps model (3D)
                                                              ! =3; nonlinear k-eps model (3D)
    integer                                    :: iveg        ! indicator vegetation with force in horizontal direction
                                                              ! =0; no vegetation
                                                              ! =1; vegetation friction included
                                                              ! =2; vegetation friction and porosity included
    integer                                    :: ivegw       ! indicator vegetation with force in vertical direction
                                                              ! =0; no vegetation
                                                              ! =1; vegetation friction included
    integer                                    :: iwind       ! indicator wind stress
                                                              ! =0; no wind stress
                                                              ! =1; constant drag coefficient
                                                              ! =2; drag coefficient based on Charnock
                                                              ! =3; drag coefficient linear dependent on wind speed
                                                              ! =4; drag coefficient based on 2nd order polynomial fit
    !
    integer                                    :: lmax        ! number of layers to be used in vegetation model
    integer                                    :: lsal        ! constituent number for salinity
    integer                                    :: lsed        ! constituent number for suspended sediment
    integer                                    :: ltemp       ! constituent number for temperature
    integer                                    :: ltrans      ! number of transport constituents employed
    integer                                    :: ltur        ! number of turbulence quantities employed
    !
    real, dimension(mbot)                      :: pbot        ! parameters for the bottom friction coefficient
                                                              ! =1; friction parameter (constant, Chezy, Manning or Colebrook-White)
                                                              ! =2; Nikuradse roughness height
    real, dimension(mcor)                      :: pcor        ! parameters for the Coriolis force
                                                              ! =1; constant Coriolis parameter
                                                              ! =2; (small) parameter to modify the weights of AB2 scheme
    real, dimension(miwg)                      :: piwg        ! parameters for internal wave generation
                                                              ! =1; centre of gravity of the source area
                                                              ! =2; width of the source area
                                                              ! =3; still water depth at the source area
                                                              ! =4; delta shape factor of the source area
                                                              ! =5; period for smoothing source function during cold start
                                                              ! =6;  boundary marker (unstructured grids)
    real, dimension(mpor)                      :: ppor        ! parameters for porous structures
                                                              ! =1; characteristic grain size
                                                              ! =2; structure height
                                                              ! =3; dimensionless constant for laminar friction factor (=alpha0)
                                                              ! =4; dimensionless constant for turbulent friction factor (=beta0)
                                                              ! =5; typical wave period
    real, dimension(msed)                      :: psed        ! parameters for suspended sediment transport
                                                              ! =1; fall velocity for suspended sediment
                                                              ! =2; median sediment diameter
                                                              ! =3; Schmidt number for sediment
                                                              ! =4; alfa = ratio between roughness height and particle size (=5.5)
                                                              ! =5; critical Shields parameter (=0.05) (noncohesive sediment)
                                                              ! =6; dimensionless parameter for sediment pickup (=0.00033)
                                                              ! =7; ek = coefficient enhancing pickup under breaking waves due to breaking-induced turbulence
                                                              ! =8; factor = (s-1) * grav * d50
                                                              ! =9; proportionality factor of pickup function
                                                              ! =10; critical bed shear stress for erosion (cohesive sediment)
                                                              ! =11; critical bed shear stress for deposition (cohesive sediment)
                                                              ! =12; entrainment rate for erosion flux (cohesive sediment)
                                                              ! =13; empirical constant used in stability factor for sediment-laden BBL (=5.5)
    real, dimension(mship)                     :: pship       ! parameters for floating rigid bodies
                                                              ! = 1; compressibility factor for flow beneath fixed floating object
                                                              ! = 2; theta; implicitness factor for time integration used beneath floating object
                                                              !      =0.5: Crank-Nicolson
                                                              !      =1.0: implicit Euler
                                                              ! = 3; beta; implicitness factor for Newmark scheme
                                                              ! = 4; gamma; implicitness factor for Newmark scheme
                                                              ! = 5; rhoinf; spectral radius (for infinite frequency) to be used in generalized-alpha schemes
                                                              ! = 6; not used
                                                              ! = 7; not used
                                                              ! = 8; maximum error for the fluid-structure coupling
                                                              ! = 9; maximum number of iterations to be performed in the fluid-structure coupling
                                                              ! =10; under-relaxation for rigid body accelerations (added-mass effect)
                                                              ! =11; theta; implicitness factor for update kinematic boundary conditions
                                                              !      =0.5: Crank-Nicolson
                                                              !      =1.0: implicit Euler
    real, dimension(msurf)                     :: psurf       ! parameters for control of wave breaking
                                                              ! =1; threshold parameter corresponding to onset of wave breaking
                                                              ! =2; threshold parameter corresponding to re-initiation of breaking in post-breaking area
                                                              !     (not used if isurf=2)
                                                              ! =3; proportionality factor for water depth based mixing length model
                                                              !     (not used if isurf=2)
    real, dimension(mturb)                     :: pturb       ! parameters for turbulence modelling
                                                              ! =1; not used
                                                              ! =2; vegetative drag-related closure constant for turbulent kinetic energy
                                                              ! =3; vegetative drag-related closure constant for dissipation rate
    real, dimension(mwnd)                      :: pwnd        ! parameters for the wind stress coefficient
                                                              ! =1; constant wind stress coefficient
                                                              ! =2; beta for Charnock formula
                                                              ! =3; height
                                                              ! =4; correction factor for wind at 10 metres height to wind at surface
                                                              ! =5; first coefficient of linear drag formula
                                                              ! =6; second coefficient of linear drag formula
                                                              ! =7; third coefficient of linear drag formula
                                                              ! =8; minimum wind speed to apply linear dependence of drag to speed
                                                              ! =9; maximum wind speed to apply linear dependence of drag to speed
                                                              ! =10; ratio of forced crest height to maximum surface elevation
    !
    real                                       :: cdcap       ! maximum wind drag coefficient
                                                              ! =cdcap; set by command SET ... [cdcap] ...
    real                                       :: cvm         ! virtual/added mass coefficient in inertia force due to vegetation
    real                                       :: prmean      ! mean atmospheric pressure
                                                              ! =prmean; set by command SET ... [prmean] ...
    real                                       :: u10         ! wind velocity
    real                                       :: wdic        ! =pi2*((wdip/pi2-nint(wdip/pi2))
    real                                       :: wdip        ! wind direction with respect to problem coordinates
    !
    logical                                    :: relwav      ! indicates whether wind stress is based on wind velocity relative to wave celerity or not
    logical                                    :: relwnd      ! indicates whether wind stress is based on wind velocity relative to flow velocity or not
    !
    ! numerical parameters
    !
    integer                                    :: dpsopt      ! method of determining the bottom levels in cell centers
                                                              ! =1; DPS is the minimum of the four surrounding bottom values at cell corners
                                                              ! =2; DPS is the mean of the four surrounding bottom values at cell corners
                                                              ! =3; DPS is the maximum of the four surrounding bottom values at cell corners
                                                              ! =4; bottom level at upper-right corner is shifted to cell center
    integer                                    :: icreep      ! indicates method of computing anti-creep terms in transport equation
                                                              ! =0; no anti-creep
                                                              ! =1; standard method as described in Technical documentation WAQUA / TRIWAQ
                                                              ! =2; method of Stelling and Van Kester (1994)
    integer                                    :: ihydro      ! indicates (non-)hydrostatic flow computation
                                                              ! =0; hydrostatic
                                                              ! =1; non-hydrostatic using box scheme for vertical pressure gradient
                                                              ! =2; non-hydrostatic using standard discretization for vertical pressure gradient
                                                              ! =3; box scheme applied for a pressure layer containing a number of velocity layers
    integer                                    :: inewt       ! method of Newton-type iteration
                                                              ! =0; no Newton
                                                              ! =1; nested Newton of Brugnano and Casulli (2008)
    integer                                    :: iproj       ! indicates type of projection method
                                                              ! =1; pressure correction method
                                                              ! =2; classical projection method
    integer                                    :: jhk         ! starting index of the time integration loop for global continuity equation
                                                              ! =1; second order RK3 scheme
                                                              ! =2; second order Hancock scheme
                                                              ! =3; first order Euler scheme
    integer                                    :: mtimei      ! method of time integration
                                                              ! =1; explicit method
                                                              != 2; implicit method
    integer                                    :: ncompt      ! number of COMPUTE commands
    integer                                    :: nstatc      ! indicates stationarity of computation:
                                                              ! =0; stationary computation
                                                              ! =1; nonstationary computation
    integer                                    :: nstatm      ! = 0; stationary mode
                                                              ! = 1; nonstationary mode
                                                              ! =-1; unknown
    integer                                    :: propsc      ! flag for indicating advection approximation
                                                              ! 1=first order upwind scheme
                                                              ! 2=not used
                                                              ! 3=linear higher order schemes (kappa-scheme)
                                                              ! 4=Sweby Phi-limiter
                                                              ! 5=R-kappa limiter
                                                              ! 6=PL-kappa limiter
    !
    integer                                    :: qlay        ! number of layers for which pressure is constant to be used in reduced pressure equation method
    integer                                    :: qmax        ! total number of layers for reduced pressure Poisson equation
                                                              ! (=kmax-qlay)
    !
    real                                       :: cflmax      ! maximum CFL value found in flow computation
    real                                       :: epsdry      ! minimal depth for drying/wetting
                                                              ! =epsdry; set by command SET ... [depmin] ...
    real                                       :: epshu       ! minimal depth at velocity face
                                                              ! =epshu; set by command SET ... [epshu] ...
    real                                       :: epsuf       ! flow velocity with absolute value smaller than epsuf is zero
                                                              ! to be used for determination of upwind depth
    real                                       :: kappa       ! controls the spatial accuracy of the higher order schemes
    real                                       :: mbound      ! gives the maximum bound of the flux limiter
    real                                       :: phieby      ! indicates a parameter of the Sweby Phi-limiter
    real, dimension(mnums)                     :: pnums       ! numerical parameters
                                                              ! = 1; theta; implicitness factor for time integration of continuity equation
                                                              !      =0.5: Crank-Nicolson
                                                              !      =1.0: implicit Euler
                                                              ! = 2; cfllow; lowest CFL number for determining time step
                                                              ! = 3; cflhig; highest CFL number for determining time step
                                                              ! = 4; theta; implicitness factor for time integration of water level gradient
                                                              !      =0.5: Crank-Nicolson
                                                              !      =1.0: implicit Euler
                                                              ! = 5; theta; implicitness factor for time integration of non-hydrostatic pressure gradient
                                                              !      =0.5: Crank-Nicolson
                                                              !      =1.0: implicit Euler
                                                              ! = 6; propsc; discretization scheme for horizontal advection in u-momentum equation
                                                              ! = 7; kappa; parameter for kappa scheme i.c. horizontal advection in u-momentum equation
                                                              ! = 8; mbound; upper bound for PL-kappa scheme i.c. horizontal advection in u-momentum equation
                                                              ! = 9; phieby; parameter for Sweby-phi limiter i.c. horizontal advection in u-momentum equation
                                                              ! =10; not used
                                                              ! =11; propsc; discretization scheme for water depth correction
                                                              ! =12; kappa; parameter for kappa scheme i.c. water depth correction
                                                              ! =13; mbound; upper bound for PL-kappa scheme i.c. water depth correction
                                                              ! =14; phieby; parameter for Sweby-phi limiter i.c. water depth correction
                                                              ! =15; not used
                                                              ! =16; propsc; discretization scheme for horizontal advection in w-momentum equation
                                                              ! =17; kappa; parameter for kappa scheme i.c. horizontal advection in w-momentum equation
                                                              ! =18; mbound; upper bound for PL-kappa scheme i.c. horizontal advection in w-momentum equation
                                                              ! =19; phieby; parameter for Sweby-phi limiter i.c. horizontal advection in w-momentum equation
                                                              ! =20; not used
                                                              ! =21; relative accuracy with respect to the initial residual for the CG method
                                                              ! =22; relative accuracy with respect to the right-hand side for the SIP solver or the BiCGSTAB method
                                                              ! =23; relative accuracy with respect to the initial residual for the BiCGSTAB method
                                                              ! =24; maximum number of iterations to be performed in the CG method
                                                              ! =25; maximum number of iterations to be performed in the SIP solver or the BiCGSTAB method
                                                              ! =26; weight factor for convex combination of ILUD and MILUD for the CG method
                                                              ! =27; relaxation parameter for the SIP solver or weight factor for convex combination
                                                              !      of ILU(D) and MILU(D) for the BiCGSTAB method
                                                              ! =28; frequency with which preconditioner is applied
                                                              ! =31; thetau; implicitness factor for time integration of the vertical terms in u-momentum equation
                                                              !      =0.5: Crank-Nicolson
                                                              !      =1.0: implicit Euler
                                                              ! =32; thetaw; implicitness factor for time integration of the vertical terms in w-momentum equation
                                                              !      =0.5: Crank-Nicolson
                                                              !      =1.0: implicit Euler
                                                              ! =33; thetat; implicitness factor for time integration of the vertical terms in transport equation
                                                              !      =0.5: Crank-Nicolson
                                                              !      =1.0: implicit Euler
                                                              ! =34; not used
                                                              ! =35; not used
                                                              ! =36; propsc; discretization scheme for vertical advection in u-momentum equation
                                                              ! =37; kappa; parameter for kappa scheme i.c. vertical advection in u-momentum equation
                                                              ! =38; mbound; upper bound for PL-kappa scheme i.c. vertical advection in u-momentum equation
                                                              ! =39; phieby; parameter for Sweby-phi limiter i.c. vertical advection in u-momentum equation
                                                              ! =40; not used
                                                              ! =41; propsc; discretization scheme for vertical advection in w-momentum equation
                                                              ! =42; kappa; parameter for kappa scheme i.c. vertical advection in w-momentum equation
                                                              ! =43; mbound; upper bound for PL-kappa scheme i.c. vertical advection in w-momentum equation
                                                              ! =44; phieby; parameter for Sweby-phi limiter i.c. vertical advection in w-momentum equation
                                                              ! =45; not used
                                                              ! =46; propsc; discretization scheme for horizontal advection in transport equation
                                                              ! =47; kappa; parameter for kappa scheme i.c. horizontal advection in transport equation
                                                              ! =48; mbound; upper bound for PL-kappa scheme i.c. horizontal advection in transport equation
                                                              ! =49; phieby; parameter for Sweby-phi limiter i.c. horizontal advection in transport equation
                                                              ! =50; not used
                                                              ! =51; propsc; discretization scheme for vertical advection in transport equation
                                                              ! =52; kappa; parameter for kappa scheme i.c. vertical advection in transport equation
                                                              ! =53; mbound; upper bound for PL-kappa scheme i.c. vertical advection in transport equation
                                                              ! =54; phieby; parameter for Sweby-phi limiter i.c. vertical advection in transport equation
                                                              ! =55; not used
                                                              ! =58; relative accuracy with respect to water level in the pressure projection method
                                                              ! =59; maximum number of iterations to be performed in the pressure projection method
                                                              ! =61; propsc; discretization scheme for horizontal U-advection with ambient current
                                                              ! =62; kappa; parameter for kappa scheme i.c. horizontal U-advection with ambient current
                                                              ! =63; mbound; upper bound for PL-kappa scheme i.c. horizontal U-advection with ambient current
                                                              ! =64; phieby; parameter for Sweby-phi limiter i.c. horizontal U-advection with ambient current
                                                              ! =65; not used
                                                              ! =66; propsc; discretization scheme for horizontal W-advection with ambient current
                                                              ! =67; kappa; parameter for kappa scheme i.c. horizontal W-advection with ambient current
                                                              ! =68; mbound; upper bound for PL-kappa scheme i.c. horizontal W-advection with ambient current
                                                              ! =69; phieby; parameter for Sweby-phi limiter i.c. horizontal W-advection with ambient current
                                                              ! =70; not used
                                                              ! =71; thetamb; implicitness factor for time integration related to ambient current in global continuity eq
                                                              !      =0.5: Crank-Nicolson
                                                              !      =1.0: implicit Euler
    real*8, dimension(300,5)                   :: rcompt      !
    !
    logical                                    :: bnaut       ! indicates whether nautical or Cartesian directions are used
    logical                                    :: ccn         ! indicates whether non-Delaunay (obtuse) cells are corrected or not
    logical                                    :: coriolis    ! indicates whether Coriolis force should be included or not
    logical                                    :: corrdep     ! water depth in velocity point is corrected to 2nd order accuracy or not
    logical                                    :: depcds      ! water depth in velocity point approximated by central differences only
    logical                                    :: horwinc     ! horizontal terms of w-momentum equation included or not
    logical                                    :: linswe      ! indicates whether the shallow equations are linearized or not
    logical                                    :: lpproj      ! indicates whether non-hydrostatic pressure part in pressure projection method is included in continuity equation or not
    logical                                    :: lprecon     ! indicates whether preconditioner should be applied or not
    logical                                    :: lreptx      ! indicates whether the domain repeats itself in x-direction or not
    logical                                    :: lrepty      ! indicates whether the domain repeats itself in y-direction or not
    logical                                    :: mimetic     ! indicates whether advection term is dicretized mimetically (in both time and space) or not
    logical                                    :: momskip     ! indicates whether the momentum equations should be solved or not
    logical                                    :: numdisp     ! indicates whether exact dispersion relation should be applied or an approximate one based on box scheme
    logical                                    :: oned        ! indicates whether the calculation should be performed in 1D mode
    logical                                    :: stricthead  ! indicates whether advection term should be strictly energy head conservative
    logical                                    :: strictmom   ! indicates whether advection term should be strictly momentum conservative
    logical                                    :: svwp        ! indicates whether space varying wind and atmospheric pressure should be included or not
    logical                                    :: verwinc     ! vertical terms of w-momentum equation included or not
!
!   Source text
!
end module SwashCommdata3
!
module SwashCommdata4
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
!    1.00, February 2010: New Module based on module SWCOMM4
!
!   Purpose
!
!   Module containing data of test output and spherical coordinates
!
!   Method
!
!   Data with respect to test output and spherical coordinates
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
!   ---
!
!   Module variables
!
    ! information for test output
    !
    integer                                    :: intes   ! testing parameter,
                                                          ! =30; for command INTE, not documented in the manual
                                                          ! =intes; for command INTE [intes], not documented
    integer                                    :: ioutes  ! minimum value for ITEST,
                                                          ! =30; for command OUTE, not documented in the manual
                                                          ! =itest; for command OUTE [itest], not documented
    integer                                    :: iptst   ! sequence number of a test point
    integer                                    :: ifpar   ! unit reference number for output of flow variables in test points
    integer                                    :: maxmes  ! set by command SET [maxmes]
    integer                                    :: nptst   ! number of test points; set by command TEST
    integer                                    :: nptsta  ! number of test points, equal to max(1,nptst)
    !
    logical                                    :: testfl  ! test output must or must not be made, mainly for test points
    !
    !$omp threadprivate(iptst,testfl)
    !
    ! information for spherical coordinates
    !
    integer                                    :: kspher  ! indicates whether spherical coordinates are used or not
                                                          ! =0; Cartesian coordinates
                                                          ! >0; spherical coordinates
    integer                                    :: mproj   ! projection method
                                                          ! =0; (quasi-)Cartesian
                                                          ! =1; uniform Mercator (only spherical coordinates)
    !
    real                                       :: lendeg  ! length of a degree of the sphere
    real                                       :: rearth  ! radius of the Earth
!
!   Source text
!
end module SwashCommdata4
!
module SwashTimecomm
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
!    1.00, February 2010: New Module based on module TIMECOMM
!
!   Purpose
!
!   Module containing data of time-related parameters and timings
!
!   Method
!
!   Data with respect to time-related parameters and timings
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
!TIMG    integer, parameter                          :: nsectm = 300 ! number of sections that may be defined for an application
!TIMG    integer, parameter                          :: mxtimr = 10  ! maximum number of simultaneous timings, dimension of work-arrays listtm, timers
!
!   Module variables
!
    ! information for time-related parameters
    !
    real*8                                      :: dt     ! time step of the computation (in seconds)
                                                          ! =deltc; set by command COMP ... [deltc] ...
                                                          ! =deltc*60; set by command COMP ... [deltc] MI ...
                                                          ! =deltc*60*60; set by command COMP ... [deltc] HR ...
                                                          ! =deltc*60*60*24; set by command COMP ... [deltc] DA ...
    real*8                                      :: rdtim  ! =0 when in stationary mode
                                                          ! =1/dt when in nonstationary mode
    real*8                                      :: timco  ! time and date of the computation during the simulation (in seconds since the reference day (REFDAY))
    real*8                                      :: tinic  ! start time and date of the computation (in seconds since the reference day (REFDAY))
                                                          ! =tbegc; set by command COMP [tbegc] ...
    real*8                                      :: tfinc  ! end time and date of the computation (in seconds since the reference day (REFDAY))
    !
    character(20)                               :: chtime ! character string representation of date-time of computation
!TIMG    !
!TIMG    ! information for cpu and wall-clock times
!TIMG    !
!TIMG    integer, dimension(nsectm), save            :: ncumtm ! for each section of the application the number of timings that contributed to time in dcputm
!TIMG    integer, dimension(mxtimr), save            :: listtm ! list of section numbers for all running/active timers in array timers.
!TIMG                                                          ! A value of -1 signals that the corresponding timer is not running
!TIMG    integer, save                               :: lasttm ! last occupied position in listtm, 0 if all positions in listtm are free
!TIMG    double precision, dimension(nsectm,2), save :: dcumtm ! cumulative time; columns 1,2: cpu-time, wall-clock time
!TIMG    double precision, dimension(mxtimr,2), save :: timers ! start-time of active timers; columns 1,2: cpu-time, wall-clock time
!TIMG    !
!TIMG    !$omp threadprivate(dcumtm,timers,ncumtm,listtm,lasttm)
!
!   Source text
!
end module SwashTimecomm
