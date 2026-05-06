module SwashFlowdata
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
!    1.00, January 2010: New Module
!
!   Purpose
!
!   Module containing data for flow computation
!
!   Method
!
!   Data with respect to flow properties
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
    integer                                       :: mf           ! =1;     index of first point in x-direction of computational grid
                                                                  !         first wl-point considering as a virtual one
                                                                  !         first  u-point in computational grid
    integer                                       :: mfu          ! =2;     index of point mf+1 in x-direction of computational grid
                                                                  !         first wl-point considering as an internal one
    integer                                       :: ml           ! =mxc-1; index of last point in x-direction of computational grid
                                                                  !         last wl-point considering as an internal one
                                                                  !         last  u-point in computational grid
    integer                                       :: mlu          ! =mxc;   index of point ml+1 in x-direction of computational grid
                                                                  !         last wl-point considering as a virtual one
    integer                                       :: nf           ! =1;     index of first point in y-direction of computational grid
                                                                  !         first wl-point considering as a virtual one
                                                                  !         first  v-point in computational grid
    integer                                       :: nfu          ! =2;     index of point nf+1 in y-direction of computational grid
                                                                  !         first wl-point considering as an internal one
    integer                                       :: nl           ! =myc-1; index of last point in y-direction of computational grid
                                                                  !         last wl-point considering as an internal one
                                                                  !         last  v-point in computational grid
    integer                                       :: nlu          ! =myc;   index of point nl+1 in y-direction of computational grid
                                                                  !         last wl-point considering as a virtual one
    !
    integer, dimension(:)  , save, allocatable    :: brks         ! mask value for wave breaking at wl-point
                                                                  ! =0; no breaking
                                                                  ! =1; breaking
    integer, dimension(:)  , save, allocatable    :: ibb          ! boundary condition type at lower boundary
                                                                  ! = 1; closed boundary
                                                                  ! = 2; water level opening
                                                                  ! = 3; velocity opening
                                                                  ! = 4; not used
                                                                  ! = 5; discharge opening
                                                                  ! = 6; Riemann invariant opening
                                                                  ! = 7; weakly reflective opening
                                                                  ! = 8; Sommerfeld radiation condition
                                                                  ! =10; outflow condition
                                                                  ! < 0; =-ibb, where both normal and tangential components are described
    integer, dimension(:)  , save, allocatable    :: ibl          ! boundary condition type at left boundary
                                                                  ! see variable ibb for meaning
    integer, dimension(:)  , save, allocatable    :: ibr          ! boundary condition type at right boundary
                                                                  ! see variable ibb for meaning
    integer, dimension(:)  , save, allocatable    :: ibt          ! boundary condition type at upper boundary
                                                                  ! see variable ibb for meaning
    integer, dimension(:)  , save, allocatable    :: lfbs         ! label of floating body in wl-point
    integer, dimension(:)  , save, allocatable    :: mcell        ! cell-based index table for each time-step level (unstructured-mesh transport module)
    integer, dimension(:)  , save, allocatable    :: mlev         ! time-step level for local time stepping (unstructured-mesh transport module)
    integer, dimension(:)  , save, allocatable    :: presp        ! mask value for free surface/pressurized flow at pressure/wl-point
                                                                  ! =0; free surface flow
                                                                  ! =1; pressurized flow
    integer, dimension(:)  , save, allocatable    :: presu        ! mask value for free surface/pressurized flow at u-point
                                                                  ! =0; free surface flow
                                                                  ! =1; pressurized flow
    integer, dimension(:)  , save, allocatable    :: presv        ! mask value for free surface/pressurized flow at v-point
                                                                  ! =0; free surface flow
                                                                  ! =1; pressurized flow
    integer, dimension(:)  , save, allocatable    :: uwetp        ! indicate a wet and valid u-point for 3D computation
                                                                  ! =0; dry
                                                                  ! =1; wet
    integer, dimension(:)  , save, allocatable    :: vwetp        ! indicate a wet and valid v-point for 3D computation
                                                                  ! =0; dry
                                                                  ! =1; wet
    integer, dimension(:)  , save, allocatable    :: wets         ! mask value for wetting and drying at wl-point
                                                                  ! =0; dry
                                                                  ! =1; wet
    integer, dimension(:)  , save, allocatable    :: wetu         ! mask value for wetting and drying at u-point
                                                                  ! =0; dry
                                                                  ! =1; wet
    integer, dimension(:)  , save, allocatable    :: wetv         ! mask value for wetting and drying at v-point
                                                                  ! =0; dry
                                                                  ! =1; wet
    !
    integer                                       :: nconct       ! number of connections in one row of the matrix
    integer, dimension(:)  , save, allocatable    :: ishif        ! the shifts of the layer index for a given point in the Keller-box scheme
    !
    real  , dimension(:)      , save, allocatable :: a            ! lower diagonal of tri-diagonal matrix
    real  , dimension(:)      , save, allocatable :: advec        ! contribution to advective term (1DH)
    real  , dimension(:)      , save, allocatable :: advecx       ! contribution to advective term in x-direction
    real  , dimension(:)      , save, allocatable :: advecy       ! contribution to advective term in y-direction
    real  , dimension(:,:)    , save, allocatable :: amat         ! coefficient matrix of the system of equations
    real  , dimension(:,:,:)  , save, allocatable :: amatc        ! tri-diagonal matrix containing vertical terms in transport equation
    real*8, dimension(:,:)    , save, allocatable :: amatn        ! coefficient matrix of the mildly nonlinear system of equations
    real  , dimension(:,:,:)  , save, allocatable :: amatp        ! coefficient matrix of Poisson equation for non-hydrostatic pressure correction
    real  , dimension(:,:,:)  , save, allocatable :: amatt        ! coefficient matrix of turbulence model equations
    real  , dimension(:,:,:)  , save, allocatable :: amatu        ! tri-diagonal matrix containing vertical terms in u-momentum equation
    real  , dimension(:,:,:)  , save, allocatable :: amatw        ! tri-diagonal matrix containing vertical terms in w-momentum equation
    real*8, dimension(:)      , save, allocatable :: an           ! lower diagonal of tri-diagonal matrix of mildly nonlinear system
    real  , dimension(:,:)    , save, allocatable :: apoks        ! laminar friction factor at center of layer interfaces
    real  , dimension(:,:)    , save, allocatable :: apomu        ! laminar friction factor in u-point based on mean porosity
    real  , dimension(:,:)    , save, allocatable :: apomv        ! laminar friction factor in v-point based on mean porosity
    real  , dimension(:)      , save, allocatable :: apors        ! laminar friction factor in wl-point
    real  , dimension(:)      , save, allocatable :: aporu        ! laminar friction factor in u-point
    real  , dimension(:)      , save, allocatable :: aporv        ! laminar friction factor in v-point
    real  , dimension(:)      , save, allocatable :: b            ! main diagonal of tri-diagonal matrix
    real  , dimension(:,:,:)  , save, allocatable :: bcrp         ! prescribed constituent at boundary face (current time level)
    real  , dimension(:,:,:)  , save, allocatable :: bcrpo        ! prescribed constituent at boundary face (previous time level)
    real  , dimension(:)      , save, allocatable :: bcs          ! prescribed water level at boundary face (current time level)
    real  , dimension(:)      , save, allocatable :: bcso         ! prescribed water level at boundary face (previous time level)
    real*8, dimension(:)      , save, allocatable :: bn           ! main diagonal of tri-diagonal matrix of mildly nonlinear system
    real  , dimension(:,:)    , save, allocatable :: bdx          ! coefficient matrix of horizontal u-momentum equation corresponding to point (m-1,n,k)
    real  , dimension(:,:)    , save, allocatable :: bdy          ! coefficient matrix of horizontal u-momentum equation corresponding to point (m,n-1,k)
    real  , dimension(:,:)    , save, allocatable :: bpoks        ! turbulent friction factor at center of layer interfaces
    real  , dimension(:,:)    , save, allocatable :: bpomu        ! turbulent friction factor in u-point based on mean porosity
    real  , dimension(:,:)    , save, allocatable :: bpomv        ! turbulent friction factor in v-point based on mean porosity
    real  , dimension(:)      , save, allocatable :: bpors        ! turbulent friction factor in wl-point
    real  , dimension(:)      , save, allocatable :: bporu        ! turbulent friction factor in u-point
    real  , dimension(:)      , save, allocatable :: bporv        ! turbulent friction factor in v-point
    real  , dimension(:,:)    , save, allocatable :: bux          ! coefficient matrix of horizontal u-momentum equation corresponding to point (m+1,n,k)
    real  , dimension(:,:)    , save, allocatable :: buy          ! coefficient matrix of horizontal u-momentum equation corresponding to point (m,n+1,k)
    real  , dimension(:)      , save, allocatable :: c            ! upper diagonal of tri-diagonal matrix
    real  , dimension(:)      , save, allocatable :: cbot         ! bottom friction term
    real  , dimension(:,:,:)  , save, allocatable :: cbndb        ! boundary condition of Dirichlet type for transport constituent at lower boundary
    real  , dimension(:,:,:)  , save, allocatable :: cbndl        ! boundary condition of Dirichlet type for transport constituent at left boundary
    real  , dimension(:,:,:)  , save, allocatable :: cbndr        ! boundary condition of Dirichlet type for transport constituent at right boundary
    real  , dimension(:,:,:)  , save, allocatable :: cbndt        ! boundary condition of Dirichlet type for transport constituent at upper boundary
    real  , dimension(:,:,:)  , save, allocatable :: cbndu        ! boundary condition of Dirichlet type for transport constituent at boundary of unstructured mesh
    real  , dimension(:)      , save, allocatable :: cfricu       ! friction coefficient in u-point
    real  , dimension(:)      , save, allocatable :: cfricv       ! friction coefficient in v-point
    real  , dimension(:,:)    , save, allocatable :: cfu          ! stored Coriolis forcing in u-point for use at two time levels with AB2 scheme
    real  , dimension(:,:)    , save, allocatable :: cfv          ! stored Coriolis forcing in v-point for use at two time levels with AB2 scheme
    real*8, dimension(:)      , save, allocatable :: cn           ! upper diagonal of tri-diagonal matrix of mildly nonlinear system
    real  , dimension(:)      , save, allocatable :: corf         ! Coriolis force term
    real  , dimension(:,:)    , save, allocatable :: coutb        ! computed salinity at lower open boundary when tidal flow returns inward
    real  , dimension(:,:)    , save, allocatable :: coutl        ! computed salinity at left open boundary when tidal flow returns inward
    real  , dimension(:,:)    , save, allocatable :: coutr        ! computed salinity at right open boundary when tidal flow returns inward
    real  , dimension(:,:)    , save, allocatable :: coutt        ! computed salinity at upper open boundary when tidal flow returns inward
    real  , dimension(:,:)    , save, allocatable :: coutu        ! computed salinity at open boundary when tidal flow returns inward (unstructured mesh)
    real  , dimension(:,:)    , save, allocatable :: cpoks        ! added mass coefficient of porous medium flow at center of layer interfaces
    real  , dimension(:,:)    , save, allocatable :: cpomu        ! added mass coefficient of porous medium flow in u-point based on mean porosity
    real  , dimension(:,:)    , save, allocatable :: cpomv        ! added mass coefficient of porous medium flow in v-point based on mean porosity
    real  , dimension(:)      , save, allocatable :: cpors        ! added mass coefficient of porous medium flow in wl-point
    real  , dimension(:)      , save, allocatable :: cporu        ! added mass coefficient of porous medium flow in u-point
    real  , dimension(:)      , save, allocatable :: cporv        ! added mass coefficient of porous medium flow in v-point
    real  , dimension(:,:,:)  , save, allocatable :: cvegu        ! vegetation coefficients in u-point
    real  , dimension(:,:,:)  , save, allocatable :: cvegv        ! vegetation coefficients in v-point
    real  , dimension(:)      , save, allocatable :: cwndu        ! wind stress coefficient in u-point
    real  , dimension(:)      , save, allocatable :: cwndv        ! wind stress coefficient in v-point
    real  , dimension(:)      , save, allocatable :: d            ! right-hand side of tri-diagonal system of equations
    real  , dimension(:,:)    , save, allocatable :: divu         ! divergence of flow velocity per cell/layer
    real  , dimension(:,:)    , save, allocatable :: divq         ! divergence of mass flux per cell/layer
    real  , dimension(:,:,:)  , save, allocatable :: dmat         ! divergence matrix (structured grids)
    real  , dimension(:,:,:,:), save, allocatable :: dmatu        ! divergence matrix (unstructured grids)
    real*8, dimension(:)      , save, allocatable :: dn           ! right-hand side of tri-diagonal matrix of mildly nonlinear system
    real  , dimension(:)      , save, allocatable :: dps          ! bottom depth in wl-point
    real  , dimension(:)      , save, allocatable :: dpu          ! bottom depth in u-point (tiled)
    real  , dimension(:)      , save, allocatable :: dpv          ! bottom depth in v-point (tiled)
    real  , dimension(:,:)    , save, allocatable :: dq           ! non-hydrostatic pressure correction
    real  , dimension(:,:)    , save, allocatable :: dqgrd        ! gradient of pressure correction to be employed in pressure projection method (1DH)
    real  , dimension(:,:)    , save, allocatable :: dqgrdu       ! gradient of pressure correction in x-direction to be employed in pressure projection method (2DH)
    real  , dimension(:,:)    , save, allocatable :: dqgrdv       ! gradient of pressure correction in y-direction to be employed in pressure projection method (2DH)
    real  , dimension(:,:)    , save, allocatable :: dqv          ! non-hydrostatic pressure correction on fine (velocity) grid
    real  , dimension(:)      , save, allocatable :: dq0          ! non-hydrostatic pressure correction for a pressure layer in layer-averaged case
    real  , dimension(:)      , save, allocatable :: ds           ! water level correction
    real  , dimension(:)      , save, allocatable :: dsrc         ! distance source to boundary for internal wave generation (unstructured grids)
    real  , dimension(:)      , save, allocatable :: dtc          ! cell-based time step for local time stepping (unstructured-mesh transport module)
    real  , dimension(:,:)    , save, allocatable :: fcor         ! Coriolis parameter
    real  , dimension(:)      , save, allocatable :: flos         ! draft of floating object in wl-point
    real  , dimension(:)      , save, allocatable :: flou         ! draft of floating object in u-point
    real  , dimension(:)      , save, allocatable :: flov         ! draft of floating object in v-point
    real  , dimension(:,:,:)  , save, allocatable :: flux         ! total flux at cell-faces for transport constituent
    real  , dimension(:,:)    , save, allocatable :: fluxt        ! total flux at unstructured mesh faces for transport constituent
    real  , dimension(:,:,:)  , save, allocatable :: gmatu        ! gradient matrix for pressure in u-point
    real  , dimension(:,:,:)  , save, allocatable :: gmatuv       ! gradient matrix for pressure in u-point in a velocity layer
    real  , dimension(:,:)    , save, allocatable :: gmatu0       ! gradient matrix for pressure in u-point in one pressure layer
    real  , dimension(:,:,:)  , save, allocatable :: gmatv        ! gradient matrix for pressure in v-point
    real  , dimension(:,:,:)  , save, allocatable :: gmatvv       ! gradient matrix for pressure in v-point in a velocity layer
    real  , dimension(:,:)    , save, allocatable :: gmatv0       ! gradient matrix for pressure in v-point in one pressure layer
    real  , dimension(:,:,:)  , save, allocatable :: gmatw        ! gradient matrix for pressure in w-point
    real  , dimension(:)      , save, allocatable :: hindun       ! represents inundation depth in wl-point
                                                                  ! =0; not inundated
                                                                  ! =1; inundated
    real  , dimension(:,:)    , save, allocatable :: hks          ! layer thicknesses in wl-point
    real  , dimension(:,:)    , save, allocatable :: hksc         ! layer thicknesses in wl-point on coarse (pressure) grid
    real  , dimension(:,:)    , save, allocatable :: hksnew       ! layer thicknesses in wl-point at new time level
    real  , dimension(:,:)    , save, allocatable :: hkso         ! layer thicknesses in wl-point at old time level
    real  , dimension(:,:)    , save, allocatable :: hku          ! layer thicknesses in u-point based on upwinding
    real  , dimension(:,:)    , save, allocatable :: hkuc         ! layer thicknesses in u-point based on upwinding on coarse (pressure) grid
    real  , dimension(:,:)    , save, allocatable :: hkum         ! layer thicknesses in u-point based on averaging
    real  , dimension(:,:)    , save, allocatable :: hkumc        ! layer thicknesses in u-point based on averaging on coarse (pressure) grid
    real  , dimension(:,:)    , save, allocatable :: hkumn        ! layer thicknesses in u-point based on averaging extrapolated in time
    real  , dimension(:,:)    , save, allocatable :: hkv          ! layer thicknesses in v-point based on upwinding
    real  , dimension(:,:)    , save, allocatable :: hkvc         ! layer thicknesses in v-point based on upwinding on coarse (pressure) grid
    real  , dimension(:,:)    , save, allocatable :: hkvm         ! layer thicknesses in v-point based on averaging
    real  , dimension(:,:)    , save, allocatable :: hkvmc        ! layer thicknesses in v-point based on averaging on coarse (pressure) grid
    real  , dimension(:,:)    , save, allocatable :: hkvmn        ! layer thicknesses in v-point based on averaging extrapolated in time
    real  , dimension(:)      , save, allocatable :: hpors        ! porous structure height in wl-point
    real  , dimension(:)      , save, allocatable :: hporu        ! porous structure height in u-point
    real  , dimension(:)      , save, allocatable :: hporv        ! porous structure height in v-point
    real  , dimension(:)      , save, allocatable :: hs           ! water depth in wl-point
    real  , dimension(:)      , save, allocatable :: hso          ! water depth in wl-point at old time level
    real  , dimension(:)      , save, allocatable :: hu           ! water depth in u-point based on upwinding
    real  , dimension(:)      , save, allocatable :: hum          ! water depth in u-point based on averaging
    real  , dimension(:)      , save, allocatable :: humn         ! water depth in u-point based on averaging extrapolated in time
    real  , dimension(:)      , save, allocatable :: humo         ! water depth in u-point based on averaging at old time level
    real  , dimension(:)      , save, allocatable :: hv           ! water depth in v-point based on upwinding
    real  , dimension(:)      , save, allocatable :: hvm          ! water depth in v-point based on averaging
    real  , dimension(:)      , save, allocatable :: hvmn         ! water depth in v-point based on averaging extrapolated in time
    real  , dimension(:)      , save, allocatable :: hvmo         ! water depth in v-point based on averaging at old time level
    real  , dimension(:,:)    , save, allocatable :: icretb       ! step counter for constituent return time at lower open boundary (countdown)
    real  , dimension(:,:)    , save, allocatable :: icretl       ! step counter for constituent return time at left open boundary (countdown)
    real  , dimension(:,:)    , save, allocatable :: icretr       ! step counter for constituent return time at right open boundary (countdown)
    real  , dimension(:,:)    , save, allocatable :: icrett       ! step counter for constituent return time at upper open boundary (countdown)
    real  , dimension(:,:)    , save, allocatable :: icretu       ! step counter for constituent return time at open boundary of unstructured mesh (countdown)
    real  , dimension(:,:)    , save, allocatable :: logfrc       ! bottom friction coefficients based on logarithmic wall-law in wl-point
    real  , dimension(:,:)    , save, allocatable :: logwlb       ! wall shear stress at lower boundary based on logarithmic wall-law
    real  , dimension(:,:)    , save, allocatable :: logwll       ! wall shear stress at left boundary based on logarithmic wall-law
    real  , dimension(:,:)    , save, allocatable :: logwlr       ! wall shear stress at right boundary based on logarithmic wall-law
    real  , dimension(:,:)    , save, allocatable :: logwlt       ! wall shear stress at upper boundary based on logarithmic wall-law
    real*8, dimension(:)      , save, allocatable :: lon          ! lower bound of solution of mildly nonlinear system
    real  , dimension(:,:)    , save, allocatable :: npoks        ! volumetric porosity at center of layer interfaces
    real  , dimension(:)      , save, allocatable :: npors        ! volumetric porosity in wl-point
    real  , dimension(:)      , save, allocatable :: nporu        ! volumetric porosity in u-point
    real  , dimension(:)      , save, allocatable :: nporv        ! volumetric porosity in v-point
    real  , dimension(:)      , save, allocatable :: patm         ! atmospheric pressure in wl-point
    real  , dimension(:)      , save, allocatable :: pfac         ! factor for hydrostatic pressure gradient at water level boundary
                                                                  ! = 0.5; take distance between boundary face and cell center twice
                                                                  ! = 1; represents distance between boundary face and cell center
    real  , dimension(:)      , save, allocatable :: pgrad        ! pressure gradient term due to density and atmospheric pressure
    real  , dimension(:,:)    , save, allocatable :: q            ! non-hydrostatic pressure
    real  , dimension(:)      , save, allocatable :: qbot         ! non-hydrostatic pressure at bottom
    real  , dimension(:)      , save, allocatable :: qgrad        ! pressure gradient term
    real  , dimension(:,:)    , save, allocatable :: qm           ! mass flux in wl-point based on averaging (1DH)
    real  , dimension(:,:)    , save, allocatable :: qn           ! mass flux at triangular face (unstructured grids)
    real  , dimension(:,:,:)  , save, allocatable :: quf          ! mass flux at cell face in flow direction (unstructured grids)
    real  , dimension(:,:)    , save, allocatable :: qv           ! non-hydrostatic pressure on fine (velocity) grid
    real  , dimension(:,:)    , save, allocatable :: qx           ! mass flux per width in x-direction
    real  , dimension(:,:)    , save, allocatable :: qxm          ! mass flux per width in x-direction based on averaging
    real  , dimension(:,:)    , save, allocatable :: qy           ! mass flux per width in y-direction
    real  , dimension(:,:)    , save, allocatable :: qym          ! mass flux per width in y-direction based on averaging
    real  , dimension(:,:)    , save, allocatable :: rho          ! water density relative to reference density
    real  , dimension(:)      , save, allocatable :: rhs          ! right-hand side of the system of equations
    real  , dimension(:,:)    , save, allocatable :: rhsc         ! right-hand side of tri-diagonal system of transport equation
    real*8, dimension(:)      , save, allocatable :: rhsn         ! right-hand side of the mildly nonlinear system of equations
    real  , dimension(:,:)    , save, allocatable :: rhsp         ! right-hand side of Poisson equation for non-hydrostatic pressure correction
    real  , dimension(:,:)    , save, allocatable :: rhst         ! right-hand side of turbulence model equations
    real  , dimension(:,:)    , save, allocatable :: rhsu         ! right-hand side of tri-diagonal system of u-momentum equation
    real  , dimension(:,:)    , save, allocatable :: rhsw         ! right-hand side of tri-diagonal system of w-momentum equation
    real  , dimension(:,:,:)  , save, allocatable :: rp           ! transport constituents in wl-point
    real  , dimension(:,:)    , save, allocatable :: rp0          ! cell-based transport constituent at previous substep (for local time stepping)
    real  , dimension(:,:)    , save, allocatable :: rp1          ! cell-based transport constituent at current substep (for local time stepping)
    real  , dimension(:,:,:)  , save, allocatable :: rpi          ! intermediate transport constituent at various time-step levels (for local time stepping)
    real  , dimension(:,:)    , save, allocatable :: rpo          ! transport constituents in wl-point at previous time level
    real  , dimension(:,:)    , save, allocatable :: rsuu         ! Reynolds stress component -u'u'
    real  , dimension(:,:)    , save, allocatable :: rsuv         ! Reynolds stress component -u'v'
    real  , dimension(:,:)    , save, allocatable :: rsuw         ! Reynolds stress component -u'w'
    real  , dimension(:,:)    , save, allocatable :: rsvu         ! Reynolds stress component -v'u'
    real  , dimension(:,:)    , save, allocatable :: rsvv         ! Reynolds stress component -v'v'
    real  , dimension(:,:)    , save, allocatable :: rsvw         ! Reynolds stress component -v'w'
    real  , dimension(:,:)    , save, allocatable :: rswu         ! Reynolds stress component -w'u'
    real  , dimension(:,:)    , save, allocatable :: rswv         ! Reynolds stress component -w'v'
    real  , dimension(:,:)    , save, allocatable :: rsww         ! Reynolds stress component -w'w'
    real  , dimension(:,:,:)  , save, allocatable :: rtur         ! turbulence quantities at center of layer interfaces
    real  , dimension(:)      , save, allocatable :: s0           ! water level at previous time level
    real  , dimension(:)      , save, allocatable :: s1           ! water level at current time level
    real  , dimension(:)      , save, allocatable :: skc          ! time rate of change of displacement of rigid body (derived from kinematic boundary conditions)
    real  , dimension(:)      , save, allocatable :: smax         ! maximum water level stored
    real  , dimension(:)      , save, allocatable :: so           ! water level at old time level
    real  , dimension(:)      , save, allocatable :: sponxl       ! left sponge layer in x-direction
    real  , dimension(:)      , save, allocatable :: sponxr       ! right sponge layer in x-direction
    real  , dimension(:)      , save, allocatable :: sponyb       ! lower sponge layer in y-direction
    real  , dimension(:)      , save, allocatable :: sponyt       ! upper sponge layer in y-direction
    real  , dimension(:)      , save, allocatable :: srcm         ! mass source for internal wave generation
    real  , dimension(:,:)    , save, allocatable :: tbndx        ! target boundary values for left/right sponge layers to be made available to other nodes
    real  , dimension(:,:)    , save, allocatable :: tbndy        ! target boundary values for lower/upper sponge layers to be made available to other nodes
    real  , dimension(:)      , save, allocatable :: teta         ! space varying implicitness factor for time integration of continuity equation (1DH)
    real  , dimension(:)      , save, allocatable :: teta2        ! space varying implicitness factor for water level gradient
    real  , dimension(:)      , save, allocatable :: tetau        ! space varying implicitness factor for time integration of u-contribution of continuity equation (2DH)
    real  , dimension(:)      , save, allocatable :: tetav        ! space varying implicitness factor for time integration of v-contribution of continuity equation (2DH)
    real  , dimension(:,:)    , save, allocatable :: u0           ! u-velocity at previous time level
    real  , dimension(:,:)    , save, allocatable :: u0p          ! u-velocity on coarse (pressure) grid at previous time level
    real  , dimension(:,:)    , save, allocatable :: u1           ! u-velocity at current time level
    real  , dimension(:,:)    , save, allocatable :: u1p          ! u-velocity on coarse (pressure) grid at current time level
    real  , dimension(:,:)    , save, allocatable :: ua           ! advective velocity based on finite differencing
    real  , dimension(:,:)    , save, allocatable :: ua2          ! another advective velocity based on finite differencing
    real  , dimension(:,:)    , save, allocatable :: Ubck         ! ambient U-velocity in wl-point or u-point
    real  , dimension(:)      , save, allocatable :: udep         ! depth-averaged u-velocity
    real  , dimension(:,:)    , save, allocatable :: ui           ! intermediate u-velocity in iterative process of pressure projection method
    real  , dimension(:,:)    , save, allocatable :: up           ! upwind velocity (for energy head computation)
    real*8, dimension(:)      , save, allocatable :: upn          ! upper bound of solution of mildly nonlinear system
    real  , dimension(:,:,:)  , save, allocatable :: uvc          ! velocity components in cell circumcenter
    real  , dimension(:,:)    , save, allocatable :: v0           ! v-velocity at previous time level
    real  , dimension(:,:)    , save, allocatable :: v0p          ! v-velocity on coarse (pressure) grid at previous time level
    real  , dimension(:,:)    , save, allocatable :: v1           ! v-velocity at current time level
    real  , dimension(:,:)    , save, allocatable :: v1p          ! v-velocity on coarse (pressure) grid at current time level
    real  , dimension(:,:)    , save, allocatable :: Vbck         ! ambient V-velocity in wl-point or v-point
    real  , dimension(:)      , save, allocatable :: vdep         ! depth-averaged v-velocity
    real  , dimension(:,:)    , save, allocatable :: vi           ! intermediate v-velocity in iterative process of pressure projection method
    real  , dimension(:)      , save, allocatable :: visc         ! horizontal viscosity terms for momentum equations
    real  , dimension(:)      , save, allocatable :: vnu2d        ! horizontal eddy viscosity coefficient in depth point (2DH), wl-point (1DH) or cell face (unstructured)
    real  , dimension(:,:)    , save, allocatable :: vnu3d        ! vertical eddy viscosity coefficient at center of layer interfaces
    real  , dimension(:,:)    , save, allocatable :: w0           ! w-velocity at previous time level
    real  , dimension(:)      , save, allocatable :: w0bot        ! w-velocity at bottom at previous time level (1DH/2DH)
    real  , dimension(:,:)    , save, allocatable :: w0p          ! w-velocity on coarse (pressure) grid at previous time level
    real  , dimension(:)      , save, allocatable :: w0top        ! w-velocity at free surface at previous time level (1DH/2DH)
    real  , dimension(:,:)    , save, allocatable :: w1           ! w-velocity at current time level
    real  , dimension(:)      , save, allocatable :: w1bot        ! w-velocity at bottom at current time level (1DH/2DH)
    real  , dimension(:,:)    , save, allocatable :: w1p          ! w-velocity on coarse (pressure) grid at current time level
    real  , dimension(:)      , save, allocatable :: w1top        ! w-velocity at free surface at current time level (1DH/2DH)
    real  , dimension(:)      , save, allocatable :: windu        ! wind stress term in u-point
    real  , dimension(:)      , save, allocatable :: windv        ! wind stress term in v-point
    real  , dimension(:)      , save, allocatable :: wndimp       ! implicit part of wind stress term
    real  , dimension(:,:)    , save, allocatable :: wom          ! relative vertical velocity
    real  , dimension(:,:)    , save, allocatable :: womp         ! relative vertical velocity on coarse (pressure) grid
    real  , dimension(:,:)    , save, allocatable :: wphy         ! physical vertical velocity (only used for output purposes)
    real  , dimension(:)      , save, allocatable :: xshear       ! magnitude of shear squared in case of horizontal mixing (only used for flexible mesh)
    real  , dimension(:,:)    , save, allocatable :: zks          ! layer interfaces in wl-point
    real  , dimension(:,:)    , save, allocatable :: zksnew       ! layer interfaces in wl-point at new time level
    real  , dimension(:,:)    , save, allocatable :: zkso         ! layer interfaces in wl-point at old time level
    real  , dimension(:,:)    , save, allocatable :: zku          ! layer interfaces in u-point based on upwinding
    real  , dimension(:,:)    , save, allocatable :: zkum         ! layer interfaces in u-point based on averaging
    real  , dimension(:,:)    , save, allocatable :: zkuo         ! layer interfaces in u-point at old time level
    real  , dimension(:,:)    , save, allocatable :: zkv          ! layer interfaces in v-point based on upwinding
    real  , dimension(:,:)    , save, allocatable :: zkvm         ! layer interfaces in v-point based on averaging
    real  , dimension(:,:)    , save, allocatable :: zkvo         ! layer interfaces in v-point at old time level
    real  , dimension(:,:)    , save, allocatable :: zshear       ! magnitude of shear squared in case of vertical mixing
    !
    type sputype                                                  ! parameters sponge layer for flow velocity
      integer :: bface                                            ! closest boundary face
      real    :: gamma                                            ! relaxation parameter
    end type sputype

    type spstype                                                  ! parameters sponge layer for water level
      integer :: bcell                                            ! closest boundary cell
      real    :: gamma                                            ! relaxation parameter
    end type spstype
    !
    type(sputype), dimension(:,:), save, allocatable :: sponu     ! sponge layer and its attributes for flow velocity
    type(spstype), dimension(:,:), save, allocatable :: spons     ! sponge layer and its attributes for water level
!
!   Source text
!
end module SwashFlowdata
