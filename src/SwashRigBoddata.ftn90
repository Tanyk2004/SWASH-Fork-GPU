module SwashRigBoddata
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
!    1.00, October 2022: New Module
!
!   Purpose
!
!   Module containing data for rigid bodies computation
!
!   Method
!
!   Data with respect to rigid body properties
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
    integer, parameter                           :: ndim = 3 ! number of dimensions related to body motions
!
!   Module variables
!
    integer                                      :: mbod     ! number of rigid bodies
    integer                                      :: mfoti    ! method of time integration (generalized-alpha method)
                                                             ! =1; Newmark scheme
                                                             ! =2; Chung-Hulbert scheme
                                                             ! =3; Hilber-Hughes-Taylor scheme
                                                             ! =4; Wood-Bossak-Zienkiewicz scheme
    integer                                      :: mxfen    ! maximum number of fenders
    integer                                      :: mxmli    ! maximum number of mooring lines
    !
    real                                         :: alfaf    ! implicitness factor w.r.t. forces to control the amplification of high-frequency modes
                                                             ! note: this factor is internally determined by generalized-alpha scheme
    real                                         :: alfam    ! implicitness factor w.r.t. inertia to control the amplification of high-frequency modes
                                                             ! note: this factor is internally determined by generalized-alpha scheme
    !
    integer, dimension(:) , save, allocatable    :: nfeb     ! number of fenders per rigid body
    integer, dimension(:) , save, allocatable    :: nmlb     ! number of mooring lines per rigid body
    !
    real, dimension(:)    , save, allocatable    :: barea    ! the water plane area of the body
    real, dimension(:)    , save, allocatable    :: bcob     ! the center of buoyancy relative to the equilibrium position
    real, dimension(:,:)  , save, allocatable    :: bcog     ! the center of gravity with respect to problem coordinates
    real, dimension(:,:,:), save, allocatable    :: bfen     ! parameters for each fender connected to each rigid body
                                                             ! bfen(:,:,1): spring coefficient
                                                             ! bfen(:,:,2): x-component of attachment point at body relative to center of gravity
                                                             ! bfen(:,:,3): y-component of attachment point at body relative to center of gravity
                                                             ! bfen(:,:,4): z-component of attachment point at body relative to center of gravity
    real, dimension(:)    , save, allocatable    :: bmass    ! the body mass
    real, dimension(:,:,:), save, allocatable    :: bmli     ! mooring line parameters for each line connected to each rigid body
                                                             ! bmli(:,:, 1): spring coefficient
                                                             ! bmli(:,:, 2): damping coefficient
                                                             ! bmli(:,:, 3): x-component of attachment point at bottom relative to center of gravity
                                                             ! bmli(:,:, 4): y-component of attachment point at bottom relative to center of gravity
                                                             ! bmli(:,:, 5): z-component of attachment point at bottom relative to center of gravity
                                                             ! bmli(:,:, 6): x-component of attachment point at body relative to center of gravity
                                                             ! bmli(:,:, 7): y-component of attachment point at body relative to center of gravity
                                                             ! bmli(:,:, 8): z-component of attachment point at body relative to center of gravity
                                                             ! bmli(:,:, 9): equilibrium length
                                                             ! bmli(:,:,10): pretension coefficient
    real, dimension(:)    , save, allocatable    :: bmoax    ! the second moment of area for rigid body with respect to x-axis
    real, dimension(:)    , save, allocatable    :: bmoay    ! the second moment of area for rigid body with respect to y-axis
    real, dimension(:,:)  , save, allocatable    :: bmoi     ! the moment of inertia relative to the rotation axis
    real, dimension(:)    , save, allocatable    :: bvol     ! the volume of submerged part of the body
    !
    real, dimension(:,:)  , save, allocatable    :: extml0   ! extension of mooring line at previous time level
    real, dimension(:,:)  , save, allocatable    :: extml1   ! extension of mooring line at current time level
    !
    real, dimension(:)    , save, allocatable    :: ptop     ! power of PTO (power take-off) system (e.g. WEC)
    !
    real, dimension(:,:)  , save, allocatable    :: afor0    ! angular acceleration of rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: afor1    ! angular acceleration of rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: afot0    ! linear acceleration of rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: afot1    ! linear acceleration of rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: fbod0    ! spring-damping forces acting on rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: fbod1    ! spring-damping forces acting on rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: fhyd0    ! hydrodynamic forces acting on rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: fhyd1    ! hydrodynamic forces acting on rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: tbod0    ! spring-damping torques acting on rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: tbod1    ! spring-damping torques acting on rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: thyd0    ! hydrodynamic torques acting on rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: thyd1    ! hydrodynamic torques acting on rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: vfor0    ! angular velocity of rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: vfor1    ! angular velocity of rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: vfot0    ! linear velocity of rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: vfot1    ! linear velocity of rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: xfor0    ! angular displacement of rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: xfor1    ! angular displacement of rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: xfori    ! intermediate angular displacements in iterative process of fluid-structure interaction
    real, dimension(:,:)  , save, allocatable    :: xfot0    ! linear displacement of rigid body at previous time level
    real, dimension(:,:)  , save, allocatable    :: xfot1    ! linear displacement of rigid body at current time level
    real, dimension(:,:)  , save, allocatable    :: xfoti    ! intermediate linear displacements in iterative process of fluid-structure interaction
    !
    logical, dimension(:,:,:), save, allocatable :: bdof     ! indicates whether a degree of freedom of a rigid body
                                                             ! is included (true) or not (false). Meaning:
                                                             ! bdof(:,1,1): surge
                                                             ! bdof(:,2,1): sway
                                                             ! bdof(:,3,1): heave
                                                             ! bdof(:,1,2): roll
                                                             ! bdof(:,2,2): pitch
                                                             ! bdof(:,3,2): yaw
    logical, dimension(:)    , save, allocatable :: cpto     ! indicates whether pretension of mooring line should be included or not
    !
    type ribdat
       integer               :: label
       integer               :: nfen, nmli
       real                  :: parm(7)
       logical               :: pret
       logical               :: dof(6)
       real, pointer         :: fen(:,:)
       real, pointer         :: mli(:,:)
       type(ribdat), pointer :: nextrib
    end type ribdat
    !
    type(ribdat), save, target :: frigbod
!
!   Source text
!
end module SwashRigBoddata
