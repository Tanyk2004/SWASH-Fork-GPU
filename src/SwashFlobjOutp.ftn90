subroutine SwashFlobjOutp ( oqproc, nvoqp, voqr, voq )
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
!    1.00: Dirk Rijnsdorp
!
!   Updates
!
!    1.00,   February 2018: New subroutine
!
!   Purpose
!
!   Requests hydrodynamic loads and/or body motions for the purpose of output
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata1, only: nmovar, outpar
    use SwashCommdata3, only: ifloat, pi
    use SwashTimeComm, only: timco
    use SwashRigBoddata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                     :: nvoqp  ! number of quantities
    integer, dimension(nmovar), intent(in)  :: voqr   ! place of each output quantity
    real   , dimension(nvoqp) , intent(out) :: voq    ! output quantity at request
    logical, dimension(nmovar), intent(in)  :: oqproc ! indicates whether or not an output quantity must be processed
!
!   Local variables
!
    integer, save                   :: ient = 0 ! number of entries in this subroutine
    integer                         :: j        ! loop counter
    integer                         :: k        ! loop counter
    integer                         :: l        ! counter
    integer                         :: nflob    ! number of floating bodies
    !
    real, dimension(:), allocatable :: fx       ! hydrodynamic force in x-direction (surge)
    real, dimension(:), allocatable :: fy       ! hydrodynamic force in y-direction (sway)
    real, dimension(:), allocatable :: fz       ! hydrodynamic force in z-direction (heave)
    !
    real, dimension(:), allocatable :: mx       ! hydrodynamic moment in x-direction (roll)
    real, dimension(:), allocatable :: my       ! hydrodynamic moment in y-direction (pitch)
    real, dimension(:), allocatable :: mz       ! hydrodynamic moment in z-direction (yaw)
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashFlobjOutp')
    !
    ! determine number of floating bodies
    !
    if ( ifloat == 1 ) then
       nflob = 1
    else if ( ifloat == 2 ) then
       nflob = mbod
    endif
    !
    ! allocate arrays for storing forces and moments
    !
    allocate (fx(nflob))
    allocate (fy(nflob))
    allocate (fz(nflob))
    !
    allocate (mx(nflob))
    allocate (my(nflob))
    allocate (mz(nflob))
    !
    ! determine hydrodynamic loads acting on ...
    !
    if ( ifloat == 1 ) then
       !
       ! a fixed body
       !
       do j = 1, ndim
          do k = 1, 2
             l = j + ndim*(k-1)
             bdof(1,j,k) = oqproc(100+l)
          enddo
       enddo
       !
       call SwashHydroLoads ( nflob, fx, fy, fz, mx, my, mz )
       !
    else if ( ifloat == 2 ) then
       !
       ! moving bodies
       !
       fx(:) = fhyd1(:,1)
       fy(:) = fhyd1(:,2)
       fz(:) = fhyd1(:,3)
       !
       mx(:) = thyd1(:,1)
       my(:) = thyd1(:,2)
       mz(:) = thyd1(:,3)
       !
    endif
    !
    ! store the requested output quantities
    !
    ! Tsec
    !
    if ( oqproc(41) ) voq(voqr(41)) = real(timco) - outpar(1)
    !
    do j = 1, nflob
       !
       ! force in x-direction (surge)
       !
       if ( oqproc(101) ) voq(voqr(101)+j-1) = fx(j) / 1000.
       !
       ! force in y-direction (sway)
       !
       if ( oqproc(102) ) voq(voqr(102)+j-1) = fy(j) / 1000.
       !
       ! force in z-direction (heave)
       !
       if ( oqproc(103) ) voq(voqr(103)+j-1) = fz(j) / 1000.
       !
       ! rotational force in x-direction (roll)
       !
       if ( oqproc(104) ) voq(voqr(104)+j-1) = mx(j) / 1000.
       !
       ! rotational force in y-direction (pitch)
       !
       if ( oqproc(105) ) voq(voqr(105)+j-1) = my(j) / 1000.
       !
       ! rotational force in z-direction (yaw)
       !
       if ( oqproc(106) ) voq(voqr(106)+j-1) = mz(j) / 1000.
       !
    enddo
    !
    ! include body motion, if appropriate
    !
    if ( ifloat == 2 ) then
       !
       do j = 1, nflob
          !
          ! translation in x-direction
          !
          if ( oqproc(107) ) voq(voqr(107)+j-1) = xfot1(j,1)
          !
          ! translation in y-direction
          !
          if ( oqproc(108) ) voq(voqr(108)+j-1) = xfot1(j,2)
          !
          ! translation in z-direction
          !
          if ( oqproc(109) ) voq(voqr(109)+j-1) = xfot1(j,3)
          !
          ! rotation around x-axis
          !
          if ( oqproc(110) ) voq(voqr(110)+j-1) = xfor1(j,1) * 180./pi
          !
          ! rotation around y-axis
          !
          if ( oqproc(111) ) voq(voqr(111)+j-1) = xfor1(j,2) * 180./pi
          !
          ! rotation around z-axis
          !
          if ( oqproc(112) ) voq(voqr(112)+j-1) = xfor1(j,3) * 180./pi
          !
          ! PTO power
          !
          if ( oqproc(113) ) voq(voqr(113)+j-1) = ptop(j) / 1000.
          !
       enddo
       !
    endif
    !
    deallocate ( fx, fy, fz, mx, my, mz )
    !
end subroutine SwashFlobjOutp
