subroutine SwashMotionRigidBod
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
!    9.01: Dirk Rijnsdorp, Marcel Zijlema
!
!   Updates
!
!    9.01, October 2022: New subroutine
!
!   Purpose
!
!   Solves the equations of motion for floating rigid bodies
!
!   Method
!
!   The motion of a rigid body is governed by the force and torque balances derived from
!   the conservation of linear and angular momentum, respectively
!
!   The ODEs are numerically solved using the generalized-alpha method (Chung and Hulbert, 1993)
!
!   Special cases are the Newmark (1959) scheme which is non-dissipative, and the
!   Chung-Hulbert (1993) scheme which is a second order dissipative scheme (default in SWASH)
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm, only: dt
    use SwashRigBoddata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    !
    real          :: accu     ! auxiliary real to store acceleration
    real          :: alfa     ! under-relaxation factor
    real          :: beta     ! implicitness factor for generalized-alpha method
    real          :: f0       ! total force at previous time level
    real          :: f1       ! total force at current time level
    real          :: gamma    ! another implicitness factor for generalized-alpha method
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashMotionRigidBod')
    !
    ! if not moving rigid bodies, return
    !
    if ( ifloat /= 2 ) return
    !
    alfa  = pship(10)
    beta  = pship( 3)
    gamma = pship( 4)
    !
    ! translational motion
    !
    do n = 1, ndim
       !
       do m = 1, mbod
          !
          if ( .not.bdof(m,n,1) ) cycle
          !
          ! get total force per unit mass for body m
          !
          f0 = ( fhyd0(m,n) + fbod0(m,n) ) / bmass(m)
          f1 = ( fhyd1(m,n) + fbod1(m,n) ) / bmass(m)
          !
          ! save last iteration for under-relaxation
          !
          accu = afot1(m,n)
          !
          ! update acceleration of body under motion
          !
          afot1(m,n) = ( -alfam * afot0(m,n) + alfaf * f0 + (1.-alfaf) * f1 ) / (1.-alfam)
          !
          ! apply under-relaxation
          !
          afot1(m,n) = alfa * afot1(m,n) + (1.-alfa) * accu
          !
          ! update body velocity and displacement based on under-relaxed acceleration
          !
          vfot1(m,n) = vfot0(m,n) + dt * ( gamma * afot1(m,n) + (1.-gamma) * afot0(m,n) )
          !
          xfot1(m,n) = xfot0(m,n) + dt * ( vfot0(m,n) + dt * ( beta * afot1(m,n) + (0.5-beta) * afot0(m,n) ) )
          !
       enddo
       !
    enddo
    !
    ! rotational motion
    !
    do n = 1, ndim
       !
       do m = 1, mbod
          !
          if ( .not.bdof(m,n,2) ) cycle
          !
          ! get total rotational force per unit angular mass for body m
          !
          f0 = ( thyd0(m,n) + tbod0(m,n) ) / bmoi(m,n)
          f1 = ( thyd1(m,n) + tbod1(m,n) ) / bmoi(m,n)
          !
          ! save last iteration for under-relaxation
          !
          accu = afor1(m,n)
          !
          ! update acceleration of body under motion
          !
          afor1(m,n) = ( -alfam * afor0(m,n) + alfaf * f0 + (1.-alfaf) * f1 ) / (1.-alfam)
          !
          ! apply under-relaxation
          !
          afor1(m,n) = alfa * afor1(m,n) + (1.-alfa) * accu
          !
          ! update body velocity and displacement based on under-relaxed acceleration
          !
          vfor1(m,n) = vfor0(m,n) + dt * ( gamma * afor1(m,n) + (1.-gamma) * afor0(m,n) )
          !
          xfor1(m,n) = xfor0(m,n) + dt * ( vfor0(m,n) + dt * ( beta * afor1(m,n) + (0.5-beta) * afor0(m,n) ) )
          !
       enddo
       !
    enddo
    !
end subroutine SwashMotionRigidBod
