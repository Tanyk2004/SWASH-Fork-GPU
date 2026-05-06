subroutine SwashForcesRigidBod
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
!   Computes forces and torques on moving rigid bodies
!
!   Method
!
!   The following forces are added to the equations of motion:
!
!   1) Froude-Krylov forces induced by pressure in the fluid
!   2) forces due to mooring lines and fenders
!   3) restoring forces due to buoyancy
!
!   Note that the restoring forces have to be added explicitly since the
!   body motions are not included in the computational schematization
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
    integer, save               :: ient = 0 ! number of entries in this subroutine
    integer                     :: l        ! loop counter over mooring lines
    integer                     :: m        ! loop counter over bodies
    integer                     :: n        ! loop counter over body dimension
    !
    real                        :: fdmp     ! damping force due to fender
    real                        :: fspr     ! spring force due to fender
    real                        :: lml      ! length of mooring line
    !
    real, dimension(mxfen,3)    :: angfe    ! vector to store angle for each fender
    real, dimension(mxmli,3)    :: angml    ! vector to store angle for each mooring line
    real, dimension(mbod,mxmli) :: demldt   ! rate of change of mooring line length
    real, dimension(mxfen,3)    :: disfe    ! displacement vector for each fender
    real, dimension(mxmli,3)    :: disml    ! displacement vector for each mooring line
    real, dimension(mxfen)      :: ffe      ! force induced by each fender
    real, dimension(mxmli)      :: fml      ! force induced by each mooring line
    real, dimension(3)          :: vec      ! auxiliary vector to store position vector
    real, dimension(3)          :: vecr     ! auxiliary vector to store rotated vector
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashForcesRigidBod')
    !
    ! if not moving rigid bodies, return
    !
    if ( ifloat /= 2 ) return
    !
    ! compute hydrodynamic loads acting on rigid bodies
    !
    call SwashHydroLoads ( mbod, fhyd1(1,1), fhyd1(1,2), fhyd1(1,3), thyd1(1,1), thyd1(1,2), thyd1(1,3) )
    !
    ! compute spring and damping forces of mooring lines at attachment point to the rigid body
    !
    do m = 1, mbod
       !
       ! mooring line modelled as a linear spring-mass-damper system
       !
       do l = 1, nmlb(m)
          !
          ! determine initial attachment point to the body
          !
          vec(1) = bmli(m,l,6)
          vec(2) = bmli(m,l,7)
          vec(3) = bmli(m,l,8)
          !
          ! rotate attachment point with angular displacements ...
          !
          call rotatep ( vec, vecr, xfor1(m,1), xfor1(m,2), xfor1(m,3) )
          !
          disml(l,1) = vecr(1)
          disml(l,2) = vecr(2)
          disml(l,3) = vecr(3)
          !
          ! ... and translate it in each direction
          !
          vec(1) = disml(l,1) + xfot1(m,1)
          vec(2) = disml(l,2) + xfot1(m,2)
          vec(3) = disml(l,3) + xfot1(m,3)
          !
          ! compute new position relative to attachment point at the bottom
          !
          vec(1) = vec(1) - bmli(m,l,3)
          vec(2) = vec(2) - bmli(m,l,4)
          vec(3) = vec(3) - bmli(m,l,5)
          !
          ! compute length of mooring line after motion
          !
          lml = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )
          !
          ! compute extension of mooring line w.r.t. its equilibrium
          !
          extml1(m,l) = lml - bmli(m,l,9)
          !
          ! ... and its rate of change
          !
          demldt(m,l) = ( extml1(m,l) - extml0(m,l) ) / dt
          !
          ! compute total force due to each mooring line
          !
          fml(l) = bmli(m,l,10) + bmli(m,l,1) * extml1(m,l) + bmli(m,l,2) * demldt(m,l)
          !
          ! compute angle of mooring line with respect to the anchor at bottom
          !
          angml(l,1) = atan2( vec(1), sqrt( vec(2)**2 + vec(3)**2 ) )
          angml(l,2) = atan2( vec(2), sqrt( vec(3)**2 + vec(1)**2 ) )
          angml(l,3) = atan2( vec(3), sqrt( vec(1)**2 + vec(2)**2 ) )
          !
       enddo
       !
       ! fender modelled as a linear spring-mass system
       !
       do l = 1, nfeb(m)
          !
          ! determine initial attachment point to the body
          !
          vec(1) = bfen(m,l,2)
          vec(2) = bfen(m,l,3)
          vec(3) = bfen(m,l,4)
          !
          ! rotate attachment point with counter angular displacements ...
          !
          call rotatep ( vec, vecr, -xfor1(m,1), -xfor1(m,2), -xfor1(m,3) )
          !
          ! ... and move backward in each direction
          !
          vecr(1) = vecr(1) - xfot1(m,1)
          vecr(2) = vecr(2) - xfot1(m,2)
          vecr(3) = vecr(3) - xfot1(m,3)
          !
          disfe(l,1) = vecr(1)
          disfe(l,2) = vecr(2)
          disfe(l,3) = vecr(3)
          !
          ! if rotated fender position is 'inside' hull then compute force
          ! note: estimated as original position of fender and only considering the y-coordinate, which should be okay if motions are small compared to body size
          !
          if ( .not. vecr(2) > vec(2) ) then
             !
             fspr = bfen(m,l,1) * ( vec(2) - vecr(2) )
             fdmp = 1.e6 * vfot1(m,1)
             !
             ! compute total force
             !
             ffe(l) = sqrt( fspr**2 + fdmp**2 )
             !
             ! compute fender angle at connection point
             !
             angfe(l,1) = asin( fdmp /ffe(l) )
             angfe(l,2) = asin( fspr /ffe(l) )
             angfe(l,3) = 0.
             !
          else
             !
             ffe  (l  ) = 0.
             angfe(l,1) = 0.
             angfe(l,2) = 0.
             angfe(l,3) = 0.
             !
          endif
          !
       enddo
       !
       ! compute total spring-damping force acting on the body
       !
       do n = 1, ndim
          !
          fbod1(m,n) = 0.
          if ( .not.bdof(m,n,1) ) cycle
          do l = 1, nmlb(m)
             fbod1(m,n) = fbod1(m,n) - fml(l) * sin(angml(l,n))
          enddo
          do l = 1, nfeb(m)
             fbod1(m,n) = fbod1(m,n) - ffe(l) * sin(angfe(l,n))
          enddo
          !
       enddo
       !
       ! compute total spring-damping rotational force acting on the body
       !
       tbod1(m,1) = 0.
       if ( .not.bdof(m,1,2) ) goto 100
       do l = 1, nmlb(m)
          tbod1(m,1) = tbod1(m,1) - fml(l) * ( sin(angml(l,3)) * disml(l,2) - sin(angml(l,2)) * disml(l,3) )
       enddo
       do l = 1, nfeb(m)
          tbod1(m,1) = tbod1(m,1) - ffe(l) * ( sin(angfe(l,3)) * disfe(l,2) - sin(angfe(l,2)) * disfe(l,3) )
       enddo
       !
 100   tbod1(m,2) = 0.
       if ( .not.bdof(m,2,2) ) goto 200
       do l = 1, nmlb(m)
          tbod1(m,2) = tbod1(m,2) + fml(l) * ( sin(angml(l,3)) * disml(l,1) - sin(angml(l,1)) * disml(l,3) )
       enddo
       do l = 1, nfeb(m)
          tbod1(m,2) = tbod1(m,2) + ffe(l) * ( sin(angfe(l,3)) * disfe(l,1) - sin(angfe(l,1)) * disfe(l,3) )
       enddo
       !
 200   tbod1(m,3) = 0.
       if ( .not.bdof(m,3,2) ) goto 300
       do l = 1, nmlb(m)
          tbod1(m,3) = tbod1(m,3) - fml(l) * ( sin(angml(l,2)) * disml(l,1) - sin(angml(l,1)) * disml(l,2) )
       enddo
       do l = 1, nfeb(m)
          tbod1(m,3) = tbod1(m,3) - ffe(l) * ( sin(angfe(l,2)) * disfe(l,1) - sin(angfe(l,1)) * disfe(l,2) )
       enddo
       !
 300   continue
       !
    enddo
    !
    ! add spring forces to restore the equilibrium between buoyancy and gravity
    ! note: they are present for the heave, roll and pitch motions only
    !
    do m = 1, mbod
       !
       if ( bdof(m,3,1) ) fbod1(m,3) = fbod1(m,3) - rhow * grav * barea(m) * xfot1(m,3)
       !
       if ( bdof(m,1,2) ) tbod1(m,1) = tbod1(m,1) - ( rhow * grav * ( bmoax(m) + bvol(m)*bcob(m) ) - grav*bmass(m)*bcog(m,3) ) * xfor1(m,1)
       !
       if ( bdof(m,2,2) ) tbod1(m,2) = tbod1(m,2) - ( rhow * grav * ( bmoay(m) + bvol(m)*bcob(m) ) - grav*bmass(m)*bcog(m,3) ) * xfor1(m,2)
       !
    enddo
    !
    ! add gravity force
    !
    do m = 1, mbod
       !
       if ( bdof(m,3,1) ) fbod1(m,3) = fbod1(m,3) - bmass(m) * grav
       !
    enddo
    !
    ! compute PTO power
    ! note: based on damping coefficient of the PTO unit
    !
    do m = 1, mbod
       !
       ptop(m) = 0.
       !
       do l = 1, nmlb(m)
          !
          ptop(m) = ptop(m) + bmli(m,l,2) * demldt(m,l) * demldt(m,l)
          !
       enddo
       !
    enddo
    !
    contains
    !
    subroutine rotatep ( v, rv, ax, ay, az )
    !
    implicit none
    !
    real,               intent(in)  :: ax ! roll angle
    real,               intent(in)  :: ay ! pitch angle
    real,               intent(in)  :: az ! yaw angle
    real, dimension(3), intent(out) :: rv ! rotated vector
    real, dimension(3), intent(in)  :: v  ! position vector
    !
    real, dimension(3,3)            :: rx ! rotation matrix that rotate vector around the x-axis (counterclockwise)
    real, dimension(3,3)            :: ry ! rotation matrix that rotate vector around the y-axis (counterclockwise)
    real, dimension(3,3)            :: rz ! rotation matrix that rotate vector around the z-axis (counterclockwise)
    !
    ! define rotation matrices using the right-hand rule
    !
    rx(1,1) =  1.
    rx(1,2) =  0.
    rx(1,3) =  0.
    rx(2,1) =  0.
    rx(2,2) =  cos(ax)
    rx(2,3) = -sin(ax)
    rx(3,1) =  0.
    rx(3,2) =  sin(ax)
    rx(3,3) =  cos(ax)
    !
    ry(1,1) =  cos(ay)
    ry(1,2) =  0.
    ry(1,3) =  sin(ay)
    ry(2,1) =  0.
    ry(2,2) =  1.
    ry(2,3) =  0.
    ry(3,1) = -sin(ay)
    ry(3,2) =  0.
    ry(3,3) =  cos(ay)
    !
    rz(1,1) =  cos(az)
    rz(1,2) = -sin(az)
    rz(1,3) =  0.
    rz(2,1) =  sin(az)
    rz(2,2) =  cos(az)
    rz(2,3) =  0.
    rz(3,1) =  0.
    rz(3,2) =  0.
    rz(3,3) =  1.
    !
    ! the order of rotation operations is from right to left, so
    !
    ! rv = rz * ry * rx * v
    !
    rv = matmul( matmul( matmul(rz,ry), rx ), v )
    !
    end subroutine rotatep
    !
end subroutine SwashForcesRigidBod
