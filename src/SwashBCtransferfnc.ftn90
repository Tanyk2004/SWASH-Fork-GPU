subroutine SwashBCtransferfnc ( bcfour, nfreq, xp, yp, ibgrpt, swd, wdir, rsgn, vdir, shape )
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
!    1.00: Panagiotis Vasarmidis
!   10.05: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, November 2023: New subroutine
!   10.05,   August 2024: bug fix: exclude very small omega3
!
!   Purpose
!
!   Computes second order sub- and super-harmonic transfer functions
!   of surface elevation and horizontal layer-averaged velocities
!
!   Method
!
!   Details on the computation of interaction coefficients of two primary wave components can be found in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_bndspec
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)         :: ibgrpt ! actual boundary grid point
    integer, intent(in)         :: nfreq  ! number of frequencies
    integer, intent(in)         :: shape  ! spectral shape
                                          ! = 1; Pierson Moskowitz
                                          ! = 2; Jonswap
                                          ! = 3; TMA
    !
    real, intent(in)            :: rsgn   ! sign for indicating in- and outflowing depending on boundary
                                          ! =+1; refers to inflowing at left and lower boundaries
                                          ! =-1; refers to outflowing at right and upper boundaries
    real, intent(in)            :: swd    ! still water depth
    real, intent(in)            :: wdir   ! incident or peak wave direction with respect to problem coordinates
    real, intent(in)            :: xp     ! x-coordinate of grid point
    real, intent(in)            :: yp     ! y-coordinate of grid point
    !
    logical, intent(in)         :: vdir   ! indicates direction of in- or outcoming velocity on boundary
                                          ! =.true.; u-velocity
                                          ! =.false.; v-velocity
    !
    type(bfsdat), intent(inout) :: bcfour ! list containing parameters for Fourier series
!
!   Parameter variables
!
    real, parameter :: khmin = 0.2  ! minimum dimensionless depth
    real, parameter :: khmax = 5.   ! maximum dimensionless depth
    real, parameter :: odmin = 0.02 ! minimum difference angular frequency
!
!   Local variables
!
    integer, save         :: ient  = 0 ! number of entries in this subroutine
    integer               :: j         ! loop counter
    integer               :: k         ! loop counter
    integer               :: l         ! frequency index of bound wave
    !
    real                  :: ampl1     ! amplitude of first primary wave component
    real                  :: ampl2     ! amplitude of second primary wave component
    real                  :: dw        ! increment in frequency space
    real                  :: eta       ! local 2nd order transfer function of surface elevation
    real                  :: kd1       ! dimensionless depth of first primary wave component
    real                  :: kd2       ! dimensionless depth of second primary wave component
    real                  :: kwav1     ! wave number of first primary wave component
    real                  :: kwav2     ! wave number of second primary wave component
    real                  :: kwav3     ! wave number of bound wave component
    real                  :: n         ! ratio of group and phase velocity
    real                  :: omega1    ! angular frequency of first primary wave component
    real                  :: omega2    ! angular frequency of second primary wave component
    real                  :: omega3    ! angular frequency of bound wave component
    real                  :: phase1    ! phase of first primary wave component
    real                  :: phase2    ! phase of second primary wave component
    real                  :: phase3    ! phase of bound wave component
    real                  :: rval      ! auxiliary real
    real                  :: s         ! sign
    real                  :: theta1    ! wave direction of first primary wave component with respect to computational coordinates
    real                  :: theta2    ! wave direction of second primary wave component with respect to computational coordinates
    real                  :: theta3    ! wave direction of bound wave component with respect to problem coordinates
    real, dimension(kmax) :: vel       ! local 2nd order transfer function of layer-dependent velocity
    real                  :: w0        ! first angular frequency
    real                  :: wn        ! last angular frequency
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashBCtransferfnc')
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
    dw = bcfour%omega(2) - bcfour%omega(1)
    w0 = bcfour%omega(1)
    wn = bcfour%omega(nfreq)
    !
    floop: do j = 1, nfreq
       !
       ! get first primary wave component
       !
       ampl1  = bcfour%ampl (j)
       omega1 = bcfour%omega(j)
       phase1 = bcfour%phase(j)
       !
       theta1 = wdir - alpc + s*bcfour%theta(j)
       !
       ! calculate wave number of this first component
       !
       call disprel ( swd, omega1, kwav1, rval, n )
       !
       kd1 = kwav1 * swd
       !
       ! correct amplitude in case of TMA spectrum for shallow water
       !
       if ( shape == 3 ) ampl1 = ampl1 * omega1 * omega1 / ( grav * kwav1 * sqrt(2.*n) )
       !
       ! in case of periodicity, wave direction must be corrected so that wave number is an integer multiple of 2pi/length with length the periodicity length
       !
       if ( bcperx ) then
          !
          rval = nint( kwav1*cos(theta1) / ( pi2/xclen ) ) * pi2/xclen / kwav1
          if ( rval > 1. ) then
             theta1 = acos ( rval - pi2/xclen / kwav1 )
          else if ( rval < -1. ) then
             theta1 = acos ( rval + pi2/xclen / kwav1 )
          else
             theta1 = acos ( rval )
          endif
          if ( rsgn == -1. ) theta1 = pi2 - theta1
          !
       else if ( bcpery ) then
          !
          rval = nint( kwav1*sin(theta1) / ( pi2/yclen ) ) * pi2/yclen / kwav1
          if ( rval > 1. ) then
             theta1 = asin ( rval - pi2/yclen / kwav1 )
          else if ( rval < -1. ) then
             theta1 = asin ( rval + pi2/yclen / kwav1 )
          else
             theta1 = asin ( rval )
          endif
          if ( rsgn == -1. ) theta1 = pi - theta1
          !
       endif
       !
       ! check this direction with respect to the normal of boundary
       ! (must be within -80 degrees to 80 degrees)
       !
       if ( vdir ) then
          if ( rsgn == 1. ) then
             if ( cos(theta1) <  0.174 ) cycle floop
          else
             if ( cos(theta1) > -0.174 ) cycle floop
          endif
       else
          if ( rsgn == 1. ) then
             if ( sin(theta1) <  0.174 ) cycle floop
          else
             if ( sin(theta1) > -0.174 ) cycle floop
          endif
       endif
       !
       if ( kd1 > khmin .and. kd1 < khmax .and. ampl1 /= 0. ) then
          !
          sloop: do k = j+1, nfreq
             !
             ! get second primary wave component
             !
             ampl2  = bcfour%ampl (k)
             omega2 = bcfour%omega(k)
             phase2 = bcfour%phase(k)
             !
             theta2 = wdir - alpc + s*bcfour%theta(k)
             !
             ! calculate wave number of this second component
             !
             call disprel ( swd, omega2, kwav2, rval, n )
             !
             kd2 = kwav2 * swd
             !
             ! correct amplitude in case of TMA spectrum for shallow water
             !
             if ( shape == 3 ) ampl2 = ampl2 * omega2 * omega2 / ( grav * kwav2 * sqrt(2.*n) )
             !
             ! in case of periodicity, wave direction must be corrected so that wave number is an integer multiple of 2pi/length with length the periodicity length
             !
             if ( bcperx ) then
                !
                rval = nint( kwav2*cos(theta2) / ( pi2/xclen ) ) * pi2/xclen / kwav2
                if ( rval > 1. ) then
                   theta2 = acos ( rval - pi2/xclen / kwav2 )
                else if ( rval < -1. ) then
                   theta2 = acos ( rval + pi2/xclen / kwav2 )
                else
                   theta2 = acos ( rval )
                endif
                if ( rsgn == -1. ) theta2 = pi2 - theta2
                !
             else if ( bcpery ) then
                !
                rval = nint( kwav2*sin(theta2) / ( pi2/yclen ) ) * pi2/yclen / kwav2
                if ( rval > 1. ) then
                   theta2 = asin ( rval - pi2/yclen / kwav2 )
                else if ( rval < -1. ) then
                   theta2 = asin ( rval + pi2/yclen / kwav2 )
                else
                   theta2 = asin ( rval )
                endif
                if ( rsgn == -1. ) theta2 = pi - theta2
                !
             endif
             !
             ! check this direction with respect to the normal of boundary
             ! (must be within -80 degrees to 80 degrees)
             !
             if ( vdir ) then
                if ( rsgn == 1. ) then
                   if ( cos(theta2) <  0.174 ) cycle sloop
                else
                   if ( cos(theta2) > -0.174 ) cycle sloop
                endif
             else
                if ( rsgn == 1. ) then
                   if ( sin(theta2) <  0.174 ) cycle sloop
                else
                   if ( sin(theta2) > -0.174 ) cycle sloop
                endif
             endif
             !
             if ( kd2 > khmin .and. kd2 < khmax .and. ampl2 /= 0. ) then
                !
                ! calculate wave frequency of bound sub-harmonic component
                !
                omega3 = omega2 - omega1
                !
                if ( nfreq == 2 .or. .not. omega3 < odmin ) then
                   !
                   ! calculate frequency index
                   !
                   l = nint(omega3/dw)
                   !
                   ! calculate wave number of bound sub-harmonic component
                   !
                   kwav3 = real(sqrt( dble(kwav1)**2 + dble(kwav2)**2 - 2.*dble(kwav1)*dble(kwav2)*cos(theta1-theta2) ))
                   !
                   ! calculate direction of bound sub-harmonic component
                   !
                   theta3 = alpc + atan( ( kwav2*sin(theta2) - kwav1*sin(theta1) ) / ( kwav2*cos(theta2) - kwav1*cos(theta1) ) )
                   !
                   ! calculate phase of bound sub-harmonic component
                   !
                   phase3 = phase2 - phase1
                   !
                   ! include phase shift related to wave direction and wave number
                   !
                   phase3 = phase3 + kwav3 * ( cos(theta3)*xp + sin(theta3)*yp )
                   !
                   ! compute sub-harmonic transfer functions
                   !
                   if ( kmax == 1 ) then
                      eta    = etasub1()
                      vel(1) = velsb11()
                   else if ( kmax == 2 ) then
                      eta    = etasub2()
                      vel(1) = velsb12()
                      vel(2) = velsb22()
                   else if ( kmax == 3 ) then
                      eta    = etasub3()
                      vel(1) = velsb13()
                      vel(2) = velsb23()
                      vel(3) = velsb33()
                   else if ( kmax == 4 ) then
                      eta    = etasub4()
                      vel(1) = velsb14()
                      vel(2) = velsb24()
                      vel(3) = velsb34()
                      vel(4) = velsb44()
                   endif
                   !
                   ! store surface elevation and velocity of sub-harmonic component
                   !
                   rval = ampl1 * ampl2
                   !
                   subzc(ibgrpt,l) = subzc(ibgrpt,l) + eta * rval * cos( phase3 )
                   subzs(ibgrpt,l) = subzs(ibgrpt,l) + eta * rval * sin( phase3 )
                   !
                   if ( oned .or. optg == 5 ) then
                      subuc(ibgrpt,l,:) = subuc(ibgrpt,l,:) + vel(:) * rval * cos( phase3 )
                      subus(ibgrpt,l,:) = subus(ibgrpt,l,:) + vel(:) * rval * sin( phase3 )
                   else
                      subuc(ibgrpt,l,:) = subuc(ibgrpt,l,:) + vel(:) * rval * cos( theta3 ) * cos( phase3 )
                      subus(ibgrpt,l,:) = subus(ibgrpt,l,:) + vel(:) * rval * cos( theta3 ) * sin( phase3 )
                      subvc(ibgrpt,l,:) = subvc(ibgrpt,l,:) + vel(:) * rval * sin( theta3 ) * cos( phase3 )
                      subvs(ibgrpt,l,:) = subvs(ibgrpt,l,:) + vel(:) * rval * sin( theta3 ) * sin( phase3 )
                   endif
                   !
                endif
                !
                ! calculate wave frequency of bound super-harmonic component
                !
                omega3 = omega1 + omega2
                !
                if ( nfreq == 2 .or. .not. omega3 > wn ) then
                   !
                   ! calculate frequency index
                   !
                   l = nint((omega3-2.*w0)/dw)
                   !
                   ! calculate wave number of bound super-harmonic component
                   !
                   kwav3 = real(sqrt( dble(kwav1)**2 + dble(kwav2)**2 + 2.*dble(kwav1)*dble(kwav2)*cos(theta1-theta2) ))
                   !
                   ! calculate direction of bound super-harmonic component
                   !
                   theta3 = alpc + atan( ( kwav1*sin(theta1) + kwav2*sin(theta2) ) / ( kwav1*cos(theta1) + kwav2*cos(theta2) ) )
                   !
                   ! calculate phase of bound super-harmonic component
                   !
                   phase3 = phase1 + phase2
                   !
                   ! include phase shift related to wave direction and wave number
                   !
                   phase3 = phase3 + kwav3 * ( cos(theta3)*xp + sin(theta3)*yp )
                   !
                   ! compute super-harmonic transfer functions
                   !
                   if ( kmax == 1 ) then
                      eta    = etasup1()
                      vel(1) = velsp11()
                   else if ( kmax == 2 ) then
                      eta    = etasup2()
                      vel(1) = velsp12()
                      vel(2) = velsp22()
                   else if ( kmax == 3 ) then
                      eta    = etasup3()
                      vel(1) = velsp13()
                      vel(2) = velsp23()
                      vel(3) = velsp33()
                   else if ( kmax == 4 ) then
                      eta    = etasup4()
                      vel(1) = velsp14()
                      vel(2) = velsp24()
                      vel(3) = velsp34()
                      vel(4) = velsp44()
                   endif
                   !
                   ! store surface elevation and velocity of super-harmonic component
                   !
                   rval = ampl1 * ampl2
                   !
                   supzc(ibgrpt,l) = supzc(ibgrpt,l) + eta * rval * cos( phase3 )
                   supzs(ibgrpt,l) = supzs(ibgrpt,l) + eta * rval * sin( phase3 )
                   !
                   if ( oned .or. optg == 5 ) then
                      supuc(ibgrpt,l,:) = supuc(ibgrpt,l,:) + vel(:) * rval * cos( phase3 )
                      supus(ibgrpt,l,:) = supus(ibgrpt,l,:) + vel(:) * rval * sin( phase3 )
                   else
                      supuc(ibgrpt,l,:) = supuc(ibgrpt,l,:) + vel(:) * rval * cos( theta3 ) * cos( phase3 )
                      supus(ibgrpt,l,:) = supus(ibgrpt,l,:) + vel(:) * rval * cos( theta3 ) * sin( phase3 )
                      supvc(ibgrpt,l,:) = supvc(ibgrpt,l,:) + vel(:) * rval * sin( theta3 ) * cos( phase3 )
                      supvs(ibgrpt,l,:) = supvs(ibgrpt,l,:) + vel(:) * rval * sin( theta3 ) * sin( phase3 )
                   endif
                   !
                endif
                !
             endif
             !
          enddo sloop
          !
       endif
       !
    enddo floop
    !
    contains
!
real function etasub1()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic surface elevation for 1 layer
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   surface transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'etasub1')
    !
    etasub1 = ((kd1 - kd2)**2 * (4. * grav * kd1 * kd2 + &
              (-(kd2 * (2. * kd1 + kd2) * omega1**2) + &
              2. * (4. + kd1**2 + kd2**2) * omega1 * omega2 - &
              kd1 * (kd1 + 2. * kd2) * omega2**2) * swd)) / &
              (2. * kd1 * kd2 * swd * (-4. * grav * (kd1 - kd2)**2 + &
              (4. + (kd1 - kd2)**2) * (omega1 - omega2)**2 * swd))
    !
end function etasub1
!
real function etasup1()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic surface elevation for 1 layer
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   surface transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'etasup1')
    !
    etasup1 = ((kd1 + kd2)**2 * (4. * grav * kd1 * kd2 + &
              (kd2 * (-2. * kd1 + kd2) * omega1**2 + &
              2. * (4. + kd1**2 + kd2**2) * omega1 * omega2 + &
              kd1 * (kd1 - 2. * kd2) * omega2**2) * swd)) / &
              (2. * kd1 * kd2 * swd * (-4. * grav * (kd1 + kd2)**2 + &
              (4. + (kd1 + kd2)**2) * (omega1 + omega2)**2 * swd))
    !
end function etasup1
!
real function velsb11()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic depth-averaged velocity for 1 layer
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb11')
    !
    velsb11 = (4. * grav * (kd1 - kd2) * ((2. * kd1 - kd2) * kd2 * omega1 + &
              kd1 * (kd1 - 2. * kd2) * omega2) + &
              (omega1 - omega2) * (kd2 * (-4. + 3. * kd1 * (-kd1 + kd2)) * omega1**2 + &
              (kd1 - kd2) * (4. + (kd1 + kd2)**2) * omega1 * omega2 + &
              kd1 * (4. + 3. * kd2 * (-kd1 + kd2)) * omega2**2) * swd) / &
              (2. * kd1 * kd2 * swd * (-4. * grav * (kd1 - kd2)**2 + &
              (4. + (kd1 - kd2)**2) * (omega1 - omega2)**2 * swd))
    !
end function velsb11
!
real function velsp11()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic depth-averaged velocity for 1 layer
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp11')
    !
    velsp11 = (4. * grav * (kd1 + kd2) * (kd2 * (2. * kd1 + kd2) * omega1 + &
              kd1 * (kd1 + 2. * kd2) * omega2) - &
              (omega1 + omega2) * (kd2 * (4. + 3. * kd1 * (kd1 + kd2)) * omega1**2 - &
              (4. + (kd1 - kd2)**2) * (kd1 + kd2) * omega1 * omega2 + &
              kd1 * (4. + 3. * kd2 * (kd1 + kd2)) * omega2**2) * swd) / &
              (2. * kd1 * kd2 * swd * (-4. * grav * (kd1 + kd2)**2 + &
              (4. + (kd1 + kd2)**2) * (omega1 + omega2)**2 * swd))
    !
end function velsp11
!
real function etasub2()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic surface elevation for 2 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   surface transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'etasub2')
    !
    etasub2 = ((kd1 - kd2)**2*(-256.*grav**2*kd1**3*(-32. + (kd1 - kd2)**2)*(16. + (kd1 - kd2)**2)*kd2**3 + 16.*grav*kd1*kd2*(kd2**2*(-32.*(256. - 32.*kd1**2 + 5.*kd1**4) - &
              4.*kd1*(512. - 16.*kd1**2 + kd1**4)*kd2 + (16. + kd1**2)*(-16. + 11.*kd1**2)*kd2**2 - 2.*kd1*(48. + 5.*kd1**2)*kd2**3 + (16. + 3.*kd1**2)*kd2**4)*omega1**2 + &
              kd1*kd2*(kd1**5*kd2 + kd1*(-8. + kd2)*kd2*(8. + kd2)*(48. + kd2**2) + 4.*kd1**3*kd2*(-4. + 3.*kd2**2) - kd1**4*(16. + 7.*kd2**2) + kd1**2*(1536. + 80.*kd2**2 - &
              7.*kd2**4) - 16.*(256. - 96.*kd2**2 + kd2**4))*omega1*omega2 + kd1**2*(kd1**4*(16. + 3.*kd2**2) - 2.*kd1**3*kd2*(48. + 5.*kd2**2) + &
              kd1**2*(16. + kd2**2)*(-16. + 11.*kd2**2) - 4.*kd1*kd2*(512. - 16.*kd2**2 + kd2**4) - 32.*(256. - 32.*kd2**2 + 5.*kd2**4))*omega2**2)*swd + &
              omega1*omega2*(-(kd2**2*(-16. + kd2**2)*(kd1**6 + 2.*kd1**5*kd2 - 256.*(-16. + kd2**2) - 384.*kd1**2*(6. + kd2**2) + 4.*kd1**3*kd2*(184. + kd2**2) - &
              kd1**4*(144. + 7.*kd2**2))*omega1**2) + 2.*kd1*kd2*(kd1**6*(-16. + kd2**2) - 4.*kd1**4*kd2**2*(-8. + kd2**2) + 112.*kd1*kd2*(16. + kd2**2)**2 - &
              16.*(-32. + kd2**2)*(16. + kd2**2)**2 + kd1**5*kd2*(112. + kd2**2) + kd1**3*kd2*(3584. - 48.*kd2**2 + kd2**4) + kd1**2*(12288. - 10240.*kd2**2 + &
              32.*kd2**4 + kd2**6))*omega1*omega2 - kd1**2*(-16. + kd1**2)*(-256.*(-16. + kd1**2) - 384.*(6. + kd1**2)*kd2**2 + 4.*kd1*(184. + kd1**2)*kd2**3 - &
              (144. + 7.*kd1**2)*kd2**4 + 2.*kd1*kd2**5 + kd2**6)*omega2**2)*swd**2))/ &
              (2.*kd1**2*(-16. + kd1**2)*kd2**2*(-16. + kd2**2)*omega1*omega2*swd**2*(-16.*grav*(16. + (kd1 - kd2)**2)*(kd1 - kd2)**2 + (256. + &
              (96. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function etasub2
!
real function etasup2()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic surface elevation for 2 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   surface transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'etasup2')
    !
    etasup2 = ((kd1 + kd2)**2*(-256.*grav**2*kd1**3*kd2**3*(-32. + (kd1 + kd2)**2)*(16. + (kd1 + kd2)**2) - 16.*grav*kd1*kd2*(-(kd2**2*(-32.*(256. - &
              32.*kd1**2 + 5.*kd1**4) + 4.*kd1*(512. - 16.*kd1**2 + kd1**4)*kd2 + (16. + kd1**2)*(-16. + 11.*kd1**2)*kd2**2 + 2.*kd1*(48. + 5.*kd1**2)*kd2**3 + &
              (16. + 3.*kd1**2)*kd2**4)*omega1**2) + kd1*kd2*(kd1**5*kd2 + kd1*(-8. + kd2)*kd2*(8. + kd2)*(48. + kd2**2) + 4.*kd1**3*kd2*(-4. + 3.*kd2**2) + &
              kd1**4*(16. + 7.*kd2**2) + 16.*(256. - 96.*kd2**2 + kd2**4) + kd1**2*(-1536. - 80.*kd2**2 + 7.*kd2**4))*omega1*omega2 - kd1**2*(kd1**4*(16. + &
              3.*kd2**2) + 2.*kd1**3*kd2*(48. + 5.*kd2**2) + kd1**2*(16. + kd2**2)*(-16. + 11.*kd2**2) + 4.*kd1*kd2*(512. - 16.*kd2**2 + kd2**4) - 32.*(256. - &
              32.*kd2**2 + 5.*kd2**4))*omega2**2)*swd + omega1*omega2*(kd2**2*(-16. + kd2**2)*(-kd1**6 + 2.*kd1**5*kd2 + 256.*(-16. + kd2**2) + &
              384.*kd1**2*(6. + kd2**2) + 4.*kd1**3*kd2*(184. + kd2**2) + kd1**4*(144. + 7.*kd2**2))*omega1**2 + 2.*kd1*kd2*(kd1**6*(-16. + kd2**2) - &
              4.*kd1**4*kd2**2*(-8. + kd2**2) - 112.*kd1*kd2*(16. + kd2**2)**2 - 16.*(-32. + kd2**2)*(16. + kd2**2)**2 - kd1**5*kd2*(112. + kd2**2) - &
              kd1**3*kd2*(3584. - 48.*kd2**2 + kd2**4) + kd1**2*(12288. - 10240.*kd2**2 + 32.*kd2**4 + kd2**6))*omega1*omega2 + &
              kd1**2*(-16. + kd1**2)*(256.*(-16. + kd1**2) + 384.*(6. + kd1**2)*kd2**2 + 4.*kd1*(184. + kd1**2)*kd2**3 + (144. + 7.*kd1**2)*kd2**4 + &
              2.*kd1*kd2**5 - kd2**6)*omega2**2)*swd**2))/ &
              (2.*kd1**2*(-16. + kd1**2)*kd2**2*(-16. + kd2**2)*omega1*omega2*swd**2*(-16.*grav*(kd1 + kd2)**2*(16. + (kd1 + kd2)**2) + &
              (256. + (kd1 + kd2)**2*(96. + (kd1 + kd2)**2))*(omega1 + omega2)**2*swd))
    !
end function etasup2
!
real function velsb12()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic velocity of 1st layer for 2 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb12')
    !
    velsb12 = -(-16.*grav*kd1**3*(kd1 - kd2)**3*kd2**3*((kd1 - 4.*kd2)*(3.*kd1 - 2.*kd2)*omega1**4 + (-5.*kd1**2 + 20.*kd1*kd2 - 24.*kd2**2)*omega1**3*omega2 + &
              4.*(5.*kd1**2 - 4.*kd1*kd2 + 5.*kd2**2)*omega1**2*omega2**2 + (-24.*kd1**2 + 20.*kd1*kd2 - 5.*kd2**2)*omega1*omega2**3 + &
              (2.*kd1 - 3.*kd2)*(4.*kd1 - kd2)*omega2**4)*swd - 65536.*grav*(kd1 - kd2)*(-(kd1*kd2**3*omega1**4) + kd2**2*(kd1**2 + kd2**2)*omega1**3*omega2 + &
              2.*kd1*kd2*(kd1**2 - 3.*kd1*kd2 + kd2**2)*omega1**2*omega2**2 + kd1**2*(kd1**2 + kd2**2)*omega1*omega2**3 - kd1**3*kd2*omega2**4)*swd - &
              256.*grav*kd1*(kd1 - kd2)*kd2*(kd2**2*(-19.*kd1**4 + 20.*kd1**3*kd2 + 12.*kd1**2*kd2**2 - 18.*kd1*kd2**3 + 4.*kd2**4)*omega1**4 - &
              kd2*(11.*kd1**5 - 49.*kd1**4*kd2 + 32.*kd1**3*kd2**2 + 18.*kd1**2*kd2**3 - 20.*kd1*kd2**4 + 7.*kd2**5)*omega1**3*omega2 + &
              (kd1 - kd2)**4*(5.*kd1**2 + 18.*kd1*kd2 + 5.*kd2**2)*omega1**2*omega2**2 - kd1*(7.*kd1**5 - 20.*kd1**4*kd2 + 18.*kd1**3*kd2**2 + &
              32.*kd1**2*kd2**3 - 49.*kd1*kd2**4 + 11.*kd2**5)*omega1*omega2**3 + kd1**2*(4.*kd1**4 - 18.*kd1**3*kd2 + 12.*kd1**2*kd2**2 + 20.*kd1*kd2**3 - &
              19.*kd2**4)*omega2**4)*swd + 4096.*grav*(kd1 - kd2)*(kd2**6*omega1**3*omega2 + kd1**6*omega1*omega2**3 - kd1**5*kd2*omega2**2*(5.*omega1**2 - &
              omega1*omega2 + omega2**2) - kd1*kd2**5*omega1**2*(omega1**2 - omega1*omega2 + 5.*omega2**2) + kd1**2*kd2**4*omega1*(8.*omega1**3 - &
              31.*omega1**2*omega2 + 26.*omega1*omega2**2 - 2.*omega2**3) + kd1**4*kd2**2*omega2*(-2.*omega1**3 + 26.*omega1**2*omega2 - 31.*omega1*omega2**2 + &
              8.*omega2**3) - kd1**3*kd2**3*(omega1**4 - 24.*omega1**3*omega2 + 40.*omega1**2*omega2**2 - 24.*omega1*omega2**3 + omega2**4))*swd - &
              65536.*kd1*(kd1 - kd2)*kd2*omega1**2*(omega1 - omega2)**2*omega2**2*swd**2 + &
              kd1**3*(kd1 - kd2)**3*kd2**3*omega1*(omega1 - omega2)**2*omega2*(kd1**2*omega2*(5.*omega1 + 3.*omega2) + kd2**2*omega1*(3.*omega1 + &
              5.*omega2) + 8.*kd1*kd2*(omega1**2 - 4.*omega1*omega2 + omega2**2))*swd**2 + 4096.*(16.*grav**2*kd1*(kd1 - kd2)*kd2*(kd2**2*(-kd1**2 - &
              2.*kd1*kd2 + kd2**2)*omega1**2 + kd1*kd2*(kd1 + kd2)**2*omega1*omega2 + kd1**2*(kd1**2 - 2.*kd1*kd2 - kd2**2)*omega2**2) + &
              omega1*(omega1 - omega2)**2*omega2*(kd2**2*(8.*kd1**3 - 27.*kd1**2*kd2 + 11.*kd1*kd2**2 - 4.*kd2**3)*omega1**2 + &
              kd1*(kd1 - 5.*kd2)*(kd1 - kd2)*(5.*kd1 - kd2)*kd2*omega1*omega2 + kd1**2*(4.*kd1**3 - 11.*kd1**2*kd2 + 27.*kd1*kd2**2 - 8.*kd2**3)*omega2**2)*swd**2) + &
              16.*kd1*(kd1 - kd2)*kd2*(32.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd2**2*(omega1 - 3.*omega2)*(omega1 - omega2) + &
              kd1*kd2*(-4.*omega1**2 + 7.*omega1*omega2 - 4.*omega2**2) + kd1**2*(3.*omega1**2 - 4.*omega1*omega2 + omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(kd1**6*omega2*(7.*omega1 + omega2) + kd2**6*omega1*(omega1 + 7.*omega2) - 2.*kd1**3*kd2**3*(omega1**2 - &
              21.*omega1*omega2 + omega2**2) - 2.*kd1*kd2**5*(41.*omega1**2 + 19.*omega1*omega2 + 4.*omega2**2) + kd1**2*kd2**4*(87.*omega1**2 + &
              omega1*omega2 + 13.*omega2**2) - 2.*kd1**5*kd2*(4.*omega1**2 + 19.*omega1*omega2 + 41.*omega2**2) + kd1**4*kd2**2*(13.*omega1**2 + &
              omega1*omega2 + 87.*omega2**2))*swd**2) + 256.*(16.*grav**2*kd1*(kd1 - kd2)**3*kd2*(kd2**2*(-9.*kd1**2 + 2.*kd1*kd2 + kd2**2)*omega1**2 + &
              kd1*kd2*(-3.*kd1**2 + 11.*kd1*kd2 - 3.*kd2**2)*omega1*omega2 + kd1**2*(kd1**2 + 2.*kd1*kd2 - 9.*kd2**2)*omega2**2) - &
              omega1*(omega1 - omega2)**2*omega2*(-4.*kd2**7*omega1**2 + kd1**2*kd2**5*omega1*(56.*omega1 - 45.*omega2) + &
              kd1**5*kd2**2*(45.*omega1 - 56.*omega2)*omega2 + 4.*kd1**7*omega2**2 + 5.*kd1*kd2**6*omega1*(2.*omega1 + omega2) - &
              5.*kd1**6*kd2*omega2*(omega1 + 2.*omega2) + kd1**3*kd2**4*(-164.*omega1**2 + 113.*omega1*omega2 - 90.*omega2**2) + &
              kd1**4*kd2**3*(90.*omega1**2 - 113.*omega1*omega2 + 164.*omega2**2))*swd**2))/ &
              (2.*kd1**2*(-16. + kd1**2)*kd2**2*(-16. + kd2**2)*omega1*(omega1 - omega2)*omega2*swd**2*(-16.*grav*(16. + (kd1 - kd2)**2)*(kd1 - kd2)**2 + &
              (256. + (96. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function velsb12
!
real function velsp12()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic velocity of 1st layer for 2 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp12')
    !
    velsp12 = -(-65536.*grav*(kd1 + kd2)*(-(kd1*kd2**3*omega1**4) + kd2**2*(kd1**2 + kd2**2)*omega1**3*omega2 + 2.*kd1*kd2*(kd1**2 + 3.*kd1*kd2 + &
              kd2**2)*omega1**2*omega2**2 + kd1**2*(kd1**2 + kd2**2)*omega1*omega2**3 - kd1**3*kd2*omega2**4)*swd - &
              256.*grav*kd1*kd2*(kd1 + kd2)*(kd2**2*(-19.*kd1**4 - 20.*kd1**3*kd2 + 12.*kd1**2*kd2**2 + 18.*kd1*kd2**3 + 4.*kd2**4)*omega1**4 + &
              kd2*(-11.*kd1**5 - 49.*kd1**4*kd2 - 32.*kd1**3*kd2**2 + 18.*kd1**2*kd2**3 + 20.*kd1*kd2**4 + 7.*kd2**5)*omega1**3*omega2 + &
              (kd1 + kd2)**4*(5.*kd1**2 - 18.*kd1*kd2 + 5.*kd2**2)*omega1**2*omega2**2 + kd1*(7.*kd1**5 + 20.*kd1**4*kd2 + 18.*kd1**3*kd2**2 - &
              32.*kd1**2*kd2**3 - 49.*kd1*kd2**4 - 11.*kd2**5)*omega1*omega2**3 + kd1**2*(4.*kd1**4 + 18.*kd1**3*kd2 + 12.*kd1**2*kd2**2 - &
              20.*kd1*kd2**3 - 19.*kd2**4)*omega2**4)*swd + 4096.*grav*(kd1 + kd2)*(kd2**6*omega1**3*omega2 + kd1**6*omega1*omega2**3 - &
              kd1**5*kd2*omega2**2*(5.*omega1**2 + omega1*omega2 + omega2**2) - kd1*kd2**5*omega1**2*(omega1**2 + omega1*omega2 + 5.*omega2**2) - &
              kd1**2*kd2**4*omega1*(8.*omega1**3 + 31.*omega1**2*omega2 + 26.*omega1*omega2**2 + 2.*omega2**3) - kd1**4*kd2**2*omega2*(2.*omega1**3 + &
              26.*omega1**2*omega2 + 31.*omega1*omega2**2 + 8.*omega2**3) - kd1**3*kd2**3*(omega1**4 + 24.*omega1**3*omega2 + 40.*omega1**2*omega2**2 + &
              24.*omega1*omega2**3 + omega2**4))*swd - 16.*grav*kd1**3*kd2**3*(kd1 + kd2)**3*(kd2**2*(2.*omega1 + 3.*omega2)*(4.*omega1**3 + &
              6.*omega1**2*omega2 + omega1*omega2**2 + omega2**3) + kd1**2*(3.*omega1 + 2.*omega2)*(omega1**3 + omega1**2*omega2 + 6.*omega1*omega2**2 + &
              4.*omega2**3) + 2.*kd1*kd2*(7.*omega1**4 + 10.*omega1**3*omega2 + 8.*omega1**2*omega2**2 + 10.*omega1*omega2**3 + 7.*omega2**4))*swd - &
              65536.*kd1*kd2*(kd1 + kd2)*omega1**2*omega2**2*(omega1 + omega2)**2*swd**2 - &
              kd1**3*kd2**3*(kd1 + kd2)**3*omega1*omega2*(omega1 + omega2)**2*(kd2**2*omega1*(3.*omega1 - 5.*omega2) + kd1**2*omega2*(-5.*omega1 + &
              3.*omega2) - 8.*kd1*kd2*(omega1**2 + 4.*omega1*omega2 + omega2**2))*swd**2 + 4096.*(16.*grav**2*kd1*kd2*(kd1 + kd2)*(kd2**4*omega1**2 + &
              kd1**4*omega2**2 - kd1**2*kd2**2*(omega1 + omega2)**2 + kd1*kd2**3*omega1*(2.*omega1 + omega2) + kd1**3*kd2*omega2*(omega1 + 2.*omega2)) + &
              omega1*omega2*(omega1 + omega2)**2*(kd2**2*(8.*kd1**3 + 27.*kd1**2*kd2 + 11.*kd1*kd2**2 + &
              4.*kd2**3)*omega1**2 + kd1*kd2*(kd1 + kd2)*(5.*kd1 + kd2)*(kd1 + 5.*kd2)*omega1*omega2 + kd1**2*(4.*kd1**3 + 11.*kd1**2*kd2 + &
              27.*kd1*kd2**2 + 8.*kd2**3)*omega2**2)*swd**2) + 256.*(16.*grav**2*kd1*kd2*(kd1 + kd2)**3*(kd2**2*(-9.*kd1**2 - 2.*kd1*kd2 + &
              kd2**2)*omega1**2 - kd1*kd2*(3.*kd1**2 + 11.*kd1*kd2 + 3.*kd2**2)*omega1*omega2 + kd1**2*(kd1**2 - 2.*kd1*kd2 - 9.*kd2**2)*omega2**2) + &
              omega1*omega2*(omega1 + omega2)**2*(2.*kd2**3*(45.*kd1**4 + 82.*kd1**3*kd2 + 28.*kd1**2*kd2**2 - 5.*kd1*kd2**3 - 2.*kd2**4)*omega1**2 + &
              kd1*kd2*(kd1 + kd2)*(5.*kd1**4 + 40.*kd1**3*kd2 + 73.*kd1**2*kd2**2 + 40.*kd1*kd2**3 + 5.*kd2**4)*omega1*omega2 + 2.*kd1**3*(-2.*kd1**4 - &
              5.*kd1**3*kd2 + 28.*kd1**2*kd2**2 + 82.*kd1*kd2**3 + 45.*kd2**4)*omega2**2)*swd**2) + &
              16.*kd1*kd2*(kd1 + kd2)*(32.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(kd1**2*(omega1 + omega2)*(3.*omega1 + omega2) + &
              kd2**2*(omega1 + omega2)*(omega1 + 3.*omega2) + kd1*kd2*(4.*omega1**2 + 7.*omega1*omega2 + 4.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(kd2**6*omega1*(omega1 - 7.*omega2) + kd1**6*omega2*(-7.*omega1 + omega2) + &
              2.*kd1**3*kd2**3*(omega1**2 + 21.*omega1*omega2 + omega2**2) + 2.*kd1*kd2**5*(41.*omega1**2 - 19.*omega1*omega2 + 4.*omega2**2) + &
              kd1**2*kd2**4*(87.*omega1**2 - omega1*omega2 + 13.*omega2**2) + 2.*kd1**5*kd2*(4.*omega1**2 - 19.*omega1*omega2 + 41.*omega2**2) + &
              kd1**4*kd2**2*(13.*omega1**2 - omega1*omega2 + 87.*omega2**2))*swd**2))/ &
              (2.*kd1**2*(-16. + kd1**2)*kd2**2*(-16. + kd2**2)*omega1*omega2*(omega1 + omega2)*swd**2*(-16.*grav*(kd1 + kd2)**2*(16. + (kd1 + kd2)**2) + &
              (256. + (kd1 + kd2)**2*(96. + (kd1 + kd2)**2))*(omega1 + omega2)**2*swd))
    !
end function velsp12
!
real function velsb22()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic velocity of 2nd layer for 2 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb22')
    !
    velsb22 = -(16.*grav*kd1**3*(kd1 - kd2)**3*kd2**3*((3.*kd1**2 - 6.*kd1*kd2 + 2.*kd2**2)*omega1**4 + (-7.*kd1**2 + 12.*kd1*kd2 - &
              12.*kd2**2)*omega1**3*omega2 + 16.*(kd1**2 - kd1*kd2 + kd2**2)*omega1**2*omega2**2 + (-12.*kd1**2 + 12.*kd1*kd2 - &
              7.*kd2**2)*omega1*omega2**3 + (2.*kd1**2 - 6.*kd1*kd2 + 3.*kd2**2)*omega2**4)*swd + 65536.*grav*(kd1 - kd2)*(kd2**4*omega1**3*omega2 + &
              kd1**4*omega1*omega2**3 + kd1**2*kd2**2*omega1*omega2*(omega1**2 - 6.*omega1*omega2 + omega2**2) + kd1**3*kd2*omega2**2*(4.*omega1**2 - &
              6.*omega1*omega2 + 3.*omega2**2) + kd1*kd2**3*omega1**2*(3.*omega1**2 - 6.*omega1*omega2 + 4.*omega2**2))*swd + &
              256.*grav*kd1*(kd1 - kd2)*kd2*(kd1**6*omega2**2*(5.*omega1**2 - 5.*omega1*omega2 + 2.*omega2**2) + kd2**6*omega1**2*(2.*omega1**2 - &
              5.*omega1*omega2 + 5.*omega2**2) - kd1**5*kd2*omega2*(7.*omega1**3 + 2.*omega1**2*omega2 - 4.*omega1*omega2**2 + 6.*omega2**3) - &
              kd1*kd2**5*omega1*(6.*omega1**3 - 4.*omega1**2*omega2 + 2.*omega1*omega2**2 + 7.*omega2**3) + kd1**4*kd2**2*(omega1**4 + &
              5.*omega1**3*omega2 - 29.*omega1**2*omega2**2 + 18.*omega1*omega2**3 - 8.*omega2**4) + kd1**2*kd2**4*(-8.*omega1**4 + &
              18.*omega1**3*omega2 - 29.*omega1**2*omega2**2 + 5.*omega1*omega2**3 + omega2**4) + 4.*kd1**3*kd2**3*(3.*omega1**4 - &
              5.*omega1**3*omega2 + 15.*omega1**2*omega2**2 - 5.*omega1*omega2**3 + 3.*omega2**4))*swd - &
              4096.*grav*(kd1 - kd2)*(kd2**6*omega1**3*omega2 + kd1**6*omega1*omega2**3 + kd1*kd2**5*omega1**2*(-3.*omega1**2 + 5.*omega1*omega2 - &
              7.*omega2**2) + kd1**5*kd2*omega2**2*(-7.*omega1**2 + 5.*omega1*omega2 - 3.*omega2**2) + kd1**4*kd2**2*omega2*(10.*omega1**3 - &
              18.*omega1**2*omega2 + 17.*omega1*omega2**2 - 8.*omega2**3) + kd1**2*kd2**4*omega1*(-8.*omega1**3 + 17.*omega1**2*omega2 - &
              18.*omega1*omega2**2 + 10.*omega2**3) + kd1**3*kd2**3*(7.*omega1**4 - 20.*omega1**3*omega2 + 32.*omega1**2*omega2**2 - &
              20.*omega1*omega2**3 + 7.*omega2**4))*swd + 65536.*omega1*(omega1 - omega2)**2*omega2*(2.*kd2**3*omega1**2 + &
              5.*kd1*kd2*(-kd1 + kd2)*omega1*omega2 - 2.*kd1**3*omega2**2)*swd**2 + &
              kd1**3*(kd1 - kd2)**3*kd2**3*omega1*(omega1 - omega2)**2*omega2*(kd2**2*omega1*(3.*omega1 - 7.*omega2) + &
              kd1**2*omega2*(-7.*omega1 + 3.*omega2) - 4.*kd1*kd2*(omega1**2 - 4.*omega1*omega2 + omega2**2))*swd**2 - &
              4096.*kd1*kd2*(16.*grav**2*(kd1 - kd2)*(kd2**4*omega1**2 + 3.*kd1**2*kd2**2*(omega1 - omega2)**2 + kd1**3*kd2*(omega1 - 2.*omega2)*omega2 + &
              kd1**4*omega2**2 + kd1*kd2**3*omega1*(-2.*omega1 + omega2)) + omega1*(omega1 - omega2)**2*omega2*(kd2*(-20.*kd1**2 + 15.*kd1*kd2 - &
              3.*kd2**2)*omega1**2 + (kd1 - kd2)*(7.*kd1**2 + 26.*kd1*kd2 + 7.*kd2**2)*omega1*omega2 + kd1*(3.*kd1**2 - 15.*kd1*kd2 + &
              20.*kd2**2)*omega2**2)*swd**2) - 256.*(16.*grav**2*kd1*(kd1 - kd2)**3*kd2*(kd2**4*omega1**2 + kd1*kd2**3*omega1*(2.*omega1 - 3.*omega2) + &
              kd1**4*omega2**2 + kd1**3*kd2*omega2*(-3.*omega1 + 2.*omega2) - 7.*kd1**2*kd2**2*(omega1**2 - omega1*omega2 + omega2**2)) - &
              omega1*(omega1 - omega2)**2*omega2*(-2.*kd2**2*(-4.*kd1**5 + 2.*kd1**4*kd2 + 26.*kd1**3*kd2**2 - 21.*kd1**2*kd2**3 + kd1*kd2**4 + &
              kd2**5)*omega1**2 + kd1*kd2*(-kd1 + kd2)*(15.*kd1**4 - 77.*kd1**2*kd2**2 + 15.*kd2**4)*omega1*omega2 + 2.*kd1**2*(kd1**5 + kd1**4*kd2 - &
              21.*kd1**3*kd2**2 + 26.*kd1**2*kd2**3 + 2.*kd1*kd2**4 - 4.*kd2**5)*omega2**2)*swd**2) - &
              16.*kd1*(kd1 - kd2)*kd2*(32.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(2.*kd1**2*omega1*(omega1 - omega2) + 2.*kd2**2*omega2*(-omega1 + omega2) + &
              kd1*kd2*(-2.*omega1**2 + 3.*omega1*omega2 - 2.*omega2**2)) - omega1*(omega1 - omega2)**2*omega2*(kd2**6*omega1*(omega1 - 5.*omega2) + &
              kd1**6*omega2*(-5.*omega1 + omega2) + 2.*kd1**3*kd2**3*(omega1**2 - 19.*omega1*omega2 + omega2**2) + 2.*kd1*kd2**5*(14.*omega1**2 + &
              9.*omega1*omega2 + 2.*omega2**2) - kd1**4*kd2**2*(11.*omega1**2 + 11.*omega1*omega2 + 7.*omega2**2) - kd1**2*kd2**4*(7.*omega1**2 + &
              11.*omega1*omega2 + 11.*omega2**2) + 2.*kd1**5*kd2*(2.*omega1**2 + 9.*omega1*omega2 + 14.*omega2**2))*swd**2))/ &
              (2.*kd1**2*(-16. + kd1**2)*kd2**2*(-16. + kd2**2)*omega1*(omega1 - omega2)*omega2*swd**2*(-16.*grav*(16. + (kd1 - kd2)**2)*(kd1 - kd2)**2 + &
              (256. + (96. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function velsb22
!
real function velsp22()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic velocity of 2nd layer for 2 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp22')
    !
    velsp22 = (-16.*grav*kd1**3*kd2**3*(kd1 + kd2)**3*((3.*kd1**2 + 6.*kd1*kd2 + 2.*kd2**2)*omega1**4 + (7.*kd1**2 + 12.*kd1*kd2 + &
              12.*kd2**2)*omega1**3*omega2 + 16.*(kd1**2 + kd1*kd2 + kd2**2)*omega1**2*omega2**2 + (12.*kd1**2 + 12.*kd1*kd2 + &
              7.*kd2**2)*omega1*omega2**3 + (2.*kd1**2 + 6.*kd1*kd2 + 3.*kd2**2)*omega2**4)*swd - 65536.*grav*(kd1 + kd2)*(kd2**4*omega1**3*omega2 + &
              kd1**4*omega1*omega2**3 + kd1**2*kd2**2*omega1*omega2*(omega1**2 + 6.*omega1*omega2 + omega2**2) + kd1**3*kd2*omega2**2*(4.*omega1**2 + &
              6.*omega1*omega2 + 3.*omega2**2) + kd1*kd2**3*omega1**2*(3.*omega1**2 + 6.*omega1*omega2 + 4.*omega2**2))*swd - &
              256.*grav*kd1*kd2*(kd1 + kd2)*(kd1**6*omega2**2*(5.*omega1**2 + 5.*omega1*omega2 + 2.*omega2**2) + kd2**6*omega1**2*(2.*omega1**2 + &
              5.*omega1*omega2 + 5.*omega2**2) + kd1*kd2**5*omega1*(6.*omega1**3 + 4.*omega1**2*omega2 + 2.*omega1*omega2**2 - 7.*omega2**3) + &
              kd1**5*kd2*omega2*(-7.*omega1**3 + 2.*omega1**2*omega2 + 4.*omega1*omega2**2 + 6.*omega2**3) + kd1**4*kd2**2*(omega1**4 - &
              5.*omega1**3*omega2 - 29.*omega1**2*omega2**2 - 18.*omega1*omega2**3 - 8.*omega2**4) + kd1**2*kd2**4*(-8.*omega1**4 - 18.*omega1**3*omega2 - &
              29.*omega1**2*omega2**2 - 5.*omega1*omega2**3 + omega2**4) - 4.*kd1**3*kd2**3*(3.*omega1**4 + 5.*omega1**3*omega2 + 15.*omega1**2*omega2**2 + &
              5.*omega1*omega2**3 + 3.*omega2**4))*swd + 4096.*grav*(kd1 + kd2)*(kd2**6*omega1**3*omega2 + kd1**6*omega1*omega2**3 - &
              kd1**5*kd2*omega2**2*(7.*omega1**2 + 5.*omega1*omega2 + 3.*omega2**2) - kd1*kd2**5*omega1**2*(3.*omega1**2 + 5.*omega1*omega2 + &
              7.*omega2**2) + kd1**4*kd2**2*omega2*(10.*omega1**3 + 18.*omega1**2*omega2 + 17.*omega1*omega2**2 + 8.*omega2**3) + &
              kd1**2*kd2**4*omega1*(8.*omega1**3 + 17.*omega1**2*omega2 + 18.*omega1*omega2**2 + 10.*omega2**3) + kd1**3*kd2**3*(7.*omega1**4 + &
              20.*omega1**3*omega2 + 32.*omega1**2*omega2**2 + 20.*omega1*omega2**3 + 7.*omega2**4))*swd + &
              65536.*omega1*omega2*(omega1 + omega2)**2*(2.*kd2**3*omega1**2 + 5.*kd1*kd2*(kd1 + kd2)*omega1*omega2 + 2.*kd1**3*omega2**2)*swd**2 + &
              kd1**3*kd2**3*(kd1 + kd2)**3*omega1*omega2*(omega1 + omega2)**2*(kd1**2*omega2*(7.*omega1 + 3.*omega2) + &
              kd2**2*omega1*(3.*omega1 + 7.*omega2) + 4.*kd1*kd2*(omega1**2 + 4.*omega1*omega2 + omega2**2))*swd**2 + &
              4096.*kd1*kd2*(16.*grav**2*(kd1 + kd2)*(kd2**4*omega1**2 + kd1**4*omega2**2 + 3.*kd1**2*kd2**2*(omega1 + omega2)**2 + &
              kd1*kd2**3*omega1*(2.*omega1 + omega2) + kd1**3*kd2*omega2*(omega1 + 2.*omega2)) - omega1*omega2*(omega1 + omega2)**2*(kd2*(20.*kd1**2 + &
              15.*kd1*kd2 + 3.*kd2**2)*omega1**2 - (kd1 + kd2)*(7.*kd1**2 - 26.*kd1*kd2 + 7.*kd2**2)*omega1*omega2 + kd1*(3.*kd1**2 + 15.*kd1*kd2 + &
              20.*kd2**2)*omega2**2)*swd**2) + 256.*(16.*grav**2*kd1*kd2*(kd1 + kd2)**3*(kd2**4*omega1**2 + kd1**4*omega2**2 - &
              kd1**3*kd2*omega2*(3.*omega1 + 2.*omega2) - kd1*kd2**3*omega1*(2.*omega1 + 3.*omega2) - 7.*kd1**2*kd2**2*(omega1**2 + omega1*omega2 + &
              omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(2.*kd2**2*(4.*kd1**5 + 2.*kd1**4*kd2 - 26.*kd1**3*kd2**2 - 21.*kd1**2*kd2**3 - &
              kd1*kd2**4 + kd2**5)*omega1**2 - kd1*kd2*(kd1 + kd2)*(15.*kd1**4 - 77.*kd1**2*kd2**2 + 15.*kd2**4)*omega1*omega2 + 2.*kd1**2*(kd1**5 - &
              kd1**4*kd2 - 21.*kd1**3*kd2**2 - 26.*kd1**2*kd2**3 + 2.*kd1*kd2**4 + 4.*kd2**5)*omega2**2)*swd**2) + &
              16.*kd1*kd2*(kd1 + kd2)*(32.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(2.*kd1**2*omega1*(omega1 + omega2) + &
              2.*kd2**2*omega2*(omega1 + omega2) + kd1*kd2*(2.*omega1**2 + 3.*omega1*omega2 + 2.*omega2**2)) + omega1*omega2*(omega1 + &
              omega2)**2*(kd1**6*omega2*(5.*omega1 + omega2) + kd2**6*omega1*(omega1 + 5.*omega2) + kd1**2*kd2**4*(-7.*omega1**2 + &
              11.*omega1*omega2 - 11.*omega2**2) + kd1**4*kd2**2*(-11.*omega1**2 + 11.*omega1*omega2 - 7.*omega2**2) - 2.*kd1**3*kd2**3*(omega1**2 + &
              19.*omega1*omega2 + omega2**2) - 2.*kd1*kd2**5*(14.*omega1**2 - 9.*omega1*omega2 + 2.*omega2**2) - 2.*kd1**5*kd2*(2.*omega1**2 - &
              9.*omega1*omega2 + 14.*omega2**2))*swd**2))/ &
              (2.*kd1**2*(-16. + kd1**2)*kd2**2*(-16. + kd2**2)*omega1*omega2*(omega1 + omega2)*swd**2*(-16.*grav*(kd1 + kd2)**2*(16. + (kd1 + kd2)**2) + &
              (256. + (kd1 + kd2)**2*(96. + (kd1 + kd2)**2))*(omega1 + omega2)**2*swd))
    !
end function velsp22
!
real function etasub3()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic surface elevation for 3 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   surface transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'etasub3')
    !
    etasub3 = ((kd1 - kd2)**2*(2592.*grav**2*kd1**3*(-72. + (kd1 - kd2)**2)*kd2**3*(kd1**6*(36. + kd2**2) + &
              36.*(-6. + kd2)*(6. + kd2)*(36. + kd2**2)**2 - 4.*kd1**5*kd2*(48. + kd2**2) - 192.*kd1*kd2**3*(54. + kd2**2) + &
              6.*kd1**4*(216. + 50.*kd2**2 + kd2**4) - 4.*kd1**3*kd2*(2592. + 72.*kd2**2 + kd2**4) + kd1**2*(-46656. + 15120.*kd2**2 + &
              300.*kd2**4 + kd2**6)) - 12.*grav*kd1*kd2*(108.*(-6. + kd1)*(6. + kd1)*(36. + kd1**2)*(-186624. - 9072.*kd1**2 - &
              504.*kd1**4 + kd1**6)*kd2**2*omega1**2 - 1728.*kd1*(-3219264. - 381024.*kd1**2 - 5076.*kd1**4 + 108.*kd1**6 + kd1**8)*kd2**3*omega1**2 + &
              9.*(40310784. - 66624768.*kd1**2 - 2301696.*kd1**4 + 107136.*kd1**6 + 1128.*kd1**8 + kd1**10)*kd2**4*omega1**2 - &
              48.*kd1*(-5458752. - 104976.*kd1**2 + 17820.*kd1**4 + 567.*kd1**6 + kd1**8)*kd2**5*omega1**2 + 108.*(-279936. + 146448.*kd1**2 - &
              432.*kd1**4 + 405.*kd1**6 + kd1**8)*kd2**6*omega1**2 - 12.*kd1*(419904. - 9072.*kd1**2 + 3996.*kd1**4 + 11.*kd1**6)*kd2**7*omega1**2 + &
              3.*(-93312. + 45360.*kd1**2 + 11520.*kd1**4 + 31.*kd1**6)*kd2**8*omega1**2 - 36.*kd1*(2160. + 384.*kd1**2 + kd1**4)*kd2**9*omega1**2 + &
              6.*(1296. + 360.*kd1**2 + kd1**4)*kd2**10*omega1**2 + kd1*kd2*(-3.*kd1**9*kd2*(-36. + kd2**2) + 18.*kd1**6*(-6. + kd2)*(6. + kd2)*(20. + &
              kd2**2)*(-72. + 13.*kd2**2) + kd1**7*kd2*(-6480. + 1368.*kd2**2 - 161.*kd2**4) + 3888.*(-6. + kd2)*(6. + kd2)*(108. + kd2**2)*(-432. + &
              168.*kd2**2 + kd2**4) + kd1**8*(3888. - 72.*kd2**2 + 47.*kd2**4) + kd1**5*kd2*(-7884864. + 497664.*kd2**2 + 7272.*kd2**4 - 161.*kd2**6) + &
              108.*kd1*(-6. + kd2)*kd2*(6. + kd2)*(-1586304. - 73872.*kd2**2 - 24.*kd2**4 + kd2**6) - 3.*kd1**3*kd2*(-38631168. + 9580032.*kd2**2 - &
              165888.*kd2**4 - 456.*kd2**6 + kd2**8) - 72.*kd1**2*(36951552. + 4712256.*kd2**2 - 296784.*kd2**4 + 2052.*kd2**6 + kd2**8) + &
              kd1**4*(30233088. + 21368448.*kd2**2 - 682992.*kd2**4 - 5040.*kd2**6 + 47.*kd2**8))*omega1*omega2 + &
              3.*kd1**2*(2592.*(-6. + kd1)*(6. + kd1)*(-72. + kd1**2)*(36. + kd1**2)**2 - 5184.*kd1*(-357696. - 16848.*kd1**2 + 324.*kd1**4 + &
              5.*kd1**6)*kd2 + 144.*(2939328. - 1388016.*kd1**2 + 36612.*kd1**4 + 315.*kd1**6 + 5.*kd1**8)*kd2**2 - 576.*kd1*(-381024. - &
              2916.*kd1**2 - 63.*kd1**4 + 8.*kd1**6)*kd2**3 + 2.*(8398080. - 3452544.*kd1**2 - 7776.*kd1**4 + 5760.*kd1**6 + kd1**8)*kd2**4 - &
              12.*kd1*(-243648. + 23760.*kd1**2 + 1332.*kd1**4 + kd1**6)*kd2**5 + (-373248. + 321408.*kd1**2 + 14580.*kd1**4 + 31.*kd1**6)*kd2**6 - &
              4.*kd1*(15552. + 2268.*kd1**2 + 11.*kd1**4)*kd2**7 + 36.*(-504. + 94.*kd1**2 + kd1**4)*kd2**8 - 16.*kd1*(36. + kd1**2)*kd2**9 + &
              3.*(12. + kd1**2)*kd2**10)*omega2**2)*swd - omega1*omega2*(-46656.*kd1**2*(-36. + kd1**2)**4*omega2**2 + &
              3888.*(-6. + kd1)*kd1*(6. + kd1)*kd2*omega2*(-((-72. + kd1**2)*(36. + kd1**2)*(1296. + 216.*kd1**2 + kd1**4)*omega1) + &
              24.*kd1**2*(-36. + kd1**2)**2*omega2) + 3.*kd1*kd2**11*omega1*(2.*(-6. + kd1)*kd1**2*(6. + kd1)*omega1 - &
              3.*(12. + kd1**2)*(36. + kd1**2)*omega2) + kd2**10*(-((46656. - 24624.*kd1**2 + 3276.*kd1**4 + 23.*kd1**6)*omega1**2) + 8.*kd1**2*(6480. + &
              504.*kd1**2 + 5.*kd1**4)*omega1*omega2 + kd1**2*(-36. + kd1**2)**2*omega2**2) + kd1*kd2**9*(8.*(11664. - 25272.*kd1**2 + 1701.*kd1**4 + &
              4.*kd1**6)*omega1**2 - (559872. + 133488.*kd1**2 + 18072.*kd1**4 + 79.*kd1**6)*omega1*omega2 + 2.*kd1**2*(-36. + kd1**2)**2*omega2**2) - &
              6.*kd2**8*(3.*(-6. + kd1)*(6. + kd1)*(10368. + 8496.*kd1**2 + 952.*kd1**4 + kd1**6)*omega1**2 - 16.*kd1**2*(124416. + 8640.*kd1**2 + &
              402.*kd1**4 + kd1**6)*omega1*omega2 + 3.*kd1**2*(-36. + kd1**2)**2*(42. + kd1**2)*omega2**2) + kd1*kd2**7*(2.*(-5038848. + 606528.*kd1**2 - &
              454896.*kd1**4 + 3024.*kd1**6 + kd1**8)*omega1**2 - (-60466176. + 34898688.*kd1**2 + 2494800.*kd1**4 + 48240.*kd1**6 + &
              79.*kd1**8)*omega1*omega2 + 32.*kd1**2*(-36. + kd1**2)**2*(261. + kd1**2)*omega2**2) - 216.*kd1*kd2**3*(-12.*(-6. + kd1)*(6. + kd1)*(46656. + &
              198288.*kd1**2 + 4212.*kd1**4 + kd1**6)*omega1**2 + 2.*(166281984. - 45349632.*kd1**2 - 917568.*kd1**4 + 80784.*kd1**6 + &
              309.*kd1**8 + kd1**10)*omega1*omega2 + kd1**2*(-36. + kd1**2)**2*(65664. + 1008.*kd1**2 + kd1**4)*omega2**2) - 3.*kd1*kd2**5*(48.*(-2519424. - &
              5132160.*kd1**2 - 79704.*kd1**4 + 3888.*kd1**6 + kd1**8)*omega1**2 + 3.*(-120932352. - 44043264.*kd1**2 - 5961600.*kd1**4 + 277200.*kd1**6 + &
              2008.*kd1**8 + kd1**10)*omega1*omega2 - 2.*kd1**2*(-36. + kd1**2)**2*(15552. + 2340.*kd1**2 + kd1**4)*omega2**2) + kd2**6*((-362797056. - &
              443418624.*kd1**2 - 17402688.*kd1**4 + 1185840.*kd1**6 + 540.*kd1**8 + kd1**10)*omega1**2 + 8.*kd1**2*(-24354432. + 233280.*kd1**2 + &
              428328.*kd1**4 + 4824.*kd1**6 + 5.*kd1**8)*omega1*omega2 - kd1**2*(-36. + kd1**2)**2*(90720. + 18144.*kd1**2 + 23.*kd1**4)*omega2**2) - &
              36.*kd2**4*(2.*(-120932352. - 92378880.*kd1**2 - 839808.*kd1**4 + 235872.*kd1**6 - 432.*kd1**8 + kd1**10)*omega1**2 - 16.*kd1**2*(-31072896. - &
              2612736.*kd1**2 + 3240.*kd1**4 + 1440.*kd1**6 + 7.*kd1**8)*omega1*omega2 + kd1**2*(-36. + kd1**2)**2*(-90720. - 6336.*kd1**2 + &
              91.*kd1**4)*omega2**2) + 1296.*kd2**2*((-6. + kd1)*(6. + kd1)*(1679616. - 933120.*kd1**2 - 116640.*kd1**4 - 720.*kd1**6 + kd1**8)*omega1**2 + &
              8.*kd1**2*(-9237888. - 1726272.*kd1**2 - 18792.*kd1**4 + 1152.*kd1**6 + 5.*kd1**8)*omega1*omega2 + kd1**2*(-36. + kd1**2)**2*(27216. + &
              19.*kd1**2*(288. + kd1**2))*omega2**2))*swd**2))/(2.*kd1**2*(-36. + kd1**2)**2*kd2**2*(-36. + kd2**2)**2*omega1*omega2*swd**2*(-36.*grav*(12. + &
              (kd1 - kd2)**2)*(108. + (kd1 - kd2)**2)*(kd1 - kd2)**2 + (36. + (kd1 - kd2)**2)*(1296. + (504. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function etasub3
!
real function etasup3()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic surface elevation for 3 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   surface transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'etasup3')
    !
    etasup3 = ((kd1 + kd2)**2*(2592.*grav**2*kd1**3*kd2**3*(-72. + (kd1 + kd2)**2)*(kd1**6*(36. + kd2**2) + 36.*(-36. + kd2**2)*(36. + kd2**2)**2 + &
              4.*kd1**5*kd2*(48. + kd2**2) + 192.*kd1*kd2**3*(54. + kd2**2) + 6.*kd1**4*(216. + 50.*kd2**2 + kd2**4) + 4.*kd1**3*kd2*(2592. + &
              72.*kd2**2 + kd2**4) + kd1**2*(-46656. + 15120.*kd2**2 + 300.*kd2**4 + kd2**6)) - &
              12.*grav*kd1*kd2*(108.*(-6. + kd1)*(6. + kd1)*(36. + kd1**2)*(-186624. - 9072.*kd1**2 - 504.*kd1**4 + kd1**6)*kd2**2*omega1**2 + &
              1728.*kd1*(-3219264. - 381024.*kd1**2 - 5076.*kd1**4 + 108.*kd1**6 + kd1**8)*kd2**3*omega1**2 + 9.*(40310784. - 66624768.*kd1**2 - &
              2301696.*kd1**4 + 107136.*kd1**6 + 1128.*kd1**8 + kd1**10)*kd2**4*omega1**2 + 48.*kd1*(-5458752. - 104976.*kd1**2 + 17820.*kd1**4 + &
              567.*kd1**6 + kd1**8)*kd2**5*omega1**2 + 108.*(-279936. + 146448.*kd1**2 - 432.*kd1**4 + 405.*kd1**6 + kd1**8)*kd2**6*omega1**2 + &
              12.*kd1*(419904. - 9072.*kd1**2 + 3996.*kd1**4 + 11.*kd1**6)*kd2**7*omega1**2 + 3.*(-93312. + 45360.*kd1**2 + &
              11520.*kd1**4 + 31.*kd1**6)*kd2**8*omega1**2 + 36.*kd1*(2160. + 384.*kd1**2 + kd1**4)*kd2**9*omega1**2 + 6.*(1296. + &
              360.*kd1**2 + kd1**4)*kd2**10*omega1**2 + kd1*kd2*(3.*kd1**9*kd2*(-36. + kd2**2) + &
              18.*kd1**6*(-6. + kd2)*(6. + kd2)*(20. + kd2**2)*(-72. + 13.*kd2**2) + &
              3888.*(-6. + kd2)*(6. + kd2)*(108. + kd2**2)*(-432. + 168.*kd2**2 + kd2**4) + kd1**8*(3888. - 72.*kd2**2 + 47.*kd2**4) + &
              kd1**7*kd2*(6480. - 1368.*kd2**2 + 161.*kd2**4) - 108.*kd1*(-6. + kd2)*kd2*(6. + kd2)*(-1586304. - 73872.*kd2**2 - 24.*kd2**4 + kd2**6) + &
              kd1**5*kd2*(7884864. - 497664.*kd2**2 - 7272.*kd2**4 + 161.*kd2**6) + 3.*kd1**3*kd2*(-38631168. + &
              9580032.*kd2**2 - 165888.*kd2**4 - 456.*kd2**6 + kd2**8) - 72.*kd1**2*(36951552. + 4712256.*kd2**2 - 296784.*kd2**4 + 2052.*kd2**6 + kd2**8) + &
              kd1**4*(30233088. + 21368448.*kd2**2 - 682992.*kd2**4 - 5040.*kd2**6 + 47.*kd2**8))*omega1*omega2 + &
              3.*kd1**2*(2592.*(-6. + kd1)*(6. + kd1)*(-72. + kd1**2)*(36. + kd1**2)**2 + 5184.*kd1*(-357696. - 16848.*kd1**2 + 324.*kd1**4 + 5.*kd1**6)*kd2 + &
              144.*(2939328. - 1388016.*kd1**2 + 36612.*kd1**4 + 315.*kd1**6 + 5.*kd1**8)*kd2**2 + 576.*kd1*(-381024. - 2916.*kd1**2 - 63.*kd1**4 + &
              8.*kd1**6)*kd2**3 + 2.*(8398080. - 3452544.*kd1**2 - 7776.*kd1**4 + 5760.*kd1**6 + kd1**8)*kd2**4 + 12.*kd1*(-243648. + 23760.*kd1**2 + &
              1332.*kd1**4 + kd1**6)*kd2**5 + (-373248. + 321408.*kd1**2 + 14580.*kd1**4 + 31.*kd1**6)*kd2**6 + 4.*kd1*(15552. + 2268.*kd1**2 + 11.*kd1**4)*kd2**7 + &
              36.*(-504. + 94.*kd1**2 + kd1**4)*kd2**8 + 16.*kd1*(36. + kd1**2)*kd2**9 + 3.*(12. + kd1**2)*kd2**10)*omega2**2)*swd + &
              omega1*omega2*(46656.*kd1**2*(-36. + kd1**2)**4*omega2**2 + 3888.*(-6. + kd1)*kd1*(6. + kd1)*kd2*omega2*((-72. + kd1**2)*(36. + kd1**2)*(1296. + &
              216.*kd1**2 + kd1**4)*omega1 + 24.*kd1**2*(-36. + kd1**2)**2*omega2) + 3.*kd1*kd2**11*omega1*(2.*(-6. + kd1)*kd1**2*(6. + kd1)*omega1 + &
              3.*(12. + kd1**2)*(36. + kd1**2)*omega2) + kd2**10*((46656. - 24624.*kd1**2 + 3276.*kd1**4 + 23.*kd1**6)*omega1**2 + 8.*kd1**2*(6480. + &
              504.*kd1**2 + 5.*kd1**4)*omega1*omega2 - kd1**2*(-36. + kd1**2)**2*omega2**2) + kd1*kd2**9*(8.*(11664. - 25272.*kd1**2 + 1701.*kd1**4 + &
              4.*kd1**6)*omega1**2 + (559872. + 133488.*kd1**2 + 18072.*kd1**4 + 79.*kd1**6)*omega1*omega2 + 2.*kd1**2*(-36. + kd1**2)**2*omega2**2) + &
              6.*kd2**8*(3.*(-6. + kd1)*(6. + kd1)*(10368. + 8496.*kd1**2 + 952.*kd1**4 + kd1**6)*omega1**2 + 16.*kd1**2*(124416. + 8640.*kd1**2 + &
              402.*kd1**4 + kd1**6)*omega1*omega2 + 3.*kd1**2*(-36. + kd1**2)**2*(42. + kd1**2)*omega2**2) + kd1*kd2**7*(2.*(-5038848. + 606528.*kd1**2 - &
              454896.*kd1**4 + 3024.*kd1**6 + kd1**8)*omega1**2 + (-60466176. + 34898688.*kd1**2 + 2494800.*kd1**4 + 48240.*kd1**6 + 79.*kd1**8)*omega1*omega2 + &
              32.*kd1**2*(-36. + kd1**2)**2*(261. + kd1**2)*omega2**2) + 216.*kd1*kd2**3*(12.*(-6. + kd1)*(6. + kd1)*(46656. + 198288.*kd1**2 + 4212.*kd1**4 + &
              kd1**6)*omega1**2 + 2.*(166281984. - 45349632.*kd1**2 - 917568.*kd1**4 + 80784.*kd1**6 + 309.*kd1**8 + kd1**10)*omega1*omega2 - kd1**2*(-36. + &
              kd1**2)**2*(65664. + 1008.*kd1**2 + kd1**4)*omega2**2) + 3.*kd1*kd2**5*(-48.*(-2519424. - 5132160.*kd1**2 - 79704.*kd1**4 + 3888.*kd1**6 + &
              kd1**8)*omega1**2 + 3.*(-120932352. - 44043264.*kd1**2 - 5961600.*kd1**4 + 277200.*kd1**6 + 2008.*kd1**8 + kd1**10)*omega1*omega2 + &
              2.*kd1**2*(-36. + kd1**2)**2*(15552. + 2340.*kd1**2 + kd1**4)*omega2**2) + kd2**6*(-((-362797056. - 443418624.*kd1**2 - 17402688.*kd1**4 + &
              1185840.*kd1**6 + 540.*kd1**8 + kd1**10)*omega1**2) + 8.*kd1**2*(-24354432. + 233280.*kd1**2 + 428328.*kd1**4 + 4824.*kd1**6 + &
              5.*kd1**8)*omega1*omega2 + kd1**2*(-36. + kd1**2)**2*(90720. + 18144.*kd1**2 + 23.*kd1**4)*omega2**2) + 36.*kd2**4*(2.*(-120932352. - &
              92378880.*kd1**2 - 839808.*kd1**4 + 235872.*kd1**6 - 432.*kd1**8 + kd1**10)*omega1**2 + 16.*kd1**2*(-31072896. - 2612736.*kd1**2 + &
              3240.*kd1**4 + 1440.*kd1**6 + 7.*kd1**8)*omega1*omega2 + kd1**2*(-36. + kd1**2)**2*(-90720. - 6336.*kd1**2 + 91.*kd1**4)*omega2**2) - &
              1296.*kd2**2*((-6. + kd1)*(6. + kd1)*(1679616. - 933120.*kd1**2 - 116640.*kd1**4 - 720.*kd1**6 + kd1**8)*omega1**2 + 8.*kd1**2*(9237888. + &
              1726272.*kd1**2 + 18792.*kd1**4 - 1152.*kd1**6 - 5.*kd1**8)*omega1*omega2 + kd1**2*(-36. + kd1**2)**2*(27216. + 19.*kd1**2*(288. + &
              kd1**2))*omega2**2))*swd**2))/(2.*kd1**2*(-36. + kd1**2)**2*kd2**2*(-36. + kd2**2)**2*omega1*omega2*swd**2*(-36.*grav*(kd1 + kd2)**2*(12. + &
              (kd1 + kd2)**2)*(108. + (kd1 + kd2)**2) + (36. + (kd1 + kd2)**2)*(1296. + (kd1 + kd2)**2*(504. + (kd1 + kd2)**2))*(omega1 + omega2)**2*swd))
    !
end function etasup3
!
real function velsb13()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic velocity of 1st layer for 3 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb13')
    !
    velsb13 = (3.*(13060694016.*grav*(kd1 - kd2)**3*omega1*omega2*(kd2**2*omega1**2 + 6.*kd1*kd2*omega1*omega2 + kd1**2*omega2**2)*swd**2 + &
              6.*grav*kd1**5*(kd1 - kd2)**5*kd2**5*(-2.*(8.*kd1**2 - 14.*kd1*kd2 + 9.*kd2**2)*omega1**4 + 2.*(18.*kd1**2 - 55.*kd1*kd2 + &
              11.*kd2**2)*omega1**3*omega2 + (-27.*kd1**2 + 170.*kd1*kd2 - 27.*kd2**2)*omega1**2*omega2**2 + 2.*(11.*kd1**2 - 55.*kd1*kd2 + &
              18.*kd2**2)*omega1*omega2**3 - 2.*(9.*kd1**2 - 14.*kd1*kd2 + 8.*kd2**2)*omega2**4)*swd**2 - &
              kd1**5*(kd1 - kd2)**5*kd2**5*omega1*(omega1 - omega2)**2*omega2*(2.*kd2*(kd1 + kd2)*omega1**2 + (-5.*kd1**2 + 2.*kd1*kd2 - &
              5.*kd2**2)*omega1*omega2 + 2.*kd1*(kd1 + kd2)*omega2**2)*swd**3 - 2176782336.*swd*(6.*grav**2*kd1*(kd1 - kd2)*kd2*(4.*kd2**4*omega1**2 + &
              kd1**3*kd2*(omega1 - 8.*omega2)*omega2 + 4.*kd1**4*omega2**2 + kd1*kd2**3*omega1*(-8.*omega1 + omega2) + &
              3.*kd1**2*kd2**2*(omega1**2 + omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(-(kd2**2*(-2.*kd1**3 + &
              9.*kd1**2*kd2 - 3.*kd1*kd2**2 + kd2**3)*omega1**2) + 3.*kd1*(kd1 - 2.*kd2)*(kd1 - kd2)*(2.*kd1 - kd2)*kd2*omega1*omega2 + &
              kd1**2*(kd1**3 - 3.*kd1**2*kd2 + 9.*kd1*kd2**2 - 2.*kd2**3)*omega2**2)*swd**2) + &
              60466176.*swd*(2.*grav**2*kd1*(kd1 - kd2)**2*kd2*(kd2**2*(-81.*kd1**3 + 55.*kd1**2*kd2 + 20.*kd1*kd2**2 + 18.*kd2**3)*omega1**2 + &
              5.*kd1*(kd1 - kd2)*kd2*(7.*kd1**2 + 3.*kd1*kd2 + 7.*kd2**2)*omega1*omega2 - kd1**2*(18.*kd1**3 + 20.*kd1**2*kd2 + 55.*kd1*kd2**2 - &
              81.*kd2**3)*omega2**2) + omega1*(omega1 - omega2)**2*omega2*(kd2**2*(4.*kd1**5 + 71.*kd1**4*kd2 - 151.*kd1**3*kd2**2 + 55.*kd1**2*kd2**3 + &
              9.*kd1*kd2**4 - 3.*kd2**5)*omega1**2 + 3.*kd1*kd2*(-kd1 + kd2)*(13.*kd1**4 - 33.*kd1**3*kd2 + 45.*kd1**2*kd2**2 - 33.*kd1*kd2**3 + &
              13.*kd2**4)*omega1*omega2 + kd1**2*(3.*kd1**5 - 9.*kd1**4*kd2 - 55.*kd1**3*kd2**2 + 151.*kd1**2*kd2**3 - 71.*kd1*kd2**4 - &
              4.*kd2**5)*omega2**2)*swd**2) + 362797056.*grav*(108.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3 + (2.*kd1*kd2**3*(15.*kd1**3 - &
              46.*kd1**2*kd2 + 41.*kd1*kd2**2 - 12.*kd2**3)*omega1**4 + kd2**2*(-6.*kd1**5 - 85.*kd1**4*kd2 + 236.*kd1**3*kd2**2 - 179.*kd1**2*kd2**3 + &
              39.*kd1*kd2**4 + 3.*kd2**5)*omega1**3*omega2 + kd1*(kd1 - kd2)*kd2*(63.*kd1**4 - 138.*kd1**3*kd2 + 172.*kd1**2*kd2**2 - 138.*kd1*kd2**3 + &
              63.*kd2**4)*omega1**2*omega2**2 + kd1**2*(-3.*kd1**5 - 39.*kd1**4*kd2 + 179.*kd1**3*kd2**2 - 236.*kd1**2*kd2**3 + 85.*kd1*kd2**4 + &
              6.*kd2**5)*omega1*omega2**3 + 2.*kd1**3*kd2*(12.*kd1**3 - 41.*kd1**2*kd2 + 46.*kd1*kd2**2 - 15.*kd2**3)*omega2**4)*swd**2) + &
              36.*kd1**3*(kd1 - kd2)**3*kd2**3*swd*(2.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(3.*kd2**2*(21.*omega1**2 - 25.*omega1*omega2 + 14.*omega2**2) - &
              2.*kd1*kd2*(51.*omega1**2 - 82.*omega1*omega2 + 51.*omega2**2) + kd1**2*(42.*omega1**2 - 75.*omega1*omega2 + 63.*omega2**2)) - &
              omega1*(omega1 - omega2)**2*omega2*(2.*kd2**6*omega1*(omega1 - 10.*omega2) + 378.*kd1**3*kd2**3*omega1*omega2 + &
              2.*kd1**6*omega2*(-10.*omega1 + omega2) - 2.*kd1*kd2**5*(57.*omega1**2 - 56.*omega1*omega2 + 2.*omega2**2) + kd1**2*kd2**4*(131.*omega1**2 - &
              300.*omega1*omega2 + 4.*omega2**2) - 2.*kd1**5*kd2*(2.*omega1**2 - 56.*omega1*omega2 + 57.*omega2**2) + kd1**4*kd2**2*(4.*omega1**2 - &
              300.*omega1*omega2 + 131.*omega2**2))*swd**2) - 1679616.*swd*(-2.*grav**2*kd1**2*(kd1 - kd2)*kd2**2*(kd1**6*(3.*omega1 - 82.*omega2)*omega2 + &
              kd2**6*omega1*(-82.*omega1 + 3.*omega2) + kd1*kd2**5*(114.*omega1**2 + 179.*omega1*omega2 - 213.*omega2**2) + 2.*kd1**4*kd2**2*(227.*omega1**2 - &
              251.*omega1*omega2 + 35.*omega2**2) + kd1**5*kd2*(-213.*omega1**2 + 179.*omega1*omega2 + 114.*omega2**2) - 2.*kd1**3*kd2**3*(167.*omega1**2 - &
              311.*omega1*omega2 + 167.*omega2**2) + 2.*kd1**2*kd2**4*(35.*omega1**2 - 251.*omega1*omega2 + 227.*omega2**2)) - omega1*(omega1 - &
              omega2)**2*omega2*(3.*kd2**9*omega1**2 - 3.*kd1**9*omega2**2 + 3.*kd1**8*kd2*omega2*(-7.*omega1 + 3.*omega2) + &
              3.*kd1*kd2**8*omega1*(-3.*omega1 + 7.*omega2) - 4.*kd1**2*kd2**7*(31.*omega1**2 - 5.*omega1*omega2 + 15.*omega2**2) + &
              4.*kd1**7*kd2**2*(15.*omega1**2 - 5.*omega1*omega2 + 31.*omega2**2) + kd1**5*kd2**4*(317.*omega1**2 + 507.*omega1*omega2 + 239.*omega2**2) + &
              kd1**3*kd2**6*(268.*omega1**2 + 71.*omega1*omega2 + 261.*omega2**2) - kd1**6*kd2**3*(261.*omega1**2 + 71.*omega1*omega2 + 268.*omega2**2) - &
              kd1**4*kd2**5*(239.*omega1**2 + 507.*omega1*omega2 + 317.*omega2**2))*swd**2) - &
              46656.*swd*(2.*grav**2*kd1*(kd1 - kd2)**2*kd2*(6.*kd2**9*omega1**2 - 6.*kd1**9*omega2**2 - 3.*kd1*kd2**8*omega1*(10.*omega1 + omega2) + &
              3.*kd1**8*kd2*omega2*(omega1 + 10.*omega2) + kd1**2*kd2**7*(-87.*omega1**2 + 82.*omega1*omega2 - 39.*omega2**2) + &
              kd1**6*kd2**3*(203.*omega1**2 + 8.*omega1*omega2 + 37.*omega2**2) + kd1**7*kd2**2*(39.*omega1**2 - 82.*omega1*omega2 + 87.*omega2**2) - &
              kd1**3*kd2**6*(37.*omega1**2 + 8.*omega1*omega2 + 203.*omega2**2) - 2.*kd1**5*kd2**4*(378.*omega1**2 - 152.*omega1*omega2 + 325.*omega2**2) + &
              2.*kd1**4*kd2**5*(325.*omega1**2 - 152.*omega1*omega2 + 378.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(kd2**11*omega1**2 - &
              3.*kd1*kd2**10*omega1*(omega1 - 5.*omega2) - kd1**11*omega2**2 + 3.*kd1**10*kd2*omega2*(-5.*omega1 + omega2) + kd1**5*kd2**6*(163.*omega1**2 - &
              764.*omega1*omega2 - 231.*omega2**2) + kd1**6*kd2**5*(231.*omega1**2 + 764.*omega1*omega2 - 163.*omega2**2) + kd1**3*kd2**8*(73.*omega1**2 + &
              596.*omega1*omega2 - 147.*omega2**2) + kd1**8*kd2**3*(147.*omega1**2 - 596.*omega1*omega2 - 73.*omega2**2) + kd1**2*kd2**9*(-47.*omega1**2 - &
              214.*omega1*omega2 + 4.*omega2**2) + kd1**9*kd2**2*(-4.*omega1**2 + 214.*omega1*omega2 + 47.*omega2**2) + kd1**7*kd2**4*(-411.*omega1**2 + &
              290.*omega1*omega2 + 145.*omega2**2) + kd1**4*kd2**7*(-145.*omega1**2 - 290.*omega1*omega2 + 411.*omega2**2))*swd**2) + &
              1296.*kd1*(kd1 - kd2)*kd2*swd*(2.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd1**5*kd2*(-246.*omega1**2 + 433.*omega1*omega2 - 420.*omega2**2) + &
              kd1*kd2**5*(-420.*omega1**2 + 433.*omega1*omega2 - 246.*omega2**2) + 3.*kd2**6*(25.*omega1**2 - 27.*omega1*omega2 + 18.*omega2**2) + &
              kd1**6*(54.*omega1**2 - 81.*omega1*omega2 + 75.*omega2**2) - 2.*kd1**3*kd2**3*(299.*omega1**2 - 272.*omega1*omega2 + 299.*omega2**2) + &
              kd1**2*kd2**4*(720.*omega1**2 - 621.*omega1*omega2 + 418.*omega2**2) + kd1**4*kd2**2*(418.*omega1**2 - 621.*omega1*omega2 + 720.*omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(3.*kd1**10*omega1*omega2 + 3.*kd2**10*omega1*omega2 - kd1*kd2**9*(13.*omega1**2 + 34.*omega1*omega2 + &
              2.*omega2**2) + kd1**2*kd2**8*(31.*omega1**2 + 188.*omega1*omega2 + 6.*omega2**2) - kd1**3*kd2**7*(112.*omega1**2 + 949.*omega1*omega2 + &
              12.*omega2**2) - kd1**9*kd2*(2.*omega1**2 + 34.*omega1*omega2 + 13.*omega2**2) + kd1**8*kd2**2*(6.*omega1**2 + 188.*omega1*omega2 + &
              31.*omega2**2) - kd1**7*kd2**3*(12.*omega1**2 + 949.*omega1*omega2 + 112.*omega2**2) + kd1**4*kd2**6*(527.*omega1**2 + 2527.*omega1*omega2 + &
              282.*omega2**2) - 2.*kd1**5*kd2**5*(355.*omega1**2 + 1732.*omega1*omega2 + 355.*omega2**2) + kd1**6*kd2**4*(282.*omega1**2 + &
              2527.*omega1*omega2 + 527.*omega2**2))*swd**2) + 279936.*grav*(12.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(7.*kd1**4 + 26.*kd1**3*kd2 - &
              47.*kd1**2*kd2**2 + 26.*kd1*kd2**3 + 7.*kd2**4) + (kd2**11*omega1**3*omega2 - 3.*kd1*kd2**10*omega1**2*(omega1 - 5.*omega2)*omega2 - &
              kd1**11*omega1*omega2**3 + 3.*kd1**10*kd2*omega1*omega2**2*(-5.*omega1 + omega2) + kd1**2*kd2**9*omega1*(-194.*omega1**3 + &
              357.*omega1**2*omega2 - 405.*omega1*omega2**2 + 10.*omega2**3) + kd1**9*kd2**2*omega2*(-10.*omega1**3 + 405.*omega1**2*omega2 - &
              357.*omega1*omega2**2 + 194.*omega2**3) + kd1**5*kd2**6*(-370.*omega1**4 + 679.*omega1**3*omega2 - 458.*omega1**2*omega2**2 + &
              545.*omega1*omega2**3 - 636.*omega2**4) + kd1**8*kd2**3*(126.*omega1**4 + 29.*omega1**3*omega2 - 1190.*omega1**2*omega2**2 + &
              1032.*omega1*omega2**3 - 574.*omega2**4) + kd1**3*kd2**8*(574.*omega1**4 - 1032.*omega1**3*omega2 + 1190.*omega1**2*omega2**2 - &
              29.*omega1*omega2**3 - 126.*omega2**4) + kd1**7*kd2**4*(-412.*omega1**4 + 17.*omega1**3*omega2 + 926.*omega1**2*omega2**2 - &
              515.*omega1*omega2**3 + 364.*omega2**4) + kd1**6*kd2**5*(636.*omega1**4 - 545.*omega1**3*omega2 + 458.*omega1**2*omega2**2 - &
              679.*omega1*omega2**3 + 370.*omega2**4) + kd1**4*kd2**7*(-364.*omega1**4 + 515.*omega1**3*omega2 - 926.*omega1**2*omega2**2 - &
              17.*omega1*omega2**3 + 412.*omega2**4))*swd**2) + 10077696.*grav*(12.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(11.*kd1**2 + 2.*kd1*kd2 + &
              11.*kd2**2) + (-3.*kd2**9*omega1**3*omega2 + 3.*kd1**9*omega1*omega2**3 + 3.*kd1**8*kd2*omega2**2*(19.*omega1**2 - 27.*omega1*omega2 + &
              12.*omega2**2) - 3.*kd1*kd2**8*omega1**2*(12.*omega1**2 - 27.*omega1*omega2 + 19.*omega2**2) + kd1**7*kd2**2*omega2*(-116.*omega1**3 + &
              59.*omega1**2*omega2 + 156.*omega1*omega2**2 - 124.*omega2**3) + kd1**2*kd2**7*omega1*(124.*omega1**3 - 156.*omega1**2*omega2 - &
              59.*omega1*omega2**2 + 116.*omega2**3) + kd1**5*kd2**4*(-564.*omega1**4 + 553.*omega1**3*omega2 - 751.*omega1**2*omega2**2 + &
              1005.*omega1*omega2**3 - 730.*omega2**4) + kd1**3*kd2**6*(-396.*omega1**4 + 503.*omega1**3*omega2 + 3.*omega1**2*omega2**2 - &
              127.*omega1*omega2**3 - 150.*omega2**4) + kd1**6*kd2**3*(150.*omega1**4 + 127.*omega1**3*omega2 - 3.*omega1**2*omega2**2 - &
              503.*omega1*omega2**3 + 396.*omega2**4) + kd1**4*kd2**5*(730.*omega1**4 - 1005.*omega1**3*omega2 + 751.*omega1**2*omega2**2 - &
              553.*omega1*omega2**3 + 564.*omega2**4))*swd**2) - 216.*grav*kd1**3*(kd1 - kd2)**3*kd2**3*(108.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2 + &
              (kd1**5*kd2*(-92.*omega1**4 + 229.*omega1**3*omega2 - 906.*omega1**2*omega2**2 + 1240.*omega1*omega2**3 - 520.*omega2**4) + &
              kd1*kd2**5*(-520.*omega1**4 + 1240.*omega1**3*omega2 - 906.*omega1**2*omega2**2 + 229.*omega1*omega2**3 - 92.*omega2**4) + &
              6.*kd2**6*(19.*omega1**4 - 45.*omega1**3*omega2 + 31.*omega1**2*omega2**2 - 4.*omega1*omega2**3 + 2.*omega2**4) + &
              6.*kd1**6*(2.*omega1**4 - 4.*omega1**3*omega2 + 31.*omega1**2*omega2**2 - 45.*omega1*omega2**3 + 19.*omega2**4) - &
              2.*kd1**3*kd2**3*(358.*omega1**4 - 685.*omega1**3*omega2 + 878.*omega1**2*omega2**2 - 685.*omega1*omega2**3 + 358.*omega2**4) + &
              kd1**2*kd2**4*(846.*omega1**4 - 1767.*omega1**3*omega2 + 1531.*omega1**2*omega2**2 - 717.*omega1*omega2**3 + 362.*omega2**4) + &
              kd1**4*kd2**2*(362.*omega1**4 - 717.*omega1**3*omega2 + 1531.*omega1**2*omega2**2 - 1767.*omega1*omega2**3 + 846.*omega2**4))*swd**2) - &
              7776.*grav*kd1*(kd1 - kd2)*kd2*(12.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(3.*kd1**4 - 14.*kd1**3*kd2 - 4.*kd1**2*kd2**2 - &
              14.*kd1*kd2**3 + 3.*kd2**4) + (3.*kd1**10*omega2**2*(5.*omega1**2 - 8.*omega1*omega2 + 4.*omega2**2) + 3.*kd2**10*omega1**2*(4.*omega1**2 - &
              8.*omega1*omega2 + 5.*omega2**2) - kd1*kd2**9*omega1*(112.*omega1**3 - 245.*omega1**2*omega2 + 124.*omega1*omega2**2 + 3.*omega2**3) - &
              kd1**9*kd2*omega2*(3.*omega1**3 + 124.*omega1**2*omega2 - 245.*omega1*omega2**2 + 112.*omega2**3) + kd1**3*kd2**7*(-186.*omega1**4 + &
              252.*omega1**3*omega2 - 1108.*omega1**2*omega2**2 + 853.*omega1*omega2**3 - 454.*omega2**4) + kd1**7*kd2**3*(-454.*omega1**4 + &
              853.*omega1**3*omega2 - 1108.*omega1**2*omega2**2 + 252.*omega1*omega2**3 - 186.*omega2**4) + kd1**2*kd2**8*(200.*omega1**4 - &
              537.*omega1**3*omega2 + 345.*omega1**2*omega2**2 - 39.*omega1*omega2**3 + 6.*omega2**4) + kd1**8*kd2**2*(6.*omega1**4 - &
              39.*omega1**3*omega2 + 345.*omega1**2*omega2**2 - 537.*omega1*omega2**3 + 200.*omega2**4) + kd1**6*kd2**4*(1552.*omega1**4 - &
              2325.*omega1**3*omega2 + 3076.*omega1**2*omega2**2 - 327.*omega1*omega2**3 + 750.*omega2**4) - 2.*kd1**5*kd2**5*(884.*omega1**4 - &
              955.*omega1**3*omega2 + 2209.*omega1**2*omega2**2 - 955.*omega1*omega2**3 + 884.*omega2**4) + kd1**4*kd2**6*(750.*omega1**4 - &
              327.*omega1**3*omega2 + 3076.*omega1**2*omega2**2 - 2325.*omega1*omega2**3 + 1552.*omega2**4))*swd**2)))/ &
              (kd1**2*(-36. + kd1**2)**2*kd2**2*(-36. + kd2**2)**2*omega1*(omega1 - omega2)*omega2*swd**3*(-36.*grav*(12. + &
              (kd1 - kd2)**2)*(108. + (kd1 - kd2)**2)*(kd1 - kd2)**2 + (36. + (kd1 - kd2)**2)*(1296. + &
              (504. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function velsb13
!
real function velsp13()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic velocity of 1st layer for 3 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp13')
    !
    velsp13 = (3.*(13060694016.*grav*(kd1 + kd2)**3*omega1*omega2*(kd2**2*omega1**2 + 6.*kd1*kd2*omega1*omega2 + kd1**2*omega2**2)*swd**2 - &
              6.*grav*kd1**5*kd2**5*(kd1 + kd2)**5*(2.*(8.*kd1**2 + 14.*kd1*kd2 + 9.*kd2**2)*omega1**4 + 2.*(18.*kd1**2 + 55.*kd1*kd2 + 11.*kd2**2)*omega1**3*omega2 + &
              (27.*kd1**2 + 170.*kd1*kd2 + 27.*kd2**2)*omega1**2*omega2**2 + 2.*(11.*kd1**2 + 55.*kd1*kd2 + 18.*kd2**2)*omega1*omega2**3 + 2.*(9.*kd1**2 + &
              14.*kd1*kd2 + 8.*kd2**2)*omega2**4)*swd**2 + kd1**5*kd2**5*(kd1 + kd2)**5*omega1*omega2*(omega1 + omega2)**2*(kd1**2*omega2*(5.*omega1 + 2.*omega2) + &
              kd2**2*omega1*(2.*omega1 + 5.*omega2) - 2.*kd1*kd2*(omega1**2 - omega1*omega2 + omega2**2))*swd**3 - &
              2176782336.*swd*(6.*grav**2*kd1*kd2*(kd1 + kd2)*(kd2**2*(kd1 + 2.*kd2)*(3.*kd1 + 2.*kd2)*omega1**2 + kd1*kd2*(kd1**2 + kd2**2)*omega1*omega2 + &
              kd1**2*(2.*kd1 + kd2)*(2.*kd1 + 3.*kd2)*omega2**2) + omega1*omega2*(omega1 + omega2)**2*(kd2**2*(2.*kd1**3 + 9.*kd1**2*kd2 + 3.*kd1*kd2**2 + &
              kd2**3)*omega1**2 + 3.*kd1*kd2*(kd1 + kd2)*(2.*kd1 + kd2)*(kd1 + 2.*kd2)*omega1*omega2 + kd1**2*(kd1**3 + 3.*kd1**2*kd2 + 9.*kd1*kd2**2 + &
              2.*kd2**3)*omega2**2)*swd**2) - 60466176.*swd*(2.*grav**2*kd1*kd2*(kd1 + kd2)**2*(kd2**2*(81.*kd1**3 + 55.*kd1**2*kd2 - 20.*kd1*kd2**2 + &
              18.*kd2**3)*omega1**2 - 5.*kd1*kd2*(kd1 + kd2)*(7.*kd1**2 - 3.*kd1*kd2 + 7.*kd2**2)*omega1*omega2 + kd1**2*(18.*kd1**3 - 20.*kd1**2*kd2 + &
              55.*kd1*kd2**2 + 81.*kd2**3)*omega2**2) - omega1*omega2*(omega1 + omega2)**2*(kd2**2*(4.*kd1**5 - 71.*kd1**4*kd2 - 151.*kd1**3*kd2**2 - &
              55.*kd1**2*kd2**3 + 9.*kd1*kd2**4 + 3.*kd2**5)*omega1**2 - 3.*kd1*kd2*(kd1 + kd2)*(13.*kd1**4 + 33.*kd1**3*kd2 + 45.*kd1**2*kd2**2 + &
              33.*kd1*kd2**3 + 13.*kd2**4)*omega1*omega2 + kd1**2*(3.*kd1**5 + 9.*kd1**4*kd2 - 55.*kd1**3*kd2**2 - 151.*kd1**2*kd2**3 - 71.*kd1*kd2**4 + &
              4.*kd2**5)*omega2**2)*swd**2) + 362797056.*grav*(108.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3 + (2.*kd1*kd2**3*(15.*kd1**3 + 46.*kd1**2*kd2 + &
              41.*kd1*kd2**2 + 12.*kd2**3)*omega1**4 + kd2**2*(-6.*kd1**5 + 85.*kd1**4*kd2 + 236.*kd1**3*kd2**2 + 179.*kd1**2*kd2**3 + 39.*kd1*kd2**4 - &
              3.*kd2**5)*omega1**3*omega2 + kd1*kd2*(kd1 + kd2)*(63.*kd1**4 + 138.*kd1**3*kd2 + 172.*kd1**2*kd2**2 + 138.*kd1*kd2**3 + &
              63.*kd2**4)*omega1**2*omega2**2 + kd1**2*(-3.*kd1**5 + 39.*kd1**4*kd2 + 179.*kd1**3*kd2**2 + 236.*kd1**2*kd2**3 + 85.*kd1*kd2**4 - &
              6.*kd2**5)*omega1*omega2**3 + 2.*kd1**3*kd2*(12.*kd1**3 + 41.*kd1**2*kd2 + 46.*kd1*kd2**2 + 15.*kd2**3)*omega2**4)*swd**2) + &
              10077696.*grav*(12.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(11.*kd1**2 - 2.*kd1*kd2 + 11.*kd2**2) + (2.*kd1*kd2**3*(75.*kd1**5 + &
              282.*kd1**4*kd2 + 365.*kd1**3*kd2**2 + 198.*kd1**2*kd2**3 + 62.*kd1*kd2**4 + 18.*kd2**5)*omega1**4 + kd2**2*(-116.*kd1**7 - 127.*kd1**6*kd2 + &
              553.*kd1**5*kd2**2 + 1005.*kd1**4*kd2**3 + 503.*kd1**3*kd2**4 + 156.*kd1**2*kd2**5 + 81.*kd1*kd2**6 + 3.*kd2**7)*omega1**3*omega2 + &
              kd1*kd2*(kd1 + kd2)*(57.*kd1**6 - 116.*kd1**5*kd2 + 113.*kd1**4*kd2**2 + 638.*kd1**3*kd2**3 + 113.*kd1**2*kd2**4 - 116.*kd1*kd2**5 + &
              57.*kd2**6)*omega1**2*omega2**2 + kd1**2*(3.*kd1**7 + 81.*kd1**6*kd2 + 156.*kd1**5*kd2**2 + 503.*kd1**4*kd2**3 + 1005.*kd1**3*kd2**4 + &
              553.*kd1**2*kd2**5 - 127.*kd1*kd2**6 - 116.*kd2**7)*omega1*omega2**3 + 2.*kd1**3*kd2*(18.*kd1**5 + 62.*kd1**4*kd2 + 198.*kd1**3*kd2**2 + &
              365.*kd1**2*kd2**3 + 282.*kd1*kd2**4 + 75.*kd2**5)*omega2**4)*swd**2) - &
              7776.*grav*kd1*kd2*(kd1 + kd2)*(12.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(3.*kd1**4 + 14.*kd1**3*kd2 - 4.*kd1**2*kd2**2 + 14.*kd1*kd2**3 + &
              3.*kd2**4) + (2.*kd2**2*(kd1 + kd2)*(3.*kd1**7 + 224.*kd1**6*kd2 + 552.*kd1**5*kd2**2 + 332.*kd1**4*kd2**3 + 43.*kd1**3*kd2**4 + 50.*kd1**2*kd2**5 + &
              50.*kd1*kd2**6 + 6.*kd2**7)*omega1**4 + kd2*(-3.*kd1**9 + 39.*kd1**8*kd2 + 853.*kd1**7*kd2**2 + 2325.*kd1**6*kd2**3 + 1910.*kd1**5*kd2**4 + &
              327.*kd1**4*kd2**5 + 252.*kd1**3*kd2**6 + 537.*kd1**2*kd2**7 + 245.*kd1*kd2**8 + 24.*kd2**9)*omega1**3*omega2 + (15.*kd1**10 + 124.*kd1**9*kd2 + &
              345.*kd1**8*kd2**2 + 1108.*kd1**7*kd2**3 + 3076.*kd1**6*kd2**4 + 4418.*kd1**5*kd2**5 + 3076.*kd1**4*kd2**6 + 1108.*kd1**3*kd2**7 + &
              345.*kd1**2*kd2**8 + 124.*kd1*kd2**9 + 15.*kd2**10)*omega1**2*omega2**2 + kd1*(24.*kd1**9 + 245.*kd1**8*kd2 + 537.*kd1**7*kd2**2 + &
              252.*kd1**6*kd2**3 + 327.*kd1**5*kd2**4 + 1910.*kd1**4*kd2**5 + 2325.*kd1**3*kd2**6 + 853.*kd1**2*kd2**7 + 39.*kd1*kd2**8 - 3.*kd2**9)*omega1*omega2**3 + &
              2.*kd1**2*(kd1 + kd2)*(6.*kd1**7 + 50.*kd1**6*kd2 + 50.*kd1**5*kd2**2 + 43.*kd1**4*kd2**3 + 332.*kd1**3*kd2**4 + 552.*kd1**2*kd2**5 + &
              224.*kd1*kd2**6 + 3.*kd2**7)*omega2**4)*swd**2) + 36.*kd1**3*kd2**3*(kd1 + kd2)**3*swd*(2.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(3.*kd2**2*(21.*omega1**2 + &
              25.*omega1*omega2 + 14.*omega2**2) + 2.*kd1*kd2*(51.*omega1**2 + 82.*omega1*omega2 + 51.*omega2**2) + kd1**2*(42.*omega1**2 + 75.*omega1*omega2 + &
              63.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(378.*kd1**3*kd2**3*omega1*omega2 + 2.*kd1**6*omega2*(10.*omega1 + omega2) + &
              2.*kd2**6*omega1*(omega1 + 10.*omega2) + 2.*kd1*kd2**5*(57.*omega1**2 + 56.*omega1*omega2 + 2.*omega2**2) + kd1**2*kd2**4*(131.*omega1**2 + &
              300.*omega1*omega2 + 4.*omega2**2) + 2.*kd1**5*kd2*(2.*omega1**2 + 56.*omega1*omega2 + 57.*omega2**2) + kd1**4*kd2**2*(4.*omega1**2 + &
              300.*omega1*omega2 + 131.*omega2**2))*swd**2) + 1679616.*swd*(2.*grav**2*kd1**2*kd2**2*(kd1 + kd2)*(kd2**6*omega1*(82.*omega1 + 3.*omega2) + &
              kd1**6*omega2*(3.*omega1 + 82.*omega2) + kd1*kd2**5*(114.*omega1**2 - 179.*omega1*omega2 - 213.*omega2**2) - 2.*kd1**4*kd2**2*(227.*omega1**2 + &
              251.*omega1*omega2 + 35.*omega2**2) + kd1**5*kd2*(-213.*omega1**2 - 179.*omega1*omega2 + 114.*omega2**2) - 2.*kd1**3*kd2**3*(167.*omega1**2 + &
              311.*omega1*omega2 + 167.*omega2**2) - 2.*kd1**2*kd2**4*(35.*omega1**2 + 251.*omega1*omega2 + 227.*omega2**2)) + omega1*omega2*(omega1 + &
              omega2)**2*(-3.*kd2**9*omega1**2 - 3.*kd1**9*omega2**2 - 3.*kd1**8*kd2*omega2*(7.*omega1 + 3.*omega2) - 3.*kd1*kd2**8*omega1*(3.*omega1 + 7.*omega2) + &
              4.*kd1**2*kd2**7*(31.*omega1**2 + 5.*omega1*omega2 + 15.*omega2**2) + 4.*kd1**7*kd2**2*(15.*omega1**2 + 5.*omega1*omega2 + 31.*omega2**2) + &
              kd1**5*kd2**4*(317.*omega1**2 - 507.*omega1*omega2 + 239.*omega2**2) + kd1**3*kd2**6*(268.*omega1**2 - 71.*omega1*omega2 + 261.*omega2**2) + &
              kd1**6*kd2**3*(261.*omega1**2 - 71.*omega1*omega2 + 268.*omega2**2) + kd1**4*kd2**5*(239.*omega1**2 - 507.*omega1*omega2 + 317.*omega2**2))*swd**2) + &
              46656.*swd*(2.*grav**2*kd1*kd2*(kd1 + kd2)**2*(6.*kd2**9*omega1**2 + 3.*kd1*kd2**8*omega1*(10.*omega1 - omega2) - 3.*kd1**8*kd2*(omega1 - &
              10.*omega2)*omega2 + 6.*kd1**9*omega2**2 + kd1**6*kd2**3*(203.*omega1**2 - 8.*omega1*omega2 + 37.*omega2**2) - kd1**2*kd2**7*(87.*omega1**2 + &
              82.*omega1*omega2 + 39.*omega2**2) - kd1**7*kd2**2*(39.*omega1**2 + 82.*omega1*omega2 + 87.*omega2**2) + kd1**3*kd2**6*(37.*omega1**2 - &
              8.*omega1*omega2 + 203.*omega2**2) + 2.*kd1**5*kd2**4*(378.*omega1**2 + 152.*omega1*omega2 + 325.*omega2**2) + 2.*kd1**4*kd2**5*(325.*omega1**2 + &
              152.*omega1*omega2 + 378.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(kd2**11*omega1**2 + kd1**11*omega2**2 + &
              3.*kd1**10*kd2*omega2*(5.*omega1 + omega2) + 3.*kd1*kd2**10*omega1*(omega1 + 5.*omega2) + kd1**6*kd2**5*(231.*omega1**2 - 764.*omega1*omega2 - &
              163.*omega2**2) + kd1**7*kd2**4*(411.*omega1**2 + 290.*omega1*omega2 - 145.*omega2**2) + kd1**8*kd2**3*(147.*omega1**2 + 596.*omega1*omega2 - &
              73.*omega2**2) + kd1**9*kd2**2*(4.*omega1**2 + 214.*omega1*omega2 - 47.*omega2**2) + kd1**2*kd2**9*(-47.*omega1**2 + 214.*omega1*omega2 + &
              4.*omega2**2) + kd1**3*kd2**8*(-73.*omega1**2 + 596.*omega1*omega2 + 147.*omega2**2) + kd1**5*kd2**6*(-163.*omega1**2 - 764.*omega1*omega2 + &
              231.*omega2**2) + kd1**4*kd2**7*(-145.*omega1**2 + 290.*omega1*omega2 + 411.*omega2**2))*swd**2) + &
              1296.*kd1*kd2*(kd1 + kd2)*swd*(2.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(3.*kd2**6*(25.*omega1**2 + 27.*omega1*omega2 + 18.*omega2**2) + &
              kd1**6*(54.*omega1**2 + 81.*omega1*omega2 + 75.*omega2**2) + kd1*kd2**5*(420.*omega1**2 + 433.*omega1*omega2 + 246.*omega2**2) + &
              2.*kd1**3*kd2**3*(299.*omega1**2 + 272.*omega1*omega2 + 299.*omega2**2) + kd1**2*kd2**4*(720.*omega1**2 + 621.*omega1*omega2 + 418.*omega2**2) + &
              kd1**5*kd2*(246.*omega1**2 + 433.*omega1*omega2 + 420.*omega2**2) + kd1**4*kd2**2*(418.*omega1**2 + 621.*omega1*omega2 + 720.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(-3.*kd1**10*omega1*omega2 - 3.*kd2**10*omega1*omega2 + kd1*kd2**9*(13.*omega1**2 - 34.*omega1*omega2 + &
              2.*omega2**2) + kd1**2*kd2**8*(31.*omega1**2 - 188.*omega1*omega2 + 6.*omega2**2) + kd1**3*kd2**7*(112.*omega1**2 - 949.*omega1*omega2 + &
              12.*omega2**2) + kd1**9*kd2*(2.*omega1**2 - 34.*omega1*omega2 + 13.*omega2**2) + kd1**8*kd2**2*(6.*omega1**2 - 188.*omega1*omega2 + &
              31.*omega2**2) + kd1**7*kd2**3*(12.*omega1**2 - 949.*omega1*omega2 + 112.*omega2**2) + kd1**4*kd2**6*(527.*omega1**2 - 2527.*omega1*omega2 + &
              282.*omega2**2) + 2.*kd1**5*kd2**5*(355.*omega1**2 - 1732.*omega1*omega2 + 355.*omega2**2) + kd1**6*kd2**4*(282.*omega1**2 - 2527.*omega1*omega2 + &
              527.*omega2**2))*swd**2) + 279936.*grav*(12.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(7.*kd1**4 - 26.*kd1**3*kd2 - 47.*kd1**2*kd2**2 - 26.*kd1*kd2**3 + &
              7.*kd2**4) - (kd2**11*omega1**3*omega2 + kd1**11*omega1*omega2**3 + 3.*kd1**10*kd2*omega1*omega2**2*(5.*omega1 + omega2) + &
              3.*kd1*kd2**10*omega1**2*omega2*(omega1 + 5.*omega2) + kd1**2*kd2**9*omega1*(194.*omega1**3 + 357.*omega1**2*omega2 + 405.*omega1*omega2**2 + &
              10.*omega2**3) + kd1**9*kd2**2*omega2*(10.*omega1**3 + 405.*omega1**2*omega2 + 357.*omega1*omega2**2 + 194.*omega2**3) + kd1**4*kd2**7*(364.*omega1**4 + &
              515.*omega1**3*omega2 + 926.*omega1**2*omega2**2 - 17.*omega1*omega2**3 - 412.*omega2**4) + kd1**3*kd2**8*(574.*omega1**4 + 1032.*omega1**3*omega2 + &
              1190.*omega1**2*omega2**2 + 29.*omega1*omega2**3 - 126.*omega2**4) + kd1**7*kd2**4*(-412.*omega1**4 - 17.*omega1**3*omega2 + 926.*omega1**2*omega2**2 + &
              515.*omega1*omega2**3 + 364.*omega2**4) - kd1**6*kd2**5*(636.*omega1**4 + 545.*omega1**3*omega2 + 458.*omega1**2*omega2**2 + 679.*omega1*omega2**3 + &
              370.*omega2**4) + kd1**8*kd2**3*(-126.*omega1**4 + 29.*omega1**3*omega2 + 1190.*omega1**2*omega2**2 + 1032.*omega1*omega2**3 + 574.*omega2**4) - &
              kd1**5*kd2**6*(370.*omega1**4 + 679.*omega1**3*omega2 + 458.*omega1**2*omega2**2 + 545.*omega1*omega2**3 + 636.*omega2**4))*swd**2) - &
              216.*grav*kd1**3*kd2**3*(kd1 + kd2)**3*(108.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4 + (6.*kd2**6*(19.*omega1**4 + 45.*omega1**3*omega2 + &
              31.*omega1**2*omega2**2 + 4.*omega1*omega2**3 + 2.*omega2**4) + 6.*kd1**6*(2.*omega1**4 + 4.*omega1**3*omega2 + 31.*omega1**2*omega2**2 + &
              45.*omega1*omega2**3 + 19.*omega2**4) + kd1*kd2**5*(520.*omega1**4 + 1240.*omega1**3*omega2 + 906.*omega1**2*omega2**2 + 229.*omega1*omega2**3 + &
              92.*omega2**4) + 2.*kd1**3*kd2**3*(358.*omega1**4 + 685.*omega1**3*omega2 + 878.*omega1**2*omega2**2 + 685.*omega1*omega2**3 + 358.*omega2**4) + &
              kd1**2*kd2**4*(846.*omega1**4 + 1767.*omega1**3*omega2 + 1531.*omega1**2*omega2**2 + 717.*omega1*omega2**3 + 362.*omega2**4) + kd1**5*kd2*(92.*omega1**4 + &
              229.*omega1**3*omega2 + 906.*omega1**2*omega2**2 + 1240.*omega1*omega2**3 + 520.*omega2**4) + kd1**4*kd2**2*(362.*omega1**4 + 717.*omega1**3*omega2 + &
              1531.*omega1**2*omega2**2 + 1767.*omega1*omega2**3 + 846.*omega2**4))*swd**2)))/ &
              (kd1**2*(-36. + kd1**2)**2*kd2**2*(-36. + kd2**2)**2*omega1*omega2*(omega1 + omega2)*swd**3*(-36.*grav*(kd1 + kd2)**2*(12. + (kd1 + kd2)**2)*(108. + &
              (kd1 + kd2)**2) + (36. + (kd1 + kd2)**2)*(1296. + (kd1 + kd2)**2*(504. + (kd1 + kd2)**2))*(omega1 + omega2)**2*swd))
    !
end function velsp13
!
real function velsb23()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic velocity of 2nd layer for 3 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb23')
    !
    velsb23 = (-3.*(52242776064.*grav*(kd1 - kd2)*(3.*kd1*kd2**3*omega1**4 - kd2**2*(kd1**2 + 4.*kd1*kd2 + kd2**2)*omega1**3*omega2 + &
              6.*kd1**2*kd2**2*omega1**2*omega2**2 - kd1**2*(kd1**2 + 4.*kd1*kd2 + kd2**2)*omega1*omega2**3 + 3.*kd1**3*kd2*omega2**4)*swd**2 - &
              12.*grav*kd1**5*(kd1 - kd2)**5*kd2**5*(3.*kd2**2*(12.*omega1**4 - 6.*omega1**3*omega2 + 6.*omega1**2*omega2**2 - 13.*omega1*omega2**3 + &
              7.*omega2**4) - 4.*kd1*kd2*(12.*omega1**4 - 37.*omega1**3*omega2 + 59.*omega1**2*omega2**2 - 37.*omega1*omega2**3 + 12.*omega2**4) + &
              3.*kd1**2*(7.*omega1**4 - 13.*omega1**3*omega2 + 6.*omega1**2*omega2**2 - 6.*omega1*omega2**3 + 12.*omega2**4))*swd**2 - &
              725594112.*grav*(kd1 - kd2)*(2.*kd2**6*omega1**3*omega2 + 2.*kd1**6*omega1*omega2**3 + kd1**5*kd2*omega2**2*(79.*omega1**2 - &
              101.*omega1*omega2 + 42.*omega2**2) + kd1*kd2**5*omega1**2*(42.*omega1**2 - 101.*omega1*omega2 + 79.*omega2**2) + &
              2.*kd1**4*kd2**2*omega2*(30.*omega1**3 - 195.*omega1**2*omega2 + 214.*omega1*omega2**2 - 76.*omega2**3) + &
              2.*kd1**2*kd2**4*omega1*(-76.*omega1**3 + 214.*omega1**2*omega2 - 195.*omega1*omega2**2 + 30.*omega2**3) + kd1**3*kd2**3*(21.*omega1**4 - &
              255.*omega1**3*omega2 + 532.*omega1**2*omega2**2 - 255.*omega1*omega2**3 + 21.*omega2**4))*swd**2 - &
              156728328192.*kd1*(kd1 - kd2)*kd2*omega1**2*(omega1 - omega2)**2*omega2**2*swd**3 + &
              kd1**5*(kd1 - kd2)**5*kd2**5*omega1*(omega1 - omega2)**2*omega2*(kd1**2*omega2*(7.*omega1 + 5.*omega2) + kd2**2*omega1*(5.*omega1 + 7.*omega2) - &
              8.*kd1*kd2*(omega1**2 + omega1*omega2 + omega2**2))*swd**3 + &
              2176782336.*swd*(24.*grav**2*kd1*(kd1 - kd2)*kd2*(kd2**2*(-3.*kd1 + kd2)*(kd1 + kd2)*omega1**2 + &
              kd1*kd2*(kd1**2 + 6.*kd1*kd2 + kd2**2)*omega1*omega2 + kd1**2*(kd1 - 3.*kd2)*(kd1 + kd2)*omega2**2) + &
              omega1*(omega1 - omega2)**2*omega2*(kd2**2*(24.*kd1**3 - 59.*kd1**2*kd2 + 35.*kd1*kd2**2 - 12.*kd2**3)*omega1**2 + &
              kd1*(kd1 - kd2)*kd2*(19.*kd1**2 - 104.*kd1*kd2 + 19.*kd2**2)*omega1*omega2 + kd1**2*(12.*kd1**3 - 35.*kd1**2*kd2 + 59.*kd1*kd2**2 - &
              24.*kd2**3)*omega2**2)*swd**2) + 120932352.*swd*(2.*grav**2*kd1*(kd1 - kd2)*kd2*(36.*kd2**6*omega1**2 - 2.*kd1**3*kd2**3*(7.*omega1 - &
              13.*omega2)*(13.*omega1 - 7.*omega2) + 2.*kd1**5*kd2*(omega1 - 70.*omega2)*omega2 + 36.*kd1**6*omega2**2 + &
              2.*kd1*kd2**5*omega1*(-70.*omega1 + omega2) + kd1**2*kd2**4*(193.*omega1**2 - 148.*omega1*omega2 + 21.*omega2**2) + kd1**4*kd2**2*(21.*omega1**2 - &
              148.*omega1*omega2 + 193.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(10.*kd2**7*omega1**2 - 10.*kd1**7*omega2**2 + &
              kd1**6*kd2*omega2*(83.*omega1 + 21.*omega2) - kd1*kd2**6*omega1*(21.*omega1 + 83.*omega2) + kd1**2*kd2**5*(-99.*omega1**2 + &
              371.*omega1*omega2 - 8.*omega2**2) + 6.*kd1**3*kd2**4*(50.*omega1**2 - 95.*omega1*omega2 + 34.*omega2**2) - 6.*kd1**4*kd2**3*(34.*omega1**2 - &
              95.*omega1*omega2 + 50.*omega2**2) + kd1**5*kd2**2*(8.*omega1**2 - 371.*omega1*omega2 + 99.*omega2**2))*swd**2) - &
              72.*kd1**3*(kd1 - kd2)**3*kd2**3*swd*(6.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd1**2*(-25.*omega1**2 + 35.*omega1*omega2 - 46.*omega2**2) + &
              kd2**2*(-46.*omega1**2 + 35.*omega1*omega2 - 25.*omega2**2) + 4.*kd1*kd2*(17.*omega1**2 - 22.*omega1*omega2 + 17.*omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(2.*kd2**6*omega1*(omega1 - 19.*omega2) + 2.*kd1**6*omega2*(-19.*omega1 + omega2) - &
              4.*kd1*kd2**5*(40.*omega1**2 - 45.*omega1*omega2 + 2.*omega2**2) + 2.*kd1**3*kd2**3*(3.*omega1**2 + 283.*omega1*omega2 + 3.*omega2**2) + &
              kd1**2*kd2**4*(105.*omega1**2 - 391.*omega1*omega2 + 21.*omega2**2) - 4.*kd1**5*kd2*(2.*omega1**2 - 45.*omega1*omega2 + 40.*omega2**2) + &
              kd1**4*kd2**2*(21.*omega1**2 - 391.*omega1*omega2 + 105.*omega2**2))*swd**2) + 93312.*swd*(2.*grav**2*kd1*(kd1 - kd2)**3*kd2*(36.*kd2**8*omega1**2 + &
              36.*kd1**8*omega2**2 - 18.*kd1*kd2**7*omega1*(4.*omega1 + omega2) - 18.*kd1**7*kd2*omega2*(omega1 + 4.*omega2) + 2.*kd1**3*kd2**5*(159.*omega1**2 - &
              86.*omega1*omega2 - 292.*omega2**2) - 2.*kd1**5*kd2**3*(292.*omega1**2 + 86.*omega1*omega2 - 159.*omega2**2) + 3.*kd1**6*kd2**2*(29.*omega1**2 - &
              2.*omega1*omega2 - 123.*omega2**2) - 3.*kd1**2*kd2**6*(123.*omega1**2 + 2.*omega1*omega2 - 29.*omega2**2) - 3.*kd1**4*kd2**4*(19.*omega1**2 - &
              564.*omega1*omega2 + 19.*omega2**2)) - omega1*(omega1 - omega2)**2*omega2*(2.*kd2**11*omega1**2 - 2.*kd1**11*omega2**2 + &
              kd1**10*kd2*omega2*(-62.*omega1 + 13.*omega2) + kd1*kd2**10*omega1*(-13.*omega1 + 62.*omega2) + kd1**5*kd2**6*(-284.*omega1**2 + 125.*omega1*omega2 - &
              980.*omega2**2) + kd1**9*kd2**2*(24.*omega1**2 + 282.*omega1*omega2 - 157.*omega2**2) + kd1**2*kd2**9*(157.*omega1**2 - 282.*omega1*omega2 - &
              24.*omega2**2) - 2.*kd1**3*kd2**8*(158.*omega1**2 - 39.*omega1*omega2 + 14.*omega2**2) - 2.*kd1**7*kd2**4*(314.*omega1**2 + 165.*omega1*omega2 + &
              34.*omega2**2) + 2.*kd1**8*kd2**3*(14.*omega1**2 - 39.*omega1*omega2 + 158.*omega2**2) + kd1**6*kd2**5*(980.*omega1**2 - 125.*omega1*omega2 + &
              284.*omega2**2) + 2.*kd1**4*kd2**7*(34.*omega1**2 + 165.*omega1*omega2 + 314.*omega2**2))*swd**2) + &
              3359232.*swd*(2.*grav**2*kd1*kd2*(-kd1 + kd2)*(-66.*kd2**8*omega1**2 - 66.*kd1**8*omega2**2 + 6.*kd1*kd2**7*omega1*(56.*omega1 + 9.*omega2) + &
              6.*kd1**7*kd2*omega2*(9.*omega1 + 56.*omega2) + kd1**4*kd2**4*(-1027.*omega1**2 + 2138.*omega1*omega2 - 1027.*omega2**2) + &
              3.*kd1**6*kd2**2*(5.*omega1**2 + 13.*omega1*omega2 - 176.*omega2**2) + 3.*kd1**2*kd2**6*(-176.*omega1**2 + 13.*omega1*omega2 + 5.*omega2**2) + &
              2.*kd1**3*kd2**5*(427.*omega1**2 - 629.*omega1*omega2 + 256.*omega2**2) + 2.*kd1**5*kd2**3*(256.*omega1**2 - 629.*omega1*omega2 + 427.*omega2**2)) - &
              omega1*(omega1 - omega2)**2*omega2*(2.*kd2**9*omega1**2 - 2.*kd1**9*omega2**2 - 2.*kd1**8*kd2*omega2*(64.*omega1 + 5.*omega2) + &
              2.*kd1*kd2**8*omega1*(5.*omega1 + 64.*omega2) - 2.*kd1**2*kd2**7*(143.*omega1**2 + 249.*omega1*omega2 + 60.*omega2**2) + kd1**3*kd2**6*(631.*omega1**2 + &
              102.*omega1*omega2 + 111.*omega2**2) + kd1**4*kd2**5*(-257.*omega1**2 + 1000.*omega1*omega2 + 139.*omega2**2) + 2.*kd1**7*kd2**2*(60.*omega1**2 + &
              249.*omega1*omega2 + 143.*omega2**2) + kd1**5*kd2**4*(-139.*omega1**2 - 1000.*omega1*omega2 + 257.*omega2**2) - kd1**6*kd2**3*(111.*omega1**2 + &
              102.*omega1*omega2 + 631.*omega2**2))*swd**2) - 1296.*kd1*(kd1 - kd2)*kd2*swd*(12.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd1**6*(-39.*omega1**2 + &
              45.*omega1*omega2 - 78.*omega2**2) + 56.*kd1**3*kd2**3*(7.*omega1**2 - 6.*omega1*omega2 + 7.*omega2**2) - 3.*kd2**6*(26.*omega1**2 - 15.*omega1*omega2 + &
              13.*omega2**2) + 2.*kd1*kd2**5*(188.*omega1**2 - 128.*omega1*omega2 + 69.*omega2**2) - 2.*kd1**2*kd2**4*(289.*omega1**2 - 188.*omega1*omega2 + &
              106.*omega2**2) + 2.*kd1**5*kd2*(69.*omega1**2 - 128.*omega1*omega2 + 188.*omega2**2) - 2.*kd1**4*kd2**2*(106.*omega1**2 - 188.*omega1*omega2 + &
              289.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(kd2**10*omega1*(omega1 - 13.*omega2) + kd1**10*omega2*(-13.*omega1 + omega2) + &
              4.*kd1*kd2**9*omega2*(23.*omega1 + 2.*omega2) + 4.*kd1**9*kd2*omega1*(2.*omega1 + 23.*omega2) + kd1**2*kd2**8*(2.*omega1**2 - 657.*omega1*omega2 - &
              37.*omega2**2) + 4.*kd1**3*kd2**7*(91.*omega1**2 + 748.*omega1*omega2 - 7.*omega2**2) + kd1**8*kd2**2*(-37.*omega1**2 - 657.*omega1*omega2 + &
              2.*omega2**2) + 4.*kd1**7*kd2**3*(-7.*omega1**2 + 748.*omega1*omega2 + 91.*omega2**2) - kd1**4*kd2**6*(1675.*omega1**2 + 7623.*omega1*omega2 + &
              270.*omega2**2) + 4.*kd1**5*kd2**5*(403.*omega1**2 + 2616.*omega1*omega2 + 403.*omega2**2) - kd1**6*kd2**4*(270.*omega1**2 + 7623.*omega1*omega2 + &
              1675.*omega2**2))*swd**2) + 20155392.*grav*(kd1 - kd2)*(24.*grav**2*kd1**3*(kd1 - kd2)**2*kd2**3*(kd1**2 + 10.*kd1*kd2 + kd2**2) + &
              (10.*kd2**8*omega1**3*omega2 + 10.*kd1**8*omega1*omega2**3 + kd1*kd2**7*omega1**2*(-126.*omega1**2 + 209.*omega1*omega2 - 223.*omega2**2) + &
              kd1**7*kd2*omega2**2*(-223.*omega1**2 + 209.*omega1*omega2 - 126.*omega2**2) + 4.*kd1**2*kd2**6*omega1*(108.*omega1**3 - 181.*omega1**2*omega2 + &
              157.*omega1*omega2**2 + 36.*omega2**3) + 4.*kd1**6*kd2**2*omega2*(36.*omega1**3 + 157.*omega1**2*omega2 - 181.*omega1*omega2**2 + 108.*omega2**3) + &
              kd1**5*kd2**3*(54.*omega1**4 + 198.*omega1**3*omega2 - 1214.*omega1**2*omega2**2 + 1309.*omega1*omega2**3 - 471.*omega2**4) + &
              kd1**3*kd2**5*(-471.*omega1**4 + 1309.*omega1**3*omega2 - 1214.*omega1**2*omega2**2 + 198.*omega1*omega2**3 + 54.*omega2**4) + &
              4.*kd1**4*kd2**4*(86.*omega1**4 - 401.*omega1**3*omega2 + 517.*omega1**2*omega2**2 - 401.*omega1*omega2**3 + 86.*omega2**4))*swd**2) + &
              559872.*grav*(kd1 - kd2)*(240.*grav**2*kd1**4*(kd1 - kd2)**2*kd2**4*(3.*kd1**2 - 4.*kd1*kd2 + 3.*kd2**2) - (6.*kd2**10*omega1**3*omega2 + &
              6.*kd1**10*omega1*omega2**3 + kd1**9*kd2*omega2**2*(229.*omega1**2 - 223.*omega1*omega2 + 102.*omega2**2) + kd1*kd2**9*omega1**2*(102.*omega1**2 - &
              223.*omega1*omega2 + 229.*omega2**2) - 2.*kd1**2*kd2**8*omega1*(112.*omega1**3 - 110.*omega1**2*omega2 + 375.*omega1*omega2**2 + 78.*omega2**3) - &
              2.*kd1**8*kd2**2*omega2*(78.*omega1**3 + 375.*omega1**2*omega2 - 110.*omega1*omega2**2 + 112.*omega2**3) + kd1**7*kd2**3*(48.*omega1**4 + &
              32.*omega1**3*omega2 - 411.*omega1**2*omega2**2 + 2284.*omega1*omega2**3 - 781.*omega2**4) + 2.*kd1**4*kd2**6*(558.*omega1**4 - 1494.*omega1**3*omega2 + &
              1077.*omega1**2*omega2**2 + 755.*omega1*omega2**3 - 290.*omega2**4) + kd1**3*kd2**7*(-781.*omega1**4 + 2284.*omega1**3*omega2 - 411.*omega1**2*omega2**2 + &
              32.*omega1*omega2**3 + 48.*omega2**4) + kd1**5*kd2**5*(260.*omega1**4 - 515.*omega1**3*omega2 - 2666.*omega1**2*omega2**2 - 515.*omega1*omega2**3 + &
              260.*omega2**4) + 2.*kd1**6*kd2**4*(-290.*omega1**4 + 755.*omega1**3*omega2 + 1077.*omega1**2*omega2**2 - 1494.*omega1*omega2**3 + 558.*omega2**4))*swd**2) - &
              15552.*grav*kd1*(kd1 - kd2)*kd2*(24.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(3.*kd1**2 - 11.*kd1*kd2 + 3.*kd2**2)*(3.*kd1**2 + 5.*kd1*kd2 + 3.*kd2**2) + &
              (3.*kd1**10*omega2**2*(15.*omega1**2 - 17.*omega1*omega2 + 8.*omega2**2) + 3.*kd2**10*omega1**2*(8.*omega1**2 - 17.*omega1*omega2 + 15.*omega2**2) - &
              6.*kd1*kd2**9*omega1*(28.*omega1**3 - 75.*omega1**2*omega2 + 34.*omega1*omega2**2 + 3.*omega2**3) - 6.*kd1**9*kd2*omega2*(3.*omega1**3 + &
              34.*omega1**2*omega2 - 75.*omega1*omega2**2 + 28.*omega2**3) + kd1**2*kd2**8*(357.*omega1**4 - 1061.*omega1**3*omega2 + 552.*omega1**2*omega2**2 - &
              298.*omega1*omega2**3 + 156.*omega2**4) + kd1**8*kd2**2*(156.*omega1**4 - 298.*omega1**3*omega2 + 552.*omega1**2*omega2**2 - 1061.*omega1*omega2**3 + &
              357.*omega2**4) - 2.*kd1**7*kd2**3*(376.*omega1**4 - 951.*omega1**3*omega2 + 1680.*omega1**2*omega2**2 - 756.*omega1*omega2**3 + 358.*omega2**4) - &
              4.*kd1**5*kd2**5*(375.*omega1**4 - 319.*omega1**3*omega2 + 2662.*omega1**2*omega2**2 - 319.*omega1*omega2**3 + 375.*omega2**4) - 2.*kd1**3*kd2**7*(358.*omega1**4 - &
              756.*omega1**3*omega2 + 1680.*omega1**2*omega2**2 - 951.*omega1*omega2**3 + 376.*omega2**4) + 2.*kd1**6*kd2**4*(664.*omega1**4 - 1285.*omega1**3*omega2 + &
              4147.*omega1**2*omega2**2 - 574.*omega1*omega2**3 + 637.*omega2**4) + 2.*kd1**4*kd2**6*(637.*omega1**4 - 574.*omega1**3*omega2 + 4147.*omega1**2*omega2**2 - &
              1285.*omega1*omega2**3 + 664.*omega2**4))*swd**2) - 432.*grav*kd1**3*(kd1 - kd2)**3*kd2**3*(360.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2 + &
              (3.*kd2**6*(60.*omega1**4 - 137.*omega1**3*omega2 + 123.*omega1**2*omega2**2 - 19.*omega1*omega2**3 + 9.*omega2**4) + 3.*kd1**6*(9.*omega1**4 - &
              19.*omega1**3*omega2 + 123.*omega1**2*omega2**2 - 137.*omega1*omega2**3 + 60.*omega2**4) - 2.*kd1*kd2**5*(408.*omega1**4 - 1031.*omega1**3*omega2 + &
              839.*omega1**2*omega2**2 - 196.*omega1*omega2**3 + 67.*omega2**4) - 4.*kd1**3*kd2**3*(266.*omega1**4 - 476.*omega1**3*omega2 + 681.*omega1**2*omega2**2 - &
              476.*omega1*omega2**3 + 266.*omega2**4) - 2.*kd1**5*kd2*(67.*omega1**4 - 196.*omega1**3*omega2 + 839.*omega1**2*omega2**2 - 1031.*omega1*omega2**3 + &
              408.*omega2**4) + kd1**2*kd2**4*(1404.*omega1**4 - 2970.*omega1**3*omega2 + 2604.*omega1**2*omega2**2 - 865.*omega1*omega2**3 + 415.*omega2**4) + &
              kd1**4*kd2**2*(415.*omega1**4 - 865.*omega1**3*omega2 + 2604.*omega1**2*omega2**2 - 2970.*omega1*omega2**3 + 1404.*omega2**4))*swd**2)))/ &
              (4.*kd1**2*(-36. + kd1**2)**2*kd2**2*(-36. + kd2**2)**2*omega1*(omega1 - omega2)*omega2*swd**3*(-36.*grav*(12. + (kd1 - kd2)**2)*(108. + &
              (kd1 - kd2)**2)*(kd1 - kd2)**2 + (36. + (kd1 - kd2)**2)*(1296. + (504. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function velsb23
!
real function velsp23()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic velocity of 2nd layer for 3 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp23')
    !
    velsp23 = (-3.*(-52242776064.*grav*(kd1 + kd2)*(kd2**4*omega1**3*omega2 + kd1**4*omega1*omega2**3 - kd1**3*kd2*omega2**3*(4.*omega1 + 3.*omega2) - &
              kd1*kd2**3*omega1**3*(3.*omega1 + 4.*omega2) + kd1**2*kd2**2*omega1*omega2*(omega1**2 + 6.*omega1*omega2 + omega2**2))*swd**2 - &
              12.*grav*kd1**5*kd2**5*(kd1 + kd2)**5*(3.*kd2**2*(12.*omega1**4 + 6.*omega1**3*omega2 + 6.*omega1**2*omega2**2 + 13.*omega1*omega2**3 + &
              7.*omega2**4) + 3.*kd1**2*(7.*omega1**4 + 13.*omega1**3*omega2 + 6.*omega1**2*omega2**2 + 6.*omega1*omega2**3 + 12.*omega2**4) + &
              4.*kd1*kd2*(12.*omega1**4 + 37.*omega1**3*omega2 + 59.*omega1**2*omega2**2 + 37.*omega1*omega2**3 + 12.*omega2**4))*swd**2 - &
              725594112.*grav*(kd1 + kd2)*(2.*kd2**6*omega1**3*omega2 + 2.*kd1**6*omega1*omega2**3 + kd1**5*kd2*omega2**2*(79.*omega1**2 + &
              101.*omega1*omega2 + 42.*omega2**2) + kd1*kd2**5*omega1**2*(42.*omega1**2 + 101.*omega1*omega2 + 79.*omega2**2) + &
              2.*kd1**2*kd2**4*omega1*(76.*omega1**3 + 214.*omega1**2*omega2 + 195.*omega1*omega2**2 + 30.*omega2**3) + &
              2.*kd1**4*kd2**2*omega2*(30.*omega1**3 + 195.*omega1**2*omega2 + 214.*omega1*omega2**2 + 76.*omega2**3) + kd1**3*kd2**3*(21.*omega1**4 + &
              255.*omega1**3*omega2 + 532.*omega1**2*omega2**2 + 255.*omega1*omega2**3 + 21.*omega2**4))*swd**2 - &
              156728328192.*kd1*kd2*(kd1 + kd2)*omega1**2*omega2**2*(omega1 + omega2)**2*swd**3 - &
              kd1**5*kd2**5*(kd1 + kd2)**5*omega1*omega2*(omega1 + omega2)**2*(kd2**2*omega1*(5.*omega1 - 7.*omega2) + kd1**2*omega2*(-7.*omega1 + 5.*omega2) + &
              8.*kd1*kd2*(omega1**2 - omega1*omega2 + omega2**2))*swd**3 + 2176782336.*swd*(24.*grav**2*kd1*kd2*(kd1 + kd2)*(kd2**4*omega1**2 + &
              kd1**4*omega2**2 - 3.*kd1**2*kd2**2*(omega1 + omega2)**2 + kd1*kd2**3*omega1*(2.*omega1 + omega2) + kd1**3*kd2*omega2*(omega1 + 2.*omega2)) + &
              omega1*omega2*(omega1 + omega2)**2*(kd2**2*(24.*kd1**3 + 59.*kd1**2*kd2 + 35.*kd1*kd2**2 + 12.*kd2**3)*omega1**2 + kd1*kd2*(kd1 + kd2)*(19.*kd1**2 + &
              104.*kd1*kd2 + 19.*kd2**2)*omega1*omega2 + kd1**2*(12.*kd1**3 + 35.*kd1**2*kd2 + 59.*kd1*kd2**2 + 24.*kd2**3)*omega2**2)*swd**2) + &
              120932352.*swd*(2.*grav**2*kd1*kd2*(kd1 + kd2)*(36.*kd2**6*omega1**2 + 36.*kd1**6*omega2**2 + 2.*kd1*kd2**5*omega1*(70.*omega1 + omega2) + &
              2.*kd1**3*kd2**3*(13.*omega1 + 7.*omega2)*(7.*omega1 + 13.*omega2) + 2.*kd1**5*kd2*omega2*(omega1 + 70.*omega2) + kd1**2*kd2**4*(193.*omega1**2 + &
              148.*omega1*omega2 + 21.*omega2**2) + kd1**4*kd2**2*(21.*omega1**2 + 148.*omega1*omega2 + 193.*omega2**2)) + &
              omega1*omega2*(omega1 + omega2)**2*(-10.*kd2**7*omega1**2 + kd1**6*kd2*(83.*omega1 - 21.*omega2)*omega2 - 10.*kd1**7*omega2**2 + &
              kd1*kd2**6*omega1*(-21.*omega1 + 83.*omega2) + kd1**2*kd2**5*(99.*omega1**2 + 371.*omega1*omega2 + 8.*omega2**2) + &
              6.*kd1**3*kd2**4*(50.*omega1**2 + 95.*omega1*omega2 + 34.*omega2**2) + 6.*kd1**4*kd2**3*(34.*omega1**2 + 95.*omega1*omega2 + 50.*omega2**2) + &
              kd1**5*kd2**2*(8.*omega1**2 + 371.*omega1*omega2 + 99.*omega2**2))*swd**2) + &
              72.*kd1**3*kd2**3*(kd1 + kd2)**3*swd*(6.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*((25.*kd1**2 + 68.*kd1*kd2 + 46.*kd2**2)*omega1**2 + &
              (35.*kd1**2 + 88.*kd1*kd2 + 35.*kd2**2)*omega1*omega2 + (46.*kd1**2 + 68.*kd1*kd2 + 25.*kd2**2)*omega2**2) + &
              omega1*omega2*(omega1 + omega2)**2*(2.*kd1**6*omega2*(19.*omega1 + omega2) + 2.*kd2**6*omega1*(omega1 + 19.*omega2) + &
              4.*kd1*kd2**5*(40.*omega1**2 + 45.*omega1*omega2 + 2.*omega2**2) - 2.*kd1**3*kd2**3*(3.*omega1**2 - 283.*omega1*omega2 + 3.*omega2**2) + &
              kd1**2*kd2**4*(105.*omega1**2 + 391.*omega1*omega2 + 21.*omega2**2) + 4.*kd1**5*kd2*(2.*omega1**2 + 45.*omega1*omega2 + 40.*omega2**2) + &
              kd1**4*kd2**2*(21.*omega1**2 + 391.*omega1*omega2 + 105.*omega2**2))*swd**2) + 3359232.*swd*(2.*grav**2*kd1*kd2*(kd1 + kd2)*(66.*kd2**8*omega1**2 + &
              6.*kd1*kd2**7*omega1*(56.*omega1 - 9.*omega2) + 66.*kd1**8*omega2**2 + 6.*kd1**7*kd2*omega2*(-9.*omega1 + 56.*omega2) + &
              3.*kd1**2*kd2**6*(176.*omega1**2 + 13.*omega1*omega2 - 5.*omega2**2) + 3.*kd1**6*kd2**2*(-5.*omega1**2 + 13.*omega1*omega2 + 176.*omega2**2) + &
              2.*kd1**3*kd2**5*(427.*omega1**2 + 629.*omega1*omega2 + 256.*omega2**2) + 2.*kd1**5*kd2**3*(256.*omega1**2 + 629.*omega1*omega2 + 427.*omega2**2) + &
              kd1**4*kd2**4*(1027.*omega1**2 + 2138.*omega1*omega2 + 1027.*omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(-2.*kd2**9*omega1**2 + &
              2.*kd1*kd2**8*omega1*(5.*omega1 - 64.*omega2) - 2.*kd1**9*omega2**2 + 2.*kd1**8*kd2*omega2*(-64.*omega1 + 5.*omega2) + kd1**4*kd2**5*(257.*omega1**2 + &
              1000.*omega1*omega2 - 139.*omega2**2) + 2.*kd1**2*kd2**7*(143.*omega1**2 - 249.*omega1*omega2 + 60.*omega2**2) + kd1**3*kd2**6*(631.*omega1**2 - &
              102.*omega1*omega2 + 111.*omega2**2) + 2.*kd1**7*kd2**2*(60.*omega1**2 - 249.*omega1*omega2 + 143.*omega2**2) + kd1**5*kd2**4*(-139.*omega1**2 + &
              1000.*omega1*omega2 + 257.*omega2**2) + kd1**6*kd2**3*(111.*omega1**2 - 102.*omega1*omega2 + 631.*omega2**2))*swd**2) + &
              93312.*swd*(2.*grav**2*kd1*kd2*(kd1 + kd2)**3*(36.*kd2**8*omega1**2 + 18.*kd1*kd2**7*omega1*(4.*omega1 - omega2) - 18.*kd1**7*kd2*(omega1 - &
              4.*omega2)*omega2 + 36.*kd1**8*omega2**2 - 2.*kd1**3*kd2**5*(159.*omega1**2 + 86.*omega1*omega2 - 292.*omega2**2) + 2.*kd1**5*kd2**3*(292.*omega1**2 - &
              86.*omega1*omega2 - 159.*omega2**2) + 3.*kd1**6*kd2**2*(29.*omega1**2 + 2.*omega1*omega2 - 123.*omega2**2) - 3.*kd1**4*kd2**4*(19.*omega1**2 + &
              564.*omega1*omega2 + 19.*omega2**2) + 3.*kd1**2*kd2**6*(-123.*omega1**2 + 2.*omega1*omega2 + 29.*omega2**2)) + &
              omega1*omega2*(omega1 + omega2)**2*(2.*kd2**11*omega1**2 + 2.*kd1**11*omega2**2 + kd1**10*kd2*omega2*(62.*omega1 + 13.*omega2) + &
              kd1*kd2**10*omega1*(13.*omega1 + 62.*omega2) + kd1**2*kd2**9*(157.*omega1**2 + 282.*omega1*omega2 - 24.*omega2**2) + 2.*kd1**3*kd2**8*(158.*omega1**2 + &
              39.*omega1*omega2 + 14.*omega2**2) + 2.*kd1**7*kd2**4*(314.*omega1**2 - 165.*omega1*omega2 + 34.*omega2**2) + kd1**9*kd2**2*(-24.*omega1**2 + &
              282.*omega1*omega2 + 157.*omega2**2) + 2.*kd1**8*kd2**3*(14.*omega1**2 + 39.*omega1*omega2 + 158.*omega2**2) + kd1**6*kd2**5*(980.*omega1**2 + &
              125.*omega1*omega2 + 284.*omega2**2) + 2.*kd1**4*kd2**7*(34.*omega1**2 - 165.*omega1*omega2 + 314.*omega2**2) + kd1**5*kd2**6*(284.*omega1**2 + &
              125.*omega1*omega2 + 980.*omega2**2))*swd**2) + 1296.*kd1*kd2*(kd1 + kd2)*swd*(12.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(56.*kd1**3*kd2**3*(7.*omega1**2 + &
              6.*omega1*omega2 + 7.*omega2**2) + 3.*kd2**6*(26.*omega1**2 + 15.*omega1*omega2 + 13.*omega2**2) + 2.*kd1*kd2**5*(188.*omega1**2 + 128.*omega1*omega2 + &
              69.*omega2**2) + kd1**6*(39.*omega1**2 + 45.*omega1*omega2 + 78.*omega2**2) + 2.*kd1**2*kd2**4*(289.*omega1**2 + 188.*omega1*omega2 + 106.*omega2**2) + &
              2.*kd1**5*kd2*(69.*omega1**2 + 128.*omega1*omega2 + 188.*omega2**2) + 2.*kd1**4*kd2**2*(106.*omega1**2 + 188.*omega1*omega2 + 289.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(4.*kd1**9*kd2*omega1*(2.*omega1 - 23.*omega2) - kd1**10*omega2*(13.*omega1 + omega2) + &
              4.*kd1*kd2**9*omega2*(-23.*omega1 + 2.*omega2) - kd2**10*omega1*(omega1 + 13.*omega2) - 4.*kd1**7*kd2**3*(7.*omega1**2 + &
              748.*omega1*omega2 - 91.*omega2**2) + 4.*kd1**3*kd2**7*(91.*omega1**2 - 748.*omega1*omega2 - 7.*omega2**2) + kd1**8*kd2**2*(37.*omega1**2 - &
              657.*omega1*omega2 - 2.*omega2**2) + kd1**2*kd2**8*(-2.*omega1**2 - 657.*omega1*omega2 + 37.*omega2**2) + kd1**4*kd2**6*(1675.*omega1**2 - &
              7623.*omega1*omega2 + 270.*omega2**2) + 4.*kd1**5*kd2**5*(403.*omega1**2 - 2616.*omega1*omega2 + 403.*omega2**2) + kd1**6*kd2**4*(270.*omega1**2 - &
              7623.*omega1*omega2 + 1675.*omega2**2))*swd**2) + 20155392.*grav*(kd1 + kd2)*(24.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**2*(kd1**2 - 10.*kd1*kd2 + kd2**2) + &
              (10.*kd2**8*omega1**3*omega2 + 10.*kd1**8*omega1*omega2**3 - kd1**7*kd2*omega2**2*(223.*omega1**2 + 209.*omega1*omega2 + 126.*omega2**2) - &
              kd1*kd2**7*omega1**2*(126.*omega1**2 + 209.*omega1*omega2 + 223.*omega2**2) - 4.*kd1**2*kd2**6*omega1*(108.*omega1**3 + 181.*omega1**2*omega2 + &
              157.*omega1*omega2**2 - 36.*omega2**3) - 4.*kd1**6*kd2**2*omega2*(-36.*omega1**3 + 157.*omega1**2*omega2 + 181.*omega1*omega2**2 + 108.*omega2**3) + &
              kd1**5*kd2**3*(54.*omega1**4 - 198.*omega1**3*omega2 - 1214.*omega1**2*omega2**2 - 1309.*omega1*omega2**3 - 471.*omega2**4) - &
              kd1**3*kd2**5*(471.*omega1**4 + 1309.*omega1**3*omega2 + 1214.*omega1**2*omega2**2 + 198.*omega1*omega2**3 - 54.*omega2**4) - &
              4.*kd1**4*kd2**4*(86.*omega1**4 + 401.*omega1**3*omega2 + 517.*omega1**2*omega2**2 + 401.*omega1*omega2**3 + 86.*omega2**4))*swd**2) - &
              559872.*grav*(kd1 + kd2)*(240.*grav**2*kd1**4*kd2**4*(kd1 + kd2)**2*(3.*kd1**2 + 4.*kd1*kd2 + 3.*kd2**2) + (6.*kd2**10*omega1**3*omega2 + &
              6.*kd1**10*omega1*omega2**3 + kd1**9*kd2*omega2**2*(229.*omega1**2 + 223.*omega1*omega2 + 102.*omega2**2) + kd1*kd2**9*omega1**2*(102.*omega1**2 + &
              223.*omega1*omega2 + 229.*omega2**2) + 2.*kd1**2*kd2**8*omega1*(112.*omega1**3 + 110.*omega1**2*omega2 + 375.*omega1*omega2**2 - 78.*omega2**3) + &
              2.*kd1**8*kd2**2*omega2*(-78.*omega1**3 + 375.*omega1**2*omega2 + 110.*omega1*omega2**2 + 112.*omega2**3) + kd1**7*kd2**3*(48.*omega1**4 - &
              32.*omega1**3*omega2 - 411.*omega1**2*omega2**2 - 2284.*omega1*omega2**3 - 781.*omega2**4) + 2.*kd1**6*kd2**4*(290.*omega1**4 + 755.*omega1**3*omega2 - &
              1077.*omega1**2*omega2**2 - 1494.*omega1*omega2**3 - 558.*omega2**4) - 2.*kd1**4*kd2**6*(558.*omega1**4 + 1494.*omega1**3*omega2 + &
              1077.*omega1**2*omega2**2 - 755.*omega1*omega2**3 - 290.*omega2**4) - kd1**3*kd2**7*(781.*omega1**4 + 2284.*omega1**3*omega2 + 411.*omega1**2*omega2**2 + &
              32.*omega1*omega2**3 - 48.*omega2**4) + kd1**5*kd2**5*(260.*omega1**4 + 515.*omega1**3*omega2 - 2666.*omega1**2*omega2**2 + 515.*omega1*omega2**3 + &
              260.*omega2**4))*swd**2) - 15552.*grav*kd1*kd2*(kd1 + kd2)*(24.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(3.*kd1**2 - 5.*kd1*kd2 + 3.*kd2**2)*(3.*kd1**2 + &
              11.*kd1*kd2 + 3.*kd2**2) + (3.*kd1**10*omega2**2*(15.*omega1**2 + 17.*omega1*omega2 + 8.*omega2**2) + 3.*kd2**10*omega1**2*(8.*omega1**2 + &
              17.*omega1*omega2 + 15.*omega2**2) + 6.*kd1*kd2**9*omega1*(28.*omega1**3 + 75.*omega1**2*omega2 + 34.*omega1*omega2**2 - 3.*omega2**3) + &
              6.*kd1**9*kd2*omega2*(-3.*omega1**3 + 34.*omega1**2*omega2 + 75.*omega1*omega2**2 + 28.*omega2**3) + kd1**2*kd2**8*(357.*omega1**4 + &
              1061.*omega1**3*omega2 + 552.*omega1**2*omega2**2 + 298.*omega1*omega2**3 + 156.*omega2**4) + kd1**8*kd2**2*(156.*omega1**4 + 298.*omega1**3*omega2 + &
              552.*omega1**2*omega2**2 + 1061.*omega1*omega2**3 + 357.*omega2**4) + 2.*kd1**7*kd2**3*(376.*omega1**4 + 951.*omega1**3*omega2 + &
              1680.*omega1**2*omega2**2 + 756.*omega1*omega2**3 + 358.*omega2**4) + 4.*kd1**5*kd2**5*(375.*omega1**4 + 319.*omega1**3*omega2 + &
              2662.*omega1**2*omega2**2 + 319.*omega1*omega2**3 + 375.*omega2**4) + 2.*kd1**3*kd2**7*(358.*omega1**4 + 756.*omega1**3*omega2 + &
              1680.*omega1**2*omega2**2 + 951.*omega1*omega2**3 + 376.*omega2**4) + 2.*kd1**6*kd2**4*(664.*omega1**4 + 1285.*omega1**3*omega2 + &
              4147.*omega1**2*omega2**2 + 574.*omega1*omega2**3 + 637.*omega2**4) + 2.*kd1**4*kd2**6*(637.*omega1**4 + 574.*omega1**3*omega2 + &
              4147.*omega1**2*omega2**2 + 1285.*omega1*omega2**3 + 664.*omega2**4))*swd**2) - &
              432.*grav*kd1**3*kd2**3*(kd1 + kd2)**3*(360.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4 + (3.*kd2**6*(60.*omega1**4 + 137.*omega1**3*omega2 + &
              123.*omega1**2*omega2**2 + 19.*omega1*omega2**3 + 9.*omega2**4) + 3.*kd1**6*(9.*omega1**4 + 19.*omega1**3*omega2 + 123.*omega1**2*omega2**2 + &
              137.*omega1*omega2**3 + 60.*omega2**4) + 2.*kd1*kd2**5*(408.*omega1**4 + 1031.*omega1**3*omega2 + 839.*omega1**2*omega2**2 + 196.*omega1*omega2**3 + &
              67.*omega2**4) + 4.*kd1**3*kd2**3*(266.*omega1**4 + 476.*omega1**3*omega2 + 681.*omega1**2*omega2**2 + 476.*omega1*omega2**3 + 266.*omega2**4) + &
              2.*kd1**5*kd2*(67.*omega1**4 + 196.*omega1**3*omega2 + 839.*omega1**2*omega2**2 + 1031.*omega1*omega2**3 + 408.*omega2**4) + kd1**2*kd2**4*(1404.*omega1**4 + &
              2970.*omega1**3*omega2 + 2604.*omega1**2*omega2**2 + 865.*omega1*omega2**3 + 415.*omega2**4) + kd1**4*kd2**2*(415.*omega1**4 + 865.*omega1**3*omega2 + &
              2604.*omega1**2*omega2**2 + 2970.*omega1*omega2**3 + 1404.*omega2**4))*swd**2)))/ &
              (4.*kd1**2*(-36. + kd1**2)**2*kd2**2*(-36. + kd2**2)**2*omega1*omega2*(omega1 + omega2)*swd**3*(-36.*grav*(kd1 + kd2)**2*(12. + (kd1 + kd2)**2)*(108. + &
              (kd1 + kd2)**2) + (36. + (kd1 + kd2)**2)*(1296. + (kd1 + kd2)**2*(504. + (kd1 + kd2)**2))*(omega1 + omega2)**2*swd))
    !
end function velsp23
!
real function velsb33()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic velocity of 3rd layer for 3 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb33')
    !
    velsb33 = (3.*(-52242776064.*grav*(kd1 - kd2)*(2.*kd2**4*omega1**3*omega2 + 2.*kd1**4*omega1*omega2**3 + 2.*kd1**2*kd2**2*omega1*omega2*(omega1**2 - &
              9.*omega1*omega2 + omega2**2) + kd1**3*kd2*omega2**2*(15.*omega1**2 - 19.*omega1*omega2 + 9.*omega2**2) + kd1*kd2**3*omega1**2*(9.*omega1**2 - &
              19.*omega1*omega2 + 15.*omega2**2))*swd**2 + 12.*grav*kd1**5*(kd1 - kd2)**5*kd2**5*(8.*kd1*kd2*(2.*omega1**4 - 5.*omega1**3*omega2 + &
              9.*omega1**2*omega2**2 - 5.*omega1*omega2**3 + 2.*omega2**4) - kd2**2*(12.*omega1**4 + 2.*omega1**3*omega2 - 9.*omega1*omega2**3 + 7.*omega2**4) - &
              kd1**2*(7.*omega1**4 - 9.*omega1**3*omega2 + 2.*omega1*omega2**3 + 12.*omega2**4))*swd**2 - &
              156728328192.*omega1*(omega1 - omega2)**2*omega2*(kd2**3*omega1**2 + 4.*kd1*kd2*(-kd1 + kd2)*omega1*omega2 - kd1**3*omega2**2)*swd**3 + &
              kd1**5*(kd1 - kd2)**5*kd2**5*omega1*(omega1 - omega2)**2*omega2*(3.*kd2**2*omega1*(omega1 + omega2) + 3.*kd1**2*omega2*(omega1 + omega2) - &
              4.*kd1*kd2*(omega1**2 + omega1*omega2 + omega2**2))*swd**3 + 2176782336.*swd*(24.*grav**2*kd1*(kd1 - kd2)*kd2*(5.*kd2**4*omega1**2 + &
              2.*kd1**3*kd2*(omega1 - 5.*omega2)*omega2 + 5.*kd1**4*omega2**2 + 2.*kd1*kd2**3*omega1*(-5.*omega1 + omega2) + 6.*kd1**2*kd2**2*(2.*omega1**2 - &
              3.*omega1*omega2 + 2.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(kd2**2*(-36.*kd1**3 + 11.*kd1**2*kd2 + 9.*kd1*kd2**2 - 8.*kd2**3)*omega1**2 + &
              kd1*(kd1 - kd2)*kd2*(83.*kd1**2 - 16.*kd1*kd2 + 83.*kd2**2)*omega1*omega2 + kd1**2*(8.*kd1**3 - 9.*kd1**2*kd2 - 11.*kd1*kd2**2 + &
              36.*kd2**3)*omega2**2)*swd**2) - 725594112.*grav*(216.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3 + (kd1*kd2**3*(123.*kd1**3 - 215.*kd1**2*kd2 + &
              186.*kd1*kd2**2 - 102.*kd2**3)*omega1**4 + kd2**2*(-48.*kd1**5 - 221.*kd1**4*kd2 + 483.*kd1**3*kd2**2 - 399.*kd1**2*kd2**3 + 197.*kd1*kd2**4 + &
              4.*kd2**5)*omega1**3*omega2 + kd1*(kd1 - kd2)*kd2*(209.*kd1**4 - 250.*kd1**3*kd2 + 456.*kd1**2*kd2**2 - 250.*kd1*kd2**3 + &
              209.*kd2**4)*omega1**2*omega2**2 + kd1**2*(-4.*kd1**5 - 197.*kd1**4*kd2 + 399.*kd1**3*kd2**2 - 483.*kd1**2*kd2**3 + 221.*kd1*kd2**4 + &
              48.*kd2**5)*omega1*omega2**3 + kd1**3*kd2*(102.*kd1**3 - 186.*kd1**2*kd2 + 215.*kd1*kd2**2 - 123.*kd2**3)*omega2**4)*swd**2) - &
              559872.*grav*(24.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(7.*kd1**4 - 4.*kd1**3*kd2 - 7.*kd1**2*kd2**2 - 4.*kd1*kd2**3 + 7.*kd2**4) + &
              (kd1*kd2**3*(252.*kd1**7 - 1028.*kd1**6*kd2 + 848.*kd1**5*kd2**2 + 1220.*kd1**4*kd2**3 - 2163.*kd1**3*kd2**4 + 811.*kd1**2*kd2**5 + &
              142.*kd1*kd2**6 - 90.*kd2**7)*omega1**4 - kd2**2*(144.*kd1**9 + 28.*kd1**8*kd2 - 1922.*kd1**7*kd2**2 + 2671.*kd1**6*kd2**3 + 1311.*kd1**5*kd2**4 - &
              4188.*kd1**4*kd2**5 + 1958.*kd1**3*kd2**6 + 169.*kd1**2*kd2**7 - 191.*kd1*kd2**8 + 4.*kd2**9)*omega1**3*omega2 + kd1*(kd1 - kd2)*kd2*(179.*kd1**8 - &
              258.*kd1**7*kd2 - 701.*kd1**6*kd2**2 + 1230.*kd1**5*kd2**3 - 946.*kd1**4*kd2**4 + 1230.*kd1**3*kd2**5 - 701.*kd1**2*kd2**6 - 258.*kd1*kd2**7 + &
              179.*kd2**8)*omega1**2*omega2**2 + kd1**2*(4.*kd1**9 - 191.*kd1**8*kd2 + 169.*kd1**7*kd2**2 + 1958.*kd1**6*kd2**3 - 4188.*kd1**5*kd2**4 + &
              1311.*kd1**4*kd2**5 + 2671.*kd1**3*kd2**6 - 1922.*kd1**2*kd2**7 + 28.*kd1*kd2**8 + 144.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(90.*kd1**7 - &
              142.*kd1**6*kd2 - 811.*kd1**5*kd2**2 + 2163.*kd1**4*kd2**3 - 1220.*kd1**3*kd2**4 - 848.*kd1**2*kd2**5 + 1028.*kd1*kd2**6 - &
              252.*kd2**7)*omega2**4)*swd**2) - 72.*kd1**3*(kd1 - kd2)**3*kd2**3*swd*(-2.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(3.*(3.*kd1 - &
              4.*kd2)**2*omega1**2 + (-27.*kd1**2 + 80.*kd1*kd2 - 27.*kd2**2)*omega1*omega2 + 3.*(4.*kd1 - 3.*kd2)**2*omega2**2) - &
              omega1*(omega1 - omega2)**2*omega2*(6.*kd1**6*omega2*(2.*omega1 + omega2) + 6.*kd2**6*omega1*(omega1 + 2.*omega2) + &
              kd1**4*kd2**2*(-11.*omega1**2 + 104.*omega1*omega2 - 14.*omega2**2) + kd1**2*kd2**4*(-14.*omega1**2 + 104.*omega1*omega2 - 11.*omega2**2) + &
              kd1*kd2**5*(41.*omega1**2 - 52.*omega1*omega2 + 4.*omega2**2) - 2.*kd1**3*kd2**3*(7.*omega1**2 + 76.*omega1*omega2 + 7.*omega2**2) + &
              kd1**5*kd2*(4.*omega1**2 - 52.*omega1*omega2 + 41.*omega2**2))*swd**2) + 120932352.*swd*(2.*grav**2*kd1*(kd1 - kd2)*kd2*(72.*kd2**6*omega1**2 + &
              72.*kd1**6*omega2**2 - 68.*kd1*kd2**5*omega1*(2.*omega1 + omega2) - 68.*kd1**5*kd2*omega2*(omega1 + 2.*omega2) - 2.*kd1**3*kd2**3*(191.*omega1**2 - &
              106.*omega1*omega2 + 191.*omega2**2) + kd1**2*kd2**4*(299.*omega1**2 - 110.*omega1*omega2 + 219.*omega2**2) + kd1**4*kd2**2*(219.*omega1**2 - &
              110.*omega1*omega2 + 299.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(10.*kd2**7*omega1**2 - 10.*kd1**7*omega2**2 + &
              kd1**6*kd2*omega2*(157.*omega1 + 13.*omega2) - kd1*kd2**6*omega1*(13.*omega1 + 157.*omega2) + 2.*kd1**4*kd2**3*(14.*omega1**2 + 192.*omega1*omega2 - &
              77.*omega2**2) + 2.*kd1**3*kd2**4*(77.*omega1**2 - 192.*omega1*omega2 - 14.*omega2**2) + kd1**2*kd2**5*(-105.*omega1**2 + &
              347.*omega1*omega2 + 56.*omega2**2) + kd1**5*kd2**2*(-56.*omega1**2 - 347.*omega1*omega2 + 105.*omega2**2))*swd**2) + &
              3359232.*swd*(2.*grav**2*kd1*(kd1 - kd2)*kd2*(66.*kd2**8*omega1**2 + 66.*kd1**8*omega2**2 - 4.*kd1*kd2**7*omega1*(43.*omega1 + 15.*omega2) - &
              4.*kd1**7*kd2*omega2*(15.*omega1 + 43.*omega2) - 2.*kd1**5*kd2**3*(386.*omega1**2 - 483.*omega1*omega2 + 173.*omega2**2) + &
              kd1**6*kd2**2*(303.*omega1**2 - 181.*omega1*omega2 + 192.*omega2**2) + kd1**2*kd2**6*(192.*omega1**2 - 181.*omega1*omega2 + 303.*omega2**2) - &
              2.*kd1**3*kd2**5*(173.*omega1**2 - 483.*omega1*omega2 + 386.*omega2**2) + kd1**4*kd2**4*(783.*omega1**2 - 1558.*omega1*omega2 + &
              783.*omega2**2)) - omega1*(omega1 - omega2)**2*omega2*(4.*kd2**9*omega1**2 - 4.*kd1**9*omega2**2 + 4.*kd1**8*kd2*omega2*(-37.*omega1 + 3.*omega2) + &
              4.*kd1*kd2**8*omega1*(-3.*omega1 + 37.*omega2) + kd1**4*kd2**5*(-231.*omega1**2 + 896.*omega1*omega2 - 337.*omega2**2) - &
              4.*kd1**2*kd2**7*(48.*omega1**2 + 101.*omega1*omega2 + 39.*omega2**2) + 4.*kd1**7*kd2**2*(39.*omega1**2 + 101.*omega1*omega2 + 48.*omega2**2) + &
              kd1**5*kd2**4*(337.*omega1**2 - 896.*omega1*omega2 + 231.*omega2**2) - kd1**6*kd2**3*(407.*omega1**2 + 98.*omega1*omega2 + 341.*omega2**2) + &
              kd1**3*kd2**6*(341.*omega1**2 + 98.*omega1*omega2 + 407.*omega2**2))*swd**2) - 93312.*swd*(2.*grav**2*kd1*(kd1 - kd2)**2*kd2*(24.*kd2**9*omega1**2 - &
              24.*kd1**9*omega2**2 - 12.*kd1*kd2**8*omega1*(4.*omega1 + omega2) + 12.*kd1**8*kd2*omega2*(omega1 + 4.*omega2) + kd1**6*kd2**3*(277.*omega1**2 + &
              126.*omega1*omega2 - 749.*omega2**2) + kd1**4*kd2**5*(-1255.*omega1**2 + 1632.*omega1*omega2 - 565.*omega2**2) + kd1**3*kd2**6*(749.*omega1**2 - &
              126.*omega1*omega2 - 277.*omega2**2) + kd1**2*kd2**7*(-159.*omega1**2 - 80.*omega1*omega2 + 129.*omega2**2) + kd1**7*kd2**2*(-129.*omega1**2 + &
              80.*omega1*omega2 + 159.*omega2**2) + kd1**5*kd2**4*(565.*omega1**2 - 1632.*omega1*omega2 + 1255.*omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(kd2**11*omega1**2 - kd1**11*omega2**2 - kd1**10*kd2*omega2*(31.*omega1 + 3.*omega2) + &
              kd1*kd2**10*omega1*(3.*omega1 + 31.*omega2) + kd1**9*kd2**2*(24.*omega1**2 + 61.*omega1*omega2 - 91.*omega2**2) + kd1**2*kd2**9*(91.*omega1**2 - &
              61.*omega1*omega2 - 24.*omega2**2) + kd1**3*kd2**8*(-286.*omega1**2 - 246.*omega1*omega2 + 61.*omega2**2) + kd1**8*kd2**3*(-61.*omega1**2 + &
              246.*omega1*omega2 + 286.*omega2**2) + kd1**4*kd2**7*(569.*omega1**2 + 398.*omega1*omega2 + 314.*omega2**2) - 2.*kd1**5*kd2**6*(479.*omega1**2 - &
              42.*omega1*omega2 + 470.*omega2**2) + 2.*kd1**6*kd2**5*(470.*omega1**2 - 42.*omega1*omega2 + 479.*omega2**2) - kd1**7*kd2**4*(314.*omega1**2 + &
              398.*omega1*omega2 + 569.*omega2**2))*swd**2) + 1296.*kd1*(kd1 - kd2)*kd2*swd*(4.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(15.*kd2**6*(8.*omega1**2 - &
              3.*omega1*omega2 + 3.*omega2**2) + 15.*kd1**6*(3.*omega1**2 - 3.*omega1*omega2 + 8.*omega2**2) + 2.*kd1**2*kd2**4*(261.*omega1**2 - &
              171.*omega1*omega2 + 14.*omega2**2) - 4.*kd1**3*kd2**3*(31.*omega1**2 - 52.*omega1*omega2 + 31.*omega2**2) - 2.*kd1*kd2**5*(240.*omega1**2 - &
              143.*omega1*omega2 + 57.*omega2**2) - 2.*kd1**5*kd2*(57.*omega1**2 - 143.*omega1*omega2 + 240.*omega2**2) + 2.*kd1**4*kd2**2*(14.*omega1**2 - &
              171.*omega1*omega2 + 261.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(kd1**10*omega2*(5.*omega1 + omega2) + &
              kd2**10*omega1*(omega1 + 5.*omega2) + 4.*kd1**7*kd2**3*(15.*omega1**2 - 157.*omega1*omega2 - 128.*omega2**2) + kd1**4*kd2**6*(853.*omega1**2 + &
              1603.*omega1*omega2 - 114.*omega2**2) - 4.*kd1**3*kd2**7*(128.*omega1**2 + 157.*omega1*omega2 - 15.*omega2**2) + 4.*kd1*kd2**9*(omega1**2 - &
              6.*omega1*omega2 - omega2**2) - 4.*kd1**9*kd2*(omega1**2 + 6.*omega1*omega2 - omega2**2) + kd1**2*kd2**8*(154.*omega1**2 + 133.*omega1*omega2 + &
              19.*omega2**2) - 4.*kd1**5*kd2**5*(121.*omega1**2 + 533.*omega1*omega2 + 121.*omega2**2) + kd1**8*kd2**2*(19.*omega1**2 + 133.*omega1*omega2 + &
              154.*omega2**2) + kd1**6*kd2**4*(-114.*omega1**2 + 1603.*omega1*omega2 + 853.*omega2**2))*swd**2) - &
              15552.*grav*kd1*(kd1 - kd2)*kd2*(24.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(6.*kd1**4 - 4.*kd1**3*kd2 - 33.*kd1**2*kd2**2 - 4.*kd1*kd2**3 + 6.*kd2**4) + &
              (3.*kd1**10*omega2**2*(7.*omega1**2 - 7.*omega1*omega2 + 4.*omega2**2) + 3.*kd2**10*omega1**2*(4.*omega1**2 - 7.*omega1*omega2 + 7.*omega2**2) - &
              4.*kd1*kd2**9*omega1*(16.*omega1**3 - 44.*omega1**2*omega2 + 13.*omega1*omega2**2 + 3.*omega2**3) - 4.*kd1**9*kd2*omega2*(3.*omega1**3 + &
              13.*omega1**2*omega2 - 44.*omega1*omega2**2 + 16.*omega2**3) - 4.*kd1**3*kd2**7*(44.*omega1**4 - 96.*omega1**3*omega2 + 176.*omega1**2*omega2**2 - &
              125.*omega1*omega2**3 + 33.*omega2**4) - 4.*kd1**7*kd2**3*(33.*omega1**4 - 125.*omega1**3*omega2 + 176.*omega1**2*omega2**2 - 96.*omega1*omega2**3 + &
              44.*omega2**4) + kd1**2*kd2**8*(167.*omega1**4 - 373.*omega1**3*omega2 - 66.*omega1**2*omega2**2 - 32.*omega1*omega2**3 + 60.*omega2**4) - &
              2.*kd1**4*kd2**6*(149.*omega1**4 - 270.*omega1**3*omega2 - 920.*omega1**2*omega2**2 + 34.*omega1*omega2**3 + 144.*omega2**4) - &
              2.*kd1**6*kd2**4*(144.*omega1**4 + 34.*omega1**3*omega2 - 920.*omega1**2*omega2**2 - 270.*omega1*omega2**3 + 149.*omega2**4) + &
              kd1**8*kd2**2*(60.*omega1**4 - 32.*omega1**3*omega2 - 66.*omega1**2*omega2**2 - 373.*omega1*omega2**3 + 167.*omega2**4) + &
              4.*kd1**5*kd2**5*(179.*omega1**4 - 270.*omega1**3*omega2 - 525.*omega1**2*omega2**2 - 270.*omega1*omega2**3 + 179.*omega2**4))*swd**2) - &
              432.*grav*kd1**3*(kd1 - kd2)**3*kd2**3*(144.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2 + (3.*kd2**6*(24.*omega1**4 - 39.*omega1**3*omega2 + &
              41.*omega1**2*omega2**2 - 5.*omega1*omega2**3 + 3.*omega2**4) - 2.*kd1*kd2**5*(152.*omega1**4 - 341.*omega1**3*omega2 + &
              271.*omega1**2*omega2**2 - 61.*omega1*omega2**3 + 17.*omega2**4) + 3.*kd1**6*(3.*omega1**4 - 5.*omega1**3*omega2 + 41.*omega1**2*omega2**2 - &
              39.*omega1*omega2**3 + 24.*omega2**4) - 4.*kd1**3*kd2**3*(70.*omega1**4 - 93.*omega1**3*omega2 + 83.*omega1**2*omega2**2 - &
              93.*omega1*omega2**3 + 70.*omega2**4) + kd1**2*kd2**4*(456.*omega1**4 - 866.*omega1**3*omega2 + 568.*omega1**2*omega2**2 - &
              161.*omega1*omega2**3 + 81.*omega2**4) - 2.*kd1**5*kd2*(17.*omega1**4 - 61.*omega1**3*omega2 + 271.*omega1**2*omega2**2 - &
              341.*omega1*omega2**3 + 152.*omega2**4) + kd1**4*kd2**2*(81.*omega1**4 - 161.*omega1**3*omega2 + 568.*omega1**2*omega2**2 - &
              866.*omega1*omega2**3 + 456.*omega2**4))*swd**2) - 20155392.*grav*(48.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(5.*kd1**2 - &
              4.*kd1*kd2 + 5.*kd2**2) + (4.*kd2**9*omega1**3*omega2 - 4.*kd1**9*omega1*omega2**3 + kd1*kd2**8*omega1**2*(-162.*omega1**2 + &
              317.*omega1*omega2 - 329.*omega2**2) + kd1**8*kd2*omega2**2*(329.*omega1**2 - 317.*omega1*omega2 + 162.*omega2**2) + &
              kd1**2*kd2**7*omega1*(458.*omega1**3 - 697.*omega1**2*omega2 + 597.*omega1*omega2**2 + 312.*omega2**3) - kd1**7*kd2**2*omega2*(312.*omega1**3 + &
              597.*omega1**2*omega2 - 697.*omega1*omega2**2 + 458.*omega2**3) + kd1**5*kd2**4*(-694.*omega1**4 + 1078.*omega1**3*omega2 - &
              1790.*omega1**2*omega2**2 + 1537.*omega1*omega2**3 - 777.*omega2**4) - kd1**3*kd2**6*(669.*omega1**4 - 1045.*omega1**3*omega2 + &
              938.*omega1**2*omega2**2 + 70.*omega1*omega2**3 + 306.*omega2**4) + kd1**6*kd2**3*(306.*omega1**4 + 70.*omega1**3*omega2 + &
              938.*omega1**2*omega2**2 - 1045.*omega1*omega2**3 + 669.*omega2**4) + kd1**4*kd2**5*(777.*omega1**4 - 1537.*omega1**3*omega2 + &
              1790.*omega1**2*omega2**2 - 1078.*omega1*omega2**3 + 694.*omega2**4))*swd**2)))/ &
              (4.*kd1**2*(-36. + kd1**2)**2*kd2**2*(-36. + kd2**2)**2*omega1*(omega1 - omega2)*omega2*swd**3*(-36.*grav*(12. + (kd1 - kd2)**2)*(108. + &
              (kd1 - kd2)**2)*(kd1 - kd2)**2 + (36. + (kd1 - kd2)**2)*(1296. + (504. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function velsb33
!
real function velsp33()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic velocity of 3rd layer for 3 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp33')
    !
    velsp33 = (-3.*(52242776064.*grav*(kd1 + kd2)*(2.*kd2**4*omega1**3*omega2 + 2.*kd1**4*omega1*omega2**3 + 2.*kd1**2*kd2**2*omega1*omega2*(omega1**2 + &
              9.*omega1*omega2 + omega2**2) + kd1**3*kd2*omega2**2*(15.*omega1**2 + 19.*omega1*omega2 + 9.*omega2**2) + kd1*kd2**3*omega1**2*(9.*omega1**2 + &
              19.*omega1*omega2 + 15.*omega2**2))*swd**2 + 12.*grav*kd1**5*kd2**5*(kd1 + kd2)**5*(8.*kd1*kd2*(2.*omega1**4 + 5.*omega1**3*omega2 + &
              9.*omega1**2*omega2**2 + 5.*omega1*omega2**3 + 2.*omega2**4) + kd2**2*(12.*omega1**4 - 2.*omega1**3*omega2 + 9.*omega1*omega2**3 + &
              7.*omega2**4) + kd1**2*(7.*omega1**4 + 9.*omega1**3*omega2 - 2.*omega1*omega2**3 + 12.*omega2**4))*swd**2 - &
              156728328192.*omega1*omega2*(omega1 + omega2)**2*(kd2**3*omega1**2 + 4.*kd1*kd2*(kd1 + kd2)*omega1*omega2 + kd1**3*omega2**2)*swd**3 + &
              kd1**5*kd2**5*(kd1 + kd2)**5*omega1*omega2*(omega1 + omega2)**2*(3.*kd2**2*omega1*(omega1 - omega2) + 3.*kd1**2*omega2*(-omega1 + omega2) + &
              4.*kd1*kd2*(omega1**2 - omega1*omega2 + omega2**2))*swd**3 - 2176782336.*swd*(24.*grav**2*kd1*kd2*(kd1 + kd2)*(5.*kd2**4*omega1**2 + &
              5.*kd1**4*omega2**2 + 2.*kd1*kd2**3*omega1*(5.*omega1 + omega2) + 2.*kd1**3*kd2*omega2*(omega1 + 5.*omega2) + 6.*kd1**2*kd2**2*(2.*omega1**2 + &
              3.*omega1*omega2 + 2.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(kd2**2*(-36.*kd1**3 - 11.*kd1**2*kd2 + 9.*kd1*kd2**2 + &
              8.*kd2**3)*omega1**2 + kd1*kd2*(kd1 + kd2)*(83.*kd1**2 + 16.*kd1*kd2 + 83.*kd2**2)*omega1*omega2 + kd1**2*(8.*kd1**3 + 9.*kd1**2*kd2 - &
              11.*kd1*kd2**2 - 36.*kd2**3)*omega2**2)*swd**2) - 120932352.*swd*(2.*grav**2*kd1*kd2*(kd1 + kd2)*(72.*kd2**6*omega1**2 + &
              68.*kd1*kd2**5*omega1*(2.*omega1 - omega2) + 72.*kd1**6*omega2**2 + 68.*kd1**5*kd2*omega2*(-omega1 + 2.*omega2) + 2.*kd1**3*kd2**3*(191.*omega1**2 + &
              106.*omega1*omega2 + 191.*omega2**2) + kd1**2*kd2**4*(299.*omega1**2 + 110.*omega1*omega2 + 219.*omega2**2) + kd1**4*kd2**2*(219.*omega1**2 + &
              110.*omega1*omega2 + 299.*omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(kd2**2*(56.*kd1**5 + 28.*kd1**4*kd2 - 154.*kd1**3*kd2**2 - &
              105.*kd1**2*kd2**3 + 13.*kd1*kd2**4 + 10.*kd2**5)*omega1**2 - kd1*kd2*(kd1 + kd2)*(157.*kd1**4 + 190.*kd1**3*kd2 + 194.*kd1**2*kd2**2 + &
              190.*kd1*kd2**3 + 157.*kd2**4)*omega1*omega2 + kd1**2*(10.*kd1**5 + 13.*kd1**4*kd2 - 105.*kd1**3*kd2**2 - 154.*kd1**2*kd2**3 + &
              28.*kd1*kd2**4 + 56.*kd2**5)*omega2**2)*swd**2) + 725594112.*grav*(216.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3 + (kd1*kd2**3*(123.*kd1**3 + &
              215.*kd1**2*kd2 + 186.*kd1*kd2**2 + 102.*kd2**3)*omega1**4 + kd2**2*(-48.*kd1**5 + 221.*kd1**4*kd2 + 483.*kd1**3*kd2**2 + &
              399.*kd1**2*kd2**3 + 197.*kd1*kd2**4 - 4.*kd2**5)*omega1**3*omega2 + kd1*kd2*(kd1 + kd2)*(209.*kd1**4 + 250.*kd1**3*kd2 + &
              456.*kd1**2*kd2**2 + 250.*kd1*kd2**3 + 209.*kd2**4)*omega1**2*omega2**2 + kd1**2*(-4.*kd1**5 + 197.*kd1**4*kd2 + 399.*kd1**3*kd2**2 + &
              483.*kd1**2*kd2**3 + 221.*kd1*kd2**4 - 48.*kd2**5)*omega1*omega2**3 + kd1**3*kd2*(102.*kd1**3 + 186.*kd1**2*kd2 + 215.*kd1*kd2**2 + &
              123.*kd2**3)*omega2**4)*swd**2) + 432.*grav*kd1**3*kd2**3*(kd1 + kd2)**3*(144.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4 + ((kd1 + kd2)*(9.*kd1**5 + &
              25.*kd1**4*kd2 + 56.*kd1**3*kd2**2 + 224.*kd1**2*kd2**3 + 232.*kd1*kd2**4 + 72.*kd2**5)*omega1**4 + (15.*kd1**6 + 122.*kd1**5*kd2 + &
              161.*kd1**4*kd2**2 + 372.*kd1**3*kd2**3 + 866.*kd1**2*kd2**4 + 682.*kd1*kd2**5 + 117.*kd2**6)*omega1**3*omega2 + (123.*kd1**6 + 542.*kd1**5*kd2 + &
              568.*kd1**4*kd2**2 + 332.*kd1**3*kd2**3 + 568.*kd1**2*kd2**4 + 542.*kd1*kd2**5 + 123.*kd2**6)*omega1**2*omega2**2 + (117.*kd1**6 + 682.*kd1**5*kd2 + &
              866.*kd1**4*kd2**2 + 372.*kd1**3*kd2**3 + 161.*kd1**2*kd2**4 + 122.*kd1*kd2**5 + 15.*kd2**6)*omega1*omega2**3 + (kd1 + kd2)*(72.*kd1**5 + &
              232.*kd1**4*kd2 + 224.*kd1**3*kd2**2 + 56.*kd1**2*kd2**3 + 25.*kd1*kd2**4 + 9.*kd2**5)*omega2**4)*swd**2) + &
              20155392.*grav*(48.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(5.*kd1**2 + 4.*kd1*kd2 + 5.*kd2**2) + (kd1*kd2**3*(306.*kd1**5 + 694.*kd1**4*kd2 + &
              777.*kd1**3*kd2**2 + 669.*kd1**2*kd2**3 + 458.*kd1*kd2**4 + 162.*kd2**5)*omega1**4 + kd2**2*(-312.*kd1**7 - 70.*kd1**6*kd2 + 1078.*kd1**5*kd2**2 + &
              1537.*kd1**4*kd2**3 + 1045.*kd1**3*kd2**4 + 697.*kd1**2*kd2**5 + 317.*kd1*kd2**6 - 4.*kd2**7)*omega1**3*omega2 + kd1*kd2*(kd1 + kd2)*(329.*kd1**6 + &
              268.*kd1**5*kd2 + 670.*kd1**4*kd2**2 + 1120.*kd1**3*kd2**3 + 670.*kd1**2*kd2**4 + 268.*kd1*kd2**5 + 329.*kd2**6)*omega1**2*omega2**2 + &
              kd1**2*(-4.*kd1**7 + 317.*kd1**6*kd2 + 697.*kd1**5*kd2**2 + 1045.*kd1**4*kd2**3 + 1537.*kd1**3*kd2**4 + 1078.*kd1**2*kd2**5 - 70.*kd1*kd2**6 - &
              312.*kd2**7)*omega1*omega2**3 + kd1**3*kd2*(162.*kd1**5 + 458.*kd1**4*kd2 + 669.*kd1**3*kd2**2 + 777.*kd1**2*kd2**3 + 694.*kd1*kd2**4 + &
              306.*kd2**5)*omega2**4)*swd**2) + 559872.*grav*(24.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(7.*kd1**4 + 4.*kd1**3*kd2 - 7.*kd1**2*kd2**2 + &
              4.*kd1*kd2**3 + 7.*kd2**4) + (kd1*kd2**3*(252.*kd1**7 + 1028.*kd1**6*kd2 + 848.*kd1**5*kd2**2 - 1220.*kd1**4*kd2**3 - 2163.*kd1**3*kd2**4 - &
              811.*kd1**2*kd2**5 + 142.*kd1*kd2**6 + 90.*kd2**7)*omega1**4 + kd2**2*(-144.*kd1**9 + 28.*kd1**8*kd2 + 1922.*kd1**7*kd2**2 + 2671.*kd1**6*kd2**3 - &
              1311.*kd1**5*kd2**4 - 4188.*kd1**4*kd2**5 - 1958.*kd1**3*kd2**6 + 169.*kd1**2*kd2**7 + 191.*kd1*kd2**8 + 4.*kd2**9)*omega1**3*omega2 + &
              kd1*kd2*(kd1 + kd2)*(179.*kd1**8 + 258.*kd1**7*kd2 - 701.*kd1**6*kd2**2 - 1230.*kd1**5*kd2**3 - 946.*kd1**4*kd2**4 - 1230.*kd1**3*kd2**5 - &
              701.*kd1**2*kd2**6 + 258.*kd1*kd2**7 + 179.*kd2**8)*omega1**2*omega2**2 + kd1**2*(4.*kd1**9 + 191.*kd1**8*kd2 + 169.*kd1**7*kd2**2 - &
              1958.*kd1**6*kd2**3 - 4188.*kd1**5*kd2**4 - 1311.*kd1**4*kd2**5 + 2671.*kd1**3*kd2**6 + 1922.*kd1**2*kd2**7 + 28.*kd1*kd2**8 - &
              144.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(90.*kd1**7 + 142.*kd1**6*kd2 - 811.*kd1**5*kd2**2 - 2163.*kd1**4*kd2**3 - 1220.*kd1**3*kd2**4 + &
              848.*kd1**2*kd2**5 + 1028.*kd1*kd2**6 + 252.*kd2**7)*omega2**4)*swd**2) - &
              72.*kd1**3*kd2**3*(kd1 + kd2)**3*swd*(2.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(3.*(3.*kd1 + 4.*kd2)**2*omega1**2 + (27.*kd1**2 + 80.*kd1*kd2 + &
              27.*kd2**2)*omega1*omega2 + 3.*(4.*kd1 + 3.*kd2)**2*omega2**2) + omega1*omega2*(omega1 + omega2)**2*(-6.*kd2**6*omega1*(omega1 - 2.*omega2) + &
              6.*kd1**6*(2.*omega1 - omega2)*omega2 + kd1*kd2**5*(41.*omega1**2 + 52.*omega1*omega2 + 4.*omega2**2) - 2.*kd1**3*kd2**3*(7.*omega1**2 - &
              76.*omega1*omega2 + 7.*omega2**2) + kd1**2*kd2**4*(14.*omega1**2 + 104.*omega1*omega2 + 11.*omega2**2) + kd1**4*kd2**2*(11.*omega1**2 + &
              104.*omega1*omega2 + 14.*omega2**2) + kd1**5*kd2*(4.*omega1**2 + 52.*omega1*omega2 + 41.*omega2**2))*swd**2) - &
              3359232.*swd*(2.*grav**2*kd1*kd2*(kd1 + kd2)*(66.*kd2**8*omega1**2 + 4.*kd1*kd2**7*omega1*(43.*omega1 - 15.*omega2) + 66.*kd1**8*omega2**2 + &
              4.*kd1**7*kd2*omega2*(-15.*omega1 + 43.*omega2) + 2.*kd1**5*kd2**3*(386.*omega1**2 + 483.*omega1*omega2 + 173.*omega2**2) + &
              kd1**6*kd2**2*(303.*omega1**2 + 181.*omega1*omega2 + 192.*omega2**2) + kd1**2*kd2**6*(192.*omega1**2 + 181.*omega1*omega2 + 303.*omega2**2) + &
              2.*kd1**3*kd2**5*(173.*omega1**2 + 483.*omega1*omega2 + 386.*omega2**2) + kd1**4*kd2**4*(783.*omega1**2 + 1558.*omega1*omega2 + 783.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(-4.*kd2**9*omega1**2 - 4.*kd1**9*omega2**2 - 4.*kd1**8*kd2*omega2*(37.*omega1 + 3.*omega2) - &
              4.*kd1*kd2**8*omega1*(3.*omega1 + 37.*omega2) + 4.*kd1**2*kd2**7*(48.*omega1**2 - 101.*omega1*omega2 + 39.*omega2**2) + &
              4.*kd1**7*kd2**2*(39.*omega1**2 - 101.*omega1*omega2 + 48.*omega2**2) + kd1**5*kd2**4*(337.*omega1**2 + 896.*omega1*omega2 + 231.*omega2**2) + &
              kd1**4*kd2**5*(231.*omega1**2 + 896.*omega1*omega2 + 337.*omega2**2) + kd1**6*kd2**3*(407.*omega1**2 - 98.*omega1*omega2 + 341.*omega2**2) + &
              kd1**3*kd2**6*(341.*omega1**2 - 98.*omega1*omega2 + 407.*omega2**2))*swd**2) - 93312.*swd*(2.*grav**2*kd1*kd2*(kd1 + kd2)**2*(24.*kd2**9*omega1**2 + &
              12.*kd1*kd2**8*omega1*(4.*omega1 - omega2) - 12.*kd1**8*kd2*(omega1 - 4.*omega2)*omega2 + 24.*kd1**9*omega2**2 + kd1**6*kd2**3*(277.*omega1**2 - &
              126.*omega1*omega2 - 749.*omega2**2) + kd1**7*kd2**2*(129.*omega1**2 + 80.*omega1*omega2 - 159.*omega2**2) + kd1**2*kd2**7*(-159.*omega1**2 + &
              80.*omega1*omega2 + 129.*omega2**2) + kd1**3*kd2**6*(-749.*omega1**2 - 126.*omega1*omega2 + 277.*omega2**2) - kd1**4*kd2**5*(1255.*omega1**2 + &
              1632.*omega1*omega2 + 565.*omega2**2) - kd1**5*kd2**4*(565.*omega1**2 + 1632.*omega1*omega2 + 1255.*omega2**2)) + &
              omega1*omega2*(omega1 + omega2)**2*(kd2**11*omega1**2 + kd1**10*kd2*(31.*omega1 - 3.*omega2)*omega2 + kd1**11*omega2**2 + &
              kd1*kd2**10*omega1*(-3.*omega1 + 31.*omega2) + kd1**3*kd2**8*(286.*omega1**2 - 246.*omega1*omega2 - 61.*omega2**2) + kd1**2*kd2**9*(91.*omega1**2 + &
              61.*omega1*omega2 - 24.*omega2**2) + kd1**9*kd2**2*(-24.*omega1**2 + 61.*omega1*omega2 + 91.*omega2**2) + kd1**8*kd2**3*(-61.*omega1**2 - &
              246.*omega1*omega2 + 286.*omega2**2) + kd1**4*kd2**7*(569.*omega1**2 - 398.*omega1*omega2 + 314.*omega2**2) + 2.*kd1**5*kd2**6*(479.*omega1**2 + &
              42.*omega1*omega2 + 470.*omega2**2) + 2.*kd1**6*kd2**5*(470.*omega1**2 + 42.*omega1*omega2 + 479.*omega2**2) + kd1**7*kd2**4*(314.*omega1**2 - &
              398.*omega1*omega2 + 569.*omega2**2))*swd**2) - 1296.*kd1*kd2*(kd1 + kd2)*swd*(4.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(15.*kd2**6*(8.*omega1**2 + &
              3.*omega1*omega2 + 3.*omega2**2) + 15.*kd1**6*(3.*omega1**2 + 3.*omega1*omega2 + 8.*omega2**2) + 2.*kd1**2*kd2**4*(261.*omega1**2 + &
              171.*omega1*omega2 + 14.*omega2**2) + 4.*kd1**3*kd2**3*(31.*omega1**2 + 52.*omega1*omega2 + 31.*omega2**2) + 2.*kd1*kd2**5*(240.*omega1**2 + &
              143.*omega1*omega2 + 57.*omega2**2) + 2.*kd1**5*kd2*(57.*omega1**2 + 143.*omega1*omega2 + 240.*omega2**2) + 2.*kd1**4*kd2**2*(14.*omega1**2 + &
              171.*omega1*omega2 + 261.*omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(kd2**10*omega1*(omega1 - 5.*omega2) + &
              kd1**10*omega2*(-5.*omega1 + omega2) - 4.*kd1**7*kd2**3*(15.*omega1**2 + 157.*omega1*omega2 - 128.*omega2**2) + kd1**4*kd2**6*(853.*omega1**2 - &
              1603.*omega1*omega2 - 114.*omega2**2) + 4.*kd1**3*kd2**7*(128.*omega1**2 - 157.*omega1*omega2 - 15.*omega2**2) + 4.*kd1**9*kd2*(omega1**2 - &
              6.*omega1*omega2 - omega2**2) - 4.*kd1*kd2**9*(omega1**2 + 6.*omega1*omega2 - omega2**2) + kd1**2*kd2**8*(154.*omega1**2 - 133.*omega1*omega2 + &
              19.*omega2**2) + 4.*kd1**5*kd2**5*(121.*omega1**2 - 533.*omega1*omega2 + 121.*omega2**2) + kd1**8*kd2**2*(19.*omega1**2 - 133.*omega1*omega2 + &
              154.*omega2**2) + kd1**6*kd2**4*(-114.*omega1**2 - 1603.*omega1*omega2 + 853.*omega2**2))*swd**2) + &
              15552.*grav*kd1*kd2*(kd1 + kd2)*(24.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(6.*kd1**4 + 4.*kd1**3*kd2 - 33.*kd1**2*kd2**2 + 4.*kd1*kd2**3 + &
              6.*kd2**4) + (3.*kd1**10*omega2**2*(7.*omega1**2 + 7.*omega1*omega2 + 4.*omega2**2) + 3.*kd2**10*omega1**2*(4.*omega1**2 + 7.*omega1*omega2 + &
              7.*omega2**2) + 4.*kd1*kd2**9*omega1*(16.*omega1**3 + 44.*omega1**2*omega2 + 13.*omega1*omega2**2 - &
              3.*omega2**3) + 4.*kd1**9*kd2*omega2*(-3.*omega1**3 + 13.*omega1**2*omega2 + 44.*omega1*omega2**2 + 16.*omega2**3) + &
              4.*kd1**3*kd2**7*(44.*omega1**4 + 96.*omega1**3*omega2 + 176.*omega1**2*omega2**2 + 125.*omega1*omega2**3 + 33.*omega2**4) + &
              4.*kd1**7*kd2**3*(33.*omega1**4 + 125.*omega1**3*omega2 + 176.*omega1**2*omega2**2 + 96.*omega1*omega2**3 + 44.*omega2**4) + &
              kd1**2*kd2**8*(167.*omega1**4 + 373.*omega1**3*omega2 - 66.*omega1**2*omega2**2 + 32.*omega1*omega2**3 + 60.*omega2**4) - &
              2.*kd1**4*kd2**6*(149.*omega1**4 + 270.*omega1**3*omega2 - 920.*omega1**2*omega2**2 - 34.*omega1*omega2**3 + 144.*omega2**4) - &
              2.*kd1**6*kd2**4*(144.*omega1**4 - 34.*omega1**3*omega2 - 920.*omega1**2*omega2**2 + 270.*omega1*omega2**3 + 149.*omega2**4) + &
              kd1**8*kd2**2*(60.*omega1**4 + 32.*omega1**3*omega2 - 66.*omega1**2*omega2**2 + 373.*omega1*omega2**3 + 167.*omega2**4) - &
              4.*kd1**5*kd2**5*(179.*omega1**4 + 270.*omega1**3*omega2 - 525.*omega1**2*omega2**2 + 270.*omega1*omega2**3 + 179.*omega2**4))*swd**2)))/ &
              (4.*kd1**2*(-36. + kd1**2)**2*kd2**2*(-36. + kd2**2)**2*omega1*omega2*(omega1 + omega2)*swd**3*(-36.*grav*(kd1 + kd2)**2*(12. + &
              (kd1 + kd2)**2)*(108. + (kd1 + kd2)**2) + (36. + (kd1 + kd2)**2)*(1296. + (kd1 + kd2)**2*(504. + (kd1 + kd2)**2))*(omega1 + omega2)**2*swd))
    !
end function velsp33
!
real function etasub4()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic surface elevation for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   surface transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'etasub4')
    !
    etasub4 = ((kd1 - kd2)**2*(-1152921504606846976.*grav*kd1*kd2*(6.*kd2**2*omega1**2 + kd1*kd2*omega1*omega2 + 6.*kd1**2*omega2**2)*swd + &
              64.*grav*kd1**7*(kd1 - kd2)**6*kd2**7*((2.*kd1**2 - 8.*kd1*kd2 + 5.*kd2**2)*omega1**2 + (kd1**2 - 21.*kd1*kd2 + kd2**2)*omega1*omega2 + &
              (5.*kd1**2 - 8.*kd1*kd2 + 2.*kd2**2)*omega2**2)*swd - 18014398509481984.*grav*kd1*kd2*(11.*kd2**4*omega1**2 + 2.*kd1*kd2**3*omega1*(34.*omega1 - &
              13.*omega2) + 11.*kd1**4*omega2**2 + 2.*kd1**3*kd2*omega2*(-13.*omega1 + 34.*omega2) + 2.*kd1**2*kd2**2*(29.*omega1**2 + 32.*omega1*omega2 + &
              29.*omega2**2))*swd + 281474976710656.*grav*kd1*kd2*(27.*kd2**6*omega1**2 + kd1**5*kd2*(omega1 - 250.*omega2)*omega2 + 27.*kd1**6*omega2**2 + &
              kd1*kd2**5*omega1*(-250.*omega1 + omega2) + kd1**3*kd2**3*(-824.*omega1**2 + 339.*omega1*omega2 - 824.*omega2**2) + kd1**2*kd2**4*(348.*omega1**2 - &
              157.*omega1*omega2 + 42.*omega2**2) + kd1**4*kd2**2*(42.*omega1**2 - 157.*omega1*omega2 + 348.*omega2**2))*swd + &
              4096.*grav*kd1**5*(kd1 - kd2)**4*kd2**5*(kd1**3*kd2**3*(-324.*omega1**2 + 35.*omega1*omega2 - 324.*omega2**2) + kd2**6*(57.*omega1**2 - &
              2.*omega1*omega2 + 12.*omega2**2) - kd1*kd2**5*(250.*omega1**2 + 11.*omega1*omega2 + 56.*omega2**2) + kd1**6*(12.*omega1**2 - 2.*omega1*omega2 + &
              57.*omega2**2) + kd1**2*kd2**4*(434.*omega1**2 - 37.*omega1*omega2 + 128.*omega2**2) - kd1**5*kd2*(56.*omega1**2 + 11.*omega1*omega2 + 250.*omega2**2) + &
              kd1**4*kd2**2*(128.*omega1**2 - 37.*omega1*omega2 + 434.*omega2**2))*swd + 4398046511104.*grav*kd1*kd2*(10.*kd2**8*omega1**2 + &
              4.*kd1*kd2**7*omega1*(60.*omega1 - 13.*omega2) + 10.*kd1**8*omega2**2 + 4.*kd1**7*kd2*omega2*(-13.*omega1 + 60.*omega2) + kd1**6*kd2**2*(22.*omega1**2 + &
              349.*omega1*omega2 - 421.*omega2**2) + kd1**2*kd2**6*(-421.*omega1**2 + 349.*omega1*omega2 + 22.*omega2**2) - kd1**3*kd2**5*(1182.*omega1**2 + &
              777.*omega1*omega2 + 380.*omega2**2) + kd1**4*kd2**4*(853.*omega1**2 + 960.*omega1*omega2 + 853.*omega2**2) - kd1**5*kd2**3*(380.*omega1**2 + &
              777.*omega1*omega2 + 1182.*omega2**2))*swd + 262144.*grav*kd1**3*(kd1 - kd2)**2*kd2**3*(kd1**5*kd2**5*(-4752.*omega1**2 + 739.*omega1*omega2 - &
              4752.*omega2**2) + kd2**10*(63.*omega1**2 + omega1*omega2 + 2.*omega2**2) - kd1*kd2**9*(468.*omega1**2 + 15.*omega1*omega2 + 32.*omega2**2) + &
              kd1**10*(2.*omega1**2 + omega1*omega2 + 63.*omega2**2) + kd1**2*kd2**8*(1283.*omega1**2 + 49.*omega1*omega2 + 357.*omega2**2) - kd1**9*kd2*(32.*omega1**2 + &
              15.*omega1*omega2 + 468.*omega2**2) - 2.*kd1**3*kd2**7*(1209.*omega1**2 - 60.*omega1*omega2 + 758.*omega2**2) - 2.*kd1**7*kd2**3*(758.*omega1**2 - &
              60.*omega1*omega2 + 1209.*omega2**2) + kd1**8*kd2**2*(357.*omega1**2 + 49.*omega1*omega2 + 1283.*omega2**2) + kd1**4*kd2**6*(4044.*omega1**2 - &
              539.*omega1*omega2 + 3438.*omega2**2) + kd1**6*kd2**4*(3438.*omega1**2 - 539.*omega1*omega2 + 4044.*omega2**2))*swd - &
              68719476736.*grav*kd1*kd2*(24.*kd2**10*omega1**2 + 24.*kd1**10*omega2**2 - kd1*kd2**9*omega1*(276.*omega1 + omega2) - kd1**9*kd2*omega2*(omega1 + &
              276.*omega2) - 2.*kd1**8*kd2**2*(47.*omega1**2 + 35.*omega1*omega2 - 430.*omega2**2) + 2.*kd1**2*kd2**8*(430.*omega1**2 - 35.*omega1*omega2 - &
              47.*omega2**2) - 2.*kd1**3*kd2**7*(790.*omega1**2 - 21.*omega1*omega2 + 648.*omega2**2) - 2.*kd1**7*kd2**3*(648.*omega1**2 - 21.*omega1*omega2 + &
              790.*omega2**2) + 2.*kd1**6*kd2**4*(2402.*omega1**2 + 242.*omega1*omega2 + 1577.*omega2**2) + 2.*kd1**4*kd2**6*(1577.*omega1**2 + &
              242.*omega1*omega2 + 2402.*omega2**2) - kd1**5*kd2**5*(5252.*omega1**2 + 883.*omega1*omega2 + 5252.*omega2**2))*swd + &
              1073741824.*grav*kd1*kd2*(kd2**12*omega1**2 + 26.*kd1**11*kd2*(omega1 - 6.*omega2)*omega2 + kd1**12*omega2**2 + &
              26.*kd1*kd2**11*omega1*(-6.*omega1 + omega2) - 2.*kd1**10*kd2**2*(47.*omega1**2 + 103.*omega1*omega2 - 252.*omega2**2) + &
              2.*kd1**2*kd2**10*(252.*omega1**2 - 103.*omega1*omega2 - 47.*omega2**2) + 4.*kd1**3*kd2**9*(376.*omega1**2 + 124.*omega1*omega2 + 89.*omega2**2) + &
              4.*kd1**9*kd2**3*(89.*omega1**2 + 124.*omega1*omega2 + 376.*omega2**2) + 2.*kd1**5*kd2**7*(4062.*omega1**2 + 953.*omega1*omega2 + 458.*omega2**2) - &
              8.*kd1**6*kd2**6*(491.*omega1**2 + 316.*omega1*omega2 + 491.*omega2**2) - kd1**4*kd2**8*(6775.*omega1**2 + 952.*omega1*omega2 + 497.*omega2**2) + &
              2.*kd1**7*kd2**5*(458.*omega1**2 + 953.*omega1*omega2 + 4062.*omega2**2) - kd1**8*kd2**4*(497.*omega1**2 + 952.*omega1*omega2 + &
              6775.*omega2**2))*swd + 16777216.*grav*kd1*kd2*(3.*kd2**14*omega1**2 + 3.*kd1**14*omega2**2 - kd1*kd2**13*omega1*(42.*omega1 + omega2) - &
              kd1**13*kd2*omega2*(omega1 + 42.*omega2) + kd1**2*kd2**12*(168.*omega1**2 + 7.*omega1*omega2 - 2.*omega2**2) + kd1**12*kd2**2*(-2.*omega1**2 + &
              7.*omega1*omega2 + 168.*omega2**2) - kd1**3*kd2**11*(740.*omega1**2 + 25.*omega1*omega2 + 344.*omega2**2) - kd1**11*kd2**3*(344.*omega1**2 + &
              25.*omega1*omega2 + 740.*omega2**2) + kd1**4*kd2**10*(3335.*omega1**2 + 363.*omega1*omega2 + 1628.*omega2**2) + kd1**6*kd2**8*(9326.*omega1**2 + &
              3749.*omega1*omega2 + 2015.*omega2**2) - kd1**5*kd2**9*(8074.*omega1**2 + 1665.*omega1*omega2 + 2370.*omega2**2) + kd1**10*kd2**4*(1628.*omega1**2 + &
              363.*omega1*omega2 + 3335.*omega2**2) - kd1**7*kd2**7*(4904.*omega1**2 + 4855.*omega1*omega2 + 4904.*omega2**2) - kd1**9*kd2**5*(2370.*omega1**2 + &
              1665.*omega1*omega2 + 8074.*omega2**2) + kd1**8*kd2**6*(2015.*omega1**2 + 3749.*omega1*omega2 + 9326.*omega2**2))*swd + &
              1152921504606846976.*omega1*omega2*(kd2**2*omega1**2 + 8.*kd1*kd2*omega1*omega2 + kd1**2*omega2**2)*swd**2 - &
              kd1**7*(kd1 - kd2)**6*kd2**7*omega1*omega2*(kd2*(kd1 + 8.*kd2)*omega1**2 - 2.*(2.*kd1 + kd2)*(kd1 + 2.*kd2)*omega1*omega2 + &
              kd1*(8.*kd1 + kd2)*omega2**2)*swd**2 + 18014398509481984.*(384.*grav**2*kd1**3*kd2**3 + omega1*omega2*(-6.*kd2**4*omega1**2 - &
              6.*kd1**4*omega2**2 + 4.*kd1**3*kd2*omega2*(27.*omega1 + omega2) + 4.*kd1*kd2**3*omega1*(omega1 + 27.*omega2) + kd1**2*kd2**2*(-37.*omega1**2 + &
              90.*omega1*omega2 - 37.*omega2**2))*swd**2) - 64.*kd1**5*(kd1 - kd2)**4*kd2**5*(192.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2 - &
              omega1*omega2*(4.*kd2**6*omega1*(4.*omega1 - 19.*omega2) + 4.*kd1**6*omega2*(-19.*omega1 + 4.*omega2) + 2.*kd1**4*kd2**2*(9.*omega1**2 - &
              170.*omega1*omega2 - 259.*omega2**2) - 2.*kd1**2*kd2**4*(259.*omega1**2 + 170.*omega1*omega2 - 9.*omega2**2) + kd1*kd2**5*(322.*omega1**2 + &
              314.*omega1*omega2 + 3.*omega2**2) - 2.*kd1**3*kd2**3*(4.*omega1**2 - 269.*omega1*omega2 + 4.*omega2**2) + kd1**5*kd2*(3.*omega1**2 + &
              314.*omega1*omega2 + 322.*omega2**2))*swd**2) + 281474976710656.*(64.*grav**2*kd1**3*kd2**3*(11.*kd1**2 + 18.*kd1*kd2 + 11.*kd2**2) + &
              omega1*omega2*(15.*kd2**6*omega1**2 + 15.*kd1**6*omega2**2 - 4.*kd1**5*kd2*omega2*(16.*omega1 + 5.*omega2) - 4.*kd1*kd2**5*omega1*(5.*omega1 + &
              16.*omega2) + kd1**4*kd2**2*(-279.*omega1**2 + 1328.*omega1*omega2 - 417.*omega2**2) + kd1**2*kd2**4*(-417.*omega1**2 + 1328.*omega1*omega2 - &
              279.*omega2**2) + 2.*kd1**3*kd2**3*(533.*omega1**2 - 34.*omega1*omega2 + 533.*omega2**2))*swd**2) - &
              4398046511104.*(64.*grav**2*kd1**3*kd2**3*(27.*kd1**4 - 152.*kd1**3*kd2 + 134.*kd1**2*kd2**2 - 152.*kd1*kd2**3 + 27.*kd2**4) + &
              omega1*omega2*(20.*kd2**8*omega1**2 + 20.*kd1**7*kd2*(11.*omega1 - 2.*omega2)*omega2 + 20.*kd1**8*omega2**2 + &
              20.*kd1*kd2**7*omega1*(-2.*omega1 + 11.*omega2) + 2.*kd1**3*kd2**5*(681.*omega1**2 + 1658.*omega1*omega2 - 89.*omega2**2) + &
              6.*kd1**4*kd2**4*(213.*omega1**2 - 2239.*omega1*omega2 + 213.*omega2**2) - kd1**2*kd2**6*(1108.*omega1**2 + 970.*omega1*omega2 + 315.*omega2**2) + &
              2.*kd1**5*kd2**3*(-89.*omega1**2 + 1658.*omega1*omega2 + 681.*omega2**2) - kd1**6*kd2**2*(315.*omega1**2 + 970.*omega1*omega2 + &
              1108.*omega2**2))*swd**2) - 68719476736.*(128.*grav**2*kd1**3*kd2**3*(5.*kd1**6 + 10.*kd1**5*kd2 - 72.*kd1**4*kd2**2 + 78.*kd1**3*kd2**3 - &
              72.*kd1**2*kd2**4 + 10.*kd1*kd2**5 + 5.*kd2**6) + omega1*omega2*(-15.*kd2**10*omega1**2 + 8.*kd1*kd2**9*omega1*(5.*omega1 - 13.*omega2) - &
              15.*kd1**10*omega2**2 + 8.*kd1**9*kd2*omega2*(-13.*omega1 + 5.*omega2) + kd1**4*kd2**6*(-5322.*omega1**2 + 3682.*omega1*omega2 - 5431.*omega2**2) + &
              kd1**6*kd2**4*(-5431.*omega1**2 + 3682.*omega1*omega2 - 5322.*omega2**2) + kd1**2*kd2**8*(546.*omega1**2 + 3280.*omega1*omega2 - 315.*omega2**2) + &
              4.*kd1**7*kd2**3*(783.*omega1**2 - 2445.*omega1*omega2 + 409.*omega2**2) + kd1**8*kd2**2*(-315.*omega1**2 + 3280.*omega1*omega2 + 546.*omega2**2) + &
              4.*kd1**3*kd2**7*(409.*omega1**2 - 2445.*omega1*omega2 + 783.*omega2**2) + 2.*kd1**5*kd2**5*(2215.*omega1**2 + 1362.*omega1*omega2 + &
              2215.*omega2**2))*swd**2) - 4096.*kd1**3*(kd1 - kd2)**2*kd2**3*(64.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(10.*kd1**4 - 28.*kd1**3*kd2 + &
              19.*kd1**2*kd2**2 - 28.*kd1*kd2**3 + 10.*kd2**4) + omega1*omega2*(4.*kd1**10*omega2*(13.*omega1 + 2.*omega2) + 4.*kd2**10*omega1*(2.*omega1 + &
              13.*omega2) + 2.*kd1**8*kd2**2*(6.*omega1**2 + 1682.*omega1*omega2 - 315.*omega2**2) + kd1*kd2**9*(32.*omega1**2 - 494.*omega1*omega2 + &
              3.*omega2**2) + 2.*kd1**2*kd2**8*(-315.*omega1**2 + 1682.*omega1*omega2 + 6.*omega2**2) + kd1**9*kd2*(3.*omega1**2 - 494.*omega1*omega2 + &
              32.*omega2**2) + kd1**3*kd2**7*(1432.*omega1**2 - 11636.*omega1*omega2 + 33.*omega2**2) - 4.*kd1**4*kd2**6*(713.*omega1**2 - 6066.*omega1*omega2 + &
              429.*omega2**2) - 4.*kd1**6*kd2**4*(429.*omega1**2 - 6066.*omega1*omega2 + 713.*omega2**2) + kd1**7*kd2**3*(33.*omega1**2 - 11636.*omega1*omega2 + &
              1432.*omega2**2) + kd1**5*kd2**5*(3957.*omega1**2 - 31658.*omega1*omega2 + 3957.*omega2**2))*swd**2) - &
              262144.*kd1*kd2*(64.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(3.*kd1**8 - 20.*kd1**7*kd2 + 24.*kd1**6*kd2**2 - 88.*kd1**5*kd2**3 + 213.*kd1**4*kd2**4 - &
              88.*kd1**3*kd2**5 + 24.*kd1**2*kd2**6 - 20.*kd1*kd2**7 + 3.*kd2**8) - omega1*omega2*(-4.*kd1**14*omega1*omega2 - 4.*kd2**14*omega1*omega2 + &
              2.*kd1**4*kd2**10*(871.*omega1**2 - 15542.*omega1*omega2 - 985.*omega2**2) - 2.*kd1**10*kd2**4*(985.*omega1**2 + 15542.*omega1*omega2 - 871.*omega2**2) + &
              2.*kd1**12*kd2**2*(omega1**2 - 462.*omega1*omega2 - 325.*omega2**2) + 2.*kd1**2*kd2**12*(-325.*omega1**2 - 462.*omega1*omega2 + omega2**2) + &
              kd1*kd2**13*(94.*omega1**2 + 78.*omega1*omega2 + omega2**2) + 6.*kd1**3*kd2**11*(39.*omega1**2 + 1347.*omega1*omega2 + 13.*omega2**2) + &
              6.*kd1**11*kd2**3*(13.*omega1**2 + 1347.*omega1*omega2 + 39.*omega2**2) + kd1**13*kd2*(omega1**2 + 78.*omega1*omega2 + 94.*omega2**2) + &
              2.*kd1**9*kd2**5*(3577.*omega1**2 + 30645.*omega1*omega2 + 124.*omega2**2) - 2.*kd1**8*kd2**6*(5163.*omega1**2 + 38446.*omega1*omega2 + 2783.*omega2**2) + &
              2.*kd1**5*kd2**9*(124.*omega1**2 + 30645.*omega1*omega2 + 3577.*omega2**2) + 2.*kd1**7*kd2**7*(4451.*omega1**2 + 39511.*omega1*omega2 + 4451.*omega2**2) - &
              2.*kd1**6*kd2**8*(2783.*omega1**2 + 38446.*omega1*omega2 + 5163.*omega2**2))*swd**2) + 1073741824.*(128.*grav**2*kd1**3*kd2**3*(12.*kd1**8 - &
              116.*kd1**7*kd2 + 316.*kd1**6*kd2**2 - 410.*kd1**5*kd2**3 + 407.*kd1**4*kd2**4 - 410.*kd1**3*kd2**5 + 316.*kd1**2*kd2**6 - 116.*kd1*kd2**7 + 12.*kd2**8) - &
              omega1*omega2*(6.*kd2**12*omega1**2 + 6.*kd1**12*omega2**2 - 4.*kd1**11*kd2*omega2*(29.*omega1 + 5.*omega2) - 4.*kd1*kd2**11*omega1*(5.*omega1 + 29.*omega2) + &
              kd1**2*kd2**10*(285.*omega1**2 + 242.*omega1*omega2 + 279.*omega2**2) + kd1**10*kd2**2*(279.*omega1**2 + 242.*omega1*omega2 + 285.*omega2**2) - &
              4.*kd1**3*kd2**9*(607.*omega1**2 - 103.*omega1*omega2 + 298.*omega2**2) - 4.*kd1**9*kd2**3*(298.*omega1**2 - 103.*omega1*omega2 + 607.*omega2**2) + &
              2.*kd1**4*kd2**8*(2746.*omega1**2 + 8662.*omega1*omega2 + 965.*omega2**2) - 4.*kd1**5*kd2**7*(2872.*omega1**2 + 16027.*omega1*omega2 + 2057.*omega2**2) + &
              2.*kd1**8*kd2**4*(965.*omega1**2 + 8662.*omega1*omega2 + 2746.*omega2**2) - 4.*kd1**7*kd2**5*(2057.*omega1**2 + 16027.*omega1*omega2 + 2872.*omega2**2) + &
              3.*kd1**6*kd2**6*(5411.*omega1**2 + 29962.*omega1*omega2 + 5411.*omega2**2))*swd**2) - 16777216.*(64.*grav**2*kd1**3*(kd1 - kd2)**2*kd2**3*(kd1**8 - &
              48.*kd1**7*kd2 + 168.*kd1**6*kd2**2 - 36.*kd1**5*kd2**3 - 217.*kd1**4*kd2**4 - 36.*kd1**3*kd2**5 + 168.*kd1**2*kd2**6 - 48.*kd1*kd2**7 + kd2**8) - &
              omega1*omega2*(kd2**14*omega1**2 + kd1**14*omega2**2 - 4.*kd1**13*kd2*omega2*(12.*omega1 + omega2) - 4.*kd1*kd2**13*omega1*(omega1 + 12.*omega2) + &
              kd1**4*kd2**10*(1509.*omega1**2 - 388.*omega1*omega2 - 1875.*omega2**2) + 8.*kd1**9*kd2**5*(137.*omega1**2 + 555.*omega1*omega2 - 1210.*omega2**2) + &
              kd1**6*kd2**8*(12705.*omega1**2 + 19896.*omega1*omega2 - 173.*omega2**2) + kd1**2*kd2**12*(83.*omega1**2 + 1056.*omega1*omega2 - 37.*omega2**2) + &
              2.*kd1**11*kd2**3*(345.*omega1**2 - 1534.*omega1*omega2 + 81.*omega2**2) + kd1**12*kd2**2*(-37.*omega1**2 + 1056.*omega1*omega2 + 83.*omega2**2) + &
              8.*kd1**5*kd2**9*(-1210.*omega1**2 + 555.*omega1*omega2 + 137.*omega2**2) - 22.*kd1**7*kd2**7*(189.*omega1**2 + 2014.*omega1*omega2 + 189.*omega2**2) + &
              2.*kd1**3*kd2**11*(81.*omega1**2 - 1534.*omega1*omega2 + 345.*omega2**2) + kd1**10*kd2**4*(-1875.*omega1**2 - 388.*omega1*omega2 + 1509.*omega2**2) + &
              kd1**8*kd2**6*(-173.*omega1**2 + 19896.*omega1*omega2 + 12705.*omega2**2))*swd**2)))/ &
              (2.*(-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*omega2*swd**2*(-64.*grav*(64. + (kd1 - kd2)**2)*(4096. + &
              (384. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2 + (16777216. + (7340032. + (286720. + (1792. + &
              (kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function etasub4
!
real function etasup4()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic surface elevation for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   surface transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'etasup4')
    !
    etasup4 = ((kd1 + kd2)**2*(-1152921504606846976.*grav*kd1*kd2*(6.*kd2**2*omega1**2 + kd1*kd2*omega1*omega2 + 6.*kd1**2*omega2**2)*swd + &
              64.*grav*kd1**7*kd2**7*(kd1 + kd2)**6*((2.*kd1**2 + 8.*kd1*kd2 + 5.*kd2**2)*omega1**2 - (kd1**2 + 21.*kd1*kd2 + kd2**2)*omega1*omega2 + &
              (5.*kd1**2 + 8.*kd1*kd2 + 2.*kd2**2)*omega2**2)*swd + 4398046511104.*grav*kd1*kd2*(kd2**2*(22.*kd1**6 + 380.*kd1**5*kd2 + 853.*kd1**4*kd2**2 + &
              1182.*kd1**3*kd2**3 - 421.*kd1**2*kd2**4 - 240.*kd1*kd2**5 + 10.*kd2**6)*omega1**2 - kd1*kd2*(kd1 + kd2)**2*(52.*kd1**4 + 245.*kd1**3*kd2 + &
              235.*kd1**2*kd2**2 + 245.*kd1*kd2**3 + 52.*kd2**4)*omega1*omega2 + kd1**2*(10.*kd1**6 - 240.*kd1**5*kd2 - 421.*kd1**4*kd2**2 + &
              1182.*kd1**3*kd2**3 + 853.*kd1**2*kd2**4 + 380.*kd1*kd2**5 + 22.*kd2**6)*omega2**2)*swd - 18014398509481984.*grav*kd1*kd2*(11.*kd2**4*omega1**2 + &
              11.*kd1**4*omega2**2 - 2.*kd1*kd2**3*omega1*(34.*omega1 + 13.*omega2) - 2.*kd1**3*kd2*omega2*(13.*omega1 + 34.*omega2) + &
              2.*kd1**2*kd2**2*(29.*omega1**2 - 32.*omega1*omega2 + 29.*omega2**2))*swd + 4096.*grav*kd1**5*kd2**5*(kd1 + kd2)**4*(kd2**6*(57.*omega1**2 + &
              2.*omega1*omega2 + 12.*omega2**2) + kd1*kd2**5*(250.*omega1**2 - 11.*omega1*omega2 + 56.*omega2**2) + kd1**6*(12.*omega1**2 + 2.*omega1*omega2 + &
              57.*omega2**2) + kd1**2*kd2**4*(434.*omega1**2 + 37.*omega1*omega2 + 128.*omega2**2) + kd1**5*kd2*(56.*omega1**2 - 11.*omega1*omega2 + 250.*omega2**2) + &
              kd1**3*kd2**3*(324.*omega1**2 + 35.*omega1*omega2 + 324.*omega2**2) + kd1**4*kd2**2*(128.*omega1**2 + 37.*omega1*omega2 + &
              434.*omega2**2))*swd + 281474976710656.*grav*kd1*kd2*(27.*kd2**6*omega1**2 + 27.*kd1**6*omega2**2 + kd1*kd2**5*omega1*(250.*omega1 + omega2) + &
              kd1**5*kd2*omega2*(omega1 + 250.*omega2) + kd1**2*kd2**4*(348.*omega1**2 + 157.*omega1*omega2 + 42.*omega2**2) + kd1**4*kd2**2*(42.*omega1**2 + &
              157.*omega1*omega2 + 348.*omega2**2) + kd1**3*kd2**3*(824.*omega1**2 + 339.*omega1*omega2 + 824.*omega2**2))*swd + &
              1073741824.*grav*kd1*kd2*(kd2**12*omega1**2 + kd1**12*omega2**2 + 26.*kd1*kd2**11*omega1*(6.*omega1 + omega2) + &
              26.*kd1**11*kd2*omega2*(omega1 + 6.*omega2) + kd1**8*kd2**4*(-497.*omega1**2 + 952.*omega1*omega2 - 6775.*omega2**2) + &
              kd1**4*kd2**8*(-6775.*omega1**2 + 952.*omega1*omega2 - 497.*omega2**2) + 2.*kd1**2*kd2**10*(252.*omega1**2 + 103.*omega1*omega2 - &
              47.*omega2**2) - 4.*kd1**3*kd2**9*(376.*omega1**2 - 124.*omega1*omega2 + 89.*omega2**2) + 2.*kd1**10*kd2**2*(-47.*omega1**2 + 103.*omega1*omega2 + &
              252.*omega2**2) - 4.*kd1**9*kd2**3*(89.*omega1**2 - 124.*omega1*omega2 + 376.*omega2**2) - 2.*kd1**5*kd2**7*(4062.*omega1**2 - 953.*omega1*omega2 + &
              458.*omega2**2) - 8.*kd1**6*kd2**6*(491.*omega1**2 - 316.*omega1*omega2 + 491.*omega2**2) - 2.*kd1**7*kd2**5*(458.*omega1**2 - 953.*omega1*omega2 + &
              4062.*omega2**2))*swd + 262144.*grav*kd1**3*kd2**3*(kd1 + kd2)**2*(kd2**10*(63.*omega1**2 - omega1*omega2 + 2.*omega2**2) + kd1*kd2**9*(468.*omega1**2 - &
              15.*omega1*omega2 + 32.*omega2**2) + kd1**10*(2.*omega1**2 - omega1*omega2 + 63.*omega2**2) + kd1**2*kd2**8*(1283.*omega1**2 - 49.*omega1*omega2 + &
              357.*omega2**2) + kd1**9*kd2*(32.*omega1**2 - 15.*omega1*omega2 + 468.*omega2**2) + 2.*kd1**3*kd2**7*(1209.*omega1**2 + 60.*omega1*omega2 + &
              758.*omega2**2) + 2.*kd1**7*kd2**3*(758.*omega1**2 + 60.*omega1*omega2 + 1209.*omega2**2) + kd1**8*kd2**2*(357.*omega1**2 - 49.*omega1*omega2 + &
              1283.*omega2**2) + kd1**4*kd2**6*(4044.*omega1**2 + 539.*omega1*omega2 + 3438.*omega2**2) + kd1**6*kd2**4*(3438.*omega1**2 + 539.*omega1*omega2 + &
              4044.*omega2**2) + kd1**5*kd2**5*(4752.*omega1**2 + 739.*omega1*omega2 + 4752.*omega2**2))*swd - 68719476736.*grav*kd1*kd2*(24.*kd2**10*omega1**2 + &
              kd1*kd2**9*omega1*(276.*omega1 - omega2) + 24.*kd1**10*omega2**2 + kd1**9*kd2*omega2*(-omega1 + 276.*omega2) + 2.*kd1**2*kd2**8*(430.*omega1**2 + &
              35.*omega1*omega2 - 47.*omega2**2) + 2.*kd1**8*kd2**2*(-47.*omega1**2 + 35.*omega1*omega2 + 430.*omega2**2) + 2.*kd1**3*kd2**7*(790.*omega1**2 + &
              21.*omega1*omega2 + 648.*omega2**2) + 2.*kd1**7*kd2**3*(648.*omega1**2 + 21.*omega1*omega2 + 790.*omega2**2) + 2.*kd1**6*kd2**4*(2402.*omega1**2 - &
              242.*omega1*omega2 + 1577.*omega2**2) + 2.*kd1**4*kd2**6*(1577.*omega1**2 - 242.*omega1*omega2 + 2402.*omega2**2) + kd1**5*kd2**5*(5252.*omega1**2 - &
              883.*omega1*omega2 + 5252.*omega2**2))*swd + 16777216.*grav*kd1*kd2*(3.*kd2**14*omega1**2 + kd1*kd2**13*omega1*(42.*omega1 - omega2) + &
              3.*kd1**14*omega2**2 + kd1**13*kd2*omega2*(-omega1 + 42.*omega2) + kd1**2*kd2**12*(168.*omega1**2 - 7.*omega1*omega2 - 2.*omega2**2) + &
              kd1**12*kd2**2*(-2.*omega1**2 - 7.*omega1*omega2 + 168.*omega2**2) + kd1**3*kd2**11*(740.*omega1**2 - 25.*omega1*omega2 + 344.*omega2**2) + &
              kd1**11*kd2**3*(344.*omega1**2 - 25.*omega1*omega2 + 740.*omega2**2) + kd1**4*kd2**10*(3335.*omega1**2 - 363.*omega1*omega2 + 1628.*omega2**2) + &
              kd1**6*kd2**8*(9326.*omega1**2 - 3749.*omega1*omega2 + 2015.*omega2**2) + kd1**5*kd2**9*(8074.*omega1**2 - 1665.*omega1*omega2 + 2370.*omega2**2) + &
              kd1**10*kd2**4*(1628.*omega1**2 - 363.*omega1*omega2 + 3335.*omega2**2) + kd1**7*kd2**7*(4904.*omega1**2 - 4855.*omega1*omega2 + 4904.*omega2**2) + &
              kd1**9*kd2**5*(2370.*omega1**2 - 1665.*omega1*omega2 + 8074.*omega2**2) + kd1**8*kd2**6*(2015.*omega1**2 - 3749.*omega1*omega2 + &
              9326.*omega2**2))*swd + 1152921504606846976.*omega1*omega2*(kd2**2*omega1**2 + 8.*kd1*kd2*omega1*omega2 + kd1**2*omega2**2)*swd**2 + &
              kd1**7*kd2**7*(kd1 + kd2)**6*omega1*omega2*(4.*kd2**2*omega1*(2.*omega1 + omega2) + 4.*kd1**2*omega2*(omega1 + 2.*omega2) - kd1*kd2*(omega1**2 + &
              10.*omega1*omega2 + omega2**2))*swd**2 - 64.*kd1**5*kd2**5*(kd1 + kd2)**4*(192.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4 + &
              omega1*omega2*(4.*kd1**6*omega2*(19.*omega1 + 4.*omega2) + 4.*kd2**6*omega1*(4.*omega1 + 19.*omega2) + kd1**5*kd2*(-3.*omega1**2 + 314.*omega1*omega2 - &
              322.*omega2**2) + 2.*kd1**4*kd2**2*(9.*omega1**2 + 170.*omega1*omega2 - 259.*omega2**2) + kd1*kd2**5*(-322.*omega1**2 + 314.*omega1*omega2 - 3.*omega2**2) + &
              2.*kd1**3*kd2**3*(4.*omega1**2 + 269.*omega1*omega2 + 4.*omega2**2) + 2.*kd1**2*kd2**4*(-259.*omega1**2 + 170.*omega1*omega2 + &
              9.*omega2**2))*swd**2) + 18014398509481984.*(384.*grav**2*kd1**3*kd2**3 - omega1*omega2*(6.*kd2**4*omega1**2 + 4.*kd1*kd2**3*omega1*(omega1 - 27.*omega2) + &
              6.*kd1**4*omega2**2 + 4.*kd1**3*kd2*omega2*(-27.*omega1 + omega2) + kd1**2*kd2**2*(37.*omega1**2 + 90.*omega1*omega2 + &
              37.*omega2**2))*swd**2) - 4398046511104.*(64.*grav**2*kd1**3*kd2**3*(27.*kd1**4 + 152.*kd1**3*kd2 + 134.*kd1**2*kd2**2 + 152.*kd1*kd2**3 + 27.*kd2**4) + &
              omega1*omega2*(20.*kd2**8*omega1**2 + 20.*kd1**8*omega2**2 + 20.*kd1**7*kd2*omega2*(11.*omega1 + 2.*omega2) + 20.*kd1*kd2**7*omega1*(2.*omega1 + 11.*omega2) + &
              kd1**6*kd2**2*(-315.*omega1**2 + 970.*omega1*omega2 - 1108.*omega2**2) + 2.*kd1**5*kd2**3*(89.*omega1**2 + 1658.*omega1*omega2 - 681.*omega2**2) + &
              kd1**2*kd2**6*(-1108.*omega1**2 + 970.*omega1*omega2 - 315.*omega2**2) + 2.*kd1**3*kd2**5*(-681.*omega1**2 + 1658.*omega1*omega2 + 89.*omega2**2) + &
              6.*kd1**4*kd2**4*(213.*omega1**2 + 2239.*omega1*omega2 + 213.*omega2**2))*swd**2) + 281474976710656.*(64.*grav**2*kd1**3*kd2**3*(11.*kd1**2 - 18.*kd1*kd2 + &
              11.*kd2**2) - omega1*omega2*(-15.*kd2**6*omega1**2 + 4.*kd1**5*kd2*(16.*omega1 - 5.*omega2)*omega2 - 15.*kd1**6*omega2**2 + &
              4.*kd1*kd2**5*omega1*(-5.*omega1 + 16.*omega2) + kd1**2*kd2**4*(417.*omega1**2 + 1328.*omega1*omega2 + 279.*omega2**2) + kd1**4*kd2**2*(279.*omega1**2 + &
              1328.*omega1*omega2 + 417.*omega2**2) + 2.*kd1**3*kd2**3*(533.*omega1**2 + 34.*omega1*omega2 + 533.*omega2**2))*swd**2) - &
              4096.*kd1**3*kd2**3*(kd1 + kd2)**2*(64.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(10.*kd1**4 + 28.*kd1**3*kd2 + 19.*kd1**2*kd2**2 + 28.*kd1*kd2**3 + &
              10.*kd2**4) + omega1*omega2*(4.*kd1**10*(13.*omega1 - 2.*omega2)*omega2 + 4.*kd2**10*omega1*(-2.*omega1 + 13.*omega2) + 2.*kd1**2*kd2**8*(315.*omega1**2 + &
              1682.*omega1*omega2 - 6.*omega2**2) + kd1*kd2**9*(32.*omega1**2 + 494.*omega1*omega2 + 3.*omega2**2) + kd1**9*kd2*(3.*omega1**2 + 494.*omega1*omega2 + &
              32.*omega2**2) + kd1**3*kd2**7*(1432.*omega1**2 + 11636.*omega1*omega2 + 33.*omega2**2) + 2.*kd1**8*kd2**2*(-6.*omega1**2 + 1682.*omega1*omega2 + &
              315.*omega2**2) + 4.*kd1**4*kd2**6*(713.*omega1**2 + 6066.*omega1*omega2 + 429.*omega2**2) + 4.*kd1**6*kd2**4*(429.*omega1**2 + 6066.*omega1*omega2 + &
              713.*omega2**2) + kd1**7*kd2**3*(33.*omega1**2 + 11636.*omega1*omega2 + 1432.*omega2**2) + kd1**5*kd2**5*(3957.*omega1**2 + 31658.*omega1*omega2 + &
              3957.*omega2**2))*swd**2) - 262144.*kd1*kd2*(64.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(3.*kd1**8 + 20.*kd1**7*kd2 + 24.*kd1**6*kd2**2 + 88.*kd1**5*kd2**3 + &
              213.*kd1**4*kd2**4 + 88.*kd1**3*kd2**5 + 24.*kd1**2*kd2**6 + 20.*kd1*kd2**7 + 3.*kd2**8) - omega1*omega2*(-4.*kd1**14*omega1*omega2 - &
              4.*kd2**14*omega1*omega2 - 2.*kd1**4*kd2**10*(871.*omega1**2 + 15542.*omega1*omega2 - 985.*omega2**2) + 2.*kd1**10*kd2**4*(985.*omega1**2 - &
              15542.*omega1*omega2 - 871.*omega2**2) - 2.*kd1**12*kd2**2*(omega1**2 + 462.*omega1*omega2 - 325.*omega2**2) + 2.*kd1**2*kd2**12*(325.*omega1**2 - &
              462.*omega1*omega2 - omega2**2) + kd1*kd2**13*(94.*omega1**2 - 78.*omega1*omega2 + omega2**2) + 6.*kd1**3*kd2**11*(39.*omega1**2 - 1347.*omega1*omega2 + &
              13.*omega2**2) + 6.*kd1**11*kd2**3*(13.*omega1**2 - 1347.*omega1*omega2 + 39.*omega2**2) + kd1**13*kd2*(omega1**2 - 78.*omega1*omega2 + 94.*omega2**2) + &
              2.*kd1**9*kd2**5*(3577.*omega1**2 - 30645.*omega1*omega2 + 124.*omega2**2) + 2.*kd1**8*kd2**6*(5163.*omega1**2 - 38446.*omega1*omega2 + 2783.*omega2**2) + &
              2.*kd1**5*kd2**9*(124.*omega1**2 - 30645.*omega1*omega2 + 3577.*omega2**2) + 2.*kd1**7*kd2**7*(4451.*omega1**2 - 39511.*omega1*omega2 + 4451.*omega2**2) + &
              2.*kd1**6*kd2**8*(2783.*omega1**2 - 38446.*omega1*omega2 + 5163.*omega2**2))*swd**2) + 1073741824.*(128.*grav**2*kd1**3*kd2**3*(12.*kd1**8 + &
              116.*kd1**7*kd2 + 316.*kd1**6*kd2**2 + 410.*kd1**5*kd2**3 + 407.*kd1**4*kd2**4 + 410.*kd1**3*kd2**5 + 316.*kd1**2*kd2**6 + 116.*kd1*kd2**7 + 12.*kd2**8) - &
              omega1*omega2*(6.*kd2**12*omega1**2 + 4.*kd1*kd2**11*omega1*(5.*omega1 - 29.*omega2) + 6.*kd1**12*omega2**2 + 4.*kd1**11*kd2*omega2*(-29.*omega1 + 5.*omega2) + &
              kd1**2*kd2**10*(285.*omega1**2 - 242.*omega1*omega2 + 279.*omega2**2) + kd1**10*kd2**2*(279.*omega1**2 - 242.*omega1*omega2 + 285.*omega2**2) + &
              4.*kd1**3*kd2**9*(607.*omega1**2 + 103.*omega1*omega2 + 298.*omega2**2) + 4.*kd1**9*kd2**3*(298.*omega1**2 + 103.*omega1*omega2 + 607.*omega2**2) + &
              2.*kd1**4*kd2**8*(2746.*omega1**2 - 8662.*omega1*omega2 + 965.*omega2**2) + 4.*kd1**5*kd2**7*(2872.*omega1**2 - 16027.*omega1*omega2 + 2057.*omega2**2) + &
              2.*kd1**8*kd2**4*(965.*omega1**2 - 8662.*omega1*omega2 + 2746.*omega2**2) + 4.*kd1**7*kd2**5*(2057.*omega1**2 - 16027.*omega1*omega2 + 2872.*omega2**2) + &
              3.*kd1**6*kd2**6*(5411.*omega1**2 - 29962.*omega1*omega2 + 5411.*omega2**2))*swd**2) - 68719476736.*(128.*grav**2*kd1**3*kd2**3*(5.*kd1**6 - 10.*kd1**5*kd2 - &
              72.*kd1**4*kd2**2 - 78.*kd1**3*kd2**3 - 72.*kd1**2*kd2**4 - 10.*kd1*kd2**5 + 5.*kd2**6) - omega1*omega2*(15.*kd2**10*omega1**2 + 15.*kd1**10*omega2**2 + &
              8.*kd1**9*kd2*omega2*(13.*omega1 + 5.*omega2) + 8.*kd1*kd2**9*omega1*(5.*omega1 + 13.*omega2) + kd1**8*kd2**2*(315.*omega1**2 + 3280.*omega1*omega2 - &
              546.*omega2**2) + kd1**2*kd2**8*(-546.*omega1**2 + 3280.*omega1*omega2 + 315.*omega2**2) + 4.*kd1**7*kd2**3*(783.*omega1**2 + 2445.*omega1*omega2 + &
              409.*omega2**2) + 4.*kd1**3*kd2**7*(409.*omega1**2 + 2445.*omega1*omega2 + 783.*omega2**2) + 2.*kd1**5*kd2**5*(2215.*omega1**2 - 1362.*omega1*omega2 + &
              2215.*omega2**2) + kd1**6*kd2**4*(5431.*omega1**2 + 3682.*omega1*omega2 + 5322.*omega2**2) + kd1**4*kd2**6*(5322.*omega1**2 + 3682.*omega1*omega2 + &
              5431.*omega2**2))*swd**2) - 16777216.*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**2*(kd1**8 + 48.*kd1**7*kd2 + 168.*kd1**6*kd2**2 + 36.*kd1**5*kd2**3 - &
              217.*kd1**4*kd2**4 + 36.*kd1**3*kd2**5 + 168.*kd1**2*kd2**6 + 48.*kd1*kd2**7 + kd2**8) - omega1*omega2*(kd2**14*omega1**2 + 4.*kd1*kd2**13*omega1*(omega1 - &
              12.*omega2) + kd1**14*omega2**2 + 4.*kd1**13*kd2*omega2*(-12.*omega1 + omega2) + kd1**4*kd2**10*(1509.*omega1**2 + 388.*omega1*omega2 - 1875.*omega2**2) + &
              kd1**6*kd2**8*(12705.*omega1**2 - 19896.*omega1*omega2 - 173.*omega2**2) + 8.*kd1**5*kd2**9*(1210.*omega1**2 + 555.*omega1*omega2 - 137.*omega2**2) + &
              kd1**2*kd2**12*(83.*omega1**2 - 1056.*omega1*omega2 - 37.*omega2**2) - 2.*kd1**11*kd2**3*(345.*omega1**2 + 1534.*omega1*omega2 + 81.*omega2**2) + &
              kd1**12*kd2**2*(-37.*omega1**2 - 1056.*omega1*omega2 + 83.*omega2**2) + 22.*kd1**7*kd2**7*(189.*omega1**2 - 2014.*omega1*omega2 + 189.*omega2**2) - &
              2.*kd1**3*kd2**11*(81.*omega1**2 + 1534.*omega1*omega2 + 345.*omega2**2) + 8.*kd1**9*kd2**5*(-137.*omega1**2 + 555.*omega1*omega2 + 1210.*omega2**2) + &
              kd1**10*kd2**4*(-1875.*omega1**2 + 388.*omega1*omega2 + 1509.*omega2**2) + kd1**8*kd2**6*(-173.*omega1**2 - 19896.*omega1*omega2 + 12705.*omega2**2))*swd**2)))/ &
              (2.*(-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*omega2*swd**2*(-64.*grav*(kd1 + kd2)**2*(64. + (kd1 + kd2)**2)*(4096. + &
              (kd1 + kd2)**2*(384. + (kd1 + kd2)**2)) + (16777216. + (kd1 + kd2)**2*(7340032. + (kd1 + kd2)**2*(286720. + (kd1 + kd2)**2*(1792. + (kd1 + kd2)**2))))*(omega1 + omega2)**2*swd))
    !
end function etasup4
!
real function velsb14()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic velocity of 1st layer for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb14')
    !
    velsb14 = (-4.*(-144115188075855872.*grav*(kd1 - kd2)**3*omega1*omega2*(kd2**2*omega1**2 + 8.*kd1*kd2*omega1*omega2 + kd1**2*omega2**2)*swd**2 - &
              16.*grav*kd1**7*(kd1 - kd2)**7*kd2**7*((9.*kd1**2 - 29.*kd1*kd2 + 16.*kd2**2)*omega1**4 - (15.*kd1**2 + 16.*kd1*kd2 + 45.*kd2**2)*omega1**3*omega2 + &
              (37.*kd1**2 + 86.*kd1*kd2 + 37.*kd2**2)*omega1**2*omega2**2 - (45.*kd1**2 + 16.*kd1*kd2 + 15.*kd2**2)*omega1*omega2**3 + (16.*kd1**2 - 29.*kd1*kd2 + &
              9.*kd2**2)*omega2**4)*swd**2 + kd1**7*(kd1 - kd2)**7*kd2**7*omega1*(omega1 - omega2)**2*omega2*(kd2**2*omega1*(3.*omega1 + omega2) + &
              kd1**2*omega2*(omega1 + 3.*omega2) + 2.*kd1*kd2*(omega1**2 - 6.*omega1*omega2 + omega2**2))*swd**3 + &
              18014398509481984.*swd*(8.*grav**2*kd1*(kd1 - kd2)*kd2*(6.*kd2**4*omega1**2 + kd1**3*kd2*(omega1 - 12.*omega2)*omega2 + 6.*kd1**4*omega2**2 + &
              kd1*kd2**3*omega1*(-12.*omega1 + omega2) + 5.*kd1**2*kd2**2*(omega1**2 + omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(-(kd2**2*(-2.*kd1**3 + &
              11.*kd1**2*kd2 - 3.*kd1*kd2**2 + kd2**3)*omega1**2) + kd1*(kd1 - kd2)*kd2*(8.*kd1**2 - 19.*kd1*kd2 + 8.*kd2**2)*omega1*omega2 + kd1**2*(kd1**3 - &
              3.*kd1**2*kd2 + 11.*kd1*kd2**2 - 2.*kd2**3)*omega2**2)*swd**2) - 4398046511104.*swd*(8.*grav**2*kd1*kd2*(-kd1 + kd2)*(-10.*kd2**8*omega1**2 + &
              2.*kd1*kd2**7*omega1*(109.*omega1 - 7.*omega2) - 10.*kd1**8*omega2**2 + 2.*kd1**7*kd2*omega2*(-7.*omega1 + 109.*omega2) + kd1**6*kd2**2*(183.*omega1**2 - &
              113.*omega1*omega2 - 376.*omega2**2) + 2.*kd1**3*kd2**5*(369.*omega1**2 - 62.*omega1*omega2 + 129.*omega2**2) + kd1**2*kd2**6*(-376.*omega1**2 - &
              113.*omega1*omega2 + 183.*omega2**2) - 5.*kd1**4*kd2**4*(203.*omega1**2 - 102.*omega1*omega2 + 203.*omega2**2) + 2.*kd1**5*kd2**3*(129.*omega1**2 - &
              62.*omega1*omega2 + 369.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(kd2**2*(246.*kd1**7 - 1081.*kd1**6*kd2 + 1107.*kd1**5*kd2**2 - &
              301.*kd1**4*kd2**3 + 415.*kd1**3*kd2**4 - 288.*kd1**2*kd2**5 - 38.*kd1*kd2**6 + 10.*kd2**7)*omega1**2 + 2.*kd1*kd2*(-kd1 + kd2)*(26.*kd1**6 + &
              334.*kd1**5*kd2 + 43.*kd1**4*kd2**2 - 861.*kd1**3*kd2**3 + 43.*kd1**2*kd2**4 + 334.*kd1*kd2**5 + 26.*kd2**6)*omega1*omega2 + kd1**2*(-10.*kd1**7 + &
              38.*kd1**6*kd2 + 288.*kd1**5*kd2**2 - 415.*kd1**4*kd2**3 + 301.*kd1**3*kd2**4 - 1107.*kd1**2*kd2**5 + 1081.*kd1*kd2**6 - &
              246.*kd2**7)*omega2**2)*swd**2) - 2251799813685248.*grav*(320.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3 + (2.*kd1*kd2**3*(28.*kd1**3 - 85.*kd1**2*kd2 + &
              78.*kd1*kd2**2 - 24.*kd2**3)*omega1**4 + kd2**2*(-18.*kd1**5 - 131.*kd1**4*kd2 + 406.*kd1**3*kd2**2 - 329.*kd1**2*kd2**3 + 79.*kd1*kd2**4 + &
              5.*kd2**5)*omega1**3*omega2 + 2.*kd1*(kd1 - kd2)*kd2*(kd1**2 + kd2**2)*(82.*kd1**2 - 155.*kd1*kd2 + 82.*kd2**2)*omega1**2*omega2**2 + kd1**2*(-5.*kd1**5 - &
              79.*kd1**4*kd2 + 329.*kd1**3*kd2**2 - 406.*kd1**2*kd2**3 + 131.*kd1*kd2**4 + 18.*kd2**5)*omega1*omega2**3 + 2.*kd1**3*kd2*(24.*kd1**3 - &
              78.*kd1**2*kd2 + 85.*kd1*kd2**2 - 28.*kd2**3)*omega2**4)*swd**2) - 35184372088832.*grav*(384.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(2.*kd1**2 + &
              3.*kd1*kd2 + 2.*kd2**2) + (2.*kd1*kd2**3*(360.*kd1**5 - 1134.*kd1**4*kd2 + 1205.*kd1**3*kd2**2 - 494.*kd1**2*kd2**3 + 143.*kd1*kd2**4 - &
              68.*kd2**5)*omega1**4 + kd2**2*(-405.*kd1**7 + 232.*kd1**6*kd2 + 2011.*kd1**5*kd2**2 - 2817.*kd1**4*kd2**3 + 969.*kd1**3*kd2**4 - 338.*kd1**2*kd2**5 + &
              310.*kd1*kd2**6 - 10.*kd2**7)*omega1**3*omega2 + 2.*kd1*(kd1 - kd2)*kd2*(94.*kd1**6 + 492.*kd1**5*kd2 + 38.*kd1**4*kd2**2 - 1333.*kd1**3*kd2**3 + &
              38.*kd1**2*kd2**4 + 492.*kd1*kd2**5 + 94.*kd2**6)*omega1**2*omega2**2 + kd1**2*(10.*kd1**7 - 310.*kd1**6*kd2 + 338.*kd1**5*kd2**2 - 969.*kd1**4*kd2**3 + &
              2817.*kd1**3*kd2**4 - 2011.*kd1**2*kd2**5 - 232.*kd1*kd2**6 + 405.*kd2**7)*omega1*omega2**3 + 2.*kd1**3*kd2*(68.*kd1**5 - 143.*kd1**4*kd2 + &
              494.*kd1**3*kd2**2 - 1205.*kd1**2*kd2**3 + 1134.*kd1*kd2**4 - 360.*kd2**5)*omega2**4)*swd**2) - &
              512.*grav*kd1**5*(kd1 - kd2)**5*kd2**5*(128.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2 + (2.*(kd1 - kd2)*(43.*kd1**5 - 177.*kd1**4*kd2 + 391.*kd1**3*kd2**2 - &
              878.*kd1**2*kd2**3 + 781.*kd1*kd2**4 - 236.*kd2**5)*omega1**4 + (-178.*kd1**6 + 721.*kd1**5*kd2 - 2391.*kd1**4*kd2**2 + 5222.*kd1**3*kd2**3 - &
              6491.*kd1**2*kd2**4 + 3406.*kd1*kd2**5 - 838.*kd2**6)*omega1**3*omega2 + 2.*(259.*kd1**6 - 1086.*kd1**5*kd2 + 2540.*kd1**4*kd2**2 - 2877.*kd1**3*kd2**3 + &
              2540.*kd1**2*kd2**4 - 1086.*kd1*kd2**5 + 259.*kd2**6)*omega1**2*omega2**2 + (-838.*kd1**6 + 3406.*kd1**5*kd2 - 6491.*kd1**4*kd2**2 + 5222.*kd1**3*kd2**3 - &
              2391.*kd1**2*kd2**4 + 721.*kd1*kd2**5 - 178.*kd2**6)*omega1*omega2**3 + 2.*(kd1 - kd2)*(236.*kd1**5 - 781.*kd1**4*kd2 + 878.*kd1**3*kd2**2 - &
              391.*kd1**2*kd2**3 + 177.*kd1*kd2**4 - 43.*kd2**5)*omega2**4)*swd**2) + 1099511627776.*grav*(64.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(2.*kd1**4 - &
              55.*kd1**3*kd2 + 50.*kd1**2*kd2**2 - 55.*kd1*kd2**3 + 2.*kd2**4) + (-(kd1*kd2**3*(148.*kd1**7 + 1113.*kd1**6*kd2 - 3702.*kd1**5*kd2**2 + &
              4145.*kd1**4*kd2**3 - 3258.*kd1**3*kd2**4 + 2424.*kd1**2*kd2**5 - 928.*kd1*kd2**6 + 40.*kd2**7)*omega1**4) + kd2**2*(38.*kd1**9 - 749.*kd1**8*kd2 + &
              5435.*kd1**7*kd2**2 - 10814.*kd1**6*kd2**3 + 9430.*kd1**5*kd2**4 - 6286.*kd1**4*kd2**5 + 4621.*kd1**3*kd2**6 - 1807.*kd1**2*kd2**7 + 101.*kd1*kd2**8 - &
              5.*kd2**9)*omega1**3*omega2 + kd1*(kd1 - kd2)*kd2*(124.*kd1**8 - 1950.*kd1**7*kd2 + 5108.*kd1**6*kd2**2 - 11657.*kd1**5*kd2**3 + 16642.*kd1**4*kd2**4 - &
              11657.*kd1**3*kd2**5 + 5108.*kd1**2*kd2**6 - 1950.*kd1*kd2**7 + 124.*kd2**8)*omega1**2*omega2**2 + kd1**2*(5.*kd1**9 - 101.*kd1**8*kd2 + &
              1807.*kd1**7*kd2**2 - 4621.*kd1**6*kd2**3 + 6286.*kd1**5*kd2**4 - 9430.*kd1**4*kd2**5 + 10814.*kd1**3*kd2**6 - 5435.*kd1**2*kd2**7 + 749.*kd1*kd2**8 - &
              38.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(40.*kd1**7 - 928.*kd1**6*kd2 + 2424.*kd1**5*kd2**2 - 3258.*kd1**4*kd2**3 + 4145.*kd1**3*kd2**4 - &
              3702.*kd1**2*kd2**5 + 1113.*kd1*kd2**6 + 148.*kd2**7)*omega2**4)*swd**2) + 8589934592.*grav*(128.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(7.*kd1**6 - &
              37.*kd1**5*kd2 + 51.*kd1**4*kd2**2 - 81.*kd1**3*kd2**3 + 51.*kd1**2*kd2**4 - 37.*kd1*kd2**5 + 7.*kd2**6) + (2.*kd1*kd2**3*(336.*kd1**9 - 5388.*kd1**8*kd2 + &
              20023.*kd1**7*kd2**2 - 36685.*kd1**6*kd2**3 + 40472.*kd1**5*kd2**4 - 27645.*kd1**4*kd2**5 + 10899.*kd1**3*kd2**6 - 2574.*kd1**2*kd2**7 + 630.*kd1*kd2**8 - &
              80.*kd2**9)*omega1**4 + kd2**2*(-391.*kd1**11 + 739.*kd1**10*kd2 + 16127.*kd1**9*kd2**2 - 65207.*kd1**8*kd2**3 + 115523.*kd1**7*kd2**4 - &
              127615.*kd1**6*kd2**5 + 94299.*kd1**5*kd2**6 - 42507.*kd1**4*kd2**7 + 11427.*kd1**3*kd2**8 - 2649.*kd1**2*kd2**9 + 297.*kd1*kd2**10 + 5.*kd2**11)*omega1**3*omega2 + &
              2.*kd1*(kd1 - kd2)*kd2*(112.*kd1**10 - 52.*kd1**9*kd2 + 726.*kd1**8*kd2**2 - 17491.*kd1**7*kd2**3 + 49247.*kd1**6*kd2**4 - 65058.*kd1**5*kd2**5 + &
              49247.*kd1**4*kd2**6 - 17491.*kd1**3*kd2**7 + 726.*kd1**2*kd2**8 - 52.*kd1*kd2**9 + 112.*kd2**10)*omega1**2*omega2**2 - kd1**2*(5.*kd1**11 + 297.*kd1**10*kd2 - &
              2649.*kd1**9*kd2**2 + 11427.*kd1**8*kd2**3 - 42507.*kd1**7*kd2**4 + 94299.*kd1**6*kd2**5 - 127615.*kd1**5*kd2**6 + 115523.*kd1**4*kd2**7 - 65207.*kd1**3*kd2**8 + &
              16127.*kd1**2*kd2**9 + 739.*kd1*kd2**10 - 391.*kd2**11)*omega1*omega2**3 + 2.*kd1**3*kd2*(80.*kd1**9 - 630.*kd1**8*kd2 + 2574.*kd1**7*kd2**2 - 10899.*kd1**6*kd2**3 + &
              27645.*kd1**5*kd2**4 - 40472.*kd1**4*kd2**5 + 36685.*kd1**3*kd2**6 - 20023.*kd1**2*kd2**7 + 5388.*kd1*kd2**8 - 336.*kd2**9)*omega2**4)*swd**2) + &
              64.*kd1**5*(kd1 - kd2)**5*kd2**5*swd*(8.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd2**2*(20.*omega1 - 27.*omega2)*(omega1 - omega2) + kd1**2*(27.*omega1 - &
              20.*omega2)*(omega1 - omega2) - 2.*kd1*kd2*(24.*omega1**2 - 41.*omega1*omega2 + 24.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(60.*kd1**6*omega1*omega2 + &
              60.*kd2**6*omega1*omega2 + kd1**2*kd2**4*(omega1 + omega2)*(407.*omega1 + 3.*omega2) + kd1**4*kd2**2*(omega1 + omega2)*(3.*omega1 + 407.*omega2) + &
              4.*kd1**3*kd2**3*(omega1**2 - 130.*omega1*omega2 + omega2**2) - 6.*kd1*kd2**5*(58.*omega1**2 + 45.*omega1*omega2 + omega2**2) - 6.*kd1**5*kd2*(omega1**2 + &
              45.*omega1*omega2 + 58.*omega2**2))*swd**2) + 281474976710656.*swd*(8.*grav**2*kd1*(kd1 - kd2)*kd2*(17.*kd2**6*omega1**2 + kd1*kd2**5*omega1*(10.*omega1 - &
              23.*omega2) + 17.*kd1**6*omega2**2 + kd1**5*kd2*omega2*(-23.*omega1 + 10.*omega2) + kd1**4*kd2**2*(94.*omega1**2 + 41.*omega1*omega2 + 29.*omega2**2) - &
              4.*kd1**3*kd2**3*(37.*omega1**2 + 10.*omega1*omega2 + 37.*omega2**2) + kd1**2*kd2**4*(29.*omega1**2 + 41.*omega1*omega2 + 94.*omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(5.*kd2**7*omega1**2 - 5.*kd1**7*omega2**2 + kd1**6*kd2*omega2*(116.*omega1 + 17.*omega2) - &
              kd1*kd2**6*omega1*(17.*omega1 + 116.*omega2) + kd1**2*kd2**5*(-151.*omega1**2 + 365.*omega1*omega2 + 8.*omega2**2) + 2.*kd1**3*kd2**4*(191.*omega1**2 - &
              271.*omega1*omega2 + 88.*omega2**2) + kd1**5*kd2**2*(-8.*omega1**2 - 365.*omega1*omega2 + 151.*omega2**2) - 2.*kd1**4*kd2**3*(88.*omega1**2 - &
              271.*omega1*omega2 + 191.*omega2**2))*swd**2) + 4096.*kd1**3*(kd1 - kd2)**3*kd2**3*swd*(8.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*((kd1 - kd2)*(84.*kd1**5 - &
              356.*kd1**4*kd2 + 520.*kd1**3*kd2**2 - 584.*kd1**2*kd2**3 + 521.*kd1*kd2**4 - 155.*kd2**5)*omega1**2 + (-156.*kd1**6 + 729.*kd1**5*kd2 - 1359.*kd1**4*kd2**2 + &
              1544.*kd1**3*kd2**3 - 1359.*kd1**2*kd2**4 + 729.*kd1*kd2**5 - 156.*kd2**6)*omega1*omega2 + (kd1 - kd2)*(155.*kd1**5 - 521.*kd1**4*kd2 + 584.*kd1**3*kd2**2 - &
              520.*kd1**2*kd2**3 + 356.*kd1*kd2**4 - 84.*kd2**5)*omega2**2) + omega1*(omega1 - omega2)**2*omega2*(-3.*kd2**10*omega1*(omega1 - 21.*omega2) + &
              3.*kd1**10*(21.*omega1 - omega2)*omega2 + kd1**8*kd2**2*(-15.*omega1**2 + 3193.*omega1*omega2 - 250.*omega2**2) + kd1**2*kd2**8*(-250.*omega1**2 + &
              3193.*omega1*omega2 - 15.*omega2**2) + kd1*kd2**9*(61.*omega1**2 - 556.*omega1*omega2 + 6.*omega2**2) + kd1**3*kd2**7*(526.*omega1**2 - 11167.*omega1*omega2 + &
              24.*omega2**2) + kd1**9*kd2*(6.*omega1**2 - 556.*omega1*omega2 + 61.*omega2**2) + 8.*kd1**5*kd2**5*(411.*omega1**2 - 3841.*omega1*omega2 + 411.*omega2**2) + &
              kd1**7*kd2**3*(24.*omega1**2 - 11167.*omega1*omega2 + 526.*omega2**2) - 2.*kd1**4*kd2**6*(1156.*omega1**2 - 11886.*omega1*omega2 + 633.*omega2**2) - &
              2.*kd1**6*kd2**4*(633.*omega1**2 - 11886.*omega1*omega2 + 1156.*omega2**2))*swd**2) + &
              262144.*kd1*(kd1 - kd2)*kd2*swd*(8.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd1**7*kd2**3*(-3332.*omega1**2 + 3735.*omega1*omega2 - 4660.*omega2**2) + &
              kd1**3*kd2**7*(-4660.*omega1**2 + 3735.*omega1*omega2 - 3332.*omega2**2) + kd2**10*(78.*omega1**2 - 53.*omega1*omega2 + 33.*omega2**2) - &
              4.*kd1*kd2**9*(163.*omega1**2 - 109.*omega1*omega2 + 65.*omega2**2) + kd1**10*(33.*omega1**2 - 53.*omega1*omega2 + 78.*omega2**2) - &
              4.*kd1**9*kd2*(65.*omega1**2 - 109.*omega1*omega2 + 163.*omega2**2) + kd1**2*kd2**8*(2101.*omega1**2 - 1333.*omega1*omega2 + 972.*omega2**2) + &
              kd1**8*kd2**2*(972.*omega1**2 - 1333.*omega1*omega2 + 2101.*omega2**2) - 2.*kd1**5*kd2**5*(5595.*omega1**2 - 6406.*omega1*omega2 + 5595.*omega2**2) + &
              kd1**4*kd2**6*(8761.*omega1**2 - 9195.*omega1*omega2 + 8150.*omega2**2) + kd1**6*kd2**4*(8150.*omega1**2 - 9195.*omega1*omega2 + 8761.*omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(4.*kd1**14*omega1*omega2 + 4.*kd2**14*omega1*omega2 + kd1**5*kd2**9*(2803.*omega1**2 - 80713.*omega1*omega2 - &
              5412.*omega2**2) + kd1**10*kd2**4*(1339.*omega1**2 + 37970.*omega1*omega2 - 2551.*omega2**2) + 2.*kd1**8*kd2**6*(3603.*omega1**2 + 56359.*omega1*omega2 - &
              643.*omega2**2) + kd1**3*kd2**11*(571.*omega1**2 - 9908.*omega1*omega2 - 36.*omega2**2) - 2.*kd1*kd2**13*(43.*omega1**2 + 35.*omega1*omega2 + omega2**2) + &
              kd1**2*kd2**12*(331.*omega1**2 + 1194.*omega1*omega2 + 9.*omega2**2) - 2.*kd1**13*kd2*(omega1**2 + 35.*omega1*omega2 + 43.*omega2**2) + &
              kd1**12*kd2**2*(9.*omega1**2 + 1194.*omega1*omega2 + 331.*omega2**2) + kd1**11*kd2**3*(-36.*omega1**2 - 9908.*omega1*omega2 + 571.*omega2**2) + &
              kd1**4*kd2**10*(-2551.*omega1**2 + 37970.*omega1*omega2 + 1339.*omega2**2) - 2.*kd1**7*kd2**7*(1441.*omega1**2 + 61199.*omega1*omega2 + 1441.*omega2**2) + &
              kd1**9*kd2**5*(-5412.*omega1**2 - 80713.*omega1*omega2 + 2803.*omega2**2) + 2.*kd1**6*kd2**8*(-643.*omega1**2 + 56359.*omega1*omega2 + &
              3603.*omega2**2))*swd**2) - 68719476736.*swd*(16.*grav**2*kd1*(kd1 - kd2)*kd2*(10.*kd2**10*omega1**2 + 10.*kd1**10*omega2**2 - &
              13.*kd1*kd2**9*omega1*(3.*omega1 + omega2) - 13.*kd1**9*kd2*omega2*(omega1 + 3.*omega2) + kd1**5*kd2**5*(-3743.*omega1**2 + 1390.*omega1*omega2 - &
              3743.*omega2**2) + kd1**7*kd2**3*(-808.*omega1**2 + 400.*omega1*omega2 - 1163.*omega2**2) + kd1**3*kd2**7*(-1163.*omega1**2 + 400.*omega1*omega2 - &
              808.*omega2**2) + kd1**2*kd2**8*(80.*omega1**2 + 69.*omega1*omega2 + 22.*omega2**2) + kd1**8*kd2**2*(22.*omega1**2 + 69.*omega1*omega2 + 80.*omega2**2) + &
              kd1**4*kd2**6*(3131.*omega1**2 - 1148.*omega1*omega2 + 2507.*omega2**2) + kd1**6*kd2**4*(2507.*omega1**2 - 1148.*omega1*omega2 + 3131.*omega2**2)) - &
              omega1*(omega1 - omega2)**2*omega2*(10.*kd2**11*omega1**2 + 6.*kd1**2*kd2**9*omega1*(11.*omega1 - 365.*omega2) - 42.*kd1*kd2**10*omega1*(omega1 - 4.*omega2) + &
              6.*kd1**9*kd2**2*(365.*omega1 - 11.*omega2)*omega2 - 10.*kd1**11*omega2**2 + 42.*kd1**10*kd2*omega2*(-4.*omega1 + omega2) + kd1**5*kd2**6*(-1079.*omega1**2 + &
              28106.*omega1*omega2 - 4409.*omega2**2) + kd1**3*kd2**8*(-659.*omega1**2 + 7565.*omega1*omega2 - 1304.*omega2**2) + kd1**7*kd2**4*(-4542.*omega1**2 + &
              17429.*omega1*omega2 - 603.*omega2**2) + kd1**8*kd2**3*(1304.*omega1**2 - 7565.*omega1*omega2 + 659.*omega2**2) + kd1**6*kd2**5*(4409.*omega1**2 - &
              28106.*omega1*omega2 + 1079.*omega2**2) + kd1**4*kd2**7*(603.*omega1**2 - 17429.*omega1*omega2 + 4542.*omega2**2))*swd**2) + &
              1073741824.*swd*(8.*grav**2*kd1*(kd1 - kd2)*kd2*(4.*kd2**12*omega1**2 + kd1**11*kd2*(13.*omega1 - 134.*omega2)*omega2 + 4.*kd1**12*omega2**2 + &
              kd1*kd2**11*omega1*(-134.*omega1 + 13.*omega2) + kd1**6*kd2**6*(-8413.*omega1**2 + 13248.*omega1*omega2 - 8413.*omega2**2) + kd1**2*kd2**10*(241.*omega1**2 + &
              286.*omega1*omega2 - 317.*omega2**2) + kd1**10*kd2**2*(-317.*omega1**2 + 286.*omega1*omega2 + 241.*omega2**2) - 2.*kd1**4*kd2**8*(3443.*omega1**2 - &
              3462.*omega1*omega2 + 1426.*omega2**2) + kd1**3*kd2**9*(1834.*omega1**2 - 2523.*omega1*omega2 + 1568.*omega2**2) + kd1**9*kd2**3*(1568.*omega1**2 - &
              2523.*omega1*omega2 + 1834.*omega2**2) - 2.*kd1**8*kd2**4*(1426.*omega1**2 - 3462.*omega1*omega2 + 3443.*omega2**2) + kd1**5*kd2**7*(10306.*omega1**2 - &
              11321.*omega1*omega2 + 4646.*omega2**2) + kd1**7*kd2**5*(4646.*omega1**2 - 11321.*omega1*omega2 + 10306.*omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(-5.*kd2**13*omega1**2 + 5.*kd1**13*omega2**2 - kd1**12*kd2*omega2*(64.*omega1 + 23.*omega2) + &
              kd1*kd2**12*omega1*(23.*omega1 + 64.*omega2) + kd1**10*kd2**3*(-1753.*omega1**2 + 1115.*omega1*omega2 - 614.*omega2**2) + kd1**2*kd2**11*(-251.*omega1**2 + &
              355.*omega1*omega2 - 246.*omega2**2) + kd1**11*kd2**2*(246.*omega1**2 - 355.*omega1*omega2 + 251.*omega2**2) + kd1**3*kd2**10*(614.*omega1**2 - &
              1115.*omega1*omega2 + 1753.*omega2**2) + kd1**9*kd2**4*(4301.*omega1**2 + 14902.*omega1*omega2 + 1785.*omega2**2) - 5.*kd1**6*kd2**7*(3859.*omega1**2 + &
              25027.*omega1*omega2 + 3372.*omega2**2) + 5.*kd1**7*kd2**6*(3372.*omega1**2 + 25027.*omega1*omega2 + 3859.*omega2**2) - kd1**4*kd2**9*(1785.*omega1**2 + &
              14902.*omega1*omega2 + 4301.*omega2**2) + kd1**5*kd2**8*(9517.*omega1**2 + 67890.*omega1*omega2 + 8507.*omega2**2) - kd1**8*kd2**5*(8507.*omega1**2 + &
              67890.*omega1*omega2 + 9517.*omega2**2))*swd**2) - 16777216.*swd*(8.*grav**2*kd1*(kd1 - kd2)**2*kd2*(3.*kd2**13*omega1**2 - 3.*kd1**13*omega2**2 - &
              kd1*kd2**12*omega1*(29.*omega1 + 3.*omega2) + kd1**12*kd2*omega2*(3.*omega1 + 29.*omega2) + kd1**5*kd2**8*(-8853.*omega1**2 + 3811.*omega1*omega2 - &
              4545.*omega2**2) + 2.*kd1**2*kd2**11*(45.*omega1**2 - 13.*omega1*omega2 + 23.*omega2**2) - 2.*kd1**11*kd2**2*(23.*omega1**2 - 13.*omega1*omega2 + 45.*omega2**2) - &
              2.*kd1**3*kd2**10*(654.*omega1**2 - 483.*omega1*omega2 + 431.*omega2**2) + 2.*kd1**10*kd2**3*(431.*omega1**2 - 483.*omega1*omega2 + 654.*omega2**2) + &
              3.*kd1**4*kd2**9*(1873.*omega1**2 - 1213.*omega1*omega2 + 1101.*omega2**2) - 3.*kd1**9*kd2**4*(1101.*omega1**2 - 1213.*omega1*omega2 + 1873.*omega2**2) + &
              2.*kd1**6*kd2**7*(3125.*omega1**2 + 42.*omega1*omega2 + 1912.*omega2**2) - 2.*kd1**7*kd2**6*(1912.*omega1**2 + 42.*omega1*omega2 + 3125.*omega2**2) + &
              kd1**8*kd2**5*(4545.*omega1**2 - 3811.*omega1*omega2 + 8853.*omega2**2)) - omega1*(omega1 - omega2)**2*omega2*(kd2**15*omega1**2 - kd1**15*omega2**2 + &
              kd1**14*kd2*omega2*(52.*omega1 + 5.*omega2) - kd1*kd2**14*omega1*(5.*omega1 + 52.*omega2) + kd1**4*kd2**11*(1367.*omega1**2 + 5454.*omega1*omega2 - &
              1770.*omega2**2) + kd1**11*kd2**4*(1770.*omega1**2 - 5454.*omega1*omega2 - 1367.*omega2**2) + kd1**2*kd2**13*(-27.*omega1**2 + 961.*omega1*omega2 - 8.*omega2**2) + &
              kd1**13*kd2**2*(8.*omega1**2 - 961.*omega1*omega2 + 27.*omega2**2) + 2.*kd1**3*kd2**12*(245.*omega1**2 - 1800.*omega1*omega2 + 228.*omega2**2) - &
              2.*kd1**12*kd2**3*(228.*omega1**2 - 1800.*omega1*omega2 + 245.*omega2**2) + kd1**5*kd2**10*(-10887.*omega1**2 - 12931.*omega1*omega2 + 721.*omega2**2) + &
              kd1**6*kd2**9*(20273.*omega1**2 + 48969.*omega1*omega2 + 5223.*omega2**2) + kd1**10*kd2**5*(-721.*omega1**2 + 12931.*omega1*omega2 + 10887.*omega2**2) - &
              kd1**7*kd2**8*(18712.*omega1**2 + 97964.*omega1*omega2 + 12129.*omega2**2) + kd1**8*kd2**7*(12129.*omega1**2 + 97964.*omega1*omega2 + 18712.*omega2**2) - &
              kd1**9*kd2**6*(5223.*omega1**2 + 48969.*omega1*omega2 + 20273.*omega2**2))*swd**2) - &
              32768.*grav*kd1**3*(kd1 - kd2)**3*kd2**3*(64.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(12.*kd1**4 - 32.*kd1**3*kd2 + 21.*kd1**2*kd2**2 - 32.*kd1*kd2**3 + &
              12.*kd2**4) + (kd1**5*kd2**5*(-41530.*omega1**4 + 92003.*omega1**3*omega2 - 127268.*omega1**2*omega2**2 + 92003.*omega1*omega2**3 - 41530.*omega2**4) + &
              kd1**9*kd2*(-282.*omega1**4 + 518.*omega1**3*omega2 - 4852.*omega1**2*omega2**2 + 7643.*omega1*omega2**3 - 3726.*omega2**4) + kd1*kd2**9*(-3726.*omega1**4 + &
              7643.*omega1**3*omega2 - 4852.*omega1**2*omega2**2 + 518.*omega1*omega2**3 - 282.*omega2**4) + 4.*kd2**10*(124.*omega1**4 - 268.*omega1**3*omega2 + &
              165.*omega1**2*omega2**2 - 12.*omega1*omega2**3 + 6.*omega2**4) + 4.*kd1**10*(6.*omega1**4 - 12.*omega1**3*omega2 + 165.*omega1**2*omega2**2 - &
              268.*omega1*omega2**3 + 124.*omega2**4) + 2.*kd1**2*kd2**8*(5583.*omega1**4 - 10995.*omega1**3*omega2 + 8498.*omega1**2*omega2**2 - 2542.*omega1*omega2**3 + &
              1347.*omega2**4) + 2.*kd1**8*kd2**2*(1347.*omega1**4 - 2542.*omega1**3*omega2 + 8498.*omega1**2*omega2**2 - 10995.*omega1*omega2**3 + 5583.*omega2**4) - &
              2.*kd1**3*kd2**7*(11124.*omega1**4 - 22636.*omega1**3*omega2 + 23816.*omega1**2*omega2**2 - 12594.*omega1*omega2**3 + 6183.*omega2**4) - &
              2.*kd1**7*kd2**3*(6183.*omega1**4 - 12594.*omega1**3*omega2 + 23816.*omega1**2*omega2**2 - 22636.*omega1*omega2**3 + 11124.*omega2**4) + &
              kd1**4*kd2**6*(35978.*omega1**4 - 78369.*omega1**3*omega2 + 98706.*omega1**2*omega2**2 - 64313.*omega1*omega2**3 + 29802.*omega2**4) + &
              kd1**6*kd2**4*(29802.*omega1**4 - 64313.*omega1**3*omega2 + 98706.*omega1**2*omega2**2 - 78369.*omega1*omega2**3 + 35978.*omega2**4))*swd**2) + &
              134217728.*grav*(64.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(kd1**8 + 54.*kd1**7*kd2 - 172.*kd1**6*kd2**2 + 66.*kd1**5*kd2**3 + 83.*kd1**4*kd2**4 + &
              66.*kd1**3*kd2**5 - 172.*kd1**2*kd2**6 + 54.*kd1*kd2**7 + kd2**8) + (-(kd2**15*omega1**3*omega2) + kd1**15*omega1*omega2**3 + kd1**14*kd2*omega2**2*(-84.*omega1**2 + &
              59.*omega1*omega2 - 32.*omega2**2) + kd1*kd2**14*omega1**2*(32.*omega1**2 - 59.*omega1*omega2 + &
              84.*omega2**2) + kd1**2*kd2**13*omega1*(-1300.*omega1**3 + 2885.*omega1**2*omega2 - 2378.*omega1*omega2**2 + 126.*omega2**3) + kd1**13*kd2**2*omega2*(-126.*omega1**3 + &
              2378.*omega1**2*omega2 - 2885.*omega1*omega2**2 + 1300.*omega2**3) + kd1**9*kd2**6*(4620.*omega1**4 - 3547.*omega1**3*omega2 - 30198.*omega1**2*omega2**2 + &
              165771.*omega1*omega2**3 - 72052.*omega2**4) + 2.*kd1**7*kd2**8*(-29152.*omega1**4 + 67283.*omega1**3*omega2 + 5991.*omega1**2*omega2**2 + 22134.*omega1*omega2**3 - &
              8095.*omega2**4) + kd1**12*kd2**3*(632.*omega1**4 + 375.*omega1**3*omega2 - 8978.*omega1**2*omega2**2 + 11180.*omega1*omega2**3 - 4918.*omega2**4) + &
              kd1**6*kd2**9*(72052.*omega1**4 - 165771.*omega1**3*omega2 + 30198.*omega1**2*omega2**2 + 3547.*omega1*omega2**3 - 4620.*omega2**4) + &
              kd1**3*kd2**12*(4918.*omega1**4 - 11180.*omega1**3*omega2 + 8978.*omega1**2*omega2**2 - 375.*omega1*omega2**3 - 632.*omega2**4) + &
              2.*kd1**4*kd2**11*(1309.*omega1**4 - 3295.*omega1**3*omega2 - 914.*omega1**2*omega2**2 + 421.*omega1*omega2**3 + 909.*omega2**4) - &
              2.*kd1**11*kd2**4*(909.*omega1**4 + 421.*omega1**3*omega2 - 914.*omega1**2*omega2**2 - 3295.*omega1*omega2**3 + 1309.*omega2**4) + kd1**5*kd2**10*(-38174.*omega1**4 + &
              89253.*omega1**3*omega2 - 29112.*omega1**2*omega2**2 - 5293.*omega1*omega2**3 + 1460.*omega2**4) + 2.*kd1**8*kd2**7*(8095.*omega1**4 - &
              22134.*omega1**3*omega2 - 5991.*omega1**2*omega2**2 - 67283.*omega1*omega2**3 + 29152.*omega2**4) + kd1**10*kd2**5*(-1460.*omega1**4 + 5293.*omega1**3*omega2 + &
              29112.*omega1**2*omega2**2 - 89253.*omega1*omega2**3 + 38174.*omega2**4))*swd**2) - &
              2097152.*grav*kd1*(kd1 - kd2)*kd2*(128.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(kd1**8 - 8.*kd1**7*kd2 + 18.*kd1**6*kd2**2 - 59.*kd1**5*kd2**3 + 119.*kd1**4*kd2**4 - &
              59.*kd1**3*kd2**5 + 18.*kd1**2*kd2**6 - 8.*kd1*kd2**7 + kd2**8) + (4.*kd1**14*omega2**2*(7.*omega1**2 - 12.*omega1*omega2 + 6.*omega2**2) + &
              4.*kd2**14*omega1**2*(6.*omega1**2 - 12.*omega1*omega2 + 7.*omega2**2) - kd1*kd2**13*omega1*(326.*omega1**3 - 612.*omega1**2*omega2 + 320.*omega1*omega2**2 + &
              19.*omega2**3) - kd1**13*kd2*omega2*(19.*omega1**3 + 320.*omega1**2*omega2 - 612.*omega1*omega2**2 + 326.*omega2**3) + kd1**7*kd2**7*(-79300.*omega1**4 + &
              135409.*omega1**3*omega2 - 234710.*omega1**2*omega2**2 + 135409.*omega1*omega2**3 - 79300.*omega2**4) + kd1**11*kd2**3*(-3628.*omega1**4 + 7772.*omega1**3*omega2 - &
              23028.*omega1**2*omega2**2 + 19929.*omega1*omega2**3 - 9018.*omega2**4) + kd1**3*kd2**11*(-9018.*omega1**4 + 19929.*omega1**3*omega2 - 23028.*omega1**2*omega2**2 + &
              7772.*omega1*omega2**3 - 3628.*omega2**4) + kd1**2*kd2**12*(1562.*omega1**4 - 3357.*omega1**3*omega2 + 2888.*omega1**2*omega2**2 - 453.*omega1*omega2**3 + &
              208.*omega2**4) + kd1**12*kd2**2*(208.*omega1**4 - 453.*omega1**3*omega2 + 2888.*omega1**2*omega2**2 - 3357.*omega1*omega2**3 + 1562.*omega2**4) - &
              2.*kd1**5*kd2**9*(43299.*omega1**4 - 85405.*omega1**3*omega2 + 91010.*omega1**2*omega2**2 - 24676.*omega1*omega2**3 + 15515.*omega2**4) + kd1**4*kd2**10*(38260.*omega1**4 - &
              79037.*omega1**3*omega2 + 89260.*omega1**2*omega2**2 - 29359.*omega1*omega2**3 + 15618.*omega2**4) + kd1**10*kd2**4*(15618.*omega1**4 - 29359.*omega1**3*omega2 + &
              89260.*omega1**2*omega2**2 - 79037.*omega1*omega2**3 + 38260.*omega2**4) - 2.*kd1**9*kd2**5*(15515.*omega1**4 - 24676.*omega1**3*omega2 + 91010.*omega1**2*omega2**2 - &
              85405.*omega1*omega2**3 + 43299.*omega2**4) + kd1**6*kd2**8*(107100.*omega1**4 - 201741.*omega1**3*omega2 + 230554.*omega1**2*omega2**2 - 69877.*omega1*omega2**3 + &
              47128.*omega2**4) + kd1**8*kd2**6*(47128.*omega1**4 - 69877.*omega1**3*omega2 + 230554.*omega1**2*omega2**2 - 201741.*omega1*omega2**3 + 107100.*omega2**4))*swd**2)))/ &
              ((-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*(omega1 - omega2)*omega2*swd**3*(-64.*grav*(64. + (kd1 - kd2)**2)*(4096. + (384. + &
              (kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2 + (16777216. + (7340032. + (286720. + (1792. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function velsb14
!
real function velsp14()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic velocity of 1st layer for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp14')
    !
    velsp14 = (-4.*(-144115188075855872.*grav*(kd1 + kd2)**3*omega1*omega2*(kd2**2*omega1**2 + 8.*kd1*kd2*omega1*omega2 + kd1**2*omega2**2)*swd**2 - &
              16.*grav*kd1**7*kd2**7*(kd1 + kd2)**7*((9.*kd1**2 + 29.*kd1*kd2 + 16.*kd2**2)*omega1**4 + (15.*kd1**2 - 16.*kd1*kd2 + 45.*kd2**2)*omega1**3*omega2 + &
              (37.*kd1**2 - 86.*kd1*kd2 + 37.*kd2**2)*omega1**2*omega2**2 + (45.*kd1**2 - 16.*kd1*kd2 + 15.*kd2**2)*omega1*omega2**3 + (16.*kd1**2 + 29.*kd1*kd2 + &
              9.*kd2**2)*omega2**4)*swd**2 - kd1**7*kd2**7*(kd1 + kd2)**7*omega1*omega2*(omega1 + omega2)**2*(kd2**2*omega1*(3.*omega1 - omega2) + &
              kd1**2*omega2*(-omega1 + 3.*omega2) - 2.*kd1*kd2*(omega1**2 + 6.*omega1*omega2 + omega2**2))*swd**3 + &
              18014398509481984.*swd*(8.*grav**2*kd1*kd2*(kd1 + kd2)*(6.*kd2**4*omega1**2 + 6.*kd1**4*omega2**2 + kd1*kd2**3*omega1*(12.*omega1 + omega2) + &
              kd1**3*kd2*omega2*(omega1 + 12.*omega2) + 5.*kd1**2*kd2**2*(omega1**2 + omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(kd2**2*(2.*kd1**3 + &
              11.*kd1**2*kd2 + 3.*kd1*kd2**2 + kd2**3)*omega1**2 + kd1*kd2*(kd1 + kd2)*(8.*kd1**2 + 19.*kd1*kd2 + 8.*kd2**2)*omega1*omega2 + kd1**2*(kd1**3 + &
              3.*kd1**2*kd2 + 11.*kd1*kd2**2 + 2.*kd2**3)*omega2**2)*swd**2) + 281474976710656.*swd*(8.*grav**2*kd1*kd2*(kd1 + kd2)*(17.*kd2**6*omega1**2 + &
              17.*kd1**6*omega2**2 - kd1**5*kd2*omega2*(23.*omega1 + 10.*omega2) - kd1*kd2**5*omega1*(10.*omega1 + 23.*omega2) + kd1**4*kd2**2*(94.*omega1**2 - &
              41.*omega1*omega2 + 29.*omega2**2) + 4.*kd1**3*kd2**3*(37.*omega1**2 - 10.*omega1*omega2 + 37.*omega2**2) + kd1**2*kd2**4*(29.*omega1**2 - &
              41.*omega1*omega2 + 94.*omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(kd2**2*(8.*kd1**5 - 176.*kd1**4*kd2 - 382.*kd1**3*kd2**2 - 151.*kd1**2*kd2**3 + &
              17.*kd1*kd2**4 + 5.*kd2**5)*omega1**2 - kd1*kd2*(kd1 + kd2)*(116.*kd1**4 + 249.*kd1**3*kd2 + 293.*kd1**2*kd2**2 + 249.*kd1*kd2**3 + &
              116.*kd2**4)*omega1*omega2 + kd1**2*(5.*kd1**5 + 17.*kd1**4*kd2 - 151.*kd1**3*kd2**2 - 382.*kd1**2*kd2**3 - &
              176.*kd1*kd2**4 + 8.*kd2**5)*omega2**2)*swd**2) - 4398046511104.*swd*(8.*grav**2*kd1*kd2*(kd1 + kd2)*(10.*kd2**8*omega1**2 + 10.*kd1**8*omega2**2 + &
              2.*kd1*kd2**7*omega1*(109.*omega1 + 7.*omega2) + 2.*kd1**7*kd2*omega2*(7.*omega1 + 109.*omega2) + kd1**2*kd2**6*(376.*omega1**2 - 113.*omega1*omega2 - &
              183.*omega2**2) + 2.*kd1**3*kd2**5*(369.*omega1**2 + 62.*omega1*omega2 + 129.*omega2**2) + 5.*kd1**4*kd2**4*(203.*omega1**2 + 102.*omega1*omega2 + &
              203.*omega2**2) + 2.*kd1**5*kd2**3*(129.*omega1**2 + 62.*omega1*omega2 + 369.*omega2**2) + kd1**6*kd2**2*(-183.*omega1**2 - 113.*omega1*omega2 + &
              376.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(kd2**2*(246.*kd1**7 + 1081.*kd1**6*kd2 + 1107.*kd1**5*kd2**2 + 301.*kd1**4*kd2**3 + &
              415.*kd1**3*kd2**4 + 288.*kd1**2*kd2**5 - 38.*kd1*kd2**6 - 10.*kd2**7)*omega1**2 - 2.*kd1*kd2*(kd1 + kd2)*(26.*kd1**6 - 334.*kd1**5*kd2 + &
              43.*kd1**4*kd2**2 + 861.*kd1**3*kd2**3 + 43.*kd1**2*kd2**4 - 334.*kd1*kd2**5 + 26.*kd2**6)*omega1*omega2 + kd1**2*(-10.*kd1**7 - 38.*kd1**6*kd2 + &
              288.*kd1**5*kd2**2 + 415.*kd1**4*kd2**3 + 301.*kd1**3*kd2**4 + 1107.*kd1**2*kd2**5 + 1081.*kd1*kd2**6 + 246.*kd2**7)*omega2**2)*swd**2) - &
              2251799813685248.*grav*(320.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3 + (2.*kd1*kd2**3*(28.*kd1**3 + 85.*kd1**2*kd2 + 78.*kd1*kd2**2 + 24.*kd2**3)*omega1**4 + &
              kd2**2*(-18.*kd1**5 + 131.*kd1**4*kd2 + 406.*kd1**3*kd2**2 + 329.*kd1**2*kd2**3 + 79.*kd1*kd2**4 - 5.*kd2**5)*omega1**3*omega2 + &
              2.*kd1*kd2*(kd1 + kd2)*(kd1**2 + kd2**2)*(82.*kd1**2 + 155.*kd1*kd2 + 82.*kd2**2)*omega1**2*omega2**2 + kd1**2*(-5.*kd1**5 + 79.*kd1**4*kd2 + &
              329.*kd1**3*kd2**2 + 406.*kd1**2*kd2**3 + 131.*kd1*kd2**4 - 18.*kd2**5)*omega1*omega2**3 + 2.*kd1**3*kd2*(24.*kd1**3 + 78.*kd1**2*kd2 + 85.*kd1*kd2**2 + &
              28.*kd2**3)*omega2**4)*swd**2) - 512.*grav*kd1**5*kd2**5*(kd1 + kd2)**5*(128.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4 + (2.*(kd1 + kd2)*(43.*kd1**5 + &
              177.*kd1**4*kd2 + 391.*kd1**3*kd2**2 + 878.*kd1**2*kd2**3 + 781.*kd1*kd2**4 + 236.*kd2**5)*omega1**4 + (178.*kd1**6 + 721.*kd1**5*kd2 + 2391.*kd1**4*kd2**2 + &
              5222.*kd1**3*kd2**3 + 6491.*kd1**2*kd2**4 + 3406.*kd1*kd2**5 + 838.*kd2**6)*omega1**3*omega2 + 2.*(259.*kd1**6 + 1086.*kd1**5*kd2 + 2540.*kd1**4*kd2**2 + &
              2877.*kd1**3*kd2**3 + 2540.*kd1**2*kd2**4 + 1086.*kd1*kd2**5 + 259.*kd2**6)*omega1**2*omega2**2 + (838.*kd1**6 + 3406.*kd1**5*kd2 + 6491.*kd1**4*kd2**2 + &
              5222.*kd1**3*kd2**3 + 2391.*kd1**2*kd2**4 + 721.*kd1*kd2**5 + 178.*kd2**6)*omega1*omega2**3 + 2.*(kd1 + kd2)*(236.*kd1**5 + 781.*kd1**4*kd2 + &
              878.*kd1**3*kd2**2 + 391.*kd1**2*kd2**3 + 177.*kd1*kd2**4 + 43.*kd2**5)*omega2**4)*swd**2) - &
              35184372088832.*grav*(384.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(2.*kd1**2 - 3.*kd1*kd2 + 2.*kd2**2) + (2.*kd1*kd2**3*(360.*kd1**5 + 1134.*kd1**4*kd2 + &
              1205.*kd1**3*kd2**2 + 494.*kd1**2*kd2**3 + 143.*kd1*kd2**4 + 68.*kd2**5)*omega1**4 + kd2**2*(-405.*kd1**7 - 232.*kd1**6*kd2 + 2011.*kd1**5*kd2**2 + &
              2817.*kd1**4*kd2**3 + 969.*kd1**3*kd2**4 + 338.*kd1**2*kd2**5 + 310.*kd1*kd2**6 + 10.*kd2**7)*omega1**3*omega2 + 2.*kd1*kd2*(kd1 + kd2)*(94.*kd1**6 - &
              492.*kd1**5*kd2 + 38.*kd1**4*kd2**2 + 1333.*kd1**3*kd2**3 + 38.*kd1**2*kd2**4 - 492.*kd1*kd2**5 + 94.*kd2**6)*omega1**2*omega2**2 + kd1**2*(10.*kd1**7 + &
              310.*kd1**6*kd2 + 338.*kd1**5*kd2**2 + 969.*kd1**4*kd2**3 + 2817.*kd1**3*kd2**4 + 2011.*kd1**2*kd2**5 - 232.*kd1*kd2**6 - 405.*kd2**7)*omega1*omega2**3 + &
              2.*kd1**3*kd2*(68.*kd1**5 + 143.*kd1**4*kd2 + 494.*kd1**3*kd2**2 + 1205.*kd1**2*kd2**3 + 1134.*kd1*kd2**4 + 360.*kd2**5)*omega2**4)*swd**2) + &
              1099511627776.*grav*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(2.*kd1**4 + 55.*kd1**3*kd2 + 50.*kd1**2*kd2**2 + 55.*kd1*kd2**3 + 2.*kd2**4) + &
              (kd1*kd2**3*(-148.*kd1**7 + 1113.*kd1**6*kd2 + 3702.*kd1**5*kd2**2 + 4145.*kd1**4*kd2**3 + 3258.*kd1**3*kd2**4 + 2424.*kd1**2*kd2**5 + 928.*kd1*kd2**6 + &
              40.*kd2**7)*omega1**4 + kd2**2*(38.*kd1**9 + 749.*kd1**8*kd2 + 5435.*kd1**7*kd2**2 + 10814.*kd1**6*kd2**3 + 9430.*kd1**5*kd2**4 + 6286.*kd1**4*kd2**5 + &
              4621.*kd1**3*kd2**6 + 1807.*kd1**2*kd2**7 + 101.*kd1*kd2**8 + 5.*kd2**9)*omega1**3*omega2 + kd1*kd2*(kd1 + kd2)*(124.*kd1**8 + 1950.*kd1**7*kd2 + &
              5108.*kd1**6*kd2**2 + 11657.*kd1**5*kd2**3 + 16642.*kd1**4*kd2**4 + 11657.*kd1**3*kd2**5 + 5108.*kd1**2*kd2**6 + 1950.*kd1*kd2**7 + &
              124.*kd2**8)*omega1**2*omega2**2 + kd1**2*(5.*kd1**9 + 101.*kd1**8*kd2 + 1807.*kd1**7*kd2**2 + 4621.*kd1**6*kd2**3 + 6286.*kd1**5*kd2**4 + &
              9430.*kd1**4*kd2**5 + 10814.*kd1**3*kd2**6 + 5435.*kd1**2*kd2**7 + 749.*kd1*kd2**8 + 38.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(40.*kd1**7 + &
              928.*kd1**6*kd2 + 2424.*kd1**5*kd2**2 + 3258.*kd1**4*kd2**3 + 4145.*kd1**3*kd2**4 + 3702.*kd1**2*kd2**5 + 1113.*kd1*kd2**6 - &
              148.*kd2**7)*omega2**4)*swd**2) + 8589934592.*grav*(128.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(7.*kd1**6 + 37.*kd1**5*kd2 + 51.*kd1**4*kd2**2 + &
              81.*kd1**3*kd2**3 + 51.*kd1**2*kd2**4 + 37.*kd1*kd2**5 + 7.*kd2**6) + (2.*kd1*kd2**3*(336.*kd1**9 + 5388.*kd1**8*kd2 + 20023.*kd1**7*kd2**2 + &
              36685.*kd1**6*kd2**3 + 40472.*kd1**5*kd2**4 + 27645.*kd1**4*kd2**5 + 10899.*kd1**3*kd2**6 + 2574.*kd1**2*kd2**7 + 630.*kd1*kd2**8 + 80.*kd2**9)*omega1**4 + &
              kd2**2*(-391.*kd1**11 - 739.*kd1**10*kd2 + 16127.*kd1**9*kd2**2 + 65207.*kd1**8*kd2**3 + 115523.*kd1**7*kd2**4 + 127615.*kd1**6*kd2**5 + 94299.*kd1**5*kd2**6 + &
              42507.*kd1**4*kd2**7 + 11427.*kd1**3*kd2**8 + 2649.*kd1**2*kd2**9 + 297.*kd1*kd2**10 - 5.*kd2**11)*omega1**3*omega2 + 2.*kd1*kd2*(kd1 + kd2)*(112.*kd1**10 + &
              52.*kd1**9*kd2 + 726.*kd1**8*kd2**2 + 17491.*kd1**7*kd2**3 + 49247.*kd1**6*kd2**4 + 65058.*kd1**5*kd2**5 + 49247.*kd1**4*kd2**6 + 17491.*kd1**3*kd2**7 + &
              726.*kd1**2*kd2**8 + 52.*kd1*kd2**9 + 112.*kd2**10)*omega1**2*omega2**2 + kd1**2*(-5.*kd1**11 + 297.*kd1**10*kd2 + 2649.*kd1**9*kd2**2 + 11427.*kd1**8*kd2**3 + &
              42507.*kd1**7*kd2**4 + 94299.*kd1**6*kd2**5 + 127615.*kd1**5*kd2**6 + 115523.*kd1**4*kd2**7 + 65207.*kd1**3*kd2**8 + 16127.*kd1**2*kd2**9 - 739.*kd1*kd2**10 - &
              391.*kd2**11)*omega1*omega2**3 + 2.*kd1**3*kd2*(80.*kd1**9 + 630.*kd1**8*kd2 + 2574.*kd1**7*kd2**2 + 10899.*kd1**6*kd2**3 + 27645.*kd1**5*kd2**4 + &
              40472.*kd1**4*kd2**5 + 36685.*kd1**3*kd2**6 + 20023.*kd1**2*kd2**7 + 5388.*kd1*kd2**8 + 336.*kd2**9)*omega2**4)*swd**2) + &
              64.*kd1**5*kd2**5*(kd1 + kd2)**5*swd*(8.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(kd1**2*(omega1 + omega2)*(27.*omega1 + 20.*omega2) + &
              kd2**2*(omega1 + omega2)*(20.*omega1 + 27.*omega2) + 2.*kd1*kd2*(24.*omega1**2 + 41.*omega1*omega2 + 24.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(kd1**4*kd2**2*(3.*omega1 - 407.*omega2)*(omega1 - omega2) + kd1**2*kd2**4*(407.*omega1 - 3.*omega2)*(omega1 - omega2) - &
              60.*kd1**6*omega1*omega2 - 60.*kd2**6*omega1*omega2 + 6.*kd1*kd2**5*(58.*omega1**2 - 45.*omega1*omega2 + omega2**2) - 4.*kd1**3*kd2**3*(omega1**2 + &
              130.*omega1*omega2 + omega2**2) + 6.*kd1**5*kd2*(omega1**2 - 45.*omega1*omega2 + 58.*omega2**2))*swd**2) + &
              4096.*kd1**3*kd2**3*(kd1 + kd2)**3*swd*(8.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*((kd1 + kd2)*(84.*kd1**5 + 356.*kd1**4*kd2 + 520.*kd1**3*kd2**2 + &
              584.*kd1**2*kd2**3 + 521.*kd1*kd2**4 + 155.*kd2**5)*omega1**2 + (156.*kd1**6 + 729.*kd1**5*kd2 + 1359.*kd1**4*kd2**2 + 1544.*kd1**3*kd2**3 + &
              1359.*kd1**2*kd2**4 + 729.*kd1*kd2**5 + 156.*kd2**6)*omega1*omega2 + (kd1 + kd2)*(155.*kd1**5 + 521.*kd1**4*kd2 + 584.*kd1**3*kd2**2 + &
              520.*kd1**2*kd2**3 + 356.*kd1*kd2**4 + 84.*kd2**5)*omega2**2) + omega1*omega2*(omega1 + omega2)**2*(3.*kd1**10*omega2*(21.*omega1 + omega2) + &
              3.*kd2**10*omega1*(omega1 + 21.*omega2) + kd1*kd2**9*(61.*omega1**2 + 556.*omega1*omega2 + 6.*omega2**2) + kd1**2*kd2**8*(250.*omega1**2 + &
              3193.*omega1*omega2 + 15.*omega2**2) + kd1**3*kd2**7*(526.*omega1**2 + 11167.*omega1*omega2 + 24.*omega2**2) + kd1**9*kd2*(6.*omega1**2 + &
              556.*omega1*omega2 + 61.*omega2**2) + kd1**8*kd2**2*(15.*omega1**2 + 3193.*omega1*omega2 + 250.*omega2**2) + 8.*kd1**5*kd2**5*(411.*omega1**2 + &
              3841.*omega1*omega2 + 411.*omega2**2) + kd1**7*kd2**3*(24.*omega1**2 + 11167.*omega1*omega2 + 526.*omega2**2) + 2.*kd1**4*kd2**6*(1156.*omega1**2 + &
              11886.*omega1*omega2 + 633.*omega2**2) + 2.*kd1**6*kd2**4*(633.*omega1**2 + 11886.*omega1*omega2 + 1156.*omega2**2))*swd**2) + &
              16777216.*swd*(8.*grav**2*kd1*kd2*(kd1 + kd2)**2*(3.*kd2**13*omega1**2 + kd1*kd2**12*omega1*(29.*omega1 - 3.*omega2) + 3.*kd1**13*omega2**2 + &
              kd1**12*kd2*omega2*(-3.*omega1 + 29.*omega2) + 2.*kd1**2*kd2**11*(45.*omega1**2 + 13.*omega1*omega2 + 23.*omega2**2) + 2.*kd1**11*kd2**2*(23.*omega1**2 + &
              13.*omega1*omega2 + 45.*omega2**2) + 2.*kd1**3*kd2**10*(654.*omega1**2 + 483.*omega1*omega2 + 431.*omega2**2) + 2.*kd1**10*kd2**3*(431.*omega1**2 + &
              483.*omega1*omega2 + 654.*omega2**2) + 3.*kd1**4*kd2**9*(1873.*omega1**2 + 1213.*omega1*omega2 + 1101.*omega2**2) + 3.*kd1**9*kd2**4*(1101.*omega1**2 + &
              1213.*omega1*omega2 + 1873.*omega2**2) + 2.*kd1**6*kd2**7*(3125.*omega1**2 - 42.*omega1*omega2 + 1912.*omega2**2) + 2.*kd1**7*kd2**6*(1912.*omega1**2 - &
              42.*omega1*omega2 + 3125.*omega2**2) + kd1**5*kd2**8*(8853.*omega1**2 + 3811.*omega1*omega2 + 4545.*omega2**2) + kd1**8*kd2**5*(4545.*omega1**2 + &
              3811.*omega1*omega2 + 8853.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(-(kd2**15*omega1**2) + kd1**14*kd2*(52.*omega1 - 5.*omega2)*omega2 - &
              kd1**15*omega2**2 + kd1*kd2**14*omega1*(-5.*omega1 + 52.*omega2) + kd1**9*kd2**6*(-5223.*omega1**2 + 48969.*omega1*omega2 - 20273.*omega2**2) + &
              kd1**8*kd2**7*(-12129.*omega1**2 + 97964.*omega1*omega2 - 18712.*omega2**2) + kd1**7*kd2**8*(-18712.*omega1**2 + 97964.*omega1*omega2 - 12129.*omega2**2) + &
              kd1**10*kd2**5*(721.*omega1**2 + 12931.*omega1*omega2 - 10887.*omega2**2) + kd1**6*kd2**9*(-20273.*omega1**2 + 48969.*omega1*omega2 - 5223.*omega2**2) + &
              kd1**11*kd2**4*(1770.*omega1**2 + 5454.*omega1*omega2 - 1367.*omega2**2) + kd1**2*kd2**13*(27.*omega1**2 + 961.*omega1*omega2 + 8.*omega2**2) + &
              kd1**13*kd2**2*(8.*omega1**2 + 961.*omega1*omega2 + 27.*omega2**2) + 2.*kd1**3*kd2**12*(245.*omega1**2 + 1800.*omega1*omega2 + 228.*omega2**2) + &
              2.*kd1**12*kd2**3*(228.*omega1**2 + 1800.*omega1*omega2 + 245.*omega2**2) + kd1**5*kd2**10*(-10887.*omega1**2 + 12931.*omega1*omega2 + 721.*omega2**2) + &
              kd1**4*kd2**11*(-1367.*omega1**2 + 5454.*omega1*omega2 + 1770.*omega2**2))*swd**2) + &
              262144.*kd1*kd2*(kd1 + kd2)*swd*(8.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(kd2**10*(78.*omega1**2 + 53.*omega1*omega2 + 33.*omega2**2) + 4.*kd1*kd2**9*(163.*omega1**2 + &
              109.*omega1*omega2 + 65.*omega2**2) + kd1**10*(33.*omega1**2 + 53.*omega1*omega2 + 78.*omega2**2) + 4.*kd1**9*kd2*(65.*omega1**2 + 109.*omega1*omega2 + &
              163.*omega2**2) + kd1**2*kd2**8*(2101.*omega1**2 + 1333.*omega1*omega2 + 972.*omega2**2) + kd1**8*kd2**2*(972.*omega1**2 + 1333.*omega1*omega2 + 2101.*omega2**2) + &
              kd1**3*kd2**7*(4660.*omega1**2 + 3735.*omega1*omega2 + 3332.*omega2**2) + kd1**7*kd2**3*(3332.*omega1**2 + 3735.*omega1*omega2 + 4660.*omega2**2) + &
              2.*kd1**5*kd2**5*(5595.*omega1**2 + 6406.*omega1*omega2 + 5595.*omega2**2) + kd1**4*kd2**6*(8761.*omega1**2 + 9195.*omega1*omega2 + 8150.*omega2**2) + &
              kd1**6*kd2**4*(8150.*omega1**2 + 9195.*omega1*omega2 + 8761.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(4.*kd1**14*omega1*omega2 + &
              4.*kd2**14*omega1*omega2 + kd1**5*kd2**9*(2803.*omega1**2 + 80713.*omega1*omega2 - 5412.*omega2**2) + 2.*kd1**6*kd2**8*(643.*omega1**2 + 56359.*omega1*omega2 - &
              3603.*omega2**2) + kd1**4*kd2**10*(2551.*omega1**2 + 37970.*omega1*omega2 - 1339.*omega2**2) + kd1**12*kd2**2*(-9.*omega1**2 + 1194.*omega1*omega2 - &
              331.*omega2**2) + kd1**3*kd2**11*(571.*omega1**2 + 9908.*omega1*omega2 - 36.*omega2**2) + kd1**2*kd2**12*(-331.*omega1**2 + 1194.*omega1*omega2 - &
              9.*omega2**2) - 2.*kd1*kd2**13*(43.*omega1**2 - 35.*omega1*omega2 + omega2**2) - 2.*kd1**13*kd2*(omega1**2 - 35.*omega1*omega2 + 43.*omega2**2) + &
              kd1**11*kd2**3*(-36.*omega1**2 + 9908.*omega1*omega2 + 571.*omega2**2) + 2.*kd1**8*kd2**6*(-3603.*omega1**2 + 56359.*omega1*omega2 + 643.*omega2**2) - &
              2.*kd1**7*kd2**7*(1441.*omega1**2 - 61199.*omega1*omega2 + 1441.*omega2**2) + kd1**10*kd2**4*(-1339.*omega1**2 + 37970.*omega1*omega2 + 2551.*omega2**2) + &
              kd1**9*kd2**5*(-5412.*omega1**2 + 80713.*omega1*omega2 + 2803.*omega2**2))*swd**2) - 68719476736.*swd*(16.*grav**2*kd1*kd2*(kd1 + kd2)*(10.*kd2**10*omega1**2 + &
              13.*kd1*kd2**9*omega1*(3.*omega1 - omega2) - 13.*kd1**9*kd2*(omega1 - 3.*omega2)*omega2 + 10.*kd1**10*omega2**2 + kd1**2*kd2**8*(80.*omega1**2 - &
              69.*omega1*omega2 + 22.*omega2**2) + kd1**8*kd2**2*(22.*omega1**2 - 69.*omega1*omega2 + 80.*omega2**2) + kd1**3*kd2**7*(1163.*omega1**2 + 400.*omega1*omega2 + &
              808.*omega2**2) + kd1**7*kd2**3*(808.*omega1**2 + 400.*omega1*omega2 + 1163.*omega2**2) + kd1**4*kd2**6*(3131.*omega1**2 + 1148.*omega1*omega2 + 2507.*omega2**2) + &
              kd1**6*kd2**4*(2507.*omega1**2 + 1148.*omega1*omega2 + 3131.*omega2**2) + kd1**5*kd2**5*(3743.*omega1**2 + 1390.*omega1*omega2 + 3743.*omega2**2)) + &
              omega1*omega2*(omega1 + omega2)**2*(10.*kd2**11*omega1**2 + 10.*kd1**11*omega2**2 + 42.*kd1**10*kd2*omega2*(4.*omega1 + omega2) + &
              42.*kd1*kd2**10*omega1*(omega1 + 4.*omega2) + 6.*kd1**9*kd2**2*omega2*(365.*omega1 + 11.*omega2) + 6.*kd1**2*kd2**9*omega1*(11.*omega1 + 365.*omega2) + &
              kd1**7*kd2**4*(4542.*omega1**2 + 17429.*omega1*omega2 + 603.*omega2**2) + kd1**8*kd2**3*(1304.*omega1**2 + 7565.*omega1*omega2 + 659.*omega2**2) + &
              kd1**6*kd2**5*(4409.*omega1**2 + 28106.*omega1*omega2 + 1079.*omega2**2) + kd1**3*kd2**8*(659.*omega1**2 + 7565.*omega1*omega2 + 1304.*omega2**2) + &
              kd1**5*kd2**6*(1079.*omega1**2 + 28106.*omega1*omega2 + 4409.*omega2**2) + kd1**4*kd2**7*(603.*omega1**2 + 17429.*omega1*omega2 + &
              4542.*omega2**2))*swd**2) + 1073741824.*swd*(8.*grav**2*kd1*kd2*(kd1 + kd2)*(4.*kd2**12*omega1**2 + 4.*kd1**12*omega2**2 + kd1*kd2**11*omega1*(134.*omega1 + &
              13.*omega2) + kd1**11*kd2*omega2*(13.*omega1 + 134.*omega2) + kd1**2*kd2**10*(241.*omega1**2 - 286.*omega1*omega2 - 317.*omega2**2) + &
              kd1**10*kd2**2*(-317.*omega1**2 - 286.*omega1*omega2 + 241.*omega2**2) - 2.*kd1**4*kd2**8*(3443.*omega1**2 + 3462.*omega1*omega2 + 1426.*omega2**2) - &
              kd1**3*kd2**9*(1834.*omega1**2 + 2523.*omega1*omega2 + 1568.*omega2**2) - kd1**9*kd2**3*(1568.*omega1**2 + 2523.*omega1*omega2 + 1834.*omega2**2) - &
              2.*kd1**8*kd2**4*(1426.*omega1**2 + 3462.*omega1*omega2 + 3443.*omega2**2) - kd1**5*kd2**7*(10306.*omega1**2 + 11321.*omega1*omega2 + 4646.*omega2**2) - &
              kd1**6*kd2**6*(8413.*omega1**2 + 13248.*omega1*omega2 + 8413.*omega2**2) - kd1**7*kd2**5*(4646.*omega1**2 + 11321.*omega1*omega2 + 10306.*omega2**2)) + &
              omega1*omega2*(omega1 + omega2)**2*(5.*kd2**13*omega1**2 + kd1*kd2**12*omega1*(23.*omega1 - 64.*omega2) + 5.*kd1**13*omega2**2 + &
              kd1**12*kd2*omega2*(-64.*omega1 + 23.*omega2) + kd1**2*kd2**11*(251.*omega1**2 + 355.*omega1*omega2 + 246.*omega2**2) + kd1**11*kd2**2*(246.*omega1**2 + &
              355.*omega1*omega2 + 251.*omega2**2) + kd1**10*kd2**3*(1753.*omega1**2 + 1115.*omega1*omega2 + 614.*omega2**2) + kd1**3*kd2**10*(614.*omega1**2 + 1115.*omega1*omega2 + &
              1753.*omega2**2) + kd1**9*kd2**4*(4301.*omega1**2 - 14902.*omega1*omega2 + 1785.*omega2**2) + 5.*kd1**6*kd2**7*(3859.*omega1**2 - 25027.*omega1*omega2 + &
              3372.*omega2**2) + 5.*kd1**7*kd2**6*(3372.*omega1**2 - 25027.*omega1*omega2 + 3859.*omega2**2) + kd1**4*kd2**9*(1785.*omega1**2 - 14902.*omega1*omega2 + &
              4301.*omega2**2) + kd1**5*kd2**8*(9517.*omega1**2 - 67890.*omega1*omega2 + 8507.*omega2**2) + kd1**8*kd2**5*(8507.*omega1**2 - 67890.*omega1*omega2 + &
              9517.*omega2**2))*swd**2) - 32768.*grav*kd1**3*kd2**3*(kd1 + kd2)**3*(64.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(12.*kd1**4 + 32.*kd1**3*kd2 + 21.*kd1**2*kd2**2 + &
              32.*kd1*kd2**3 + 12.*kd2**4) + (4.*kd2**10*(124.*omega1**4 + 268.*omega1**3*omega2 + 165.*omega1**2*omega2**2 + 12.*omega1*omega2**3 + 6.*omega2**4) + &
              4.*kd1**10*(6.*omega1**4 + 12.*omega1**3*omega2 + 165.*omega1**2*omega2**2 + 268.*omega1*omega2**3 + 124.*omega2**4) + kd1*kd2**9*(3726.*omega1**4 + &
              7643.*omega1**3*omega2 + 4852.*omega1**2*omega2**2 + 518.*omega1*omega2**3 + 282.*omega2**4) + 2.*kd1**2*kd2**8*(5583.*omega1**4 + 10995.*omega1**3*omega2 + &
              8498.*omega1**2*omega2**2 + 2542.*omega1*omega2**3 + 1347.*omega2**4) + kd1**9*kd2*(282.*omega1**4 + 518.*omega1**3*omega2 + 4852.*omega1**2*omega2**2 + &
              7643.*omega1*omega2**3 + 3726.*omega2**4) + 2.*kd1**8*kd2**2*(1347.*omega1**4 + 2542.*omega1**3*omega2 + 8498.*omega1**2*omega2**2 + 10995.*omega1*omega2**3 + &
              5583.*omega2**4) + 2.*kd1**3*kd2**7*(11124.*omega1**4 + 22636.*omega1**3*omega2 + 23816.*omega1**2*omega2**2 + 12594.*omega1*omega2**3 + 6183.*omega2**4) + &
              2.*kd1**7*kd2**3*(6183.*omega1**4 + 12594.*omega1**3*omega2 + 23816.*omega1**2*omega2**2 + 22636.*omega1*omega2**3 + 11124.*omega2**4) + &
              kd1**4*kd2**6*(35978.*omega1**4 + 78369.*omega1**3*omega2 + 98706.*omega1**2*omega2**2 + 64313.*omega1*omega2**3 + 29802.*omega2**4) + kd1**6*kd2**4*(29802.*omega1**4 + &
              64313.*omega1**3*omega2 + 98706.*omega1**2*omega2**2 + 78369.*omega1*omega2**3 + 35978.*omega2**4) + kd1**5*kd2**5*(41530.*omega1**4 + 92003.*omega1**3*omega2 + &
              127268.*omega1**2*omega2**2 + 92003.*omega1*omega2**3 + 41530.*omega2**4))*swd**2) + 134217728.*grav*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(kd1**8 - &
              54.*kd1**7*kd2 - 172.*kd1**6*kd2**2 - 66.*kd1**5*kd2**3 + 83.*kd1**4*kd2**4 - 66.*kd1**3*kd2**5 - 172.*kd1**2*kd2**6 - 54.*kd1*kd2**7 + kd2**8) + &
              (kd2**15*omega1**3*omega2 + kd1**15*omega1*omega2**3 - kd1**14*kd2*omega2**2*(84.*omega1**2 + 59.*omega1*omega2 + 32.*omega2**2) - kd1*kd2**14*omega1**2*(32.*omega1**2 + &
              59.*omega1*omega2 + 84.*omega2**2) - kd1**2*kd2**13*omega1*(1300.*omega1**3 + 2885.*omega1**2*omega2 + 2378.*omega1*omega2**2 + 126.*omega2**3) - &
              kd1**13*kd2**2*omega2*(126.*omega1**3 + 2378.*omega1**2*omega2 + 2885.*omega1*omega2**2 + 1300.*omega2**3) + kd1**12*kd2**3*(632.*omega1**4 - 375.*omega1**3*omega2 - &
              8978.*omega1**2*omega2**2 - 11180.*omega1*omega2**3 - 4918.*omega2**4) + kd1**6*kd2**9*(72052.*omega1**4 + 165771.*omega1**3*omega2 + 30198.*omega1**2*omega2**2 - &
              3547.*omega1*omega2**3 - 4620.*omega2**4) + kd1**5*kd2**10*(38174.*omega1**4 + 89253.*omega1**3*omega2 + 29112.*omega1**2*omega2**2 - 5293.*omega1*omega2**3 - &
              1460.*omega2**4) - kd1**3*kd2**12*(4918.*omega1**4 + 11180.*omega1**3*omega2 + 8978.*omega1**2*omega2**2 + 375.*omega1*omega2**3 - 632.*omega2**4) + &
              2.*kd1**4*kd2**11*(1309.*omega1**4 + 3295.*omega1**3*omega2 - 914.*omega1**2*omega2**2 - 421.*omega1*omega2**3 + 909.*omega2**4) + 2.*kd1**11*kd2**4*(909.*omega1**4 - &
              421.*omega1**3*omega2 - 914.*omega1**2*omega2**2 + 3295.*omega1*omega2**3 + 1309.*omega2**4) + 2.*kd1**7*kd2**8*(29152.*omega1**4 + 67283.*omega1**3*omega2 - &
              5991.*omega1**2*omega2**2 + 22134.*omega1*omega2**3 + 8095.*omega2**4) + 2.*kd1**8*kd2**7*(8095.*omega1**4 + 22134.*omega1**3*omega2 - 5991.*omega1**2*omega2**2 + &
              67283.*omega1*omega2**3 + 29152.*omega2**4) + kd1**10*kd2**5*(-1460.*omega1**4 - 5293.*omega1**3*omega2 + 29112.*omega1**2*omega2**2 + 89253.*omega1*omega2**3 + &
              38174.*omega2**4) + kd1**9*kd2**6*(-4620.*omega1**4 - 3547.*omega1**3*omega2 + 30198.*omega1**2*omega2**2 + 165771.*omega1*omega2**3 + 72052.*omega2**4))*swd**2) - &
              2097152.*grav*kd1*kd2*(kd1 + kd2)*(128.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(kd1**8 + 8.*kd1**7*kd2 + 18.*kd1**6*kd2**2 + 59.*kd1**5*kd2**3 + 119.*kd1**4*kd2**4 + &
              59.*kd1**3*kd2**5 + 18.*kd1**2*kd2**6 + 8.*kd1*kd2**7 + kd2**8) + (4.*kd1**14*omega2**2*(7.*omega1**2 + 12.*omega1*omega2 + 6.*omega2**2) + &
              4.*kd2**14*omega1**2*(6.*omega1**2 + 12.*omega1*omega2 + 7.*omega2**2) + kd1*kd2**13*omega1*(326.*omega1**3 + 612.*omega1**2*omega2 + 320.*omega1*omega2**2 - &
              19.*omega2**3) + kd1**13*kd2*omega2*(-19.*omega1**3 + 320.*omega1**2*omega2 + 612.*omega1*omega2**2 + 326.*omega2**3) + kd1**2*kd2**12*(1562.*omega1**4 + &
              3357.*omega1**3*omega2 + 2888.*omega1**2*omega2**2 + 453.*omega1*omega2**3 + 208.*omega2**4) + kd1**12*kd2**2*(208.*omega1**4 + 453.*omega1**3*omega2 + &
              2888.*omega1**2*omega2**2 + 3357.*omega1*omega2**3 + 1562.*omega2**4) + kd1**3*kd2**11*(9018.*omega1**4 + 19929.*omega1**3*omega2 + 23028.*omega1**2*omega2**2 + &
              7772.*omega1*omega2**3 + 3628.*omega2**4) + kd1**11*kd2**3*(3628.*omega1**4 + 7772.*omega1**3*omega2 + 23028.*omega1**2*omega2**2 + 19929.*omega1*omega2**3 + &
              9018.*omega2**4) + 2.*kd1**5*kd2**9*(43299.*omega1**4 + 85405.*omega1**3*omega2 + 91010.*omega1**2*omega2**2 + 24676.*omega1*omega2**3 + 15515.*omega2**4) + &
              kd1**4*kd2**10*(38260.*omega1**4 + 79037.*omega1**3*omega2 + 89260.*omega1**2*omega2**2 + 29359.*omega1*omega2**3 + 15618.*omega2**4) + &
              kd1**10*kd2**4*(15618.*omega1**4 + 29359.*omega1**3*omega2 + 89260.*omega1**2*omega2**2 + 79037.*omega1*omega2**3 + 38260.*omega2**4) + &
              2.*kd1**9*kd2**5*(15515.*omega1**4 + 24676.*omega1**3*omega2 + 91010.*omega1**2*omega2**2 + 85405.*omega1*omega2**3 + 43299.*omega2**4) + &
              kd1**6*kd2**8*(107100.*omega1**4 + 201741.*omega1**3*omega2 + 230554.*omega1**2*omega2**2 + 69877.*omega1*omega2**3 + 47128.*omega2**4) + &
              kd1**7*kd2**7*(79300.*omega1**4 + 135409.*omega1**3*omega2 + 234710.*omega1**2*omega2**2 + 135409.*omega1*omega2**3 + 79300.*omega2**4) + &
              kd1**8*kd2**6*(47128.*omega1**4 + 69877.*omega1**3*omega2 + 230554.*omega1**2*omega2**2 + 201741.*omega1*omega2**3 + 107100.*omega2**4))*swd**2)))/ &
              ((-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*omega2*(omega1 + omega2)*swd**3*(-64.*grav*(kd1 + kd2)**2*(64. + (kd1 + kd2)**2)*(4096. + &
              (kd1 + kd2)**2*(384. + (kd1 + kd2)**2)) + (16777216. + (kd1 + kd2)**2*(7340032. + (kd1 + kd2)**2*(286720. + (kd1 + kd2)**2*(1792. + (kd1 + kd2)**2))))*(omega1 + omega2)**2*swd))
    !
end function velsp14
!
real function velsb24()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic velocity of 2nd layer for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb24')
    !
    velsb24 = (-4.*(-144115188075855872.*grav*(kd1 - kd2)**3*omega1*omega2*(kd2**2*omega1**2 + 8.*kd1*kd2*omega1*omega2 + kd1**2*omega2**2)*swd**2 + &
              16.*grav*kd1**7*(kd1 - kd2)**7*kd2**7*((11.*kd1**2 - 31.*kd1*kd2 + 16.*kd2**2)*omega1**4 + (-21.*kd1**2 + 6.*kd1*kd2 - 55.*kd2**2)*omega1**3*omega2 + &
              (55.*kd1**2 + 38.*kd1*kd2 + 55.*kd2**2)*omega1**2*omega2**2 + (-55.*kd1**2 + 6.*kd1*kd2 - 21.*kd2**2)*omega1*omega2**3 + (16.*kd1**2 - 31.*kd1*kd2 + &
              11.*kd2**2)*omega2**4)*swd**2 + kd1**7*(kd1 - kd2)**7*kd2**7*omega1*(omega1 - omega2)**2*omega2*(kd2**2*omega1*(omega1 - 5.*omega2) + &
              kd1**2*omega2*(-5.*omega1 + omega2) - 2.*kd1*kd2*(omega1**2 - 6.*omega1*omega2 + omega2**2))*swd**3 + &
              18014398509481984.*swd*(8.*grav**2*kd1*(kd1 - kd2)*kd2*(6.*kd2**4*omega1**2 + kd1**3*kd2*(omega1 - 12.*omega2)*omega2 + 6.*kd1**4*omega2**2 + &
              kd1*kd2**3*omega1*(-12.*omega1 + omega2) + 5.*kd1**2*kd2**2*(omega1**2 + omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(kd2**2*(6.*kd1**3 - &
              17.*kd1**2*kd2 + 9.*kd1*kd2**2 - 3.*kd2**3)*omega1**2 + 3.*kd1*(kd1 - kd2)*kd2*(8.*kd1**2 - 19.*kd1*kd2 + 8.*kd2**2)*omega1*omega2 + kd1**2*(3.*kd1**3 - &
              9.*kd1**2*kd2 + 17.*kd1*kd2**2 - 6.*kd2**3)*omega2**2)*swd**2) + 281474976710656.*swd*(8.*grav**2*kd1*(kd1 - kd2)*kd2*(41.*kd2**6*omega1**2 + &
              41.*kd1**6*omega2**2 - kd1*kd2**5*omega1*(118.*omega1 + 11.*omega2) - kd1**5*kd2*omega2*(11.*omega1 + 118.*omega2) - 12.*kd1**3*kd2**3*(43.*omega1**2 - &
              40.*omega1*omega2 + 43.*omega2**2) + kd1**2*kd2**4*(353.*omega1**2 - 215.*omega1*omega2 + 226.*omega2**2) + kd1**4*kd2**2*(226.*omega1**2 - &
              215.*omega1*omega2 + 353.*omega2**2)) - omega1*(omega1 - omega2)**2*omega2*(kd2**2*(8.*kd1**5 + 192.*kd1**4*kd2 - 334.*kd1**3*kd2**2 + 109.*kd1**2*kd2**3 + &
              31.*kd1*kd2**4 - 11.*kd2**5)*omega1**2 + kd1*kd2*(-kd1 + kd2)*(380.*kd1**4 - 923.*kd1**3*kd2 + 1099.*kd1**2*kd2**2 - 923.*kd1*kd2**3 + &
              380.*kd2**4)*omega1*omega2 + kd1**2*(11.*kd1**5 - 31.*kd1**4*kd2 - 109.*kd1**3*kd2**2 + 334.*kd1**2*kd2**3 - 192.*kd1*kd2**4 - &
              8.*kd2**5)*omega2**2)*swd**2) + 68719476736.*swd*(16.*grav**2*kd1*kd2*(-kd1 + kd2)*(-16.*kd2**10*omega1**2 + kd1**9*kd2*(31.*omega1 - 71.*omega2)*omega2 - &
              16.*kd1**10*omega2**2 + kd1*kd2**9*omega1*(-71.*omega1 + 31.*omega2) + kd1**6*kd2**4*(-2779.*omega1**2 + 966.*omega1*omega2 - 3985.*omega2**2) + &
              kd1**4*kd2**6*(-3985.*omega1**2 + 966.*omega1*omega2 - 2779.*omega2**2) + kd1**2*kd2**8*(-340.*omega1**2 + 199.*omega1*omega2 - 414.*omega2**2) + &
              kd1**8*kd2**2*(-414.*omega1**2 + 199.*omega1*omega2 - 340.*omega2**2) + 5.*kd1**5*kd2**5*(917.*omega1**2 - 374.*omega1*omega2 + 917.*omega2**2) + &
              kd1**3*kd2**7*(1845.*omega1**2 - 266.*omega1*omega2 + 1180.*omega2**2) + kd1**7*kd2**3*(1180.*omega1**2 - 266.*omega1*omega2 + 1845.*omega2**2)) - &
              omega1*(omega1 - omega2)**2*omega2*(kd2**2*(448.*kd1**9 - 3064.*kd1**8*kd2 + 5274.*kd1**7*kd2**2 - 2347.*kd1**6*kd2**3 - 1723.*kd1**5*kd2**4 + &
              3107.*kd1**4*kd2**5 - 2351.*kd1**3*kd2**6 + 726.*kd1**2*kd2**7 + 6.*kd1*kd2**8 - 6.*kd2**9)*omega1**2 + 5.*kd1*kd2*(-kd1 + kd2)*(40.*kd1**8 + &
              354.*kd1**7*kd2 + 651.*kd1**6*kd2**2 - 4170.*kd1**5*kd2**3 + 6292.*kd1**4*kd2**4 - 4170.*kd1**3*kd2**5 + 651.*kd1**2*kd2**6 + 354.*kd1*kd2**7 + &
              40.*kd2**8)*omega1*omega2 + kd1**2*(6.*kd1**9 - 6.*kd1**8*kd2 - 726.*kd1**7*kd2**2 + 2351.*kd1**6*kd2**3 - 3107.*kd1**5*kd2**4 + 1723.*kd1**4*kd2**5 + &
              2347.*kd1**3*kd2**6 - 5274.*kd1**2*kd2**7 + 3064.*kd1*kd2**8 - 448.*kd2**9)*omega2**2)*swd**2) - &
              2251799813685248.*grav*(320.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3 + (2.*kd1*kd2**3*(84.*kd1**3 - 257.*kd1**2*kd2 + 238.*kd1*kd2**2 - 72.*kd2**3)*omega1**4 + &
              kd2**2*(26.*kd1**5 - 503.*kd1**4*kd2 + 1246.*kd1**3*kd2**2 - 1033.*kd1**2*kd2**3 + &
              291.*kd1*kd2**4 + kd2**5)*omega1**3*omega2 + 2.*kd1*(kd1 - kd2)*kd2*(146.*kd1**4 - 379.*kd1**3*kd2 + 496.*kd1**2*kd2**2 - 379.*kd1*kd2**3 + &
              146.*kd2**4)*omega1**2*omega2**2 - kd1**2*(kd1**5 + 291.*kd1**4*kd2 - 1033.*kd1**3*kd2**2 + 1246.*kd1**2*kd2**3 - 503.*kd1*kd2**4 + 26.*kd2**5)*omega1*omega2**3 + &
              2.*kd1**3*kd2*(72.*kd1**3 - 238.*kd1**2*kd2 + 257.*kd1*kd2**2 - 84.*kd2**3)*omega2**4)*swd**2) - &
              35184372088832.*grav*(128.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(16.*kd1**2 - 19.*kd1*kd2 + 16.*kd2**2) + (2.*kd1*kd2**3*(1192.*kd1**5 - 4114.*kd1**4*kd2 + &
              5045.*kd1**3*kd2**2 - 2874.*kd1**2*kd2**3 + 1071.*kd1*kd2**4 - 300.*kd2**5)*omega1**4 + kd2**2*(-473.*kd1**7 - 3756.*kd1**6*kd2 + 16039.*kd1**5*kd2**2 - &
              20077.*kd1**4*kd2**3 + 10893.*kd1**3*kd2**4 - 3870.*kd1**2*kd2**5 + 1158.*kd1*kd2**6 + 6.*kd2**7)*omega1**3*omega2 + 2.*kd1*(kd1 - kd2)*kd2*(574.*kd1**6 - &
              1360.*kd1**5*kd2 + 5458.*kd1**4*kd2**2 - 9625.*kd1**3*kd2**3 + 5458.*kd1**2*kd2**4 - 1360.*kd1*kd2**5 + 574.*kd2**6)*omega1**2*omega2**2 + kd1**2*(-6.*kd1**7 - &
              1158.*kd1**6*kd2 + 3870.*kd1**5*kd2**2 - 10893.*kd1**4*kd2**3 + 20077.*kd1**3*kd2**4 - 16039.*kd1**2*kd2**5 + 3756.*kd1*kd2**6 + 473.*kd2**7)*omega1*omega2**3 + &
              2.*kd1**3*kd2*(300.*kd1**5 - 1071.*kd1**4*kd2 + 2874.*kd1**3*kd2**2 - 5045.*kd1**2*kd2**3 + 4114.*kd1*kd2**4 - 1192.*kd2**5)*omega2**4)*swd**2) - &
              1099511627776.*grav*(64.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(32.*kd1**4 - 65.*kd1**3*kd2 + 114.*kd1**2*kd2**2 - 65.*kd1*kd2**3 + 32.*kd2**4) + &
              (kd1*kd2**3*(1996.*kd1**7 - 7663.*kd1**6*kd2 + 15670.*kd1**5*kd2**2 - 19603.*kd1**4*kd2**3 + 12326.*kd1**3*kd2**4 - 2540.*kd1**2*kd2**5 + 44.*kd1*kd2**6 - &
              248.*kd2**7)*omega1**4 + kd2**2*(-538.*kd1**9 - 1231.*kd1**8*kd2 + 10269.*kd1**7*kd2**2 - 26478.*kd1**6*kd2**3 + 35260.*kd1**5*kd2**4 - 21344.*kd1**4*kd2**5 + &
              3393.*kd1**3*kd2**6 + 177.*kd1**2*kd2**7 + 535.*kd1*kd2**8 - 7.*kd2**9)*omega1**3*omega2 + kd1*(kd1 - kd2)*kd2*(516.*kd1**8 + 90.*kd1**7*kd2 + 7252.*kd1**6*kd2**2 - &
              27643.*kd1**5*kd2**3 + 39422.*kd1**4*kd2**4 - 27643.*kd1**3*kd2**5 + 7252.*kd1**2*kd2**6 + 90.*kd1*kd2**7 + 516.*kd2**8)*omega1**2*omega2**2 + kd1**2*(7.*kd1**9 - &
              535.*kd1**8*kd2 - 177.*kd1**7*kd2**2 - 3393.*kd1**6*kd2**3 + 21344.*kd1**5*kd2**4 - 35260.*kd1**4*kd2**5 + 26478.*kd1**3*kd2**6 - 10269.*kd1**2*kd2**7 + &
              1231.*kd1*kd2**8 + 538.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(248.*kd1**7 - 44.*kd1**6*kd2 + 2540.*kd1**5*kd2**2 - 12326.*kd1**4*kd2**3 + 19603.*kd1**3*kd2**4 - &
              15670.*kd1**2*kd2**5 + 7663.*kd1*kd2**6 - 1996.*kd2**7)*omega2**4)*swd**2) - 64.*kd1**5*(kd1 - kd2)**5*kd2**5*swd*(8.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*((37.*kd1**2 - &
              72.*kd1*kd2 + 32.*kd2**2)*omega1**2 + (-53.*kd1**2 + 82.*kd1*kd2 - 53.*kd2**2)*omega1*omega2 + (32.*kd1**2 - 72.*kd1*kd2 + 37.*kd2**2)*omega2**2) + &
              omega1*(omega1 - omega2)**2*omega2*(12.*kd1**6*omega2*(4.*omega1 + omega2) + 12.*kd2**6*omega1*(omega1 + 4.*omega2) - 4.*kd1**3*kd2**3*(omega1**2 + &
              172.*omega1*omega2 + omega2**2) - 2.*kd1*kd2**5*(164.*omega1**2 + 135.*omega1*omega2 + 3.*omega2**2) + kd1**2*kd2**4*(275.*omega1**2 + 602.*omega1*omega2 + &
              15.*omega2**2) - 2.*kd1**5*kd2*(3.*omega1**2 + 135.*omega1*omega2 + 164.*omega2**2) + kd1**4*kd2**2*(15.*omega1**2 + 602.*omega1*omega2 + 275.*omega2**2))*swd**2) - &
              4096.*kd1**3*(kd1 - kd2)**3*kd2**3*swd*(8.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd2**6*(301.*omega1**2 - 140.*omega1*omega2 + 88.*omega2**2) - &
              3.*kd1*kd2**5*(444.*omega1**2 - 213.*omega1*omega2 + 184.*omega2**2) + kd1**6*(88.*omega1**2 - 140.*omega1*omega2 + 301.*omega2**2) - 8.*kd1**3*kd2**3*(305.*omega1**2 - &
              279.*omega1*omega2 + 305.*omega2**2) - 3.*kd1**5*kd2*(184.*omega1**2 - 213.*omega1*omega2 + 444.*omega2**2) + kd1**2*kd2**4*(2427.*omega1**2 - 1633.*omega1*omega2 + &
              1504.*omega2**2) + kd1**4*kd2**2*(1504.*omega1**2 - 1633.*omega1*omega2 + 2427.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(5.*kd1**10*omega2*(11.*omega1 + &
              omega2) + 5.*kd2**10*omega1*(omega1 + 11.*omega2) + 3.*kd1**7*kd2**3*(8.*omega1**2 - 4117.*omega1*omega2 - 170.*omega2**2) + kd1**9*kd2*(6.*omega1**2 - &
              516.*omega1*omega2 - 31.*omega2**2) + kd1**2*kd2**8*(298.*omega1**2 + 3041.*omega1*omega2 - 27.*omega2**2) - 3.*kd1**3*kd2**7*(170.*omega1**2 + 4117.*omega1*omega2 - &
              8.*omega2**2) + kd1*kd2**9*(-31.*omega1**2 - 516.*omega1*omega2 + 6.*omega2**2) - 6.*kd1**4*kd2**6*(266.*omega1**2 - 4712.*omega1*omega2 + 135.*omega2**2) - &
              6.*kd1**6*kd2**4*(135.*omega1**2 - 4712.*omega1*omega2 + 266.*omega2**2) + kd1**8*kd2**2*(-27.*omega1**2 + 3041.*omega1*omega2 + 298.*omega2**2) + &
              4.*kd1**5*kd2**5*(641.*omega1**2 - 9212.*omega1*omega2 + 641.*omega2**2))*swd**2) + 4398046511104.*swd*(8.*grav**2*kd1*(kd1 - kd2)*kd2*(82.*kd2**8*omega1**2 + &
              82.*kd1**8*omega2**2 - 2.*kd1*kd2**7*omega1*(139.*omega1 + 23.*omega2) - 2.*kd1**7*kd2*omega2*(23.*omega1 + 139.*omega2) + kd1**2*kd2**6*(1580.*omega1**2 - &
              841.*omega1*omega2 + 931.*omega2**2) - 2.*kd1**3*kd2**5*(1927.*omega1**2 - 1360.*omega1*omega2 + 1511.*omega2**2) + kd1**6*kd2**2*(931.*omega1**2 - &
              841.*omega1*omega2 + 1580.*omega2**2) - 2.*kd1**5*kd2**3*(1511.*omega1**2 - 1360.*omega1*omega2 + 1927.*omega2**2) + kd1**4*kd2**4*(4541.*omega1**2 - &
              3626.*omega1*omega2 + 4541.*omega2**2)) - omega1*(omega1 - omega2)**2*omega2*(14.*kd2**9*omega1**2 - 14.*kd1**9*omega2**2 + 2.*kd1**8*kd2*omega2*(-326.*omega1 + &
              17.*omega2) + 2.*kd1*kd2**8*omega1*(-17.*omega1 + 326.*omega2) - 2.*kd1**2*kd2**7*(270.*omega1**2 + 908.*omega1*omega2 + 269.*omega2**2) + &
              2.*kd1**7*kd2**2*(269.*omega1**2 + 908.*omega1*omega2 + 270.*omega2**2) - kd1**6*kd2**3*(1707.*omega1**2 + 7830.*omega1*omega2 + 1649.*omega2**2) + &
              kd1**3*kd2**6*(1649.*omega1**2 + 7830.*omega1*omega2 + 1707.*omega2**2) + kd1**5*kd2**4*(2121.*omega1**2 + 19440.*omega1*omega2 + 2071.*omega2**2) - &
              kd1**4*kd2**5*(2071.*omega1**2 + 19440.*omega1*omega2 + 2121.*omega2**2))*swd**2) - 262144.*kd1*(kd1 - kd2)*kd2*swd*(8.*grav**2*kd1**2*(kd1 - &
              kd2)**2*kd2**2*(kd1**7*kd2**3*(-7660.*omega1**2 + 6921.*omega1*omega2 - 12524.*omega2**2) + kd1**3*kd2**7*(-12524.*omega1**2 + 6921.*omega1*omega2 - &
              7660.*omega2**2) + kd2**10*(170.*omega1**2 - 47.*omega1*omega2 + 35.*omega2**2) - 4.*kd1*kd2**9*(347.*omega1**2 - 95.*omega1*omega2 + 67.*omega2**2) + &
              kd1**10*(35.*omega1**2 - 47.*omega1*omega2 + 170.*omega2**2) - 4.*kd1**9*kd2*(67.*omega1**2 - 95.*omega1*omega2 + 347.*omega2**2) + kd1**2*kd2**8*(4955.*omega1**2 - &
              1723.*omega1*omega2 + 1580.*omega2**2) + kd1**8*kd2**2*(1580.*omega1**2 - 1723.*omega1*omega2 + 4955.*omega2**2) - 2.*kd1**5*kd2**5*(13873.*omega1**2 - &
              11288.*omega1*omega2 + 13873.*omega2**2) + kd1**4*kd2**6*(23139.*omega1**2 - 16821.*omega1*omega2 + 19706.*omega2**2) + kd1**6*kd2**4*(19706.*omega1**2 - &
              16821.*omega1*omega2 + 23139.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(4.*kd1**14*omega1*omega2 + 4.*kd2**14*omega1*omega2 + kd1**5*kd2**9*(1395.*omega1**2 - &
              125589.*omega1*omega2 - 4016.*omega2**2) + 2.*kd1**8*kd2**6*(2421.*omega1**2 + 110385.*omega1*omega2 - 3247.*omega2**2) - 2.*kd1*kd2**13*(35.*omega1**2 + &
              27.*omega1*omega2 + omega2**2) - kd1**3*kd2**11*(1369.*omega1**2 + 10424.*omega1*omega2 + 12.*omega2**2) + kd1**2*kd2**12*(447.*omega1**2 + 1074.*omega1*omega2 + &
              13.*omega2**2) - 2.*kd1**13*kd2*(omega1**2 + 27.*omega1*omega2 + 35.*omega2**2) + kd1**12*kd2**2*(13.*omega1**2 + 1074.*omega1*omega2 + 447.*omega2**2) + &
              kd1**4*kd2**10*(2029.*omega1**2 + 46930.*omega1*omega2 + 775.*omega2**2) + 2.*kd1**7*kd2**7*(1227.*omega1**2 - 132703.*omega1*omega2 + 1227.*omega2**2) - &
              kd1**11*kd2**3*(12.*omega1**2 + 10424.*omega1*omega2 + 1369.*omega2**2) + kd1**9*kd2**5*(-4016.*omega1**2 - 125589.*omega1*omega2 + 1395.*omega2**2) + &
              kd1**10*kd2**4*(775.*omega1**2 + 46930.*omega1*omega2 + 2029.*omega2**2) + 2.*kd1**6*kd2**8*(-3247.*omega1**2 + 110385.*omega1*omega2 + &
              2421.*omega2**2))*swd**2) + 16777216.*swd*(-8.*grav**2*kd1*(kd1 - kd2)**2*kd2*(-9.*kd2**13*omega1**2 + 9.*kd1**13*omega2**2 + kd1*kd2**12*omega1*(71.*omega1 + &
              9.*omega2) - kd1**12*kd2*omega2*(9.*omega1 + 71.*omega2) + kd1**8*kd2**5*(-15047.*omega1**2 + 10913.*omega1*omega2 - 29631.*omega2**2) + kd1**4*kd2**9*(-14625.*omega1**2 + &
              4217.*omega1*omega2 - 6065.*omega2**2) - 14.*kd1**2*kd2**11*(21.*omega1**2 - 5.*omega1*omega2 + 9.*omega2**2) + 14.*kd1**11*kd2**2*(9.*omega1**2 - 5.*omega1*omega2 + &
              21.*omega2**2) + 2.*kd1**3*kd2**10*(1862.*omega1**2 - 715.*omega1*omega2 + 871.*omega2**2) - 2.*kd1**10*kd2**3*(871.*omega1**2 - 715.*omega1*omega2 + &
              1862.*omega2**2) - 2.*kd1**6*kd2**7*(18425.*omega1**2 - 10516.*omega1*omega2 + 14474.*omega2**2) + kd1**9*kd2**4*(6065.*omega1**2 - 4217.*omega1*omega2 + &
              14625.*omega2**2) + kd1**5*kd2**8*(29631.*omega1**2 - 10913.*omega1*omega2 + 15047.*omega2**2) + 2.*kd1**7*kd2**6*(14474.*omega1**2 - 10516.*omega1*omega2 + &
              18425.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(-(kd2**15*omega1**2) + kd1**15*omega2**2 - kd1**14*kd2*omega2*(68.*omega1 + 5.*omega2) + &
              kd1*kd2**14*omega1*(5.*omega1 + 68.*omega2) + kd1**9*kd2**6*(8791.*omega1**2 + 230693.*omega1*omega2 - 13087.*omega2**2) + kd1**6*kd2**9*(13087.*omega1**2 - &
              230693.*omega1*omega2 - 8791.*omega2**2) + kd1**11*kd2**4*(-1470.*omega1**2 + 27574.*omega1*omega2 - 4229.*omega2**2) - 3.*kd1**8*kd2**7*(4111.*omega1**2 + &
              107400.*omega1*omega2 - 428.*omega2**2) + kd1**5*kd2**10*(-10605.*omega1**2 + 106859.*omega1*omega2 - 197.*omega2**2) + kd1**13*kd2**2*(8.*omega1**2 + &
              757.*omega1*omega2 - 99.*omega2**2) + kd1**2*kd2**13*(99.*omega1**2 - 757.*omega1*omega2 - 8.*omega2**2) - 2.*kd1**3*kd2**12*(485.*omega1**2 - 1998.*omega1*omega2 + &
              116.*omega2**2) + 2.*kd1**12*kd2**3*(116.*omega1**2 - 1998.*omega1*omega2 + 485.*omega2**2) + kd1**4*kd2**11*(4229.*omega1**2 - 27574.*omega1*omega2 + &
              1470.*omega2**2) + 3.*kd1**7*kd2**8*(-428.*omega1**2 + 107400.*omega1*omega2 + 4111.*omega2**2) + kd1**10*kd2**5*(197.*omega1**2 - 106859.*omega1*omega2 + &
              10605.*omega2**2))*swd**2) - 1073741824.*swd*(8.*grav**2*kd1*(kd1 - kd2)*kd2*(24.*kd2**12*omega1**2 + kd1**11*kd2*(19.*omega1 - 378.*omega2)*omega2 + &
              24.*kd1**12*omega2**2 + kd1*kd2**11*omega1*(-378.*omega1 + 19.*omega2) + kd1**5*kd2**7*(-8314.*omega1**2 + 3005.*omega1*omega2 - 11014.*omega2**2) + &
              kd1**7*kd2**5*(-11014.*omega1**2 + 3005.*omega1*omega2 - 8314.*omega2**2) + kd1**9*kd2**3*(-1656.*omega1**2 + 2303.*omega1*omega2 - 1946.*omega2**2) + &
              kd1**3*kd2**9*(-1946.*omega1**2 + 2303.*omega1*omega2 - 1656.*omega2**2) + kd1**2*kd2**10*(739.*omega1**2 - 66.*omega1*omega2 - 219.*omega2**2) + &
              kd1**10*kd2**2*(-219.*omega1**2 - 66.*omega1*omega2 + 739.*omega2**2) + 2.*kd1**8*kd2**4*(3534.*omega1**2 - 2520.*omega1*omega2 + 2557.*omega2**2) + &
              2.*kd1**4*kd2**8*(2557.*omega1**2 - 2520.*omega1*omega2 + 3534.*omega2**2) + kd1**6*kd2**6*(10601.*omega1**2 - 480.*omega1*omega2 + 10601.*omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(-(kd2**13*omega1**2) + kd1*kd2**12*omega1*(11.*omega1 - 160.*omega2) + kd1**12*kd2*(160.*omega1 - 11.*omega2)*omega2 + &
              kd1**13*omega2**2 + kd1**4*kd2**9*(1139.*omega1**2 + 6478.*omega1*omega2 - 4041.*omega2**2) + kd1**9*kd2**4*(4041.*omega1**2 - 6478.*omega1*omega2 - &
              1139.*omega2**2) + kd1**10*kd2**3*(-1365.*omega1**2 + 3679.*omega1*omega2 - 594.*omega2**2) + kd1**2*kd2**11*(-283.*omega1**2 + 2351.*omega1*omega2 - &
              46.*omega2**2) + kd1**11*kd2**2*(46.*omega1**2 - 2351.*omega1*omega2 + 283.*omega2**2) - kd1**5*kd2**8*(7195.*omega1**2 + 26790.*omega1*omega2 + &
              433.*omega2**2) + kd1**3*kd2**10*(594.*omega1**2 - 3679.*omega1*omega2 + 1365.*omega2**2) + kd1**8*kd2**5*(433.*omega1**2 + 26790.*omega1*omega2 + &
              7195.*omega2**2) + kd1**6*kd2**7*(14617.*omega1**2 + 54009.*omega1*omega2 + 12092.*omega2**2) - kd1**7*kd2**6*(12092.*omega1**2 + 54009.*omega1*omega2 + &
              14617.*omega2**2))*swd**2) + 512.*grav*kd1**5*(kd1 - kd2)**5*kd2**5*(640.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2 + (kd1**5*kd2*(-432.*omega1**4 + &
              615.*omega1**3*omega2 - 2772.*omega1**2*omega2**2 + 3134.*omega1*omega2**3 - 2182.*omega2**4) + kd1*kd2**5*(-2182.*omega1**4 + 3134.*omega1**3*omega2 - &
              2772.*omega1**2*omega2**2 + 615.*omega1*omega2**3 - 432.*omega2**4) + 2.*kd2**6*(252.*omega1**4 - 409.*omega1**3*omega2 + 273.*omega1**2*omega2**2 - &
              59.*omega1*omega2**3 + 33.*omega2**4) + 2.*kd1**6*(33.*omega1**4 - 59.*omega1**3*omega2 + 273.*omega1**2*omega2**2 - 409.*omega1*omega2**3 + 252.*omega2**4) + &
              kd1**2*kd2**4*(3730.*omega1**4 - 6785.*omega1**3*omega2 + 7304.*omega1**2*omega2**2 - 2917.*omega1*omega2**3 + 1472.*omega2**4) - &
              2.*kd1**3*kd2**3*(1595.*omega1**4 - 3291.*omega1**3*omega2 + 4739.*omega1**2*omega2**2 - 3291.*omega1*omega2**3 + 1595.*omega2**4) + &
              kd1**4*kd2**2*(1472.*omega1**4 - 2917.*omega1**3*omega2 + 7304.*omega1**2*omega2**2 - 6785.*omega1*omega2**3 + 3730.*omega2**4))*swd**2) - &
              8589934592.*grav*(128.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(19.*kd1**6 + 5.*kd1**5*kd2 + 13.*kd1**4*kd2**2 - 43.*kd1**3*kd2**3 + 13.*kd1**2*kd2**4 + &
              5.*kd1*kd2**5 + 19.*kd2**6) + (11.*kd2**13*omega1**3*omega2 - 11.*kd1**13*omega1*omega2**3 + kd1**12*kd2*omega2**2*(-96.*omega1**2 + 185.*omega1*omega2 - &
              64.*omega2**2) + kd1*kd2**12*omega1**2*(64.*omega1**2 - 185.*omega1*omega2 + 96.*omega2**2) + kd1**2*kd2**11*omega1*(-2420.*omega1**3 + 4821.*omega1**2*omega2 - &
              5936.*omega1*omega2**2 + 609.*omega2**3) + kd1**11*kd2**2*omega2*(-609.*omega1**3 + 5936.*omega1**2*omega2 - 4821.*omega1*omega2**2 + 2420.*omega2**3) + &
              kd1**7*kd2**6*(-4566.*omega1**4 + 19533.*omega1**3*omega2 - 152462.*omega1**2*omega2**2 + 67537.*omega1*omega2**3 - 24008.*omega2**4) + &
              kd1**3*kd2**10*(1844.*omega1**4 - 4111.*omega1**3*omega2 + 12316.*omega1**2*omega2**2 - 1445.*omega1*omega2**3 - 2272.*omega2**4) + &
              kd1**10*kd2**3*(2272.*omega1**4 + 1445.*omega1**3*omega2 - 12316.*omega1**2*omega2**2 + 4111.*omega1*omega2**3 - 1844.*omega2**4) + &
              kd1**4*kd2**9*(10658.*omega1**4 - 21885.*omega1**3*omega2 + 7342.*omega1**2*omega2**2 + 12143.*omega1*omega2**3 + 1464.*omega2**4) + &
              kd1**5*kd2**8*(-25998.*omega1**4 + 61553.*omega1**3*omega2 - 74332.*omega1**2*omega2**2 - 19091.*omega1*omega2**3 + 4390.*omega2**4) + &
              kd1**6*kd2**7*(24008.*omega1**4 - 67537.*omega1**3*omega2 + 152462.*omega1**2*omega2**2 - 19533.*omega1*omega2**3 + 4566.*omega2**4) - &
              kd1**9*kd2**4*(1464.*omega1**4 + 12143.*omega1**3*omega2 + 7342.*omega1**2*omega2**2 - 21885.*omega1*omega2**3 + 10658.*omega2**4) + &
              kd1**8*kd2**5*(-4390.*omega1**4 + 19091.*omega1**3*omega2 + 74332.*omega1**2*omega2**2 - 61553.*omega1*omega2**3 + 25998.*omega2**4))*swd**2) + &
              32768.*grav*kd1**3*(kd1 - kd2)**3*kd2**3*(64.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(32.*kd1**4 - 96.*kd1**3*kd2 + 95.*kd1**2*kd2**2 - 96.*kd1*kd2**3 + &
              32.*kd2**4) + (kd1**5*kd2**5*(-59542.*omega1**4 + 130733.*omega1**3*omega2 - 237196.*omega1**2*omega2**2 + 130733.*omega1*omega2**3 - 59542.*omega2**4) + &
              kd1**9*kd2*(-246.*omega1**4 + 466.*omega1**3*omega2 - 5332.*omega1**2*omega2**2 + 6625.*omega1*omega2**3 - 3546.*omega2**4) + kd1*kd2**9*(-3546.*omega1**4 + &
              6625.*omega1**3*omega2 - 5332.*omega1**2*omega2**2 + 466.*omega1*omega2**3 - 246.*omega2**4) + 4.*kd2**10*(116.*omega1**4 - 232.*omega1**3*omega2 + &
              167.*omega1**2*omega2**2 - 12.*omega1*omega2**3 + 6.*omega2**4) + 4.*kd1**10*(6.*omega1**4 - 12.*omega1**3*omega2 + 167.*omega1**2*omega2**2 - &
              232.*omega1*omega2**3 + 116.*omega2**4) + 2.*kd1**2*kd2**8*(6685.*omega1**4 - 12259.*omega1**3*omega2 + 11794.*omega1**2*omega2**2 - 1948.*omega1*omega2**3 + &
              1121.*omega2**4) - 2.*kd1**3*kd2**7*(16648.*omega1**4 - 32338.*omega1**3*omega2 + 39848.*omega1**2*omega2**2 - 12522.*omega1*omega2**3 + 6681.*omega2**4) + &
              2.*kd1**8*kd2**2*(1121.*omega1**4 - 1948.*omega1**3*omega2 + 11794.*omega1**2*omega2**2 - 12259.*omega1*omega2**3 + 6685.*omega2**4) - &
              2.*kd1**7*kd2**3*(6681.*omega1**4 - 12522.*omega1**3*omega2 + 39848.*omega1**2*omega2**2 - 32338.*omega1*omega2**3 + 16648.*omega2**4) + &
              kd1**4*kd2**6*(55726.*omega1**4 - 118435.*omega1**3*omega2 + 179446.*omega1**2*omega2**2 - 79771.*omega1*omega2**3 + 38142.*omega2**4) + &
              kd1**6*kd2**4*(38142.*omega1**4 - 79771.*omega1**3*omega2 + 179446.*omega1**2*omega2**2 - 118435.*omega1*omega2**3 + 55726.*omega2**4))*swd**2) + &
              134217728.*grav*(64.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(5.*kd1**8 - 170.*kd1**7*kd2 + 456.*kd1**6*kd2**2 - 662.*kd1**5*kd2**3 + 739.*kd1**4*kd2**4 - &
              662.*kd1**3*kd2**5 + 456.*kd1**2*kd2**6 - 170.*kd1*kd2**7 + 5.*kd2**8) + (3.*kd2**15*omega1**3*omega2 - 3.*kd1**15*omega1*omega2**3 + &
              kd1*kd2**14*omega1**2*(-128.*omega1**2 + 241.*omega1*omega2 - 300.*omega2**2) + kd1**14*kd2*omega2**2*(300.*omega1**2 - 241.*omega1*omega2 + &
              128.*omega2**2) + kd1**2*kd2**13*omega1*(1508.*omega1**3 - 3303.*omega1**2*omega2 + 3614.*omega1*omega2**2 + 58.*omega2**3) - &
              kd1**13*kd2**2*omega2*(58.*omega1**3 + 3614.*omega1**2*omega2 - 3303.*omega1*omega2**2 + 1508.*omega2**3) + kd1**9*kd2**6*(-80420.*omega1**4 + &
              129537.*omega1**3*omega2 - 506278.*omega1**2*omega2**2 + 292383.*omega1*omega2**3 - 159804.*omega2**4) + kd1**5*kd2**10*(-81942.*omega1**4 + &
              166165.*omega1**3*omega2 - 234512.*omega1**2*omega2**2 + 75199.*omega1*omega2**3 - 36956.*omega2**4) + kd1**3*kd2**12*(-3786.*omega1**4 + &
              10136.*omega1**3*omega2 - 12366.*omega1**2*omega2**2 + 1621.*omega1*omega2**3 - 296.*omega2**4) + kd1**12*kd2**3*(296.*omega1**4 - 1621.*omega1**3*omega2 + &
              12366.*omega1**2*omega2**2 - 10136.*omega1*omega2**3 + 3786.*omega2**4) + 2.*kd1**4*kd2**11*(10305.*omega1**4 - 24161.*omega1**3*omega2 + &
              31434.*omega1**2*omega2**2 - 11707.*omega1*omega2**3 + 4465.*omega2**4) - 2.*kd1**11*kd2**4*(4465.*omega1**4 - 11707.*omega1**3*omega2 + &
              31434.*omega1**2*omega2**2 - 24161.*omega1*omega2**3 + 10305.*omega2**4) - 2.*kd1**7*kd2**8*(88804.*omega1**4 - 144871.*omega1**3*omega2 + &
              355553.*omega1**2*omega2**2 - 99174.*omega1*omega2**3 + 66819.*omega2**4) + kd1**6*kd2**9*(159804.*omega1**4 - 292383.*omega1**3*omega2 + &
              506278.*omega1**2*omega2**2 - 129537.*omega1*omega2**3 + 80420.*omega2**4) + kd1**10*kd2**5*(36956.*omega1**4 - 75199.*omega1**3*omega2 + &
              234512.*omega1**2*omega2**2 - 166165.*omega1*omega2**3 + 81942.*omega2**4) + 2.*kd1**8*kd2**7*(66819.*omega1**4 - 99174.*omega1**3*omega2 + &
              355553.*omega1**2*omega2**2 - 144871.*omega1*omega2**3 + 88804.*omega2**4))*swd**2) + 2097152.*grav*kd1*(kd1 - kd2)*kd2*(128.*grav**2*kd1**2*(kd1 - &
              kd2)**4*kd2**2*(3.*kd1**8 - 20.*kd1**7*kd2 + 48.*kd1**6*kd2**2 - 229.*kd1**5*kd2**3 + 397.*kd1**4*kd2**4 - 229.*kd1**3*kd2**5 + 48.*kd1**2*kd2**6 - &
              20.*kd1*kd2**7 + 3.*kd2**8) + (12.*kd1**14*omega2**2*(3.*omega1**2 - 4.*omega1*omega2 + 2.*omega2**2) + 12.*kd2**14*omega1**2*(2.*omega1**2 - &
              4.*omega1*omega2 + 3.*omega2**2) - kd1*kd2**13*omega1*(290.*omega1**3 - 476.*omega1**2*omega2 + 304.*omega1*omega2**2 + 25.*omega2**3) - &
              kd1**13*kd2*omega2*(25.*omega1**3 + 304.*omega1**2*omega2 - 476.*omega1*omega2**2 + 290.*omega2**3) + kd1**7*kd2**7*(-240292.*omega1**4 + &
              494923.*omega1**3*omega2 - 966418.*omega1**2*omega2**2 + 494923.*omega1*omega2**3 - 240292.*omega2**4) + kd1**11*kd2**3*(-3692.*omega1**4 + &
              7752.*omega1**3*omega2 - 46292.*omega1**2*omega2**2 + 36947.*omega1*omega2**3 - 18262.*omega2**4) + kd1**3*kd2**11*(-18262.*omega1**4 + &
              36947.*omega1**3*omega2 - 46292.*omega1**2*omega2**2 + 7752.*omega1*omega2**3 - 3692.*omega2**4) + kd1**2*kd2**12*(2398.*omega1**4 - 4299.*omega1**3*omega2 + &
              4592.*omega1**2*omega2**2 - 563.*omega1*omega2**3 + 304.*omega2**4) + kd1**12*kd2**2*(304.*omega1**4 - 563.*omega1**3*omega2 + 4592.*omega1**2*omega2**2 - &
              4299.*omega1*omega2**3 + 2398.*omega2**4) + kd1**4*kd2**10*(70076.*omega1**4 - 145863.*omega1**3*omega2 + 194628.*omega1**2*omega2**2 - 34821.*omega1*omega2**3 + &
              18326.*omega2**4) - 2.*kd1**5*kd2**9*(79465.*omega1**4 - 167361.*omega1**3*omega2 + 240694.*omega1**2*omega2**2 - 60812.*omega1*omega2**3 + 32749.*omega2**4) + &
              kd1**10*kd2**4*(18326.*omega1**4 - 34821.*omega1**3*omega2 + 194628.*omega1**2*omega2**2 - 145863.*omega1*omega2**3 + 70076.*omega2**4) - &
              2.*kd1**9*kd2**5*(32749.*omega1**4 - 60812.*omega1**3*omega2 + 240694.*omega1**2*omega2**2 - 167361.*omega1*omega2**3 + 79465.*omega2**4) + &
              kd1**6*kd2**8*(238420.*omega1**4 - 502403.*omega1**3*omega2 + 811942.*omega1**2*omega2**2 - 308427.*omega1*omega2**3 + 157416.*omega2**4) + &
              kd1**8*kd2**6*(157416.*omega1**4 - 308427.*omega1**3*omega2 + 811942.*omega1**2*omega2**2 - 502403.*omega1*omega2**3 + 238420.*omega2**4))*swd**2)))/ &
              ((-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*(omega1 - omega2)*omega2*swd**3*(-64.*grav*(64. + (kd1 - kd2)**2)*(4096. + &
              (384. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2 + (16777216. + (7340032. + (286720. + (1792. + (kd1 - kd2)**2)*(kd1 - &
              kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function velsb24
!
real function velsp24()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic velocity of 2nd layer for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp24')
    !
    velsp24 = (4.*(144115188075855872.*grav*(kd1 + kd2)**3*omega1*omega2*(kd2**2*omega1**2 + 8.*kd1*kd2*omega1*omega2 + kd1**2*omega2**2)*swd**2 - &
              16.*grav*kd1**7*kd2**7*(kd1 + kd2)**7*((11.*kd1**2 + 31.*kd1*kd2 + 16.*kd2**2)*omega1**4 + (21.*kd1**2 + 6.*kd1*kd2 + 55.*kd2**2)*omega1**3*omega2 + &
              (55.*kd1**2 - 38.*kd1*kd2 + 55.*kd2**2)*omega1**2*omega2**2 + (55.*kd1**2 + 6.*kd1*kd2 + 21.*kd2**2)*omega1*omega2**3 + (16.*kd1**2 + 31.*kd1*kd2 + &
              11.*kd2**2)*omega2**4)*swd**2 + kd1**7*kd2**7*(kd1 + kd2)**7*omega1*omega2*(omega1 + omega2)**2*(kd1**2*omega2*(5.*omega1 + omega2) + &
              kd2**2*omega1*(omega1 + 5.*omega2) + 2.*kd1*kd2*(omega1**2 + 6.*omega1*omega2 + omega2**2))*swd**3 - 18014398509481984.*swd*(8.*grav**2*kd1*kd2*(kd1 + &
              kd2)*(6.*kd2**4*omega1**2 + 6.*kd1**4*omega2**2 + kd1*kd2**3*omega1*(12.*omega1 + omega2) + kd1**3*kd2*omega2*(omega1 + 12.*omega2) + &
              5.*kd1**2*kd2**2*(omega1**2 + omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(kd2**2*(6.*kd1**3 + 17.*kd1**2*kd2 + 9.*kd1*kd2**2 + 3.*kd2**3)*omega1**2 + &
              3.*kd1*kd2*(kd1 + kd2)*(8.*kd1**2 + 19.*kd1*kd2 + 8.*kd2**2)*omega1*omega2 + kd1**2*(3.*kd1**3 + 9.*kd1**2*kd2 + 17.*kd1*kd2**2 + &
              6.*kd2**3)*omega2**2)*swd**2) - 281474976710656.*swd*(8.*grav**2*kd1*kd2*(kd1 + kd2)*(41.*kd2**6*omega1**2 + kd1*kd2**5*omega1*(118.*omega1 - 11.*omega2) + &
              41.*kd1**6*omega2**2 + kd1**5*kd2*omega2*(-11.*omega1 + 118.*omega2) + 12.*kd1**3*kd2**3*(43.*omega1**2 + 40.*omega1*omega2 + 43.*omega2**2) + &
              kd1**2*kd2**4*(353.*omega1**2 + 215.*omega1*omega2 + 226.*omega2**2) + kd1**4*kd2**2*(226.*omega1**2 + 215.*omega1*omega2 + 353.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(kd2**2*(8.*kd1**5 - 192.*kd1**4*kd2 - 334.*kd1**3*kd2**2 - 109.*kd1**2*kd2**3 + 31.*kd1*kd2**4 + 11.*kd2**5)*omega1**2 - &
              kd1*kd2*(kd1 + kd2)*(380.*kd1**4 + 923.*kd1**3*kd2 + 1099.*kd1**2*kd2**2 + 923.*kd1*kd2**3 + 380.*kd2**4)*omega1*omega2 + kd1**2*(11.*kd1**5 + 31.*kd1**4*kd2 - &
              109.*kd1**3*kd2**2 - 334.*kd1**2*kd2**3 - 192.*kd1*kd2**4 + 8.*kd2**5)*omega2**2)*swd**2) - 4398046511104.*swd*(8.*grav**2*kd1*kd2*(kd1 + &
              kd2)*(82.*kd2**8*omega1**2 + 2.*kd1*kd2**7*omega1*(139.*omega1 - 23.*omega2) + 82.*kd1**8*omega2**2 + 2.*kd1**7*kd2*omega2*(-23.*omega1 + 139.*omega2) + &
              kd1**2*kd2**6*(1580.*omega1**2 + 841.*omega1*omega2 + 931.*omega2**2) + 2.*kd1**3*kd2**5*(1927.*omega1**2 + 1360.*omega1*omega2 + 1511.*omega2**2) + &
              kd1**6*kd2**2*(931.*omega1**2 + 841.*omega1*omega2 + 1580.*omega2**2) + 2.*kd1**5*kd2**3*(1511.*omega1**2 + 1360.*omega1*omega2 + 1927.*omega2**2) + &
              kd1**4*kd2**4*(4541.*omega1**2 + 3626.*omega1*omega2 + 4541.*omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(kd2**2*(538.*kd1**7 + 1707.*kd1**6*kd2 + &
              2121.*kd1**5*kd2**2 + 2071.*kd1**4*kd2**3 + 1649.*kd1**3*kd2**4 + 540.*kd1**2*kd2**5 - 34.*kd1*kd2**6 - 14.*kd2**7)*omega1**2 - 2.*kd1*kd2*(kd1 + kd2)*(326.*kd1**6 + &
              582.*kd1**5*kd2 + 3333.*kd1**4*kd2**2 + 6387.*kd1**3*kd2**3 + 3333.*kd1**2*kd2**4 + 582.*kd1*kd2**5 + 326.*kd2**6)*omega1*omega2 + kd1**2*(-14.*kd1**7 - &
              34.*kd1**6*kd2 + 540.*kd1**5*kd2**2 + 1649.*kd1**4*kd2**3 + 2071.*kd1**3*kd2**4 + 2121.*kd1**2*kd2**5 + 1707.*kd1*kd2**6 + &
              538.*kd2**7)*omega2**2)*swd**2) - 68719476736.*swd*(16.*grav**2*kd1*kd2*(kd1 + kd2)*(16.*kd2**10*omega1**2 + 16.*kd1**10*omega2**2 - kd1*kd2**9*omega1*(71.*omega1 + &
              31.*omega2) - kd1**9*kd2*omega2*(31.*omega1 + 71.*omega2) + kd1**8*kd2**2*(414.*omega1**2 + 199.*omega1*omega2 + 340.*omega2**2) + kd1**2*kd2**8*(340.*omega1**2 + &
              199.*omega1*omega2 + 414.*omega2**2) + 5.*kd1**5*kd2**5*(917.*omega1**2 + 374.*omega1*omega2 + 917.*omega2**2) + kd1**3*kd2**7*(1845.*omega1**2 + &
              266.*omega1*omega2 + 1180.*omega2**2) + kd1**7*kd2**3*(1180.*omega1**2 + 266.*omega1*omega2 + 1845.*omega2**2) + kd1**4*kd2**6*(3985.*omega1**2 + &
              966.*omega1*omega2 + 2779.*omega2**2) + kd1**6*kd2**4*(2779.*omega1**2 + 966.*omega1*omega2 + 3985.*omega2**2)) - omega1*omega2*(omega1 + &
              omega2)**2*(kd2**2*(448.*kd1**9 + 3064.*kd1**8*kd2 + 5274.*kd1**7*kd2**2 + 2347.*kd1**6*kd2**3 - 1723.*kd1**5*kd2**4 - 3107.*kd1**4*kd2**5 - 2351.*kd1**3*kd2**6 - &
              726.*kd1**2*kd2**7 + 6.*kd1*kd2**8 + 6.*kd2**9)*omega1**2 - 5.*kd1*kd2*(kd1 + kd2)*(40.*kd1**8 - 354.*kd1**7*kd2 + 651.*kd1**6*kd2**2 + 4170.*kd1**5*kd2**3 + &
              6292.*kd1**4*kd2**4 + 4170.*kd1**3*kd2**5 + 651.*kd1**2*kd2**6 - 354.*kd1*kd2**7 + 40.*kd2**8)*omega1*omega2 + kd1**2*(6.*kd1**9 + 6.*kd1**8*kd2 - &
              726.*kd1**7*kd2**2 - 2351.*kd1**6*kd2**3 - 3107.*kd1**5*kd2**4 - 1723.*kd1**4*kd2**5 + 2347.*kd1**3*kd2**6 + 5274.*kd1**2*kd2**7 + 3064.*kd1*kd2**8 + &
              448.*kd2**9)*omega2**2)*swd**2) + 2251799813685248.*grav*(320.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3 + (2.*kd1*kd2**3*(84.*kd1**3 + 257.*kd1**2*kd2 + &
              238.*kd1*kd2**2 + 72.*kd2**3)*omega1**4 + kd2**2*(26.*kd1**5 + 503.*kd1**4*kd2 + 1246.*kd1**3*kd2**2 + 1033.*kd1**2*kd2**3 + 291.*kd1*kd2**4 - &
              kd2**5)*omega1**3*omega2 + 2.*kd1*kd2*(kd1 + kd2)*(146.*kd1**4 + 379.*kd1**3*kd2 + 496.*kd1**2*kd2**2 + 379.*kd1*kd2**3 + 146.*kd2**4)*omega1**2*omega2**2 + &
              kd1**2*(-kd1**5 + 291.*kd1**4*kd2 + 1033.*kd1**3*kd2**2 + 1246.*kd1**2*kd2**3 + 503.*kd1*kd2**4 + 26.*kd2**5)*omega1*omega2**3 + 2.*kd1**3*kd2*(72.*kd1**3 + &
              238.*kd1**2*kd2 + 257.*kd1*kd2**2 + 84.*kd2**3)*omega2**4)*swd**2) + 35184372088832.*grav*(128.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(16.*kd1**2 + &
              19.*kd1*kd2 + 16.*kd2**2) + (2.*kd1*kd2**3*(1192.*kd1**5 + 4114.*kd1**4*kd2 + 5045.*kd1**3*kd2**2 + 2874.*kd1**2*kd2**3 + 1071.*kd1*kd2**4 + &
              300.*kd2**5)*omega1**4 + kd2**2*(-473.*kd1**7 + 3756.*kd1**6*kd2 + 16039.*kd1**5*kd2**2 + 20077.*kd1**4*kd2**3 + 10893.*kd1**3*kd2**4 + 3870.*kd1**2*kd2**5 + &
              1158.*kd1*kd2**6 - 6.*kd2**7)*omega1**3*omega2 + 2.*kd1*kd2*(kd1 + kd2)*(574.*kd1**6 + 1360.*kd1**5*kd2 + 5458.*kd1**4*kd2**2 + 9625.*kd1**3*kd2**3 + &
              5458.*kd1**2*kd2**4 + 1360.*kd1*kd2**5 + 574.*kd2**6)*omega1**2*omega2**2 + kd1**2*(-6.*kd1**7 + 1158.*kd1**6*kd2 + 3870.*kd1**5*kd2**2 + 10893.*kd1**4*kd2**3 + &
              20077.*kd1**3*kd2**4 + 16039.*kd1**2*kd2**5 + 3756.*kd1*kd2**6 - 473.*kd2**7)*omega1*omega2**3 + 2.*kd1**3*kd2*(300.*kd1**5 + 1071.*kd1**4*kd2 + 2874.*kd1**3*kd2**2 + &
              5045.*kd1**2*kd2**3 + 4114.*kd1*kd2**4 + 1192.*kd2**5)*omega2**4)*swd**2) + 1099511627776.*grav*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(32.*kd1**4 + &
              65.*kd1**3*kd2 + 114.*kd1**2*kd2**2 + 65.*kd1*kd2**3 + 32.*kd2**4) + (kd1*kd2**3*(1996.*kd1**7 + 7663.*kd1**6*kd2 + 15670.*kd1**5*kd2**2 + 19603.*kd1**4*kd2**3 + &
              12326.*kd1**3*kd2**4 + 2540.*kd1**2*kd2**5 + 44.*kd1*kd2**6 + 248.*kd2**7)*omega1**4 + kd2**2*(-538.*kd1**9 + 1231.*kd1**8*kd2 + 10269.*kd1**7*kd2**2 + &
              26478.*kd1**6*kd2**3 + 35260.*kd1**5*kd2**4 + 21344.*kd1**4*kd2**5 + 3393.*kd1**3*kd2**6 - 177.*kd1**2*kd2**7 + 535.*kd1*kd2**8 + 7.*kd2**9)*omega1**3*omega2 + &
              kd1*kd2*(kd1 + kd2)*(516.*kd1**8 - 90.*kd1**7*kd2 + 7252.*kd1**6*kd2**2 + 27643.*kd1**5*kd2**3 + 39422.*kd1**4*kd2**4 + 27643.*kd1**3*kd2**5 + 7252.*kd1**2*kd2**6 - &
              90.*kd1*kd2**7 + 516.*kd2**8)*omega1**2*omega2**2 + kd1**2*(7.*kd1**9 + 535.*kd1**8*kd2 - 177.*kd1**7*kd2**2 + 3393.*kd1**6*kd2**3 + 21344.*kd1**5*kd2**4 + &
              35260.*kd1**4*kd2**5 + 26478.*kd1**3*kd2**6 + 10269.*kd1**2*kd2**7 + 1231.*kd1*kd2**8 - 538.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(248.*kd1**7 + 44.*kd1**6*kd2 + &
              2540.*kd1**5*kd2**2 + 12326.*kd1**4*kd2**3 + 19603.*kd1**3*kd2**4 + 15670.*kd1**2*kd2**5 + 7663.*kd1*kd2**6 + 1996.*kd2**7)*omega2**4)*swd**2) + &
              64.*kd1**5*kd2**5*(kd1 + kd2)**5*swd*(8.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*((37.*kd1**2 + 72.*kd1*kd2 + 32.*kd2**2)*omega1**2 + (53.*kd1**2 + 82.*kd1*kd2 + &
              53.*kd2**2)*omega1*omega2 + (32.*kd1**2 + 72.*kd1*kd2 + 37.*kd2**2)*omega2**2) - omega1*omega2*(omega1 + omega2)**2*(12.*kd2**6*omega1*(omega1 - 4.*omega2) + &
              12.*kd1**6*omega2*(-4.*omega1 + omega2) + 4.*kd1**3*kd2**3*(omega1**2 - 172.*omega1*omega2 + omega2**2) + 2.*kd1*kd2**5*(164.*omega1**2 - 135.*omega1*omega2 + &
              3.*omega2**2) + kd1**2*kd2**4*(275.*omega1**2 - 602.*omega1*omega2 + 15.*omega2**2) + 2.*kd1**5*kd2*(3.*omega1**2 - 135.*omega1*omega2 + 164.*omega2**2) + &
              kd1**4*kd2**2*(15.*omega1**2 - 602.*omega1*omega2 + 275.*omega2**2))*swd**2) + 4096.*kd1**3*kd2**3*(kd1 + kd2)**3*swd*(8.*grav**2*kd1**2*kd2**2*(kd1 + &
              kd2)**2*(kd2**6*(301.*omega1**2 + 140.*omega1*omega2 + 88.*omega2**2) + 3.*kd1*kd2**5*(444.*omega1**2 + 213.*omega1*omega2 + 184.*omega2**2) + &
              kd1**6*(88.*omega1**2 + 140.*omega1*omega2 + 301.*omega2**2) + 8.*kd1**3*kd2**3*(305.*omega1**2 + 279.*omega1*omega2 + 305.*omega2**2) + &
              3.*kd1**5*kd2*(184.*omega1**2 + 213.*omega1*omega2 + 444.*omega2**2) + kd1**2*kd2**4*(2427.*omega1**2 + 1633.*omega1*omega2 + 1504.*omega2**2) + &
              kd1**4*kd2**2*(1504.*omega1**2 + 1633.*omega1*omega2 + 2427.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(-5.*kd2**10*omega1*(omega1 - 11.*omega2) + &
              5.*kd1**10*(11.*omega1 - omega2)*omega2 + kd1**8*kd2**2*(27.*omega1**2 + 3041.*omega1*omega2 - 298.*omega2**2) + 3.*kd1**7*kd2**3*(8.*omega1**2 + &
              4117.*omega1*omega2 - 170.*omega2**2) + kd1**9*kd2*(6.*omega1**2 + 516.*omega1*omega2 - 31.*omega2**2) + kd1*kd2**9*(-31.*omega1**2 + 516.*omega1*omega2 + &
              6.*omega2**2) + 3.*kd1**3*kd2**7*(-170.*omega1**2 + 4117.*omega1*omega2 + 8.*omega2**2) + kd1**2*kd2**8*(-298.*omega1**2 + 3041.*omega1*omega2 + 27.*omega2**2) + &
              6.*kd1**4*kd2**6*(266.*omega1**2 + 4712.*omega1*omega2 + 135.*omega2**2) + 6.*kd1**6*kd2**4*(135.*omega1**2 + 4712.*omega1*omega2 + 266.*omega2**2) + &
              4.*kd1**5*kd2**5*(641.*omega1**2 + 9212.*omega1*omega2 + 641.*omega2**2))*swd**2) + 262144.*kd1*kd2*(kd1 + kd2)*swd*(8.*grav**2*kd1**2*kd2**2*(kd1 + &
              kd2)**2*(kd2**10*(170.*omega1**2 + 47.*omega1*omega2 + 35.*omega2**2) + 4.*kd1*kd2**9*(347.*omega1**2 + 95.*omega1*omega2 + 67.*omega2**2) + kd1**10*(35.*omega1**2 + &
              47.*omega1*omega2 + 170.*omega2**2) + 4.*kd1**9*kd2*(67.*omega1**2 + 95.*omega1*omega2 + 347.*omega2**2) + kd1**2*kd2**8*(4955.*omega1**2 + 1723.*omega1*omega2 + &
              1580.*omega2**2) + kd1**8*kd2**2*(1580.*omega1**2 + 1723.*omega1*omega2 + 4955.*omega2**2) + kd1**3*kd2**7*(12524.*omega1**2 + 6921.*omega1*omega2 + &
              7660.*omega2**2) + kd1**7*kd2**3*(7660.*omega1**2 + 6921.*omega1*omega2 + 12524.*omega2**2) + 2.*kd1**5*kd2**5*(13873.*omega1**2 + 11288.*omega1*omega2 + &
              13873.*omega2**2) + kd1**4*kd2**6*(23139.*omega1**2 + 16821.*omega1*omega2 + 19706.*omega2**2) + kd1**6*kd2**4*(19706.*omega1**2 + 16821.*omega1*omega2 + &
              23139.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(4.*kd1**14*omega1*omega2 + 4.*kd2**14*omega1*omega2 + kd1**5*kd2**9*(1395.*omega1**2 + 125589.*omega1*omega2 - &
              4016.*omega2**2) + 2.*kd1**6*kd2**8*(3247.*omega1**2 + 110385.*omega1*omega2 - 2421.*omega2**2) + kd1**10*kd2**4*(-775.*omega1**2 + 46930.*omega1*omega2 - &
              2029.*omega2**2) + kd1**11*kd2**3*(-12.*omega1**2 + 10424.*omega1*omega2 - 1369.*omega2**2) + kd1**4*kd2**10*(-2029.*omega1**2 + 46930.*omega1*omega2 - &
              775.*omega2**2) + kd1**12*kd2**2*(-13.*omega1**2 + 1074.*omega1*omega2 - 447.*omega2**2) + kd1**2*kd2**12*(-447.*omega1**2 + 1074.*omega1*omega2 - 13.*omega2**2) + &
              kd1**3*kd2**11*(-1369.*omega1**2 + 10424.*omega1*omega2 - 12.*omega2**2) - 2.*kd1*kd2**13*(35.*omega1**2 - 27.*omega1*omega2 + omega2**2) - 2.*kd1**13*kd2*(omega1**2 - &
              27.*omega1*omega2 + 35.*omega2**2) + 2.*kd1**7*kd2**7*(1227.*omega1**2 + 132703.*omega1*omega2 + 1227.*omega2**2) + kd1**9*kd2**5*(-4016.*omega1**2 + &
              125589.*omega1*omega2 + 1395.*omega2**2) + 2.*kd1**8*kd2**6*(-2421.*omega1**2 + 110385.*omega1*omega2 + 3247.*omega2**2))*swd**2) + &
              1073741824.*swd*(8.*grav**2*kd1*kd2*(kd1 + kd2)*(24.*kd2**12*omega1**2 + 24.*kd1**12*omega2**2 + kd1*kd2**11*omega1*(378.*omega1 + 19.*omega2) + &
              kd1**11*kd2*omega2*(19.*omega1 + 378.*omega2) + kd1**2*kd2**10*(739.*omega1**2 + 66.*omega1*omega2 - 219.*omega2**2) + kd1**10*kd2**2*(-219.*omega1**2 + &
              66.*omega1*omega2 + 739.*omega2**2) + kd1**3*kd2**9*(1946.*omega1**2 + 2303.*omega1*omega2 + 1656.*omega2**2) + kd1**9*kd2**3*(1656.*omega1**2 + 2303.*omega1*omega2 + &
              1946.*omega2**2) + 2.*kd1**8*kd2**4*(3534.*omega1**2 + 2520.*omega1*omega2 + 2557.*omega2**2) + 2.*kd1**4*kd2**8*(2557.*omega1**2 + 2520.*omega1*omega2 + &
              3534.*omega2**2) + kd1**7*kd2**5*(11014.*omega1**2 + 3005.*omega1*omega2 + 8314.*omega2**2) + kd1**6*kd2**6*(10601.*omega1**2 + 480.*omega1*omega2 + &
              10601.*omega2**2) + kd1**5*kd2**7*(8314.*omega1**2 + 3005.*omega1*omega2 + 11014.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(kd2**13*omega1**2 + &
              kd1**13*omega2**2 + kd1**12*kd2*omega2*(160.*omega1 + 11.*omega2) + kd1*kd2**12*omega1*(11.*omega1 + 160.*omega2) + kd1**7*kd2**6*(-12092.*omega1**2 + &
              54009.*omega1*omega2 - 14617.*omega2**2) + kd1**6*kd2**7*(-14617.*omega1**2 + 54009.*omega1*omega2 - 12092.*omega2**2) + kd1**8*kd2**5*(-433.*omega1**2 + &
              26790.*omega1*omega2 - 7195.*omega2**2) + kd1**9*kd2**4*(4041.*omega1**2 + 6478.*omega1*omega2 - 1139.*omega2**2) + kd1**5*kd2**8*(-7195.*omega1**2 + &
              26790.*omega1*omega2 - 433.*omega2**2) + kd1**2*kd2**11*(283.*omega1**2 + 2351.*omega1*omega2 + 46.*omega2**2) + kd1**11*kd2**2*(46.*omega1**2 + 2351.*omega1*omega2 + &
              283.*omega2**2) + kd1**10*kd2**3*(1365.*omega1**2 + 3679.*omega1*omega2 + 594.*omega2**2) + kd1**3*kd2**10*(594.*omega1**2 + 3679.*omega1*omega2 + &
              1365.*omega2**2) + kd1**4*kd2**9*(-1139.*omega1**2 + 6478.*omega1*omega2 + 4041.*omega2**2))*swd**2) + 16777216.*swd*(8.*grav**2*kd1*kd2*(kd1 + &
              kd2)**2*(9.*kd2**13*omega1**2 + kd1*kd2**12*omega1*(71.*omega1 - 9.*omega2) + 9.*kd1**13*omega2**2 + kd1**12*kd2*omega2*(-9.*omega1 + 71.*omega2) + &
              14.*kd1**2*kd2**11*(21.*omega1**2 + 5.*omega1*omega2 + 9.*omega2**2) + 14.*kd1**11*kd2**2*(9.*omega1**2 + 5.*omega1*omega2 + 21.*omega2**2) + &
              2.*kd1**3*kd2**10*(1862.*omega1**2 + 715.*omega1*omega2 + 871.*omega2**2) + 2.*kd1**10*kd2**3*(871.*omega1**2 + 715.*omega1*omega2 + 1862.*omega2**2) + &
              kd1**4*kd2**9*(14625.*omega1**2 + 4217.*omega1*omega2 + 6065.*omega2**2) + 2.*kd1**6*kd2**7*(18425.*omega1**2 + 10516.*omega1*omega2 + 14474.*omega2**2) + &
              kd1**9*kd2**4*(6065.*omega1**2 + 4217.*omega1*omega2 + 14625.*omega2**2) + kd1**5*kd2**8*(29631.*omega1**2 + 10913.*omega1*omega2 + 15047.*omega2**2) + &
              2.*kd1**7*kd2**6*(14474.*omega1**2 + 10516.*omega1*omega2 + 18425.*omega2**2) + kd1**8*kd2**5*(15047.*omega1**2 + 10913.*omega1*omega2 + 29631.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(kd2**15*omega1**2 + kd1*kd2**14*omega1*(5.*omega1 - 68.*omega2) + kd1**15*omega2**2 + kd1**14*kd2*omega2*(-68.*omega1 + &
              5.*omega2) + kd1**9*kd2**6*(8791.*omega1**2 - 230693.*omega1*omega2 - 13087.*omega2**2) - 3.*kd1**7*kd2**8*(428.*omega1**2 + 107400.*omega1*omega2 - &
              4111.*omega2**2) + 3.*kd1**8*kd2**7*(4111.*omega1**2 - 107400.*omega1*omega2 - 428.*omega2**2) + kd1**13*kd2**2*(8.*omega1**2 - 757.*omega1*omega2 - 99.*omega2**2) + &
              kd1**2*kd2**13*(-99.*omega1**2 - 757.*omega1*omega2 + 8.*omega2**2) - 2.*kd1**3*kd2**12*(485.*omega1**2 + 1998.*omega1*omega2 + 116.*omega2**2) - &
              kd1**5*kd2**10*(10605.*omega1**2 + 106859.*omega1*omega2 + 197.*omega2**2) - 2.*kd1**12*kd2**3*(116.*omega1**2 + 1998.*omega1*omega2 + 485.*omega2**2) - &
              kd1**4*kd2**11*(4229.*omega1**2 + 27574.*omega1*omega2 + 1470.*omega2**2) - kd1**11*kd2**4*(1470.*omega1**2 + 27574.*omega1*omega2 + 4229.*omega2**2) + &
              kd1**6*kd2**9*(-13087.*omega1**2 - 230693.*omega1*omega2 + 8791.*omega2**2) - kd1**10*kd2**5*(197.*omega1**2 + 106859.*omega1*omega2 + 10605.*omega2**2))*swd**2) - &
              512.*grav*kd1**5*kd2**5*(kd1 + kd2)**5*(640.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4 + (2.*kd2**6*(252.*omega1**4 + 409.*omega1**3*omega2 + 273.*omega1**2*omega2**2 + &
              59.*omega1*omega2**3 + 33.*omega2**4) + 2.*kd1**6*(33.*omega1**4 + 59.*omega1**3*omega2 + 273.*omega1**2*omega2**2 + 409.*omega1*omega2**3 + 252.*omega2**4) + &
              kd1*kd2**5*(2182.*omega1**4 + 3134.*omega1**3*omega2 + 2772.*omega1**2*omega2**2 + 615.*omega1*omega2**3 + 432.*omega2**4) + kd1**2*kd2**4*(3730.*omega1**4 + &
              6785.*omega1**3*omega2 + 7304.*omega1**2*omega2**2 + 2917.*omega1*omega2**3 + 1472.*omega2**4) + 2.*kd1**3*kd2**3*(1595.*omega1**4 + 3291.*omega1**3*omega2 + &
              4739.*omega1**2*omega2**2 + 3291.*omega1*omega2**3 + 1595.*omega2**4) + kd1**5*kd2*(432.*omega1**4 + 615.*omega1**3*omega2 + 2772.*omega1**2*omega2**2 + &
              3134.*omega1*omega2**3 + 2182.*omega2**4) + kd1**4*kd2**2*(1472.*omega1**4 + 2917.*omega1**3*omega2 + 7304.*omega1**2*omega2**2 + 6785.*omega1*omega2**3 + &
              3730.*omega2**4))*swd**2) + 8589934592.*grav*(128.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(19.*kd1**6 - 5.*kd1**5*kd2 + 13.*kd1**4*kd2**2 + 43.*kd1**3*kd2**3 + &
              13.*kd1**2*kd2**4 - 5.*kd1*kd2**5 + 19.*kd2**6) - (11.*kd2**13*omega1**3*omega2 + 11.*kd1**13*omega1*omega2**3 + kd1**12*kd2*omega2**2*(96.*omega1**2 + &
              185.*omega1*omega2 + 64.*omega2**2) + kd1*kd2**12*omega1**2*(64.*omega1**2 + 185.*omega1*omega2 + 96.*omega2**2) + kd1**2*kd2**11*omega1*(2420.*omega1**3 + &
              4821.*omega1**2*omega2 + 5936.*omega1*omega2**2 + 609.*omega2**3) + kd1**11*kd2**2*omega2*(609.*omega1**3 + 5936.*omega1**2*omega2 + 4821.*omega1*omega2**2 + &
              2420.*omega2**3) + kd1**8*kd2**5*(4390.*omega1**4 + 19091.*omega1**3*omega2 - 74332.*omega1**2*omega2**2 - 61553.*omega1*omega2**3 - 25998.*omega2**4) + &
              kd1**3*kd2**10*(1844.*omega1**4 + 4111.*omega1**3*omega2 + 12316.*omega1**2*omega2**2 + 1445.*omega1*omega2**3 - 2272.*omega2**4) - kd1**4*kd2**9*(10658.*omega1**4 + &
              21885.*omega1**3*omega2 + 7342.*omega1**2*omega2**2 - 12143.*omega1*omega2**3 + 1464.*omega2**4) + kd1**10*kd2**3*(-2272.*omega1**4 + 1445.*omega1**3*omega2 + &
              12316.*omega1**2*omega2**2 + 4111.*omega1*omega2**3 + 1844.*omega2**4) + kd1**5*kd2**8*(-25998.*omega1**4 - 61553.*omega1**3*omega2 - 74332.*omega1**2*omega2**2 + &
              19091.*omega1*omega2**3 + 4390.*omega2**4) - kd1**6*kd2**7*(24008.*omega1**4 + 67537.*omega1**3*omega2 + 152462.*omega1**2*omega2**2 + 19533.*omega1*omega2**3 + &
              4566.*omega2**4) - kd1**9*kd2**4*(1464.*omega1**4 - 12143.*omega1**3*omega2 + 7342.*omega1**2*omega2**2 + 21885.*omega1*omega2**3 + 10658.*omega2**4) - &
              kd1**7*kd2**6*(4566.*omega1**4 + 19533.*omega1**3*omega2 + 152462.*omega1**2*omega2**2 + 67537.*omega1*omega2**3 + 24008.*omega2**4))*swd**2) - &
              32768.*grav*kd1**3*kd2**3*(kd1 + kd2)**3*(64.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(32.*kd1**4 + 96.*kd1**3*kd2 + 95.*kd1**2*kd2**2 + 96.*kd1*kd2**3 + 32.*kd2**4) + &
              (4.*kd2**10*(116.*omega1**4 + 232.*omega1**3*omega2 + 167.*omega1**2*omega2**2 + 12.*omega1*omega2**3 + 6.*omega2**4) + 4.*kd1**10*(6.*omega1**4 + 12.*omega1**3*omega2 + &
              167.*omega1**2*omega2**2 + 232.*omega1*omega2**3 + 116.*omega2**4) + kd1*kd2**9*(3546.*omega1**4 + 6625.*omega1**3*omega2 + 5332.*omega1**2*omega2**2 + &
              466.*omega1*omega2**3 + 246.*omega2**4) + 2.*kd1**2*kd2**8*(6685.*omega1**4 + 12259.*omega1**3*omega2 + 11794.*omega1**2*omega2**2 + 1948.*omega1*omega2**3 + &
              1121.*omega2**4) + kd1**9*kd2*(246.*omega1**4 + 466.*omega1**3*omega2 + 5332.*omega1**2*omega2**2 + 6625.*omega1*omega2**3 + 3546.*omega2**4) + &
              2.*kd1**3*kd2**7*(16648.*omega1**4 + 32338.*omega1**3*omega2 + 39848.*omega1**2*omega2**2 + 12522.*omega1*omega2**3 + 6681.*omega2**4) + &
              2.*kd1**8*kd2**2*(1121.*omega1**4 + 1948.*omega1**3*omega2 + 11794.*omega1**2*omega2**2 + 12259.*omega1*omega2**3 + 6685.*omega2**4) + 2.*kd1**7*kd2**3*(6681.*omega1**4 + &
              12522.*omega1**3*omega2 + 39848.*omega1**2*omega2**2 + 32338.*omega1*omega2**3 + 16648.*omega2**4) + kd1**4*kd2**6*(55726.*omega1**4 + 118435.*omega1**3*omega2 + &
              179446.*omega1**2*omega2**2 + 79771.*omega1*omega2**3 + 38142.*omega2**4) + kd1**6*kd2**4*(38142.*omega1**4 + 79771.*omega1**3*omega2 + 179446.*omega1**2*omega2**2 + &
              118435.*omega1*omega2**3 + 55726.*omega2**4) + kd1**5*kd2**5*(59542.*omega1**4 + 130733.*omega1**3*omega2 + 237196.*omega1**2*omega2**2 + 130733.*omega1*omega2**3 + &
              59542.*omega2**4))*swd**2) - 134217728.*grav*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(5.*kd1**8 + 170.*kd1**7*kd2 + 456.*kd1**6*kd2**2 + 662.*kd1**5*kd2**3 + &
              739.*kd1**4*kd2**4 + 662.*kd1**3*kd2**5 + 456.*kd1**2*kd2**6 + 170.*kd1*kd2**7 + 5.*kd2**8) + (-3.*kd2**15*omega1**3*omega2 - 3.*kd1**15*omega1*omega2**3 + &
              kd1**14*kd2*omega2**2*(300.*omega1**2 + 241.*omega1*omega2 + 128.*omega2**2) + kd1*kd2**14*omega1**2*(128.*omega1**2 + 241.*omega1*omega2 + 300.*omega2**2) + &
              kd1**2*kd2**13*omega1*(1508.*omega1**3 + 3303.*omega1**2*omega2 + 3614.*omega1*omega2**2 - 58.*omega2**3) + kd1**13*kd2**2*omega2*(-58.*omega1**3 + 3614.*omega1**2*omega2 + &
              3303.*omega1*omega2**2 + 1508.*omega2**3) + kd1**3*kd2**12*(3786.*omega1**4 + 10136.*omega1**3*omega2 + 12366.*omega1**2*omega2**2 + 1621.*omega1*omega2**3 + &
              296.*omega2**4) + kd1**12*kd2**3*(296.*omega1**4 + 1621.*omega1**3*omega2 + 12366.*omega1**2*omega2**2 + 10136.*omega1*omega2**3 + 3786.*omega2**4) + &
              2.*kd1**4*kd2**11*(10305.*omega1**4 + 24161.*omega1**3*omega2 + 31434.*omega1**2*omega2**2 + 11707.*omega1*omega2**3 + 4465.*omega2**4) + &
              2.*kd1**11*kd2**4*(4465.*omega1**4 + 11707.*omega1**3*omega2 + 31434.*omega1**2*omega2**2 + 24161.*omega1*omega2**3 + 10305.*omega2**4) + &
              kd1**5*kd2**10*(81942.*omega1**4 + 166165.*omega1**3*omega2 + 234512.*omega1**2*omega2**2 + 75199.*omega1*omega2**3 + 36956.*omega2**4) + &
              2.*kd1**7*kd2**8*(88804.*omega1**4 + 144871.*omega1**3*omega2 + 355553.*omega1**2*omega2**2 + 99174.*omega1*omega2**3 + 66819.*omega2**4) + &
              kd1**6*kd2**9*(159804.*omega1**4 + 292383.*omega1**3*omega2 + 506278.*omega1**2*omega2**2 + 129537.*omega1*omega2**3 + 80420.*omega2**4) + &
              kd1**10*kd2**5*(36956.*omega1**4 + 75199.*omega1**3*omega2 + 234512.*omega1**2*omega2**2 + 166165.*omega1*omega2**3 + 81942.*omega2**4) + &
              2.*kd1**8*kd2**7*(66819.*omega1**4 + 99174.*omega1**3*omega2 + 355553.*omega1**2*omega2**2 + 144871.*omega1*omega2**3 + 88804.*omega2**4) + &
              kd1**9*kd2**6*(80420.*omega1**4 + 129537.*omega1**3*omega2 + 506278.*omega1**2*omega2**2 + 292383.*omega1*omega2**3 + 159804.*omega2**4))*swd**2) - &
              2097152.*grav*kd1*kd2*(kd1 + kd2)*(128.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(3.*kd1**8 + 20.*kd1**7*kd2 + 48.*kd1**6*kd2**2 + 229.*kd1**5*kd2**3 + &
              397.*kd1**4*kd2**4 + 229.*kd1**3*kd2**5 + 48.*kd1**2*kd2**6 + 20.*kd1*kd2**7 + 3.*kd2**8) + (12.*kd1**14*omega2**2*(3.*omega1**2 + 4.*omega1*omega2 + &
              2.*omega2**2) + 12.*kd2**14*omega1**2*(2.*omega1**2 + 4.*omega1*omega2 + 3.*omega2**2) + kd1*kd2**13*omega1*(290.*omega1**3 + 476.*omega1**2*omega2 + &
              304.*omega1*omega2**2 - 25.*omega2**3) + kd1**13*kd2*omega2*(-25.*omega1**3 + 304.*omega1**2*omega2 + 476.*omega1*omega2**2 + 290.*omega2**3) + &
              kd1**2*kd2**12*(2398.*omega1**4 + 4299.*omega1**3*omega2 + 4592.*omega1**2*omega2**2 + 563.*omega1*omega2**3 + 304.*omega2**4) + kd1**12*kd2**2*(304.*omega1**4 + &
              563.*omega1**3*omega2 + 4592.*omega1**2*omega2**2 + 4299.*omega1*omega2**3 + 2398.*omega2**4) + kd1**3*kd2**11*(18262.*omega1**4 + 36947.*omega1**3*omega2 + &
              46292.*omega1**2*omega2**2 + 7752.*omega1*omega2**3 + 3692.*omega2**4) + kd1**11*kd2**3*(3692.*omega1**4 + 7752.*omega1**3*omega2 + 46292.*omega1**2*omega2**2 + &
              36947.*omega1*omega2**3 + 18262.*omega2**4) + kd1**4*kd2**10*(70076.*omega1**4 + 145863.*omega1**3*omega2 + 194628.*omega1**2*omega2**2 + 34821.*omega1*omega2**3 + &
              18326.*omega2**4) + 2.*kd1**5*kd2**9*(79465.*omega1**4 + 167361.*omega1**3*omega2 + 240694.*omega1**2*omega2**2 + 60812.*omega1*omega2**3 + 32749.*omega2**4) + &
              kd1**10*kd2**4*(18326.*omega1**4 + 34821.*omega1**3*omega2 + 194628.*omega1**2*omega2**2 + 145863.*omega1*omega2**3 + 70076.*omega2**4) + &
              2.*kd1**9*kd2**5*(32749.*omega1**4 + 60812.*omega1**3*omega2 + 240694.*omega1**2*omega2**2 + 167361.*omega1*omega2**3 + 79465.*omega2**4) + &
              kd1**6*kd2**8*(238420.*omega1**4 + 502403.*omega1**3*omega2 + 811942.*omega1**2*omega2**2 + 308427.*omega1*omega2**3 + 157416.*omega2**4) + &
              kd1**8*kd2**6*(157416.*omega1**4 + 308427.*omega1**3*omega2 + 811942.*omega1**2*omega2**2 + 502403.*omega1*omega2**3 + 238420.*omega2**4) + &
              kd1**7*kd2**7*(240292.*omega1**4 + 494923.*omega1**3*omega2 + 966418.*omega1**2*omega2**2 + 494923.*omega1*omega2**3 + 240292.*omega2**4))*swd**2)))/ &
              ((-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*omega2*(omega1 + omega2)*swd**3*(-64.*grav*(kd1 + kd2)**2*(64. + &
              (kd1 + kd2)**2)*(4096. + (kd1 + kd2)**2*(384. + (kd1 + kd2)**2)) + (16777216. + (kd1 + kd2)**2*(7340032. + (kd1 + kd2)**2*(286720. + (kd1 + kd2)**2*(1792. + &
              (kd1 + kd2)**2))))*(omega1 + omega2)**2*swd))
    !
end function velsp24
!
real function velsb34()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic velocity of 3rd layer for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb34')
    !
    velsb34 = (576460752303423488.*grav*(kd1 - kd2)*(-6.*kd1*kd2**3*omega1**4 + kd2**2*(kd1**2 + 10.*kd1*kd2 + kd2**2)*omega1**3*omega2 - 4.*kd1*kd2*(kd1**2 + &
              kd1*kd2 + kd2**2)*omega1**2*omega2**2 + kd1**2*(kd1**2 + 10.*kd1*kd2 + kd2**2)*omega1*omega2**3 - 6.*kd1**3*kd2*omega2**4)*swd**2 + 64.*grav*kd1**7*(kd1 - &
              kd2)**7*kd2**7*((14.*kd1**2 - 33.*kd1*kd2 + 16.*kd2**2)*omega1**4 + (-30.*kd1**2 + 38.*kd1*kd2 - 51.*kd2**2)*omega1**3*omega2 + (61.*kd1**2 - 30.*kd1*kd2 + &
              61.*kd2**2)*omega1**2*omega2**2 + (-51.*kd1**2 + 38.*kd1*kd2 - 30.*kd2**2)*omega1*omega2**3 + (16.*kd1**2 - 33.*kd1*kd2 + 14.*kd2**2)*omega2**4)*swd**2 + &
              3458764513820540928.*kd1*(kd1 - kd2)*kd2*omega1**2*(omega1 - omega2)**2*omega2**2*swd**3 + kd1**7*(kd1 - kd2)**7*kd2**7*omega1*(omega1 - &
              omega2)**2*omega2*(kd2**2*omega1*(5.*omega1 - 21.*omega2) + kd1**2*omega2*(-21.*omega1 + 5.*omega2) - 8.*kd1*kd2*(omega1**2 - 6.*omega1*omega2 + &
              omega2**2))*swd**3 - 18014398509481984.*swd*(32.*grav**2*kd1**2*(kd1 - kd2)*kd2**2*(-7.*kd1*kd2*omega1**2 + (kd1**2 + 12.*kd1*kd2 + kd2**2)*omega1*omega2 - &
              7.*kd1*kd2*omega2**2) + omega1*(omega1 - omega2)**2*omega2*(kd2**2*(40.*kd1**3 - 91.*kd1**2*kd2 + 59.*kd1*kd2**2 - 20.*kd2**3)*omega1**2 + kd1*(kd1 - &
              kd2)*kd2*(9.*kd1**2 - 182.*kd1*kd2 + 9.*kd2**2)*omega1*omega2 + kd1**2*(20.*kd1**3 - 59.*kd1**2*kd2 + 91.*kd1*kd2**2 - 40.*kd2**3)*omega2**2)*swd**2) - &
              281474976710656.*swd*(32.*grav**2*kd1*(kd1 - kd2)*kd2*(21.*kd2**6*omega1**2 + kd1**5*kd2*(13.*omega1 - 98.*omega2)*omega2 + 21.*kd1**6*omega2**2 + &
              kd1*kd2**5*omega1*(-98.*omega1 + 13.*omega2) + kd1**2*kd2**4*(135.*omega1**2 - 95.*omega1*omega2 + 8.*omega2**2) - 4.*kd1**3*kd2**3*(31.*omega1**2 - &
              70.*omega1*omega2 + 31.*omega2**2) + kd1**4*kd2**2*(8.*omega1**2 - 95.*omega1*omega2 + 135.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(2.*kd2**2*(80.*kd1**5 - &
              691.*kd1**4*kd2 + 938.*kd1**3*kd2**2 - 396.*kd1**2*kd2**3 + 17.*kd1*kd2**4 + 10.*kd2**5)*omega1**2 + kd1*(kd1 - kd2)*kd2*(753.*kd1**4 - 2600.*kd1**3*kd2 + &
              1793.*kd1**2*kd2**2 - 2600.*kd1*kd2**3 + 753.*kd2**4)*omega1*omega2 - 2.*kd1**2*(10.*kd1**5 + 17.*kd1**4*kd2 - 396.*kd1**3*kd2**2 + 938.*kd1**2*kd2**3 - &
              691.*kd1*kd2**4 + 80.*kd2**5)*omega2**2)*swd**2) + 140737488355328.*grav*(128.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(2.*kd1**2 - 13.*kd1*kd2 + 2.*kd2**2) + &
              (2.*kd1*kd2**3*(153.*kd1**5 - 1143.*kd1**4*kd2 + 1843.*kd1**3*kd2**2 - 1684.*kd1**2*kd2**3 + 1049.*kd1*kd2**4 - 230.*kd2**5)*omega1**4 + kd2**2*(15.*kd1**7 - &
              2208.*kd1**6*kd2 + 9199.*kd1**5*kd2**2 - 12159.*kd1**4*kd2**3 + 8583.*kd1**3*kd2**4 - 4224.*kd1**2*kd2**5 + 820.*kd1*kd2**6 + 22.*kd2**7)*omega1**3*omega2 + &
              2.*kd1*(kd1 - kd2)*kd2*(453.*kd1**6 - 1870.*kd1**5*kd2 + 3441.*kd1**4*kd2**2 - 5059.*kd1**3*kd2**3 + 3441.*kd1**2*kd2**4 - 1870.*kd1*kd2**5 + &
              453.*kd2**6)*omega1**2*omega2**2 - kd1**2*(22.*kd1**7 + 820.*kd1**6*kd2 - 4224.*kd1**5*kd2**2 + 8583.*kd1**4*kd2**3 - 12159.*kd1**3*kd2**4 + 9199.*kd1**2*kd2**5 - &
              2208.*kd1*kd2**6 + 15.*kd2**7)*omega1*omega2**3 + 2.*kd1**3*kd2*(230.*kd1**5 - 1049.*kd1**4*kd2 + 1684.*kd1**3*kd2**2 - 1843.*kd1**2*kd2**3 + 1143.*kd1*kd2**4 - &
              153.*kd2**5)*omega2**4)*swd**2) + 4398046511104.*grav*(64.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(12.*kd1**4 - 67.*kd1**3*kd2 + 54.*kd1**2*kd2**2 - 67.*kd1*kd2**3 + &
              12.*kd2**4) + (kd1*kd2**3*(657.*kd1**7 - 5042.*kd1**6*kd2 + 10123.*kd1**5*kd2**2 - 9344.*kd1**4*kd2**3 + 6066.*kd1**3*kd2**4 - 4790.*kd1**2*kd2**5 + 2914.*kd1*kd2**6 - &
              566.*kd2**7)*omega1**4 - kd2**2*(386.*kd1**9 + 1382.*kd1**8*kd2 - 13262.*kd1**7*kd2**2 + 25243.*kd1**6*kd2**3 - 20999.*kd1**5*kd2**4 + 11102.*kd1**4*kd2**5 - &
              7915.*kd1**3*kd2**6 + 5237.*kd1**2*kd2**7 - 1145.*kd1*kd2**8 + 7.*kd2**9)*omega1**3*omega2 + kd1*(kd1 - kd2)*kd2*(1146.*kd1**8 - 5022.*kd1**7*kd2 + 8808.*kd1**6*kd2**2 - &
              15973.*kd1**5*kd2**3 + 19914.*kd1**4*kd2**4 - 15973.*kd1**3*kd2**5 + 8808.*kd1**2*kd2**6 - 5022.*kd1*kd2**7 + 1146.*kd2**8)*omega1**2*omega2**2 + kd1**2*(7.*kd1**9 - &
              1145.*kd1**8*kd2 + 5237.*kd1**7*kd2**2 - 7915.*kd1**6*kd2**3 + 11102.*kd1**5*kd2**4 - 20999.*kd1**4*kd2**5 + 25243.*kd1**3*kd2**6 - 13262.*kd1**2*kd2**7 + &
              1382.*kd1*kd2**8 + 386.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(566.*kd1**7 - 2914.*kd1**6*kd2 + 4790.*kd1**5*kd2**2 - 6066.*kd1**4*kd2**3 + 9344.*kd1**3*kd2**4 - &
              10123.*kd1**2*kd2**5 + 5042.*kd1*kd2**6 - 657.*kd2**7)*omega2**4)*swd**2) + 34359738368.*grav*(128.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(25.*kd1**6 - 129.*kd1**5*kd2 + &
              59.*kd1**4*kd2**2 + 39.*kd1**3*kd2**3 + 59.*kd1**2*kd2**4 - 129.*kd1*kd2**5 + 25.*kd2**6) + (2.*kd1*kd2**3*(1607.*kd1**9 - 11915.*kd1**8*kd2 + 27193.*kd1**7*kd2**2 - &
              30775.*kd1**6*kd2**3 + 23145.*kd1**5*kd2**4 - 14306.*kd1**4*kd2**5 + 8375.*kd1**3*kd2**6 - 6142.*kd1**2*kd2**7 + 3351.*kd1*kd2**8 - 545.*kd2**9)*omega1**4 - &
              kd2**2*(913.*kd1**11 + 4951.*kd1**10*kd2 - 52389.*kd1**9*kd2**2 + 112641.*kd1**8*kd2**3 - 101465.*kd1**7*kd2**4 + 55949.*kd1**6*kd2**5 - 39157.*kd1**5*kd2**6 + &
              33005.*kd1**4*kd2**7 - 25969.*kd1**3*kd2**8 + 13715.*kd1**2*kd2**9 - 2247.*kd1*kd2**10 + 5.*kd2**11)*omega1**3*omega2 + 2.*kd1*(kd1 - kd2)*kd2*(1080.*kd1**10 - &
              5358.*kd1**9*kd2 + 9574.*kd1**8*kd2**2 - 25595.*kd1**7*kd2**3 + 42837.*kd1**6*kd2**4 - 45840.*kd1**5*kd2**5 + 42837.*kd1**4*kd2**6 - 25595.*kd1**3*kd2**7 + &
              9574.*kd1**2*kd2**8 - 5358.*kd1*kd2**9 + 1080.*kd2**10)*omega1**2*omega2**2 + kd1**2*(5.*kd1**11 - 2247.*kd1**10*kd2 + 13715.*kd1**9*kd2**2 - 25969.*kd1**8*kd2**3 + &
              33005.*kd1**7*kd2**4 - 39157.*kd1**6*kd2**5 + 55949.*kd1**5*kd2**6 - 101465.*kd1**4*kd2**7 + 112641.*kd1**3*kd2**8 - 52389.*kd1**2*kd2**9 + 4951.*kd1*kd2**10 + &
              913.*kd2**11)*omega1*omega2**3 + 2.*kd1**3*kd2*(545.*kd1**9 - 3351.*kd1**8*kd2 + 6142.*kd1**7*kd2**2 - 8375.*kd1**6*kd2**3 + 14306.*kd1**5*kd2**4 - 23145.*kd1**4*kd2**5 + &
              30775.*kd1**3*kd2**6 - 27193.*kd1**2*kd2**7 + 11915.*kd1*kd2**8 - 1607.*kd2**9)*omega2**4)*swd**2) - 64.*kd1**5*(kd1 - kd2)**5*kd2**5*swd*(32.*grav**2*kd1**2*(kd1 - &
              kd2)**2*kd2**2*((45.*kd1**2 - 88.*kd1*kd2 + 40.*kd2**2)*omega1**2 + (-53.*kd1**2 + 86.*kd1*kd2 - 53.*kd2**2)*omega1*omega2 + (40.*kd1**2 - 88.*kd1*kd2 + &
              45.*kd2**2)*omega2**2) - omega1*(omega1 - omega2)**2*omega2*(kd2**6*omega1*(13.*omega1 - 253.*omega2) + kd1**6*omega2*(-253.*omega1 + 13.*omega2) + &
              2.*kd1*kd2**5*(387.*omega1**2 + 601.*omega1*omega2 + 12.*omega2**2) - 2.*kd1**3*kd2**3*(53.*omega1**2 - 1525.*omega1*omega2 + 53.*omega2**2) - &
              kd1**2*kd2**4*(565.*omega1**2 + 2551.*omega1*omega2 + 63.*omega2**2) + 2.*kd1**5*kd2*(12.*omega1**2 + 601.*omega1*omega2 + 387.*omega2**2) - kd1**4*kd2**2*(63.*omega1**2 + &
              2551.*omega1*omega2 + 565.*omega2**2))*swd**2) + 4096.*kd1**3*(kd1 - kd2)**3*kd2**3*swd*(32.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd1**4*kd2**2*(-1538.*omega1**2 + &
              1479.*omega1*omega2 - 2555.*omega2**2) + kd1**2*kd2**4*(-2555.*omega1**2 + 1479.*omega1*omega2 - 1538.*omega2**2) + kd1**6*(-112.*omega1**2 + 132.*omega1*omega2 - &
              339.*omega2**2) + kd2**6*(-339.*omega1**2 + 132.*omega1*omega2 - 112.*omega2**2) + 71.*kd1*kd2**5*(20.*omega1**2 - 7.*omega1*omega2 + 8.*omega2**2) + &
              71.*kd1**5*kd2*(8.*omega1**2 - 7.*omega1*omega2 + 20.*omega2**2) + 8.*kd1**3*kd2**3*(320.*omega1**2 - 277.*omega1*omega2 + 320.*omega2**2)) - &
              omega1*(omega1 - omega2)**2*omega2*(kd1**10*omega2*(223.*omega1 + 17.*omega2) + kd2**10*omega1*(17.*omega1 + 223.*omega2) + kd1**6*kd2**4*(-876.*omega1**2 + &
              101101.*omega1*omega2 - 6145.*omega2**2) + kd1**4*kd2**6*(-6145.*omega1**2 + 101101.*omega1*omega2 - 876.*omega2**2) - 2.*kd1**7*kd2**3*(138.*omega1**2 + &
              24111.*omega1*omega2 - 245.*omega2**2) + 2.*kd1**3*kd2**7*(245.*omega1**2 - 24111.*omega1*omega2 - 138.*omega2**2) + kd1**2*kd2**8*(196.*omega1**2 + &
              12378.*omega1*omega2 - 111.*omega2**2) + 4.*kd1**9*kd2*(6.*omega1**2 - 439.*omega1*omega2 - 68.*omega2**2) - 4.*kd1*kd2**9*(68.*omega1**2 + 439.*omega1*omega2 - &
              6.*omega2**2) + kd1**8*kd2**2*(-111.*omega1**2 + 12378.*omega1*omega2 + 196.*omega2**2) + 2.*kd1**5*kd2**5*(3513.*omega1**2 - 63797.*omega1*omega2 + &
              3513.*omega2**2))*swd**2) - 4398046511104.*swd*(32.*grav**2*kd1*(kd1 - kd2)*kd2*(106.*kd2**8*omega1**2 + 106.*kd1**8*omega2**2 - 2.*kd1*kd2**7*omega1*(269.*omega1 + &
              15.*omega2) - 2.*kd1**7*kd2*omega2*(15.*omega1 + 269.*omega2) + kd1**2*kd2**6*(1060.*omega1**2 - 279.*omega1*omega2 + 201.*omega2**2) - 2.*kd1**3*kd2**5*(937.*omega1**2 - &
              850.*omega1*omega2 + 685.*omega2**2) - 2.*kd1**5*kd2**3*(685.*omega1**2 - 850.*omega1*omega2 + 937.*omega2**2) + kd1**6*kd2**2*(201.*omega1**2 - 279.*omega1*omega2 + &
              1060.*omega2**2) + kd1**4*kd2**4*(2187.*omega1**2 - 2326.*omega1*omega2 + 2187.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(56.*kd2**9*omega1**2 - &
              56.*kd1**9*omega2**2 + kd1**8*kd2*omega2*(3181.*omega1 + 407.*omega2) - kd1*kd2**8*omega1*(407.*omega1 + 3181.*omega2) + kd1**4*kd2**5*(-6897.*omega1**2 + &
              28903.*omega1*omega2 - 10423.*omega2**2) + kd1**2*kd2**7*(1693.*omega1**2 + 15523.*omega1*omega2 + 536.*omega2**2) + kd1**6*kd2**3*(-3293.*omega1**2 + &
              25563.*omega1*omega2 + 775.*omega2**2) - kd1**7*kd2**2*(536.*omega1**2 + 15523.*omega1*omega2 + 1693.*omega2**2) + kd1**3*kd2**6*(-775.*omega1**2 - &
              25563.*omega1*omega2 + 3293.*omega2**2) + kd1**5*kd2**4*(10423.*omega1**2 - 28903.*omega1*omega2 + 6897.*omega2**2))*swd**2) - &
              68719476736.*swd*(64.*grav**2*kd1*kd2*(-kd1 + kd2)*(-78.*kd2**10*omega1**2 - 78.*kd1**10*omega2**2 + kd1*kd2**9*omega1*(419.*omega1 + 31.*omega2) + &
              kd1**9*kd2*omega2*(31.*omega1 + 419.*omega2) + kd1**4*kd2**6*(-2455.*omega1**2 + 1834.*omega1*omega2 - 3593.*omega2**2) + kd1**6*kd2**4*(-3593.*omega1**2 + &
              1834.*omega1*omega2 - 2455.*omega2**2) + kd1**8*kd2**2*(-268.*omega1**2 + 223.*omega1*omega2 - 744.*omega2**2) + kd1**2*kd2**8*(-744.*omega1**2 + 223.*omega1*omega2 - &
              268.*omega2**2) + kd1**7*kd2**3*(1940.*omega1**2 - 1638.*omega1*omega2 + 1451.*omega2**2) + kd1**3*kd2**7*(1451.*omega1**2 - 1638.*omega1*omega2 + 1940.*omega2**2) + &
              kd1**5*kd2**5*(3415.*omega1**2 - 1074.*omega1*omega2 + 3415.*omega2**2)) - omega1*(omega1 - omega2)**2*omega2*(88.*kd2**11*omega1**2 - 88.*kd1**11*omega2**2 + &
              kd1**10*kd2*omega2*(-5079.*omega1 + 428.*omega2) + kd1*kd2**10*omega1*(-428.*omega1 + 5079.*omega2) + kd1**5*kd2**6*(3871.*omega1**2 + 95625.*omega1*omega2 - &
              24995.*omega2**2) + kd1**8*kd2**3*(820.*omega1**2 - 49565.*omega1*omega2 - 7621.*omega2**2) + kd1**6*kd2**5*(24995.*omega1**2 - 95625.*omega1*omega2 - &
              3871.*omega2**2) + kd1**3*kd2**8*(7621.*omega1**2 + 49565.*omega1*omega2 - 820.*omega2**2) + kd1**9*kd2**2*(1792.*omega1**2 + 28567.*omega1*omega2 + 360.*omega2**2) - &
              kd1**2*kd2**9*(360.*omega1**2 + 28567.*omega1*omega2 + 1792.*omega2**2) + kd1**7*kd2**4*(-19896.*omega1**2 + 64727.*omega1*omega2 + 18335.*omega2**2) + &
              kd1**4*kd2**7*(-18335.*omega1**2 - 64727.*omega1*omega2 + 19896.*omega2**2))*swd**2) - 262144.*kd1*(kd1 - kd2)*kd2*swd*(32.*grav**2*kd1**2*(kd1 - &
              kd2)**2*kd2**2*(kd1**7*kd2**3*(-8164.*omega1**2 + 7215.*omega1*omega2 - 11724.*omega2**2) + kd1**3*kd2**7*(-11724.*omega1**2 + 7215.*omega1*omega2 - 8164.*omega2**2) + &
              kd2**10*(230.*omega1**2 - 39.*omega1*omega2 + 35.*omega2**2) - 4.*kd1*kd2**9*(381.*omega1**2 - 67.*omega1*omega2 + 55.*omega2**2) + kd1**10*(35.*omega1**2 - &
              39.*omega1*omega2 + 230.*omega2**2) - 4.*kd1**9*kd2*(55.*omega1**2 - 67.*omega1*omega2 + 381.*omega2**2) + kd1**2*kd2**8*(5005.*omega1**2 - 1935.*omega1*omega2 + &
              1816.*omega2**2) + kd1**8*kd2**2*(1816.*omega1**2 - 1935.*omega1*omega2 + 5005.*omega2**2) - 2.*kd1**5*kd2**5*(11433.*omega1**2 - 8192.*omega1*omega2 + &
              11433.*omega2**2) + kd1**4*kd2**6*(19591.*omega1**2 - 13701.*omega1*omega2 + 17820.*omega2**2) + kd1**6*kd2**4*(17820.*omega1**2 - 13701.*omega1*omega2 + &
              19591.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(kd1**14*omega2*(15.*omega1 + omega2) + kd2**14*omega1*(omega1 + 15.*omega2) + 2.*kd1**11*kd2**3*(165.*omega1**2 - &
              20913.*omega1*omega2 - 4121.*omega2**2) - 2.*kd1**3*kd2**11*(4121.*omega1**2 + 20913.*omega1*omega2 - 165.*omega2**2) + kd1**4*kd2**10*(17005.*omega1**2 + &
              164123.*omega1*omega2 - 35.*omega2**2) - 2.*kd1*kd2**13*(111.*omega1**2 + 73.*omega1*omega2 + 4.*omega2**2) + kd1**2*kd2**12*(1437.*omega1**2 + 5367.*omega1*omega2 + &
              53.*omega2**2) - 2.*kd1**13*kd2*(4.*omega1**2 + 73.*omega1*omega2 + 111.*omega2**2) + kd1**12*kd2**2*(53.*omega1**2 + 5367.*omega1*omega2 + 1437.*omega2**2) + &
              kd1**8*kd2**6*(27777.*omega1**2 + 743187.*omega1*omega2 + 1951.*omega2**2) - 2.*kd1**9*kd2**5*(5831.*omega1**2 + 210135.*omega1*omega2 + 5411.*omega2**2) - &
              2.*kd1**5*kd2**9*(5411.*omega1**2 + 210135.*omega1*omega2 + 5831.*omega2**2) - 2.*kd1**7*kd2**7*(8765.*omega1**2 + 450483.*omega1*omega2 + 8765.*omega2**2) + &
              kd1**10*kd2**4*(-35.*omega1**2 + 164123.*omega1*omega2 + 17005.*omega2**2) + kd1**6*kd2**8*(1951.*omega1**2 + 743187.*omega1*omega2 + &
              27777.*omega2**2))*swd**2) + 16777216.*swd*(-32.*grav**2*kd1*(kd1 - kd2)**2*kd2*(-15.*kd2**13*omega1**2 + 15.*kd1**13*omega2**2 + 3.*kd1*kd2**12*omega1*(27.*omega1 + &
              5.*omega2) - 3.*kd1**12*kd2*omega2*(5.*omega1 + 27.*omega2) + kd1**8*kd2**5*(-10229.*omega1**2 + 6815.*omega1*omega2 - 18857.*omega2**2) + kd1**4*kd2**9*(-12055.*omega1**2 + &
              2839.*omega1*omega2 - 4955.*omega2**2) - 2.*kd1**2*kd2**11*(406.*omega1**2 - 141.*omega1*omega2 + 176.*omega2**2) + 2.*kd1**11*kd2**2*(176.*omega1**2 - 141.*omega1*omega2 + &
              406.*omega2**2) + 2.*kd1**3*kd2**10*(2467.*omega1**2 - 877.*omega1*omega2 + 1110.*omega2**2) - 2.*kd1**10*kd2**3*(1110.*omega1**2 - 877.*omega1*omega2 + &
              2467.*omega2**2) - 2.*kd1**6*kd2**7*(11312.*omega1**2 - 7444.*omega1*omega2 + 9391.*omega2**2) + kd1**5*kd2**8*(18857.*omega1**2 - 6815.*omega1*omega2 + &
              10229.*omega2**2) + 2.*kd1**7*kd2**6*(9391.*omega1**2 - 7444.*omega1*omega2 + 11312.*omega2**2) + kd1**9*kd2**4*(4955.*omega1**2 - 2839.*omega1*omega2 + &
              12055.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(-4.*kd2**15*omega1**2 + 4.*kd1**15*omega2**2 - kd1**14*kd2*omega2*(379.*omega1 + 46.*omega2) + &
              kd1*kd2**14*omega1*(46.*omega1 + 379.*omega2) + kd1**9*kd2**6*(-16812.*omega1**2 + 642647.*omega1*omega2 - 92864.*omega2**2) + kd1**5*kd2**10*(-71086.*omega1**2 + &
              476883.*omega1*omega2 - 23582.*omega2**2) + kd1**11*kd2**4*(-7324.*omega1**2 + 175655.*omega1*omega2 - 19932.*omega2**2) + kd1**7*kd2**8*(-41300.*omega1**2 + &
              553955.*omega1*omega2 - 678.*omega2**2) - kd1**12*kd2**3*(250.*omega1**2 + 26907.*omega1*omega2 + 14.*omega2**2) - 5.*kd1**2*kd2**13*(112.*omega1**2 + &
              503.*omega1*omega2 + 32.*omega2**2) + 5.*kd1**13*kd2**2*(32.*omega1**2 + 503.*omega1*omega2 + 112.*omega2**2) + kd1**3*kd2**12*(14.*omega1**2 + 26907.*omega1*omega2 + &
              250.*omega2**2) + kd1**4*kd2**11*(19932.*omega1**2 - 175655.*omega1*omega2 + 7324.*omega2**2) + kd1**6*kd2**9*(92864.*omega1**2 - 642647.*omega1*omega2 + &
              16812.*omega2**2) + kd1**8*kd2**7*(678.*omega1**2 - 553955.*omega1*omega2 + 41300.*omega2**2) + kd1**10*kd2**5*(23582.*omega1**2 - 476883.*omega1*omega2 + &
              71086.*omega2**2))*swd**2) - 1073741824.*swd*(32.*grav**2*kd1*(kd1 - kd2)*kd2*(86.*kd2**12*omega1**2 + 86.*kd1**12*omega2**2 - kd1*kd2**11*omega1*(542.*omega1 + &
              35.*omega2) - kd1**11*kd2*omega2*(35.*omega1 + 542.*omega2) + kd1**7*kd2**5*(-7782.*omega1**2 + 5579.*omega1*omega2 - 10982.*omega2**2) + kd1**5*kd2**7*(-10982.*omega1**2 + &
              5579.*omega1*omega2 - 7782.*omega2**2) + kd1**9*kd2**3*(-6100.*omega1**2 + 6985.*omega1*omega2 - 6906.*omega2**2) + kd1**3*kd2**9*(-6906.*omega1**2 + 6985.*omega1*omega2 - &
              6100.*omega2**2) + kd1**2*kd2**10*(1769.*omega1**2 - 934.*omega1*omega2 + 923.*omega2**2) + kd1**10*kd2**2*(923.*omega1**2 - 934.*omega1*omega2 + 1769.*omega2**2) + &
              kd1**6*kd2**6*(5019.*omega1**2 + 1004.*omega1*omega2 + 5019.*omega2**2) + 2.*kd1**4*kd2**8*(6674.*omega1**2 - 6026.*omega1*omega2 + 5561.*omega2**2) + &
              2.*kd1**8*kd2**4*(5561.*omega1**2 - 6026.*omega1*omega2 + 6674.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(28.*kd2**13*omega1**2 - 28.*kd1**13*omega2**2 + &
              kd1**12*kd2*omega2*(2875.*omega1 + 67.*omega2) - kd1*kd2**12*omega1*(67.*omega1 + 2875.*omega2) + kd1**10*kd2**3*(5275.*omega1**2 + 33053.*omega1*omega2 - &
              8170.*omega2**2) + kd1**3*kd2**10*(8170.*omega1**2 - 33053.*omega1*omega2 - 5275.*omega2**2) + kd1**4*kd2**9*(-126.*omega1**2 + 88855.*omega1*omega2 - &
              2105.*omega2**2) - 3.*kd1**11*kd2**2*(600.*omega1**2 + 5759.*omega1*omega2 - 651.*omega2**2) + kd1**9*kd2**4*(2105.*omega1**2 - 88855.*omega1*omega2 + &
              126.*omega2**2) + 3.*kd1**2*kd2**11*(-651.*omega1**2 + 5759.*omega1*omega2 + 600.*omega2**2) + kd1**5*kd2**8*(-51501.*omega1**2 - 266911.*omega1*omega2 + &
              1534.*omega2**2) + kd1**6*kd2**7*(91685.*omega1**2 + 465929.*omega1*omega2 + 50342.*omega2**2) + kd1**8*kd2**5*(-1534.*omega1**2 + 266911.*omega1*omega2 + &
              51501.*omega2**2) - kd1**7*kd2**6*(50342.*omega1**2 + 465929.*omega1*omega2 + 91685.*omega2**2))*swd**2) - 9007199254740992.*grav*(64.*grav**2*kd1**3*(kd1 - &
              kd2)**3*kd2**3 + (kd1*kd2**6*omega1**2*(46.*omega1 - 41.*omega2)*(omega1 - 2.*omega2) + 7.*kd2**7*omega1**3*omega2 - kd1**6*kd2*(41.*omega1 - 46.*omega2)*(2.*omega1 - &
              omega2)*omega2**2 - 7.*kd1**7*omega1*omega2**3 + kd1**2*kd2**5*omega1*(-242.*omega1**3 + 671.*omega1**2*omega2 - 576.*omega1*omega2**2 + 82.*omega2**3) + &
              kd1**5*kd2**2*omega2*(-82.*omega1**3 + 576.*omega1**2*omega2 - 671.*omega1*omega2**2 + 242.*omega2**3) + kd1**4*kd2**3*(22.*omega1**4 + 285.*omega1**3*omega2 - &
              1046.*omega1**2*omega2**2 + 736.*omega1*omega2**3 - 168.*omega2**4) + kd1**3*kd2**4*(168.*omega1**4 - 736.*omega1**3*omega2 + 1046.*omega1**2*omega2**2 - &
              285.*omega1*omega2**3 - 22.*omega2**4))*swd**2) + 2048.*grav*kd1**5*(kd1 - kd2)**5*kd2**5*(896.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2 + (kd1**5*kd2*(-416.*omega1**4 + &
              617.*omega1**3*omega2 - 3376.*omega1**2*omega2**2 + 2942.*omega1*omega2**3 - 1850.*omega2**4) + kd1*kd2**5*(-1850.*omega1**4 + 2942.*omega1**3*omega2 - &
              3376.*omega1**2*omega2**2 + 617.*omega1*omega2**3 - 416.*omega2**4) + 2.*kd2**6*(220.*omega1**4 - 428.*omega1**3*omega2 + 396.*omega1**2*omega2**2 - &
              77.*omega1*omega2**3 + 39.*omega2**4) + 2.*kd1**6*(39.*omega1**4 - 77.*omega1**3*omega2 + 396.*omega1**2*omega2**2 - 428.*omega1*omega2**3 + 220.*omega2**4) + &
              kd1**2*kd2**4*(3298.*omega1**4 - 6289.*omega1**3*omega2 + 8070.*omega1**2*omega2**2 - 3103.*omega1*omega2**3 + 1546.*omega2**4) - 2.*kd1**3*kd2**3*(1567.*omega1**4 - &
              3403.*omega1**3*omega2 + 5411.*omega1**2*omega2**2 - 3403.*omega1*omega2**3 + 1567.*omega2**4) + kd1**4*kd2**2*(1546.*omega1**4 - 3103.*omega1**3*omega2 + &
              8070.*omega1**2*omega2**2 - 6289.*omega1*omega2**3 + 3298.*omega2**4))*swd**2) + 131072.*grav*kd1**3*(kd1 - kd2)**3*kd2**3*(64.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(40.*kd1**4 - &
              96.*kd1**3*kd2 + 83.*kd1**2*kd2**2 - 96.*kd1*kd2**3 + 40.*kd2**4) + (kd1**5*kd2**5*(-50326.*omega1**4 + 106739.*omega1**3*omega2 - 211276.*omega1**2*omega2**2 + &
              106739.*omega1*omega2**3 - 50326.*omega2**4) + kd1**9*kd2*(-178.*omega1**4 + 342.*omega1**3*omega2 - 4892.*omega1**2*omega2**2 + 4887.*omega1*omega2**3 - &
              2790.*omega2**4) + kd1*kd2**9*(-2790.*omega1**4 + 4887.*omega1**3*omega2 - 4892.*omega1**2*omega2**2 + 342.*omega1*omega2**3 - 178.*omega2**4) + 2.*kd2**10*(200.*omega1**4 - &
              394.*omega1**3*omega2 + 354.*omega1**2*omega2**2 - 21.*omega1*omega2**3 + 11.*omega2**4) + kd1**10*(22.*omega1**4 - 42.*omega1**3*omega2 + 708.*omega1**2*omega2**2 - &
              788.*omega1*omega2**3 + 400.*omega2**4) + 2.*kd1**2*kd2**8*(6023.*omega1**4 - 11153.*omega1**3*omega2 + 12283.*omega1**2*omega2**2 - 2420.*omega1*omega2**3 + &
              1241.*omega2**4) + 2.*kd1**8*kd2**2*(1241.*omega1**4 - 2420.*omega1**3*omega2 + 12283.*omega1**2*omega2**2 - 11153.*omega1*omega2**3 + 6023.*omega2**4) - &
              2.*kd1**3*kd2**7*(16006.*omega1**4 - 31142.*omega1**3*omega2 + 39680.*omega1**2*omega2**2 - 11806.*omega1*omega2**3 + 6449.*omega2**4) - 2.*kd1**7*kd2**3*(6449.*omega1**4 - &
              11806.*omega1**3*omega2 + 39680.*omega1**2*omega2**2 - 31142.*omega1*omega2**3 + 16006.*omega2**4) + kd1**4*kd2**6*(50566.*omega1**4 - 104853.*omega1**3*omega2 + &
              164636.*omega1**2*omega2**2 - 65037.*omega1*omega2**3 + 32670.*omega2**4) + kd1**6*kd2**4*(32670.*omega1**4 - 65037.*omega1**3*omega2 + 164636.*omega1**2*omega2**2 - &
              104853.*omega1*omega2**3 + 50566.*omega2**4))*swd**2) + 536870912.*grav*(64.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(41.*kd1**8 - 226.*kd1**7*kd2 + 304.*kd1**6*kd2**2 - &
              334.*kd1**5*kd2**3 + 407.*kd1**4*kd2**4 - 334.*kd1**3*kd2**5 + 304.*kd1**2*kd2**6 - 226.*kd1*kd2**7 + 41.*kd2**8) + (5.*kd2**15*omega1**3*omega2 - &
              5.*kd1**15*omega1*omega2**3 + kd1*kd2**14*omega1**2*(-326.*omega1**2 + 617.*omega1*omega2 - 666.*omega2**2) + kd1**14*kd2*omega2**2*(666.*omega1**2 - 617.*omega1*omega2 + &
              326.*omega2**2) + kd1**2*kd2**13*omega1*(1786.*omega1**3 - 3383.*omega1**2*omega2 + 4096.*omega1*omega2**2 + 430.*omega2**3) - kd1**13*kd2**2*omega2*(430.*omega1**3 + &
              4096.*omega1**2*omega2 - 3383.*omega1*omega2**2 + 1786.*omega2**3) + kd1**9*kd2**6*(-41614.*omega1**4 + 91797.*omega1**3*omega2 - 417434.*omega1**2*omega2**2 + &
              327243.*omega1*omega2**3 - 159854.*omega2**4) + kd1**5*kd2**10*(-115824.*omega1**4 + 257033.*omega1**3*omega2 - 314300.*omega1**2*omega2**2 + 110651.*omega1*omega2**3 - &
              45614.*omega2**4) + kd1**3*kd2**12*(-6316.*omega1**4 + 13486.*omega1**3*omega2 - 20074.*omega1**2*omega2**2 + 4413.*omega1*omega2**3 - 2910.*omega2**4) + &
              4.*kd1**4*kd2**11*(9657.*omega1**4 - 22078.*omega1**3*omega2 + 28800.*omega1**2*omega2**2 - 11668.*omega1*omega2**3 + 5232.*omega2**4) + kd1**12*kd2**3*(2910.*omega1**4 - &
              4413.*omega1**3*omega2 + 20074.*omega1**2*omega2**2 - 13486.*omega1*omega2**3 + 6316.*omega2**4) - 4.*kd1**11*kd2**4*(5232.*omega1**4 - 11668.*omega1**3*omega2 + &
              28800.*omega1**2*omega2**2 - 22078.*omega1*omega2**3 + 9657.*omega2**4) - 2.*kd1**7*kd2**8*(52523.*omega1**4 - 87878.*omega1**3*omega2 + 171459.*omega1**2*omega2**2 - &
              25473.*omega1*omega2**3 + 20634.*omega2**4) + kd1**6*kd2**9*(159854.*omega1**4 - 327243.*omega1**3*omega2 + 417434.*omega1**2*omega2**2 - 91797.*omega1*omega2**3 + &
              41614.*omega2**4) + 2.*kd1**8*kd2**7*(20634.*omega1**4 - 25473.*omega1**3*omega2 + 171459.*omega1**2*omega2**2 - 87878.*omega1*omega2**3 + 52523.*omega2**4) + &
              kd1**10*kd2**5*(45614.*omega1**4 - 110651.*omega1**3*omega2 + 314300.*omega1**2*omega2**2 - 257033.*omega1*omega2**3 + 115824.*omega2**4))*swd**2) + &
              8388608.*grav*kd1*(kd1 - kd2)*kd2*(128.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(5.*kd1**8 - 20.*kd1**7*kd2 + 52.*kd1**6*kd2**2 - 201.*kd1**5*kd2**3 + 271.*kd1**4*kd2**4 - &
              201.*kd1**3*kd2**5 + 52.*kd1**2*kd2**6 - 20.*kd1*kd2**7 + 5.*kd2**8) + (2.*kd1**14*omega2**2*(21.*omega1**2 - 23.*omega1*omega2 + 12.*omega2**2) + &
              2.*kd2**14*omega1**2*(12.*omega1**2 - 23.*omega1*omega2 + 21.*omega2**2) - kd1*kd2**13*omega1*(222.*omega1**3 - 264.*omega1**2*omega2 + 244.*omega1*omega2**2 + &
              31.*omega2**3) - kd1**13*kd2*omega2*(31.*omega1**3 + 244.*omega1**2*omega2 - 264.*omega1*omega2**2 + 222.*omega2**3) + kd1**7*kd2**7*(-180196.*omega1**4 + &
              360577.*omega1**3*omega2 - 688730.*omega1**2*omega2**2 + 360577.*omega1*omega2**3 - 180196.*omega2**4) + kd1**11*kd2**3*(-3452.*omega1**4 + 7156.*omega1**3*omega2 - &
              57132.*omega1**2*omega2**2 + 46313.*omega1*omega2**3 - 23862.*omega2**4) + kd1**3*kd2**11*(-23862.*omega1**4 + 46313.*omega1**3*omega2 - 57132.*omega1**2*omega2**2 + &
              7156.*omega1*omega2**3 - 3452.*omega2**4) + kd1**2*kd2**12*(4014.*omega1**4 - 7047.*omega1**3*omega2 + 7574.*omega1**2*omega2**2 - 737.*omega1*omega2**3 + &
              454.*omega2**4) + kd1**12*kd2**2*(454.*omega1**4 - 737.*omega1**3*omega2 + 7574.*omega1**2*omega2**2 - 7047.*omega1*omega2**3 + 4014.*omega2**4) + &
              kd1**4*kd2**10*(62720.*omega1**4 - 128719.*omega1**3*omega2 + 188228.*omega1**2*omega2**2 - 34681.*omega1*omega2**3 + 16670.*omega2**4) - &
              2.*kd1**5*kd2**9*(53335.*omega1**4 - 112609.*omega1**3*omega2 + 194244.*omega1**2*omega2**2 - 61562.*omega1*omega2**3 + 31331.*omega2**4) - &
              2.*kd1**9*kd2**5*(31331.*omega1**4 - 61562.*omega1**3*omega2 + 194244.*omega1**2*omega2**2 - 112609.*omega1*omega2**3 + 53335.*omega2**4) + &
              kd1**10*kd2**4*(16670.*omega1**4 - 34681.*omega1**3*omega2 + 188228.*omega1**2*omega2**2 - 128719.*omega1*omega2**3 + 62720.*omega2**4) + &
              kd1**6*kd2**8*(153980.*omega1**4 - 319963.*omega1**3*omega2 + 594382.*omega1**2*omega2**2 - 271423.*omega1*omega2**3 + 139200.*omega2**4) + &
              kd1**8*kd2**6*(139200.*omega1**4 - 271423.*omega1**3*omega2 + 594382.*omega1**2*omega2**2 - 319963.*omega1*omega2**3 + 153980.*omega2**4))*swd**2))/ &
              ((-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*(omega1 - omega2)*omega2*swd**3*(-64.*grav*(64. + (kd1 - kd2)**2)*(4096. + &
              (384. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2 + (16777216. + (7340032. + (286720. + (1792. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - &
              kd2)**2)*(kd1 - kd2)**2)*(omega1 - omega2)**2*swd))
    !
end function velsb34
!
real function velsp34()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic velocity of 3rd layer for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp34')
    !
    velsp34 = -((-64.*grav*kd1**7*kd2**7*(kd1 + kd2)**7*((14.*kd1**2 + 33.*kd1*kd2 + 16.*kd2**2)*omega1**4 + (30.*kd1**2 + 38.*kd1*kd2 + 51.*kd2**2)*omega1**3*omega2 + &
              (61.*kd1**2 + 30.*kd1*kd2 + 61.*kd2**2)*omega1**2*omega2**2 + (51.*kd1**2 + 38.*kd1*kd2 + 30.*kd2**2)*omega1*omega2**3 + (16.*kd1**2 + 33.*kd1*kd2 + &
              14.*kd2**2)*omega2**4)*swd**2 - 576460752303423488.*grav*(kd1 + kd2)*(kd2**4*omega1**3*omega2 + kd1**4*omega1*omega2**3 - 2.*kd1*kd2**3*omega1**2*(omega1 + &
              omega2)*(3.*omega1 + 2.*omega2) - 2.*kd1**3*kd2*omega2**2*(omega1 + omega2)*(2.*omega1 + 3.*omega2) + kd1**2*kd2**2*omega1*omega2*(omega1**2 + 4.*omega1*omega2 + &
              omega2**2))*swd**2 - 3458764513820540928.*kd1*kd2*(kd1 + kd2)*omega1**2*omega2**2*(omega1 + omega2)**2*swd**3 + kd1**7*kd2**7*(kd1 + kd2)**7*omega1*omega2*(omega1 + &
              omega2)**2*(kd1**2*omega2*(21.*omega1 + 5.*omega2) + kd2**2*omega1*(5.*omega1 + 21.*omega2) + 8.*kd1*kd2*(omega1**2 + 6.*omega1*omega2 + &
              omega2**2))*swd**3 + 18014398509481984.*swd*(32.*grav**2*kd1**2*kd2**2*(kd1 + kd2)*(-7.*kd1*kd2*omega1**2 + (kd1**2 - 12.*kd1*kd2 + kd2**2)*omega1*omega2 - &
              7.*kd1*kd2*omega2**2) + omega1*omega2*(omega1 + omega2)**2*(kd2**2*(40.*kd1**3 + 91.*kd1**2*kd2 + 59.*kd1*kd2**2 + 20.*kd2**3)*omega1**2 + kd1*kd2*(kd1 + kd2)*(9.*kd1**2 + &
              182.*kd1*kd2 + 9.*kd2**2)*omega1*omega2 + kd1**2*(20.*kd1**3 + 59.*kd1**2*kd2 + 91.*kd1*kd2**2 + 40.*kd2**3)*omega2**2)*swd**2) + &
              281474976710656.*swd*(32.*grav**2*kd1*kd2*(kd1 + kd2)*(21.*kd2**6*omega1**2 + 21.*kd1**6*omega2**2 + kd1*kd2**5*omega1*(98.*omega1 + 13.*omega2) + &
              kd1**5*kd2*omega2*(13.*omega1 + 98.*omega2) + kd1**2*kd2**4*(135.*omega1**2 + 95.*omega1*omega2 + 8.*omega2**2) + 4.*kd1**3*kd2**3*(31.*omega1**2 + 70.*omega1*omega2 + &
              31.*omega2**2) + kd1**4*kd2**2*(8.*omega1**2 + 95.*omega1*omega2 + 135.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(2.*kd2**2*(80.*kd1**5 + 691.*kd1**4*kd2 + &
              938.*kd1**3*kd2**2 + 396.*kd1**2*kd2**3 + 17.*kd1*kd2**4 - 10.*kd2**5)*omega1**2 + kd1*kd2*(kd1 + kd2)*(753.*kd1**4 + 2600.*kd1**3*kd2 + 1793.*kd1**2*kd2**2 + &
              2600.*kd1*kd2**3 + 753.*kd2**4)*omega1*omega2 + 2.*kd1**2*(-10.*kd1**5 + 17.*kd1**4*kd2 + 396.*kd1**3*kd2**2 + 938.*kd1**2*kd2**3 + 691.*kd1*kd2**4 + &
              80.*kd2**5)*omega2**2)*swd**2) + 4398046511104.*swd*(32.*grav**2*kd1*kd2*(kd1 + kd2)*(106.*kd2**8*omega1**2 + 2.*kd1*kd2**7*omega1*(269.*omega1 - 15.*omega2) + &
              106.*kd1**8*omega2**2 + 2.*kd1**7*kd2*omega2*(-15.*omega1 + 269.*omega2) + kd1**2*kd2**6*(1060.*omega1**2 + 279.*omega1*omega2 + 201.*omega2**2) + &
              2.*kd1**3*kd2**5*(937.*omega1**2 + 850.*omega1*omega2 + 685.*omega2**2) + 2.*kd1**5*kd2**3*(685.*omega1**2 + 850.*omega1*omega2 + 937.*omega2**2) + &
              kd1**6*kd2**2*(201.*omega1**2 + 279.*omega1*omega2 + 1060.*omega2**2) + kd1**4*kd2**4*(2187.*omega1**2 + 2326.*omega1*omega2 + 2187.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(kd2**2*(536.*kd1**7 - 3293.*kd1**6*kd2 - 10423.*kd1**5*kd2**2 - 6897.*kd1**4*kd2**3 + 775.*kd1**3*kd2**4 + 1693.*kd1**2*kd2**5 + &
              407.*kd1*kd2**6 + 56.*kd2**7)*omega1**2 - kd1*kd2*(kd1 + kd2)*(3181.*kd1**6 + 12342.*kd1**5*kd2 + 13221.*kd1**4*kd2**2 + 15682.*kd1**3*kd2**3 + 13221.*kd1**2*kd2**4 + &
              12342.*kd1*kd2**5 + 3181.*kd2**6)*omega1*omega2 + kd1**2*(56.*kd1**7 + 407.*kd1**6*kd2 + 1693.*kd1**5*kd2**2 + 775.*kd1**4*kd2**3 - 6897.*kd1**3*kd2**4 - &
              10423.*kd1**2*kd2**5 - 3293.*kd1*kd2**6 + 536.*kd2**7)*omega2**2)*swd**2) - 140737488355328.*grav*(128.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(2.*kd1**2 + 13.*kd1*kd2 + &
              2.*kd2**2) + (2.*kd1*kd2**3*(153.*kd1**5 + 1143.*kd1**4*kd2 + 1843.*kd1**3*kd2**2 + 1684.*kd1**2*kd2**3 + 1049.*kd1*kd2**4 + 230.*kd2**5)*omega1**4 + kd2**2*(15.*kd1**7 + &
              2208.*kd1**6*kd2 + 9199.*kd1**5*kd2**2 + 12159.*kd1**4*kd2**3 + 8583.*kd1**3*kd2**4 + 4224.*kd1**2*kd2**5 + 820.*kd1*kd2**6 - 22.*kd2**7)*omega1**3*omega2 + &
              2.*kd1*kd2*(kd1 + kd2)*(453.*kd1**6 + 1870.*kd1**5*kd2 + 3441.*kd1**4*kd2**2 + 5059.*kd1**3*kd2**3 + 3441.*kd1**2*kd2**4 + 1870.*kd1*kd2**5 + &
              453.*kd2**6)*omega1**2*omega2**2 + kd1**2*(-22.*kd1**7 + 820.*kd1**6*kd2 + 4224.*kd1**5*kd2**2 + 8583.*kd1**4*kd2**3 + 12159.*kd1**3*kd2**4 + 9199.*kd1**2*kd2**5 + &
              2208.*kd1*kd2**6 + 15.*kd2**7)*omega1*omega2**3 + 2.*kd1**3*kd2*(230.*kd1**5 + 1049.*kd1**4*kd2 + 1684.*kd1**3*kd2**2 + 1843.*kd1**2*kd2**3 + 1143.*kd1*kd2**4 + &
              153.*kd2**5)*omega2**4)*swd**2) - 4398046511104.*grav*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(12.*kd1**4 + 67.*kd1**3*kd2 + 54.*kd1**2*kd2**2 + 67.*kd1*kd2**3 + &
              12.*kd2**4) + (kd1*kd2**3*(657.*kd1**7 + 5042.*kd1**6*kd2 + 10123.*kd1**5*kd2**2 + 9344.*kd1**4*kd2**3 + 6066.*kd1**3*kd2**4 + 4790.*kd1**2*kd2**5 + 2914.*kd1*kd2**6 + &
              566.*kd2**7)*omega1**4 + kd2**2*(-386.*kd1**9 + 1382.*kd1**8*kd2 + 13262.*kd1**7*kd2**2 + 25243.*kd1**6*kd2**3 + 20999.*kd1**5*kd2**4 + 11102.*kd1**4*kd2**5 + &
              7915.*kd1**3*kd2**6 + 5237.*kd1**2*kd2**7 + 1145.*kd1*kd2**8 + 7.*kd2**9)*omega1**3*omega2 + kd1*kd2*(kd1 + kd2)*(1146.*kd1**8 + 5022.*kd1**7*kd2 + 8808.*kd1**6*kd2**2 + &
              15973.*kd1**5*kd2**3 + 19914.*kd1**4*kd2**4 + 15973.*kd1**3*kd2**5 + 8808.*kd1**2*kd2**6 + 5022.*kd1*kd2**7 + 1146.*kd2**8)*omega1**2*omega2**2 + kd1**2*(7.*kd1**9 + &
              1145.*kd1**8*kd2 + 5237.*kd1**7*kd2**2 + 7915.*kd1**6*kd2**3 + 11102.*kd1**5*kd2**4 + 20999.*kd1**4*kd2**5 + 25243.*kd1**3*kd2**6 + 13262.*kd1**2*kd2**7 + &
              1382.*kd1*kd2**8 - 386.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(566.*kd1**7 + 2914.*kd1**6*kd2 + 4790.*kd1**5*kd2**2 + 6066.*kd1**4*kd2**3 + 9344.*kd1**3*kd2**4 + &
              10123.*kd1**2*kd2**5 + 5042.*kd1*kd2**6 + 657.*kd2**7)*omega2**4)*swd**2) - 34359738368.*grav*(128.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(25.*kd1**6 + 129.*kd1**5*kd2 + &
              59.*kd1**4*kd2**2 - 39.*kd1**3*kd2**3 + 59.*kd1**2*kd2**4 + 129.*kd1*kd2**5 + 25.*kd2**6) + (2.*kd1*kd2**3*(1607.*kd1**9 + 11915.*kd1**8*kd2 + 27193.*kd1**7*kd2**2 + &
              30775.*kd1**6*kd2**3 + 23145.*kd1**5*kd2**4 + 14306.*kd1**4*kd2**5 + 8375.*kd1**3*kd2**6 + 6142.*kd1**2*kd2**7 + 3351.*kd1*kd2**8 + 545.*kd2**9)*omega1**4 + &
              kd2**2*(-913.*kd1**11 + 4951.*kd1**10*kd2 + 52389.*kd1**9*kd2**2 + 112641.*kd1**8*kd2**3 + 101465.*kd1**7*kd2**4 + 55949.*kd1**6*kd2**5 + 39157.*kd1**5*kd2**6 + &
              33005.*kd1**4*kd2**7 + 25969.*kd1**3*kd2**8 + 13715.*kd1**2*kd2**9 + 2247.*kd1*kd2**10 + 5.*kd2**11)*omega1**3*omega2 + 2.*kd1*kd2*(kd1 + kd2)*(1080.*kd1**10 + &
              5358.*kd1**9*kd2 + 9574.*kd1**8*kd2**2 + 25595.*kd1**7*kd2**3 + 42837.*kd1**6*kd2**4 + 45840.*kd1**5*kd2**5 + 42837.*kd1**4*kd2**6 + 25595.*kd1**3*kd2**7 + &
              9574.*kd1**2*kd2**8 + 5358.*kd1*kd2**9 + 1080.*kd2**10)*omega1**2*omega2**2 + kd1**2*(5.*kd1**11 + 2247.*kd1**10*kd2 + 13715.*kd1**9*kd2**2 + 25969.*kd1**8*kd2**3 + &
              33005.*kd1**7*kd2**4 + 39157.*kd1**6*kd2**5 + 55949.*kd1**5*kd2**6 + 101465.*kd1**4*kd2**7 + 112641.*kd1**3*kd2**8 + 52389.*kd1**2*kd2**9 + 4951.*kd1*kd2**10 - &
              913.*kd2**11)*omega1*omega2**3 + 2.*kd1**3*kd2*(545.*kd1**9 + 3351.*kd1**8*kd2 + 6142.*kd1**7*kd2**2 + 8375.*kd1**6*kd2**3 + 14306.*kd1**5*kd2**4 + 23145.*kd1**4*kd2**5 + &
              30775.*kd1**3*kd2**6 + 27193.*kd1**2*kd2**7 + 11915.*kd1*kd2**8 + 1607.*kd2**9)*omega2**4)*swd**2) + 64.*kd1**5*kd2**5*(kd1 + kd2)**5*swd*(32.*grav**2*kd1**2*kd2**2*(kd1 + &
              kd2)**2*((45.*kd1**2 + 88.*kd1*kd2 + 40.*kd2**2)*omega1**2 + (53.*kd1**2 + 86.*kd1*kd2 + 53.*kd2**2)*omega1*omega2 + (40.*kd1**2 + 88.*kd1*kd2 + 45.*kd2**2)*omega2**2) + &
              omega1*omega2*(omega1 + omega2)**2*(kd1**6*omega2*(253.*omega1 + 13.*omega2) + kd2**6*omega1*(13.*omega1 + 253.*omega2) + kd1**4*kd2**2*(-63.*omega1**2 + &
              2551.*omega1*omega2 - 565.*omega2**2) + kd1**2*kd2**4*(-565.*omega1**2 + 2551.*omega1*omega2 - 63.*omega2**2) - 2.*kd1*kd2**5*(387.*omega1**2 - 601.*omega1*omega2 + &
              12.*omega2**2) + 2.*kd1**3*kd2**3*(53.*omega1**2 + 1525.*omega1*omega2 + 53.*omega2**2) - 2.*kd1**5*kd2*(12.*omega1**2 - 601.*omega1*omega2 + 387.*omega2**2))*swd**2) + &
              4096.*kd1**3*kd2**3*(kd1 + kd2)**3*swd*(32.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(71.*kd1*kd2**5*(20.*omega1**2 + 7.*omega1*omega2 + 8.*omega2**2) + &
              71.*kd1**5*kd2*(8.*omega1**2 + 7.*omega1*omega2 + 20.*omega2**2) + kd2**6*(339.*omega1**2 + 132.*omega1*omega2 + 112.*omega2**2) + 8.*kd1**3*kd2**3*(320.*omega1**2 + &
              277.*omega1*omega2 + 320.*omega2**2) + kd1**6*(112.*omega1**2 + 132.*omega1*omega2 + 339.*omega2**2) + kd1**2*kd2**4*(2555.*omega1**2 + 1479.*omega1*omega2 + &
              1538.*omega2**2) + kd1**4*kd2**2*(1538.*omega1**2 + 1479.*omega1*omega2 + 2555.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(kd1**10*(223.*omega1 - &
              17.*omega2)*omega2 + kd2**10*omega1*(-17.*omega1 + 223.*omega2) + kd1**8*kd2**2*(111.*omega1**2 + 12378.*omega1*omega2 - 196.*omega2**2) + &
              2.*kd1**3*kd2**7*(245.*omega1**2 + 24111.*omega1*omega2 - 138.*omega2**2) + 4.*kd1**9*kd2*(6.*omega1**2 + 439.*omega1*omega2 - 68.*omega2**2) + &
              4.*kd1*kd2**9*(-68.*omega1**2 + 439.*omega1*omega2 + 6.*omega2**2) + kd1**2*kd2**8*(-196.*omega1**2 + 12378.*omega1*omega2 + 111.*omega2**2) + &
              2.*kd1**7*kd2**3*(-138.*omega1**2 + 24111.*omega1*omega2 + 245.*omega2**2) + kd1**4*kd2**6*(6145.*omega1**2 + 101101.*omega1*omega2 + 876.*omega2**2) + &
              2.*kd1**5*kd2**5*(3513.*omega1**2 + 63797.*omega1*omega2 + 3513.*omega2**2) + kd1**6*kd2**4*(876.*omega1**2 + 101101.*omega1*omega2 + &
              6145.*omega2**2))*swd**2) + 68719476736.*swd*(64.*grav**2*kd1*kd2*(kd1 + kd2)*(78.*kd2**10*omega1**2 + kd1*kd2**9*omega1*(419.*omega1 - 31.*omega2) + &
              78.*kd1**10*omega2**2 + kd1**9*kd2*omega2*(-31.*omega1 + 419.*omega2) + kd1**2*kd2**8*(744.*omega1**2 + 223.*omega1*omega2 + 268.*omega2**2) + &
              kd1**8*kd2**2*(268.*omega1**2 + 223.*omega1*omega2 + 744.*omega2**2) + kd1**7*kd2**3*(1940.*omega1**2 + 1638.*omega1*omega2 + 1451.*omega2**2) + &
              kd1**3*kd2**7*(1451.*omega1**2 + 1638.*omega1*omega2 + 1940.*omega2**2) + kd1**6*kd2**4*(3593.*omega1**2 + 1834.*omega1*omega2 + 2455.*omega2**2) + &
              kd1**5*kd2**5*(3415.*omega1**2 + 1074.*omega1*omega2 + 3415.*omega2**2) + kd1**4*kd2**6*(2455.*omega1**2 + 1834.*omega1*omega2 + 3593.*omega2**2)) + &
              omega1*omega2*(omega1 + omega2)**2*(88.*kd2**11*omega1**2 + 88.*kd1**11*omega2**2 + kd1**10*kd2*omega2*(5079.*omega1 + 428.*omega2) + kd1*kd2**10*omega1*(428.*omega1 + &
              5079.*omega2) + kd1**7*kd2**4*(19896.*omega1**2 + 64727.*omega1*omega2 - 18335.*omega2**2) + kd1**8*kd2**3*(820.*omega1**2 + 49565.*omega1*omega2 - &
              7621.*omega2**2) + kd1**6*kd2**5*(24995.*omega1**2 + 95625.*omega1*omega2 - 3871.*omega2**2) + kd1**2*kd2**9*(-360.*omega1**2 + 28567.*omega1*omega2 - &
              1792.*omega2**2) + kd1**9*kd2**2*(-1792.*omega1**2 + 28567.*omega1*omega2 - 360.*omega2**2) + kd1**3*kd2**8*(-7621.*omega1**2 + 49565.*omega1*omega2 + &
              820.*omega2**2) + kd1**4*kd2**7*(-18335.*omega1**2 + 64727.*omega1*omega2 + 19896.*omega2**2) + kd1**5*kd2**6*(-3871.*omega1**2 + 95625.*omega1*omega2 + &
              24995.*omega2**2))*swd**2) + 262144.*kd1*kd2*(kd1 + kd2)*swd*(32.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(kd2**10*(230.*omega1**2 + 39.*omega1*omega2 + 35.*omega2**2) + &
              4.*kd1*kd2**9*(381.*omega1**2 + 67.*omega1*omega2 + 55.*omega2**2) + kd1**10*(35.*omega1**2 + 39.*omega1*omega2 + 230.*omega2**2) + 4.*kd1**9*kd2*(55.*omega1**2 + &
              67.*omega1*omega2 + 381.*omega2**2) + kd1**2*kd2**8*(5005.*omega1**2 + 1935.*omega1*omega2 + 1816.*omega2**2) + kd1**8*kd2**2*(1816.*omega1**2 + 1935.*omega1*omega2 + &
              5005.*omega2**2) + kd1**3*kd2**7*(11724.*omega1**2 + 7215.*omega1*omega2 + 8164.*omega2**2) + 2.*kd1**5*kd2**5*(11433.*omega1**2 + 8192.*omega1*omega2 + &
              11433.*omega2**2) + kd1**7*kd2**3*(8164.*omega1**2 + 7215.*omega1*omega2 + 11724.*omega2**2) + kd1**4*kd2**6*(19591.*omega1**2 + 13701.*omega1*omega2 + &
              17820.*omega2**2) + kd1**6*kd2**4*(17820.*omega1**2 + 13701.*omega1*omega2 + 19591.*omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(kd2**14*omega1*(omega1 - &
              15.*omega2) + kd1**14*omega2*(-15.*omega1 + omega2) - 2.*kd1**11*kd2**3*(165.*omega1**2 + 20913.*omega1*omega2 - 4121.*omega2**2) + 2.*kd1**3*kd2**11*(4121.*omega1**2 - &
              20913.*omega1*omega2 - 165.*omega2**2) + kd1**4*kd2**10*(17005.*omega1**2 - 164123.*omega1*omega2 - 35.*omega2**2) + 2.*kd1*kd2**13*(111.*omega1**2 - 73.*omega1*omega2 + &
              4.*omega2**2) + kd1**2*kd2**12*(1437.*omega1**2 - 5367.*omega1*omega2 + 53.*omega2**2) + 2.*kd1**13*kd2*(4.*omega1**2 - 73.*omega1*omega2 + 111.*omega2**2) + &
              kd1**12*kd2**2*(53.*omega1**2 - 5367.*omega1*omega2 + 1437.*omega2**2) + kd1**8*kd2**6*(27777.*omega1**2 - 743187.*omega1*omega2 + 1951.*omega2**2) + &
              2.*kd1**9*kd2**5*(5831.*omega1**2 - 210135.*omega1*omega2 + 5411.*omega2**2) + 2.*kd1**5*kd2**9*(5411.*omega1**2 - 210135.*omega1*omega2 + 5831.*omega2**2) + &
              2.*kd1**7*kd2**7*(8765.*omega1**2 - 450483.*omega1*omega2 + 8765.*omega2**2) + kd1**10*kd2**4*(-35.*omega1**2 - 164123.*omega1*omega2 + 17005.*omega2**2) + &
              kd1**6*kd2**8*(1951.*omega1**2 - 743187.*omega1*omega2 + 27777.*omega2**2))*swd**2) + 1073741824.*swd*(32.*grav**2*kd1*kd2*(kd1 + kd2)*(86.*kd2**12*omega1**2 + &
              kd1*kd2**11*omega1*(542.*omega1 - 35.*omega2) + 86.*kd1**12*omega2**2 + kd1**11*kd2*omega2*(-35.*omega1 + 542.*omega2) + kd1**2*kd2**10*(1769.*omega1**2 + &
              934.*omega1*omega2 + 923.*omega2**2) + kd1**10*kd2**2*(923.*omega1**2 + 934.*omega1*omega2 + 1769.*omega2**2) + kd1**6*kd2**6*(5019.*omega1**2 - 1004.*omega1*omega2 + &
              5019.*omega2**2) + 2.*kd1**4*kd2**8*(6674.*omega1**2 + 6026.*omega1*omega2 + 5561.*omega2**2) + kd1**3*kd2**9*(6906.*omega1**2 + 6985.*omega1*omega2 + &
              6100.*omega2**2) + 2.*kd1**8*kd2**4*(5561.*omega1**2 + 6026.*omega1*omega2 + 6674.*omega2**2) + kd1**9*kd2**3*(6100.*omega1**2 + 6985.*omega1*omega2 + &
              6906.*omega2**2) + kd1**5*kd2**7*(10982.*omega1**2 + 5579.*omega1*omega2 + 7782.*omega2**2) + kd1**7*kd2**5*(7782.*omega1**2 + 5579.*omega1*omega2 + 10982.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(28.*kd2**13*omega1**2 + kd1*kd2**12*omega1*(67.*omega1 - 2875.*omega2) + 28.*kd1**13*omega2**2 + kd1**12*kd2*omega2*(-2875.*omega1 + &
              67.*omega2) + kd1**10*kd2**3*(5275.*omega1**2 - 33053.*omega1*omega2 - 8170.*omega2**2) + kd1**5*kd2**8*(51501.*omega1**2 - 266911.*omega1*omega2 - &
              1534.*omega2**2) + 3.*kd1**11*kd2**2*(600.*omega1**2 - 5759.*omega1*omega2 - 651.*omega2**2) - 3.*kd1**2*kd2**11*(651.*omega1**2 + 5759.*omega1*omega2 - &
              600.*omega2**2) - kd1**9*kd2**4*(2105.*omega1**2 + 88855.*omega1*omega2 + 126.*omega2**2) - kd1**4*kd2**9*(126.*omega1**2 + 88855.*omega1*omega2 + &
              2105.*omega2**2) + kd1**3*kd2**10*(-8170.*omega1**2 - 33053.*omega1*omega2 + 5275.*omega2**2) + kd1**6*kd2**7*(91685.*omega1**2 - 465929.*omega1*omega2 + &
              50342.*omega2**2) + kd1**8*kd2**5*(-1534.*omega1**2 - 266911.*omega1*omega2 + 51501.*omega2**2) + kd1**7*kd2**6*(50342.*omega1**2 - 465929.*omega1*omega2 + &
              91685.*omega2**2))*swd**2) + 16777216.*swd*(32.*grav**2*kd1*kd2*(kd1 + kd2)**2*(15.*kd2**13*omega1**2 + 3.*kd1*kd2**12*omega1*(27.*omega1 - 5.*omega2) + &
              15.*kd1**13*omega2**2 + 3.*kd1**12*kd2*omega2*(-5.*omega1 + 27.*omega2) + 2.*kd1**2*kd2**11*(406.*omega1**2 + 141.*omega1*omega2 + 176.*omega2**2) + &
              2.*kd1**11*kd2**2*(176.*omega1**2 + 141.*omega1*omega2 + 406.*omega2**2) + 2.*kd1**3*kd2**10*(2467.*omega1**2 + 877.*omega1*omega2 + 1110.*omega2**2) + &
              2.*kd1**10*kd2**3*(1110.*omega1**2 + 877.*omega1*omega2 + 2467.*omega2**2) + kd1**4*kd2**9*(12055.*omega1**2 + 2839.*omega1*omega2 + 4955.*omega2**2) + &
              2.*kd1**6*kd2**7*(11312.*omega1**2 + 7444.*omega1*omega2 + 9391.*omega2**2) + kd1**5*kd2**8*(18857.*omega1**2 + 6815.*omega1*omega2 + 10229.*omega2**2) + &
              2.*kd1**7*kd2**6*(9391.*omega1**2 + 7444.*omega1*omega2 + 11312.*omega2**2) + kd1**9*kd2**4*(4955.*omega1**2 + 2839.*omega1*omega2 + 12055.*omega2**2) + &
              kd1**8*kd2**5*(10229.*omega1**2 + 6815.*omega1*omega2 + 18857.*omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(4.*kd2**15*omega1**2 + kd1*kd2**14*omega1*(46.*omega1 - &
              379.*omega2) + 4.*kd1**15*omega2**2 + kd1**14*kd2*omega2*(-379.*omega1 + 46.*omega2) + kd1**12*kd2**3*(250.*omega1**2 - 26907.*omega1*omega2 + 14.*omega2**2) + &
              5.*kd1**2*kd2**13*(112.*omega1**2 - 503.*omega1*omega2 + 32.*omega2**2) + 5.*kd1**13*kd2**2*(32.*omega1**2 - 503.*omega1*omega2 + 112.*omega2**2) + &
              kd1**3*kd2**12*(14.*omega1**2 - 26907.*omega1*omega2 + 250.*omega2**2) - kd1**7*kd2**8*(41300.*omega1**2 + 553955.*omega1*omega2 + 678.*omega2**2) - &
              kd1**4*kd2**11*(19932.*omega1**2 + 175655.*omega1*omega2 + 7324.*omega2**2) - kd1**6*kd2**9*(92864.*omega1**2 + 642647.*omega1*omega2 + 16812.*omega2**2) - &
              kd1**11*kd2**4*(7324.*omega1**2 + 175655.*omega1*omega2 + 19932.*omega2**2) - kd1**5*kd2**10*(71086.*omega1**2 + 476883.*omega1*omega2 + 23582.*omega2**2) - &
              kd1**8*kd2**7*(678.*omega1**2 + 553955.*omega1*omega2 + 41300.*omega2**2) - kd1**10*kd2**5*(23582.*omega1**2 + 476883.*omega1*omega2 + 71086.*omega2**2) - &
              kd1**9*kd2**6*(16812.*omega1**2 + 642647.*omega1*omega2 + 92864.*omega2**2))*swd**2) + 9007199254740992.*grav*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3 - &
              (7.*kd2**7*omega1**3*omega2 + 7.*kd1**7*omega1*omega2**3 + kd1*kd2**6*omega1**2*(omega1 + 2.*omega2)*(46.*omega1 + 41.*omega2) + kd1**6*kd2*omega2**2*(2.*omega1 + &
              omega2)*(41.*omega1 + 46.*omega2) + kd1**2*kd2**5*omega1*(242.*omega1**3 + 671.*omega1**2*omega2 + 576.*omega1*omega2**2 + 82.*omega2**3) + &
              kd1**5*kd2**2*omega2*(82.*omega1**3 + 576.*omega1**2*omega2 + 671.*omega1*omega2**2 + 242.*omega2**3) + kd1**3*kd2**4*(168.*omega1**4 + 736.*omega1**3*omega2 + &
              1046.*omega1**2*omega2**2 + 285.*omega1*omega2**3 - 22.*omega2**4) + kd1**4*kd2**3*(-22.*omega1**4 + 285.*omega1**3*omega2 + 1046.*omega1**2*omega2**2 + &
              736.*omega1*omega2**3 + 168.*omega2**4))*swd**2) - 2048.*grav*kd1**5*kd2**5*(kd1 + kd2)**5*(896.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4 + (2.*kd2**6*(220.*omega1**4 + &
              428.*omega1**3*omega2 + 396.*omega1**2*omega2**2 + 77.*omega1*omega2**3 + 39.*omega2**4) + 2.*kd1**6*(39.*omega1**4 + 77.*omega1**3*omega2 + 396.*omega1**2*omega2**2 + &
              428.*omega1*omega2**3 + 220.*omega2**4) + kd1*kd2**5*(1850.*omega1**4 + 2942.*omega1**3*omega2 + 3376.*omega1**2*omega2**2 + 617.*omega1*omega2**3 + &
              416.*omega2**4) + kd1**2*kd2**4*(3298.*omega1**4 + 6289.*omega1**3*omega2 + 8070.*omega1**2*omega2**2 + 3103.*omega1*omega2**3 + 1546.*omega2**4) + &
              2.*kd1**3*kd2**3*(1567.*omega1**4 + 3403.*omega1**3*omega2 + 5411.*omega1**2*omega2**2 + 3403.*omega1*omega2**3 + 1567.*omega2**4) + kd1**5*kd2*(416.*omega1**4 + &
              617.*omega1**3*omega2 + 3376.*omega1**2*omega2**2 + 2942.*omega1*omega2**3 + 1850.*omega2**4) + kd1**4*kd2**2*(1546.*omega1**4 + 3103.*omega1**3*omega2 + &
              8070.*omega1**2*omega2**2 + 6289.*omega1*omega2**3 + 3298.*omega2**4))*swd**2) - 131072.*grav*kd1**3*kd2**3*(kd1 + kd2)**3*(64.*grav**2*kd1**2*kd2**2*(kd1 + &
              kd2)**4*(40.*kd1**4 + 96.*kd1**3*kd2 + 83.*kd1**2*kd2**2 + 96.*kd1*kd2**3 + 40.*kd2**4) + (2.*kd2**10*(200.*omega1**4 + 394.*omega1**3*omega2 + 354.*omega1**2*omega2**2 + &
              21.*omega1*omega2**3 + 11.*omega2**4) + kd1*kd2**9*(2790.*omega1**4 + 4887.*omega1**3*omega2 + 4892.*omega1**2*omega2**2 + 342.*omega1*omega2**3 + 178.*omega2**4) + &
              kd1**10*(22.*omega1**4 + 42.*omega1**3*omega2 + 708.*omega1**2*omega2**2 + 788.*omega1*omega2**3 + 400.*omega2**4) + 2.*kd1**2*kd2**8*(6023.*omega1**4 + &
              11153.*omega1**3*omega2 + 12283.*omega1**2*omega2**2 + 2420.*omega1*omega2**3 + 1241.*omega2**4) + kd1**9*kd2*(178.*omega1**4 + 342.*omega1**3*omega2 + &
              4892.*omega1**2*omega2**2 + 4887.*omega1*omega2**3 + 2790.*omega2**4) + 2.*kd1**8*kd2**2*(1241.*omega1**4 + 2420.*omega1**3*omega2 + 12283.*omega1**2*omega2**2 + &
              11153.*omega1*omega2**3 + 6023.*omega2**4) + 2.*kd1**3*kd2**7*(16006.*omega1**4 + 31142.*omega1**3*omega2 + 39680.*omega1**2*omega2**2 + 11806.*omega1*omega2**3 + &
              6449.*omega2**4) + 2.*kd1**7*kd2**3*(6449.*omega1**4 + 11806.*omega1**3*omega2 + 39680.*omega1**2*omega2**2 + 31142.*omega1*omega2**3 + 16006.*omega2**4) + &
              kd1**4*kd2**6*(50566.*omega1**4 + 104853.*omega1**3*omega2 + 164636.*omega1**2*omega2**2 + 65037.*omega1*omega2**3 + 32670.*omega2**4) + kd1**5*kd2**5*(50326.*omega1**4 + &
              106739.*omega1**3*omega2 + 211276.*omega1**2*omega2**2 + 106739.*omega1*omega2**3 + 50326.*omega2**4) + kd1**6*kd2**4*(32670.*omega1**4 + 65037.*omega1**3*omega2 + &
              164636.*omega1**2*omega2**2 + 104853.*omega1*omega2**3 + 50566.*omega2**4))*swd**2) - 536870912.*grav*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(41.*kd1**8 + 226.*kd1**7*kd2 + &
              304.*kd1**6*kd2**2 + 334.*kd1**5*kd2**3 + 407.*kd1**4*kd2**4 + 334.*kd1**3*kd2**5 + 304.*kd1**2*kd2**6 + 226.*kd1*kd2**7 + 41.*kd2**8) + (-5.*kd2**15*omega1**3*omega2 - &
              5.*kd1**15*omega1*omega2**3 + kd1**14*kd2*omega2**2*(666.*omega1**2 + 617.*omega1*omega2 + 326.*omega2**2) + kd1*kd2**14*omega1**2*(326.*omega1**2 + 617.*omega1*omega2 + &
              666.*omega2**2) + kd1**2*kd2**13*omega1*(1786.*omega1**3 + 3383.*omega1**2*omega2 + 4096.*omega1*omega2**2 - 430.*omega2**3) + kd1**13*kd2**2*omega2*(-430.*omega1**3 + &
              4096.*omega1**2*omega2 + 3383.*omega1*omega2**2 + 1786.*omega2**3) + kd1**3*kd2**12*(6316.*omega1**4 + 13486.*omega1**3*omega2 + 20074.*omega1**2*omega2**2 + &
              4413.*omega1*omega2**3 + 2910.*omega2**4) + 4.*kd1**4*kd2**11*(9657.*omega1**4 + 22078.*omega1**3*omega2 + 28800.*omega1**2*omega2**2 + 11668.*omega1*omega2**3 + &
              5232.*omega2**4) + kd1**12*kd2**3*(2910.*omega1**4 + 4413.*omega1**3*omega2 + 20074.*omega1**2*omega2**2 + 13486.*omega1*omega2**3 + 6316.*omega2**4) + &
              4.*kd1**11*kd2**4*(5232.*omega1**4 + 11668.*omega1**3*omega2 + 28800.*omega1**2*omega2**2 + 22078.*omega1*omega2**3 + 9657.*omega2**4) + &
              2.*kd1**7*kd2**8*(52523.*omega1**4 + 87878.*omega1**3*omega2 + 171459.*omega1**2*omega2**2 + 25473.*omega1*omega2**3 + 20634.*omega2**4) + &
              kd1**6*kd2**9*(159854.*omega1**4 + 327243.*omega1**3*omega2 + 417434.*omega1**2*omega2**2 + 91797.*omega1*omega2**3 + 41614.*omega2**4) + &
              kd1**5*kd2**10*(115824.*omega1**4 + 257033.*omega1**3*omega2 + 314300.*omega1**2*omega2**2 + 110651.*omega1*omega2**3 + 45614.*omega2**4) + &
              2.*kd1**8*kd2**7*(20634.*omega1**4 + 25473.*omega1**3*omega2 + 171459.*omega1**2*omega2**2 + 87878.*omega1*omega2**3 + 52523.*omega2**4) + &
              kd1**10*kd2**5*(45614.*omega1**4 + 110651.*omega1**3*omega2 + 314300.*omega1**2*omega2**2 + 257033.*omega1*omega2**3 + 115824.*omega2**4) + &
              kd1**9*kd2**6*(41614.*omega1**4 + 91797.*omega1**3*omega2 + 417434.*omega1**2*omega2**2 + 327243.*omega1*omega2**3 + 159854.*omega2**4))*swd**2) - &
              8388608.*grav*kd1*kd2*(kd1 + kd2)*(128.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(5.*kd1**8 + 20.*kd1**7*kd2 + 52.*kd1**6*kd2**2 + 201.*kd1**5*kd2**3 + 271.*kd1**4*kd2**4 + &
              201.*kd1**3*kd2**5 + 52.*kd1**2*kd2**6 + 20.*kd1*kd2**7 + 5.*kd2**8) + (2.*kd1**14*omega2**2*(21.*omega1**2 + 23.*omega1*omega2 + 12.*omega2**2) + &
              2.*kd2**14*omega1**2*(12.*omega1**2 + 23.*omega1*omega2 + 21.*omega2**2) + kd1*kd2**13*omega1*(222.*omega1**3 + 264.*omega1**2*omega2 + 244.*omega1*omega2**2 - &
              31.*omega2**3) + kd1**13*kd2*omega2*(-31.*omega1**3 + 244.*omega1**2*omega2 + 264.*omega1*omega2**2 + 222.*omega2**3) + kd1**2*kd2**12*(4014.*omega1**4 + &
              7047.*omega1**3*omega2 + 7574.*omega1**2*omega2**2 + 737.*omega1*omega2**3 + 454.*omega2**4) + kd1**3*kd2**11*(23862.*omega1**4 + 46313.*omega1**3*omega2 + &
              57132.*omega1**2*omega2**2 + 7156.*omega1*omega2**3 + 3452.*omega2**4) + kd1**12*kd2**2*(454.*omega1**4 + 737.*omega1**3*omega2 + 7574.*omega1**2*omega2**2 + &
              7047.*omega1*omega2**3 + 4014.*omega2**4) + kd1**4*kd2**10*(62720.*omega1**4 + 128719.*omega1**3*omega2 + 188228.*omega1**2*omega2**2 + 34681.*omega1*omega2**3 + &
              16670.*omega2**4) + kd1**11*kd2**3*(3452.*omega1**4 + 7156.*omega1**3*omega2 + 57132.*omega1**2*omega2**2 + 46313.*omega1*omega2**3 + 23862.*omega2**4) + &
              2.*kd1**5*kd2**9*(53335.*omega1**4 + 112609.*omega1**3*omega2 + 194244.*omega1**2*omega2**2 + 61562.*omega1*omega2**3 + 31331.*omega2**4) + &
              2.*kd1**9*kd2**5*(31331.*omega1**4 + 61562.*omega1**3*omega2 + 194244.*omega1**2*omega2**2 + 112609.*omega1*omega2**3 + 53335.*omega2**4) + &
              kd1**10*kd2**4*(16670.*omega1**4 + 34681.*omega1**3*omega2 + 188228.*omega1**2*omega2**2 + 128719.*omega1*omega2**3 + 62720.*omega2**4) + &
              kd1**6*kd2**8*(153980.*omega1**4 + 319963.*omega1**3*omega2 + 594382.*omega1**2*omega2**2 + 271423.*omega1*omega2**3 + 139200.*omega2**4) + &
              kd1**8*kd2**6*(139200.*omega1**4 + 271423.*omega1**3*omega2 + 594382.*omega1**2*omega2**2 + 319963.*omega1*omega2**3 + 153980.*omega2**4) + &
              kd1**7*kd2**7*(180196.*omega1**4 + 360577.*omega1**3*omega2 + 688730.*omega1**2*omega2**2 + 360577.*omega1*omega2**3 + 180196.*omega2**4))*swd**2))/ &
              ((-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*omega2*(omega1 + omega2)*swd**3*(-64.*grav*(kd1 + kd2)**2*(64. + &
              (kd1 + kd2)**2)*(4096. + (kd1 + kd2)**2*(384. + (kd1 + kd2)**2)) + (16777216. + (kd1 + kd2)**2*(7340032. + (kd1 + kd2)**2*(286720. + (kd1 + kd2)**2*(1792. + &
              (kd1 + kd2)**2))))*(omega1 + omega2)**2*swd)))
    !
end function velsp34
!
real function velsb44()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order sub-harmonic velocity of 4th layer for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsb44')
    !
    velsb44 = -((576460752303423488.*grav*(kd1 - kd2)*(3.*kd2**4*omega1**3*omega2 + 3.*kd1**4*omega1*omega2**3 + 3.*kd1**2*kd2**2*omega1*omega2*(omega1**2 - 12.*omega1*omega2 + &
              omega2**2) + 2.*kd1**3*kd2*omega2**2*(16.*omega1**2 - 19.*omega1*omega2 + 9.*omega2**2) + 2.*kd1*kd2**3*omega1**2*(9.*omega1**2 - 19.*omega1*omega2 + 16.*omega2**2))*swd**2 - &
              64.*grav*kd1**7*(kd1 - kd2)**7*kd2**7*(3.*kd1*kd2*(omega1**2 + omega2**2)*(5.*omega1**2 - 8.*omega1*omega2 + 5.*omega2**2) + kd2**2*(-6.*omega1**4 + 21.*omega1**3*omega2 - &
              31.*omega1**2*omega2**2 + 18.*omega1*omega2**3 - 8.*omega2**4) + kd1**2*(-8.*omega1**4 + 18.*omega1**3*omega2 - 31.*omega1**2*omega2**2 + 21.*omega1*omega2**3 - &
              6.*omega2**4))*swd**2 + 1152921504606846976.*omega1*(omega1 - omega2)**2*omega2*(2.*kd2**3*omega1**2 + 11.*kd1*kd2*(-kd1 + kd2)*omega1*omega2 - 2.*kd1**3*omega2**2)*swd**3 + &
              kd1**7*(kd1 - kd2)**7*kd2**7*omega1*(omega1 - omega2)**2*omega2*(kd2**2*omega1*(3.*omega1 - 11.*omega2) + kd1**2*omega2*(-11.*omega1 + 3.*omega2) - 4.*kd1*kd2*(omega1**2 - &
              6.*omega1*omega2 + omega2**2))*swd**3 - 18014398509481984.*swd*(96.*grav**2*kd1*(kd1 - kd2)*kd2*(4.*kd2**4*omega1**2 + kd1**3*kd2*(omega1 - 8.*omega2)*omega2 + &
              4.*kd1**4*omega2**2 + kd1*kd2**3*omega1*(-8.*omega1 + omega2) + 3.*kd1**2*kd2**2*(3.*omega1**2 - 4.*omega1*omega2 + 3.*omega2**2)) - omega1*(omega1 - &
              omega2)**2*omega2*(kd2**2*(52.*kd1**3 + 9.*kd1**2*kd2 - 37.*kd1*kd2**2 + 24.*kd2**3)*omega1**2 + kd1*kd2*(-kd1 + kd2)*(303.*kd1**2 - 194.*kd1*kd2 + 303.*kd2**2)*omega1*omega2 - &
              kd1**2*(24.*kd1**3 - 37.*kd1**2*kd2 + 9.*kd1*kd2**2 + 52.*kd2**3)*omega2**2)*swd**2) - 281474976710656.*swd*(32.*grav**2*kd1*(kd1 - kd2)*kd2*(kd2**2*(372.*kd1**4 - &
              716.*kd1**3*kd2 + 561.*kd1**2*kd2**2 - 206.*kd1*kd2**3 + 79.*kd2**4)*omega1**2 - 3.*kd1*kd2*(7.*kd1**4 + 119.*kd1**3*kd2 - 192.*kd1**2*kd2**2 + 119.*kd1*kd2**3 + &
              7.*kd2**4)*omega1*omega2 + kd1**2*(79.*kd1**4 - 206.*kd1**3*kd2 + 561.*kd1**2*kd2**2 - 716.*kd1*kd2**3 + 372.*kd2**4)*omega2**2) - omega1*(omega1 - &
              omega2)**2*omega2*(2.*kd2**2*(220.*kd1**5 - 30.*kd1**4*kd2 - 626.*kd1**3*kd2**2 + 407.*kd1**2*kd2**3 + 33.*kd1*kd2**4 - 27.*kd2**5)*omega1**2 + kd1*kd2*(-kd1 + &
              kd2)*(2631.*kd1**4 - 4408.*kd1**3*kd2 + 6703.*kd1**2*kd2**2 - 4408.*kd1*kd2**3 + 2631.*kd2**4)*omega1*omega2 + 2.*kd1**2*(27.*kd1**5 - 33.*kd1**4*kd2 - 407.*kd1**3*kd2**2 + &
              626.*kd1**2*kd2**3 + 30.*kd1*kd2**4 - 220.*kd2**5)*omega2**2)*swd**2) + 9007199254740992.*grav*(576.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3 + (2.*kd1*kd2**3*(217.*kd1**3 - &
              406.*kd1**2*kd2 + 323.*kd1*kd2**2 - 141.*kd2**3)*omega1**4 - kd2**2*(30.*kd1**5 + 935.*kd1**4*kd2 - 1804.*kd1**3*kd2**2 + 1385.*kd1**2*kd2**3 - 575.*kd1*kd2**4 + &
              kd2**5)*omega1**3*omega2 + 2.*kd1*(kd1 - kd2)*kd2*(283.*kd1**4 - 497.*kd1**3*kd2 + 840.*kd1**2*kd2**2 - 497.*kd1*kd2**3 + 283.*kd2**4)*omega1**2*omega2**2 + kd1**2*(kd1**5 - &
              575.*kd1**4*kd2 + 1385.*kd1**3*kd2**2 - 1804.*kd1**2*kd2**3 + 935.*kd1*kd2**4 + 30.*kd2**5)*omega1*omega2**3 + 2.*kd1**3*kd2*(141.*kd1**3 - 323.*kd1**2*kd2 + &
              406.*kd1*kd2**2 - 217.*kd2**3)*omega2**4)*swd**2) + 140737488355328.*grav*(128.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(24.*kd1**2 - 23.*kd1*kd2 + 24.*kd2**2) + &
              (2.*kd1*kd2**3*(1621.*kd1**5 - 4659.*kd1**4*kd2 + 5749.*kd1**3*kd2**2 - 3856.*kd1**2*kd2**3 + 1709.*kd1*kd2**4 - 544.*kd2**5)*omega1**4 + kd2**2*(-823.*kd1**7 - &
              4908.*kd1**6*kd2 + 18457.*kd1**5*kd2**2 - 23733.*kd1**4*kd2**3 + 15149.*kd1**3*kd2**4 - 6356.*kd1**2*kd2**5 + 2116.*kd1*kd2**6 + 18.*kd2**7)*omega1**3*omega2 + &
              2.*kd1*(kd1 - kd2)*kd2*(1089.*kd1**6 - 2306.*kd1**5*kd2 + 7589.*kd1**4*kd2**2 - 11385.*kd1**3*kd2**3 + 7589.*kd1**2*kd2**4 - 2306.*kd1*kd2**5 + &
              1089.*kd2**6)*omega1**2*omega2**2 + kd1**2*(-18.*kd1**7 - 2116.*kd1**6*kd2 + 6356.*kd1**5*kd2**2 - 15149.*kd1**4*kd2**3 + 23733.*kd1**3*kd2**4 - 18457.*kd1**2*kd2**5 + &
              4908.*kd1*kd2**6 + 823.*kd2**7)*omega1*omega2**3 + 2.*kd1**3*kd2*(544.*kd1**5 - 1709.*kd1**4*kd2 + 3856.*kd1**3*kd2**2 - 5749.*kd1**2*kd2**3 + 4659.*kd1*kd2**4 - &
              1621.*kd2**5)*omega2**4)*swd**2) + 4398046511104.*grav*(64.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(42.*kd1**4 - 77.*kd1**3*kd2 + 118.*kd1**2*kd2**2 - 77.*kd1*kd2**3 + &
              42.*kd2**4) + (kd1*kd2**3*(2757.*kd1**7 - 10788.*kd1**6*kd2 + 19625.*kd1**5*kd2**2 - 20732.*kd1**4*kd2**3 + 13612.*kd1**3*kd2**4 - 6228.*kd1**2*kd2**5 + 2490.*kd1*kd2**6 - &
              754.*kd2**7)*omega1**4 - kd2**2*(858.*kd1**9 + 2650.*kd1**8*kd2 - 18998.*kd1**7*kd2**2 + 39793.*kd1**6*kd2**3 - 42479.*kd1**5*kd2**4 + 25712.*kd1**4*kd2**5 - &
              10391.*kd1**3*kd2**6 + 4349.*kd1**2*kd2**7 - 1539.*kd1*kd2**8 + 9.*kd2**9)*omega1**3*omega2 + kd1*(kd1 - kd2)*kd2*(1518.*kd1**8 - 3742.*kd1**7*kd2 + 13404.*kd1**6*kd2**2 - &
              32287.*kd1**5*kd2**3 + 43438.*kd1**4*kd2**4 - 32287.*kd1**3*kd2**5 + 13404.*kd1**2*kd2**6 - 3742.*kd1*kd2**7 + 1518.*kd2**8)*omega1**2*omega2**2 + kd1**2*(9.*kd1**9 - &
              1539.*kd1**8*kd2 + 4349.*kd1**7*kd2**2 - 10391.*kd1**6*kd2**3 + 25712.*kd1**5*kd2**4 - 42479.*kd1**4*kd2**5 + 39793.*kd1**3*kd2**6 - 18998.*kd1**2*kd2**7 + &
              2650.*kd1*kd2**8 + 858.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(754.*kd1**7 - 2490.*kd1**6*kd2 + 6228.*kd1**5*kd2**2 - 13612.*kd1**4*kd2**3 + 20732.*kd1**3*kd2**4 - &
              19625.*kd1**2*kd2**5 + 10788.*kd1*kd2**6 - 2757.*kd2**7)*omega2**4)*swd**2) + 34359738368.*grav*(128.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(37.*kd1**6 - 87.*kd1**5*kd2 + &
              21.*kd1**4*kd2**2 + 77.*kd1**3*kd2**3 + 21.*kd1**2*kd2**4 - 87.*kd1*kd2**5 + 37.*kd2**6) + (2.*kd1*kd2**3*(2219.*kd1**9 - 9663.*kd1**8*kd2 + 17175.*kd1**7*kd2**2 - &
              16485.*kd1**6*kd2**3 + 11489.*kd1**5*kd2**4 - 9128.*kd1**4*kd2**5 + 7685.*kd1**3*kd2**6 - 4918.*kd1**2*kd2**7 + 2111.*kd1*kd2**8 - 481.*kd2**9)*omega1**4 + &
              kd2**2*(-1179.*kd1**11 - 3473.*kd1**10*kd2 + 33427.*kd1**9*kd2**2 - 74047.*kd1**8*kd2**3 + 78251.*kd1**7*kd2**4 - 54891.*kd1**6*kd2**5 + 40311.*kd1**5*kd2**6 - &
              31023.*kd1**4*kd2**7 + 19039.*kd1**3*kd2**8 - 8345.*kd1**2*kd2**9 + 1913.*kd1*kd2**10 + kd2**11)*omega1**3*omega2 + 2.*kd1*(kd1 - kd2)*kd2*(946.*kd1**10 - 2786.*kd1**9*kd2 + &
              7044.*kd1**8*kd2**2 - 20269.*kd1**7*kd2**3 + 36725.*kd1**6*kd2**4 - 42852.*kd1**5*kd2**5 + 36725.*kd1**4*kd2**6 - 20269.*kd1**3*kd2**7 + 7044.*kd1**2*kd2**8 - &
              2786.*kd1*kd2**9 + 946.*kd2**10)*omega1**2*omega2**2 - kd1**2*(kd1**11 + 1913.*kd1**10*kd2 - 8345.*kd1**9*kd2**2 + 19039.*kd1**8*kd2**3 - 31023.*kd1**7*kd2**4 + &
              40311.*kd1**6*kd2**5 - 54891.*kd1**5*kd2**6 + 78251.*kd1**4*kd2**7 - 74047.*kd1**3*kd2**8 + 33427.*kd1**2*kd2**9 - 3473.*kd1*kd2**10 - 1179.*kd2**11)*omega1*omega2**3 + &
              2.*kd1**3*kd2*(481.*kd1**9 - 2111.*kd1**8*kd2 + 4918.*kd1**7*kd2**2 - 7685.*kd1**6*kd2**3 + 9128.*kd1**5*kd2**4 - 11489.*kd1**4*kd2**5 + 16485.*kd1**3*kd2**6 - &
              17175.*kd1**2*kd2**7 + 9663.*kd1*kd2**8 - 2219.*kd2**9)*omega2**4)*swd**2) - 64.*kd1**5*(kd1 - kd2)**5*kd2**5*swd*(32.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*((23.*kd1**2 - &
              40.*kd1*kd2 + 16.*kd2**2)*omega1**2 + (-23.*kd1**2 + 38.*kd1*kd2 - 23.*kd2**2)*omega1*omega2 + (16.*kd1**2 - 40.*kd1*kd2 + 23.*kd2**2)*omega2**2) - &
              omega1*(omega1 - omega2)**2*omega2*(5.*kd2**6*omega1*(7.*omega1 - 31.*omega2) + 5.*kd1**6*omega2*(-31.*omega1 + 7.*omega2) + 2.*kd1*kd2**5*(96.*omega1**2 + 299.*omega1*omega2 + &
              6.*omega2**2) - kd1**2*kd2**4*(89.*omega1**2 + 1089.*omega1*omega2 + 33.*omega2**2) - 2.*kd1**3*kd2**3*(45.*omega1**2 - 619.*omega1*omega2 + 45.*omega2**2) - &
              kd1**4*kd2**2*(33.*omega1**2 + 1089.*omega1*omega2 + 89.*omega2**2) + 2.*kd1**5*kd2*(6.*omega1**2 + 299.*omega1*omega2 + 96.*omega2**2))*swd**2) - &
              4096.*kd1**3*(kd1 - kd2)**3*kd2**3*swd*(32.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd1**5*kd2*(-264.*omega1**2 + 203.*omega1*omega2 - 572.*omega2**2) + &
              kd1*kd2**5*(-572.*omega1**2 + 203.*omega1*omega2 - 264.*omega2**2) + 17.*kd2**6*(9.*omega1**2 - 4.*omega1*omega2 + 4.*omega2**2) + 17.*kd1**6*(4.*omega1**2 - 4.*omega1*omega2 + &
              9.*omega2**2) - 8.*kd1**3*kd2**3*(106.*omega1**2 - 97.*omega1*omega2 + 106.*omega2**2) + kd1**2*kd2**4*(893.*omega1**2 - 525.*omega1*omega2 + 570.*omega2**2) + &
              kd1**4*kd2**2*(570.*omega1**2 - 525.*omega1*omega2 + 893.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(5.*kd1**10*(29.*omega1 - 5.*omega2)*omega2 + &
              5.*kd2**10*omega1*(-5.*omega1 + 29.*omega2) + kd1**6*kd2**4*(378.*omega1**2 + 35243.*omega1*omega2 - 3621.*omega2**2) + kd1**8*kd2**2*(-57.*omega1**2 + &
              6318.*omega1*omega2 - 826.*omega2**2) - 6.*kd1**7*kd2**3*(42.*omega1**2 + 3443.*omega1*omega2 - 354.*omega2**2) + kd1**2*kd2**8*(-826.*omega1**2 + 6318.*omega1*omega2 - &
              57.*omega2**2) + 6.*kd1**3*kd2**7*(354.*omega1**2 - 3443.*omega1*omega2 - 42.*omega2**2) + 2.*kd1**9*kd2*(6.*omega1**2 - 446.*omega1*omega2 + omega2**2) + &
              2.*kd1*kd2**9*(omega1**2 - 446.*omega1*omega2 + 6.*omega2**2) + kd1**4*kd2**6*(-3621.*omega1**2 + 35243.*omega1*omega2 + 378.*omega2**2) + 2.*kd1**5*kd2**5*(1162.*omega1**2 - &
              20215.*omega1*omega2 + 1162.*omega2**2))*swd**2) - 262144.*kd1*(kd1 - kd2)*kd2*swd*(32.*grav**2*kd1**2*(kd1 - kd2)**2*kd2**2*(kd1**7*kd2**3*(-3212.*omega1**2 + &
              2781.*omega1*omega2 - 3236.*omega2**2) + kd1**3*kd2**7*(-3236.*omega1**2 + 2781.*omega1*omega2 - 3212.*omega2**2) + 21.*kd2**10*(6.*omega1**2 - omega1*omega2 + omega2**2) + &
              21.*kd1**10*(omega1**2 - omega1*omega2 + 6.*omega2**2) - 4.*kd1*kd2**9*(171.*omega1**2 - 29.*omega1*omega2 + 27.*omega2**2) - 4.*kd1**9*kd2*(27.*omega1**2 - 29.*omega1*omega2 + &
              171.*omega2**2) + kd1**2*kd2**8*(1883.*omega1**2 - 1009.*omega1*omega2 + 940.*omega2**2) + kd1**8*kd2**2*(940.*omega1**2 - 1009.*omega1*omega2 + &
              1883.*omega2**2) - 2.*kd1**5*kd2**5*(1951.*omega1**2 - 902.*omega1*omega2 + 1951.*omega2**2) + kd1**6*kd2**4*(4612.*omega1**2 - 2771.*omega1*omega2 + &
              3561.*omega2**2) + kd1**4*kd2**6*(3561.*omega1**2 - 2771.*omega1*omega2 + 4612.*omega2**2)) - omega1*(omega1 - omega2)**2*omega2*(kd2**14*omega1*(omega1 - 9.*omega2) + &
              kd1**14*omega2*(-9.*omega1 + omega2) + kd1**10*kd2**4*(1075.*omega1**2 - 64173.*omega1*omega2 - 3467.*omega2**2) + 2.*kd1**3*kd2**11*(329.*omega1**2 + 11423.*omega1*omega2 - &
              141.*omega2**2) + kd1**2*kd2**12*(219.*omega1**2 - 3905.*omega1*omega2 - 27.*omega2**2) + 2.*kd1*kd2**13*(42.*omega1**2 + 35.*omega1*omega2 + 2.*omega2**2) + &
              2.*kd1**13*kd2*(2.*omega1**2 + 35.*omega1*omega2 + 42.*omega2**2) + kd1**12*kd2**2*(-27.*omega1**2 - 3905.*omega1*omega2 + 219.*omega2**2) + &
              2.*kd1**11*kd2**3*(-141.*omega1**2 + 11423.*omega1*omega2 + 329.*omega2**2) + kd1**4*kd2**10*(-3467.*omega1**2 - 64173.*omega1*omega2 + 1075.*omega2**2) + &
              2.*kd1**5*kd2**9*(3418.*omega1**2 + 57373.*omega1*omega2 + 2118.*omega2**2) + 10.*kd1**7*kd2**7*(2447.*omega1**2 + 16409.*omega1*omega2 + 2447.*omega2**2) + &
              2.*kd1**9*kd2**5*(2118.*omega1**2 + 57373.*omega1*omega2 + 3418.*omega2**2) - 3.*kd1**8*kd2**6*(6241.*omega1**2 + 50551.*omega1*omega2 + 5017.*omega2**2) - &
              3.*kd1**6*kd2**8*(5017.*omega1**2 + 50551.*omega1*omega2 + 6241.*omega2**2))*swd**2) - 4398046511104.*kd1*kd2*swd*(32.*grav**2*(kd1 - kd2)*(178.*kd2**8*omega1**2 + &
              178.*kd1**8*omega2**2 - 2.*kd1*kd2**7*omega1*(299.*omega1 + 45.*omega2) - 2.*kd1**7*kd2*omega2*(45.*omega1 + 299.*omega2) + kd1**2*kd2**6*(2156.*omega1**2 - &
              1017.*omega1*omega2 + 1207.*omega2**2) - 2.*kd1**3*kd2**5*(2191.*omega1**2 - 1540.*omega1*omega2 + 1763.*omega2**2) + kd1**6*kd2**2*(1207.*omega1**2 - 1017.*omega1*omega2 + &
              2156.*omega2**2) - 2.*kd1**5*kd2**3*(1763.*omega1**2 - 1540.*omega1*omega2 + 2191.*omega2**2) + kd1**4*kd2**4*(5177.*omega1**2 - 4370.*omega1*omega2 + 5177.*omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(kd1**7*omega2*(5755.*omega1 + 41.*omega2) - kd2**7*omega1*(41.*omega1 + 5755.*omega2) - kd1**4*kd2**3*(3955.*omega1**2 + &
              83073.*omega1*omega2 + 1191.*omega2**2) - kd1**6*kd2*(2844.*omega1**2 + 19149.*omega1*omega2 + 1861.*omega2**2) + kd1**5*kd2**2*(6577.*omega1**2 + 48437.*omega1*omega2 + &
              2645.*omega2**2) + kd1*kd2**6*(1861.*omega1**2 + 19149.*omega1*omega2 + 2844.*omega2**2) + kd1**3*kd2**4*(1191.*omega1**2 + 83073.*omega1*omega2 + 3955.*omega2**2) - &
              kd1**2*kd2**5*(2645.*omega1**2 + 48437.*omega1*omega2 + 6577.*omega2**2))*swd**2) - 68719476736.*swd*(64.*grav**2*kd1*kd2*(-kd1 + kd2)*(-84.*kd2**10*omega1**2 - &
              84.*kd1**10*omega2**2 + kd1*kd2**9*omega1*(309.*omega1 + 49.*omega2) + kd1**9*kd2*omega2*(49.*omega1 + 309.*omega2) + kd1**4*kd2**6*(-3597.*omega1**2 + 2228.*omega1*omega2 - &
              4153.*omega2**2) + kd1**6*kd2**4*(-4153.*omega1**2 + 2228.*omega1*omega2 - 3597.*omega2**2) + kd1**8*kd2**2*(-640.*omega1**2 + 451.*omega1*omega2 - 984.*omega2**2) + &
              kd1**2*kd2**8*(-984.*omega1**2 + 451.*omega1*omega2 - 640.*omega2**2) + 3.*kd1**5*kd2**5*(1523.*omega1**2 - 726.*omega1*omega2 + 1523.*omega2**2) + &
              kd1**7*kd2**3*(2352.*omega1**2 - 1584.*omega1*omega2 + 2173.*omega2**2) + kd1**3*kd2**7*(2173.*omega1**2 - 1584.*omega1*omega2 + 2352.*omega2**2)) - &
              omega1*(omega1 - omega2)**2*omega2*(54.*kd2**11*omega1**2 - 54.*kd1**11*omega2**2 + 37.*kd1**10*kd2*omega2*(-141.*omega1 + 4.*omega2) + 37.*kd1*kd2**10*omega1*(-4.*omega1 + &
              141.*omega2) + kd1**5*kd2**6*(-15699.*omega1**2 + 190391.*omega1*omega2 - 16237.*omega2**2) + kd1**7*kd2**4*(2836.*omega1**2 + 116305.*omega1*omega2 - &
              4141.*omega2**2) + kd1**2*kd2**9*(692.*omega1**2 - 19521.*omega1*omega2 - 3152.*omega2**2) + kd1**4*kd2**7*(4141.*omega1**2 - 116305.*omega1*omega2 - &
              2836.*omega2**2) + kd1**9*kd2**2*(3152.*omega1**2 + 19521.*omega1*omega2 - 692.*omega2**2) + kd1**8*kd2**3*(-10334.*omega1**2 - 50195.*omega1*omega2 + &
              733.*omega2**2) + kd1**3*kd2**8*(-733.*omega1**2 + 50195.*omega1*omega2 + 10334.*omega2**2) + kd1**6*kd2**5*(16237.*omega1**2 - 190391.*omega1*omega2 + &
              15699.*omega2**2))*swd**2) - 1073741824.*swd*(32.*grav**2*kd1*(kd1 - kd2)*kd2*(66.*kd2**12*omega1**2 + 66.*kd1**12*omega2**2 - kd1*kd2**11*omega1*(298.*omega1 + &
              41.*omega2) - kd1**11*kd2*omega2*(41.*omega1 + 298.*omega2) + kd1**6*kd2**6*(-10739.*omega1**2 + 8220.*omega1*omega2 - 10739.*omega2**2) + kd1**9*kd2**3*(-3804.*omega1**2 + &
              4015.*omega1*omega2 - 4054.*omega2**2) + kd1**3*kd2**9*(-4054.*omega1**2 + 4015.*omega1*omega2 - 3804.*omega2**2) + kd1**2*kd2**10*(1367.*omega1**2 - 774.*omega1*omega2 + &
              921.*omega2**2) + kd1**10*kd2**2*(921.*omega1**2 - 774.*omega1*omega2 + 1367.*omega2**2) + 2.*kd1**4*kd2**8*(1938.*omega1**2 - 2572.*omega1*omega2 + &
              1865.*omega2**2) + 2.*kd1**8*kd2**4*(1865.*omega1**2 - 2572.*omega1*omega2 + 1938.*omega2**2) + kd1**7*kd2**5*(4598.*omega1**2 - 2187.*omega1*omega2 + &
              4358.*omega2**2) + kd1**5*kd2**7*(4358.*omega1**2 - 2187.*omega1*omega2 + 4598.*omega2**2)) + omega1*(omega1 - omega2)**2*omega2*(24.*kd2**13*omega1**2 - &
              24.*kd1**13*omega2**2 + kd1**12*kd2*omega2*(2189.*omega1 + 93.*omega2) - kd1*kd2**12*omega1*(93.*omega1 + 2189.*omega2) + kd1**6*kd2**7*(13831.*omega1**2 + &
              53623.*omega1*omega2 - 15178.*omega2**2) + kd1**7*kd2**6*(15178.*omega1**2 - 53623.*omega1*omega2 - 13831.*omega2**2) + kd1**3*kd2**10*(1906.*omega1**2 - &
              21139.*omega1*omega2 - 7065.*omega2**2) + kd1**10*kd2**3*(7065.*omega1**2 + 21139.*omega1*omega2 - 1906.*omega2**2) + kd1**11*kd2**2*(-1580.*omega1**2 - &
              9587.*omega1*omega2 + 815.*omega2**2) + kd1**2*kd2**11*(-815.*omega1**2 + 9587.*omega1*omega2 + 1580.*omega2**2) + kd1**4*kd2**9*(4754.*omega1**2 + 35505.*omega1*omega2 + &
              4017.*omega2**2) - kd1**9*kd2**4*(4017.*omega1**2 + 35505.*omega1*omega2 + 4754.*omega2**2) + kd1**5*kd2**8*(-19979.*omega1**2 - 48953.*omega1*omega2 + &
              16242.*omega2**2) + kd1**8*kd2**5*(-16242.*omega1**2 + 48953.*omega1*omega2 + 19979.*omega2**2))*swd**2) + 16777216.*swd*(-32.*grav**2*kd1*(kd1 - &
              kd2)**2*kd2*(-9.*kd2**13*omega1**2 + 9.*kd1**13*omega2**2 + 3.*kd1*kd2**12*omega1*(13.*omega1 + 3.*omega2) - 3.*kd1**12*kd2*omega2*(3.*omega1 + 13.*omega2) + &
              kd1**4*kd2**9*(-2185.*omega1**2 + 533.*omega1*omega2 - 1329.*omega2**2) + kd1**5*kd2**8*(-2737.*omega1**2 + 1919.*omega1*omega2 - 1089.*omega2**2) - &
              2.*kd1**2*kd2**11*(302.*omega1**2 - 115.*omega1*omega2 + 134.*omega2**2) + 2.*kd1**11*kd2**2*(134.*omega1**2 - 115.*omega1*omega2 + 302.*omega2**2) + &
              2.*kd1**3*kd2**10*(1161.*omega1**2 - 449.*omega1*omega2 + 572.*omega2**2) - 2.*kd1**10*kd2**3*(572.*omega1**2 - 449.*omega1*omega2 + 1161.*omega2**2) + &
              kd1**9*kd2**4*(1329.*omega1**2 - 533.*omega1*omega2 + 2185.*omega2**2) + kd1**8*kd2**5*(1089.*omega1**2 - 1919.*omega1*omega2 + 2737.*omega2**2) + &
              2.*kd1**6*kd2**7*(3626.*omega1**2 - 2390.*omega1*omega2 + 2809.*omega2**2) - 2.*kd1**7*kd2**6*(2809.*omega1**2 - 2390.*omega1*omega2 + 3626.*omega2**2)) + &
              omega1*(omega1 - omega2)**2*omega2*(-2.*kd2**15*omega1**2 + 2.*kd1**15*omega2**2 + kd1**14*kd2*omega2*(-269.*omega1 + 14.*omega2) + kd1*kd2**14*omega1*(-14.*omega1 + &
              269.*omega2) + kd1**5*kd2**10*(-11858.*omega1**2 + 117173.*omega1*omega2 - 21720.*omega2**2) + kd1**11*kd2**4*(-4132.*omega1**2 + 84537.*omega1*omega2 - &
              2226.*omega2**2) - kd1**2*kd2**13*(386.*omega1**2 + 1461.*omega1*omega2 + 120.*omega2**2) + kd1**13*kd2**2*(120.*omega1**2 + 1461.*omega1*omega2 + 386.*omega2**2) + &
              kd1**3*kd2**12*(1354.*omega1**2 + 18101.*omega1*omega2 + 520.*omega2**2) + kd1**8*kd2**7*(15128.*omega1**2 + 460163.*omega1*omega2 + 612.*omega2**2) - &
              kd1**12*kd2**3*(520.*omega1**2 + 18101.*omega1*omega2 + 1354.*omega2**2) + kd1**4*kd2**11*(2226.*omega1**2 - 84537.*omega1*omega2 + 4132.*omega2**2) - &
              kd1**9*kd2**6*(32908.*omega1**2 + 105135.*omega1*omega2 + 9858.*omega2**2) + kd1**10*kd2**5*(21720.*omega1**2 - 117173.*omega1*omega2 + 11858.*omega2**2) - &
              kd1**7*kd2**8*(612.*omega1**2 + 460163.*omega1*omega2 + 15128.*omega2**2) + kd1**6*kd2**9*(9858.*omega1**2 + 105135.*omega1*omega2 + 32908.*omega2**2))*swd**2) + &
              2048.*grav*kd1**5*(kd1 - kd2)**5*kd2**5*(384.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2 + (kd1**5*kd2*(-200.*omega1**4 + 331.*omega1**3*omega2 - 1688.*omega1**2*omega2**2 + &
              1294.*omega1*omega2**3 - 702.*omega2**4) + kd1*kd2**5*(-702.*omega1**4 + 1294.*omega1**3*omega2 - 1688.*omega1**2*omega2**2 + 331.*omega1*omega2**3 - 200.*omega2**4) + &
              2.*kd2**6*(90.*omega1**4 - 212.*omega1**3*omega2 + 242.*omega1**2*omega2**2 - 55.*omega1*omega2**3 + 25.*omega2**4) + 2.*kd1**6*(25.*omega1**4 - 55.*omega1**3*omega2 + &
              242.*omega1**2*omega2**2 - 212.*omega1*omega2**3 + 90.*omega2**4) - 2.*kd1**3*kd2**3*(593.*omega1**4 - 1367.*omega1**3*omega2 + 2133.*omega1**2*omega2**2 - &
              1367.*omega1*omega2**3 + 593.*omega2**4) + kd1**2*kd2**4*(1150.*omega1**4 - 2395.*omega1**3*omega2 + 3358.*omega1**2*omega2**2 - 1441.*omega1*omega2**3 + &
              698.*omega2**4) + kd1**4*kd2**2*(698.*omega1**4 - 1441.*omega1**3*omega2 + 3358.*omega1**2*omega2**2 - 2395.*omega1*omega2**3 + 1150.*omega2**4))*swd**2) + &
              131072.*grav*kd1**3*(kd1 - kd2)**3*kd2**3*(64.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(20.*kd1**4 - 32.*kd1**3*kd2 + 9.*kd1**2*kd2**2 - 32.*kd1*kd2**3 + &
              20.*kd2**4) + (kd1**5*kd2**5*(-13306.*omega1**4 + 26913.*omega1**3*omega2 - 57172.*omega1**2*omega2**2 + 26913.*omega1*omega2**3 - 13306.*omega2**4) + &
              kd1**9*kd2*(-86.*omega1**4 + 186.*omega1**3*omega2 - 2460.*omega1**2*omega2**2 + 2161.*omega1*omega2**3 - 1098.*omega2**4) + kd1*kd2**9*(-1098.*omega1**4 + &
              2161.*omega1**3*omega2 - 2460.*omega1**2*omega2**2 + 186.*omega1*omega2**3 - 86.*omega2**4) + 2.*kd2**10*(90.*omega1**4 - 210.*omega1**3*omega2 + 218.*omega1**2*omega2**2 - &
              15.*omega1*omega2**3 + 7.*omega2**4) + 2.*kd1**10*(7.*omega1**4 - 15.*omega1**3*omega2 + 218.*omega1**2*omega2**2 - 210.*omega1*omega2**3 + 90.*omega2**4) + &
              2.*kd1**2*kd2**8*(2355.*omega1**4 - 4819.*omega1**3*omega2 + 5837.*omega1**2*omega2**2 - 1654.*omega1*omega2**3 + 753.*omega2**4) + 2.*kd1**8*kd2**2*(753.*omega1**4 - &
              1654.*omega1**3*omega2 + 5837.*omega1**2*omega2**2 - 4819.*omega1*omega2**3 + 2355.*omega2**4) - 2.*kd1**3*kd2**7*(5646.*omega1**4 - 11504.*omega1**3*omega2 + &
              15240.*omega1**2*omega2**2 - 5538.*omega1*omega2**3 + 2919.*omega2**4) - 2.*kd1**7*kd2**3*(2919.*omega1**4 - 5538.*omega1**3*omega2 + 15240.*omega1**2*omega2**2 - &
              11504.*omega1*omega2**3 + 5646.*omega2**4) + kd1**4*kd2**6*(14642.*omega1**4 - 30155.*omega1**3*omega2 + 49484.*omega1**2*omega2**2 - 19871.*omega1*omega2**3 + &
              10578.*omega2**4) + kd1**6*kd2**4*(10578.*omega1**4 - 19871.*omega1**3*omega2 + 49484.*omega1**2*omega2**2 - 30155.*omega1*omega2**3 + 14642.*omega2**4))*swd**2) + &
              8388608.*grav*kd1*(kd1 - kd2)*kd2*(128.*grav**2*kd1**2*(kd1 - kd2)**4*kd2**2*(3.*kd1**8 - 8.*kd1**7*kd2 + 22.*kd1**6*kd2**2 - 31.*kd1**5*kd2**3 - 7.*kd1**4*kd2**4 - &
              31.*kd1**3*kd2**5 + 22.*kd1**2*kd2**6 - 8.*kd1*kd2**7 + 3.*kd2**8) + (2.*kd1**14*omega2**2*(13.*omega1**2 - 13.*omega1*omega2 + 6.*omega2**2) + &
              2.*kd2**14*omega1**2*(6.*omega1**2 - 13.*omega1*omega2 + 13.*omega2**2) - kd1*kd2**13*omega1*(90.*omega1**3 - 96.*omega1**2*omega2 + 132.*omega1*omega2**2 + &
              17.*omega2**3) - kd1**13*kd2*omega2*(17.*omega1**3 + 132.*omega1**2*omega2 - 96.*omega1*omega2**2 + 90.*omega2**3) + kd1**8*kd2**6*(20852.*omega1**4 - &
              31497.*omega1**3*omega2 - 3274.*omega1**2*omega2**2 + 40955.*omega1*omega2**3 - 14644.*omega2**4) + kd1**11*kd2**3*(-2012.*omega1**4 + 4560.*omega1**3*omega2 - &
              29740.*omega1**2*omega2**2 + 23447.*omega1*omega2**3 - 11658.*omega2**4) + kd1**3*kd2**11*(-11658.*omega1**4 + 23447.*omega1**3*omega2 - 29740.*omega1**2*omega2**2 + &
              4560.*omega1*omega2**3 - 2012.*omega2**4) + kd1**2*kd2**12*(2506.*omega1**4 - 4825.*omega1**3*omega2 + 5326.*omega1**2*omega2**2 - 699.*omega1*omega2**3 + &
              366.*omega2**4) + kd1**7*kd2**7*(412.*omega1**4 - 19393.*omega1**3*omega2 + 44658.*omega1**2*omega2**2 - 19393.*omega1*omega2**3 + 412.*omega2**4) - &
              2.*kd1**9*kd2**5*(9357.*omega1**4 - 19238.*omega1**3*omega2 + 30168.*omega1**2*omega2**2 - 1561.*omega1*omega2**3 + 1021.*omega2**4) + kd1**12*kd2**2*(366.*omega1**4 - &
              699.*omega1**3*omega2 + 5326.*omega1**2*omega2**2 - 4825.*omega1*omega2**3 + 2506.*omega2**4) + kd1**4*kd2**10*(17564.*omega1**4 - 36589.*omega1**3*omega2 + &
              65808.*omega1**2*omega2**2 - 17619.*omega1*omega2**3 + 7450.*omega2**4) - 2.*kd1**5*kd2**9*(1021.*omega1**4 - 1561.*omega1**3*omega2 + 30168.*omega1**2*omega2**2 - &
              19238.*omega1*omega2**3 + 9357.*omega2**4) + kd1**10*kd2**4*(7450.*omega1**4 - 17619.*omega1**3*omega2 + 65808.*omega1**2*omega2**2 - 36589.*omega1*omega2**3 + &
              17564.*omega2**4) + kd1**6*kd2**8*(-14644.*omega1**4 + 40955.*omega1**3*omega2 - 3274.*omega1**2*omega2**2 - 31497.*omega1*omega2**3 + &
              20852.*omega2**4))*swd**2) + 536870912.*grav*(64.*grav**2*kd1**3*(kd1 - kd2)**3*kd2**3*(35.*kd1**8 - 110.*kd1**7*kd2 + 20.*kd1**6*kd2**2 + 262.*kd1**5*kd2**3 - &
              415.*kd1**4*kd2**4 + 262.*kd1**3*kd2**5 + 20.*kd1**2*kd2**6 - 110.*kd1*kd2**7 + 35.*kd2**8) + (3.*kd2**15*omega1**3*omega2 - 3.*kd1**15*omega1*omega2**3 + &
              kd1*kd2**14*omega1**2*(-226.*omega1**2 + 443.*omega1*omega2 - 462.*omega2**2) + kd1**14*kd2*omega2**2*(462.*omega1**2 - 443.*omega1*omega2 + 226.*omega2**2) + &
              kd1**2*kd2**13*omega1*(950.*omega1**3 - 1677.*omega1**2*omega2 + 2112.*omega1*omega2**2 + 334.*omega2**3) - kd1**13*kd2**2*omega2*(334.*omega1**3 + 2112.*omega1**2*omega2 - &
              1677.*omega1*omega2**2 + 950.*omega2**3) + kd1**6*kd2**9*(-12406.*omega1**4 + 20535.*omega1**3*omega2 - 72194.*omega1**2*omega2**2 + 32473.*omega1*omega2**3 - &
              28534.*omega2**4) + kd1**5*kd2**10*(-28824.*omega1**4 + 63903.*omega1**3*omega2 - 78244.*omega1**2*omega2**2 + 42553.*omega1*omega2**3 - 13530.*omega2**4) + &
              kd1**3*kd2**12*(-4808.*omega1**4 + 8538.*omega1**3*omega2 - 12854.*omega1**2*omega2**2 + 3063.*omega1*omega2**3 - 2358.*omega2**4) + 4.*kd1**4*kd2**11*(4850.*omega1**4 - &
              9842.*omega1**3*omega2 + 13834.*omega1**2*omega2**2 - 6277.*omega1*omega2**3 + 2995.*omega2**4) + kd1**12*kd2**3*(2358.*omega1**4 - 3063.*omega1**3*omega2 + &
              12854.*omega1**2*omega2**2 - 8538.*omega1*omega2**3 + 4808.*omega2**4) - 4.*kd1**11*kd2**4*(2995.*omega1**4 - 6277.*omega1**3*omega2 + 13834.*omega1**2*omega2**2 - &
              9842.*omega1*omega2**3 + 4850.*omega2**4) + kd1**9*kd2**6*(28534.*omega1**4 - 32473.*omega1**3*omega2 + 72194.*omega1**2*omega2**2 - 20535.*omega1*omega2**3 + &
              12406.*omega2**4) + kd1**10*kd2**5*(13530.*omega1**4 - 42553.*omega1**3*omega2 + 78244.*omega1**2*omega2**2 - 63903.*omega1*omega2**3 + 28824.*omega2**4) - &
              2.*kd1**8*kd2**7*(44592.*omega1**4 - 84019.*omega1**3*omega2 + 159355.*omega1**2*omega2**2 - 83552.*omega1*omega2**3 + 41329.*omega2**4) + &
              2.*kd1**7*kd2**8*(41329.*omega1**4 - 83552.*omega1**3*omega2 + 159355.*omega1**2*omega2**2 - 84019.*omega1*omega2**3 + 44592.*omega2**4))*swd**2))/ &
              ((-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*(omega1 - omega2)*omega2*swd**3*(-64.*grav*(64. + (kd1 - kd2)**2)*(4096. + &
              (384. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2 + (16777216. + (7340032. + (286720. + (1792. + (kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - kd2)**2)*(kd1 - &
              kd2)**2)*(omega1 - omega2)**2*swd)))
    !
end function velsb44
!
real function velsp44()
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, August 2023: New subroutine
!
!   Purpose
!
!   Computes second order super-harmonic velocity of 4th layer for 4 layers
!
!   Method
!
!   The surface elevation and layer-averaged velocities of the
!   second order wave group are computed by means of the SWASH
!   transfer functions for second order sub- and super-harmonic
!   interactions between two primary wave components
!   as described in
!
!   P. Vasarmidis et al
!   A study of the non-linear properties and wave generation of the multi-layer non-hydrostatic wave model SWASH
!   Ocean Engineering, 2024
!
!   Note:
!   velocity transfer function needs to be multiplied with a1 * a2
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Source text
!
    if (ltrace) call strace (ient,'velsp44')
    !
    velsp44 = (-576460752303423488.*grav*(kd1 + kd2)*(3.*kd2**4*omega1**3*omega2 + 3.*kd1**4*omega1*omega2**3 + 3.*kd1**2*kd2**2*omega1*omega2*(omega1**2 + 12.*omega1*omega2 + &
              omega2**2) + 2.*kd1**3*kd2*omega2**2*(16.*omega1**2 + 19.*omega1*omega2 + 9.*omega2**2) + 2.*kd1*kd2**3*omega1**2*(9.*omega1**2 + 19.*omega1*omega2 + &
              16.*omega2**2))*swd**2 - 64.*grav*kd1**7*kd2**7*(kd1 + kd2)**7*(3.*kd1*kd2*(omega1**2 + omega2**2)*(5.*omega1**2 + 8.*omega1*omega2 + 5.*omega2**2) + kd1**2*(8.*omega1**4 + &
              18.*omega1**3*omega2 + 31.*omega1**2*omega2**2 + 21.*omega1*omega2**3 + 6.*omega2**4) + kd2**2*(6.*omega1**4 + 21.*omega1**3*omega2 + 31.*omega1**2*omega2**2 + &
              18.*omega1*omega2**3 + 8.*omega2**4))*swd**2 + 1152921504606846976.*omega1*omega2*(omega1 + omega2)**2*(2.*kd2**3*omega1**2 + 11.*kd1*kd2*(kd1 + kd2)*omega1*omega2 + &
              2.*kd1**3*omega2**2)*swd**3 + kd1**7*kd2**7*(kd1 + kd2)**7*omega1*omega2*(omega1 + omega2)**2*(kd1**2*omega2*(11.*omega1 + 3.*omega2) + kd2**2*omega1*(3.*omega1 + &
              11.*omega2) + 4.*kd1*kd2*(omega1**2 + 6.*omega1*omega2 + omega2**2))*swd**3 + 281474976710656.*swd*(32.*grav**2*kd1*kd2*(kd1 + kd2)*(79.*kd2**6*omega1**2 + &
              kd1*kd2**5*omega1*(206.*omega1 - 21.*omega2) + 79.*kd1**6*omega2**2 + kd1**5*kd2*omega2*(-21.*omega1 + 206.*omega2) + 3.*kd1**2*kd2**4*(187.*omega1**2 + 119.*omega1*omega2 + &
              124.*omega2**2) + 4.*kd1**3*kd2**3*(179.*omega1**2 + 144.*omega1*omega2 + 179.*omega2**2) + 3.*kd1**4*kd2**2*(124.*omega1**2 + 119.*omega1*omega2 + 187.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(2.*kd2**2*(220.*kd1**5 + 30.*kd1**4*kd2 - 626.*kd1**3*kd2**2 - 407.*kd1**2*kd2**3 + 33.*kd1*kd2**4 + 27.*kd2**5)*omega1**2 - &
              kd1*kd2*(kd1 + kd2)*(2631.*kd1**4 + 4408.*kd1**3*kd2 + 6703.*kd1**2*kd2**2 + 4408.*kd1*kd2**3 + 2631.*kd2**4)*omega1*omega2 + 2.*kd1**2*(27.*kd1**5 + 33.*kd1**4*kd2 - &
              407.*kd1**3*kd2**2 - 626.*kd1**2*kd2**3 + 30.*kd1*kd2**4 + 220.*kd2**5)*omega2**2)*swd**2) + &
              4398046511104.*kd1*kd2*swd*(32.*grav**2*(kd1 + kd2)*(178.*kd2**8*omega1**2 + 2.*kd1*kd2**7*omega1*(299.*omega1 - 45.*omega2) + 178.*kd1**8*omega2**2 + &
              2.*kd1**7*kd2*omega2*(-45.*omega1 + 299.*omega2) + kd1**2*kd2**6*(2156.*omega1**2 + 1017.*omega1*omega2 + 1207.*omega2**2) + 2.*kd1**3*kd2**5*(2191.*omega1**2 + &
              1540.*omega1*omega2 + 1763.*omega2**2) + kd1**6*kd2**2*(1207.*omega1**2 + 1017.*omega1*omega2 + 2156.*omega2**2) + 2.*kd1**5*kd2**3*(1763.*omega1**2 + 1540.*omega1*omega2 + &
              2191.*omega2**2) + kd1**4*kd2**4*(5177.*omega1**2 + 4370.*omega1*omega2 + 5177.*omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(kd2*(2844.*kd1**6 + 6577.*kd1**5*kd2 + &
              3955.*kd1**4*kd2**2 + 1191.*kd1**3*kd2**3 + 2645.*kd1**2*kd2**4 + 1861.*kd1*kd2**5 + 41.*kd2**6)*omega1**2 - (kd1 + kd2)*(5755.*kd1**6 + 13394.*kd1**5*kd2 + &
              35043.*kd1**4*kd2**2 + 48030.*kd1**3*kd2**3 + 35043.*kd1**2*kd2**4 + 13394.*kd1*kd2**5 + 5755.*kd2**6)*omega1*omega2 + kd1*(41.*kd1**6 + 1861.*kd1**5*kd2 + &
              2645.*kd1**4*kd2**2 + 1191.*kd1**3*kd2**3 + 3955.*kd1**2*kd2**4 + 6577.*kd1*kd2**5 + 2844.*kd2**6)*omega2**2)*swd**2) - 9007199254740992.*grav*(576.*grav**2*kd1**3*kd2**3*(kd1 + &
              kd2)**3 + (2.*kd1*kd2**3*(217.*kd1**3 + 406.*kd1**2*kd2 + 323.*kd1*kd2**2 + 141.*kd2**3)*omega1**4 + kd2**2*(-30.*kd1**5 + 935.*kd1**4*kd2 + 1804.*kd1**3*kd2**2 + &
              1385.*kd1**2*kd2**3 + 575.*kd1*kd2**4 + kd2**5)*omega1**3*omega2 + 2.*kd1*kd2*(kd1 + kd2)*(283.*kd1**4 + 497.*kd1**3*kd2 + 840.*kd1**2*kd2**2 + 497.*kd1*kd2**3 + &
              283.*kd2**4)*omega1**2*omega2**2 + kd1**2*(kd1**5 + 575.*kd1**4*kd2 + 1385.*kd1**3*kd2**2 + 1804.*kd1**2*kd2**3 + 935.*kd1*kd2**4 - 30.*kd2**5)*omega1*omega2**3 + &
              2.*kd1**3*kd2*(141.*kd1**3 + 323.*kd1**2*kd2 + 406.*kd1*kd2**2 + 217.*kd2**3)*omega2**4)*swd**2) - 140737488355328.*grav*(128.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(24.*kd1**2 + &
              23.*kd1*kd2 + 24.*kd2**2) + (2.*kd1*kd2**3*(1621.*kd1**5 + 4659.*kd1**4*kd2 + 5749.*kd1**3*kd2**2 + 3856.*kd1**2*kd2**3 + 1709.*kd1*kd2**4 + 544.*kd2**5)*omega1**4 + &
              kd2**2*(-823.*kd1**7 + 4908.*kd1**6*kd2 + 18457.*kd1**5*kd2**2 + 23733.*kd1**4*kd2**3 + 15149.*kd1**3*kd2**4 + 6356.*kd1**2*kd2**5 + 2116.*kd1*kd2**6 - &
              18.*kd2**7)*omega1**3*omega2 + 2.*kd1*kd2*(kd1 + kd2)*(1089.*kd1**6 + 2306.*kd1**5*kd2 + 7589.*kd1**4*kd2**2 + 11385.*kd1**3*kd2**3 + 7589.*kd1**2*kd2**4 + &
              2306.*kd1*kd2**5 + 1089.*kd2**6)*omega1**2*omega2**2 + kd1**2*(-18.*kd1**7 + 2116.*kd1**6*kd2 + 6356.*kd1**5*kd2**2 + 15149.*kd1**4*kd2**3 + 23733.*kd1**3*kd2**4 + &
              18457.*kd1**2*kd2**5 + 4908.*kd1*kd2**6 - 823.*kd2**7)*omega1*omega2**3 + 2.*kd1**3*kd2*(544.*kd1**5 + 1709.*kd1**4*kd2 + 3856.*kd1**3*kd2**2 + 5749.*kd1**2*kd2**3 + &
              4659.*kd1*kd2**4 + 1621.*kd2**5)*omega2**4)*swd**2) - 4398046511104.*grav*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(42.*kd1**4 + 77.*kd1**3*kd2 + 118.*kd1**2*kd2**2 + &
              77.*kd1*kd2**3 + 42.*kd2**4) + (kd1*kd2**3*(2757.*kd1**7 + 10788.*kd1**6*kd2 + 19625.*kd1**5*kd2**2 + 20732.*kd1**4*kd2**3 + 13612.*kd1**3*kd2**4 + 6228.*kd1**2*kd2**5 + &
              2490.*kd1*kd2**6 + 754.*kd2**7)*omega1**4 + kd2**2*(-858.*kd1**9 + 2650.*kd1**8*kd2 + 18998.*kd1**7*kd2**2 + 39793.*kd1**6*kd2**3 + 42479.*kd1**5*kd2**4 + &
              25712.*kd1**4*kd2**5 + 10391.*kd1**3*kd2**6 + 4349.*kd1**2*kd2**7 + 1539.*kd1*kd2**8 + 9.*kd2**9)*omega1**3*omega2 + kd1*kd2*(kd1 + kd2)*(1518.*kd1**8 + 3742.*kd1**7*kd2 + &
              13404.*kd1**6*kd2**2 + 32287.*kd1**5*kd2**3 + 43438.*kd1**4*kd2**4 + 32287.*kd1**3*kd2**5 + 13404.*kd1**2*kd2**6 + 3742.*kd1*kd2**7 + 1518.*kd2**8)*omega1**2*omega2**2 + &
              kd1**2*(9.*kd1**9 + 1539.*kd1**8*kd2 + 4349.*kd1**7*kd2**2 + 10391.*kd1**6*kd2**3 + 25712.*kd1**5*kd2**4 + 42479.*kd1**4*kd2**5 + 39793.*kd1**3*kd2**6 + 18998.*kd1**2*kd2**7 + &
              2650.*kd1*kd2**8 - 858.*kd2**9)*omega1*omega2**3 + kd1**3*kd2*(754.*kd1**7 + 2490.*kd1**6*kd2 + 6228.*kd1**5*kd2**2 + 13612.*kd1**4*kd2**3 + 20732.*kd1**3*kd2**4 + &
              19625.*kd1**2*kd2**5 + 10788.*kd1*kd2**6 + 2757.*kd2**7)*omega2**4)*swd**2) - 34359738368.*grav*(128.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(37.*kd1**6 + 87.*kd1**5*kd2 + &
              21.*kd1**4*kd2**2 - 77.*kd1**3*kd2**3 + 21.*kd1**2*kd2**4 + 87.*kd1*kd2**5 + 37.*kd2**6) + (2.*kd1*kd2**3*(2219.*kd1**9 + 9663.*kd1**8*kd2 + 17175.*kd1**7*kd2**2 + &
              16485.*kd1**6*kd2**3 + 11489.*kd1**5*kd2**4 + 9128.*kd1**4*kd2**5 + 7685.*kd1**3*kd2**6 + 4918.*kd1**2*kd2**7 + 2111.*kd1*kd2**8 + 481.*kd2**9)*omega1**4 + &
              kd2**2*(-1179.*kd1**11 + 3473.*kd1**10*kd2 + 33427.*kd1**9*kd2**2 + 74047.*kd1**8*kd2**3 + 78251.*kd1**7*kd2**4 + 54891.*kd1**6*kd2**5 + 40311.*kd1**5*kd2**6 + &
              31023.*kd1**4*kd2**7 + 19039.*kd1**3*kd2**8 + 8345.*kd1**2*kd2**9 + 1913.*kd1*kd2**10 - kd2**11)*omega1**3*omega2 + 2.*kd1*kd2*(kd1 + kd2)*(946.*kd1**10 + 2786.*kd1**9*kd2 + &
              7044.*kd1**8*kd2**2 + 20269.*kd1**7*kd2**3 + 36725.*kd1**6*kd2**4 + 42852.*kd1**5*kd2**5 + 36725.*kd1**4*kd2**6 + 20269.*kd1**3*kd2**7 + 7044.*kd1**2*kd2**8 + &
              2786.*kd1*kd2**9 + 946.*kd2**10)*omega1**2*omega2**2 + kd1**2*(-kd1**11 + 1913.*kd1**10*kd2 + 8345.*kd1**9*kd2**2 + 19039.*kd1**8*kd2**3 + 31023.*kd1**7*kd2**4 + &
              40311.*kd1**6*kd2**5 + 54891.*kd1**5*kd2**6 + 78251.*kd1**4*kd2**7 + 74047.*kd1**3*kd2**8 + 33427.*kd1**2*kd2**9 + 3473.*kd1*kd2**10 - 1179.*kd2**11)*omega1*omega2**3 + &
              2.*kd1**3*kd2*(481.*kd1**9 + 2111.*kd1**8*kd2 + 4918.*kd1**7*kd2**2 + 7685.*kd1**6*kd2**3 + 9128.*kd1**5*kd2**4 + 11489.*kd1**4*kd2**5 + 16485.*kd1**3*kd2**6 + &
              17175.*kd1**2*kd2**7 + 9663.*kd1*kd2**8 + 2219.*kd2**9)*omega2**4)*swd**2) + 18014398509481984.*swd*(96.*grav**2*kd1*kd2*(kd1 + kd2)*(4.*kd2**4*omega1**2 + &
              4.*kd1**4*omega2**2 + kd1*kd2**3*omega1*(8.*omega1 + omega2) + kd1**3*kd2*omega2*(omega1 + 8.*omega2) + 3.*kd1**2*kd2**2*(3.*omega1**2 + 4.*omega1*omega2 + 3.*omega2**2)) + &
              omega1*omega2*(omega1 + omega2)**2*(24.*kd2**5*omega1**2 + 24.*kd1**5*omega2**2 + kd1**4*kd2*omega2*(303.*omega1 + 37.*omega2) + kd1*kd2**4*omega1*(37.*omega1 + &
              303.*omega2) + kd1**2*kd2**3*(9.*omega1**2 + 497.*omega1*omega2 - 52.*omega2**2) + kd1**3*kd2**2*(-52.*omega1**2 + 497.*omega1*omega2 + 9.*omega2**2))*swd**2) + &
              64.*kd1**5*kd2**5*(kd1 + kd2)**5*swd*(32.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*((23.*kd1**2 + 40.*kd1*kd2 + 16.*kd2**2)*omega1**2 + (23.*kd1**2 + 38.*kd1*kd2 + &
              23.*kd2**2)*omega1*omega2 + (16.*kd1**2 + 40.*kd1*kd2 + 23.*kd2**2)*omega2**2) + omega1*omega2*(omega1 + omega2)**2*(5.*kd1**6*omega2*(31.*omega1 + 7.*omega2) + &
              5.*kd2**6*omega1*(7.*omega1 + 31.*omega2) + kd1**4*kd2**2*(-33.*omega1**2 + 1089.*omega1*omega2 - 89.*omega2**2) + kd1**2*kd2**4*(-89.*omega1**2 + 1089.*omega1*omega2 - &
              33.*omega2**2) - 2.*kd1*kd2**5*(96.*omega1**2 - 299.*omega1*omega2 + 6.*omega2**2) + 2.*kd1**3*kd2**3*(45.*omega1**2 + 619.*omega1*omega2 + 45.*omega2**2) - &
              2.*kd1**5*kd2*(6.*omega1**2 - 299.*omega1*omega2 + 96.*omega2**2))*swd**2) + 4096.*kd1**3*kd2**3*(kd1 + kd2)**3*swd*(32.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*((kd1 + &
              kd2)*(68.*kd1**5 + 196.*kd1**4*kd2 + 374.*kd1**3*kd2**2 + 474.*kd1**2*kd2**3 + 419.*kd1*kd2**4 + 153.*kd2**5)*omega1**2 + (68.*kd1**6 + 203.*kd1**5*kd2 + 525.*kd1**4*kd2**2 + &
              776.*kd1**3*kd2**3 + 525.*kd1**2*kd2**4 + 203.*kd1*kd2**5 + 68.*kd2**6)*omega1*omega2 + (kd1 + kd2)*(153.*kd1**5 + 419.*kd1**4*kd2 + 474.*kd1**3*kd2**2 + 374.*kd1**2*kd2**3 + &
              196.*kd1*kd2**4 + 68.*kd2**5)*omega2**2) + omega1*omega2*(omega1 + omega2)**2*(5.*kd1**10*omega2*(29.*omega1 + 5.*omega2) + 5.*kd2**10*omega1*(5.*omega1 + &
              29.*omega2) + kd1**4*kd2**6*(3621.*omega1**2 + 35243.*omega1*omega2 - 378.*omega2**2) + 6.*kd1**3*kd2**7*(354.*omega1**2 + 3443.*omega1*omega2 - 42.*omega2**2) + &
              2.*kd1**9*kd2*(6.*omega1**2 + 446.*omega1*omega2 + omega2**2) + 2.*kd1*kd2**9*(omega1**2 + 446.*omega1*omega2 + 6.*omega2**2) + kd1**2*kd2**8*(826.*omega1**2 + &
              6318.*omega1*omega2 + 57.*omega2**2) + 6.*kd1**7*kd2**3*(-42.*omega1**2 + 3443.*omega1*omega2 + 354.*omega2**2) + kd1**8*kd2**2*(57.*omega1**2 + 6318.*omega1*omega2 + &
              826.*omega2**2) + 2.*kd1**5*kd2**5*(1162.*omega1**2 + 20215.*omega1*omega2 + 1162.*omega2**2) + kd1**6*kd2**4*(-378.*omega1**2 + 35243.*omega1*omega2 + &
              3621.*omega2**2))*swd**2) + 262144.*kd1*kd2*(kd1 + kd2)*swd*(32.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**2*(21.*kd2**10*(6.*omega1**2 + omega1*omega2 + omega2**2) + &
              21.*kd1**10*(omega1**2 + omega1*omega2 + 6.*omega2**2) + 4.*kd1*kd2**9*(171.*omega1**2 + 29.*omega1*omega2 + 27.*omega2**2) + 4.*kd1**9*kd2*(27.*omega1**2 + 29.*omega1*omega2 + &
              171.*omega2**2) + kd1**2*kd2**8*(1883.*omega1**2 + 1009.*omega1*omega2 + 940.*omega2**2) + kd1**8*kd2**2*(940.*omega1**2 + 1009.*omega1*omega2 + &
              1883.*omega2**2) + 2.*kd1**5*kd2**5*(1951.*omega1**2 + 902.*omega1*omega2 + 1951.*omega2**2) + kd1**3*kd2**7*(3236.*omega1**2 + 2781.*omega1*omega2 + &
              3212.*omega2**2) + kd1**7*kd2**3*(3212.*omega1**2 + 2781.*omega1*omega2 + 3236.*omega2**2) + kd1**6*kd2**4*(4612.*omega1**2 + 2771.*omega1*omega2 + &
              3561.*omega2**2) + kd1**4*kd2**6*(3561.*omega1**2 + 2771.*omega1*omega2 + 4612.*omega2**2)) + omega1*omega2*(omega1 + omega2)**2*(kd1**14*omega2*(9.*omega1 + omega2) + &
              kd2**14*omega1*(omega1 + 9.*omega2) + kd1**10*kd2**4*(1075.*omega1**2 + 64173.*omega1*omega2 - 3467.*omega2**2) + 2.*kd1**11*kd2**3*(141.*omega1**2 + 11423.*omega1*omega2 - &
              329.*omega2**2) + kd1**2*kd2**12*(219.*omega1**2 + 3905.*omega1*omega2 - 27.*omega2**2) - 2.*kd1*kd2**13*(42.*omega1**2 - 35.*omega1*omega2 + 2.*omega2**2) - &
              2.*kd1**13*kd2*(2.*omega1**2 - 35.*omega1*omega2 + 42.*omega2**2) + 2.*kd1**3*kd2**11*(-329.*omega1**2 + 11423.*omega1*omega2 + 141.*omega2**2) + &
              kd1**12*kd2**2*(-27.*omega1**2 + 3905.*omega1*omega2 + 219.*omega2**2) + kd1**4*kd2**10*(-3467.*omega1**2 + 64173.*omega1*omega2 + 1075.*omega2**2) - &
              2.*kd1**5*kd2**9*(3418.*omega1**2 - 57373.*omega1*omega2 + 2118.*omega2**2) - 10.*kd1**7*kd2**7*(2447.*omega1**2 - 16409.*omega1*omega2 + 2447.*omega2**2) - &
              2.*kd1**9*kd2**5*(2118.*omega1**2 - 57373.*omega1*omega2 + 3418.*omega2**2) - 3.*kd1**8*kd2**6*(6241.*omega1**2 - 50551.*omega1*omega2 + 5017.*omega2**2) - &
              3.*kd1**6*kd2**8*(5017.*omega1**2 - 50551.*omega1*omega2 + 6241.*omega2**2))*swd**2) + 68719476736.*swd*(64.*grav**2*kd1*kd2*(kd1 + kd2)*(84.*kd2**10*omega1**2 + &
              kd1*kd2**9*omega1*(309.*omega1 - 49.*omega2) + 84.*kd1**10*omega2**2 + kd1**9*kd2*omega2*(-49.*omega1 + 309.*omega2) + kd1**2*kd2**8*(984.*omega1**2 + 451.*omega1*omega2 + &
              640.*omega2**2) + kd1**8*kd2**2*(640.*omega1**2 + 451.*omega1*omega2 + 984.*omega2**2) + 3.*kd1**5*kd2**5*(1523.*omega1**2 + 726.*omega1*omega2 + 1523.*omega2**2) + &
              kd1**7*kd2**3*(2352.*omega1**2 + 1584.*omega1*omega2 + 2173.*omega2**2) + kd1**3*kd2**7*(2173.*omega1**2 + 1584.*omega1*omega2 + 2352.*omega2**2) + &
              kd1**6*kd2**4*(4153.*omega1**2 + 2228.*omega1*omega2 + 3597.*omega2**2) + kd1**4*kd2**6*(3597.*omega1**2 + 2228.*omega1*omega2 + 4153.*omega2**2)) + &
              omega1*omega2*(omega1 + omega2)**2*(54.*kd2**11*omega1**2 + 54.*kd1**11*omega2**2 + 37.*kd1**10*kd2*omega2*(141.*omega1 + 4.*omega2) + 37.*kd1*kd2**10*omega1*(4.*omega1 + &
              141.*omega2) + kd1**3*kd2**8*(733.*omega1**2 + 50195.*omega1*omega2 - 10334.*omega2**2) + kd1**2*kd2**9*(692.*omega1**2 + 19521.*omega1*omega2 - 3152.*omega2**2) + &
              kd1**4*kd2**7*(4141.*omega1**2 + 116305.*omega1*omega2 - 2836.*omega2**2) + kd1**9*kd2**2*(-3152.*omega1**2 + 19521.*omega1*omega2 + 692.*omega2**2) + &
              kd1**8*kd2**3*(-10334.*omega1**2 + 50195.*omega1*omega2 + 733.*omega2**2) + kd1**7*kd2**4*(-2836.*omega1**2 + 116305.*omega1*omega2 + 4141.*omega2**2) + &
              kd1**6*kd2**5*(16237.*omega1**2 + 190391.*omega1*omega2 + 15699.*omega2**2) + kd1**5*kd2**6*(15699.*omega1**2 + 190391.*omega1*omega2 + &
              16237.*omega2**2))*swd**2) + 1073741824.*swd*(32.*grav**2*kd1*kd2*(kd1 + kd2)*(66.*kd2**12*omega1**2 + kd1*kd2**11*omega1*(298.*omega1 - 41.*omega2) + &
              66.*kd1**12*omega2**2 + kd1**11*kd2*omega2*(-41.*omega1 + 298.*omega2) + kd1**2*kd2**10*(1367.*omega1**2 + 774.*omega1*omega2 + 921.*omega2**2) + &
              kd1**10*kd2**2*(921.*omega1**2 + 774.*omega1*omega2 + 1367.*omega2**2) + 2.*kd1**4*kd2**8*(1938.*omega1**2 + 2572.*omega1*omega2 + 1865.*omega2**2) + &
              2.*kd1**8*kd2**4*(1865.*omega1**2 + 2572.*omega1*omega2 + 1938.*omega2**2) + kd1**3*kd2**9*(4054.*omega1**2 + 4015.*omega1*omega2 + 3804.*omega2**2) + &
              kd1**9*kd2**3*(3804.*omega1**2 + 4015.*omega1*omega2 + 4054.*omega2**2) - kd1**7*kd2**5*(4598.*omega1**2 + 2187.*omega1*omega2 + 4358.*omega2**2) - &
              kd1**5*kd2**7*(4358.*omega1**2 + 2187.*omega1*omega2 + 4598.*omega2**2) - kd1**6*kd2**6*(10739.*omega1**2 + 8220.*omega1*omega2 + 10739.*omega2**2)) - &
              omega1*omega2*(omega1 + omega2)**2*(24.*kd2**13*omega1**2 + kd1*kd2**12*omega1*(93.*omega1 - 2189.*omega2) + 24.*kd1**13*omega2**2 + kd1**12*kd2*omega2*(-2189.*omega1 + &
              93.*omega2) + kd1**5*kd2**8*(19979.*omega1**2 - 48953.*omega1*omega2 - 16242.*omega2**2) + kd1**6*kd2**7*(13831.*omega1**2 - 53623.*omega1*omega2 - &
              15178.*omega2**2) + kd1**10*kd2**3*(7065.*omega1**2 - 21139.*omega1*omega2 - 1906.*omega2**2) + kd1**11*kd2**2*(1580.*omega1**2 - 9587.*omega1*omega2 - &
              815.*omega2**2) + kd1**2*kd2**11*(-815.*omega1**2 - 9587.*omega1*omega2 + 1580.*omega2**2) + kd1**4*kd2**9*(4754.*omega1**2 - 35505.*omega1*omega2 + &
              4017.*omega2**2) + kd1**9*kd2**4*(4017.*omega1**2 - 35505.*omega1*omega2 + 4754.*omega2**2) + kd1**3*kd2**10*(-1906.*omega1**2 - 21139.*omega1*omega2 + &
              7065.*omega2**2) + kd1**7*kd2**6*(-15178.*omega1**2 - 53623.*omega1*omega2 + 13831.*omega2**2) + kd1**8*kd2**5*(-16242.*omega1**2 - 48953.*omega1*omega2 + &
              19979.*omega2**2))*swd**2) + 16777216.*swd*(32.*grav**2*kd1*kd2*(kd1 + kd2)**2*(9.*kd2**13*omega1**2 + 3.*kd1*kd2**12*omega1*(13.*omega1 - 3.*omega2) + &
              9.*kd1**13*omega2**2 + 3.*kd1**12*kd2*omega2*(-3.*omega1 + 13.*omega2) + 2.*kd1**2*kd2**11*(302.*omega1**2 + 115.*omega1*omega2 + 134.*omega2**2) + &
              2.*kd1**11*kd2**2*(134.*omega1**2 + 115.*omega1*omega2 + 302.*omega2**2) + 2.*kd1**3*kd2**10*(1161.*omega1**2 + 449.*omega1*omega2 + 572.*omega2**2) - &
              kd1**5*kd2**8*(2737.*omega1**2 + 1919.*omega1*omega2 + 1089.*omega2**2) + 2.*kd1**10*kd2**3*(572.*omega1**2 + 449.*omega1*omega2 + 1161.*omega2**2) + &
              kd1**4*kd2**9*(2185.*omega1**2 + 533.*omega1*omega2 + 1329.*omega2**2) + kd1**9*kd2**4*(1329.*omega1**2 + 533.*omega1*omega2 + 2185.*omega2**2) - &
              kd1**8*kd2**5*(1089.*omega1**2 + 1919.*omega1*omega2 + 2737.*omega2**2) - 2.*kd1**6*kd2**7*(3626.*omega1**2 + 2390.*omega1*omega2 + 2809.*omega2**2) - &
              2.*kd1**7*kd2**6*(2809.*omega1**2 + 2390.*omega1*omega2 + 3626.*omega2**2)) - omega1*omega2*(omega1 + omega2)**2*(2.*kd2**15*omega1**2 + 2.*kd1**15*omega2**2 - &
              kd1**14*kd2*omega2*(269.*omega1 + 14.*omega2) - kd1*kd2**14*omega1*(14.*omega1 + 269.*omega2) + kd1**6*kd2**9*(-9858.*omega1**2 + 105135.*omega1*omega2 - &
              32908.*omega2**2) + kd1**7*kd2**8*(-612.*omega1**2 + 460163.*omega1*omega2 - 15128.*omega2**2) + kd1**9*kd2**6*(-32908.*omega1**2 + 105135.*omega1*omega2 - &
              9858.*omega2**2) + kd1**8*kd2**7*(-15128.*omega1**2 + 460163.*omega1*omega2 - 612.*omega2**2) + kd1**2*kd2**13*(386.*omega1**2 - 1461.*omega1*omega2 + &
              120.*omega2**2) + kd1**13*kd2**2*(120.*omega1**2 - 1461.*omega1*omega2 + 386.*omega2**2) + kd1**3*kd2**12*(1354.*omega1**2 - 18101.*omega1*omega2 + 520.*omega2**2) + &
              kd1**12*kd2**3*(520.*omega1**2 - 18101.*omega1*omega2 + 1354.*omega2**2) - kd1**11*kd2**4*(4132.*omega1**2 + 84537.*omega1*omega2 + 2226.*omega2**2) - &
              kd1**4*kd2**11*(2226.*omega1**2 + 84537.*omega1*omega2 + 4132.*omega2**2) - kd1**10*kd2**5*(21720.*omega1**2 + 117173.*omega1*omega2 + 11858.*omega2**2) - &
              kd1**5*kd2**10*(11858.*omega1**2 + 117173.*omega1*omega2 + 21720.*omega2**2))*swd**2) - 2048.*grav*kd1**5*kd2**5*(kd1 + kd2)**5*(384.*grav**2*kd1**2*kd2**2*(kd1 + &
              kd2)**4 + (2.*kd2**6*(90.*omega1**4 + 212.*omega1**3*omega2 + 242.*omega1**2*omega2**2 + 55.*omega1*omega2**3 + 25.*omega2**4) + 2.*kd1**6*(25.*omega1**4 + &
              55.*omega1**3*omega2 + 242.*omega1**2*omega2**2 + 212.*omega1*omega2**3 + 90.*omega2**4) + kd1*kd2**5*(702.*omega1**4 + 1294.*omega1**3*omega2 + 1688.*omega1**2*omega2**2 + &
              331.*omega1*omega2**3 + 200.*omega2**4) + 2.*kd1**3*kd2**3*(593.*omega1**4 + 1367.*omega1**3*omega2 + 2133.*omega1**2*omega2**2 + 1367.*omega1*omega2**3 + &
              593.*omega2**4) + kd1**2*kd2**4*(1150.*omega1**4 + 2395.*omega1**3*omega2 + 3358.*omega1**2*omega2**2 + 1441.*omega1*omega2**3 + 698.*omega2**4) + &
              kd1**5*kd2*(200.*omega1**4 + 331.*omega1**3*omega2 + 1688.*omega1**2*omega2**2 + 1294.*omega1*omega2**3 + 702.*omega2**4) + kd1**4*kd2**2*(698.*omega1**4 + &
              1441.*omega1**3*omega2 + 3358.*omega1**2*omega2**2 + 2395.*omega1*omega2**3 + 1150.*omega2**4))*swd**2) - 131072.*grav*kd1**3*kd2**3*(kd1 + &
              kd2)**3*(64.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(20.*kd1**4 + 32.*kd1**3*kd2 + 9.*kd1**2*kd2**2 + 32.*kd1*kd2**3 + 20.*kd2**4) + (2.*kd2**10*(90.*omega1**4 + &
              210.*omega1**3*omega2 + 218.*omega1**2*omega2**2 + 15.*omega1*omega2**3 + 7.*omega2**4) + kd1*kd2**9*(1098.*omega1**4 + 2161.*omega1**3*omega2 + 2460.*omega1**2*omega2**2 + &
              186.*omega1*omega2**3 + 86.*omega2**4) + 2.*kd1**10*(7.*omega1**4 + 15.*omega1**3*omega2 + 218.*omega1**2*omega2**2 + 210.*omega1*omega2**3 + 90.*omega2**4) + &
              2.*kd1**2*kd2**8*(2355.*omega1**4 + 4819.*omega1**3*omega2 + 5837.*omega1**2*omega2**2 + 1654.*omega1*omega2**3 + 753.*omega2**4) + kd1**9*kd2*(86.*omega1**4 + &
              186.*omega1**3*omega2 + 2460.*omega1**2*omega2**2 + 2161.*omega1*omega2**3 + 1098.*omega2**4) + 2.*kd1**8*kd2**2*(753.*omega1**4 + 1654.*omega1**3*omega2 + &
              5837.*omega1**2*omega2**2 + 4819.*omega1*omega2**3 + 2355.*omega2**4) + 2.*kd1**3*kd2**7*(5646.*omega1**4 + 11504.*omega1**3*omega2 + 15240.*omega1**2*omega2**2 + &
              5538.*omega1*omega2**3 + 2919.*omega2**4) + 2.*kd1**7*kd2**3*(2919.*omega1**4 + 5538.*omega1**3*omega2 + 15240.*omega1**2*omega2**2 + 11504.*omega1*omega2**3 + &
              5646.*omega2**4) + kd1**4*kd2**6*(14642.*omega1**4 + 30155.*omega1**3*omega2 + 49484.*omega1**2*omega2**2 + 19871.*omega1*omega2**3 + 10578.*omega2**4) + &
              kd1**5*kd2**5*(13306.*omega1**4 + 26913.*omega1**3*omega2 + 57172.*omega1**2*omega2**2 + 26913.*omega1*omega2**3 + 13306.*omega2**4) + kd1**6*kd2**4*(10578.*omega1**4 + &
              19871.*omega1**3*omega2 + 49484.*omega1**2*omega2**2 + 30155.*omega1*omega2**3 + 14642.*omega2**4))*swd**2) - 8388608.*grav*kd1*kd2*(kd1 + &
              kd2)*(128.*grav**2*kd1**2*kd2**2*(kd1 + kd2)**4*(3.*kd1**8 + 8.*kd1**7*kd2 + 22.*kd1**6*kd2**2 + 31.*kd1**5*kd2**3 - 7.*kd1**4*kd2**4 + 31.*kd1**3*kd2**5 + 22.*kd1**2*kd2**6 + &
              8.*kd1*kd2**7 + 3.*kd2**8) + (2.*kd1**14*omega2**2*(13.*omega1**2 + 13.*omega1*omega2 + 6.*omega2**2) + 2.*kd2**14*omega1**2*(6.*omega1**2 + 13.*omega1*omega2 + &
              13.*omega2**2) + kd1*kd2**13*omega1*(90.*omega1**3 + 96.*omega1**2*omega2 + 132.*omega1*omega2**2 - 17.*omega2**3) + kd1**13*kd2*omega2*(-17.*omega1**3 + &
              132.*omega1**2*omega2 + 96.*omega1*omega2**2 + 90.*omega2**3) + kd1**8*kd2**6*(20852.*omega1**4 + 31497.*omega1**3*omega2 - 3274.*omega1**2*omega2**2 - 40955.*omega1*omega2**3 - &
              14644.*omega2**4) + kd1**2*kd2**12*(2506.*omega1**4 + 4825.*omega1**3*omega2 + 5326.*omega1**2*omega2**2 + 699.*omega1*omega2**3 + 366.*omega2**4) - &
              kd1**7*kd2**7*(412.*omega1**4 + 19393.*omega1**3*omega2 + 44658.*omega1**2*omega2**2 + 19393.*omega1*omega2**3 + 412.*omega2**4) + 2.*kd1**9*kd2**5*(9357.*omega1**4 + &
              19238.*omega1**3*omega2 + 30168.*omega1**2*omega2**2 + 1561.*omega1*omega2**3 + 1021.*omega2**4) + kd1**3*kd2**11*(11658.*omega1**4 + 23447.*omega1**3*omega2 + &
              29740.*omega1**2*omega2**2 + 4560.*omega1*omega2**3 + 2012.*omega2**4) + kd1**12*kd2**2*(366.*omega1**4 + 699.*omega1**3*omega2 + 5326.*omega1**2*omega2**2 + &
              4825.*omega1*omega2**3 + 2506.*omega2**4) + kd1**4*kd2**10*(17564.*omega1**4 + 36589.*omega1**3*omega2 + 65808.*omega1**2*omega2**2 + 17619.*omega1*omega2**3 + &
              7450.*omega2**4) + 2.*kd1**5*kd2**9*(1021.*omega1**4 + 1561.*omega1**3*omega2 + 30168.*omega1**2*omega2**2 + 19238.*omega1*omega2**3 + 9357.*omega2**4) + &
              kd1**11*kd2**3*(2012.*omega1**4 + 4560.*omega1**3*omega2 + 29740.*omega1**2*omega2**2 + 23447.*omega1*omega2**3 + 11658.*omega2**4) + kd1**10*kd2**4*(7450.*omega1**4 + &
              17619.*omega1**3*omega2 + 65808.*omega1**2*omega2**2 + 36589.*omega1*omega2**3 + 17564.*omega2**4) + kd1**6*kd2**8*(-14644.*omega1**4 - 40955.*omega1**3*omega2 - &
              3274.*omega1**2*omega2**2 + 31497.*omega1*omega2**3 + 20852.*omega2**4))*swd**2) - 536870912.*grav*(64.*grav**2*kd1**3*kd2**3*(kd1 + kd2)**3*(35.*kd1**8 + 110.*kd1**7*kd2 + &
              20.*kd1**6*kd2**2 - 262.*kd1**5*kd2**3 - 415.*kd1**4*kd2**4 - 262.*kd1**3*kd2**5 + 20.*kd1**2*kd2**6 + 110.*kd1*kd2**7 + 35.*kd2**8) + (-3.*kd2**15*omega1**3*omega2 - &
              3.*kd1**15*omega1*omega2**3 + kd1**14*kd2*omega2**2*(462.*omega1**2 + 443.*omega1*omega2 + 226.*omega2**2) + kd1*kd2**14*omega1**2*(226.*omega1**2 + 443.*omega1*omega2 + &
              462.*omega2**2) + kd1**2*kd2**13*omega1*(950.*omega1**3 + 1677.*omega1**2*omega2 + 2112.*omega1*omega2**2 - 334.*omega2**3) + kd1**13*kd2**2*omega2*(-334.*omega1**3 + &
              2112.*omega1**2*omega2 + 1677.*omega1*omega2**2 + 950.*omega2**3) + kd1**3*kd2**12*(4808.*omega1**4 + 8538.*omega1**3*omega2 + 12854.*omega1**2*omega2**2 + &
              3063.*omega1*omega2**3 + 2358.*omega2**4) + 4.*kd1**4*kd2**11*(4850.*omega1**4 + 9842.*omega1**3*omega2 + 13834.*omega1**2*omega2**2 + 6277.*omega1*omega2**3 + &
              2995.*omega2**4) + kd1**12*kd2**3*(2358.*omega1**4 + 3063.*omega1**3*omega2 + 12854.*omega1**2*omega2**2 + 8538.*omega1*omega2**3 + 4808.*omega2**4) + &
              4.*kd1**11*kd2**4*(2995.*omega1**4 + 6277.*omega1**3*omega2 + 13834.*omega1**2*omega2**2 + 9842.*omega1*omega2**3 + 4850.*omega2**4) - kd1**9*kd2**6*(28534.*omega1**4 + &
              32473.*omega1**3*omega2 + 72194.*omega1**2*omega2**2 + 20535.*omega1*omega2**3 + 12406.*omega2**4) + kd1**5*kd2**10*(28824.*omega1**4 + 63903.*omega1**3*omega2 + &
              78244.*omega1**2*omega2**2 + 42553.*omega1*omega2**3 + 13530.*omega2**4) - kd1**6*kd2**9*(12406.*omega1**4 + 20535.*omega1**3*omega2 + 72194.*omega1**2*omega2**2 + &
              32473.*omega1*omega2**3 + 28534.*omega2**4) + kd1**10*kd2**5*(13530.*omega1**4 + 42553.*omega1**3*omega2 + 78244.*omega1**2*omega2**2 + 63903.*omega1*omega2**3 + &
              28824.*omega2**4) - 2.*kd1**8*kd2**7*(44592.*omega1**4 + 84019.*omega1**3*omega2 + 159355.*omega1**2*omega2**2 + 83552.*omega1*omega2**3 + 41329.*omega2**4) - &
              2.*kd1**7*kd2**8*(41329.*omega1**4 + 83552.*omega1**3*omega2 + 159355.*omega1**2*omega2**2 + 84019.*omega1*omega2**3 + 44592.*omega2**4))*swd**2))/ &
              ((-8. + kd1)**3*kd1**2*(8. + kd1)**3*(-8. + kd2)**3*kd2**2*(8. + kd2)**3*omega1*omega2*(omega1 + omega2)*swd**3*(-64.*grav*(kd1 + kd2)**2*(64. + (kd1 + kd2)**2)*(4096. + &
              (kd1 + kd2)**2*(384. + (kd1 + kd2)**2)) + (16777216. + (kd1 + kd2)**2*(7340032. + (kd1 + kd2)**2*(286720. + (kd1 + kd2)**2*(1792. + (kd1 + kd2)**2))))*(omega1 + omega2)**2*swd))
    !
end function velsp44
    !
end subroutine SwashBCtransferfnc
