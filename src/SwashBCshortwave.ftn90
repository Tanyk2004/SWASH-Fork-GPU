subroutine SwashBCshortwave ( bcfour, nfreq, xp, yp, ibgrpt, swd, wdir, rsgn, vdir, shape )
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
!    1.00, October 2012: New subroutine
!
!   Purpose
!
!   Computes first order free short wave components for synthesizing time series along open boundaries
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
!   Local variables
!
    integer        :: icemds              ! counter for number of evanescent modes
    integer, save  :: ient  = 0           ! number of entries in this subroutine
    integer        :: j                   ! loop counter
    !
    real           :: ampl                ! amplitude of a Fourier component
    real           :: kwav                ! wave number of a Fourier component
    real           :: n                   ! ratio of group and phase velocity
    real           :: omega               ! angular frequency of a Fourier component
    real           :: omegcf              ! cut-off frequency
    real           :: phase               ! phase of a Fourier component
    real           :: rval                ! auxiliary real
    real           :: s                   ! sign
    real           :: theta               ! wave direction of a Fourier component with respect to computational coordinates
    real           :: urmax               ! maximum Ursell number
    !
    logical, save  :: nowarn  = .false.   ! give no warning again
    logical, save  :: nowarn2 = .false.   ! give no other warning again
    !
    character(120) :: msgstr              ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashBCshortwave')
    !
    ! determine the cut-off frequency above which are the evanescent modes to be filtered out
    ! (based on numerical dispersion relation; to be saved a slightly smaller cut-off is chosen)
    !
    omegcf = 0.9 * 2. * real(kmax) * sqrt(grav/swd)
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
    urmax  = -999.
    icemds = 0
    !
    do j = 1, nfreq
       !
       ampl  = bcfour%ampl (j)
       omega = bcfour%omega(j)
       phase = bcfour%phase(j)
       !
       theta = wdir - alpc + s*bcfour%theta(j)
       !
       ! filter out the evanescent modes, if appropriate
       !
       if ( numdisp ) then
          !
          if ( omega > omegcf ) then
             icemds         = icemds + 1
             ampl           = 0.
             bcfour%ampl(j) = 0.
          endif
          !
       endif
       !
       ! calculate wave number and store it
       !
       call disprel ( swd, omega, kwav, rval, n )
       !
       kwave(ibgrpt,j) = kwav
       !
       ! correct amplitude in case of TMA spectrum for shallow water
       !
       if ( shape == 3 ) ampl = ampl * omega * omega / ( grav * kwav * sqrt(2.*n) )
       !
       ! compute the maximum Ursell number
       !
       if ( ampl /= 0. .and. kwav /= 0. ) then
          !
          ! this is the classical definition of the Ursell number
          !rval = ampl / swd / swd / swd / kwav / kwav
          ! this is another variant of the Ursell number more appropriate for large kd > 1 or intermediate depths (see Beji, 1995)
          rval = ampl * kwav / ( tanh(kwav*swd)*tanh(kwav*swd)*tanh(kwav*swd) )
          if ( rval > urmax ) urmax = rval
          !
       endif
       !
       ! in case of periodicity, wave direction must be corrected so that wave number is an integer multiple of 2pi/length with length the periodicity length
       !
       if ( bcperx ) then
          !
          rval = nint( kwav*cos(theta) / ( pi2/xclen ) ) * pi2/xclen / kwav
          if ( rval > 1. ) then
             theta = acos ( rval - pi2/xclen / kwav )
          else if ( rval < -1. ) then
             theta = acos ( rval + pi2/xclen / kwav )
          else
             theta = acos ( rval )
          endif
          if ( rsgn == -1. ) theta = pi2 - theta
          !
       else if ( bcpery ) then
          !
          rval = nint( kwav*sin(theta) / ( pi2/yclen ) ) * pi2/yclen / kwav
          if ( rval > 1. ) then
             theta = asin ( rval - pi2/yclen / kwav )
          else if ( rval < -1. ) then
             theta = asin ( rval + pi2/yclen / kwav )
          else
             theta = asin ( rval )
          endif
          if ( rsgn == -1. ) theta = pi - theta
          !
       endif
       !
       ! check this direction with respect to the normal of boundary
       ! (must be within -80 degrees to 80 degrees)
       !
       if ( vdir ) then
          if ( rsgn == 1. ) then
             if ( cos(theta) <  0.174 ) cycle
          else
             if ( cos(theta) > -0.174 ) cycle
          endif
       else
          if ( rsgn == 1. ) then
             if ( sin(theta) <  0.174 ) cycle
          else
             if ( sin(theta) > -0.174 ) cycle
          endif
       endif
       !
       ! include phase shift related to wave direction and wave number
       !
       phase = phase + kwav * ( cos(theta+alpc)*xp + sin(theta+alpc)*yp )
       !
       ! store free wave components
       !
       comp1(ibgrpt,j) = ampl * cos( phase )
       comp2(ibgrpt,j) = ampl * sin( phase )
       !
    enddo
    !
    ! give warning if maximum Ursell number > 0.2
    !
    if ( urmax > 0.2 .and. .not.nowarn2 ) then
       !
       write (msgstr,'(a,f5.2,a)') 'the Ursell number associated with the wavemaker-generated wave field = ',urmax, ' > 0.2'
       call msgerr (1, trim(msgstr) )
       write (PRINTF,'(a)') '                       (linear wave theory and possible second order bound waves not valid)'
       !
       nowarn2 = .true.
       !
    endif
    !
    ! give warning for filtering out the evanescent modes
    !
    if ( icemds > 0 .and. .not.nowarn ) then
       !
       rval = 100.*real(icemds)/real(nfreq)
       !
       if ( .not. rval < 10. ) then
          write (msgstr,'(a,f5.1,a,f5.2,a)') 'percentage of wave components on boundary that have been filtered out (evanescent modes) = ',rval, ' (cut-off = ',omegcf/pi2,' Hz)'
          call msgerr (1, trim(msgstr) )
       endif
       !
       nowarn = .true.
       !
    endif
    !
end subroutine SwashBCshortwave
