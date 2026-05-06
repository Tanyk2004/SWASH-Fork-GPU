subroutine SwashIntWavgen ( igser, nfreq, wdir, shape )
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
!    1.00, June 2019: New subroutine
!
!   Purpose
!
!   Computes the source function amplitude and shape factor of source area
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
    integer, intent(in)         :: nfreq  ! number of frequencies
    integer, intent(in)         :: shape  ! spectral shape
                                          ! = 1; Pierson Moskowitz
                                          ! = 2; Jonswap
                                          ! = 3; TMA
    !
    real, intent(in)            :: wdir   ! incident or peak wave direction with respect to problem coordinates
    !
    !
    type(bfsdat), intent(inout) :: igser  ! list containing parameters for internal-generated series of wave components
!
!   Local variables
!
    integer        :: icemds              ! counter for number of evanescent modes
    integer, save  :: ient  = 0           ! number of entries in this subroutine
    integer        :: j                   ! loop counter
    !
    real           :: ampl                ! amplitude of an internal-generated wave component
    real           :: beta                ! shape factor beta of the source area
    real           :: cdn                 ! part of formula for computing energy velocity (denominator)
    real           :: cen                 ! energy velocity
    real           :: cnu                 ! part of formula for computing energy velocity (numerator)
    real           :: fac1                ! auxiliary factor
    real           :: fac2                ! auxiliary factor
    real           :: ishap               ! shape factor I of the source area
    real           :: kwav                ! wave number of an internal-generated wave component
    real           :: n                   ! ratio of group and phase velocity
    real           :: omega               ! angular frequency of an internal-generated wave component
    real           :: omegcf              ! cut-off frequency
    real           :: rval                ! auxiliary real
    real           :: swd                 ! still water depth
    real           :: theta               ! wave direction of an internal-generated component with respect to computational coordinates
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
    if (ltrace) call strace (ient,'SwashIntWavgen')
    !
    ! get still water depth
    !
    swd = piwg(3)
    !
    ! determine the cut-off frequency above which are the evanescent modes to be filtered out
    ! (based on numerical dispersion relation; to be saved a slightly smaller cut-off is chosen)
    !
    omegcf = 0.9 * 2. * real(kmax) * sqrt(grav/swd)
    !
    allocate(igser%sfamp(nfreq))
    allocate(igser%bshap(nfreq))
    !
    igser%sfamp = 0.
    igser%bshap = 0.
    !
    urmax  = -999.
    icemds = 0
    !
    do j = 1, nfreq
       !
       ampl  = igser%ampl (j)
       omega = igser%omega(j)
       !
       theta = wdir - alpc + igser%theta(j)
       !
       ! filter out the evanescent modes, if appropriate
       !
       if ( numdisp ) then
          !
          if ( omega > omegcf ) then
             icemds         = icemds + 1
             ampl           = 0.
             igser%ampl(j)  = 0.
          endif
          !
       endif
       !
       ! calculate wave number and store it
       !
       call disprel ( swd, omega, kwav, rval, n )
       !
       kwave(j,1) = kwav
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
          !
       endif
       !
       ! calculate shape factor beta
       !
       beta = 80. / piwg(4) ** 2 / ( pi2/kwav ) ** 2
       !
       ! calculate energy velocity
       !
       fac1 = kwav * swd
       fac2 = sqrt( grav * swd )
       !
       if ( kpmax == 1 ) then
          !
          cnu = 8. * fac2
          cdn = (4. + fac1 ** 2.) ** 1.5
          !
       else if ( kpmax == 2 ) then
          !
          cnu = 64. * fac2 * (256. + 32. * fac1 ** 2. + 5. * fac1 ** 4.)
          cdn = sqrt((16. + fac1 ** 2.) * (256. + 96. * fac1 ** 2. + fac1 ** 4.) ** 3.)
          !
       else if ( kpmax == 3 ) then
          !
          cnu = 72. * fac2 * (5038848. + 933120. * fac1 ** 2. + 147744. * fac1 ** 4. + 3024. * fac1 ** 6. + 35. * fac1 ** 8.)
          cdn = sqrt((1296. + 120. * fac1 ** 2. + fac1 ** 4.) * (46656. + 19440. * fac1 ** 2. + 540. * fac1 ** 4. + fac1 ** 6.) ** 3.)
          !
       else
          !
          cnu = (0.5 + fac1 / sinh( 2. * fac1 )) * (omega / kwav)
          cdn = 1.
          !
       endif
       !
       cen = cnu / cdn
       !
       ! calculate shape factor I
       !
       if ( optg /= 5 ) then
          !
          if ( lsrcfy ) then
             !
             ishap = sqrt( pi/beta ) * exp( -( kwav * cos(theta) ) ** 2 / 4. / beta ) / cos(theta)
             !
          else if ( lsrcfx ) then
             !
             ishap = sqrt( pi/beta ) * exp( -( kwav * sin(theta) ) ** 2 / 4. / beta ) / sin(theta)
             !
          endif
          !
       else
          !
          ishap = sqrt( pi/beta ) * exp( -( kwav ) ** 2 / 4. / beta )
          !
       endif
       !
       ! store parameters for internal wave generation
       !
       igser%theta(j) = theta
       igser%bshap(j) = beta
       igser%sfamp(j) = 2. * ampl * cen / ishap
       !
    enddo
    !
    ! give warning if maximum Ursell number > 0.2
    !
    if ( urmax > 0.2 .and. .not.nowarn2 ) then
       !
       write (msgstr,'(a,f5.2,a)') 'the Ursell number associated with the internal-generated wave field = ',urmax, ' > 0.2'
       call msgerr (1, trim(msgstr) )
       write (PRINTF,'(a)') '                       (linear wave theory is not valid)'
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
          write (msgstr,'(a,f5.1,a,f5.2,a)') 'percentage of internal-generated components that have been filtered out (evanescent modes) = ',rval, ' (cut-off = ',omegcf/pi2,' Hz)'
          call msgerr (1, trim(msgstr) )
       endif
       !
       nowarn = .true.
       !
    endif
    !
end subroutine SwashIntWavgen
