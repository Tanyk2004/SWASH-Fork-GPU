subroutine SwashDensity
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
!    1.00, May 2013: New subroutine
!
!   Purpose
!
!   Calculates density of water relative to reference density based on temperature, salinity or sediment
!
!   Method
!
!   The Eckart's formula is used here
!
!   C. Eckart
!   Properties of water, part II. The equation of state of water and sea water at low temperatures and pressures
!   Amer. J. of Sci., vol. 256, 225-240, 1958
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata, only: rho, rp
    use SwanGriddata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter over vertical layers
    integer       :: n        ! loop counter over grid elements
    integer       :: ntot     ! number of grid elements
    !
    real*8        :: c        ! local sediment volume concentration
    real*8        :: fac      ! a factor
    real*8        :: lambda   ! polynomial function of temperature and salinity (for sea water) in Eckart's formula
    real*8        :: p0       ! another polynomial function of temperature and salinity (for sea water) in Eckart's formula
    real*8        :: rval1    ! auxiliary real
    real*8        :: rval2    ! auxiliary real
    real*8        :: rval3    ! auxiliary real
    real*8        :: s        ! local salinity
    real*8        :: t        ! local temperature
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashDensity')
    !
    if ( optg /= 5 ) then
       ntot = mcgrd
    else
       ntot = ncells
    endif
    !
    ! initialize local salinity and temperature
    !
    s = dble(salw )
    t = dble(tempw)
    !
    if ( lsal > 0 .and. ltemp > 0 ) then
       !
       ! both space varying salinity and temperature are taken into account
       !
       do k = 1, kmax
          do n = 1, ntot
             !
             s = dble(rp(n,k,lsal ))
             t = dble(rp(n,k,ltemp))
             !
             p0        =  5890d0  + (  38d0  - 375d-3*t)*t +              3d0*s
             lambda    = 17795d-1 + (1125d-2 - 745d-4*t)*t - (38d-1 + 1d-2*t)*s
             !
             rho(n,k) = p0/(lambda+0.698*p0)
             !
          enddo
       enddo
       !
    else if ( lsal > 0 ) then
       !
       ! space varying salinity is taken into account, temperature is constant
       !
       rval1 =  5890d0  + (  38d0  - 375d-3*t)*t
       rval2 = 17795d-1 + (1125d-2 - 745d-4*t)*t
       rval3 =    38d-1 +                 1d-2*t
       !
       do k = 1, kmax
          do n = 1, ntot
             !
             s = dble(rp(n,k,lsal))
             !
             p0        = rval1 +   3d0*s
             lambda    = rval2 - rval3*s
             !
             rho(n,k) = p0/(lambda+0.698*p0)
             !
          enddo
       enddo
       !
    else if ( ltemp > 0 ) then
       !
       ! space varying temperature is taken into account, salinity is constant
       !
       rval1 =  5890d0  +   3d0*s
       rval2 = 17795d-1 - 38d-1*s
       !
       do k = 1, kmax
          do n = 1, ntot
             !
             t = dble(rp(n,k,ltemp))
             !
             p0        = rval1 + (  38d0           - 375d-3*t)*t
             lambda    = rval2 + (1125d-2 - 1d-2*s - 745d-4*t)*t
             !
             rho(n,k) = p0/(lambda+0.698*p0)
             !
          enddo
       enddo
       !
    else
       !
       rho = 0.
       !
    endif
    !
    ! density based on Eckart is g/ml, hence multiply by 1000 to get g/l
    ! also set density relative to reference density
    !
    if ( lsal > 0 .or. ltemp > 0 ) rho = 1000.*rho - rhow
    !
    if ( lsed > 0 .and. lmixt ) then
       !
       ! sediment is taken into account
       !
       fac = dble(rhos-rhow)
       !
       do k = 1, kmax
          do n = 1, ntot
             !
             c = dble(rp(n,k,lsed))
             !
             rho(n,k) = c * fac + rho(n,k) * (1d0 - c)
             !
          enddo
       enddo
       !
    endif
    !
end subroutine SwashDensity
