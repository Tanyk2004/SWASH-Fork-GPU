subroutine SwashLayUIntfaces
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
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
!    1.00, April 2020: New subroutine
!
!   Purpose
!
!   Determines position of the layer interfaces and layer thicknesses at waterlevel and velocity points on unstructured mesh
!
!   Method
!
!   In the vertical direction the computational domain is divided into a fixed number
!   of layers of which the k-index has range (0,1,...,kmax) and counting downward from
!   the free surface.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwanGriddata, only: ncells, nfaces
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter over velocity layers
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashLayUIntfaces')
    !
    ! determine layer interfaces at wl-point
    !
    zks(:,0) = -dps(:) + hs(:)
    !
    call sigmacoor( zks, hs, ncells )
    !
    ! determine layer thicknesses at wl-point
    !
    do k = 1, kmax
       !
       hks(:,k) = zks(:,k-1) - zks(:,k)
       !
    enddo
    !
    ! first, calculate layer thicknesses based on humn (extrapolated in time; see subroutine SwashUpdateUDepths)
    ! note: use zkum temporarily!
    !
    ! --- determine layer interfaces at u-point
    !
    zkum(:,0) = zkum(:,kmax) + humn(:)
    !
    call sigmacoor( zkum, humn, nfaces )
    !
    ! --- determine layer thicknesses at u-point
    !
    do k = 1, kmax
       !
       hkumn(:,k) = zkum(:,k-1) - zkum(:,k)
       !
    enddo
    !
    ! next, calculate the layer thicknesses based on hu and hum
    !
    ! --- determine layer interfaces at u-point
    !
    zku (:,0) = zku (:,kmax) + hu (:)
    zkum(:,0) = zkum(:,kmax) + hum(:)
    !
    call sigmacoor( zku , hu , nfaces )
    call sigmacoor( zkum, hum, nfaces )
    !
    ! --- determine layer thicknesses at u-point
    !
    do k = 1, kmax
       !
       hku (:,k) = zku (:,k-1) - zku (:,k)
       hkum(:,k) = zkum(:,k-1) - zkum(:,k)
       !
    enddo
    !
end subroutine SwashLayUIntfaces
