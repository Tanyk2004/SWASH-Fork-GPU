subroutine SwashLayerIntfaces
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
!    4.01: Dirk Rijnsdorp, Pieter Smit
!
!   Updates
!
!    1.00,      July 2010: New subroutine
!    4.01, September 2016: subgrid approach
!
!   Purpose
!
!   Determines position of the layer interfaces and layer thicknesses at waterlevel and velocity points
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
    use m_genarr, only: npu
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! loop counter over velocity layers in a pressure layer
    integer       :: k        ! loop counter over velocity layers
    integer       :: kp       ! loop counter over pressure layers
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashLayerIntfaces')
    !
    ! determine layer interfaces at wl-point
    !
    zks(:,0) = -dps(:) + hs(:)
    !
    call sigmacoor( zks, hs, mcgrd )
    !
    ! determine layer thicknesses at wl-point
    !
    do k = 1, kmax
       !
       hks(:,k) = zks(:,k-1) - zks(:,k)
       !
    enddo
    !
    ! first, calculate layer thicknesses based on humn (extrapolated in time; see subroutine SwashUpdateDepths)
    ! note: use zkum temporarily!
    !
    ! --- determine layer interfaces at u-point
    !
    zkum(:,0) = zkum(:,kmax) + humn(:)
    !
    call sigmacoor( zkum, humn, mcgrd )
    !
    ! --- determine layer thicknesses at u-point
    !
    do k = 1, kmax
       !
       hkumn(:,k) = zkum(:,k-1) - zkum(:,k)
       !
    enddo
    !
    ! next, calculate the usual layer thicknesses based on hu and hum
    !
    ! --- determine layer interfaces at u-point
    !
    zku (:,0) = zku (:,kmax) + hu (:)
    zkum(:,0) = zkum(:,kmax) + hum(:)
    !
    call sigmacoor( zku , hu , mcgrd )
    call sigmacoor( zkum, hum, mcgrd )
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
    ! determine layer interfaces and layer thicknesses at v-point, if appropriate
    !
    if ( .not.oned ) then
       !
       zkvm(:,0) = zkvm(:,kmax) + hvmn(:)
       !
       call sigmacoor( zkvm, hvmn, mcgrd )
       !
       do k = 1, kmax
          !
          hkvmn(:,k) = zkvm(:,k-1) - zkvm(:,k)
          !
       enddo
       !
       zkv (:,0) = zkv (:,kmax) + hv (:)
       zkvm(:,0) = zkvm(:,kmax) + hvm(:)
       !
       call sigmacoor( zkv , hv , mcgrd )
       call sigmacoor( zkvm, hvm, mcgrd )
       !
       do k = 1, kmax
          !
          hkv (:,k) = zkv (:,k-1) - zkv (:,k)
          hkvm(:,k) = zkvm(:,k-1) - zkvm(:,k)
          !
       enddo
       !
    endif
    !
    ! determine layer thicknesses of pressure layers, if appropriate
    !
    if ( lsubg ) then
       !
       hksc  = 0.
       hkuc  = 0.
       hkumc = 0.
       !
       if ( .not.oned ) then
          !
          hkvc  = 0.
          hkvmc = 0.
          !
       endif
       !
       k = 0
       !
       do kp = 1, kpmax
          !
          do j = 1, npu(kp)
             !
             k = k + 1
             !
             hksc (:,kp) = hksc (:,kp) + hks (:,k)
             !
             hkuc (:,kp) = hkuc (:,kp) + hku (:,k)
             hkumc(:,kp) = hkumc(:,kp) + hkum(:,k)
             !
             if ( .not.oned ) then
                !
                hkvc (:,kp) = hkvc (:,kp) + hkv (:,k)
                hkvmc(:,kp) = hkvmc(:,kp) + hkvm(:,k)
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
end subroutine SwashLayerIntfaces
