subroutine SwashUBreakPoint
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
!    1.00,   March 2020: New subroutine
!
!   Purpose
!
!   Determines cells of unstructured mesh where steep bore-like wave front occurs for wave breaking
!
!   Method
!
!   Update mask array at those cells where the vertical speed of the free surface
!   exceeds a fraction of the wave phase speed. At those cells, hydrostatic pressure
!   is assumed and remains so at the front face of the breaking wave.
!
!   This approach combined with a proper momentum conservation leads to a correct
!   amount of energy dissipation on the front face of the breaking wave. Moreover,
!   nonlinear wave properties such as asymmetry and skewness are preserved as well.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm, only: dt
    use m_genarr, only: iwrk
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Local variables
!
    integer       :: icell    ! loop counter over cells
    integer       :: icell1   ! sequence number of cell 1 adjacent to present face
    integer       :: icell2   ! sequence number of cell 2 adjacent to present face
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: iface    ! face index
    integer       :: jf       ! loop counter
    !
    real          :: alpha    ! threshold parameter corresponding to onset of wave breaking
    real          :: beta     ! threshold parameter corresponding to re-initiation of breaking in post-breaking area
    real          :: dsdt     ! vertical velocity of the free surface
    real          :: rootgh   ! shallow water celerity
    !
    type(celltype), dimension(:), pointer  :: cell ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer  :: face ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUBreakPoint')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    alpha = psurf(1)
    beta  = psurf(2)
    !
    ! store array brks temporarily
    !
    iwrk(:,1) = brks(:)
    !
    do icell = 1, ncells
       !
       if ( wets(icell) == 1 ) then
          !
          rootgh = sqrt( grav * hs(icell) )
          !
          dsdt = ( s1(icell) - s0(icell) ) / dt
          !
          if ( dsdt > alpha * rootgh ) then
             !
             brks(icell) = 1
             !
             q(icell,:) = 0.
             !
          else if ( dsdt > beta * rootgh ) then
             !
             ! loop over faces of the cell
             !
             do jf = 1, cell(icell)%nof
                !
                ! face identifier
                !
                iface = cell(icell)%face(jf)%atti(FACEID)
                !
                ! consider adjacent cells of current face
                !
                icell1 = face(iface)%atti(FACEC1)
                icell2 = face(iface)%atti(FACEC2)
                !
                if ( icell2 /= 0 ) then
                   !
                   if ( (icell == icell1 .and. iwrk(icell2,1) == 1) .or. (icell == icell2 .and. iwrk(icell1,1) == 1) ) then
                      !
                      brks(icell) = 1
                      !
                      q(icell,:) = 0.
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          else if ( dsdt > 0. .and. brks(icell) == 1 ) then
             !
             q(icell,:) = 0.
             !
          else
             !
             brks(icell) = 0
             !
          endif
          !
       else
          !
          brks(icell) = 0
          !
       endif
       !
    enddo
    !
    ! re-update mask array for wetting and drying at cells by taking into account the breaking points
    !
    do icell = 1, ncells
       !
       if ( wets(icell) == 1 .and. brks(icell) == 1 ) wets(icell) = 0
       !
    enddo
    !
end subroutine SwashUBreakPoint
