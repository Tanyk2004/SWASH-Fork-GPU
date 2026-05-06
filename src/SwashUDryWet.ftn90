subroutine SwashUDryWet
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
!    1.00, February 2020: New subroutine
!
!   Purpose
!
!   Performs wetting and drying checks in cells and faces of triangular mesh
!
!   Method
!
!   Update mask arrays
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: icell    ! loop counter over cells
    integer       :: iface    ! face index / loop counter over faces
    integer       :: ival     ! integer value
    integer       :: jf       ! loop counter
    !
    type(celltype), dimension(:), pointer  :: cell     ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer  :: face     ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUDryWet')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    ! loop over faces
    !
    do iface = 1, nfaces
       !
       if ( hu(iface) > epshu ) then
          !
          wetu(iface) = 1
          !
       else
          !
          wetu(iface) = 0
          !
       endif
       !
    enddo
    !
    ! loop over cells
    !
    do icell = 1, ncells
       !
       ival = 1
       !
       do jf = 1, cell(icell)%nof
          !
          ! face identifier
          !
          iface = cell(icell)%face(jf)%atti(FACEID)
          !
          ival = ival * ( 1 - wetu(iface) )
          !
       enddo
       !
       if ( ival == 0 ) then
          !
          wets(icell) = 1
          !
       else
          !
          wets(icell) = 0
          !
       endif
       !
    enddo
    !
end subroutine SwashUDryWet
