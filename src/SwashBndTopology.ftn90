subroutine SwashBndTopology
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
!    1.00, July 2022: New subroutine
!
!   Purpose
!
!   Setups the boundary topology
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Local variables
!
    integer                               :: icell    ! loop counter over cells
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: iface    ! loop counter over faces
    integer                               :: jb       ! loop counter over boundary faces/cells
    !
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face     ! datastructure for faces with their attributes
    !
    type bcpt                                         ! linked list for boundary cells
       integer             :: jb
       type(bcpt), pointer :: nextbc
    end type bcpt
    type(bcpt), target     :: frstbc
    type(bcpt), pointer    :: currbc, tmpbc
    !
    type bfpt                                         ! linked list for boundary faces
       integer             :: jb
       type(bfpt), pointer :: nextbf
    end type bfpt
    type(bfpt), target     :: frstbf
    type(bfpt), pointer    :: currbf, tmpbf
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashBndTopology')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    ! determine number of boundary cells and their indices
    !
    ncellsb   = 0
    frstbc%jb = 0
    nullify(frstbc%nextbc)
    currbc => frstbc
    do icell = 1, ncells
       if ( cell(icell)%atti(CMARKER) /= 0 ) then   ! boundary cell
          ncellsb = ncellsb + 1
          allocate(tmpbc)
          tmpbc%jb = icell
          nullify(tmpbc%nextbc)
          currbc%nextbc => tmpbc
          currbc => tmpbc
       endif
    enddo
    !
    if (.not.allocated(jbcell)) allocate(jbcell(ncellsb))
    !
    if ( ncellsb > 0 ) then
       currbc => frstbc%nextbc
       do jb = 1, ncellsb
          jbcell(jb) = currbc%jb
          currbc => currbc%nextbc
       enddo
       deallocate(tmpbc)
    endif
    !
    ! determine number of boundary faces and their indices
    !
    nfacesb   = 0
    frstbf%jb = 0
    nullify(frstbf%nextbf)
    currbf => frstbf
    do iface = 1, nfaces
       if ( face(iface)%atti(FMARKER) /= 0 ) then   ! boundary face
          nfacesb = nfacesb + 1
          allocate(tmpbf)
          tmpbf%jb = iface
          nullify(tmpbf%nextbf)
          currbf%nextbf => tmpbf
          currbf => tmpbf
       endif
    enddo
    !
    if (.not.allocated(jbface)) allocate(jbface(nfacesb))
    !
    if ( nfacesb > 0 ) then
       currbf => frstbf%nextbf
       do jb = 1, nfacesb
          jbface(jb) = currbf%jb
          currbf => currbf%nextbf
       enddo
       deallocate(tmpbf)
    endif
    !
end subroutine SwashBndTopology
