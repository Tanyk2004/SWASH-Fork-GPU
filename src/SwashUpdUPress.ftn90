subroutine SwashUpdUPress
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
!    1.00, March 2023: New subroutine
!
!   Purpose
!
!   Initializes / updates atmospheric pressure on triangular mesh based on space varying input field
!   Also correct water level on open boundaries, if appropriate
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use m_genarr
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Local variables
!
    integer               :: btype    ! boundary type (see SwashUpdateUData.f90)
    integer               :: icell    ! cell index / loop counter over cells
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer               :: iface    ! face index
    integer               :: jf       ! loop counter over boundary faces
    integer, dimension(3) :: v        ! vertices of present cell
    !
    type(celltype), dimension(:), pointer :: cell ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdUPress')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    if ( ifldyn(13) == 1 ) then
       !
       ! determine pressure in cells
       !
       do icell = 1, ncells
          !
          ! get vertices of present cell
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          patm(icell) = ( presf(v(1),2) + presf(v(2),2) + presf(v(3),2) )/ 3.
          !
       enddo
       !
    endif
    !
    ! correct water level on boundary for local atmospheric pressure, if appropriate
    !
    if ( prmean > 0. ) then
       !
       do jf = 1, nfacesb
          !
          iface = jbface(jf)
          icell = face(iface)%atti(FACEC1)
          !
          btype = face(iface)%atti(FBTYPE)
          !
          if ( btype == 2 ) bcs(iface) = bcs(iface) + ( prmean - patm(icell) ) / (rhow*grav)
          !
       enddo
       !
    endif
    !
end subroutine SwashUpdUPress
