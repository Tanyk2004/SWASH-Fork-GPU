subroutine SwashUpdUFlowFlds
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
!   Initializes / updates flow variables on triangular mesh based on space varying input fields
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use m_genarr
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Local variables
!
    integer               :: icell    ! loop counter over cells
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer               :: iface    ! loop counter of faces
    integer, dimension(3) :: v        ! vertices of present cell/face
    !
    real                  :: nx       ! x-component of normal to face
    real                  :: ny       ! y-component of normal to face
    real                  :: ux       ! x-component of velocity at midface
    real                  :: uy       ! y-component of velocity at midface
    !
    type(celltype), dimension(:), pointer :: cell      ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face      ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdUFlowFlds')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    if ( initsf .or. ifldyn(7) == 1 ) then
       !
       ! determine water level in cells
       !
       do icell = 1, ncells
          !
          ! get vertices of present cell
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          s1(icell) = ( wlevf(v(1),2) + wlevf(v(2),2) + wlevf(v(3),2) )/ 3.
          !
       enddo
       !
       initsf = .false.
       !
    endif
    !
    if ( inituf .or. ifldyn(2) == 1 ) then
       !
       ! determine flow velocity at faces
       !
       do iface = 1, nfaces
          !
          v(1) = face(iface)%atti(FACEV1)
          v(2) = face(iface)%atti(FACEV2)
          !
          nx = face(iface)%attr(FACENORMX)
          ny = face(iface)%attr(FACENORMY)
          !
          ux = 0.5 * ( uxf(v(1),2) + uxf(v(2),2) )
          uy = 0.5 * ( uyf(v(1),2) + uyf(v(2),2) )
          !
          u1(iface,1) = nx * ux + ny * uy
          !
       enddo
       !
       inituf = .false.
       initvf = .false.
       !
    endif
    !
end subroutine SwashUpdUFlowFlds
