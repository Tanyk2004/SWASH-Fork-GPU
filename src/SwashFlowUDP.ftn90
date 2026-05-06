subroutine SwashFlowUDP
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
!    1.00, January 2017: New subroutine
!
!   Purpose
!
!   Determines bottom values in centroids and faces
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Local variables
!
    integer                               :: icell     ! loop counter over cells
    integer                               :: icell1    ! sequence number of cell 1 adjacent to present face
    integer                               :: icell2    ! sequence number of cell 2 adjacent to present face
    integer                               :: iface     ! face index / loop counter over faces
    integer, save                         :: ient = 0  ! number of entries in this subroutine
    integer                               :: jf        ! loop counter
    integer, dimension(3)                 :: v         ! vertices of present cell
    !
    real                                  :: dep       ! local depth
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
    if (ltrace) call strace (ient,'SwashFlowUDP')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    if ( dpsopt == 1 ) then
       !
       ! first determine bottom values at faces by means of used-defined input and the minimum option
       !
       do iface = 1, nfaces
          !
          ! get vertices of present face
          !
          v(1) = face(iface)%atti(FACEV1)
          v(2) = face(iface)%atti(FACEV2)
          !
          dpu(iface) = min( depf(v(1)), depf(v(2)) )
          !
       enddo
       !
       ! next determine bottom values in centroids using the tiled approach
       !
       do icell = 1, ncells
          !
          ! loop over faces of the cell
          !
          dep = 99999999.
          !
          do jf = 1, cell(icell)%nof
             !
             ! face identifier
             !
             iface = cell(icell)%face(jf)%atti(FACEID)
             !
             dep = min( dep, dpu(iface) )
             !
          enddo
          !
          dps(icell) = dep
          !
          if ( s1(icell) < epsdry - dps(icell) ) then
             !
             s1(icell) = 0.99*epsdry - dps(icell)
             !
          endif
          !
       enddo
       !
    else if ( dpsopt == 2 ) then
       !
       ! first determine bottom values in centroids by means of user-defined input and the mean option
       ! note: bottom level is constant within cell
       !
       do icell = 1, ncells
          !
          ! get vertices of present cell
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          dps(icell) = ( depf(v(1)) + depf(v(2)) + depf(v(3)) )/ 3.
          !
          if ( s1(icell) < epsdry - dps(icell) ) then
             !
             s1(icell) = 0.99*epsdry - dps(icell)
             !
          endif
          !
       enddo
       !
       ! next determine bottom values at faces using the tiled approach
       !
       do iface = 1, nfaces
          !
          icell1 = face(iface)%atti(FACEC1)
          icell2 = face(iface)%atti(FACEC2)
          !
          if ( icell2 /= 0 ) then
             !
             dpu(iface) = min( dps(icell1), dps(icell2) )
             !
          else
             !
             dpu(iface) = dps(icell1)
             !
          endif
          !
       enddo
       !
    else
       !
       call msgerr(2, 'option BOTCEL = MAX and SHIFT not supported')
       !
    endif
    !
end subroutine SwashFlowUDP
