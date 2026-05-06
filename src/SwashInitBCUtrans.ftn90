subroutine SwashInitBCUtrans
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
!    1.00, January 2023: New subroutine
!
!   Purpose
!
!   Initializes transport constituents based on input fields and stores boundary values as Dirichlet type condition in case of unstructured grid
!
!   Modules used
!
    use ocpcomm4
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
    integer                               :: icell    ! loop counter over cells
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: iface    ! face index
    integer                               :: jf       ! loop counter
    integer, dimension(3)                 :: v        ! vertices of present cell
    integer                               :: vf1      ! first vertex of present face
    integer                               :: vf2      ! second vertex of present face
    !
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face     ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashInitBCUtrans')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    if ( lsal > 0 ) then
       !
       ! specify salinity in centroids based on input field
       !
       do icell = 1, ncells
          !
          ! get vertices of present cell
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          rp(icell,:,lsal) = ( salf(v(1),:) + salf(v(2),:) + salf(v(3),:) ) / 3.
          !
       enddo
       !
       ! store boundary values at real boundaries
       !
       do jf = 1, nfacesb
          !
          iface = jbface(jf)
          !
          ! get vertices of current face
          !
          vf1 = face(iface)%atti(FACEV1)
          vf2 = face(iface)%atti(FACEV2)
          !
          cbndu(jf,:,lsal) = 0.5 * ( salf(vf1,:) + salf(vf2,:) )
          bcrp (jf,:,lsal) = cbndu(jf,:,lsal)
          !
       enddo
       !
    endif
    !
    if ( ltemp > 0 ) then
       !
       ! specify temperature in centroids based on input field
       !
       do icell = 1, ncells
          !
          ! get vertices of present cell
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          rp(icell,:,ltemp) = ( tempf(v(1),:) + tempf(v(2),:) + tempf(v(3),:) ) / 3.
          !
       enddo
       !
       ! store boundary values at real boundaries
       !
       do jf = 1, nfacesb
          !
          iface = jbface(jf)
          !
          ! get vertices of current face
          !
          vf1 = face(iface)%atti(FACEV1)
          vf2 = face(iface)%atti(FACEV2)
          !
          cbndu(jf,:,ltemp) = 0.5 * ( tempf(vf1,:) + tempf(vf2,:) )
          bcrp (jf,:,ltemp) = cbndu(jf,:,ltemp)
          !
       enddo
       !
    endif
    !
    if ( lsed > 0 ) then
       !
       ! specify volumetric sediment in centroids based on input field
       !
       do icell = 1, ncells
          !
          ! get vertices of present cell
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          rp(icell,:,lsed) = ( sedf(v(1),:) + sedf(v(2),:) + sedf(v(3),:) ) / 3. / rhos
          !
       enddo
       !
       ! store boundary values at real boundaries
       !
       do jf = 1, nfacesb
          !
          iface = jbface(jf)
          !
          ! get vertices of current face
          !
          vf1 = face(iface)%atti(FACEV1)
          vf2 = face(iface)%atti(FACEV2)
          !
          cbndu(jf,:,lsed) = 0.5 * ( sedf(vf1,:) + sedf(vf2,:) ) / rhos
          bcrp (jf,:,lsed) = cbndu(jf,:,lsed)
          !
       enddo
       !
    endif
    !
end subroutine SwashInitBCUtrans
