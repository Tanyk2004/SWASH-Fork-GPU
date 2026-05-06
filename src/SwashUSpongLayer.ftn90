subroutine SwashUSpongLayer
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
!    1.00, May 2020: New subroutine
!
!   Purpose
!
!   Determines damping function for flow variables due to sponge layers on unstructured mesh
!
!   Method
!
!   The sponge layer formula as described in Mayer et al. (1998), Eqs (43) and (44), is employed
!
!   S. Mayer, A. Garapon and L.S. Sorensen
!   A fractional step method for unsteady free-surface flow with applications to nonlinear wave dynamics
!   IJNMF, vol. 28, 293-315, 1998
!
!   Modules used
!
    use ocpcomm4
    use m_genarr, only: vmspon, wdspon
    use SwashCommdata3
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Parameter variables
!
    real, parameter :: grt = 0.5  ! growth rate of relaxation
!
!   Local variables
!
    integer               :: bcell      ! boundary cell
    integer               :: bface      ! boundary face
    integer               :: i          ! loop counter
    integer               :: icell      ! loop counter over cells
    integer, save         :: ient = 0   ! number of entries in this subroutine
    integer               :: iface      ! loop counter over faces
    integer               :: j          ! loop counter
    integer               :: jb         ! loop counter over boundary faces/cells
    integer               :: k          ! loop counter
    integer, dimension(3) :: v          ! vertices in present cell
    integer, dimension(2) :: vb         ! vertices of boundary face
    integer               :: vm         ! boundary marker
    !
    real                  :: dist       ! distance to boundary
    real                  :: dmin       ! shortest distance to boundary
    real                  :: dnb        ! normalised distance to boundary
    real                  :: width      ! width of sponge layer
    real                  :: x          ! x-coordinate of actual element
    real                  :: xb         ! x-coordinate of boundary element
    real                  :: y          ! y-coordinate of actual element
    real                  :: yb         ! y-coordinate of boundary element
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
    if (ltrace) call strace (ient,'SwashUSpongLayer')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    do i = 1, nspl
       !
       vm    = vmspon(i)
       width = wdspon(i)
       !
       do iface = 1, nfaces
          !
          x = face(iface)%attr(FACEMX)
          y = face(iface)%attr(FACEMY)
          !
          dmin = 100.*max(xclen,yclen)
          !
          do jb = 1, nfacesb
             !
             j = jbface(jb)
             !
             ! get vertices of present face
             !
             vb(1) = face(j)%atti(FACEV1)
             vb(2) = face(j)%atti(FACEV2)
             !
             ! check boundary marker
             !
             if ( vmark(vb(1)) == vm .and. vmark(vb(2)) == vm ) then
                !
                xb = face(j)%attr(FACEMX)
                yb = face(j)%attr(FACEMY)
                !
                dist = sqrt ( ( x - xb )**2 + ( y - yb )**2 )
                !
                if ( dist < dmin ) then
                   dmin  = dist
                   bface = j
                endif
                !
             endif
             !
          enddo
          !
          if ( .not. dmin > width ) then
             !
             dnb = 1. - dmin / width
             !
             if ( dnb > 0. ) then
                !
                sponu(iface,i)%gamma = grt * dnb**3 + (1.-grt) * dnb**6
                sponu(iface,i)%bface = bface
                !
             else
                !
                sponu(iface,i)%gamma = 0.
                sponu(iface,i)%bface = 1
                !
             endif
             !
          endif
          !
       enddo
       !
       do icell = 1, ncells
          !
          x = cell(icell)%attr(CELLCX)
          y = cell(icell)%attr(CELLCY)
          !
          dmin = 100.*max(xclen,yclen)
          !
          do jb = 1, ncellsb
             !
             j = jbcell(jb)
             !
             ! get vertices of present cell
             !
             v(1) = cell(j)%atti(CELLV1)
             v(2) = cell(j)%atti(CELLV2)
             v(3) = cell(j)%atti(CELLV3)
             !
             ! get two boundary vertices
             !
             do k = 1, 3
                if ( vmark(v(k)) == 0 ) then
                   vb(1) = v(mod(k  ,3)+1)
                   vb(2) = v(mod(k+1,3)+1)
                   exit
                endif
             enddo
             !
             ! check boundary marker
             !
             if ( vmark(vb(1)) == vm .and. vmark(vb(2)) == vm ) then
                !
                xb = cell(j)%attr(CELLCX)
                yb = cell(j)%attr(CELLCY)
                !
                dist = sqrt ( ( x - xb )**2 + ( y - yb )**2 )
                !
                if ( dist < dmin ) then
                   dmin  = dist
                   bcell = j
                endif
                !
             endif
             !
          enddo
          !
          if ( .not. dmin > width ) then
             !
             dnb = 1. - dmin / width
             !
             if ( dnb > 0. ) then
                !
                spons(icell,i)%gamma = grt * dnb**3 + (1.-grt) * dnb**6
                spons(icell,i)%bcell = bcell
                !
             else
                !
                spons(icell,i)%gamma = 0.
                spons(icell,i)%bcell = 1
                !
             endif
             !
          endif
          !
       enddo
       !
    enddo
    !
end subroutine SwashUSpongLayer
