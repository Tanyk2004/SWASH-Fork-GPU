subroutine SwanGridCell ( ncells, nverts, xcugrd, ycugrd, kvertc )
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
!   40.80: Marcel Zijlema
!
!   Updates
!
!   40.80, July 2007: New subroutine
!
!   Purpose
!
!   Fills cell-based data structure
!
!   Method
!
!   Based on unstructured grid
!   Note: we restrict ourselves to triangles only!
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                       :: ncells  ! number of cells in grid
    integer, intent(in)                       :: nverts  ! number of vertices in grid
    !
    integer, dimension(3, ncells), intent(in) :: kvertc  ! vertices of the cell
                                                         ! (must be filled by a gridgenerator!)
    !
    real, dimension(nverts), intent(in)       :: xcugrd  ! the x-coordinates of the grid vertices
    real, dimension(nverts), intent(in)       :: ycugrd  ! the y-coordinates of the grid vertices
!
!   Local variables
!
    integer                               :: icell    ! loop counter over cells / index of present cell
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: ivert    ! loop counter over vertices / index of present vertex
    integer                               :: j        ! loop counter
    integer                               :: jcell    ! index of next cell
    integer                               :: k        ! auxiliary integer / loop counter
    integer                               :: l        ! loop counter
    integer                               :: ncc      ! number of cells with circumcenter outside the triangle
    integer                               :: numcrs   ! number of crossings with faces of cell
    integer                               :: noc      ! number of cells around considered vertex
    integer, parameter                    :: nov = 3  ! number of vertices in present cell (triangles only)
    integer, dimension(3)                 :: v        ! vertices in present cell
    integer                               :: vc       ! considered vertex
    integer                               :: vn       ! first upwave vertex of next cell
    integer                               :: vp       ! last upwave vertex of present cell
    !
    integer, dimension(:), allocatable    :: ivlist   ! list of index vertices
    !
    real                                  :: carea    ! area of the present cell (triangles only)
    real                                  :: dx1      ! first component of covariant base vector a_(1)
    real                                  :: dx2      ! second component of covariant base vector a_(1)
    real                                  :: dy1      ! first component of covariant base vector a_(2)
    real                                  :: dy2      ! second component of covariant base vector a_(2)
    real                                  :: fac      ! factor
    real                                  :: rdet     ! reciprocal of determinant
    real                                  :: rproc    ! auxiliary variable with percentage of cells which are non-Delaunay
    real                                  :: th1      ! direction of one face pointing to present vertex
    real                                  :: th2      ! direction of another face pointing to present vertex
    real                                  :: thdiff   ! difference between th2 and th1
    real, dimension(2)                    :: vec12    ! translation vector of coordinates: vertex2 - vertex1
    real, dimension(2)                    :: vec13    ! translation vector of coordinates: vertex3 - vertex1
    real                                  :: x1       ! x-coordinate of first vertex of present cell
    real                                  :: x2       ! x-coordinate of second vertex of present cell
    real                                  :: x3       ! x-coordinate of third vertex of present cell
    real                                  :: xc       ! x-coordinate of the cell centroid / circumcenter
    real                                  :: y1       ! y-coordinate of first vertex of present cell
    real                                  :: y2       ! y-coordinate of second vertex of present cell
    real                                  :: y3       ! y-coordinate of third vertex of present cell
    real                                  :: yc       ! y-coordinate of the cell centroid / circumcenter
    !
    logical                               :: nxtcell  ! indicate whether there is next cell in counterclockwise direction
    !
    character(120)                        :: msgstr   ! string to pass message
    !
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanGridCell')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ncc = 0
    !
    ! loop over all cells
    !
    do icell = 1, ncells
       !
       ! determine number of vertices and faces (triangles only!)
       !
       cell(icell)%nov = nov
       cell(icell)%nof = nov
       !
       ! identification number
       !
       cell(icell)%atti(CELLID) = icell
       !
       ! cell is triangle and initiatively active
       !
       cell(icell)%atti(CELLRECT) = 0
       cell(icell)%active         = .true.
       !
       ! store vertices of the cell
       !
       v(1) = kvertc(1,icell)
       v(2) = kvertc(2,icell)
       v(3) = kvertc(3,icell)
       cell(icell)%atti(CELLV1) = v(1)
       cell(icell)%atti(CELLV2) = v(2)
       cell(icell)%atti(CELLV3) = v(3)
       !
       ! store area of the cell in carea
       !
       vec12(1) = xcugrd(v(2)) - xcugrd(v(1))
       vec12(2) = ycugrd(v(2)) - ycugrd(v(1))
       vec13(1) = xcugrd(v(3)) - xcugrd(v(1))
       vec13(2) = ycugrd(v(3)) - ycugrd(v(1))
       carea = 0.5*abs(vec12(1)*vec13(2) - vec13(1)*vec12(2))
       cell(icell)%attr(CELLAREA) = carea
       !
       ! store local covariant and contravariant base vectors at each vertex
       ! next, store directions of faces pointing to present vertex
       !
       dx1 = -vec12(1)
       dy1 = -vec12(2)
       dx2 = -vec13(1)
       dy2 = -vec13(2)
       !
       do j = 1, nov
          !
          rdet = 1./(dy2*dx1 - dy1*dx2)
          !
          cell(icell)%geom(j)%det  =  dy2*dx1 - dy1*dx2
          cell(icell)%geom(j)%dx1  =  dx1
          cell(icell)%geom(j)%dy1  =  dy1
          cell(icell)%geom(j)%dx2  =  dx2
          cell(icell)%geom(j)%dy2  =  dy2
          cell(icell)%geom(j)%rdx1 =  dy2*rdet
          cell(icell)%geom(j)%rdy1 = -dx2*rdet
          cell(icell)%geom(j)%rdx2 = -dy1*rdet
          cell(icell)%geom(j)%rdy2 =  dx1*rdet
          !
          th1 = atan2(dy1,dx1)
          th2 = atan2(dy2,dx2)
          !
          thdiff = th1 - th2
          do
             if ( abs(thdiff) <= PI ) exit
             th1 = th1 - sign (2., thdiff) * PI
             thdiff = th1 - th2
          enddo
          !
          cell(icell)%geom(j)%th1 = th1
          cell(icell)%geom(j)%th2 = th2
          !
          dx1 = dx2 - dx1
          dy1 = dy2 - dy1
          dx2 = dx1 - dx2
          dy2 = dy1 - dy2
          !
       enddo
       !
       ! determine orientation of the mesh
       !
       if ( icell == 1 ) then
          !
          if ( dy2*dx1 > dy1*dx2 ) then
             CVLEFT = .false.             ! right-handed
          else
             CVLEFT = .true.              ! left-handed
          endif
          !
       endif
       !
       ! store coordinates of centroid
       !
       xc = 0.
       yc = 0.
       do j = 1, nov
          xc = xc + xcugrd(kvertc(j,icell))
          yc = yc + ycugrd(kvertc(j,icell))
       enddo
       xc = xc / real(nov)
       yc = yc / real(nov)
       !
       cell(icell)%attr(CELLCX) = xc
       cell(icell)%attr(CELLCY) = yc
       !
       ! store coordinates of circumcenter
       !
       x1 = xcugrd(v(1))
       y1 = ycugrd(v(1))
       x2 = xcugrd(v(2))
       y2 = ycugrd(v(2))
       x3 = xcugrd(v(3))
       y3 = ycugrd(v(3))
       fac = ( (x2-x3)*(x3-x1)+(y2-y3)*(y3-y1) )/( (x1-x2)*(y3-y1)-(y1-y2)*(x3-x1) )
       xc = 0.5 * ( fac*(y1-y2) + x1+x2 )
       yc = 0.5 * ( fac*(x2-x1) + y1+y2 )
       !
       if ( ccn ) then
          !
          ! check if circumcenter is in the cell
          !
          numcrs = 0
          !
          ! loop over faces of the cell
          !
          do j = 1, nov
             !
             k = mod(j,nov)+1
             !
             x1 = xcugrd(kvertc(j,icell))
             y1 = ycugrd(kvertc(j,icell))
             x2 = xcugrd(kvertc(k,icell))
             y2 = ycugrd(kvertc(k,icell))
             !
             if ( ( (x1 > xc) .and. (.not. x2 > xc) ) .or. ( (x2 > xc) .and. (.not. x1 > xc) ) ) then
                !
                if ( y1 > yc .or. y2 > yc ) then
                   !
                   fac = y1 + (xc-x1) * (y2-y1) / (x2-x1)
                   if ( fac > yc ) numcrs = numcrs + 1
                   !
                endif
                !
             endif
             !
          enddo
          !
       else
          !
          numcrs = 1
          !
       endif
       !
       ! if number of crossings is odd then circumcenter is inside the cell
       !
       if ( mod(numcrs,2) == 1 ) then
          cell(icell)%attr(CELLCCX) = xc
          cell(icell)%attr(CELLCCY) = yc
       else
          ! current cell is obtuse and will be corrected by replacing circumcenter by centroid
          ncc = ncc + 1
          cell(icell)%attr(CELLCCX) = cell(icell)%attr(CELLCX)
          cell(icell)%attr(CELLCCY) = cell(icell)%attr(CELLCY)
       endif
       !
       ! store vertices of each face of the cell
       !
       do j = 1, nov
          k = mod(j,nov)+1
          cell(icell)%face(j)%atti(FACEV1) = kvertc(j,icell)
          cell(icell)%face(j)%atti(FACEV2) = kvertc(k,icell)
       enddo
       !
    enddo
    !
    ! give info in case of non-Delaunay (obtuse) cells
    !
    if ( ncc > 0 .and. ITEST >= 30 ) then
       !
       rproc = 100.*real(ncc)/real(ncells)
       !
       write (msgstr,'(a,i6,a,i6,a,f5.1,a)') 'corrected ',ncc,' of ',ncells,' cells with angles > 90 degrees (',rproc,'%)'
       call msgerr (0, trim(msgstr) )
       !
    endif
    !
    allocate(ivlist(nverts))
    !
    ! loop over all vertices
    !
    do ivert = 1, nverts
       !
       ! identify the considered vertex and store index
       !
       vc = vert(ivert)%atti(VERTID)
       ivlist(vc) = ivert
       !
       ! initialize number of cells around vertex
       !
       vert(ivert)%noc = 0
       !
    enddo
    !
    do icell = 1, ncells
       !
       v(1) = cell(icell)%atti(CELLV1)
       v(2) = cell(icell)%atti(CELLV2)
       v(3) = cell(icell)%atti(CELLV3)
       !
       ! add present cell to each of these vertices
       !
       do j = 1, nov
          !
          ivert = ivlist(v(j))
          noc   = vert(ivert)%noc +1
          !
          vert(ivert)%noc = noc
          vert(ivert)%cell(noc)%atti(CELLID) = icell
          !
       enddo
       !
    enddo
    !
    deallocate(ivlist)
    !
    do ivert = 1, nverts
       !
       noc = vert(ivert)%noc
       !
       ! loop over cells around considered vertex
       !
       do j = 1, noc
          !
          ! get cell and its vertices
          !
          icell = vert(ivert)%cell(j)%atti(CELLID)
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          ! pick up last upwave vertex (counterclockwise counting of vertices is assumed)
          !
          do k = 1, nov
             if ( v(k) == ivert ) then
                vp = v(mod(k+1,nov)+1)
                exit
             endif
          enddo
          !
          ! search for next cell in counterclockwise direction
          !
          nxtcell = .false.
          !
          do l = 1, noc
             !
             ! get a cell and its vertices
             !
             jcell = vert(ivert)%cell(l)%atti(CELLID)
             !
             v(1) = cell(jcell)%atti(CELLV1)
             v(2) = cell(jcell)%atti(CELLV2)
             v(3) = cell(jcell)%atti(CELLV3)
             !
             ! pick up first upwave vertex (counterclockwise counting of vertices is assumed)
             !
             do k = 1, nov
                if ( v(k) == ivert ) then
                   vn = v(mod(k,nov)+1)
                   exit
                endif
             enddo
             !
             ! check whether first upwave vertex of next cell equals last upwave vertex of present cell
             !
             if ( vn == vp ) then
                nxtcell = .true.
                exit
             endif
             !
          enddo
          !
          if ( .not.nxtcell ) jcell = 0
          vert(ivert)%cell(j)%atti(NEXTCELL) = jcell
          !
       enddo
       !
    enddo

end subroutine SwanGridCell
