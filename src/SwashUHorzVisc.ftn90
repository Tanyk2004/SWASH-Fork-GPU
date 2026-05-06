subroutine SwashUHorzVisc ( u )
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
!   Calculates horizontal eddy viscosity coefficient in case of flexible mesh
!
!   Method
!
!   The mixing length model is employed to account for wave breaking
!   The Smagorinsky model is utilized to account for subgrid turbulent mixing
!
!   Note: if eddy viscosity is larger than maximum, which is based on stability criterion, apply clipping
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata, only: brks, hs, vnu2d, uvc, xshear
    use SwashTimeComm, only: dt
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    real, dimension(nfaces), intent(in)   :: u        ! flow velocity at faces
!
!   Local variables
!
    integer                               :: icell    ! index of present cell
    integer                               :: icelll   ! left cell of present face
    integer                               :: icellr   ! right cell of present face
    integer                               :: icistb   ! counter for number of instable faces
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: iface    ! loop counter over faces
    integer                               :: ivert    ! loop counter over vertices
    integer                               :: jc       ! loop counter
    integer                               :: jcell    ! index of next cell
    integer, dimension(2)                 :: v        ! vertices of present face
    !
    real                                  :: area     ! twices the area of centroid dual around present vertex
    real                                  :: areal    ! area of left cell of present face
    real                                  :: arear    ! area of right cell of present face
    real                                  :: dudx     ! x-component of gradient of u in vertex
    real                                  :: dvdx     ! x-component of gradient of v in vertex
    real                                  :: dudy     ! y-component of gradient of u in vertex
    real                                  :: dvdy     ! y-component of gradient of v in vertex
    real                                  :: finp     ! interpolation factor
    real                                  :: lf       ! length of face
    real                                  :: rdx      ! reciprocal of distance between circumcenters adjacent to face
    real                                  :: rproc    ! auxiliary variable with percentage of instable faces
    real                                  :: shrm     ! shear squared at mid face
    real                                  :: stabmx   ! auxiliary variable with maximum viscosity based stability criterion
    real                                  :: u0       ! u-velocity component in centroid of present cell
    real                                  :: u1       ! u-velocity component in centroid of next cell
    real                                  :: v0       ! v-velocity component in centroid of present cell
    real                                  :: v1       ! v-velocity component in centroid of next cell
    real                                  :: x0       ! x-coordinate of the centroid of present cell
    real                                  :: x1       ! x-coordinate of the centroid of next cell
    real                                  :: y0       ! y-coordinate of the centroid of present cell
    real                                  :: y1       ! y-coordinate of the centroid of next cell
    !
    character(80)                         :: msgstr   ! string to pass message
    !
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face     ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUHorzVisc')
    !
    ! point to vertex, cell and face objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    icistb = 0
    !
    ! determine cell-based velocity components
    !
    call perot ( u, 1, 1 )
    !
    ! first, compute magnitude of shear squared at internal vertices, if appropriate
    !
    if ( ihvisc > 1 ) then
       !
       xshear = 0.
       !
       vertexloop : do ivert = 1, nverts
          !
          if ( vert(ivert)%atti(VMARKER) == 1 ) cycle vertexloop    ! boundary vertex
          !
          ! compute gradients of velocities in vertices
          !
          area = 0.
          dudx = 0.
          dvdx = 0.
          dudy = 0.
          dvdy = 0.
          !
          ! loop over cells around considered vertex
          !
          do jc = 1, vert(ivert)%noc
             !
             ! get present cell
             !
             icell = vert(ivert)%cell(jc)%atti(CELLID)
             !
             ! determine centroid of present cell
             !
             x0 = cell(icell)%attr(CELLCX)
             y0 = cell(icell)%attr(CELLCY)
             !
             ! get velocity components in centroid of present cell
             !
             u0 = uvc(icell,1,1)
             v0 = uvc(icell,1,2)
             !
             ! get next cell in counterclockwise direction
             !
             jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
             !
             ! determine centroid of next cell
             !
             x1 = cell(jcell)%attr(CELLCX)
             y1 = cell(jcell)%attr(CELLCY)
             !
             ! get velocity components in centroid of next cell
             !
             u1 = uvc(jcell,1,1)
             v1 = uvc(jcell,1,2)
             !
             ! compute contribution to area of centroid dual
             !
             area = area + x0*y1 - x1*y0
             !
             ! compute contribution to x-gradient of u- and v-velocity components
             !
             dudx = dudx + ( u0 + u1 ) * ( y1 - y0 )
             dvdx = dvdx + ( v0 + v1 ) * ( y1 - y0 )
             !
             ! compute contribution to y-gradient of u- and v-velocity components
             !
             dudy = dudy + ( u0 + u1 ) * ( x1 - x0 )
             dvdy = dvdy + ( v0 + v1 ) * ( x1 - x0 )
             !
          enddo
          !
          ! if area is non-positive, give error and go to next vertex
          !
          if ( .not. area > 0. ) then
             write (msgstr, '(a,i8)') ' Area of centroid dual is negative or zero in vertex ', ivert
             call msgerr( 2, trim(msgstr) )
             cycle vertexloop
          endif
          !
          dudx =  dudx / area
          dvdx =  dvdx / area
          dudy = -dudy / area
          dvdy = -dvdy / area
          !
          ! compute magnitude of shear squared
          !
          xshear(ivert) = 2. * ( dudx*dudx + dvdy*dvdy ) + ( dudy + dvdx )*( dudy + dvdx )
          !
       enddo vertexloop
       !
    endif
    !
    ! compute eddy viscosity coefficients at internal faces only
    !
    if ( ihvisc == 1 .and. hvisc > 0. ) then
       !
       ! constant eddy viscosity
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             rdx = face(iface)%attr(FACEDISTC)
             lf  = face(iface)%attr(FACELEN)
             !
             areal = cell(icelll)%attr(CELLAREA)
             arear = cell(icellr)%attr(CELLAREA)
             !
             vnu2d(iface) = hvisc
             !
             stabmx = ( finp*areal + (1.-finp)*arear ) / ( 3. * rdx * lf * dt )
             !
             if ( .not. vnu2d(iface) < stabmx ) then
                vnu2d(iface) = stabmx
                icistb       = icistb + 1
             endif
             !
          endif
          !
       enddo
       !
    else if ( ihvisc == 2 ) then
       !
       ! Smagorinsky model
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             v(1) = face(iface)%atti(FACEV1)
             v(2) = face(iface)%atti(FACEV2)
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             rdx = face(iface)%attr(FACEDISTC)
             lf  = face(iface)%attr(FACELEN)
             !
             areal = cell(icelll)%attr(CELLAREA)
             arear = cell(icellr)%attr(CELLAREA)
             !
             shrm = 0.5 * ( xshear(v(1)) + xshear(v(2)) )
             !
             vnu2d(iface) = csmag * csmag * lf * sqrt(shrm) / rdx
             !
             stabmx = ( finp*areal + (1.-finp)*arear ) / ( 3. * rdx * lf * dt )
             !
             if ( .not. vnu2d(iface) < stabmx ) then
                vnu2d(iface) = stabmx
                icistb       = icistb + 1
             endif
             !
          endif
          !
       enddo
       !
    else if ( ihvisc == 3 ) then
       !
       ! mixing length model
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             v(1) = face(iface)%atti(FACEV1)
             v(2) = face(iface)%atti(FACEV2)
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             rdx = face(iface)%attr(FACEDISTC)
             lf  = face(iface)%attr(FACELEN)
             !
             areal = cell(icelll)%attr(CELLAREA)
             arear = cell(icellr)%attr(CELLAREA)
             !
             shrm = 0.5 * ( xshear(v(1)) + xshear(v(2)) )
             !
             vnu2d(iface) = lmix * lmix * sqrt(shrm)
             !
             stabmx = ( finp*areal + (1.-finp)*arear ) / ( 3. * rdx * lf * dt )
             !
             if ( .not. vnu2d(iface) < stabmx ) then
                vnu2d(iface) = stabmx
                icistb       = icistb + 1
             endif
             !
          endif
          !
       enddo
       !
    else if ( ihvisc == 4 ) then
       !
       ! HFA model
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             v(1) = face(iface)%atti(FACEV1)
             v(2) = face(iface)%atti(FACEV2)
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             rdx = face(iface)%attr(FACEDISTC)
             lf  = face(iface)%attr(FACELEN)
             !
             areal = cell(icelll)%attr(CELLAREA)
             arear = cell(icellr)%attr(CELLAREA)
             !
             if ( brks(icelll) + brks(icellr) /= 0 ) then
                !
                shrm = 0.5 * ( xshear(v(1)) + xshear(v(2)) )
                !
                vnu2d(iface) = psurf(3) * (finp*hs(icelll) + (1.-finp)*hs(icellr))**2 * sqrt(shrm)
                !
             else
                !
                vnu2d(iface) = 0.
                !
             endif
             !
             stabmx = ( finp*areal + (1.-finp)*arear ) / ( 3. * rdx * lf * dt )
             !
             if ( .not. vnu2d(iface) < stabmx ) then
                vnu2d(iface) = stabmx
                icistb       = icistb + 1
             endif
             !
          endif
          !
       enddo
       !
    endif
    !
    ! give warning for instable faces
    !
    if ( icistb > 0 ) then
       !
       rproc = 100.*real(icistb)/real(nfaces - nfacesb)
       !
       if ( .not. rproc < 1. ) then
          write (msgstr,'(a,f5.1)') 'percentage of instable points for computing horizontal eddy viscosity = ',rproc
          call msgerr (1, trim(msgstr) )
       endif
       !
    endif
    !
end subroutine SwashUHorzVisc
