subroutine SwashQuanOutp ( rtype, oqproc, mip, nvoqp, nvoqk, xc, yc, voqr, voq, voqk, ionod )
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
!    1.00,   March 2010: New subroutine
!    1.05, January 2012: extension logarithmic wall-law
!    1.10,   March 2012: extension time-averaged quantities
!    7.01, January 2020: extension to unstructured grids
!
!   Purpose
!
!   Calculates requested output quantities in the output points
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use SwashTimeComm
    use outp_data
    use m_genarr
    use m_parall
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                            :: mip      ! number of output points
    integer, intent(in)                            :: nvoqk    ! number of layer-dependent quantities per output point
    integer, intent(in)                            :: nvoqp    ! number of quantities per output point
    !
    integer, dimension(mip), intent(out)           :: ionod    ! array indicating in which subdomain output points are located
    integer, dimension(nmovar), intent(in)         :: voqr     ! place of each output quantity
    !
    real, dimension(mip,nvoqp), intent(out)        :: voq      ! interpolated output quantity at request
    real, dimension(mip,0:kmax,nvoqk), intent(out) :: voqk     ! interpolated layer-dependent output quantity at request
    real, dimension(mip), intent(in)               :: xc       ! computational coordinates (or broken coordinates) of output point in x-direction
    real, dimension(mip), intent(in)               :: yc       ! computational coordinates (or broken coordinates) of output point in y-direction
    !
    character(4), intent(in)                       :: rtype    ! type of output request
    !
    logical, dimension(nmovar), intent(in)         :: oqproc   ! indicates whether or not an output quantity must be processed
!
!   Local variables
!
    integer                               :: i        ! loop counter
    integer                               :: icell    ! cell index / loop counter over cells
    integer                               :: icelll   ! left cell of present face
    integer                               :: icellr   ! right cell of present face
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: iface    ! face index
    integer                               :: ihlx     ! width of halo area in x-direction
    integer                               :: ihly     ! width of halo area in y-direction
    integer                               :: indx     ! point index
    integer                               :: indxb    ! index of point below from point of consideration
    integer                               :: indxl    ! index of point left from point of consideration
    integer                               :: indxr    ! index of point right from point of consideration
    integer                               :: indxs    ! index of point right-up from point of consideration
    integer                               :: indxu    ! index of point up from point of consderation
    integer                               :: ivert    ! loop counter over vertices
    integer                               :: ix       ! point index in x-direction
    integer                               :: ixb      ! start point index in x-direction
    integer                               :: ixe      ! end point index in x-direction
    integer                               :: iy       ! point index in y-direction
    integer                               :: iyb      ! start point index in y-direction
    integer                               :: iye      ! end point index in y-direction
    integer                               :: j        ! loop counter
    integer                               :: jc       ! loop counter
    integer                               :: jf       ! loop counter
    integer                               :: jvx      ! index of u-velocity component
    integer                               :: jvy      ! index of v-velocity component
    integer                               :: k        ! loop counter over vertical layers
    integer, dimension(:), allocatable    :: kvert    ! closest vertex indices of output points
    !
    real                                  :: area     ! area of present cell
    real                                  :: DEGCNV   ! indicates Cartesian or nautical degrees
    real                                  :: dirdeg   ! direction in degrees
    real                                  :: dxdetab  ! derivative of x to eta below from point of consideration
    real                                  :: dxdetau  ! derivative of x to eta up from point of consideration
    real                                  :: dxdksil  ! derivative of x to ksi left from point of consideration
    real                                  :: dxdksir  ! derivative of x to ksi right from point of consideration
    real                                  :: dydetab  ! derivative of y to eta below from point of consideration
    real                                  :: dydetau  ! derivative of y to eta up from point of consideration
    real                                  :: dydksil  ! derivative of y to ksi left from point of consideration
    real                                  :: dydksir  ! derivative of y to ksi right from point of consideration
    real                                  :: dxl      ! local mesh size in x-direction
    real                                  :: dyl      ! local mesh size in y-direction
    real                                  :: fac      ! a factor
    real                                  :: fac1     ! another factor
    real                                  :: fac2     ! some other factor
    real                                  :: lf       ! length of face
    real                                  :: mub      ! time-averaged u-velocity below from point of consideration
    real                                  :: muu      ! time-averaged u-velocity up from point of consideration
    real                                  :: mvl      ! time-averaged v-velocity left from point of consideration
    real                                  :: mvr      ! time-averaged v-velocity right from point of consideration
    real                                  :: nx       ! x-component of normal to face
    real                                  :: ny       ! y-component of normal to face
    real                                  :: qf       ! sum of mass fluxes over faces of present cell
    real                                  :: qxb      ! u-discharge below from point of consideration
    real                                  :: qxu      ! u-discharge up from point of consideration
    real                                  :: qyl      ! v-discharge left from point of consideration
    real                                  :: qyr      ! v-discharge right from point of consideration
    real                                  :: rdist    ! a distance
    real                                  :: rdx      ! increment in x-direction
    real                                  :: rdy      ! increment in y-direction
    real                                  :: rsgn     ! sign for indicating face orientation
    real                                  :: rval1    ! a value
    real                                  :: rval2    ! another value
    real                                  :: totarea  ! total area of all cells around vertex
    real                                  :: u        ! u-velocity at point different from its point of definition
    real                                  :: ub       ! u-velocity below from point of consideration
    real                                  :: ufb      ! u-component of friction velocity below from point of consideration
    real                                  :: ufu      ! u-component of friction velocity up from point of consideration
    real                                  :: ugradz   ! u times grad z
    real                                  :: uu       ! u-velocity up from point of consideration
    real                                  :: v        ! v-velocity at point different from its point of definition
    real                                  :: vfl      ! v-component of friction velocity left from point of consideration
    real                                  :: vfr      ! v-component of friction velocity right from point of consideration
    real                                  :: vl       ! v-velocity left from point of consideration
    real                                  :: vr       ! v-velocity right from point of consideration
    real                                  :: xp       ! local x-coordinate
    real                                  :: yp       ! local y-coordinate
    !
    logical                               :: EQREAL   ! compares two reals
    logical, dimension(:), allocatable    :: ltmp     ! a temporary array of type logical
    logical                               :: reqvel   ! velocity/discharge components are requested
    logical                               :: reqzsc   ! layer-dependent scalars are requested
    logical                               :: reqzvl   ! layer-dependent velocity/discharge components are requested
    !
    type(facetype), dimension(:), pointer :: face     ! datastructure for faces with their attributes
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashQuanOutp')
    !
    ! point to face, cell and vertex objects
    !
    face => gridobject%face_grid
    cell => gridobject%cell_grid
    vert => gridobject%vert_grid
    !
    if ( ITEST >= 100 .or. IOUTES >= 10 ) write (PRTEST, 201) ( oqproc(i), i=1,10 ), mip
    !
    ! initialize help arrays for outputting
    !
    oparr = 0.
    if ( kmax > 1 ) oparrk = 0.
    !
    ! asked for layer-dependent scalars (e.g. vertical velocities and pressure) or not
    !
    if ( oqproc(52) .or. oqproc(53) .or. oqproc(54) .or. oqproc(55) .or. oqproc(56) .or.  &
         oqproc(57) .or. oqproc(58) .or. oqproc(59) .or. oqproc(82) .or. oqproc(83) .or.  &
         oqproc(89) .or. oqproc(90) .or. oqproc(91) .or. oqproc(92) .or. oqproc(93) .or.  &
         oqproc(94) ) then
       reqzsc = .true.
    else
       reqzsc = .false.
    endif
    !
    ! asked for layer-dependent velocity/discharge components or not
    !
    if ( oqproc(72) .or. oqproc(73) .or. oqproc(74) .or. oqproc(75) .or. oqproc(76) .or.  &
         oqproc(77) .or. oqproc(78) .or. oqproc(79) .or. oqproc(80) .or. oqproc(81) .or.  &
         oqproc(84) .or. oqproc(85) .or. oqproc(86) .or. oqproc(87) .or. oqproc(88) ) then
       reqzvl = .true.
    else
       reqzvl = .false.
    endif
    !
    ! asked for velocity/discharge components or not
    !
    if ( oqproc( 7) .or. oqproc( 8) .or. oqproc( 9) .or. oqproc(10) .or. oqproc(11) .or.  &
         oqproc(14) .or. oqproc(15) .or. oqproc(16) .or. oqproc(17) .or. oqproc(18) .or.  &
         oqproc(31) .or. oqproc(32) .or. oqproc(33) .or. oqproc(34) .or. oqproc(35) .or.  &
         oqproc(36) .or. oqproc(37) .or. oqproc(45) .or.                                  &
         reqzvl ) then
       reqvel = .true.
    else
       reqvel = .false.
    endif
    !
    ! compute the depth-averaged velocity components and specific discharges, if requested
    !
    if ( reqvel .and. optg /= 5 ) then
       !
       if ( oned ) then
          !
          do i = 1, mxc
             !
             indx = kgrpnt(i,1)
             !
             if ( kmax == 1 ) then
                !
                oparr(indx,3) = u1(indx,1)
                !
             else
                !
                if ( hu(indx) > 0. ) then
                   !
                   do k = 1, kmax
                      !
                      oparr(indx,3) = oparr(indx,3) + hku(indx,k)*u1(indx,k)
                      !
                   enddo
                   !
                   oparr(indx,3) = oparr(indx,3) / hu(indx)
                   !
                endif
                !
             endif
             !
          enddo
          !
          oparr(:,7) = oparr(:,3)*hu(:)
          !
          if ( lcuroutp ) then
             !
             do i = 1, mxc
                !
                indx = kgrpnt(i,1)
                !
                if ( kmax == 1 ) then
                   !
                   oparr(indx,20) = mvelu(indx,1)
                   !
                else
                   !
                   if ( hu(indx) > 0. ) then
                      !
                      do k = 1, kmax
                         !
                         oparr(indx,20) = oparr(indx,20) + hku(indx,k)*mvelu(indx,k)
                         !
                      enddo
                      !
                      oparr(indx,20) = oparr(indx,20) / hu(indx)
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          endif
          !
       else
          !
          do j = 1, myc
             do i = 1, mxc
                !
                indx = kgrpnt(i,j)
                !
                if ( kmax == 1 ) then
                   !
                   oparr(indx,3) = u1(indx,1)
                   !
                else
                   !
                   if ( hu(indx) > 0. ) then
                      !
                      do k = 1, kmax
                         !
                         oparr(indx,3) = oparr(indx,3) + hku(indx,k)*u1(indx,k)
                         !
                      enddo
                      !
                      oparr(indx,3) = oparr(indx,3) / hu(indx)
                      !
                   endif
                   !
                endif
                !
             enddo
          enddo
          !
          oparr(:,7) = oparr(:,3)*hu(:)
          !
          do j = 1, myc
             do i = 1, mxc
                !
                indx = kgrpnt(i,j)
                !
                if ( kmax == 1 ) then
                   !
                   oparr(indx,4) = v1(indx,1)
                   !
                else
                   !
                   if ( hv(indx) > 0. ) then
                      !
                      do k = 1, kmax
                         !
                         oparr(indx,4) = oparr(indx,4) + hkv(indx,k)*v1(indx,k)
                         !
                      enddo
                      !
                      oparr(indx,4) = oparr(indx,4) / hv(indx)
                      !
                   endif
                   !
                endif
                !
             enddo
          enddo
          !
          oparr(:,8) = oparr(:,4)*hv(:)
          !
          if ( lcuroutp ) then
             !
             do j = 1, myc
                do i = 1, mxc
                   !
                   indx = kgrpnt(i,j)
                   !
                   if ( kmax == 1 ) then
                      !
                      oparr(indx,20) = mvelu(indx,1)
                      !
                   else
                      !
                      if ( hu(indx) > 0. ) then
                         !
                         do k = 1, kmax
                            !
                            oparr(indx,20) = oparr(indx,20) + hku(indx,k)*mvelu(indx,k)
                            !
                         enddo
                         !
                         oparr(indx,20) = oparr(indx,20) / hu(indx)
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
             do j = 1, myc
                do i = 1, mxc
                   !
                   indx = kgrpnt(i,j)
                   !
                   if ( kmax == 1 ) then
                      !
                      oparr(indx,21) = mvelv(indx,1)
                      !
                   else
                      !
                      if ( hv(indx) > 0. ) then
                         !
                         do k = 1, kmax
                            !
                            oparr(indx,21) = oparr(indx,21) + hkv(indx,k)*mvelv(indx,k)
                            !
                         enddo
                         !
                         oparr(indx,21) = oparr(indx,21) / hv(indx)
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
          endif
          !
       endif
       !
    endif
    !
    ! compute physical vertical velocity in case of hydrostatic computation
    !
    if ( oqproc(52) .and. (ihydro == 0 .or. ihydro == 3) ) then
       !
       if ( optg /= 5 ) then
          !
          if ( oned ) then
             !
             do i = 2, mxc - 1
                !
                indx  = kgrpnt(i  ,1)
                indxl = kgrpnt(i-1,1)
                !
                if ( wets(indx) == 1 ) then
                   !
                   u = 0.5 * ( u1(indx,kmax) + u1(indxl,kmax) )
                   !
                   wphy(indx,kmax) = u * ( zku(indx,kmax) - zku(indxl,kmax) ) / dx
                   !
                   fac = - ( hku(indx,kmax)*u1(indx,kmax) - hku(indxl,kmax)*u1(indxl,kmax) ) / dx
                   !
                   do k = kmax-1, 1, -1
                      !
                      fac1 = 0.5 * ( u1(indx,k  ) + u1(indxl,k  ) )
                      fac2 = 0.5 * ( u1(indx,k+1) + u1(indxl,k+1) )
                      !
                      u = ( fac2 * hks(indx,k) + fac1 * hks(indx,k+1) ) / ( hks(indx,k) + hks(indx,k+1) )
                      !
                      wphy(indx,k) = fac + u * ( zku(indx,k) - zku(indxl,k) ) / dx
                      !
                      fac = fac - ( hku(indx,k)*u1(indx,k) - hku(indxl,k)*u1(indxl,k) ) / dx
                      !
                   enddo
                   !
                   u = 0.5 * ( u1(indx,1) + u1(indxl,1) )
                   !
                   wphy(indx,0) = fac + u * ( zku(indx,0) - zku(indxl,0) ) / dx
                   !
                else
                   !
                   wphy(indx,:) = 0.
                   !
                endif
                !
             enddo
             !
             wphy(kgrpnt(1  ,1),:) = wphy(kgrpnt(2    ,1),:)
             wphy(kgrpnt(mxc,1),:) = wphy(kgrpnt(mxc-1,1),:)
             !
          else
             !
             do j = 2, myc - 1
                !
                do i = 2, mxc - 1
                   !
                   indx  = kgrpnt(i  ,j  )
                   indxl = kgrpnt(i-1,j  )
                   indxb = kgrpnt(i  ,j-1)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( indxl == 1 ) indxl = indx
                   if ( indxb == 1 ) indxb = indx
                   !
                   if ( wets(indx) == 1 ) then
                      !
                      dxl = 0.5 * ( gvv(indx) + gvv(indxb) )
                      dyl = 0.5 * ( guu(indx) + guu(indxl) )
                      !
                      u = 0.5 * ( u1(indx,kmax) + u1(indxl,kmax) )
                      v = 0.5 * ( v1(indx,kmax) + v1(indxb,kmax) )
                      !
                      wphy(indx,kmax) = u * ( zku(indx,kmax) - zku(indxl,kmax) ) / dxl + v * ( zkv(indx,kmax) - zkv(indxb,kmax) ) / dyl
                      !
                      fac = - ( hku(indx,kmax)*guu(indx)*u1(indx,kmax) - hku(indxl,kmax)*guu(indxl)*u1(indxl,kmax) +  &
                                hkv(indx,kmax)*gvv(indx)*v1(indx,kmax) - hkv(indxb,kmax)*gvv(indxb)*v1(indxb,kmax) ) / gsqs(indx)
                      !
                      do k = kmax-1, 1, -1
                         !
                         fac1 = 0.5 * ( u1(indx,k  ) + u1(indxl,k  ) )
                         fac2 = 0.5 * ( u1(indx,k+1) + u1(indxl,k+1) )
                         !
                         u = ( fac2 * hks(indx,k) + fac1 * hks(indx,k+1) ) / ( hks(indx,k) + hks(indx,k+1) )
                         !
                         fac1 = 0.5 * ( v1(indx,k  ) + v1(indxb,k  ) )
                         fac2 = 0.5 * ( v1(indx,k+1) + v1(indxb,k+1) )
                         !
                         v = ( fac2 * hks(indx,k) + fac1 * hks(indx,k+1) ) / ( hks(indx,k) + hks(indx,k+1) )
                         !
                         wphy(indx,k) = fac + u * ( zku(indx,k) - zku(indxl,k) ) / dxl + v * ( zkv(indx,k) - zkv(indxb,k) ) / dyl
                         !
                         fac = fac - ( hku(indx,k)*guu(indx)*u1(indx,k) - hku(indxl,k)*guu(indxl)*u1(indxl,k) +  &
                                       hkv(indx,k)*gvv(indx)*v1(indx,k) - hkv(indxb,k)*gvv(indxb)*v1(indxb,k) ) / gsqs(indx)
                         !
                      enddo
                      !
                      u = 0.5 * ( u1(indx,1) + u1(indxl,1) )
                      v = 0.5 * ( v1(indx,1) + v1(indxb,1) )
                      !
                      wphy(indx,0) = fac + u * ( zku(indx,0) - zku(indxl,0) ) / dxl + v * ( zkv(indx,0) - zkv(indxb,0) ) / dyl
                      !
                   else
                      !
                      wphy(indx,:) = 0.
                      !
                   endif
                   !
                enddo
                !
                if ( .not.lreptx ) then
                   !
                   wphy(kgrpnt(1  ,j),:) = wphy(kgrpnt(2    ,j),:)
                   wphy(kgrpnt(mxc,j),:) = wphy(kgrpnt(mxc-1,j),:)
                   !
                else
                   !
                   call periodic ( wphy, kgrpnt, 0, kmax )
                   !
                endif
                !
             enddo
             !
             if ( .not.lrepty ) then
                !
                do i = 1, mxc
                   !
                   wphy(kgrpnt(i,1  ),:) = wphy(kgrpnt(i,2    ),:)
                   wphy(kgrpnt(i,myc),:) = wphy(kgrpnt(i,myc-1),:)
                   !
                enddo
                !
             else
                !
                call periodic ( wphy, kgrpnt, 0, kmax )
                !
             endif
             !
          endif
          !
       else
          !
          ! compute cell-based velocity vector
          !
          call perot ( u1, 1, kmax )
          !
          do icell = 1, ncells
             !
             if ( wets(icell) == 1 ) then
                !
                u = uvc(icell,kmax,1)
                v = uvc(icell,kmax,2)
                !
                ugradz = 0.
                qf     = 0.
                !
                do jf = 1, cell(icell)%nof
                   !
                   ! face identifier
                   !
                   iface = cell(icell)%face(jf)%atti(FACEID)
                   !
                   ! get length of current face
                   !
                   lf = face(iface)%attr(FACELEN)
                   !
                   nx = face(iface)%attr(FACENORMX)
                   ny = face(iface)%attr(FACENORMY)
                   !
                   ! consider left and right cells of current face
                   !
                   icelll = face(iface)%atti(FACECL)
                   icellr = face(iface)%atti(FACECR)
                   !
                   ! take into account orientation of the current face
                   !
                   if ( icell == icelll ) then
                      rsgn =  1.
                   else if ( icell == icellr ) then
                      rsgn = -1.
                   endif
                   !
                   ugradz = ugradz + rsgn * lf * zku(iface,kmax) * ( nx * u + ny * v )
                   qf     = qf     + rsgn * lf * hku(iface,kmax) * u1(iface,kmax)
                   !
                enddo
                !
                ! area of cell
                !
                area = cell(icell)%attr(CELLAREA)
                !
                wphy(icell,kmax) = ugradz / area
                !
                fac = -qf / area
                !
                do k = kmax-1, 1, -1
                   !
                   u = ( uvc(icell,k+1,1) * hks(icell,k) + uvc(icell,k,1) * hks(icell,k+1) ) / ( hks(icell,k) + hks(icell,k+1) )
                   v = ( uvc(icell,k+1,2) * hks(icell,k) + uvc(icell,k,2) * hks(icell,k+1) ) / ( hks(icell,k) + hks(icell,k+1) )
                   !
                   ugradz = 0.
                   qf     = 0.
                   !
                   do jf = 1, cell(icell)%nof
                      !
                      ! face identifier
                      !
                      iface = cell(icell)%face(jf)%atti(FACEID)
                      !
                      ! get length of current face
                      !
                      lf = face(iface)%attr(FACELEN)
                      !
                      nx = face(iface)%attr(FACENORMX)
                      ny = face(iface)%attr(FACENORMY)
                      !
                      ! consider left and right cells of current face
                      !
                      icelll = face(iface)%atti(FACECL)
                      icellr = face(iface)%atti(FACECR)
                      !
                      ! take into account orientation of the current face
                      !
                      if ( icell == icelll ) then
                         rsgn =  1.
                      else if ( icell == icellr ) then
                         rsgn = -1.
                      endif
                      !
                      ugradz = ugradz + rsgn * lf * zku(iface,k) * ( nx * u + ny * v )
                      qf     = qf     + rsgn * lf * hku(iface,k) * u1(iface,k)
                      !
                   enddo
                   !
                   ! area of cell
                   !
                   area = cell(icell)%attr(CELLAREA)
                   !
                   wphy(icell,k) = fac + ugradz / area
                   !
                   fac = fac - qf / area
                   !
                enddo
                !
                u = uvc(icell,1,1)
                v = uvc(icell,1,2)
                !
                ugradz = 0.
                !
                do jf = 1, cell(icell)%nof
                   !
                   ! face identifier
                   !
                   iface = cell(icell)%face(jf)%atti(FACEID)
                   !
                   ! get length of current face
                   !
                   lf = face(iface)%attr(FACELEN)
                   !
                   nx = face(iface)%attr(FACENORMX)
                   ny = face(iface)%attr(FACENORMY)
                   !
                   ! consider left and right cells of current face
                   !
                   icelll = face(iface)%atti(FACECL)
                   icellr = face(iface)%atti(FACECR)
                   !
                   ! take into account orientation of the current face
                   !
                   if ( icell == icelll ) then
                      rsgn =  1.
                   else if ( icell == icellr ) then
                      rsgn = -1.
                   endif
                   !
                   ugradz = ugradz + rsgn * lf * zku(iface,0) * ( nx * u + ny * v )
                   !
                enddo
                !
                ! area of cell
                !
                area = cell(icell)%attr(CELLAREA)
                !
                wphy(icell,0) = fac + ugradz / area
                !
             else
                !
                wphy(icell,:) = 0.
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
    ! determine output quantities in corner points of computational grid or no interpolation
    !
    if ( optg /= 5 ) then
       !
       if ( oned ) then
          !
          if ( EQREAL(outpar(3),0.) ) then
             !
             do i = 1, mxc - 1
                !
                indx  = kgrpnt(i  ,1)
                indxr = kgrpnt(i+1,1)
                !
                ! water level, water depth and inundation depth
                !
                oparr(indx,2) = 0.5 * ( s1(indx) + s1(indxr) )
                !
                oparr(indx,1) = oparr(indx,2) + depf(indx)
                if ( ifloat /= 0 ) oparr(indx,1) = min( depf(indx)-flobjf(indx), oparr(indx,1) )
                !
                oparr(indx,16) = 0.5 * ( hindun(indx) + hindun(indxr) )
                !
                ! friction velocity based on logarithmic wall-law
                !
                if ( irough == 4 .and. kmax > 1 ) then
                   !
                   oparr(indx,11) = 0.5 * ( logfrc(indx,1)*logfrc(indx,2) + logfrc(indxr,1)*logfrc(indxr,2) )
                   !
                endif
                !
                ! non-hydrostatic pressure at bottom
                !
                if ( ihydro /= 0 ) then
                   !
                   if ( kmax == 1 ) then
                      !
                      oparr(indx,13) = 0.5 * ( q(indx,1) + q(indxr,1) )
                      !
                   else
                      !
                      if ( ihydro == 3 ) then
                         !
                         oparr(indx,13) = 0.5 * ( qbot(indx) + qbot(indxr) )
                         !
                      else
                         !
                         oparr(indx,13) = 0.5 * ( q(indx,kpmax) + q(indxr,kpmax) )
                         !
                      endif
                      !
                   endif
                   !
                endif
                !
                ! wave parameters
                !
                if ( lwavoutp ) then
                   !
                   oparr(indx,17) = 0.5 * ( setup(indx) + setup(indxr) )
                   oparr(indx,18) = 0.5 * (  hsig(indx) +  hsig(indxr) )
                   !
                endif
                !
                oparr(indx,24) = 0.5 * ( real(brks(indx)) + real(brks(indxr)) )
                !
                ! salinity, temperature and sediment
                !
                if ( lsal  > 0 ) oparr(indx,25) = 0.5 * ( rp(indx,1,lsal ) + rp(indxr,1,lsal ) )
                if ( ltemp > 0 ) oparr(indx,26) = 0.5 * ( rp(indx,1,ltemp) + rp(indxr,1,ltemp) )
                if ( lsed  > 0 ) oparr(indx,27) = 0.5 * ( rp(indx,1,lsed ) + rp(indxr,1,lsed ) )
                !
                ! time-averaged or mean salinity, temperature and sediment
                !
                if ( ltraoutp ) then
                   !
                   if ( lsal  > 0 ) oparr(indx,28) = 0.5 * ( mcons(indx,1,lsal ) + mcons(indxr,1,lsal ) )
                   if ( ltemp > 0 ) oparr(indx,29) = 0.5 * ( mcons(indx,1,ltemp) + mcons(indxr,1,ltemp) )
                   if ( lsed  > 0 ) oparr(indx,30) = 0.5 * ( mcons(indx,1,lsed ) + mcons(indxr,1,lsed ) )
                   !
                endif
                !
                ! layer-dependent scalars, if requested
                !
                if ( reqzsc ) then
                   !
                   oparrk(indx,0:kmax,11) = 0.5 * ( wom(indx,0:kmax) + wom(indxr,0:kmax) )
                   !
                   if ( ihydro == 0 .or. ihydro == 3 ) then
                      !
                      oparrk(indx,0:kmax,12) = 0.5 * ( wphy(indx,0:kmax) + wphy(indxr,0:kmax) )
                      !
                   else if ( ihydro == 1 .or. ihydro == 2 ) then
                      !
                      oparrk(indx,0:kmax,12) = 0.5 * ( w1(indx,0:kmax) + w1(indxr,0:kmax) )
                      if ( .not.lsubg ) then
                         oparrk(indx,1:kmax,13) = 0.5 * (  q(indx,1:kmax) +  q(indxr,1:kmax) )
                      else
                         oparrk(indx,1:kmax,13) = 0.5 * (  qv(indx,1:kmax) +  qv(indxr,1:kmax) )
                      endif
                      !
                   endif
                   !
                   if ( iturb /= 0 ) then
                      !
                      oparrk(indx,0:kmax,15) = 0.5 * ( rtur(indx,0:kmax,1) + rtur(indxr,0:kmax,1) )
                      oparrk(indx,0:kmax,16) = 0.5 * ( rtur(indx,0:kmax,2) + rtur(indxr,0:kmax,2) )
                      !
                   endif
                   !
                   oparrk(indx,0:kmax,17) = 0.5 * ( vnu3d(indx,0:kmax) + vnu3d(indxr,0:kmax) )
                   !
                   if ( lturoutp ) then
                      !
                      oparrk(indx,0:kmax,22) = 0.5 * ( mtke(indx,0:kmax) + mtke(indxr,0:kmax) )
                      oparrk(indx,0:kmax,23) = 0.5 * ( meps(indx,0:kmax) + meps(indxr,0:kmax) )
                      !
                   endif
                   !
                   if ( lsal  > 0 ) oparrk(indx,1:kmax,25) = 0.5 * ( rp(indx,1:kmax,lsal ) + rp(indxr,1:kmax,lsal ) )
                   if ( ltemp > 0 ) oparrk(indx,1:kmax,26) = 0.5 * ( rp(indx,1:kmax,ltemp) + rp(indxr,1:kmax,ltemp) )
                   if ( lsed  > 0 ) oparrk(indx,1:kmax,27) = 0.5 * ( rp(indx,1:kmax,lsed ) + rp(indxr,1:kmax,lsed ) )
                   !
                   if ( ltraoutp ) then
                      !
                      if ( lsal  > 0 ) oparrk(indx,1:kmax,28) = 0.5 * ( mcons(indx,1:kmax,lsal ) + mcons(indxr,1:kmax,lsal ) )
                      if ( ltemp > 0 ) oparrk(indx,1:kmax,29) = 0.5 * ( mcons(indx,1:kmax,ltemp) + mcons(indxr,1:kmax,ltemp) )
                      if ( lsed  > 0 ) oparrk(indx,1:kmax,30) = 0.5 * ( mcons(indx,1:kmax,lsed ) + mcons(indxr,1:kmax,lsed ) )
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          else
             !
             ! no interpolation (only in case of GROUP or subgrid)
             !
             do i = 1, mxc - 1
                !
                indx  = kgrpnt(i  ,1)
                indxr = kgrpnt(i+1,1)
                !
                ! water level, water depth and inundation depth
                !
                oparr(indx,2) = s1(indxr)
                !
                oparr(indx,1) = oparr(indx,2) + dps(indxr)
                if ( ifloat /= 0 ) oparr(indx,1) = min( dps(indxr)-flos(indxr), oparr(indx,1) )
                !
                oparr(indx,16) = hindun(indxr)
                !
                ! friction velocity based on logarithmic wall-law
                !
                if ( irough == 4 .and. kmax > 1 ) then
                   !
                   oparr(indx,11) = logfrc(indxr,1)*logfrc(indxr,2)
                   !
                endif
                !
                ! non-hydrostatic pressure at bottom
                !
                if ( ihydro /= 0 ) then
                   !
                   if ( kmax == 1 ) then
                      !
                      oparr(indx,13) = q(indxr,1)
                      !
                   else
                      !
                      if ( ihydro == 3 ) then
                         !
                         oparr(indx,13) = qbot(indxr)
                         !
                      else
                         !
                         oparr(indx,13) = q(indxr,kpmax)
                         !
                      endif
                      !
                   endif
                   !
                endif
                !
                ! wave parameters
                !
                if ( lwavoutp ) then
                   !
                   oparr(indx,17) = setup(indxr)
                   oparr(indx,18) =  hsig(indxr)
                   !
                endif
                !
                oparr(indx,24) = real(brks(indxr))
                !
                ! salinity, temperature and sediment
                !
                if ( lsal  > 0 ) oparr(indx,25) = rp(indxr,1,lsal )
                if ( ltemp > 0 ) oparr(indx,26) = rp(indxr,1,ltemp)
                if ( lsed  > 0 ) oparr(indx,27) = rp(indxr,1,lsed )
                !
                ! time-averaged or mean salinity, temperature and sediment
                !
                if ( ltraoutp ) then
                   !
                   if ( lsal  > 0 ) oparr(indx,28) = mcons(indxr,1,lsal )
                   if ( ltemp > 0 ) oparr(indx,29) = mcons(indxr,1,ltemp)
                   if ( lsed  > 0 ) oparr(indx,30) = mcons(indxr,1,lsed )
                   !
                endif
                !
                ! layer-dependent scalars, if requested
                !
                if ( reqzsc ) then
                   !
                   oparrk(indx,0:kmax,11) = wom(indxr,0:kmax)
                   !
                   if ( ihydro == 0 .or. ihydro == 3 ) then
                      !
                      oparrk(indx,0:kmax,12) = wphy(indxr,0:kmax)
                      !
                   else if ( ihydro == 1 .or. ihydro == 2 ) then
                      !
                      oparrk(indx,0:kmax,12) = w1(indxr,0:kmax)
                      if ( .not.lsubg ) then
                         oparrk(indx,1:kmax,13) =  q(indxr,1:kmax)
                      else
                         oparrk(indx,1:kmax,13) =  qv(indxr,1:kmax)
                      endif
                      !
                   endif
                   !
                   if ( iturb /= 0 ) then
                      !
                      oparrk(indx,0:kmax,15) = rtur(indxr,0:kmax,1)
                      oparrk(indx,0:kmax,16) = rtur(indxr,0:kmax,2)
                      !
                   endif
                   !
                   oparrk(indx,0:kmax,17) = vnu3d(indxr,0:kmax)
                   !
                   if ( lturoutp ) then
                      !
                      oparrk(indx,0:kmax,22) = mtke(indxr,0:kmax)
                      oparrk(indx,0:kmax,23) = meps(indxr,0:kmax)
                      !
                   endif
                   !
                   if ( lsal  > 0 ) oparrk(indx,1:kmax,25) = rp(indxr,1:kmax,lsal )
                   if ( ltemp > 0 ) oparrk(indx,1:kmax,26) = rp(indxr,1:kmax,ltemp)
                   if ( lsed  > 0 ) oparrk(indx,1:kmax,27) = rp(indxr,1:kmax,lsed )
                   !
                   if ( ltraoutp ) then
                      !
                      if ( lsal  > 0 ) oparrk(indx,1:kmax,28) = mcons(indxr,1:kmax,lsal )
                      if ( ltemp > 0 ) oparrk(indx,1:kmax,29) = mcons(indxr,1:kmax,ltemp)
                      if ( lsed  > 0 ) oparrk(indx,1:kmax,30) = mcons(indxr,1:kmax,lsed )
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          endif
          !
          if ( reqvel ) then
             !
             oparr(:, 5) = oparr(:,3)
             oparr(:, 9) = oparr(:,7)
             if ( irough /= 4 .or. kmax == 1 ) oparr(:,11) = sqrt(cfricu(:)) * oparr(:,3)
             oparr(:,22) = oparr(:,20)
             !
          endif
          !
          if ( reqzvl ) then
             !
             oparrk(:,1:kmax,3) = u1(:,1:kmax)
             oparrk(:,1:kmax,5) = u1(:,1:kmax)
             oparrk(:,1:kmax,7) = u1(:,:)*hku(:,1:kmax)
             oparrk(:,1:kmax,9) = oparrk(:,1:kmax,7)
             if ( lcuroutp ) then
                oparrk(:,1:kmax,18) = mvelu(:,1:kmax)
                oparrk(:,1:kmax,20) = mvelu(:,1:kmax)
             endif
             !
          endif
          !
       else
          !
          ! compute Cartesian components of flow velocities, discharges, friction velocities and vorticity
          !
          if ( reqvel ) then
             !
             do j = 1, myc - 1
                do i = 1, mxc - 1
                   !
                   indx  = kgrpnt(i  ,j  )
                   indxr = kgrpnt(i+1,j  )
                   indxu = kgrpnt(i  ,j+1)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( indxr == 1 ) indxr = indx
                   if ( indxu == 1 ) indxu = indx
                   !
                   if ( j == 1 ) then
                      yp = 2.*ycgrid(i,1) - ycgrid(i,2)
                   else
                      yp = ycgrid(i,j-1)
                   endif
                   !
                   dydetab = ( ycgrid(i,j  ) - yp          ) / guu(indx )
                   dydetau = ( ycgrid(i,j+1) - ycgrid(i,j) ) / guu(indxu)
                   !
                   ub  = oparr(indx ,3) * dydetab
                   uu  = oparr(indxu,3) * dydetau
                   !
                   qxb = oparr(indx ,7) * dydetab
                   qxu = oparr(indxu,7) * dydetau
                   !
                   ufb = sqrt(cfricu(indx )) * ub
                   ufu = sqrt(cfricu(indxu)) * uu
                   !
                   mub = oparr(indx ,20) * dydetab
                   muu = oparr(indxu,20) * dydetau
                   !
                   if ( i == 1 ) then
                      yp = 2.*ycgrid(1,j) - ycgrid(2,j)
                   else
                      yp = ycgrid(i-1,j)
                   endif
                   !
                   dydksil = ( ycgrid(i  ,j) - yp          ) / gvv(indx )
                   dydksir = ( ycgrid(i+1,j) - ycgrid(i,j) ) / gvv(indxr)
                   !
                   vl  = oparr(indx ,4) * dydksil
                   vr  = oparr(indxr,4) * dydksir
                   !
                   qyl = oparr(indx ,8) * dydksil
                   qyr = oparr(indxr,8) * dydksir
                   !
                   vfl = sqrt(cfricv(indx )) * vl
                   vfr = sqrt(cfricv(indxr)) * vr
                   !
                   mvl = oparr(indx ,21) * dydksil
                   mvr = oparr(indxr,21) * dydksir
                   !
                   oparr(indx, 5) = 0.5 * (  uu +  ub ) - 0.5 * (  vl +  vr )
                   oparr(indx, 9) = 0.5 * ( qxu + qxb ) - 0.5 * ( qyl + qyr )
                   if ( irough /= 4 .or. kmax == 1 ) oparr(indx,11) = 0.5 * ( ufu + ufb ) - 0.5 * ( vfl + vfr )
                   oparr(indx,22) = 0.5 * ( muu + mub ) - 0.5 * ( mvl + mvr )
                   !
                   if ( reqzvl ) then
                      !
                      oparrk(indx,1:kmax,5) = 0.5 * ( u1(indxu,1:kmax)*dydetau + u1(indx ,1:kmax)*dydetab ) -  &
                                              0.5 * ( v1(indx ,1:kmax)*dydksil + v1(indxr,1:kmax)*dydksir )
                      oparrk(indx,1:kmax,9) = 0.5 * ( u1(indxu,1:kmax)*hku(indxu,1:kmax)*dydetau + u1(indx ,1:kmax)*hku(indx ,1:kmax)*dydetab ) -  &
                                              0.5 * ( v1(indx ,1:kmax)*hkv(indx ,1:kmax)*dydksil + v1(indxr,1:kmax)*hkv(indxr,1:kmax)*dydksir )
                      !
                      if ( lcuroutp ) then
                         oparrk(indx,1:kmax,20) = 0.5 * ( mvelu(indxu,1:kmax)*dydetau + mvelu(indx ,1:kmax)*dydetab ) -  &
                                                  0.5 * ( mvelv(indx ,1:kmax)*dydksil + mvelv(indxr,1:kmax)*dydksir )
                      endif
                      !
                   endif
                   !
                   if ( j == 1 ) then
                      xp = 2.*xcgrid(i,1) - xcgrid(i,2)
                   else
                      xp = xcgrid(i,j-1)
                   endif
                   !
                   dxdetab = ( xcgrid(i,j  ) - xp          ) / guu(indx )
                   dxdetau = ( xcgrid(i,j+1) - xcgrid(i,j) ) / guu(indxu)
                   !
                   ub  = oparr(indx ,3) * dxdetab
                   uu  = oparr(indxu,3) * dxdetau
                   !
                   qxb = oparr(indx ,7) * dxdetab
                   qxu = oparr(indxu,7) * dxdetau
                   !
                   ufb = sqrt(cfricu(indx )) * ub
                   ufu = sqrt(cfricu(indxu)) * uu
                   !
                   mub = oparr(indx ,20) * dxdetab
                   muu = oparr(indxu,20) * dxdetau
                   !
                   if ( i == 1 ) then
                      xp = 2.*xcgrid(1,j) - xcgrid(2,j)
                   else
                      xp = xcgrid(i-1,j)
                   endif
                   !
                   dxdksil = ( xcgrid(i  ,j) - xp          ) / gvv(indx )
                   dxdksir = ( xcgrid(i+1,j) - xcgrid(i,j) ) / gvv(indxr)
                   !
                   vl  = oparr(indx ,4) * dxdksil
                   vr  = oparr(indxr,4) * dxdksir
                   !
                   qyl = oparr(indx ,8) * dxdksil
                   qyr = oparr(indxr,8) * dxdksir
                   !
                   vfl = sqrt(cfricv(indx )) * vl
                   vfr = sqrt(cfricv(indxr)) * vr
                   !
                   mvl = oparr(indx ,21) * dxdksil
                   mvr = oparr(indxr,21) * dxdksir
                   !
                   oparr(indx, 6) = 0.5 * (  vl +  vr ) - 0.5 * (  uu +  ub )
                   oparr(indx,10) = 0.5 * ( qyl + qyr ) - 0.5 * ( qxu + qxb )
                   if ( irough /= 4 .or. kmax == 1 ) oparr(indx,12) = 0.5 * ( vfl + vfr ) - 0.5 * ( ufu + ufb )
                   oparr(indx,23) = 0.5 * ( mvl + mvr ) - 0.5 * ( muu + mub )
                   !
                   if ( reqzvl ) then
                      !
                      oparrk(indx,1:kmax, 6) = 0.5 * ( v1(indx ,1:kmax)*dxdksil + v1(indxr,1:kmax)*dxdksir ) -  &
                                               0.5 * ( u1(indxu,1:kmax)*dxdetau + u1(indx ,1:kmax)*dxdetab )
                      oparrk(indx,1:kmax,10) = 0.5 * ( v1(indx ,1:kmax)*hkv(indx ,1:kmax)*dxdksil + v1(indxr,1:kmax)*hkv(indxr,1:kmax)*dxdksir ) -  &
                                               0.5 * ( u1(indxu,1:kmax)*hku(indxu,1:kmax)*dxdetau + u1(indx ,1:kmax)*hku(indx ,1:kmax)*dxdetab )
                      !
                      if ( lcuroutp ) then
                         oparrk(indx,1:kmax,21) = 0.5 * ( mvelv(indx ,1:kmax)*dxdksil + mvelv(indxr,1:kmax)*dxdksir ) -  &
                                                  0.5 * ( mvelu(indxu,1:kmax)*dxdetau + mvelu(indx ,1:kmax)*dxdetab )
                      endif
                      !
                   endif
                   !
                   ! vorticity (see Eq. 7.56.10 of Aris (1962) with i = 3, j = 1, k = 2, i.e. even permutation)
                   !
                   dxl = 0.5 * ( gvv(indxr) + gvv(indx) )
                   dyl = 0.5 * ( guu(indxu) + guu(indx) )
                   !
                   ub  = oparr(indx ,3)
                   uu  = oparr(indxu,3)
                   !
                   vl  = oparr(indx ,4)
                   vr  = oparr(indxr,4)
                   !
                   oparr(indx,31) = ( vr - vl ) / dxl - ( uu - ub ) / dyl + 0.5 * ( (vl+vr)*( guv(indxr) - guv(indx) ) - (ub+uu)*( gvu(indxu) - gvu(indx) ) ) / ( dxl*dyl )
                   !
                enddo
             enddo
             !
          endif
          !
          if ( EQREAL(outpar(3),0.) ) then
             !
             ! water level, water depth and inundation depth
             !
             do j = 1, myc - 1
                do i = 1, mxc - 1
                   !
                   indx  = kgrpnt(i  ,j  )
                   indxr = kgrpnt(i+1,j  )
                   indxu = kgrpnt(i  ,j+1)
                   indxs = kgrpnt(i+1,j+1)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( indxr == 1 ) indxr = indx
                   if ( indxu == 1 ) indxu = indx
                   if ( indxs == 1 ) indxs = indx
                   !
                   oparr(indx,2) = 0.25 * ( s1(indx) + s1(indxr) + s1(indxu) + s1(indxs) )
                   !
                   oparr(indx,1) = oparr(indx,2) + depf(indx)
                   if ( ifloat /= 0 ) oparr(indx,1) = min( depf(indx)-flobjf(indx), oparr(indx,1) )
                   !
                   oparr(indx,16) = 0.25 * ( hindun(indx) + hindun(indxr) + hindun(indxu) + hindun(indxs) )
                   !
                enddo
             enddo
             !
             ! friction velocity based on logarithmic wall-law
             !
             if ( irough == 4 .and. kmax > 1 ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxr = kgrpnt(i+1,j  )
                      indxu = kgrpnt(i  ,j+1)
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxr == 1 ) indxr = indx
                      if ( indxu == 1 ) indxu = indx
                      if ( indxs == 1 ) indxs = indx
                      !
                      fac  = 0.25 * ( logfrc(indx ,1)*logfrc(indx ,2) + logfrc(indxr,1)*logfrc(indxr,2) + &
                                      logfrc(indxu,1)*logfrc(indxu,2) + logfrc(indxs,1)*logfrc(indxs,2) )
                      !
                      fac1 = 0.5 * ( u1(indx,kmax) + u1(indxu,kmax) )
                      fac2 = 0.5 * ( v1(indx,kmax) + v1(indxr,kmax) )
                      !
                      if ( fac1 /= 0. .or. fac2 /= 0. ) then
                         oparr(indx,11) = fac * fac1 / sqrt( fac1*fac1 + fac2*fac2 )
                         oparr(indx,12) = fac * fac2 / sqrt( fac1*fac1 + fac2*fac2 )
                      else
                         oparr(indx,11) = 0.
                         oparr(indx,12) = 0.
                      endif
                      !
                   enddo
                enddo
                !
             endif
             !
             ! non-hydrostatic pressure at bottom
             !
             if ( ihydro /= 0 ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxr = kgrpnt(i+1,j  )
                      indxu = kgrpnt(i  ,j+1)
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxr == 1 ) indxr = indx
                      if ( indxu == 1 ) indxu = indx
                      if ( indxs == 1 ) indxs = indx
                      !
                      if ( kmax == 1 ) then
                         !
                         oparr(indx,13) = 0.25 * ( q(indx,1) + q(indxr,1) + q(indxu,1) + q(indxs,1) )
                         !
                      else
                         !
                         if ( ihydro == 3 ) then
                            !
                            oparr(indx,13) = 0.25 * ( qbot(indx) + qbot(indxr) + qbot(indxu) + qbot(indxs) )
                            !
                         else
                            !
                            oparr(indx,13) = 0.25 * ( q(indx,kpmax) + q(indxr,kpmax) + q(indxu,kpmax) + q(indxs,kpmax) )
                            !
                         endif
                         !
                      endif
                      !
                   enddo
                enddo
                !
             endif
             !
             ! wave parameters
             !
             if ( lwavoutp ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxr = kgrpnt(i+1,j  )
                      indxu = kgrpnt(i  ,j+1)
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxr == 1 ) indxr = indx
                      if ( indxu == 1 ) indxu = indx
                      if ( indxs == 1 ) indxs = indx
                      !
                      oparr(indx,17) = 0.25 * ( setup(indx) + setup(indxr) + setup(indxu) + setup(indxs) )
                      oparr(indx,18) = 0.25 * (  hsig(indx) +  hsig(indxr) +  hsig(indxu) +  hsig(indxs) )
                      !
                   enddo
                enddo
                !
             endif
             !
             do j = 1, myc - 1
                do i = 1, mxc - 1
                   !
                   indx  = kgrpnt(i  ,j  )
                   indxr = kgrpnt(i+1,j  )
                   indxu = kgrpnt(i  ,j+1)
                   indxs = kgrpnt(i+1,j+1)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( indxr == 1 ) indxr = indx
                   if ( indxu == 1 ) indxu = indx
                   if ( indxs == 1 ) indxs = indx
                   !
                   oparr(indx,24) = 0.25 * ( real(brks(indx)) + real(brks(indxr)) + real(brks(indxu)) + real(brks(indxs)) )
                   !
                enddo
             enddo
             !
             ! salinity, temperature and sediment
             !
             if ( itrans /= 0 ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxr = kgrpnt(i+1,j  )
                      indxu = kgrpnt(i  ,j+1)
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxr == 1 ) indxr = indx
                      if ( indxu == 1 ) indxu = indx
                      if ( indxs == 1 ) indxs = indx
                      !
                      if ( lsal  > 0 ) oparr(indx,25) = 0.25 * ( rp(indx,1,lsal ) + rp(indxr,1,lsal ) + rp(indxu,1,lsal ) + rp(indxs,1,lsal ) )
                      if ( ltemp > 0 ) oparr(indx,26) = 0.25 * ( rp(indx,1,ltemp) + rp(indxr,1,ltemp) + rp(indxu,1,ltemp) + rp(indxs,1,ltemp) )
                      if ( lsed  > 0 ) oparr(indx,27) = 0.25 * ( rp(indx,1,lsed ) + rp(indxr,1,lsed ) + rp(indxu,1,lsed ) + rp(indxs,1,lsed ) )
                      !
                   enddo
                enddo
                !
             endif
             !
             ! time-averaged or mean salinity, temperature and sediment
             !
             if ( ltraoutp ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxr = kgrpnt(i+1,j  )
                      indxu = kgrpnt(i  ,j+1)
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxr == 1 ) indxr = indx
                      if ( indxu == 1 ) indxu = indx
                      if ( indxs == 1 ) indxs = indx
                      !
                      if ( lsal  > 0 ) oparr(indx,28) = 0.25 * ( mcons(indx,1,lsal ) + mcons(indxr,1,lsal ) + mcons(indxu,1,lsal ) + mcons(indxs,1,lsal ) )
                      if ( ltemp > 0 ) oparr(indx,29) = 0.25 * ( mcons(indx,1,ltemp) + mcons(indxr,1,ltemp) + mcons(indxu,1,ltemp) + mcons(indxs,1,ltemp) )
                      if ( lsed  > 0 ) oparr(indx,30) = 0.25 * ( mcons(indx,1,lsed ) + mcons(indxr,1,lsed ) + mcons(indxu,1,lsed ) + mcons(indxs,1,lsed ) )
                      !
                   enddo
                enddo
                !
             endif
             !
             ! layer-dependent scalars, if requested
             !
             if ( reqzsc ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxr = kgrpnt(i+1,j  )
                      indxu = kgrpnt(i  ,j+1)
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxr == 1 ) indxr = indx
                      if ( indxu == 1 ) indxu = indx
                      if ( indxs == 1 ) indxs = indx
                      !
                      oparrk(indx,0:kmax,11) = 0.25 * ( wom(indx,0:kmax) + wom(indxr,0:kmax) + wom(indxu,0:kmax) + wom(indxs,0:kmax) )
                      !
                      if ( ihydro == 0 .or. ihydro == 3 ) then
                         !
                         oparrk(indx,0:kmax,12) = 0.25 * ( wphy(indx,0:kmax) + wphy(indxr,0:kmax) + wphy(indxu,0:kmax) + wphy(indxs,0:kmax) )
                         !
                      else if ( ihydro == 1 .or. ihydro == 2 ) then
                         !
                         oparrk(indx,0:kmax,12) = 0.25 * ( w1(indx,0:kmax) + w1(indxr,0:kmax) + w1(indxu,0:kmax) + w1(indxs,0:kmax) )
                         if ( .not.lsubg ) then
                            oparrk(indx,1:kmax,13) = 0.25 * (  q(indx,1:kmax) +  q(indxr,1:kmax) +  q(indxu,1:kmax) +  q(indxs,1:kmax) )
                         else
                            oparrk(indx,1:kmax,13) = 0.25 * (  qv(indx,1:kmax) +  qv(indxr,1:kmax) +  qv(indxu,1:kmax) +  qv(indxs,1:kmax) )
                         endif
                         !
                      endif
                      !
                      if ( iturb /= 0 ) then
                         !
                         oparrk(indx,0:kmax,15) = 0.25 * ( rtur(indx,0:kmax,1) + rtur(indxr,0:kmax,1) + rtur(indxu,0:kmax,1) + rtur(indxs,0:kmax,1) )
                         oparrk(indx,0:kmax,16) = 0.25 * ( rtur(indx,0:kmax,2) + rtur(indxr,0:kmax,2) + rtur(indxu,0:kmax,2) + rtur(indxs,0:kmax,2) )
                         !
                      endif
                      !
                      oparrk(indx,0:kmax,17) = 0.25 * ( vnu3d(indx,0:kmax) + vnu3d(indxr,0:kmax) + vnu3d(indxu,0:kmax) + vnu3d(indxs,0:kmax) )
                      !
                      if ( lturoutp ) then
                         !
                         oparrk(indx,0:kmax,22) = 0.25 * ( mtke(indx,0:kmax) + mtke(indxr,0:kmax) + mtke(indxu,0:kmax) + mtke(indxs,0:kmax) )
                         oparrk(indx,0:kmax,23) = 0.25 * ( meps(indx,0:kmax) + meps(indxr,0:kmax) + meps(indxu,0:kmax) + meps(indxs,0:kmax) )
                         !
                      endif
                      !
                      if ( lsal  > 0 ) oparrk(indx,1:kmax,25) = 0.25 * ( rp(indx,1:kmax,lsal ) + rp(indxr,1:kmax,lsal ) + rp(indxu,1:kmax,lsal ) + rp(indxs,1:kmax,lsal ) )
                      if ( ltemp > 0 ) oparrk(indx,1:kmax,26) = 0.25 * ( rp(indx,1:kmax,ltemp) + rp(indxr,1:kmax,ltemp) + rp(indxu,1:kmax,ltemp) + rp(indxs,1:kmax,ltemp) )
                      if ( lsed  > 0 ) oparrk(indx,1:kmax,27) = 0.25 * ( rp(indx,1:kmax,lsed ) + rp(indxr,1:kmax,lsed ) + rp(indxu,1:kmax,lsed ) + rp(indxs,1:kmax,lsed ) )
                      !
                      if ( ltraoutp ) then
                         !
                         if ( lsal  > 0 ) oparrk(indx,1:kmax,28) = 0.25 * ( mcons(indx,1:kmax,lsal ) + mcons(indxr,1:kmax,lsal ) + mcons(indxu,1:kmax,lsal ) + mcons(indxs,1:kmax,lsal ) )
                         if ( ltemp > 0 ) oparrk(indx,1:kmax,29) = 0.25 * ( mcons(indx,1:kmax,ltemp) + mcons(indxr,1:kmax,ltemp) + mcons(indxu,1:kmax,ltemp) + mcons(indxs,1:kmax,ltemp) )
                         if ( lsed  > 0 ) oparrk(indx,1:kmax,30) = 0.25 * ( mcons(indx,1:kmax,lsed ) + mcons(indxr,1:kmax,lsed ) + mcons(indxu,1:kmax,lsed ) + mcons(indxs,1:kmax,lsed ) )
                         !
                      endif
                      !
                   enddo
                enddo
                !
             endif
             !
             ! grid-oriented velocities and discharges, if requested
             !
             if ( reqvel ) then
                !
                oparr(:,15) = oparr(:,3)
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i,j  )
                      indxu = kgrpnt(i,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxu == 1 ) indxu = indx
                      !
                      oparr(indx,3) = 0.5 * ( oparr(indx,15)          + oparr(indxu,15)           )
                      oparr(indx,7) = 0.5 * ( oparr(indx,15)*hu(indx) + oparr(indxu,15)*hu(indxu) )
                      !
                      if ( reqzvl ) then
                         !
                         oparrk(indx,1:kmax,3) = 0.5 * ( u1(indx,1:kmax)                  + u1(indxu,1:kmax)                   )
                         oparrk(indx,1:kmax,7) = 0.5 * ( u1(indx,1:kmax)*hku(indx,1:kmax) + u1(indxu,1:kmax)*hku(indxu,1:kmax) )
                         !
                      endif
                      !
                   enddo
                enddo
                !
                oparr(:,15) = oparr(:,4)
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j)
                      indxr = kgrpnt(i+1,j)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxr == 1 ) indxr = indx
                      !
                      oparr(indx,4) = 0.5 * ( oparr(indx,15)          + oparr(indxr,15)           )
                      oparr(indx,8) = 0.5 * ( oparr(indx,15)*hv(indx) + oparr(indxr,15)*hv(indxr) )
                      !
                      if ( reqzvl ) then
                         !
                         oparrk(indx,1:kmax,4) = 0.5 * ( v1(indx,1:kmax)                  + v1(indxr,1:kmax)                   )
                         oparrk(indx,1:kmax,8) = 0.5 * ( v1(indx,1:kmax)*hkv(indx,1:kmax) + v1(indxr,1:kmax)*hkv(indxr,1:kmax) )
                         !
                      endif
                      !
                   enddo
                enddo
                !
                if ( lcuroutp ) then
                   !
                   oparr(:,15) = oparr(:,20)
                   !
                   do j = 1, myc - 1
                      do i = 1, mxc - 1
                         !
                         indx  = kgrpnt(i,j  )
                         indxu = kgrpnt(i,j+1)
                         !
                         ! for permanently dry neighbours, corresponding values will be mirrored
                         !
                         if ( indxu == 1 ) indxu = indx
                         !
                         oparr(indx,20) = 0.5 * ( oparr(indx,15) + oparr(indxu,15) )
                         !
                         if ( reqzvl ) then
                            !
                            oparrk(indx,1:kmax,18) = 0.5 * ( mvelu(indx,1:kmax) + mvelu(indxu,1:kmax) )
                            !
                         endif
                         !
                      enddo
                   enddo
                   !
                   oparr(:,15) = oparr(:,21)
                   !
                   do j = 1, myc - 1
                      do i = 1, mxc - 1
                         !
                         indx  = kgrpnt(i  ,j)
                         indxr = kgrpnt(i+1,j)
                         !
                         ! for permanently dry neighbours, corresponding values will be mirrored
                         !
                         if ( indxr == 1 ) indxr = indx
                         !
                         oparr(indx,21) = 0.5 * ( oparr(indx,15) + oparr(indxr,15) )
                         !
                         if ( reqzvl ) then
                            !
                            oparrk(indx,1:kmax,19) = 0.5 * ( mvelv(indx,1:kmax) + mvelv(indxr,1:kmax) )
                            !
                         endif
                         !
                      enddo
                   enddo
                   !
                endif
                !
             endif
             !
          else
             !
             ! no interpolation (only in case of GROUP or subgrid)
             !
             ! water level, water depth and inundation depth
             !
             do j = 1, myc - 1
                do i = 1, mxc - 1
                   !
                   indx  = kgrpnt(i  ,j  )
                   indxs = kgrpnt(i+1,j+1)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( indxs == 1 ) indxs = indx
                   !
                   oparr(indx,2) = s1(indxs)
                   !
                   oparr(indx,1) = oparr(indx,2) + dps(indxs)
                   if ( ifloat /= 0 ) oparr(indx,1) = min( dps(indxs)-flos(indxs), oparr(indx,1) )
                   !
                   oparr(indx,16) = hindun(indxs)
                   !
                enddo
             enddo
             !
             ! friction velocity based on logarithmic wall-law
             !
             if ( irough == 4 .and. kmax > 1 ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxr = kgrpnt(i+1,j  )
                      indxu = kgrpnt(i  ,j+1)
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxr == 1 ) indxr = indx
                      if ( indxu == 1 ) indxu = indx
                      if ( indxs == 1 ) indxs = indx
                      !
                      fac  = logfrc(indxs,1)*logfrc(indxs,2)
                      !
                      fac1 = 0.5 * ( u1(indxs,kmax) + u1(indxu,kmax) )
                      fac2 = 0.5 * ( v1(indxs,kmax) + v1(indxr,kmax) )
                      !
                      if ( fac1 /= 0. .or. fac2 /= 0. ) then
                         oparr(indx,11) = fac * fac1 / sqrt( fac1*fac1 + fac2*fac2 )
                         oparr(indx,12) = fac * fac2 / sqrt( fac1*fac1 + fac2*fac2 )
                      else
                         oparr(indx,11) = 0.
                         oparr(indx,12) = 0.
                      endif
                      !
                   enddo
                enddo
                !
             endif
             !
             ! non-hydrostatic pressure at bottom
             !
             if ( ihydro /= 0 ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxs == 1 ) indxs = indx
                      !
                      if ( kmax == 1 ) then
                         !
                         oparr(indx,13) = q(indxs,1)
                         !
                      else
                         !
                         if ( ihydro == 3 ) then
                            !
                            oparr(indx,13) = qbot(indxs)
                            !
                         else
                            !
                            oparr(indx,13) = q(indxs,kpmax)
                            !
                         endif
                         !
                      endif
                      !
                   enddo
                enddo
                !
             endif
             !
             ! wave parameters
             !
             if ( lwavoutp ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxs == 1 ) indxs = indx
                      !
                      oparr(indx,17) = setup(indxs)
                      oparr(indx,18) =  hsig(indxs)
                      !
                   enddo
                enddo
                !
             endif
             !
             do j = 1, myc - 1
                do i = 1, mxc - 1
                   !
                   indx  = kgrpnt(i  ,j  )
                   indxs = kgrpnt(i+1,j+1)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( indxs == 1 ) indxs = indx
                   !
                   oparr(indx,24) = real(brks(indxs))
                   !
                enddo
             enddo
             !
             ! salinity, temperature and sediment
             !
             if ( itrans /= 0 ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxs == 1 ) indxs = indx
                      !
                      if ( lsal  > 0 ) oparr(indx,25) = rp(indxs,1,lsal )
                      if ( ltemp > 0 ) oparr(indx,26) = rp(indxs,1,ltemp)
                      if ( lsed  > 0 ) oparr(indx,27) = rp(indxs,1,lsed )
                      !
                   enddo
                enddo
                !
             endif
             !
             ! time-averaged or mean salinity, temperature and sediment
             !
             if ( ltraoutp ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxs == 1 ) indxs = indx
                      !
                      if ( lsal  > 0 ) oparr(indx,28) = mcons(indxs,1,lsal )
                      if ( ltemp > 0 ) oparr(indx,29) = mcons(indxs,1,ltemp)
                      if ( lsed  > 0 ) oparr(indx,30) = mcons(indxs,1,lsed )
                      !
                   enddo
                enddo
                !
             endif
             !
             ! layer-dependent scalars, if requested
             !
             if ( reqzsc ) then
                !
                do j = 1, myc - 1
                   do i = 1, mxc - 1
                      !
                      indx  = kgrpnt(i  ,j  )
                      indxs = kgrpnt(i+1,j+1)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( indxs == 1 ) indxs = indx
                      !
                      oparrk(indx,0:kmax,11) = wom(indxs,0:kmax)
                      !
                      if ( ihydro == 0 .or. ihydro == 3 ) then
                         !
                         oparrk(indx,0:kmax,12) = wphy(indxs,0:kmax)
                         !
                      else if ( ihydro == 1 .or. ihydro == 2 ) then
                         !
                         oparrk(indx,0:kmax,12) = w1(indxs,0:kmax)
                         if ( .not.lsubg ) then
                            oparrk(indx,1:kmax,13) =  q(indxs,1:kmax)
                         else
                            oparrk(indx,1:kmax,13) =  qv(indxs,1:kmax)
                         endif
                         !
                      endif
                      !
                      if ( iturb /= 0 ) then
                         !
                         oparrk(indx,0:kmax,15) = rtur(indxs,0:kmax,1)
                         oparrk(indx,0:kmax,16) = rtur(indxs,0:kmax,2)
                         !
                      endif
                      !
                      oparrk(indx,0:kmax,17) = vnu3d(indxs,0:kmax)
                      !
                      if ( lturoutp ) then
                         !
                         oparrk(indx,0:kmax,22) = mtke(indxs,0:kmax)
                         oparrk(indx,0:kmax,23) = meps(indxs,0:kmax)
                         !
                      endif
                      !
                      if ( lsal  > 0 ) oparrk(indx,1:kmax,25) = rp(indxs,1:kmax,lsal )
                      if ( ltemp > 0 ) oparrk(indx,1:kmax,26) = rp(indxs,1:kmax,ltemp)
                      if ( lsed  > 0 ) oparrk(indx,1:kmax,27) = rp(indxs,1:kmax,lsed )
                      !
                      if ( ltraoutp ) then
                         !
                         if ( lsal  > 0 ) oparrk(indx,1:kmax,28) = mcons(indxs,1:kmax,lsal )
                         if ( ltemp > 0 ) oparrk(indx,1:kmax,29) = mcons(indxs,1:kmax,ltemp)
                         if ( lsed  > 0 ) oparrk(indx,1:kmax,30) = mcons(indxs,1:kmax,lsed )
                         !
                      endif
                      !
                   enddo
                enddo
                !
             endif
             !
             ! grid-oriented layer-dependent velocities and discharges, if requested
             !
             if ( reqzvl ) then
                !
                oparrk(:,1:kmax,3) = u1(:,1:kmax)
                oparrk(:,1:kmax,7) = u1(:,1:kmax)*hku(:,1:kmax)
                oparrk(:,1:kmax,4) = v1(:,1:kmax)
                oparrk(:,1:kmax,8) = v1(:,1:kmax)*hkv(:,1:kmax)
                if ( lcuroutp ) then
                   oparrk(:,1:kmax,18) = mvelu(:,1:kmax)
                   oparrk(:,1:kmax,19) = mvelv(:,1:kmax)
                endif
                !
             endif
             !
          endif
          !
       endif
       !
    else
       !
       if ( irough == 4 .and. kmax > 1 ) then
          !
          ! compute cell-based velocity vector
          !
          call perot ( u1, 1, kmax )
          !
       endif
       !
       if ( EQREAL(outpar(3),0.) ) then
          !
          do ivert = 1, nverts
             !
             totarea = 0.
             !
             fac  = 0.
             fac1 = 0.
             fac2 = 0.
             !
             do jc = 1, vert(ivert)%noc
                !
                icell = vert(ivert)%cell(jc)%atti(CELLID)
                !
                area = cell(icell)%attr(CELLAREA)
                !
                oparr(ivert, 2) = oparr(ivert, 2) + area *     s1(icell)
                oparr(ivert,16) = oparr(ivert,16) + area * hindun(icell)
                !
                if ( irough == 4 .and. kmax > 1 ) then
                   !
                   fac  = fac  + area * logfrc(icell,1)*logfrc(icell,2)
                   fac1 = fac1 + area * uvc(icell,kmax,1)
                   fac2 = fac2 + area * uvc(icell,kmax,2)
                   !
                endif
                !
                if ( ihydro /= 0 ) then
                   !
                   oparr(ivert,13) = oparr(ivert,13) + area * q(icell,kmax)
                   !
                endif
                !
                if ( lwavoutp ) then
                   !
                   oparr(ivert,17) = oparr(ivert,17) + area * setup(icell)
                   oparr(ivert,18) = oparr(ivert,18) + area *  hsig(icell)
                   !
                endif
                !
                oparr(ivert,24) = oparr(ivert,24) + area * real(brks(icell))
                !
                ! salinity, temperature and sediment
                !
                if ( itrans /= 0 ) then
                   !
                   if ( lsal  > 0 ) oparr(ivert,25) = oparr(ivert,25) + area * rp(icell,1,lsal )
                   if ( ltemp > 0 ) oparr(ivert,26) = oparr(ivert,26) + area * rp(icell,1,ltemp)
                   if ( lsed  > 0 ) oparr(ivert,27) = oparr(ivert,27) + area * rp(icell,1,lsed )
                   !
                endif
                !
                ! time-averaged or mean salinity, temperature and sediment
                !
                if ( ltraoutp ) then
                   !
                   if ( lsal  > 0 ) oparr(ivert,28) = oparr(ivert,28) + area * mcons(icell,1,lsal )
                   if ( ltemp > 0 ) oparr(ivert,29) = oparr(ivert,29) + area * mcons(icell,1,ltemp)
                   if ( lsed  > 0 ) oparr(ivert,30) = oparr(ivert,30) + area * mcons(icell,1,lsed )
                   !
                endif
                !
                ! layer-dependent scalars, if requested
                !
                if ( reqzsc ) then
                   !
                   oparrk(ivert,0:kmax,11) = oparrk(ivert,0:kmax,11) + area * wom(icell,0:kmax)
                   !
                   if ( ihydro == 0 .or. ihydro == 3 ) then
                      !
                      oparrk(ivert,0:kmax,12) = oparrk(ivert,0:kmax,12) + area * wphy(icell,0:kmax)
                      !
                   else if ( ihydro == 1 .or. ihydro == 2 ) then
                      !
                      oparrk(ivert,0:kmax,12) = oparrk(ivert,0:kmax,12) + area * w1(icell,0:kmax)
                      oparrk(ivert,1:kmax,13) = oparrk(ivert,1:kmax,13) + area *  q(icell,1:kmax)
                      !
                   endif
                   !
                   if ( iturb /= 0 ) then
                      !
                      oparrk(ivert,0:kmax,15) = oparrk(ivert,0:kmax,15) + area * rtur(icell,0:kmax,1)
                      oparrk(ivert,0:kmax,16) = oparrk(ivert,0:kmax,16) + area * rtur(icell,0:kmax,2)
                      !
                   endif
                   !
                   oparrk(ivert,0:kmax,17) = oparrk(ivert,0:kmax,17) + area * vnu3d(icell,0:kmax)
                   !
                   if ( lturoutp ) then
                      !
                      oparrk(ivert,0:kmax,22) = oparrk(ivert,0:kmax,22) + area * mtke(icell,0:kmax)
                      oparrk(ivert,0:kmax,23) = oparrk(ivert,0:kmax,23) + area * meps(icell,0:kmax)
                      !
                   endif
                   !
                   if ( itrans /= 0 ) then
                      !
                      if ( lsal  > 0 ) oparrk(ivert,1:kmax,25) = oparrk(ivert,1:kmax,25) + area * rp(icell,1:kmax,lsal )
                      if ( ltemp > 0 ) oparrk(ivert,1:kmax,26) = oparrk(ivert,1:kmax,26) + area * rp(icell,1:kmax,ltemp)
                      if ( lsed  > 0 ) oparrk(ivert,1:kmax,27) = oparrk(ivert,1:kmax,27) + area * rp(icell,1:kmax,lsed )
                      !
                   endif
                   !
                   if ( ltraoutp ) then
                      !
                      if ( lsal  > 0 ) oparrk(ivert,1:kmax,28) = oparrk(ivert,1:kmax,28) + area * mcons(icell,1:kmax,lsal )
                      if ( ltemp > 0 ) oparrk(ivert,1:kmax,29) = oparrk(ivert,1:kmax,29) + area * mcons(icell,1:kmax,ltemp)
                      if ( lsed  > 0 ) oparrk(ivert,1:kmax,30) = oparrk(ivert,1:kmax,30) + area * mcons(icell,1:kmax,lsed )
                      !
                   endif
                   !
                endif
                !
                totarea = totarea + area
                !
             enddo
             !
             ! water level, water depth and inundation depth
             !
             oparr(ivert, 2) = oparr(ivert, 2) / totarea
             oparr(ivert, 1) = oparr(ivert, 2) + depf(ivert)
             oparr(ivert,16) = oparr(ivert,16) / totarea
             !
             ! friction velocity based on logarithmic wall-law
             !
             if ( irough == 4 .and. kmax > 1 ) then
                !
                fac  = fac  / totarea
                fac1 = fac1 / totarea
                fac2 = fac2 / totarea
                !
                if ( fac1 /= 0. .or. fac2 /= 0. ) then
                   oparr(ivert,11) = fac * fac1 / sqrt( fac1*fac1 + fac2*fac2 )
                   oparr(ivert,12) = fac * fac2 / sqrt( fac1*fac1 + fac2*fac2 )
                else
                   oparr(ivert,11) = 0.
                   oparr(ivert,12) = 0.
                endif
                !
             endif
             !
             ! non-hydrostatic pressure at bottom
             !
             oparr(ivert,13) = oparr(ivert,13) / totarea
             !
             ! wave parameters
             !
             oparr(ivert,17) = oparr(ivert,17) / totarea
             oparr(ivert,18) = oparr(ivert,18) / totarea
             !
             oparr(ivert,24) = oparr(ivert,24) / totarea
             !
             ! constituents
             !
             oparr(ivert,25) = oparr(ivert,25) / totarea
             oparr(ivert,26) = oparr(ivert,26) / totarea
             oparr(ivert,27) = oparr(ivert,27) / totarea
             oparr(ivert,28) = oparr(ivert,28) / totarea
             oparr(ivert,29) = oparr(ivert,29) / totarea
             oparr(ivert,30) = oparr(ivert,30) / totarea
             !
             ! layer-dependent scalars, if requested
             !
             if ( reqzsc ) then
                !
                oparrk(ivert,0:kmax,11) = oparrk(ivert,0:kmax,11) / totarea
                oparrk(ivert,0:kmax,12) = oparrk(ivert,0:kmax,12) / totarea
                oparrk(ivert,1:kmax,13) = oparrk(ivert,1:kmax,13) / totarea
                !
                oparrk(ivert,0:kmax,15) = oparrk(ivert,0:kmax,15) / totarea
                oparrk(ivert,0:kmax,16) = oparrk(ivert,0:kmax,16) / totarea
                oparrk(ivert,0:kmax,17) = oparrk(ivert,0:kmax,17) / totarea
                oparrk(ivert,0:kmax,22) = oparrk(ivert,0:kmax,22) / totarea
                oparrk(ivert,0:kmax,23) = oparrk(ivert,0:kmax,23) / totarea
                !
                oparrk(ivert,1:kmax,25) = oparrk(ivert,1:kmax,25) / totarea
                oparrk(ivert,1:kmax,26) = oparrk(ivert,1:kmax,26) / totarea
                oparrk(ivert,1:kmax,27) = oparrk(ivert,1:kmax,27) / totarea
                oparrk(ivert,1:kmax,28) = oparrk(ivert,1:kmax,28) / totarea
                oparrk(ivert,1:kmax,29) = oparrk(ivert,1:kmax,29) / totarea
                oparrk(ivert,1:kmax,30) = oparrk(ivert,1:kmax,30) / totarea
                !
             endif
             !
          enddo
          !
       endif
       !
       if ( reqvel ) then
          !
          ! compute the depth-averaged flow velocity
          !
          if ( kmax == 1 ) then
             !
             uwrk(:) = u1(:,1)
             !
          else
             !
             uwrk = 0.
             !
             do iface = 1, nfaces
                !
                if ( hu(iface) > 0. ) then
                   !
                   do k = 1, kmax
                      !
                      uwrk(iface) = uwrk(iface) + hku(iface,k) * u1(iface,k)
                      !
                   enddo
                   !
                   uwrk(iface) = uwrk(iface) / hu(iface)
                   !
                endif
                !
             enddo
             !
          endif
          !
          ! compute cell-based velocity vector
          !
          call perot ( uwrk, 1, 1 )
          !
          if ( EQREAL(outpar(3),0.) ) then
             !
             do ivert = 1, nverts
                !
                totarea = 0.
                !
                do jc = 1, vert(ivert)%noc
                   !
                   icell = vert(ivert)%cell(jc)%atti(CELLID)
                   !
                   area = cell(icell)%attr(CELLAREA)
                   !
                   oparr(ivert, 5) = oparr(ivert, 5) + area * uvc(icell,1,1)
                   oparr(ivert, 6) = oparr(ivert, 6) + area * uvc(icell,1,2)
                   oparr(ivert, 9) = oparr(ivert, 9) + area * uvc(icell,1,1) * hs(icell)
                   oparr(ivert,10) = oparr(ivert,10) + area * uvc(icell,1,2) * hs(icell)
                   !
                   totarea = totarea + area
                   !
                enddo
                !
                oparr(ivert, 5) = oparr(ivert, 5) / totarea
                oparr(ivert, 6) = oparr(ivert, 6) / totarea
                oparr(ivert, 9) = oparr(ivert, 9) / totarea
                oparr(ivert,10) = oparr(ivert,10) / totarea
                !
             enddo
             !
          endif
          !
          ! compute friction velocity
          !
          if ( irough /= 4 .or. kmax == 1 ) then
             !
             if ( EQREAL(outpar(3),0.) ) then
                !
                do ivert = 1, nverts
                   !
                   totarea = 0.
                   !
                   do jc = 1, vert(ivert)%noc
                      !
                      icell = vert(ivert)%cell(jc)%atti(CELLID)
                      !
                      area = cell(icell)%attr(CELLAREA)
                      !
                      oparr(ivert,11) = oparr(ivert,11) + area * sqrt(cfricu(icell)) * uvc(icell,1,1)
                      oparr(ivert,12) = oparr(ivert,12) + area * sqrt(cfricu(icell)) * uvc(icell,1,2)
                      !
                      totarea = totarea + area
                      !
                   enddo
                   !
                   oparr(ivert,11) = oparr(ivert,11) / totarea
                   oparr(ivert,12) = oparr(ivert,12) / totarea
                   !
                enddo
                !
             endif
             !
          endif
          !
          ! compute depth-averaged mean velocity
          !
          if ( lcuroutp ) then
             !
             if ( kmax == 1 ) then
                !
                uwrk(:) = mvelu(:,1)
                !
             else
                !
                uwrk = 0.
                !
                do iface = 1, nfaces
                   !
                   if ( hu(iface) > 0. ) then
                      !
                      do k = 1, kmax
                         !
                         uwrk(iface) = uwrk(iface) + hku(iface,k) * mvelu(iface,k)
                         !
                      enddo
                      !
                      uwrk(iface) = uwrk(iface) / hu(iface)
                      !
                   endif
                   !
                enddo
                !
             endif
             !
             ! compute cell-based velocity vector
             !
             call perot ( uwrk, 1, 1 )
             !
             if ( EQREAL(outpar(3),0.) ) then
                !
                do ivert = 1, nverts
                   !
                   totarea = 0.
                   !
                   do jc = 1, vert(ivert)%noc
                      !
                      icell = vert(ivert)%cell(jc)%atti(CELLID)
                      !
                      area = cell(icell)%attr(CELLAREA)
                      !
                      oparr(ivert,22) = oparr(ivert,22) + area * uvc(icell,1,1)
                      oparr(ivert,23) = oparr(ivert,23) + area * uvc(icell,1,2)
                      !
                      totarea = totarea + area
                      !
                   enddo
                   !
                   oparr(ivert,22) = oparr(ivert,22) / totarea
                   oparr(ivert,23) = oparr(ivert,23) / totarea
                   !
                enddo
                !
             endif
             !
          endif
          !
       endif
       !
       if ( reqzvl ) then
          !
          ! compute cell-based velocity vector
          !
          call perot ( u1, 1, kmax )
          !
          if ( EQREAL(outpar(3),0.) ) then
             !
             do ivert = 1, nverts
                !
                totarea = 0.
                !
                do jc = 1, vert(ivert)%noc
                   !
                   icell = vert(ivert)%cell(jc)%atti(CELLID)
                   !
                   area = cell(icell)%attr(CELLAREA)
                   !
                   oparrk(ivert,1:kmax, 5) = oparrk(ivert,1:kmax, 5) + area * uvc(icell,1:kmax,1)
                   oparrk(ivert,1:kmax, 6) = oparrk(ivert,1:kmax, 6) + area * uvc(icell,1:kmax,2)
                   oparrk(ivert,1:kmax, 9) = oparrk(ivert,1:kmax, 9) + area * uvc(icell,1:kmax,1) * hks(icell,1:kmax)
                   oparrk(ivert,1:kmax,10) = oparrk(ivert,1:kmax,10) + area * uvc(icell,1:kmax,2) * hks(icell,1:kmax)
                   !
                   totarea = totarea + area
                   !
                enddo
                !
                oparrk(ivert,1:kmax, 5) = oparrk(ivert,1:kmax, 5) / totarea
                oparrk(ivert,1:kmax, 6) = oparrk(ivert,1:kmax, 6) / totarea
                oparrk(ivert,1:kmax, 9) = oparrk(ivert,1:kmax, 9) / totarea
                oparrk(ivert,1:kmax,10) = oparrk(ivert,1:kmax,10) / totarea
                !
             enddo
             !
          endif
          !
          if ( lcuroutp ) then
             !
             ! compute cell-based mean velocity vector
             !
             call perot ( mvelu, 1, kmax )
             !
             if ( EQREAL(outpar(3),0.) ) then
                !
                do ivert = 1, nverts
                   !
                   totarea = 0.
                   !
                   do jc = 1, vert(ivert)%noc
                      !
                      icell = vert(ivert)%cell(jc)%atti(CELLID)
                      !
                      area = cell(icell)%attr(CELLAREA)
                      !
                      oparrk(ivert,1:kmax,20) = oparrk(ivert,1:kmax,20) + area * uvc(icell,1:kmax,1)
                      oparrk(ivert,1:kmax,21) = oparrk(ivert,1:kmax,21) + area * uvc(icell,1:kmax,2)
                      !
                      totarea = totarea + area
                      !
                   enddo
                   !
                   oparrk(ivert,1:kmax,20) = oparrk(ivert,1:kmax,20) / totarea
                   oparrk(ivert,1:kmax,21) = oparrk(ivert,1:kmax,21) / totarea
                   !
                enddo
                !
             endif
             !
          endif
          !
       endif
       !
    endif
    !
    ! determine layer interfaces and layer thicknesses, if requested
    !
    if ( oqproc(51) .or. oqproc(71) .or. oqproc(83) ) then
       !
       oparrk(:,   0,2) = oparr(:,2)
       oparrk(:,kmax,2) = oparr(:,2) - oparr(:,1)
       !
       call sigmacoor ( oparrk(1,0,2), oparr(1,1), mcgrd )
       !
       do k = 1, kmax
          !
          oparrk(:,k,1) = oparrk(:,k-1,2) - oparrk(:,k,2)
          !
       enddo
       !
    endif
    !
    ! compute pressure at bottom
    !
    if ( oqproc(13) ) then
       !
       if ( ihydro == 0 ) then
          !
          oparr(:,14) = 0.01 * rhow * grav * oparr(:,1)
          !
       else
          !
          oparr(:,14) = 0.01 * rhow * ( grav * oparr(:,1) + oparr(:,13) )
          !
       endif
       !
    endif
    !
    ! compute pressure
    !
    if ( oqproc(83) ) then
       !
       if ( ihydro == 0 ) then
          !
          do k = 1, kmax
             !
             oparrk(:,k,14) = 0.01 * rhow * grav * ( oparr(:,2) - oparrk(:,k,2) )
             !
          enddo
          !
       else if ( ihydro /= 3 ) then
          !
          do k = 1, kmax
             !
             oparrk(:,k,14) = 0.01 * rhow * ( grav * ( oparr(:,2) - oparrk(:,k,2) ) + oparrk(:,k,13) )
             !
          enddo
          !
       endif
       !
    endif
    !
    ! compute RMS wave height
    !
    if ( oqproc(24) ) then
       !
       oparr(:,19) = 0.5 * sqrt(2.) * oparr(:,18)
       !
    endif
    !
    ! compute time-averaged vertical eddy viscosity
    !
    if ( oqproc(59) ) then
       !
       if ( iturb == 0 ) then
          oparrk(:,:,24) = vnu3d(:,:)
       else
          oparrk(:,:,24) = 0.09 * oparrk(:,:,22) * oparrk(:,:,22) / ( oparrk(:,:,23) + 1.e-8 )
       endif
       !
    endif
    !
    ! compute suspended sediment
    !
    if ( oqproc(39) ) then
       !
       oparr(:,27) = rhos * oparr(:,27)
       !
    endif
    !
    if ( oqproc(44) ) then
       !
       oparr(:,30) = rhos * oparr(:,30)
       !
    endif
    !
    if ( oqproc(91) ) then
       !
       oparrk(:,:,27) = rhos * oparrk(:,:,27)
       !
    endif
    !
    if ( oqproc(94) ) then
       !
       oparrk(:,:,30) = rhos * oparrk(:,:,30)
       !
    endif
    !
    ! find closest vertex for given point in case of unstructured grid
    !
    if ( optg == 5 ) then
       !
       allocate(kvert(mip))
       if ( .not.lcompgrd ) then
          do i = 1, mip
             call SwanFindPoint ( voq(i,1), voq(i,2), kvert(i) )
          enddo
       else
          do i = 1, mip
             kvert(i) = i
          enddo
       endif
       !
    endif
    !
    ! Tsec
    !
    if ( oqproc(41) ) voq(:,voqr(41)) = real(timco) - outpar(1)
    !
    ! distance
    !
    if ( oqproc(3) ) then
       do i = 1, mip
          if ( i == 1 ) then
             rdist = 0.
          else
             rdx = voq(i,1) - voq(i-1,1)
             rdy = voq(i,2) - voq(i-1,2)
             if ( kspher > 0 ) then
                rdx = rdx * lendeg * cos(degrad*(yoffs + 0.5*(voq(i,2) + voq(i-1,2))))
                rdy = rdy * lendeg
             endif
             rdist = rdist + sqrt( rdx*rdx + rdy*rdy )
          endif
          voq(i,voqr(3)) = rdist
       enddo
    endif
    !
    ! water depth
    !
    if ( oqproc(4) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,1), ovexcv(4), xc, yc, mip, voq(1,voqr(4)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(4)), voq(1,1), voq(1,2), oparr(1,1), mip, kvert, ovexcv(4) )
       endif
    endif
    !
    ! bottom level
    !
    if ( oqproc(5) ) then
       if ( optg /= 5 ) then
          ! interpolation done in all active and non-active points
          rval1  = epsdry
          epsdry = -99999.
          call SWIPOL ( depf, ovexcv(5), xc, yc, mip, voq(1,voqr(5)), kgrpnt, oparr(1,1) )
          epsdry = rval1
       else
          ! interpolation done in all active and non-active points
          allocate(ltmp(nverts))
          ltmp(:)        = vert(:)%active
          vert(:)%active = .true.
          call SwanInterpolateOutput ( voq(1,voqr(5)), voq(1,1), voq(1,2), depf, mip, kvert, ovexcv(5) )
          vert(:)%active = ltmp(:)
          deallocate(ltmp)
       endif
    endif
    !
    ! water level
    !
    if ( oqproc(6) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,2), ovexcv(6), xc, yc, mip, voq(1,voqr(6)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(6)), voq(1,1), voq(1,2), oparr(1,2), mip, kvert, ovexcv(6) )
       endif
    endif
    !
    ! flow velocities in x- and y-directions
    !
    if ( oqproc(9) ) then
       jvx = voqr(9)
       jvy = jvx + 1
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,5), ovexcv(9), xc, yc, mip, voq(1,jvx), kgrpnt, oparr(1,1) )
          call SWIPOL ( oparr(1,6), ovexcv(9), xc, yc, mip, voq(1,jvy), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,jvx), voq(1,1), voq(1,2), oparr(1,5), mip, kvert, ovexcv(9) )
          call SwanInterpolateOutput ( voq(1,jvy), voq(1,1), voq(1,2), oparr(1,6), mip, kvert, ovexcv(9) )
       endif
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,jvx) = coscq*rval1 - sincq*rval2
          voq(i,jvy) = sincq*rval1 + coscq*rval2
       enddo
    endif
    !
    ! velocity magnitude
    !
    if ( oqproc(7) ) then
       jvx = voqr(9)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,voqr(7)) = sqrt(rval1*rval1 + rval2*rval2)
       enddo
    endif
    !
    ! velocity direction
    !
    if ( oqproc(8) ) then
       jvx = voqr(9)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          if ( rval1 /= 0. .or. rval2 /= 0. ) then
             if ( bnaut ) then
                dirdeg = atan2(rval2,rval1) * 180./pi
             else
                dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
             endif
             if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
             !
             ! convert (if necessary) from nautical degrees to Cartesian degrees
             !
             voq(i,voqr(8)) = DEGCNV( dirdeg )
             !
          else
             voq(i,voqr(8)) = ovexcv(8)
          endif
       enddo
    endif
    !
    ! grid-oriented U-velocity
    !
    if ( oqproc(10) ) call SWIPOL ( oparr(1,3), ovexcv(10), xc, yc, mip, voq(1,voqr(10)), kgrpnt, oparr(1,1) )
    !
    ! grid-oriented V-velocity
    !
    if ( oqproc(11) ) then
       if ( .not.oned ) then
          call SWIPOL ( oparr(1,4), ovexcv(11), xc, yc, mip, voq(1,voqr(11)), kgrpnt, oparr(1,1) )
       else
          voq(:,voqr(11)) = ovexcv(11)
       endif
    endif
    !
    ! non-hydrostatic pressure at bottom
    !
    if ( oqproc(12) ) then
       if ( ihydro /= 0 ) then
          if ( optg /= 5 ) then
             call SWIPOL ( oparr(1,13), ovexcv(12), xc, yc, mip, voq(1,voqr(12)), kgrpnt, oparr(1,1) )
          else
             call SwanInterpolateOutput ( voq(1,voqr(12)), voq(1,1), voq(1,2), oparr(1,13), mip, kvert, ovexcv(12) )
          endif
       else
          voq(:,voqr(12)) = ovexcv(12)
       endif
    endif
    !
    ! pressure at bottom
    !
    if ( oqproc(13) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,14), ovexcv(13), xc, yc, mip, voq(1,voqr(13)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(13)), voq(1,1), voq(1,2), oparr(1,14), mip, kvert, ovexcv(13) )
       endif
    endif
    !
    ! specific discharges in x- and y-directions
    !
    if ( oqproc(16) ) then
       jvx = voqr(16)
       jvy = jvx + 1
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1, 9), ovexcv(16), xc, yc, mip, voq(1,jvx), kgrpnt, oparr(1,1) )
          call SWIPOL ( oparr(1,10), ovexcv(16), xc, yc, mip, voq(1,jvy), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,jvx), voq(1,1), voq(1,2), oparr(1, 9), mip, kvert, ovexcv(16) )
          call SwanInterpolateOutput ( voq(1,jvy), voq(1,1), voq(1,2), oparr(1,10), mip, kvert, ovexcv(16) )
       endif
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,jvx) = coscq*rval1 - sincq*rval2
          voq(i,jvy) = sincq*rval1 + coscq*rval2
       enddo
    endif
    !
    ! magnitude of specific discharge
    !
    if ( oqproc(14) ) then
       jvx = voqr(16)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,voqr(14)) = sqrt(rval1*rval1 + rval2*rval2)
       enddo
    endif
    !
    ! direction of specific discharge
    !
    if ( oqproc(15) ) then
       jvx = voqr(16)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          if ( rval1 /= 0. .or. rval2 /= 0. ) then
             if ( bnaut ) then
                dirdeg = atan2(rval2,rval1) * 180./pi
             else
                dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
             endif
             if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
             !
             ! convert (if necessary) from nautical degrees to Cartesian degrees
             !
             voq(i,voqr(15)) = DEGCNV( dirdeg )
             !
          else
             voq(i,voqr(15)) = ovexcv(15)
          endif
       enddo
    endif
    !
    ! grid-oriented U-discharge per unit width
    !
    if ( oqproc(17) ) call SWIPOL ( oparr(1,7), ovexcv(17), xc, yc, mip, voq(1,voqr(17)), kgrpnt, oparr(1,1) )
    !
    ! grid-oriented V-discharge per unit width
    !
    if ( oqproc(18) ) then
       if ( .not.oned ) then
          call SWIPOL ( oparr(1,8), ovexcv(18), xc, yc, mip, voq(1,voqr(18)), kgrpnt, oparr(1,1) )
       else
          voq(:,voqr(18)) = ovexcv(18)
       endif
    endif
    !
    ! salinity
    !
    if ( oqproc(19) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,25), ovexcv(19), xc, yc, mip, voq(1,voqr(19)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(19)), voq(1,1), voq(1,2), oparr(1,25), mip, kvert, ovexcv(19) )
       endif
    endif
    !
    ! temperature
    !
    if ( oqproc(20) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,26), ovexcv(20), xc, yc, mip, voq(1,voqr(20)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(20)), voq(1,1), voq(1,2), oparr(1,26), mip, kvert, ovexcv(20) )
       endif
    endif
    !
    ! suspended sediment
    !
    if ( oqproc(39) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,27), ovexcv(39), xc, yc, mip, voq(1,voqr(39)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(39)), voq(1,1), voq(1,2), oparr(1,27), mip, kvert, ovexcv(39) )
       endif
    endif
    !
    ! inundation depth / maximum horizontal runup
    !
    if ( oqproc(21) ) then
       if ( optg /= 5 ) then
          ! interpolation done in all active and non-active points
          rval1  = epsdry
          epsdry = -99999.
          call SWIPOL ( oparr(1,16), ovexcv(21), xc, yc, mip, voq(1,voqr(21)), kgrpnt, oparr(1,1) )
          epsdry = rval1
       else
          ! interpolation done in all active and non-active points
          allocate(ltmp(nverts))
          ltmp(:)        = vert(:)%active
          vert(:)%active = .true.
          call SwanInterpolateOutput ( voq(1,voqr(21)), voq(1,1), voq(1,2), oparr(1,16), mip, kvert, ovexcv(21) )
          vert(:)%active = ltmp(:)
          deallocate(ltmp)
       endif
    endif
    !
    ! wave breaking points
    !
    if ( oqproc(38) ) then
       if ( optg /= 5 ) then
          ! interpolation done in all active and non-active points
          rval1  = epsdry
          epsdry = -99999.
          call SWIPOL ( oparr(1,24), ovexcv(38), xc, yc, mip, voq(1,voqr(38)), kgrpnt, oparr(1,1) )
          epsdry = rval1
       else
          ! interpolation done in all active and non-active points
          allocate(ltmp(nverts))
          ltmp(:)        = vert(:)%active
          vert(:)%active = .true.
          call SwanInterpolateOutput ( voq(1,voqr(38)), voq(1,1), voq(1,2), oparr(1,24), mip, kvert, ovexcv(38) )
          vert(:)%active = ltmp(:)
          deallocate(ltmp)
       endif
    endif
    !
    ! wave-induced setup
    !
    if ( oqproc(22) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,17), ovexcv(22), xc, yc, mip, voq(1,voqr(22)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(22)), voq(1,1), voq(1,2), oparr(1,17), mip, kvert, ovexcv(22) )
       endif
    endif
    !
    ! significant wave height
    !
    if ( oqproc(23) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,18), ovexcv(23), xc, yc, mip, voq(1,voqr(23)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(23)), voq(1,1), voq(1,2), oparr(1,18), mip, kvert, ovexcv(23) )
       endif
    endif
    !
    ! RMS wave height
    !
    if ( oqproc(24) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,19), ovexcv(24), xc, yc, mip, voq(1,voqr(24)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(24)), voq(1,1), voq(1,2), oparr(1,19), mip, kvert, ovexcv(24) )
       endif
    endif
    !
    ! wind velocity
    !
    if ( oqproc(29) ) then
       jvx = voqr(29)
       jvy = jvx + 1
       if ( varwi ) then
          if ( optg /= 5 ) then
             ! interpolation done in all active and non-active points
             rval1  = epsdry
             epsdry = -99999.
             call SWIPOL ( wxf(1,2), ovexcv(29), xc, yc, mip, voq(1,jvx), kgrpnt, oparr(1,1) )
             call SWIPOL ( wyf(1,2), ovexcv(29), xc, yc, mip, voq(1,jvy), kgrpnt, oparr(1,1) )
             epsdry = rval1
          else
             ! interpolation done in all active and non-active points
             allocate(ltmp(nverts))
             ltmp(:)        = vert(:)%active
             vert(:)%active = .true.
             call SwanInterpolateOutput ( voq(1,jvx), voq(1,1), voq(1,2), wxf(1,2), mip, kvert, ovexcv(29) )
             call SwanInterpolateOutput ( voq(1,jvy), voq(1,1), voq(1,2), wyf(1,2), mip, kvert, ovexcv(29) )
             vert(:)%active = ltmp(:)
             deallocate(ltmp)
          endif
          do i = 1, mip
             rval1 = voq(i,jvx)
             rval2 = voq(i,jvy)
             voq(i,jvx) = coscq*rval1 - sincq*rval2
             voq(i,jvy) = sincq*rval1 + coscq*rval2
          enddo
       else
          rval1 = u10 * cos(wdip)
          rval2 = u10 * sin(wdip)
          voq(:,jvx) = coscq*rval1 - sincq*rval2
          voq(:,jvy) = sincq*rval1 + coscq*rval2
       endif
    endif
    !
    ! magnitude of wind velocity
    !
    if ( oqproc(27) ) then
       jvx = voqr(29)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,voqr(27)) = sqrt(rval1*rval1 + rval2*rval2)
       enddo
    endif
    !
    ! direction of wind velocity
    !
    if ( oqproc(28) ) then
       jvx = voqr(29)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          if ( rval1 /= 0. .or. rval2 /= 0. ) then
             if ( bnaut ) then
                dirdeg = atan2(rval2,rval1) * 180./pi
             else
                dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
             endif
             if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
             !
             ! convert (if necessary) from nautical degrees to Cartesian degrees
             !
             voq(i,voqr(28)) = DEGCNV( dirdeg )
             !
          else
             voq(i,voqr(28)) = ovexcv(28)
          endif
       enddo
    endif
    !
    ! bottom friction coefficient
    !
    if ( oqproc(30) ) then
       if ( varfr ) then
          if ( optg /= 5 ) then
             ! interpolation done in all active and non-active points
             rval1  = epsdry
             epsdry = -99999.
             call SWIPOL ( fricf(1,2), ovexcv(30), xc, yc, mip, voq(1,voqr(30)), kgrpnt, oparr(1,1) )
             epsdry = rval1
          else
             ! interpolation done in all active and non-active points
             allocate(ltmp(nverts))
             ltmp(:)        = vert(:)%active
             vert(:)%active = .true.
             call SwanInterpolateOutput ( voq(1,voqr(30)), voq(1,1), voq(1,2), fricf(1,2), mip, kvert, ovexcv(30) )
             vert(:)%active = ltmp(:)
             deallocate(ltmp)
          endif
       else
          rval1 = pbot(1)
          do i = 1, mip
             if ( .not.EQREAL(rval1,ovexcv(30)) ) voq(i,voqr(30)) = rval1
          enddo
       endif
    endif
    !
    ! friction velocities in x- and y-directions
    !
    if ( oqproc(32) ) then
       jvx = voqr(32)
       jvy = jvx + 1
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,11), ovexcv(32), xc, yc, mip, voq(1,jvx), kgrpnt, oparr(1,1) )
          call SWIPOL ( oparr(1,12), ovexcv(32), xc, yc, mip, voq(1,jvy), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,jvx), voq(1,1), voq(1,2), oparr(1,11), mip, kvert, ovexcv(32) )
          call SwanInterpolateOutput ( voq(1,jvy), voq(1,1), voq(1,2), oparr(1,12), mip, kvert, ovexcv(32) )
       endif
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,jvx) = coscq*rval1 - sincq*rval2
          voq(i,jvy) = sincq*rval1 + coscq*rval2
       enddo
    endif
    !
    ! magnitude of friction velocity
    !
    if ( oqproc(31) ) then
       jvx = voqr(32)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,voqr(31)) = sqrt(rval1*rval1 + rval2*rval2)
       enddo
    endif
    !
    ! time-averaged or mean velocities in x- and y-directions
    !
    if ( oqproc(35) ) then
       jvx = voqr(35)
       jvy = jvx + 1
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,22), ovexcv(35), xc, yc, mip, voq(1,jvx), kgrpnt, oparr(1,1) )
          call SWIPOL ( oparr(1,23), ovexcv(35), xc, yc, mip, voq(1,jvy), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,jvx), voq(1,1), voq(1,2), oparr(1,22), mip, kvert, ovexcv(35) )
          call SwanInterpolateOutput ( voq(1,jvy), voq(1,1), voq(1,2), oparr(1,23), mip, kvert, ovexcv(35) )
       endif
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,jvx) = coscq*rval1 - sincq*rval2
          voq(i,jvy) = sincq*rval1 + coscq*rval2
       enddo
    endif
    !
    ! time-averaged or mean velocity magnitude
    !
    if ( oqproc(33) ) then
       jvx = voqr(35)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,voqr(33)) = sqrt(rval1*rval1 + rval2*rval2)
       enddo
    endif
    !
    ! time-averaged or mean velocity direction
    !
    if ( oqproc(34) ) then
       jvx = voqr(35)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          if ( rval1 /= 0. .or. rval2 /= 0. ) then
             if ( bnaut ) then
                dirdeg = atan2(rval2,rval1) * 180./pi
             else
                dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
             endif
             if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
             !
             ! convert (if necessary) from nautical degrees to Cartesian degrees
             !
             voq(i,voqr(34)) = DEGCNV( dirdeg )
             !
          else
             voq(i,voqr(34)) = ovexcv(34)
          endif
       enddo
    endif
    !
    ! time-averaged or mean grid-oriented U-velocity
    !
    if ( oqproc(36) ) call SWIPOL ( oparr(1,20), ovexcv(36), xc, yc, mip, voq(1,voqr(36)), kgrpnt, oparr(1,1) )
    !
    ! time-averaged or mean grid-oriented V-velocity
    !
    if ( oqproc(37) ) then
       if ( .not.oned ) then
          call SWIPOL ( oparr(1,21), ovexcv(37), xc, yc, mip, voq(1,voqr(37)), kgrpnt, oparr(1,1) )
       else
          voq(:,voqr(37)) = ovexcv(37)
       endif
    endif
    !
    ! time-averaged or mean salinity
    !
    if ( oqproc(42) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,28), ovexcv(42), xc, yc, mip, voq(1,voqr(42)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(42)), voq(1,1), voq(1,2), oparr(1,28), mip, kvert, ovexcv(42) )
       endif
    endif
    !
    ! time-averaged or mean temperature
    !
    if ( oqproc(43) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,29), ovexcv(43), xc, yc, mip, voq(1,voqr(43)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(43)), voq(1,1), voq(1,2), oparr(1,29), mip, kvert, ovexcv(43) )
       endif
    endif
    !
    ! time-averaged or mean suspended sediment
    !
    if ( oqproc(44) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,30), ovexcv(44), xc, yc, mip, voq(1,voqr(44)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(44)), voq(1,1), voq(1,2), oparr(1,30), mip, kvert, ovexcv(44) )
       endif
    endif
    !
    ! vorticity
    !
    if ( oqproc(45) ) then
       if ( optg /= 5 ) then
          call SWIPOL ( oparr(1,31), ovexcv(45), xc, yc, mip, voq(1,voqr(45)), kgrpnt, oparr(1,1) )
       else
          call SwanInterpolateOutput ( voq(1,voqr(45)), voq(1,1), voq(1,2), oparr(1,31), mip, kvert, ovexcv(45) )
       endif
    endif
    !
    ! draft of floating object
    !
    if ( oqproc(46) ) then
       if ( optg /= 5 ) then
          ! interpolation done in all active and non-active points
          rval1  = epsdry
          epsdry = -99999.
          call SWIPOL ( flobjf, ovexcv(46), xc, yc, mip, voq(1,voqr(46)), kgrpnt, oparr(1,1) )
          epsdry = rval1
       else
          ! interpolation done in all active and non-active points
          allocate(ltmp(nverts))
          ltmp(:)        = vert(:)%active
          vert(:)%active = .true.
          call SwanInterpolateOutput ( voq(1,voqr(46)), voq(1,1), voq(1,2), flobjf, mip, kvert, ovexcv(46) )
          vert(:)%active = ltmp(:)
          deallocate(ltmp)
       endif
    endif
    !
    ! layer interfaces
    !
    if ( oqproc(51) ) then
       if ( optg /= 5 ) then
          ! interpolation done in all active and non-active points
          rval1  = epsdry
          epsdry = -99999.
          do k = 0, kmax
             call SWIPOL ( oparrk(1,k,2), ovexcv(51), xc, yc, mip, voqk(1,k,voqr(51)), kgrpnt, oparr(1,1) )
          enddo
          epsdry = rval1
       else
          ! interpolation done in all active and non-active points
          allocate(ltmp(nverts))
          ltmp(:)        = vert(:)%active
          vert(:)%active = .true.
          do k = 0, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(51)), voq(1,1), voq(1,2), oparrk(1,k,2), mip, kvert, ovexcv(51) )
          enddo
          vert(:)%active = ltmp(:)
          deallocate(ltmp)
       endif
    endif
    !
    ! physical vertical velocity
    !
    if ( oqproc(52) ) then
       if ( optg /= 5 ) then
          do k = 0, kmax
             call SWIPOL ( oparrk(1,k,12), ovexcv(52), xc, yc, mip, voqk(1,k,voqr(52)), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 0, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(52)), voq(1,1), voq(1,2), oparrk(1,k,12), mip, kvert, ovexcv(52) )
          enddo
       endif
    endif
    !
    ! relative vertical velocity
    !
    if ( oqproc(53) ) then
       if ( optg /= 5 ) then
          do k = 0, kmax
             call SWIPOL ( oparrk(1,k,11), ovexcv(53), xc, yc, mip, voqk(1,k,voqr(53)), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 0, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(53)), voq(1,1), voq(1,2), oparrk(1,k,11), mip, kvert, ovexcv(53) )
          enddo
       endif
    endif
    !
    ! turbulent kinetic energy
    !
    if ( oqproc(54) ) then
       if ( iturb /= 0 ) then
          if ( optg /= 5 ) then
             do k = 0, kmax
                call SWIPOL ( oparrk(1,k,15), ovexcv(54), xc, yc, mip, voqk(1,k,voqr(54)), kgrpnt, oparr(1,1) )
             enddo
          else
             do k = 0, kmax
                call SwanInterpolateOutput ( voqk(1,k,voqr(54)), voq(1,1), voq(1,2), oparrk(1,k,15), mip, kvert, ovexcv(54) )
             enddo
          endif
       else
          voqk(:,:,voqr(54)) = ovexcv(54)
       endif
    endif
    !
    ! dissipation rate of turbulent kinetic energy
    !
    if ( oqproc(55) ) then
       if ( iturb /= 0 ) then
          if ( optg /= 5 ) then
             do k = 0, kmax
                call SWIPOL ( oparrk(1,k,16), ovexcv(55), xc, yc, mip, voqk(1,k,voqr(55)), kgrpnt, oparr(1,1) )
             enddo
          else
             do k = 0, kmax
                call SwanInterpolateOutput ( voqk(1,k,voqr(55)), voq(1,1), voq(1,2), oparrk(1,k,16), mip, kvert, ovexcv(55) )
             enddo
          endif
       else
          voqk(:,:,voqr(55)) = ovexcv(55)
       endif
    endif
    !
    ! vertical eddy viscosity
    !
    if ( oqproc(56) ) then
       if ( optg /= 5 ) then
          do k = 0, kmax
             call SWIPOL ( oparrk(1,k,17), ovexcv(56), xc, yc, mip, voqk(1,k,voqr(56)), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 0, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(56)), voq(1,1), voq(1,2), oparrk(1,k,17), mip, kvert, ovexcv(56) )
          enddo
       endif
    endif
    !
    ! mean turbulent kinetic energy
    !
    if ( oqproc(57) ) then
       if ( iturb /= 0 ) then
          if ( optg /= 5 ) then
             do k = 0, kmax
                call SWIPOL ( oparrk(1,k,22), ovexcv(57), xc, yc, mip, voqk(1,k,voqr(57)), kgrpnt, oparr(1,1) )
             enddo
          else
             do k = 0, kmax
                call SwanInterpolateOutput ( voqk(1,k,voqr(57)), voq(1,1), voq(1,2), oparrk(1,k,22), mip, kvert, ovexcv(57) )
             enddo
          endif
       else
          voqk(:,:,voqr(57)) = ovexcv(57)
       endif
    endif
    !
    ! mean dissipation rate of turbulent kinetic energy
    !
    if ( oqproc(58) ) then
       if ( iturb /= 0 ) then
          if ( optg /= 5 ) then
             do k = 0, kmax
                call SWIPOL ( oparrk(1,k,23), ovexcv(58), xc, yc, mip, voqk(1,k,voqr(58)), kgrpnt, oparr(1,1) )
             enddo
          else
             do k = 0, kmax
                call SwanInterpolateOutput ( voqk(1,k,voqr(58)), voq(1,1), voq(1,2), oparrk(1,k,23), mip, kvert, ovexcv(58) )
             enddo
          endif
       else
          voqk(:,:,voqr(58)) = ovexcv(58)
       endif
    endif
    !
    ! mean vertical eddy viscosity
    !
    if ( oqproc(59) ) then
       if ( optg /= 5 ) then
          do k = 0, kmax
             call SWIPOL ( oparrk(1,k,24), ovexcv(59), xc, yc, mip, voqk(1,k,voqr(59)), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 0, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(59)), voq(1,1), voq(1,2), oparrk(1,k,24), mip, kvert, ovexcv(59) )
          enddo
       endif
    endif
    !
    ! layer thicknesses
    !
    if ( oqproc(71) ) then
       if ( optg /= 5 ) then
          ! interpolation done in all active and non-active points
          rval1  = epsdry
          epsdry = -99999.
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,1), ovexcv(71), xc, yc, mip, voqk(1,k,voqr(71)), kgrpnt, oparr(1,1) )
          enddo
          epsdry = rval1
       else
          ! interpolation done in all active and non-active points
          allocate(ltmp(nverts))
          ltmp(:)        = vert(:)%active
          vert(:)%active = .true.
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(71)), voq(1,1), voq(1,2), oparrk(1,k,1), mip, kvert, ovexcv(71) )
          enddo
          vert(:)%active = ltmp(:)
          deallocate(ltmp)
       endif
    endif
    !
    ! layer-dependent flow velocities in x- and y-directions
    !
    if ( oqproc(74) ) then
       jvx = voqr(74)
       jvy = jvx + 1
       if ( optg /= 5 ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,5), ovexcv(74), xc, yc, mip, voqk(1,k,jvx), kgrpnt, oparr(1,1) )
             call SWIPOL ( oparrk(1,k,6), ovexcv(74), xc, yc, mip, voqk(1,k,jvy), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,jvx), voq(1,1), voq(1,2), oparrk(1,k,5), mip, kvert, ovexcv(74) )
             call SwanInterpolateOutput ( voqk(1,k,jvy), voq(1,1), voq(1,2), oparrk(1,k,6), mip, kvert, ovexcv(74) )
          enddo
       endif
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,jvx) = coscq*rval1 - sincq*rval2
             voqk(i,k,jvy) = sincq*rval1 + coscq*rval2
          enddo
       enddo
    endif
    !
    ! layer-dependent velocity magnitude
    !
    if ( oqproc(72) ) then
       jvx = voqr(74)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,voqr(72)) = sqrt(rval1*rval1 + rval2*rval2)
          enddo
       enddo
    endif
    !
    ! layer-dependent velocity direction
    !
    if ( oqproc(73) ) then
       jvx = voqr(74)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             if ( rval1 /= 0. .or. rval2 /= 0. ) then
                if ( bnaut ) then
                   dirdeg = atan2(rval2,rval1) * 180./pi
                else
                   dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
                endif
                if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
                !
                ! convert (if necessary) from nautical degrees to Cartesian degrees
                !
                voqk(i,k,voqr(73)) = DEGCNV( dirdeg )
                !
             else
                voqk(i,k,voqr(73)) = ovexcv(73)
             endif
          enddo
       enddo
    endif
    !
    ! layer-dependent grid-oriented U-velocity
    !
    if ( oqproc(75) ) then
       do k = 1, kmax
          call SWIPOL ( oparrk(1,k,3), ovexcv(75), xc, yc, mip, voqk(1,k,voqr(75)), kgrpnt, oparr(1,1) )
       enddo
    endif
    !
    ! layer-dependent grid-oriented V-velocity
    !
    if ( oqproc(76) ) then
       if ( .not.oned ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,4), ovexcv(76), xc, yc, mip, voqk(1,k,voqr(76)), kgrpnt, oparr(1,1) )
          enddo
       else
          voqk(:,:,voqr(76)) = ovexcv(76)
       endif
    endif
    !
    ! layer-dependent specific discharges in x- and y-directions
    !
    if ( oqproc(79) ) then
       jvx = voqr(79)
       jvy = jvx + 1
       if ( optg /= 5 ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k, 9), ovexcv(79), xc, yc, mip, voqk(1,k,jvx), kgrpnt, oparr(1,1) )
             call SWIPOL ( oparrk(1,k,10), ovexcv(79), xc, yc, mip, voqk(1,k,jvy), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,jvx), voq(1,1), voq(1,2), oparrk(1,k, 9), mip, kvert, ovexcv(79) )
             call SwanInterpolateOutput ( voqk(1,k,jvy), voq(1,1), voq(1,2), oparrk(1,k,10), mip, kvert, ovexcv(79) )
          enddo
       endif
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,jvx) = coscq*rval1 - sincq*rval2
             voqk(i,k,jvy) = sincq*rval1 + coscq*rval2
          enddo
       enddo
    endif
    !
    ! layer-dependent magnitude of specific discharge
    !
    if ( oqproc(77) ) then
       jvx = voqr(79)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,voqr(77)) = sqrt(rval1*rval1 + rval2*rval2)
          enddo
       enddo
    endif
    !
    ! layer-dependent direction of specific discharge
    !
    if ( oqproc(78) ) then
       jvx = voqr(79)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             if ( rval1 /= 0. .or. rval2 /= 0. ) then
                if ( bnaut ) then
                   dirdeg = atan2(rval2,rval1) * 180./pi
                else
                   dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
                endif
                if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
                !
                ! convert (if necessary) from nautical degrees to Cartesian degrees
                !
                voqk(i,k,voqr(78)) = DEGCNV( dirdeg )
                !
             else
                voqk(i,k,voqr(78)) = ovexcv(78)
             endif
          enddo
       enddo
    endif
    !
    ! layer-dependent grid-oriented U-discharge per unit width
    !
    if ( oqproc(80) ) then
       do k = 1, kmax
          call SWIPOL ( oparrk(1,k,7), ovexcv(80), xc, yc, mip, voqk(1,k,voqr(80)), kgrpnt, oparr(1,1) )
       enddo
    endif
    !
    ! layer-dependent grid-oriented V-discharge per unit width
    !
    if ( oqproc(81) ) then
       if ( .not.oned ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,8), ovexcv(81), xc, yc, mip, voqk(1,k,voqr(81)), kgrpnt, oparr(1,1) )
          enddo
       else
          voqk(:,:,voqr(81)) = ovexcv(81)
       endif
    endif
    !
    ! layer-dependent non-hydrostatic pressure
    !
    if ( oqproc(82) ) then
       if ( ihydro == 1 .or. ihydro == 2 ) then
          if ( optg /= 5 ) then
             do k = 1, kmax
                call SWIPOL ( oparrk(1,k,13), ovexcv(82), xc, yc, mip, voqk(1,k,voqr(82)), kgrpnt, oparr(1,1) )
             enddo
          else
             do k = 1, kmax
                call SwanInterpolateOutput ( voqk(1,k,voqr(82)), voq(1,1), voq(1,2), oparrk(1,k,13), mip, kvert, ovexcv(82) )
             enddo
          endif
       else
          voqk(:,:,voqr(82)) = ovexcv(82)
       endif
    endif
    !
    ! layer-dependent pressure
    !
    if ( oqproc(83) ) then
       if ( optg /= 5 ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,14), ovexcv(83), xc, yc, mip, voqk(1,k,voqr(83)), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(83)), voq(1,1), voq(1,2), oparrk(1,k,14), mip, kvert, ovexcv(83) )
          enddo
       endif
    endif
    !
    ! time-averaged or mean layer-dependent flow velocities in x- and y-directions
    !
    if ( oqproc(86) ) then
       jvx = voqr(86)
       jvy = jvx + 1
       if ( optg /= 5 ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,20), ovexcv(86), xc, yc, mip, voqk(1,k,jvx), kgrpnt, oparr(1,1) )
             call SWIPOL ( oparrk(1,k,21), ovexcv(86), xc, yc, mip, voqk(1,k,jvy), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,jvx), voq(1,1), voq(1,2), oparrk(1,k,20), mip, kvert, ovexcv(86) )
             call SwanInterpolateOutput ( voqk(1,k,jvy), voq(1,1), voq(1,2), oparrk(1,k,21), mip, kvert, ovexcv(86) )
          enddo
       endif
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,jvx) = coscq*rval1 - sincq*rval2
             voqk(i,k,jvy) = sincq*rval1 + coscq*rval2
          enddo
       enddo
    endif
    !
    ! time-averaged or mean layer-dependent velocity magnitude
    !
    if ( oqproc(84) ) then
       jvx = voqr(86)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,voqr(84)) = sqrt(rval1*rval1 + rval2*rval2)
          enddo
       enddo
    endif
    !
    ! time-averaged or mean layer-dependent velocity direction
    !
    if ( oqproc(85) ) then
       jvx = voqr(86)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             if ( rval1 /= 0. .or. rval2 /= 0. ) then
                if ( bnaut ) then
                   dirdeg = atan2(rval2,rval1) * 180./pi
                else
                   dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
                endif
                if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
                !
                ! convert (if necessary) from nautical degrees to Cartesian degrees
                !
                voqk(i,k,voqr(85)) = DEGCNV( dirdeg )
                !
             else
                voqk(i,k,voqr(85)) = ovexcv(85)
             endif
          enddo
       enddo
    endif
    !
    ! time-averaged or mean layer-dependent grid-oriented U-velocity
    !
    if ( oqproc(87) ) then
       do k = 1, kmax
          call SWIPOL ( oparrk(1,k,18), ovexcv(87), xc, yc, mip, voqk(1,k,voqr(87)), kgrpnt, oparr(1,1) )
       enddo
    endif
    !
    ! time-averaged or mean layer-dependent grid-oriented V-velocity
    !
    if ( oqproc(88) ) then
       if ( .not.oned ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,19), ovexcv(88), xc, yc, mip, voqk(1,k,voqr(88)), kgrpnt, oparr(1,1) )
          enddo
       else
          voqk(:,:,voqr(88)) = ovexcv(88)
       endif
    endif
    !
    ! layer-dependent salinity
    !
    if ( oqproc(89) ) then
       if ( optg /= 5 ) then
          ! interpolation done in all active and non-active points
          rval1  = epsdry
          epsdry = -99999.
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,25), ovexcv(89), xc, yc, mip, voqk(1,k,voqr(89)), kgrpnt, oparr(1,1) )
          enddo
          epsdry = rval1
       else
          ! interpolation done in all active and non-active points
          allocate(ltmp(nverts))
          ltmp(:)        = vert(:)%active
          vert(:)%active = .true.
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(89)), voq(1,1), voq(1,2), oparrk(1,k,25), mip, kvert, ovexcv(89) )
          enddo
          vert(:)%active = ltmp(:)
          deallocate(ltmp)
       endif
    endif
    !
    ! layer-dependent temperature
    !
    if ( oqproc(90) ) then
       if ( optg /= 5 ) then
          ! interpolation done in all active and non-active points
          rval1  = epsdry
          epsdry = -99999.
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,26), ovexcv(90), xc, yc, mip, voqk(1,k,voqr(90)), kgrpnt, oparr(1,1) )
          enddo
          epsdry = rval1
       else
          ! interpolation done in all active and non-active points
          allocate(ltmp(nverts))
          ltmp(:)        = vert(:)%active
          vert(:)%active = .true.
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(90)), voq(1,1), voq(1,2), oparrk(1,k,26), mip, kvert, ovexcv(90) )
          enddo
          vert(:)%active = ltmp(:)
          deallocate(ltmp)
       endif
    endif
    !
    ! layer-dependent suspended sediment
    !
    if ( oqproc(91) ) then
       if ( optg /= 5 ) then
          ! interpolation done in all active and non-active points
          rval1  = epsdry
          epsdry = -99999.
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,27), ovexcv(91), xc, yc, mip, voqk(1,k,voqr(91)), kgrpnt, oparr(1,1) )
          enddo
          epsdry = rval1
       else
          ! interpolation done in all active and non-active points
          allocate(ltmp(nverts))
          ltmp(:)        = vert(:)%active
          vert(:)%active = .true.
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(91)), voq(1,1), voq(1,2), oparrk(1,k,27), mip, kvert, ovexcv(91) )
          enddo
          vert(:)%active = ltmp(:)
          deallocate(ltmp)
       endif
    endif
    !
    ! time-averaged or mean layer-dependent salinity
    !
    if ( oqproc(92) ) then
       if ( optg /= 5 ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,28), ovexcv(92), xc, yc, mip, voqk(1,k,voqr(92)), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(92)), voq(1,1), voq(1,2), oparrk(1,k,28), mip, kvert, ovexcv(92) )
          enddo
       endif
    endif
    !
    ! time-averaged or mean layer-dependent temperature
    !
    if ( oqproc(93) ) then
       if ( optg /= 5 ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,29), ovexcv(93), xc, yc, mip, voqk(1,k,voqr(93)), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(93)), voq(1,1), voq(1,2), oparrk(1,k,29), mip, kvert, ovexcv(93) )
          enddo
       endif
    endif
    !
    ! time-averaged or mean layer-dependent suspended sediment
    !
    if ( oqproc(94) ) then
       if ( optg /= 5 ) then
          do k = 1, kmax
             call SWIPOL ( oparrk(1,k,30), ovexcv(94), xc, yc, mip, voqk(1,k,voqr(94)), kgrpnt, oparr(1,1) )
          enddo
       else
          do k = 1, kmax
             call SwanInterpolateOutput ( voqk(1,k,voqr(94)), voq(1,1), voq(1,2), oparrk(1,k,30), mip, kvert, ovexcv(94) )
          enddo
       endif
    endif
    !
    ! correct problem coordinates with offset values
    !
    do i = 1, mip
       xp = voq(i,1)
       if ( .not.EQREAL(xp,ovexcv(1)) ) voq(i,1) = xp + xoffs
       yp = voq(i,2)
       if ( .not.EQREAL(yp,ovexcv(2)) ) voq(i,2) = yp + yoffs
    enddo
    !
    ! in case of parallel run, mark location points inside own subdomain
    !
    if ( PARLL ) then
       if ( rtype /= 'BLKV' ) then
          ihlx = IHALOX
          ihly = IHALOY
       else
          ihlx = 2
          ihly = 2
       endif
       ixb = 1 + ihlx
       if ( LMXF ) ixb = 1
       ixe = mxc - ihlx
       if ( LMXL ) ixe = mxc
       iyb = 1 + ihly
       if ( LMYF ) iyb = 1
       iye = myc - ihly
       if ( LMYL ) iye = myc
       do i = 1, mip
          ix = int(xc(i))
          iy = int(yc(i))
          rval1 = float(ix)
          rval2 = float(iy)
          if ( .not.EQREAL(xc(i),rval1) .or. EQREAL(xc(i),0.) .or. rtype == 'BLKV' ) ix = ix + 1
          if ( .not.EQREAL(yc(i),rval2) .or. EQREAL(yc(i),0.) .or. rtype == 'BLKV' ) iy = iy + 1
          if ( ix >= ixb .and. ix <= ixe .and. iy >= iyb .and. iy <= iye ) ionod(i) = INODE
       enddo
    endif
    !
    if (allocated(kvert)) deallocate(kvert)
    !
 201 format (' entry SwashQuanOutp ', 10l2, i8)
    !
end subroutine SwashQuanOutp
