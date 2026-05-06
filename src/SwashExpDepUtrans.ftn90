subroutine SwashExpDepUtrans ( rp, rpo, rp1, rp0, rpi, u1, qn, fluxt )
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
!   Performs the time integration for the depth-averaged transport equations on triangular mesh
!
!   Method
!
!   The time integration is fully explicit. Local time stepping is employed.
!
!   The space discretization is a cell-centered finite volume discretization and
!   is consistent with the discretization of the global continuity equation.
!
!   The concentration is located at the centroid of the triangular cell.
!
!   The discretization is mass conservative and complies with the discrete maximum principle.
!
!   The advective term is approximated by either first order upwind or higher order (flux-limited) scheme
!   (CDS, Fromm, BDF, QUICK, MUSCL, Koren, etc.). For the r-ratio the most upwave vertex of upwind cell is used.
!
!   The Thatcher-Harleman boundary condition is imposed at sea side for unsteady salt intrusion.
!   The constituent return time is given by the user.
!
!   M.L. Thatcher and R.F. Harleman
!   A mathematical model for the prediction of unsteady salinity intrusion in estuaries
!   Technical report 144, MIT, Massachusetts, USA, 1972
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
    use SwashFlowdata, fluxttmp => fluxt, &
                       qntmp    => qn   , &
                       rptmp    => rp   , &
                       rpotmp   => rpo  , &
                       rp1tmp   => rp1  , &
                       rp0tmp   => rp0  , &
                       rpitmp   => rpi  , &
                       u1tmp    => u1
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    real, dimension(nfaces)             , intent(out)   :: fluxt            ! total flux at faces
    real, dimension(nfaces)             , intent(inout) :: qn               ! mass flux at triangular face
    real, dimension(ncells,ltrans)      , intent(inout) :: rp               ! concentration at current time level
    real, dimension(ncells)             , intent(out)   :: rp0              ! concentration at previous substep for local time stepping
    real, dimension(ncells)             , intent(out)   :: rp1              ! concentration at current substep for local time stepping
    real, dimension(ncells,0:2**mtamx-1), intent(out)   :: rpi              ! intermediate transport constituent at various time-step levels
    real, dimension(ncells)             , intent(out)   :: rpo              ! concentration at previous time level
    real, dimension(nfaces)             , intent(in)    :: u1               ! face velocity at current time level
!
!   Parameter variables
!
    real, parameter                                     :: nuval = -999999. ! a special value that marks it as non-used
!
!   Local variables
!
    integer                                             :: btype            ! boundary type (see SwashUpdateUData.f90)
    integer                                             :: icell            ! cell index / loop counter over cells
    integer                                             :: icelll           ! left cell of present face
    integer                                             :: icellr           ! right cell of present face
    integer                                             :: icistb           ! counter for number of instable points
    integer, save                                       :: ient = 0         ! number of entries in this subroutine
    integer                                             :: iface            ! face index / loop counter over faces
    integer                                             :: j                ! loop counter
    integer                                             :: jc               ! loop counter over cells of vertex
    integer                                             :: jf               ! loop counter over faces (of cell)
    integer                                             :: jj               ! offset index in range of substeps
    integer                                             :: l                ! loop counter over constituents
    integer                                             :: m                ! loop counter over time-step levels
    integer                                             :: maxsb            ! maximum number of substeps
    integer                                             :: mtmax            ! maximum time-step level used in local time stepping
    integer                                             :: n                ! loop counter over substeps
    integer                                             :: nc               ! number of cells corresponding to specific time-step level
    integer                                             :: nsub             ! number of substeps in local time stepping
    integer, dimension(3)                               :: v                ! vertices of present cell
    integer                                             :: vf1              ! first vertex of present face
    integer                                             :: vf2              ! second vertex of present face
    integer                                             :: vu               ! upwind vertex
    !
    real                                                :: area             ! area of present cell
    real                                                :: areal            ! area of left cell of present face
    real                                                :: arear            ! area of right cell of present face
    real                                                :: cfl              ! cell-based CFL number
    real                                                :: contrib          ! total contribution of transport flux per cell
    real                                                :: dh               ! local increment in water depth
    real                                                :: dif2d            ! horizontal eddy diffusivity coefficient in velocity point
    real                                                :: fac              ! a factor
    real                                                :: finp             ! interpolation factor
    real                                                :: fluxlim          ! flux limiter
    real                                                :: grad1            ! solution gradient
    real                                                :: grad2            ! another solution gradient
    real                                                :: h0               ! water depth at previous substep
    real                                                :: h1               ! water depth at current substep
    real                                                :: lf               ! length of present face
    real                                                :: mass             ! total mass
    real                                                :: psm              ! Prandtl-Schmidt number
    real                                                :: qf               ! mass flux
    real                                                :: rdx              ! reciprocal of distance between circumcenters adjacent to face
    real                                                :: rnsb             ! reciprocal of number of substeps
    real                                                :: rproc            ! auxiliary variable with percentage of instable points
    real                                                :: rpu              ! averaged concentration in upwind vertex
    real                                                :: rsgn             ! sign for indicating face orientation
    real                                                :: rval             ! auxiliary real
    real                                                :: stabmx           ! auxiliary variable with maximum diffusivity based stability criterion
    real                                                :: sumqf            ! sum of outgoing mass fluxes per cell
    real                                                :: totarea          ! total area of all cells around vertex
    !
    character(120)                                      :: msgstr           ! string to pass message
    !
    type(verttype), dimension(:), pointer               :: vert             ! datastructure for vertices with their attributes
    type(celltype), dimension(:), pointer               :: cell             ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer               :: face             ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashExpDepUtrans')
    !
    ! point to vertex, cell and face objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    ! if momentum equation has been skipped, compute the mass flux
    !
    if ( momskip ) then
       !
       do iface = 1, nfaces
          !
          lf = face(iface)%attr(FACELEN)
          !
          qn(iface) = lf * hu(iface) * u1(iface)
          !
       enddo
       !
    endif
    !
    ! compute the time-step level for each cell
    !
    mlev = 0
    !
    do icell = 1, ncells
       !
       ! compute the sum of mass fluxes leaving the cell
       !
       sumqf = 0.
       !
       do jf = 1, cell(icell)%nof
          !
          ! face identifier
          !
          iface = cell(icell)%face(jf)%atti(FACEID)
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
          ! compute outgoing mass flux at current face
          !
          qf = rsgn * qn(iface)
          !
          if ( qf > 0. ) sumqf = sumqf + qf
          !
       enddo
       !
       ! compute the "flow" Courant number
       ! note: it must not be larger than 0.5 to fulfill the max-min property
       !
       if ( hso(icell) > epsdry ) then
          !
          area = cell(icell)%attr(CELLAREA)
          !
          cfl = sumqf * dt / hso(icell) / area
          !
          mlev(icell) = max( 0, 1+floor( log(2.*cfl)/log(2.) ) )
          !
       endif
       !
    enddo
    !
    ! smooth transition between various time-step levels
    !
    do icell = 1, ncells
       !
       mtmax = 0
       !
       ! loop over faces of the cell
       !
       do jf = 1, cell(icell)%nof
          !
          ! face identifier
          !
          iface = cell(icell)%face(jf)%atti(FACEID)
          !
          ! consider left and right cells of current face
          !
          icelll = face(iface)%atti(FACECL)
          icellr = face(iface)%atti(FACECR)
          !
          ! take maximum time-step level of neighbouring cell
          !
          if ( icell == icelll .and. icellr /= 0 ) then
             !
             if ( mlev(icellr) > mtmax ) mtmax = mlev(icellr)
             !
          elseif ( icell == icellr .and. icelll /= 0 ) then
             !
             if ( mlev(icelll) > mtmax ) mtmax = mlev(icelll)
             !
          endif
          !
       enddo
       !
       mlev(icell) = max( mtmax, mlev(icell) )
       !
    enddo
    !
    ! compute maximum time-step level
    !
    mtmax = maxval( mlev )
    if ( mtmax > mtamx ) then
       mtmax = mtamx
       call msgerr ( 1, 'The maximum number of time-step level is exceeded!')
       call msgerr ( 0, 'It is advised to reduce the time step!' )
    endif
    if ( ITEST >= 30 ) write (PRINTF,101) mtmax
    !
    ! compute actual time step for each cell
    !
    dtc = dt / 2.**mlev
    !
    ! loop over constituents
    !
    do l = 1, ltrans
       !
       ! determine Prandtl-Schmidt number
       !
       if ( l == lsal .or. l == ltemp ) then
          psm = 0.7
       else if ( l == lsed ) then
          psm = 1.0
       endif
       !
       icistb = 0
       !
       ! store old concentrations
       !
       rpo(:) = rp(:,l)
       !
       ! compute advective flux and concentration at open boundaries using boundary conditions
       !
       do jf = 1, nfacesb
          !
          iface = jbface(jf)
          icell = face(iface)%atti(FACEC1)
          !
          ! get length of boundary face
          !
          lf = face(iface)%attr(FACELEN)
          !
          ! get area of boundary cell
          !
          area = cell(icell)%attr(CELLAREA)
          !
          btype = face(iface)%atti(FBTYPE)
          !
          if ( btype /= 1 ) then
             !
             ! --- boundary face open
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             ! take into account orientation of the current face
             !
             if ( icell == icelll ) then
                ! leaving the cell
                rsgn =  1.
             else if ( icell == icellr ) then
                ! entering the cell
                rsgn = -1.
             endif
             !
             ! compute mass flux at boundary face
             !
             qf = qn(iface)
             !
             if ( rsgn * qf > 0. ) then
                !
                ! outflow
                !
                fluxt(iface) = qf * rpo(icell)
                !
                bcrp(jf,1,l) = bcrpo(jf,1,l) - dt * rsgn * lf * u1(iface) * ( bcrpo(jf,1,l) - rpo(icell) ) / area
                !
                if ( l == lsal ) then
                   coutu (jf,1) = bcrp(jf,1,l)
                   icretu(jf,1) = tcret
                endif
                !
             else
                !
                ! inflow
                !
                fluxt(iface) = qf * bcrpo(jf,1,l)
                !
                if ( l == lsal ) then
                   fac = max(icretu(jf,1),0.) / max(tcret,dt)
                   bcrp(jf,1,l) = coutu(jf,1) + 0.5 * ( cbndu(jf,1,l) - coutu(jf,1) ) * ( 1. + cos(fac*pi) )
                   if ( .not. icretu(jf,1) < 0. ) icretu(jf,1) = icretu(jf,1) - dt
                else
                   bcrp(jf,1,l) = cbndu(jf,1,l)
                endif
                !
             endif
             !
          else
             !
             ! --- boundary face closed: no advective flux
             !
             fluxt(iface    ) = 0.
             bcrp (jf   ,1,l) = rp(icell,l)
             !
          endif
          !
       enddo
       !
       ! initialize intermediate solution at all time-step levels
       !
       rpi = nuval
       !
       ! determine maximum number of substeps of the highest time-step level
       !
       maxsb = 2**mtmax
       !
       ! loop over time-step levels, starting with the coarsest one
       !
       do m = 0, mtmax
          !
          ! determine number of substeps of current time-step level
          !
          nsub = 2**m
          rnsb = 1./real(nsub)
          !
          jj = maxsb / 2
          !
          ! construct cell-based index table for time-step level m
          !
          nc = 0
          !
          do icell = 1, ncells
             !
             if ( mlev(icell) == m ) then
                nc = nc + 1
                mcell(nc) = icell
             endif
             !
          enddo
          !
          ! initialize solution to the concentration at previous time step n
          !
          rp1 = rpo
          !
          ! start local time stepping
          !
          do n = 1, nsub
             !
             ! store concentration to previous substep
             !
             rp0 = rp1
             !
             ! compute advective and diffusive fluxes at internal faces
             !
             do iface = 1, nfaces
                !
                if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
                   !
                   finp = face(iface)%attr(FACELINPF)
                   !
                   ! get vertices of current face
                   !
                   vf1 = face(iface)%atti(FACEV1)
                   vf2 = face(iface)%atti(FACEV2)
                   !
                   ! get length of current face
                   !
                   lf = face(iface)%attr(FACELEN)
                   !
                   ! get reciprocal distance between the centroids adjacent to current face
                   !
                   rdx = face(iface)%attr(FACEDISTG)
                   !
                   ! consider left and right cells of current face
                   !
                   icelll = face(iface)%atti(FACECL)
                   icellr = face(iface)%atti(FACECR)
                   !
                   areal = cell(icelll)%attr(CELLAREA)
                   arear = cell(icellr)%attr(CELLAREA)
                   !
                   propsc = nint(pnums(46))
                   kappa  = pnums(47)
                   mbound = pnums(48)
                   phieby = pnums(49)
                   !
                   ! compute mass flux at current face
                   !
                   qf = qn(iface)
                   !
                   ! compute the advective flux based on first order upwind
                   !
                   if ( qf > 0. ) then
                      !
                      fluxt(iface) = qf * rp0(icelll)
                      !
                   else
                      !
                      fluxt(iface) = qf * rp0(icellr)
                      !
                   endif
                   !
                   ! add second order approximation
                   !
                   if ( propsc /= 1 ) then
                      !
                      if ( qf > 0. ) then
                         !
                         ! get vertices of upwind cell
                         !
                         v(1) = cell(icelll)%atti(CELLV1)
                         v(2) = cell(icelll)%atti(CELLV2)
                         v(3) = cell(icelll)%atti(CELLV3)
                         !
                         ! search for most upwave vertex
                         !
                         do j = 1, 3
                            if ( v(j) /= vf1 .and. v(j) /= vf2 ) then
                               vu = v(j)
                               exit
                            endif
                         enddo
                         !
                         ! compute area-weighted averaged concentration at upwave vertex
                         !
                         rpu     = 0.
                         totarea = 0.
                         !
                         do jc = 1, vert(vu)%noc
                            !
                            icell = vert(vu)%cell(jc)%atti(CELLID)
                            !
                            area = cell(icell)%attr(CELLAREA)
                            !
                            rpu = rpu + area * rp0(icell)
                            !
                            totarea = totarea + area
                            !
                         enddo
                         !
                         rpu = rpu / totarea
                         !
                         ! compute solution gradients
                         !
                         grad1 = rp0(icellr) - rp0(icelll)
                         grad2 = rp0(icelll) - rpu
                         !
                         ! update flux
                         !
                         fluxt(iface) = fluxt(iface) + 0.5 * qf * fluxlim(grad1,grad2)
                         !
                      else
                         !
                         ! get vertices of upwind cell
                         !
                         v(1) = cell(icellr)%atti(CELLV1)
                         v(2) = cell(icellr)%atti(CELLV2)
                         v(3) = cell(icellr)%atti(CELLV3)
                         !
                         ! search for most upwave vertex
                         !
                         do j = 1, 3
                            if ( v(j) /= vf1 .and. v(j) /= vf2 ) then
                               vu = v(j)
                               exit
                            endif
                         enddo
                         !
                         ! compute area-weighted averaged concentration at upwave vertex
                         !
                         rpu     = 0.
                         totarea = 0.
                         !
                         do jc = 1, vert(vu)%noc
                            !
                            icell = vert(vu)%cell(jc)%atti(CELLID)
                            !
                            area = cell(icell)%attr(CELLAREA)
                            !
                            rpu = rpu + area * rp0(icell)
                            !
                            totarea = totarea + area
                            !
                         enddo
                         !
                         rpu = rpu / totarea
                         !
                         ! compute solution gradients
                         !
                         grad1 = rp0(icelll) - rp0(icellr)
                         grad2 = rp0(icellr) - rpu
                         !
                         ! update flux
                         !
                         fluxt(iface) = fluxt(iface) + 0.5 * qf * fluxlim(grad1,grad2)
                         !
                      endif
                      !
                   endif
                   !
                   ! compute effective horizontal diffusivity coefficient at current face
                   !
                   if ( hdiff > 0. ) then
                      !
                      dif2d = hdiff
                      !
                   else
                      !
                      if ( ihvisc == 2 .or. ihvisc == 3 ) then
                         !
                         dif2d = vnu2d(iface) / psm
                         !
                      else
                         !
                         dif2d = 0.
                         !
                      endif
                      !
                   endif
                   !
                   rval   = max( dtc(icelll), dtc(icellr) )
                   stabmx = 0.5 * min( areal, arear ) / rval
                   !
                   ! check stability
                   !
                   if ( .not. dif2d < stabmx ) then
                      dif2d  = stabmx
                      icistb = icistb + 1
                   endif
                   !
                   ! compute the diffusive flux and update total flux
                   !
                   fluxt(iface) = fluxt(iface) - ( finp*hso(icelll)+(1.-finp)*hso(icellr) ) * dif2d * lf * rdx * ( rp0(icellr) - rp0(icelll) )
                   !
                else if ( face(iface)%atti(FMARKER) == 0 ) then
                   !
                   fluxt(iface) = 0.
                   !
                endif
                !
             enddo
             !
             !  compute concentration in cells (based on finite volume approach)
             !
             do j = 1, nc
                !
                ! consider cells with time-step level m
                !
                icell = mcell(j)
                !
                contrib = 0.
                !
                ! loop over faces of the cell
                !
                do jf = 1, cell(icell)%nof
                   !
                   ! face identifier
                   !
                   iface = cell(icell)%face(jf)%atti(FACEID)
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
                   ! get advective/diffusive flux at current face and add to other faces of the cell
                   !
                   contrib = contrib + rsgn * fluxt(iface)
                   !
                enddo
                !
                ! update in wet cell
                !
                if ( hs(icell) > epsdry ) then
                   !
                   area = cell(icell)%attr(CELLAREA)
                   !
                   ! compute water depths at current and previous substeps
                   !
                   dh = rnsb * ( hs(icell) - hso(icell) )
                   !
                   h0 = hso(icell) + real(n-1) * dh
                   h1 = h0         + dh
                   !
                   rp1(icell) = h0 * rp0(icell) - dtc(icell) * contrib / area
                   !
                   ! compute depth-averaged concentration
                   !
                   rp1(icell) = rp1(icell) / h1
                   !
                endif
                !
             enddo
             !
             ! also use intermediate solution of coarser time-step levels
             !
             do icell = 1, ncells
                !
                if ( mlev(icell) < m ) then
                   !
                   rval = rpi(icell,n*maxsb)
                   if (rval > nuval ) rp1(icell) = rval
                   !
                endif
                !
             enddo
             !
             ! finally, estimate intermediate concentration to be used at finer time-step levels
             !
             do j = 1, nc
                !
                ! consider cells with time-step level m
                !
                icell = mcell(j)
                !
                rpi(icell,jj+(n-1)*maxsb) = 0.5 * ( rp0(icell) + rp1(icell) )
                !
             enddo
             !
          enddo
          !
          ! store final solution at time step n+1
          !
          do j = 1, nc
             !
             ! consider cells with time-step level m
             !
             icell = mcell(j)
             !
             rp(icell,l) = rp1(icell)
             !
          enddo
          !
          ! compute maximum number of substeps with respect to the next level
          !
          maxsb = maxsb / 2
          !
       enddo
       !
       ! calculate total mass
       !
       if ( ITEST >= 30 ) then
          !
          mass = 0.
          !
          do icell = 1, ncells
             !
             ! area of cell
             !
             area = cell(icell)%attr(CELLAREA)
             !
             mass = mass + hs(icell) * area * rp(icell,l)
             !
          enddo
          !
          if ( l == lsal  ) write(PRINTF,102) mass
          if ( l == ltemp ) write(PRINTF,103) mass
          if ( l == lsed  ) write(PRINTF,104) mass*rhos
          !
       endif
       !
       ! give warning for instable points
       !
       if ( icistb > 0 ) then
          !
          rproc = 100.*real(icistb)/real(nfaces - nfacesb)
          !
          if ( .not. rproc < 1. ) then
             if ( l == lsal  ) write (msgstr,'(a,f5.1)') 'percentage of instable points for computing horizontal eddy diffusivity of salinity transport = ',rproc
             if ( l == ltemp ) write (msgstr,'(a,f5.1)') 'percentage of instable points for computing horizontal eddy diffusivity of heat transport = ',rproc
             if ( l == lsed  ) write (msgstr,'(a,f5.1)') 'percentage of instable points for computing horizontal eddy diffusivity of sediment transport = ',rproc
             call msgerr (1, trim(msgstr) )
          endif
          !
          if ( rproc > 30. ) then
             call msgerr ( 4, 'INSTABLE: unable to solve the transport equation!' )
             call msgerr ( 0, '          Please reduce the horizontal diffusivity coefficient!' )
             return
          endif
          !
       endif
       !
    enddo
    !
 101 format (' maximum level used in local time stepping = ',i2)
 102 format (2x,'the total mass associated with saline water is ',e14.8e2)
 103 format (2x,'the total mass associated with heat is ',e14.8e2)
 104 format (2x,'the total mass associated with suspended sediment ',e14.8e2)
    !
end subroutine SwashExpDepUtrans
