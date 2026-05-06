subroutine SwashExpLayUtrans
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
!   Performs the time integration for the layer-averaged transport equations on triangular mesh
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
!   The space discretization of the vertical advective and diffusivity terms is based on higher order
!   (flux-limited) schemes and central differences, respectively, in a finite volume fashion.
!
!   These vertical terms are treated semi-implicit. This results in a tri-diagonal system.
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
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Parameter variables
!
    real, parameter                               :: nuval = -999999. ! a special value that marks it as non-used
!
!   Local variables
!
    integer                                       :: btype            ! boundary type (see SwashUpdateUData.f90)
    integer                                       :: icell            ! cell index / loop counter over cells
    integer                                       :: icelll           ! left cell of present face
    integer                                       :: icellr           ! right cell of present face
    integer                                       :: icistb           ! counter for number of instable points
    integer, save                                 :: ient = 0         ! number of entries in this subroutine
    integer                                       :: iface            ! face index / loop counter over faces
    integer                                       :: j                ! loop counter
    integer                                       :: jc               ! loop counter over cells of vertex
    integer                                       :: jf               ! loop counter over faces (of cell)
    integer                                       :: jj               ! offset index in range of substeps
    integer                                       :: k                ! loop counter over vertical layers
    integer                                       :: kdd              ! index of layer k-2
    integer                                       :: ku               ! index of layer k+1
    integer                                       :: l                ! loop counter over constituents
    integer                                       :: m                ! loop counter over time-step levels
    integer                                       :: maxsb            ! maximum number of substeps
    integer                                       :: mtmax            ! maximum time-step level used in local time stepping
    integer                                       :: n                ! loop counter over substeps
    integer                                       :: nc               ! number of cells corresponding to specific time-step level
    integer                                       :: nsub             ! number of substeps for each time-step level
    integer, dimension(3)                         :: v                ! vertices of present cell
    integer                                       :: vf1              ! first vertex of present face
    integer                                       :: vf2              ! second vertex of present face
    integer                                       :: vu               ! upwind vertex
    !
    real                                          :: area             ! area of present cell
    real                                          :: areal            ! area of left cell of present face
    real                                          :: arear            ! area of right cell of present face
    real                                          :: bi               ! inverse of main diagonal of the matrix
    real                                          :: cfl              ! cell-based CFL number
    real                                          :: contrib          ! total contribution of transport flux per cell
    real                                          :: dhk              ! local increment in layer thickness
    real                                          :: dif2d            ! horizontal eddy diffusivity coefficient in velocity point
    real                                          :: fac              ! a factor
    real                                          :: fac1             ! another factor
    real                                          :: fac2             ! some other factor
    real                                          :: finp             ! interpolation factor
    real                                          :: fluxlim          ! flux limiter
    real                                          :: grad1            ! solution gradient
    real                                          :: grad2            ! another solution gradient
    real                                          :: hk0              ! layer thickness at previous substep
    real                                          :: hk1              ! layer thickness at current substep
    real                                          :: kwu              ! =1. if layer k+1 exists otherwise 0.
    real                                          :: lf               ! length of present face
    real                                          :: mass             ! total mass
    real                                          :: psm              ! Prandtl-Schmidt number
    real                                          :: qf               ! mass flux
    real                                          :: rdx              ! reciprocal of distance between circumcenters adjacent to face
    real                                          :: rnsb             ! reciprocal of number of substeps
    real                                          :: rproc            ! auxiliary variable with percentage of instable points
    real                                          :: rpu              ! averaged concentration in upwind vertex
    real                                          :: rsgn             ! sign for indicating face orientation
    real                                          :: rval             ! auxiliary real
    real                                          :: stabmx           ! auxiliary variable with maximum diffusivity based stability criterion
    real                                          :: sumqf            ! sum of outgoing mass fluxes per cell
    real                                          :: theta            ! implicitness factor for vertical terms
    real                                          :: totarea          ! total area of all cells around vertex
    real                                          :: w                ! local vertical velocity
    !
    character(120)                                :: msgstr           ! string to pass message
    !
    type(verttype), dimension(:), pointer         :: vert             ! datastructure for vertices with their attributes
    type(celltype), dimension(:), pointer         :: cell             ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer         :: face             ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashExpLayUtrans')
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
          do k = 1, kmax
             !
             qn(iface,k) = lf * hku(iface,k) * u1(iface,k)
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! implicitness factor for vertical terms
    !
    theta = pnums(33)
    !
    ! compute the time-step level for each cell
    !
    mlev = 0
    !
    do icell = 1, ncells
       !
       area = cell(icell)%attr(CELLAREA)
       !
       cflmax = -999.
       !
       do k = 1, kmax
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
             qf = rsgn * qn(iface,k)
             !
             if ( qf > 0. ) sumqf = sumqf + qf
             !
          enddo
          !
          ! compute the "flow" Courant number
          !
          if ( hso(icell) > epsdry ) then
             !
             cfl = sumqf * dt / hkso(icell,k) / area
             if ( cfl > cflmax ) cflmax = cfl
             !
          endif
          !
       enddo
       !
       ! note: CFL must not be larger than 0.5 to fulfill the max-min property
       if ( hso(icell) > epsdry ) mlev(icell) = max( 0, 1+floor( log(2.*cflmax)/log(2.) ) )
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
       rpo(:,:) = rp(:,:,l)
       !
       ! initialize system of equations in dry points
       !
       do icell = 1, ncells
          !
          if ( .not. hs(icell) > epsdry ) then
             !
             amatc(icell,:,1) = 1.
             amatc(icell,:,2) = 0.
             amatc(icell,:,3) = 0.
             rhsc (icell,:  ) = rpo(icell,:)
             !
          endif
          !
       enddo
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
             do k = 1, kmax
                !
                ! get mass flux at boundary face
                !
                qf = qn(iface,k)
                !
                if ( rsgn * qf > 0. ) then
                   !
                   ! outflow
                   !
                   fluxt(iface,k) = qf * rpo(icell,k)
                   !
                   bcrp(jf,k,l) = bcrpo(jf,k,l) - dt * rsgn * lf * u1(iface,k) * ( bcrpo(jf,k,l) - rpo(icell,k) ) / area
                   !
                   if ( l == lsal ) then
                      coutu (jf,k) = bcrp(jf,k,l)
                      icretu(jf,k) = tcret
                   endif
                   !
                else
                   !
                   ! inflow
                   !
                   fluxt(iface,k) = qf * bcrpo(jf,k,l)
                   !
                   if ( l == lsal ) then
                      fac = max(icretu(jf,k),0.) / max(tcret,dt)
                      bcrp(jf,k,l) = coutu(jf,k) + 0.5 * ( cbndu(jf,k,l) - coutu(jf,k) ) * ( 1. + cos(fac*pi) )
                      if ( .not. icretu(jf,k) < 0. ) icretu(jf,k) = icretu(jf,k) - dt
                   else
                      bcrp(jf,k,l) = cbndu(jf,k,l)
                   endif
                   !
                endif
                !
             enddo
             !
          else
             !
             ! --- boundary face closed: no advective flux
             !
             fluxt(iface,:)   = 0.
             bcrp (jf   ,:,l) = rp(icell,:,l)
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
                   do k = 1, kmax
                      !
                      ! compute mass flux at current face
                      !
                      qf = qn(iface,k)
                      !
                      ! compute the advective flux based on first order upwind
                      !
                      if ( qf > 0. ) then
                         !
                         fluxt(iface,k) = qf * rp0(icelll,k)
                         !
                      else
                         !
                         fluxt(iface,k) = qf * rp0(icellr,k)
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
                               rpu = rpu + area * rp0(icell,k)
                               !
                               totarea = totarea + area
                               !
                            enddo
                            !
                            rpu = rpu / totarea
                            !
                            ! compute solution gradients
                            !
                            grad1 = rp0(icellr,k) - rp0(icelll,k)
                            grad2 = rp0(icelll,k) - rpu
                            !
                            ! update flux
                            !
                            fluxt(iface,k) = fluxt(iface,k) + 0.5 * qf * fluxlim(grad1,grad2)
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
                               rpu = rpu + area * rp0(icell,k)
                               !
                               totarea = totarea + area
                               !
                            enddo
                            !
                            rpu = rpu / totarea
                            !
                            ! compute solution gradients
                            !
                            grad1 = rp0(icelll,k) - rp0(icellr,k)
                            grad2 = rp0(icellr,k) - rpu
                            !
                            ! update flux
                            !
                            fluxt(iface,k) = fluxt(iface,k) + 0.5 * qf * fluxlim(grad1,grad2)
                            !
                         endif
                         !
                      endif
                      !
                   enddo
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
                   do k = 1, kmax
                      !
                      fluxt(iface,k) = fluxt(iface,k) - ( finp*hkso(icelll,k) + (1.-finp)*hkso(icellr,k) ) * dif2d * lf * rdx * ( rp0(icellr,k) - rp0(icelll,k) )
                      !
                   enddo
                   !
                else if ( face(iface)%atti(FMARKER) == 0 ) then
                   !
                   fluxt(iface,:) = 0.
                   !
                endif
                !
             enddo
             !
             ! compute the time derivative
             !
             do j = 1, nc
                !
                ! consider cells with time-step level m
                !
                icell = mcell(j)
                !
                if ( hs(icell) > epsdry ) then
                   !
                   do k = 1, kmax
                      !
                      ! compute layer thicknesses at current and previous substeps
                      !
                      dhk = rnsb * ( hks(icell,k) - hkso(icell,k) )
                      !
                      hk0 = hkso(icell,k) + real(n-1) * dhk
                      hk1 = hk0           + dhk
                      !
                      amatc(icell,k,1) = hk1 / dtc(icell)
                      rhsc (icell,k  ) = hk0 * rp0(icell,k) / dtc(icell)
                   enddo
                   !
                endif
                !
             enddo
             !
             !  compute layer-averaged concentration in cells (based on finite volume approach)
             !
             do j = 1, nc
                !
                ! consider cells with time-step level m
                !
                icell = mcell(j)
                !
                do k = 1, kmax
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
                      contrib = contrib + rsgn * fluxt(iface,k)
                      !
                   enddo
                   !
                   ! update in wet cell
                   !
                   if ( hs(icell) > epsdry ) then
                      !
                      area = cell(icell)%attr(CELLAREA)
                      !
                      rhsc(icell,k) = rhsc(icell,k) - contrib / area
                      !
                   endif
                   !
                enddo
                !
             enddo
             !
             ! add mass exchange at bed in case of sediment transport in buoyancy flow
             !
             if ( l == lsed ) then
                !
                if ( psed(9) > 0. ) then
                   !
                   ! noncohesive sediment sand
                   ! (for meaning of different sediment parameters psed(), see module SwashCommdata3)
                   !
                   do j = 1, nc
                      !
                      ! consider cells with time-step level m
                      !
                      icell = mcell(j)
                      !
                      if ( hs(icell) > epsdry ) then
                         !
                         ! friction velocity computed from turbulent kinetic energy at bed
                         ! may include wave breaking-induced turbulence
                         !
                         fac = ( 0.3 + psed(7) ) * rtur(icell,kmax,1)
                         !
                         ! include upward sediment flux (erosion) by means of pickup function
                         !
                         if ( fac > psed(5) ) rhsc(icell,kmax) = rhsc(icell,kmax) + psed(9) * ( ( fac - psed(5) ) / psed(5) )**1.5
                         !
                         ! include downward sediment flux (deposition) by means of fall velocity
                         !
                         amatc(icell,kmax,1) = amatc(icell,kmax,1) + psed(1)
                         !
                      endif
                      !
                   enddo
                   !
                endif
                !
                if ( psed(10) > 0. ) then
                   !
                   ! cohesive sediment mud
                   ! (for meaning of different sediment parameters psed(), see module SwashCommdata3)
                   !
                   do j = 1, nc
                      !
                      ! consider cells with time-step level m
                      !
                      icell = mcell(j)
                      !
                      if ( hs(icell) > epsdry ) then
                         !
                         ! bed shear stress computed from turbulent kinetic energy at bed
                         !
                         fac = 0.3 * rhow * rtur(icell,kmax,1)
                         !
                         ! include upward sediment flux (erosion)
                         !
                         if( fac > psed(10) ) rhsc(icell,kmax) = rhsc(icell,kmax) + psed(12) * ( fac/psed(10) - 1. )
                         !
                         ! include downward sediment flux (deposition)
                         !
                         if( fac < psed(11) ) amatc(icell,kmax,1) = amatc(icell,kmax,1) + psed(1) * ( 1. - fac/psed(11) )
                         !
                      endif
                      !
                   enddo
                   !
                endif
                !
             endif
             !
             ! compute implicit part of vertical terms
             !
             propsc = nint(pnums(51))
             kappa  = pnums(52)
             mbound = pnums(53)
             phieby = pnums(54)
             !
             do j = 1, nc
                !
                ! consider cells with time-step level m
                !
                icell = mcell(j)
                !
                if ( hs(icell) > epsdry ) then
                   !
                   do k = 2, kmax
                      !
                      kdd = max(k-2,1   )
                      ku  = min(k+1,kmax)
                      !
                      kwu = 1.
                      if ( k == kmax ) kwu = 0.
                      !
                      w   = wom(icell,k-1)
                      fac = hks(icell,k-1) + hks(icell,k)
                      !
                      ! advection term
                      !
                      if ( propsc == 3 .and. kappa == 1. ) then
                         !
                         ! central differences
                         !
                         fac1 = theta * w * hks(icell,k-1) / fac
                         fac2 = theta * w * hks(icell,k  ) / fac
                         !
                         if ( k == 2 .or. k == kmax ) then
                            !
                            if ( w > 0. ) then
                               !
                               fac1 = theta * w
                               fac2 = 0.
                               !
                            else
                               !
                               fac1 = 0.
                               fac2 = theta * w
                               !
                            endif
                            !
                         endif
                         !
                         amatc(icell,k  ,1) = amatc(icell,k  ,1) + fac1
                         amatc(icell,k-1,1) = amatc(icell,k-1,1) - fac2
                         amatc(icell,k  ,2) =                      fac2
                         amatc(icell,k-1,3) =                    - fac1
                         !
                      else
                         !
                         ! first order upwind scheme
                         !
                         fac1 = theta * w
                         !
                         amatc(icell,k  ,1) = amatc(icell,k  ,1) + max(fac1,0.)
                         amatc(icell,k-1,1) = amatc(icell,k-1,1) - min(fac1,0.)
                         amatc(icell,k  ,2) =                      min(fac1,0.)
                         amatc(icell,k-1,3) =                    - max(fac1,0.)
                         !
                         ! add higher order (flux-limited) correction, if appropriate
                         !
                         if ( propsc /= 1 ) then
                            !
                            if ( w > 0. ) then
                               !
                               grad1 = rp0(icell,k-1) - rp0(icell,k )
                               grad2 = rp0(icell,k  ) - rp0(icell,ku)
                               !
                               fac2 = 0.5 * w * fluxlim(grad1,grad2)
                               !
                               rhsc(icell,k  ) = rhsc(icell,k  ) - fac2
                               rhsc(icell,k-1) = rhsc(icell,k-1) + fac2
                               !
                            else if ( w < 0. ) then
                               !
                               grad1 = rp0(icell,k  -1) - rp0(icell,k  )
                               grad2 = rp0(icell,kdd  ) - rp0(icell,k-1)
                               !
                               fac2 = 0.5 * w * fluxlim(grad1,grad2)
                               !
                               rhsc(icell,k  ) = rhsc(icell,k  ) + fac2
                               rhsc(icell,k-1) = rhsc(icell,k-1) - fac2
                               !
                            endif
                            !
                         endif
                         !
                      endif
                      !
                      ! diffusivity term
                      !
                      fac1 = 2. * theta * vnu3d(icell,k-1) / psm / fac
                      !
                      amatc(icell,k  ,1) = amatc(icell,k  ,1) + fac1
                      amatc(icell,k-1,1) = amatc(icell,k-1,1) + fac1
                      amatc(icell,k  ,2) = amatc(icell,k  ,2) - fac1
                      amatc(icell,k-1,3) = amatc(icell,k-1,3) - fac1
                      !
                      ! include fall velocity in case of sediment transport
                      !
                      if ( l == lsed .and. psed(1) > 0. ) then
                         !
                         amatc(icell,k-1,1) = amatc(icell,k-1,1) +     psed(1)
                         amatc(icell,k  ,2) = amatc(icell,k  ,2) - kwu*psed(1)
                         !
                      endif
                      !
                   enddo
                   !
                endif
                !
             enddo
             !
             ! compute explicit part of vertical terms, if appropriate
             !
             if ( theta /= 1. ) then
                !
                do j = 1, nc
                   !
                   ! consider cells with time-step level m
                   !
                   icell = mcell(j)
                   !
                   if ( hs(icell) > epsdry ) then
                      !
                      do k = 2, kmax
                         !
                         w   = wom (icell,k-1)
                         fac = hkso(icell,k-1) + hkso(icell,k)
                         !
                         ! advection term
                         !
                         if ( propsc == 3 .and. kappa == 1. ) then
                            !
                            ! central differences
                            !
                            fac1 = (1. - theta) * w * rp0(icell,k  ) * hkso(icell,k-1) / fac
                            fac2 = (1. - theta) * w * rp0(icell,k-1) * hkso(icell,k  ) / fac
                            !
                            if ( k == 2 .or. k == kmax ) then
                               !
                               if ( w > 0. ) then
                                  !
                                  fac1 = (1. - theta) * w * rp0(icell,k)
                                  fac2 = 0.
                                  !
                               else
                                  !
                                  fac1 = 0.
                                  fac2 = (1. - theta) * w * rp0(icell,k-1)
                                  !
                               endif
                               !
                            endif
                            !
                            rhsc(icell,k  ) = rhsc(icell,k  ) - fac1  - fac2
                            rhsc(icell,k-1) = rhsc(icell,k-1) + fac1  + fac2
                            !
                         else
                            !
                            ! first order upwind scheme
                            !
                            fac1 = (1. - theta) * max(w,0.) * rp0(icell,k  )
                            fac2 = (1. - theta) * min(w,0.) * rp0(icell,k-1)
                            !
                            rhsc(icell,k  ) = rhsc(icell,k  ) - fac1 - fac2
                            rhsc(icell,k-1) = rhsc(icell,k-1) + fac1 + fac2
                            !
                         endif
                         !
                         ! diffusivity term
                         !
                         fac1 = 2. * (1. - theta) * ( rp0(icell,k-1) - rp0(icell,k) ) * vnu3d(icell,k-1) / psm / fac
                         !
                         rhsc(icell,k  ) = rhsc(icell,k  ) + fac1
                         rhsc(icell,k-1) = rhsc(icell,k-1) - fac1
                         !
                      enddo
                      !
                   endif
                   !
                enddo
                !
             endif
             !
             ! solve the transport equation
             !
             do j = 1, nc
                !
                ! consider cells with time-step level m
                !
                icell = mcell(j)
                !
                bi = 1./amatc(icell,1,1)
                !
                amatc(icell,1,1) = bi
                amatc(icell,1,3) = amatc(icell,1,3)*bi
                rhsc (icell,1  ) = rhsc (icell,1  )*bi
                !
                do k = 2, kmax
                   !
                   bi = 1./(amatc(icell,k,1) - amatc(icell,k,2)*amatc(icell,k-1,3))
                   amatc(icell,k,1) = bi
                   amatc(icell,k,3) = amatc(icell,k,3)*bi
                   rhsc (icell,k  ) = (rhsc(icell,k) - amatc(icell,k,2)*rhsc(icell,k-1))*bi
                   !
                enddo
                !
                rp1(icell,kmax) = rhsc(icell,kmax)
                do k = kmax-1, 1, -1
                   rp1(icell,k) = rhsc(icell,k) - amatc(icell,k,3)*rp1(icell,k+1)
                enddo
                !
             enddo
             !
             ! also use intermediate solution of coarser time-step levels
             !
             do k = 1, kmax
                !
                do icell = 1, ncells
                   !
                   if ( mlev(icell) < m ) then
                      !
                      rval = rpi(icell,k,n*maxsb)
                      if (rval > nuval ) rp1(icell,k) = rval
                      !
                   endif
                   !
                enddo
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
                do k = 1, kmax
                   rpi(icell,k,jj+(n-1)*maxsb) = 0.5 * ( rp0(icell,k) + rp1(icell,k) )
                enddo
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
             do k = 1, kmax
                !
                rp(icell,k,l) = rp1(icell,k)
                !
             enddo
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
             do k = 1, kmax
                !
                mass = mass + hks(icell,k) * area * rp(icell,k,l)
                !
             enddo
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
end subroutine SwashExpLayUtrans
