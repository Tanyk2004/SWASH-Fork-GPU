subroutine SwashKepsMod2DH
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
!    1.00,  February 2011: New subroutine
!    1.05,   January 2012: extension logarithmic wall-law
!    1.15,   October 2012: extension repeating grid
!    1.21,       May 2013: buoyancy term included
!    1.22, September 2013: vegetation term included
!
!   Purpose
!
!   Builds and solves standard k-eps model for vertical mixing in two dimensions
!
!   Method
!
!   The time integration is fully implicit. First order timesplitting is used.
!
!   In the first stage the horizontal advection terms are integrated implicitly.
!   In space, they are discretized with first order upwinding.
!   The systems of equations are decoupled in the vertical.
!   A four-direction Gauss seidel iterative solution method is employed in which
!   four sweeps, one for each distinct direction depending on the flow velocity,
!   are carried out per iteration.
!
!   In the second stage the vertical transport terms are integrated implicitly in
!   conjunction with the source and sink terms. Sinks are treated implicitly, while
!   sources are taken explicitly.
!   The vertical advection term is approximated by means of first order upwinding.
!   The systems of equations are decoupled in the horizontal. In the vertical,
!   the tri-diagonal systems of equations are solved by a double sweep.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
    use m_genarr, only: kgrpnt, guu, gvv
    use m_parall
    use SwashFlowdata
!
    implicit none
!
!   Parameter variables
!
    real,    parameter :: eps    = 1.e-7 ! small number
    real,    parameter :: hmin   = 1.e-3 ! minimal water depth in meters
    !
    real,    parameter :: ceps1  = 1.44  ! closure constant for standard k-eps model
    real,    parameter :: ceps2  = 1.92  ! closure constant for standard k-eps model
    real,    parameter :: cmu    = 0.09  ! closure constant for standard k-eps model
    real,    parameter :: sigmak = 1.0   ! closure constant for standard k-eps model
    real,    parameter :: sigmae = 1.3   ! closure constant for standard k-eps model
    real,    parameter :: sigrho = 0.5   ! closure constant for buoyancy term
    !
    integer, parameter :: maxit = 100    ! maximum number of iterations
    real,    parameter :: reps  = 0.5e-4 ! accuracy of the final approximation
    !
    real,    parameter :: gamma = 0.34   ! empirical coefficient for added mass
!
!   Local variables
!
    integer, save      :: ient = 0       ! number of entries in this subroutine
    integer            :: j              ! iteration counter
    integer            :: k              ! loop counter over vertical layer interfaces
    integer            :: kend           ! end index of loop over vertical layer interfaces
    integer            :: kinf           ! k-index of layer interface with largest error in solution
    integer            :: ku             ! index of layer k+1
    integer            :: l              ! loop variable
    integer            :: m              ! loop counter
    integer            :: md             ! index of point m-1
    integer            :: mend           ! end index of loop over wl-points in x-direction
    integer            :: minf           ! m-index of point with largest error in solution
    integer            :: msta           ! start index of loop over wl-points in x-direction
    integer            :: mu             ! index of point m+1
    integer            :: n              ! loop counter
    integer            :: nd             ! index of point n-1
    integer            :: ndm            ! pointer to m,n-1
    integer            :: nend           ! end index of loop over wl-points in y-direction
    integer            :: nfm            ! pointer to m,nf
    integer            :: ninf           ! n-index of point with largest error in solution
    integer            :: nlum           ! pointer to m,nlu
    integer            :: nm             ! pointer to m,n
    integer            :: nmd            ! pointer to m-1,n
    integer            :: nmf            ! pointer to mf,n
    integer            :: nmlu           ! pointer to mlu,n
    integer            :: nmu            ! pointer to m+1,n
    integer            :: nsta           ! start index of loop over wl-points in y-direction
    integer            :: nu             ! index of point n+1
    integer            :: num            ! pointer to m,n+1
    !
    real               :: alpha          ! dimensionless constant for laminar friction factor
    real               :: bi             ! inverse of main diagonal of the matrix
    real               :: buotke         ! local buoyancy term representing exchange between turbulent energy and potential energy
    real               :: cfe            ! vegetative drag-related closure constant for dissipation rate
    real               :: cfk            ! vegetative drag-related closure constant for turbulent kinetic energy
    real               :: ctrkb          ! contribution of vertical advection term below considered point
    real               :: ctrkt          ! contribution of vertical advection term above considered point
    real               :: ctrml          ! contribution of horizontal advection term left from considered point
    real               :: ctrmr          ! contribution of horizontal advection term right from considered point
    real               :: ctrnd          ! contribution of horizontal advection term rear from considered point
    real               :: ctrnu          ! contribution of horizontal advection term front from considered point
    real               :: dinv           ! inverse of grain size
    real               :: drhodz         ! vertical gradient of density
    real               :: dudz           ! vertical gradient of u-velocity
    real               :: dvdz           ! vertical gradient of v-velocity
    real               :: dxl            ! local mesh size in x-direction
    real               :: dyl            ! local mesh size in y-direction
    real               :: dz             ! local layer thickness
    real               :: epsinf         ! intrinsic volume average dissipation rate in a
                                         ! periodically fully-developed flow through porous medium
    real               :: epslin         ! required accuracy in the linear solver
    real               :: fac            ! auxiliary factor
    real               :: tkeinf         ! intrinsic volume average turbulent kinetic energy in a
                                         ! periodically fully-developed flow through porous medium
    real               :: prodveg        ! turbulence production due to vegetation
    real               :: psm            ! Prandtl-Schmidt number
    real               :: res            ! residual
    real               :: resm           ! maximum error
    real               :: rturi          ! intermediate turbulence quantity
    real               :: rturmx         ! maximum value of turbulence quantity
    real               :: rval           ! auxiliary real
    real               :: u              ! local u-velocity
    real               :: ustar2         ! friction velocity squared
    real               :: utot           ! velocity magnitude
    real               :: v              ! local v-velocity
    real               :: w              ! local relative vertical velocity
    real               :: zs             ! distance of half layer thickness to bottom / free surface in waterlevel point
    !
    logical            :: EQREAL         ! compares two reals
    logical            :: STPNOW         ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashKepsMod2DH')
    !
    ! first part of time splitting method (horizontal advection terms)
    !
    ! initialize system of equations in dry points
    !
    do n = nf, nlu
       do m = mf, mlu
          !
          nm = kgrpnt(m,n)
          !
          if ( .not. hs(nm) > hmin ) then
             !
             amatt(nm,:,1) = 1.
             amatt(nm,:,2) = 0.
             amatt(nm,:,3) = 0.
             amatt(nm,:,4) = 0.
             amatt(nm,:,5) = 0.
             rhst (nm,:  ) = eps
             !
          endif
          !
       enddo
    enddo
    !
    ! interior part
    !
    if ( irough == 4 ) then
       !
       ! in case of logarithmic wall-law an equation for turbulent kinetic energy is build up at bottom
       !
       kend = kmax
    else
       kend = kmax-1
    endif
    !
    ! build time derivative and horizontal transport (implicit)
    !
    do k = 1, kend
       !
       ku = min(k+1,kmax)
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmd == 1 ) nmd = nm
             if ( ndm == 1 ) ndm = nm
             !
             if ( hs(nm) > hmin ) then
                !
                dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
                dyl = 0.5 * ( guu(nm) + guu(nmd) )
                !
                ! time derivative
                !
                amatt(nm,k,1) = 1. / dt
                !
                ! horizontal advective terms (first order upwinding)
                !
                u = 0.25 * ( u1(nm,k) + u1(nmd,k) + u1(nm,ku) + u1(nmd,ku) )
                v = 0.25 * ( v1(nm,k) + v1(ndm,k) + v1(nm,ku) + v1(ndm,ku) )
                !
                ctrml = 0.5 * (u + abs(u)) / dxl
                ctrmr = 0.5 * (u - abs(u)) / dxl
                ctrnd = 0.5 * (v + abs(v)) / dyl
                ctrnu = 0.5 * (v - abs(v)) / dyl
                !
                amatt(nm,k,1) = amatt(nm,k,1) + ctrml - ctrmr + ctrnd - ctrnu
                amatt(nm,k,2) = -ctrml
                amatt(nm,k,3) =  ctrmr
                amatt(nm,k,4) = -ctrnd
                amatt(nm,k,5) =  ctrnu
                !
             endif
             !
          enddo
       enddo
       !
    enddo
    !
    ! impose homogeneous Neumann condition at boundaries, if appropriate
    !
    if ( .not.lreptx ) then
       !
       do n = nfu, nl
          !
          nmf  = kgrpnt(mf ,n)
          nmlu = kgrpnt(mlu,n)
          !
          if ( LMXF ) then
             do k = 1, kend
                amatt(nmf,k,1) =  1.0
                amatt(nmf,k,3) = -1.0
             enddo
          endif
          !
          if ( LMXL ) then
             do k = 1, kend
                amatt(nmlu,k,1) =  1.0
                amatt(nmlu,k,2) = -1.0
             enddo
          endif
          !
       enddo
       !
    endif
    !
    if ( .not.lrepty ) then
       !
       do m = mfu, ml
          !
          nfm  = kgrpnt(m,nf )
          nlum = kgrpnt(m,nlu)
          !
          if ( LMYF ) then
             do k = 1, kend
                amatt(nfm,k,1) =  1.0
                amatt(nfm,k,5) = -1.0
             enddo
          endif
          !
          if ( LMYL ) then
             do k = 1, kend
                amatt(nlum,k,1) =  1.0
                amatt(nlum,k,4) = -1.0
             enddo
          endif
          !
       enddo
       !
    endif
    !
    ! In case of a subdomain interface not coinciding with a real boundary, the begin point and/or
    ! end point of the loop becomes mfu and ml, respectively, as there is no boundary condition.
    ! When the grid is repeated in x-direction, the loop starts at mfu and ends at ml.
    ! The same for y-direction.
    !
    msta = mf                               ! first wl-point in x-direction
    if ( .not.LMXF .or. lreptx ) msta = mfu
    mend = mlu                              ! last  wl-point in x-direction
    if ( .not.LMXL .or. lreptx ) mend = ml
    nsta = nf                               ! first wl-point in y-direction
    if ( .not.LMYF .or. lrepty ) nsta = nfu
    nend = nlu                              ! last  wl-point in y-direction
    if ( .not.LMYL .or. lrepty ) nend = nl
    !
    ! solve each turbulence equation by means of four-direction Gauss-Seidel iteration
    !
    do l = 1, ltur
       !
       ! build right-hand side first
       !
       do k = 1, kend
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > hmin ) then
                   !
                   rhst(nm,k) = rtur(nm,k,l) / dt
                   !
                else
                   !
                   rhst(nm,k) = rtur(nm,k,l)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! determine maximum of turbulence quantity and required accuracy
       !
       rturmx = maxval(abs(rtur(:,:,l)))
       epslin = max(reps,reps*rturmx)
       call SWREDUCE ( epslin, 1, SWREAL, SWMAX )
       !
       j    = 0
       resm = epslin + 1.
       !
 10    if ( resm > epslin .and. j < maxit ) then
          !
          j = j + 1
          !
          resm = 0.
          minf = 0
          ninf = 0
          kinf = 0
          !
          do k = 1, kend
             !
             ku = min(k+1,kmax)
             !
             do n = nsta, nend
                do m = msta, mend
                   !
                   md = m - 1
                   nd = n - 1
                   if ( LMXF .and. md < mf ) md = mf
                   if ( LMYF .and. nd < nf ) nd = nf
                   !
                   nm  = kgrpnt(m ,n )
                   nmd = kgrpnt(md,n )
                   ndm = kgrpnt(m ,nd)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd == 1 ) nmd = nm
                   if ( ndm == 1 ) ndm = nm
                   !
                   u = 0.25 * ( u1(nm,k) + u1(nmd,k) + u1(nm,ku) + u1(nmd,ku) )
                   v = 0.25 * ( v1(nm,k) + v1(ndm,k) + v1(nm,ku) + v1(ndm,ku) )
                   !
                   if ( .not. u < 0. .and. .not. v < 0. ) then
                      !
                      rturi = rhst(nm,k) - amatt(nm,k,2)*rtur(nmd,k,l) - amatt(nm,k,4)*rtur(ndm,k,l)
                      rturi = rturi / amatt(nm,k,1)
                      !
                      res = abs(rtur(nm,k,l) - rturi)
                      if ( res > resm ) then
                         resm = res
                         minf = m
                         ninf = n
                         kinf = k
                      endif
                      rtur(nm,k,l) = rturi
                      !
                   endif
                   !
                enddo
             enddo
             !
             do n = nsta, nend
                do m = mend, msta, -1
                   !
                   md = m - 1
                   mu = m + 1
                   nd = n - 1
                   if ( LMXF .and. md < mf  ) md = mf
                   if ( LMXL .and. mu > mlu ) mu = mlu
                   if ( LMYF .and. nd < nf  ) nd = nf
                   !
                   nm  = kgrpnt(m ,n )
                   nmd = kgrpnt(md,n )
                   nmu = kgrpnt(mu,n )
                   ndm = kgrpnt(m ,nd)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd == 1 ) nmd = nm
                   if ( nmu == 1 ) nmu = nm
                   if ( ndm == 1 ) ndm = nm
                   !
                   u = 0.25 * ( u1(nm,k) + u1(nmd,k) + u1(nm,ku) + u1(nmd,ku) )
                   v = 0.25 * ( v1(nm,k) + v1(ndm,k) + v1(nm,ku) + v1(ndm,ku) )
                   !
                   if ( u < 0. .and. .not. v < 0. ) then
                      !
                      rturi = rhst(nm,k) - amatt(nm,k,3)*rtur(nmu,k,l) - amatt(nm,k,4)*rtur(ndm,k,l)
                      rturi = rturi / amatt(nm,k,1)
                      !
                      res = abs(rtur(nm,k,l) - rturi)
                      if ( res > resm ) then
                         resm = res
                         minf = m
                         ninf = n
                         kinf = k
                      endif
                      rtur(nm,k,l) = rturi
                      !
                   endif
                   !
                enddo
             enddo
             !
             do n = nend, nsta, -1
                do m = msta, mend
                   !
                   md = m - 1
                   nd = n - 1
                   nu = n + 1
                   if ( LMXF .and. md < mf  ) md = mf
                   if ( LMYF .and. nd < nf  ) nd = nf
                   if ( LMYL .and. nu > nlu ) nu = nlu
                   !
                   nm  = kgrpnt(m ,n )
                   nmd = kgrpnt(md,n )
                   ndm = kgrpnt(m ,nd)
                   num = kgrpnt(m ,nu)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd == 1 ) nmd = nm
                   if ( ndm == 1 ) ndm = nm
                   if ( num == 1 ) num = nm
                   !
                   u = 0.25 * ( u1(nm,k) + u1(nmd,k) + u1(nm,ku) + u1(nmd,ku) )
                   v = 0.25 * ( v1(nm,k) + v1(ndm,k) + v1(nm,ku) + v1(ndm,ku) )
                   !
                   if ( .not. u < 0. .and. v < 0. ) then
                      !
                      rturi = rhst(nm,k) - amatt(nm,k,2)*rtur(nmd,k,l) - amatt(nm,k,5)*rtur(num,k,l)
                      rturi = rturi / amatt(nm,k,1)
                      !
                      res = abs(rtur(nm,k,l) - rturi)
                      if ( res > resm ) then
                         resm = res
                         minf = m
                         ninf = n
                         kinf = k
                      endif
                      rtur(nm,k,l) = rturi
                      !
                   endif
                   !
                enddo
             enddo
             !
             do n = nend, nsta, -1
                do m = mend, msta, -1
                   !
                   md = m - 1
                   mu = m + 1
                   nd = n - 1
                   nu = n + 1
                   if ( LMXF .and. md < mf  ) md = mf
                   if ( LMXL .and. mu > mlu ) mu = mlu
                   if ( LMYF .and. nd < nf  ) nd = nf
                   if ( LMYL .and. nu > nlu ) nu = nlu
                   !
                   nm  = kgrpnt(m ,n )
                   nmd = kgrpnt(md,n )
                   nmu = kgrpnt(mu,n )
                   ndm = kgrpnt(m ,nd)
                   num = kgrpnt(m ,nu)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd == 1 ) nmd = nm
                   if ( nmu == 1 ) nmu = nm
                   if ( ndm == 1 ) ndm = nm
                   if ( num == 1 ) num = nm
                   !
                   u = 0.25 * ( u1(nm,k) + u1(nmd,k) + u1(nm,ku) + u1(nmd,ku) )
                   v = 0.25 * ( v1(nm,k) + v1(ndm,k) + v1(nm,ku) + v1(ndm,ku) )
                   !
                   if ( u < 0. .and. v < 0. ) then
                      !
                      rturi = rhst(nm,k) - amatt(nm,k,3)*rtur(nmu,k,l) - amatt(nm,k,5)*rtur(num,k,l)
                      rturi = rturi / amatt(nm,k,1)
                      !
                      res = abs(rtur(nm,k,l) - rturi)
                      if ( res > resm ) then
                         resm = res
                         minf = m
                         ninf = n
                         kinf = k
                      endif
                      rtur(nm,k,l) = rturi
                      !
                   endif
                   !
                enddo
             enddo
             !
          enddo
          !
          rval = resm
          call SWREDUCE ( resm, 1, SWREAL, SWMAX )
          !
          if ( ITEST >= 30 .and. EQREAL(resm,rval) ) then
             !
             write (PRINTF,'(a,i3,a,e12.6,a,i1,a,i4,a,i4,a,i3)') ' ++ keps: iter = ',j,' res = ',resm,' of l=',l,' in (m,n)=(',minf+MXF-2,',',ninf+MYF-2,'), layer=',kinf
             !
          endif
          !
          ! exchange turbulence quantity with neighbouring subdomains, if appropriate
          !
          call SWEXCHG ( rtur(1,0,l), kgrpnt, 0, kmax )
          if (STPNOW()) return
          !
          ! synchronize turbulence quantity at appropriate boundaries in case of repeating grid, if appropriate
          !
          call periodic ( rtur(1,0,l), kgrpnt, 0, kmax )
          !
          goto 10
          !
       endif
       !
    enddo
    !
    ! second part of time splitting method (vertical advection and diffusion, production, buoyancy, vegetation and dissipation rates)
    !
    cfk = pturb(2)
    cfe = pturb(3)
    !
    ! compute magnitude of shear squared in case of vertical mixing
    !
    do k = 1, kmax-1
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmd == 1 ) nmd = nm
             if ( ndm == 1 ) ndm = nm
             !
             if ( hs(nm) > hmin ) then
                !
                dudz = ( u1(nm,k) + u1(nmd,k) - u1(nm,k+1) - u1(nmd,k+1) ) / ( hks(nm,k) + hks(nm,k+1) )
                !
                dvdz = ( v1(nm,k) + v1(ndm,k) - v1(nm,k+1) - v1(ndm,k+1) ) / ( hks(nm,k) + hks(nm,k+1) )
                !
                zshear(nm,k) = dudz*dudz + dvdz*dvdz
                !
             endif
             !
          enddo
       enddo
       !
    enddo
    !
    do l = 1, ltur
       !
       if ( l == 1 ) psm = sigmak
       if ( l == 2 ) psm = sigmae
       !
       if ( irough == 4 .and. l == 1 ) then
          kend = kmax
       else
          kend = kmax-1
       endif
       !
       ! build time-derivative
       !
       do k = 1, kend
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > hmin ) then
                   !
                   ! time derivative
                   !
                   amatt(nm,k,1) = 1. / dt
                   amatt(nm,k,2) = 0.
                   amatt(nm,k,3) = 0.
                   rhst (nm,k  ) = rtur(nm,k,l) / dt
                   !
                else
                   !
                   rhst(nm,k) = rtur(nm,k,l)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! build vertical transport (implicit)
       !
       do k = 1, kmax-1
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > hmin ) then
                   !
                   ! vertical advective term (first order upwinding)
                   !
                   w = wom(nm,k)
                   !
                   ctrkb = 0.5 * (w + abs(w)) / hks(nm,k+1)
                   ctrkt = 0.5 * (w - abs(w)) / hks(nm,k  )
                   !
                   if ( k == kmax-1 ) ctrkb = 0.
                   !
                   amatt(nm,k,1) = amatt(nm,k,1) + ctrkb - ctrkt
                   amatt(nm,k,2) =  ctrkt
                   amatt(nm,k,3) = -ctrkb
                   !
                   ! vertical diffusive term
                   !
                   dz = 0.5 * ( hks(nm,k) + hks(nm,k+1) )
                   !
                   ctrkb = 0.5 * ( vnu3d(nm,k  ) + vnu3d(nm,k+1) ) / ( psm * dz * hks(nm,k+1) )
                   ctrkt = 0.5 * ( vnu3d(nm,k-1) + vnu3d(nm,k  ) ) / ( psm * dz * hks(nm,k  ) )
                   !
                   amatt(nm,k,1) = amatt(nm,k,1) + ctrkb + ctrkt
                   amatt(nm,k,2) = amatt(nm,k,2) - ctrkt
                   amatt(nm,k,3) = amatt(nm,k,3) - ctrkb
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! build source and sink terms
       !
       if ( l == 1 ) then
          !
          do k = 1, kmax-1
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( hs(nm) > hmin ) then
                      !
                      ! production rate term
                      !
                      rhst(nm,k) = rhst(nm,k) + vnu3d(nm,k) * zshear(nm,k)
                      !
                      ! dissipation rate term (Newton linearization)
                      !
                      amatt(nm,k,1) = amatt(nm,k,1) + 2. * rtur(nm,k,2) / ( rtur(nm,k,1)+eps )
                      rhst (nm,k  ) = rhst (nm,k  ) + rtur(nm,k,2)
                      !
                   endif
                   !
                enddo
             enddo
             !
          enddo
          !
          if ( idens /= 0 ) then
             !
             ! buoyancy term
             !
             do k = 1, kmax-1
                !
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      nm = kgrpnt(m,n)
                      !
                      if ( hs(nm) > hmin ) then
                         !
                         drhodz = 2. * ( rho(nm,k) - rho(nm,k+1) ) / ( hks(nm,k) + hks(nm,k+1) )
                         !
                         buotke = vnu3d(nm,k) * grav * drhodz / ( rhow*sigrho )
                         !
                         if ( buotke > 0. ) then
                            !
                            rhst(nm,k) = rhst(nm,k) + buotke
                            !
                         else
                            !
                            amatt(nm,k,1) = amatt(nm,k,1) - buotke / ( rtur(nm,k,1)+eps )
                            !
                         endif
                         !
                      endif
                      !
                   enddo
                enddo
                !
             enddo
             !
          endif
          !
          if ( iveg /= 0 ) then
             !
             ! vegetation term
             !
             do k = 1, kmax-1
                !
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      md = m - 1
                      nd = n - 1
                      !
                      nm  = kgrpnt(m ,n )
                      nmd = kgrpnt(md,n )
                      ndm = kgrpnt(m ,nd)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( nmd == 1 ) nmd = nm
                      if ( ndm == 1 ) ndm = nm
                      !
                      if ( hs(nm) > hmin ) then
                         !
                         u = 0.25 * ( u1(nm,k) + u1(nmd,k) + u1(nm,k+1) + u1(nmd,k+1) )
                         v = 0.25 * ( v1(nm,k) + v1(ndm,k) + v1(nm,k+1) + v1(ndm,k+1) )
                         !
                         utot = sqrt( u*u + v*v )
                         !
                         prodveg = 0.25 * utot * ( cvegu(nm,k  ,1) * u1(nm,k  )**2 + cvegu(nmd,k  ,1) * u1(nmd,k  )**2 + &
                                                   cvegu(nm,k+1,1) * u1(nm,k+1)**2 + cvegu(nmd,k+1,1) * u1(nmd,k+1)**2 + &
                                                   cvegv(nm,k  ,1) * v1(nm,k  )**2 + cvegv(ndm,k  ,1) * v1(ndm,k  )**2 + &
                                                   cvegv(nm,k+1,1) * v1(nm,k+1)**2 + cvegv(ndm,k+1,1) * v1(ndm,k+1)**2 )
                         !
                         rhst(nm,k) = rhst(nm,k) + cfk * prodveg
                         !
                      endif
                      !
                   enddo
                enddo
                !
             enddo
             !
          endif
          !
          if ( iporos == 1 ) then
             !
             ! closure term for porous medium flow
             !
             alpha = ppor(3)
             fac   = gamma / sqrt(kinvis*alpha)
             !
             do k = 1, kmax-1
                !
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      md = m - 1
                      nd = n - 1
                      !
                      nm  = kgrpnt(m ,n )
                      nmd = kgrpnt(md,n )
                      ndm = kgrpnt(m ,nd)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( nmd == 1 ) nmd = nm
                      if ( ndm == 1 ) ndm = nm
                      !
                      if ( npoks(nm,k) > 0.1 .and. npoks(nm,k) < 0.99 ) then
                         !
                         u = 0.25 * ( u1(nm,k) + u1(nmd,k) + u1(nm,k+1) + u1(nmd,k+1) )
                         v = 0.25 * ( v1(nm,k) + v1(ndm,k) + v1(nm,k+1) + v1(ndm,k+1) )
                         !
                         utot = sqrt( u*u + v*v + w1(nm,k)*w1(nm,k) )
                         !
                         ! this factor is the reciprocal of grain size
                         dinv = fac * sqrt(apoks(nm,k)) / cpoks(nm,k)
                         !
                         epsinf = 39. * ( 1. - npoks(nm,k) )**2.5 * npoks(nm,k)**2 * utot**3 * dinv
                         !
                         rhst(nm,k) = rhst(nm,k) + epsinf
                         !
                      endif
                      !
                   enddo
                enddo
                !
             enddo
             !
          endif
          !
       else if ( l == 2 ) then
          !
          do k = 1, kmax-1
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( hs(nm) > hmin ) then
                      !
                      ! production rate term
                      !
                      rhst(nm,k) = rhst(nm,k) + ceps1 * cmu * rtur(nm,k,1) * zshear(nm,k)
                      !
                      ! dissipation rate term (Newton linearization)
                      !
                      fac = ceps2 * rtur(nm,k,2) / ( rtur(nm,k,1)+eps )
                      !
                      amatt(nm,k,1) = amatt(nm,k,1) + 2. * fac
                      rhst (nm,k  ) = rhst (nm,k  ) + fac * rtur(nm,k,2)
                      !
                   endif
                   !
                enddo
             enddo
             !
          enddo
          !
          if ( idens /= 0 ) then
             !
             ! buoyancy term
             !
             do k = 1, kmax-1
                !
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      nm = kgrpnt(m,n)
                      !
                      if ( hs(nm) > hmin ) then
                         !
                         drhodz = 2. * ( rho(nm,k) - rho(nm,k+1) ) / ( hks(nm,k) + hks(nm,k+1) )
                         !
                         if ( drhodz > 0. ) then
                            !
                            rhst(nm,k) = rhst(nm,k) + ceps1 * cmu * rtur(nm,k,1) * grav * drhodz / ( rhow*sigrho )
                            !
                         endif
                         !
                      endif
                      !
                   enddo
                enddo
                !
             enddo
             !
          endif
          !
          if ( iveg /= 0 ) then
             !
             ! vegetation term
             !
             do k = 1, kmax-1
                !
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      md = m - 1
                      nd = n - 1
                      !
                      nm  = kgrpnt(m ,n )
                      nmd = kgrpnt(md,n )
                      ndm = kgrpnt(m ,nd)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( nmd == 1 ) nmd = nm
                      if ( ndm == 1 ) ndm = nm
                      !
                      if ( hs(nm) > hmin ) then
                         !
                         u = 0.25 * ( u1(nm,k) + u1(nmd,k) + u1(nm,k+1) + u1(nmd,k+1) )
                         v = 0.25 * ( v1(nm,k) + v1(ndm,k) + v1(nm,k+1) + v1(ndm,k+1) )
                         !
                         utot = sqrt( u*u + v*v )
                         !
                         prodveg = 0.25 * utot * ( cvegu(nm,k  ,1) * u1(nm,k  )**2 + cvegu(nmd,k  ,1) * u1(nmd,k  )**2 + &
                                                   cvegu(nm,k+1,1) * u1(nm,k+1)**2 + cvegu(nmd,k+1,1) * u1(nmd,k+1)**2 + &
                                                   cvegv(nm,k  ,1) * v1(nm,k  )**2 + cvegv(ndm,k  ,1) * v1(ndm,k  )**2 + &
                                                   cvegv(nm,k+1,1) * v1(nm,k+1)**2 + cvegv(ndm,k+1,1) * v1(ndm,k+1)**2 )
                         !
                         rhst(nm,k) = rhst(nm,k) + ceps1 * cfe * prodveg * rtur(nm,k,2) / ( rtur(nm,k,1)+eps )
                         !
                      endif
                      !
                   enddo
                enddo
                !
             enddo
             !
          endif
          !
          if ( iporos == 1 ) then
             !
             ! closure term for porous medium flow
             !
             alpha = ppor(3)
             fac   = gamma / sqrt(kinvis*alpha)
             !
             do k = 1, kmax-1
                !
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      md = m - 1
                      nd = n - 1
                      !
                      nm  = kgrpnt(m ,n )
                      nmd = kgrpnt(md,n )
                      ndm = kgrpnt(m ,nd)
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( nmd == 1 ) nmd = nm
                      if ( ndm == 1 ) ndm = nm
                      !
                      if ( npoks(nm,k) > 0.1 .and. npoks(nm,k) < 0.99 ) then
                         !
                         u = 0.25 * ( u1(nm,k) + u1(nmd,k) + u1(nm,k+1) + u1(nmd,k+1) )
                         v = 0.25 * ( v1(nm,k) + v1(ndm,k) + v1(nm,k+1) + v1(ndm,k+1) )
                         !
                         utot = sqrt( u*u + v*v + w1(nm,k)*w1(nm,k) )
                         !
                         ! this factor is the reciprocal of grain size
                         dinv = fac * sqrt(apoks(nm,k)) / cpoks(nm,k)
                         !
                         tkeinf = 3.7 * ( 1. - npoks(nm,k) )      * npoks(nm,k)**1.5 * utot**2
                         epsinf = 39. * ( 1. - npoks(nm,k) )**2.5 * npoks(nm,k)**2   * utot**3 * dinv
                         !
                         rhst(nm,k) = rhst(nm,k) + ceps2 * epsinf**2 / ( tkeinf+eps )
                         !
                      endif
                      !
                   enddo
                enddo
                !
             enddo
             !
          endif
          !
       endif
       !
       ! free surface
       !
       if ( l == 1 ) then
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > hmin ) then
                   !
                   amatt(nm,0,1) =  1. + rtur(nm,2,1) / ( rtur(nm,0,1)+eps )
                   amatt(nm,0,2) =  0.
                   amatt(nm,0,3) = -2.
                   rhst (nm,0  ) =  0.
                   !
                else
                   !
                   rhst(nm,0) = rtur(nm,0,1)
                   !
                endif
                !
             enddo
          enddo
          !
       else if ( l == 2 ) then
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > hmin ) then
                   !
                   zs = 0.5 * ( zks(nm,0) - zks(nm,1) )
                   !
                   amatt(nm,0,1) = 1.
                   rhst (nm,0  ) = (cmu**.75) * (rtur(nm,0,1)**1.5) / ( vonkar*zs )
                   !
                else
                   !
                   rhst(nm,0) = rtur(nm,0,2)
                   !
                endif
                !
             enddo
          enddo
          !
       endif
       !
       ! bottom
       !
       if ( l == 1 ) then
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd == 1 ) nmd = nm
                if ( ndm == 1 ) ndm = nm
                !
                if ( hs(nm) > hmin ) then
                   !
                   ! generation of turbulence due to bottom friction
                   !
                   if ( irough /= 4 ) then
                      !
                      ustar2 = 0.25 * (cfricu(nm)+cfricu(nmd)) * (udep(nm)*udep(nm) + udep(nmd)*udep(nmd)) +  &
                               0.25 * (cfricv(nm)+cfricv(ndm)) * (vdep(nm)*vdep(nm) + vdep(ndm)*vdep(ndm))
                      !
                      amatt(nm,kmax,1) = 1.
                      rhst (nm,kmax  ) = ustar2 / sqrt(cmu)
                      !
                   else
                      !
                      zs = 0.5 * ( zks(nm,kmax-1) - zks(nm,kmax) )
                      !
                      ! production rate term
                      !
                      rhst(nm,kmax) = rhst(nm,kmax) + 0.25 * logfrc(nm,1) * ((u1(nm,kmax) + u1(nmd,kmax))**2 + (v1(nm,kmax) + v1(ndm,kmax))**2) / zs
                      !
                      ! dissipation rate term
                      !
                      amatt(nm,kmax,1) = amatt(nm,kmax,1) + (cmu**.75) * sqrt(rtur(nm,kmax,1)) * logfrc(nm,2) / zs
                      !
                   endif
                   !
                else
                   !
                   rhst(nm,kmax) = rtur(nm,kmax,1)
                   !
                endif
                !
             enddo
          enddo
          !
       else if ( l == 2 ) then
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > hmin ) then
                   !
                   zs = 0.5 * ( zks(nm,kmax-1) - zks(nm,kmax) )
                   !
                   amatt(nm,kmax,1) = 1.
                   rhst (nm,kmax  ) = (cmu**.75) * (rtur(nm,kmax,1)**1.5) / ( vonkar*zs )
                   !
                else
                   !
                   rhst(nm,kmax) = rtur(nm,kmax,2)
                   !
                endif
                !
             enddo
          enddo
          !
       endif
       !
       ! solve the equation
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             bi = 1./amatt(nm,0,1)
             !
             amatt(nm,0,1) = bi
             amatt(nm,0,3) = amatt(nm,0,3)*bi
             rhst (nm,0  ) = rhst (nm,0  )*bi
             !
             do k = 1, kmax
                !
                bi = 1./(amatt(nm,k,1) - amatt(nm,k,2)*amatt(nm,k-1,3))
                amatt(nm,k,1) = bi
                amatt(nm,k,3) = amatt(nm,k,3)*bi
                rhst (nm,k  ) = (rhst(nm,k) - amatt(nm,k,2)*rhst(nm,k-1))*bi
                !
             enddo
             !
             rtur(nm,kmax,l) = rhst(nm,kmax)
             do k = kmax-1, 0, -1
                rtur(nm,k,l) = rhst(nm,k) - amatt(nm,k,3)*rtur(nm,k+1,l)
             enddo
             !
          enddo
       enddo
       !
       ! if negative values occur, apply clipping or give warning
       !
       if ( ITEST < 30 ) then
          !
          do k = 0, kmax
             do n = nfu, nl
                do m = mfu, ml
                   nm = kgrpnt(m,n)
                   if ( rtur(nm,k,l) < 0. ) rtur(nm,k,l) = eps
                enddo
             enddo
          enddo
          !
       else
          !
          do k = 0, kmax
             !
             do n = nfu, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( rtur(nm,k,l) < 0. ) then
                      write (PRINTF,'(a,i1)') ' ++ keps: negative value found for l=',l
                      write (PRINTF,'(a,i4,a,i4,a)') '          in (m,n)  = (',m+MXF-2,',',n+MYF-2,')'
                      write (PRINTF,'(a,i3)') '          and layer = ',k
                      write (PRINTF,'(a,f14.8)') '          value     = ',rtur(nm,k,l)
                      rtur(nm,k,l) = eps
                   endif
                   !
                enddo
             enddo
             !
          enddo
          !
       endif
       !
    enddo
    !
    ! exchange turbulence quantities with neighbouring subdomains, if appropriate
    !
    call SWEXCHG ( rtur(1,0,1), kgrpnt, 0, kmax )
    call SWEXCHG ( rtur(1,0,2), kgrpnt, 0, kmax )
    if (STPNOW()) return
    !
    ! synchronize turbulence quantities at appropriate boundaries in case of repeating grid, if appropriate
    !
    call periodic ( rtur(1,0,1), kgrpnt, 0, kmax )
    call periodic ( rtur(1,0,2), kgrpnt, 0, kmax )
    !
end subroutine SwashKepsMod2DH
