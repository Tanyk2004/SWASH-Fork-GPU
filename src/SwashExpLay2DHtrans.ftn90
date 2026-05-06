subroutine SwashExpLay2DHtrans
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
!    1.00, November 2013: New subroutine
!    1.31,     July 2014: extension sediment transport
!
!   Purpose
!
!   Performs the time integration for the layer-averaged 2D transport equations
!
!   Method
!
!   The time integration is fully explicit, except for the vertical terms (see below).
!
!   The time integration with respect to the horizontal advective term is based on the predictor-corrector
!   scheme of MacCormack, while that for the horizontal diffusive term is based on Euler explicit.
!
!   The space discretization of the horizontal terms is based on the finite volume approach. This
!   discretization is such that it is consistent with the discretization of the global continuity equation.
!
!   The horizontal advective term is approximated by either first order upwind or higher order (flux-limited)
!   scheme (CDS, Fromm, BDF, QUICK, MUSCL, Koren, etc.). The higher order scheme is treated by defect
!   correction consistent with the MacCormack scheme.
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
    use m_genarr, only: kgrpnt, guu, guv, gvu, gvv, gsqs
    use m_parall
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer        :: icistb   ! counter for number of instable points
    integer, save  :: ient = 0 ! number of entries in this subroutine
    integer        :: k        ! loop counter over vertical layers
    integer        :: kdd      ! index of layer k-2
    integer        :: ku       ! index of layer k+1
    integer        :: l        ! loop counter over constituents
    integer        :: m        ! loop counter
    integer        :: md       ! index of point m-1
    integer        :: mend     ! end index of loop over wl-points
    integer        :: msta     ! start index of loop over wl-points
    integer        :: mu       ! index of point m+1
    integer        :: muu      ! index of point m+2
    integer        :: n        ! loop counter
    integer        :: nd       ! index of point n-1
    integer        :: ndm      ! pointer to m,n-1
    integer        :: nend     ! end index of loop over wl-points in y-direction
    integer        :: nfm      ! pointer to m,nf
    integer        :: nfum     ! pointer to m,nfu
    integer        :: nlm      ! pointer to m,nl
    integer        :: nlum     ! pointer to m,nlu
    integer        :: nm       ! pointer to m,n
    integer        :: nmd      ! pointer to m-1,n
    integer        :: nmf      ! pointer to mf,n
    integer        :: nmfu     ! pointer to mfu,n
    integer        :: nml      ! pointer to ml,n
    integer        :: nmlu     ! pointer to mlu,n
    integer        :: nmu      ! pointer to m+1,n
    integer        :: nmuu     ! pointer to m+2,n
    integer        :: nsta     ! start index of loop over wl-points in y-direction
    integer        :: nu       ! index of point n+1
    integer        :: num      ! pointer to m,n+1
    integer        :: nuu      ! index of point n+2
    integer        :: nuum     ! pointer to m,n+2
    !
    real           :: bi       ! inverse of main diagonal of the matrix
    real           :: dc       ! higher order correction to advective flux
    real           :: dif2d    ! horizontal eddy diffusivity coefficient in velocity point
    real           :: dxl      ! local mesh size in x-direction
    real           :: dyl      ! local mesh size in y-direction
    real           :: fac      ! a factor
    real           :: fac1     ! another factor
    real           :: fac2     ! some other factor
    real           :: fluxlim  ! flux limiter
    real           :: grad1    ! solution gradient
    real           :: grad2    ! another solution gradient
    real           :: kwu      ! =1. if layer k+1 exists otherwise 0.
    real           :: mass     ! total mass
    real           :: psm      ! Prandtl-Schmidt number
    real           :: rproc    ! auxiliary variable with percentage of instable points
    real           :: stabmx   ! auxiliary variable with maximum diffusivity based stability criterion
    real           :: theta    ! implicitness factor for vertical terms
    real           :: w        ! local vertical velocity
    !
    logical        :: STPNOW   ! indicates that program must stop
    !
    character(120) :: msgstr   ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashExpLay2DHtrans')
    !
    theta = pnums(33)
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
          psm = psed(3)
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
       ! not computed for end point at subdomain interface since, this end point is owned by the neighbouring subdomain
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       nend = nl - 1
       if ( LMYL ) nend = nl
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             if ( .not. hs(nm) > epsdry ) then
                !
                amatc(nm,:,1) = 1.
                amatc(nm,:,2) = 0.
                amatc(nm,:,3) = 0.
                rhsc (nm,:  ) = rpo(nm,:)
                !
             endif
             !
          enddo
       enddo
       !
       ! compute the time derivative
       !
       do k = 1, kmax
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > epsdry ) then
                   !
                   amatc(nm,k,1) = hks (nm,k) / dt
                   rhsc (nm,k  ) = hkso(nm,k) * rpo(nm,k) / dt
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! compute advective flux and concentration at left and right boundaries using boundary conditions
       !
       if ( .not.lreptx ) then
          !
          if ( LMXF ) then
             !
             ! left boundary is ...
             !
             do n = nfu, nl
                !
                nmf  = kgrpnt(mf ,n)
                nmfu = kgrpnt(mfu,n)
                !
                if ( ibl(n) /= 1 ) then
                   !
                   ! ... open
                   !
                   do k = 1, kmax
                      !
                      if ( qx(nmf,k) > 0. ) then
                         !
                         ! inflow
                         !
                         flux(nmf,k,1) = qx(nmf,k) * rpo(nmf,k)
                         !
                         if ( l == lsal ) then
                            fac = max(icretl(n,k),0.) / max(tcret,dt)
                            rp(nmf,k,l) = coutl(n,k) + 0.5 * ( cbndl(n,k,l) - coutl(n,k) ) * ( 1. + cos(fac*pi) )
                            if ( .not. icretl(n,k) < 0. ) icretl(n,k) = icretl(n,k) - dt
                         else
                            rp(nmf,k,l) = cbndl(n,k,l)
                         endif
                         !
                      else
                         !
                         ! outflow
                         !
                         flux(nmf,k,1) = qx(nmf,k) * rpo(nmfu,k)
                         !
                         rp(nmf,k,l) = rpo(nmf,k) - dt * guu(nmf) * u1(nmf,k) * ( rpo(nmfu,k) - rpo(nmf,k) ) / gsqs(nmfu)
                         !
                         if ( l == lsal ) then
                            coutl (n,k) = rp(nmf,k,l)
                            icretl(n,k) = tcret
                         endif
                         !
                      endif
                      !
                   enddo
                   !
                else
                   !
                   ! ... closed
                   !
                   flux(nmf,:,1) = 0.            ! no advective and diffusive fluxes at closed boundary
                   rp  (nmf,:,l) = rp(nmfu,:,l)
                   !
                endif
                !
             enddo
             !
          endif
          !
          if ( LMXL ) then
             !
             ! right boundary is ...
             !
             do n = nfu, nl
                !
                nml  = kgrpnt(ml ,n)
                nmlu = kgrpnt(mlu,n)
                !
                if ( ibr(n) /= 1 ) then
                   !
                   ! ... open
                   !
                   do k = 1, kmax
                      !
                      if ( qx(nml,k) > 0. ) then
                         !
                         ! outflow
                         !
                         flux(nml,k,1) = qx(nml,k) * rpo(nml,k)
                         !
                         rp(nmlu,k,l) = rpo(nmlu,k) - dt * guu(nml) * u1(nml,k) * ( rpo(nmlu,k) - rpo(nml,k) ) / gsqs(nml)
                         !
                         if ( l == lsal ) then
                            coutr (n,k) = rp(nmlu,k,l)
                            icretr(n,k) = tcret
                         endif
                         !
                      else
                         !
                         ! inflow
                         !
                         flux(nml,k,1) = qx(nml,k) * rpo(nmlu,k)
                         !
                         if ( l == lsal ) then
                            fac = max(icretr(n,k),0.) / max(tcret,dt)
                            rp(nmlu,k,l) = coutr(n,k) + 0.5 * ( cbndr(n,k,l) - coutr(n,k) ) * ( 1. + cos(fac*pi) )
                            if ( .not. icretr(n,k) < 0. ) icretr(n,k) = icretr(n,k) - dt
                         else
                            rp(nmlu,k,l) = cbndr(n,k,l)
                         endif
                         !
                      endif
                      !
                   enddo
                   !
                else
                   !
                   ! ... closed
                   !
                   flux(nml ,:,1) = 0.           ! no advective and diffusive fluxes at closed boundary
                   rp  (nmlu,:,l) = rp(nml,:,l)
                   !
                endif
                !
             enddo
             !
          endif
          !
       endif
       !
       if ( .not.lrepty ) then
          !
          if ( LMYF ) then
             !
             ! lower boundary is ...
             !
             do m = mfu, ml
                !
                nfm  = kgrpnt(m,nf )
                nfum = kgrpnt(m,nfu)
                !
                if ( ibb(m) /= 1 ) then
                   !
                   ! ... open
                   !
                   do k = 1, kmax
                      !
                      if ( qy(nfm,k) > 0. ) then
                         !
                         ! inflow
                         !
                         flux(nfm,k,2) = qy(nfm,k) * rpo(nfm,k)
                         !
                         if ( l == lsal ) then
                            fac = max(icretb(m,k),0.) / max(tcret,dt)
                            rp(nfm,k,l) = coutb(m,k) + 0.5 * ( cbndb(m,k,l) - coutb(m,k) ) * ( 1. + cos(fac*pi) )
                            if ( .not. icretb(m,k) < 0. ) icretb(m,k) = icretb(m,k) - dt
                         else
                            rp(nfm,k,l) = cbndb(m,k,l)
                         endif
                         !
                      else
                         !
                         ! outflow
                         !
                         flux(nfm,k,2) = qy(nfm,k) * rpo(nfum,k)
                         !
                         rp(nfm,k,l) = rpo(nfm,k) - dt * gvv(nfm) * v1(nfm,k) * ( rpo(nfum,k) - rpo(nfm,k) ) / gsqs(nfum)
                         !
                         if ( l == lsal ) then
                            coutb (m,k) = rp(nfm,k,l)
                            icretb(m,k) = tcret
                         endif
                         !
                      endif
                      !
                   enddo
                   !
                else
                   !
                   ! ... closed
                   !
                   flux(nfm,:,2) = 0.            ! no advective and diffusive fluxes at closed boundary
                   rp  (nfm,:,l) = rp(nfum,:,l)
                   !
                endif
                !
             enddo
             !
          endif
          !
          if ( LMYL ) then
             !
             ! upper boundary is ...
             !
             do m = mfu, ml
                !
                nlm  = kgrpnt(m,nl )
                nlum = kgrpnt(m,nlu)
                !
                if ( ibt(m) /= 1 ) then
                   !
                   ! ... open
                   !
                   do k = 1, kmax
                      !
                      if ( qy(nlm,k) > 0. ) then
                         !
                         ! outflow
                         !
                         flux(nlm,k,2) = qy(nlm,k) * rpo(nlm,k)
                         !
                         rp(nlum,k,l) = rpo(nlum,k) - dt * gvv(nlm) * v1(nlm,k) * ( rpo(nlum,k) - rpo(nlm,k) ) / gsqs(nlm)
                         !
                         if ( l == lsal ) then
                            coutt (m,k) = rp(nlum,k,l)
                            icrett(m,k) = tcret
                         endif
                         !
                      else
                         !
                         ! inflow
                         !
                         flux(nlm,k,2) = qy(nlm,k) * rpo(nlum,k)
                         !
                         if ( l == lsal ) then
                            fac = max(icrett(m,k),0.) / max(tcret,dt)
                            rp(nlum,k,l) = coutt(m,k) + 0.5 * ( cbndt(m,k,l) - coutt(m,k) ) * ( 1. + cos(fac*pi) )
                            if ( .not. icrett(m,k) < 0. ) icrett(m,k) = icrett(m,k) - dt
                         else
                            rp(nlum,k,l) = cbndt(m,k,l)
                         endif
                         !
                      endif
                      !
                   enddo
                   !
                else
                   !
                   ! ... closed
                   !
                   flux(nlm ,:,2) = 0.           ! no advective and diffusive fluxes at closed boundary
                   rp  (nlum,:,l) = rp(nlm,:,l)
                   !
                endif
                !
             enddo
             !
          endif
          !
       endif
       !
       ! compute advective and diffusive flux in x-direction at internal cell-faces
       !
       if ( LMXF ) then
          msta = mf + 1         ! first internal u-point
       else
          msta = mf             ! first internal u-point at subdomain interface
       endif
       if ( lreptx ) then
          mend = ml             ! last  internal u-point in case of repeating grid
       else
          mend = ml - 1         ! last  internal u-point
       endif
       !
       do n = nfu, nl
          !
          do m = msta, mend
             !
             mu = m + 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmu = kgrpnt(mu,n )
             ndm = kgrpnt(m ,nd)
             !
             ! compute the advective flux based on first order upwind
             !
             do k = 1, kmax
                !
                if ( qx(nm,k) > 0. ) then
                   !
                   flux(nm,k,1) = qx(nm,k) * rpo(nm,k)
                   !
                else
                   !
                   flux(nm,k,1) = qx(nm,k) * rpo(nmu,k)
                   !
                endif
                !
             enddo
             !
             ! compute effective horizontal diffusivity coefficient in u-point
             !
             if ( icreep == 2 ) then
                !
                dif2d = 0.
                !
             else if ( hdiff > 0. ) then
                !
                dif2d = hdiff
                !
             else
                !
                if ( ihvisc == 2 .or. ihvisc == 3 ) then
                   !
                   dif2d = 0.5 * ( vnu2d(nm) + vnu2d(ndm) ) / psm
                   !
                else
                   !
                   dif2d = 0.
                   !
                endif
                !
             endif
             !
             ! check stability
             !
             dxl = gvu(nm)
             dyl = guu(nm)
             !
             stabmx = 1./(2.*dt*(1./(dxl*dxl) + 1./(dyl*dyl)))
             !
             if ( .not. dif2d < stabmx ) then
                dif2d  = stabmx
                icistb = icistb + 1
             endif
             !
             ! compute the diffusive flux along sigma planes and update total flux
             !
             do k = 1, kmax
                !
                flux(nm,k,1) = flux(nm,k,1) - 0.5 * ( hkso(nm,k) + hkso(nmu,k) ) * dif2d * guu(nm) * ( rpo(nmu,k) - rpo(nm,k) ) / gvu(nm)
                !
             enddo
             !
          enddo
          !
       enddo
       !
       ! synchronize fluxes at appropriate boundaries in case of repeating grid
       !
       call periodic ( flux(1,1,1), kgrpnt, 1, kmax )
       !
       ! compute advective and diffusive flux in y-direction at internal cell-faces
       !
       if ( LMYF ) then
          nsta = nf + 1         ! first internal v-point
       else
          nsta = nf             ! first internal v-point at subdomain interface
       endif
       if ( lrepty ) then
          nend = nl             ! last  internal v-point in case of repeating grid
       else
          nend = nl - 1         ! last  internal v-point
       endif
       !
       do m = mfu, ml
          !
          do n = nsta, nend
             !
             md = m - 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             num = kgrpnt(m ,nu)
             !
             ! compute the advective flux based on first order upwind
             !
             do k = 1, kmax
                !
                if ( qy(nm,k) > 0. ) then
                   !
                   flux(nm,k,2) = qy(nm,k) * rpo(nm,k)
                   !
                else
                   !
                   flux(nm,k,2) = qy(nm,k) * rpo(num,k)
                   !
                endif
                !
             enddo
             !
             ! compute effective horizontal diffusivity coefficient in v-point
             !
             if ( icreep == 2 ) then
                !
                dif2d = 0.
                !
             else if ( hdiff > 0. ) then
                !
                dif2d = hdiff
                !
             else
                !
                if ( ihvisc == 2 .or. ihvisc == 3 ) then
                   !
                   dif2d = 0.5 * ( vnu2d(nm) + vnu2d(nmd) ) / psm
                   !
                else
                   !
                   dif2d = 0.
                   !
                endif
                !
             endif
             !
             ! check stability
             !
             dxl = gvv(nm)
             dyl = guv(nm)
             !
             stabmx = 1./(2.*dt*(1./(dxl*dxl) + 1./(dyl*dyl)))
             !
             if ( .not. dif2d < stabmx ) then
                dif2d  = stabmx
                icistb = icistb + 1
             endif
             !
             ! compute the diffusive flux along sigma planes and update total flux
             !
             do k = 1, kmax
                !
                flux(nm,k,2) = flux(nm,k,2) - 0.5 * ( hkso(nm,k) + hkso(num,k) ) * dif2d * gvv(nm) * ( rpo(num,k) - rpo(nm,k) ) / guv(nm)
                !
             enddo
             !
          enddo
          !
       enddo
       !
       ! synchronize fluxes at appropriate boundaries in case of repeating grid
       !
       call periodic ( flux(1,1,2), kgrpnt, 1, kmax )
       !
       ! compute concentration in internal wl-points (based on finite volume approach)
       !
       ! not computed for end point at subdomain interface since, this end point is owned by the neighbouring subdomain
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       nend = nl - 1
       if ( LMYL ) nend = nl
       !
       do k = 1, kmax
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
                !
                if ( hs(nm) > epsdry ) then
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd == 1 ) nmd = nm
                   if ( ndm == 1 ) ndm = nm
                   !
                   rhsc(nm,k) = rhsc(nm,k) - ( flux(nm,k,1) - flux(nmd,k,1) + flux(nm,k,2) - flux(ndm,k,2) ) / gsqs(nm)
                   !
                endif
                !
             enddo
          enddo
          !
       enddo
       !
       ! compute the diffusive flux along strictly horizontal planes, if appropriate
       !
       if ( icreep == 2 ) call SwashHDiffZplane2DH ( rhsc, rpo, psm, icistb )
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
             do n = nfu, nend
                do m = mfu, mend
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( hs(nm) > epsdry ) then
                      !
                      ! friction velocity computed from turbulent kinetic energy at bed
                      ! may include wave breaking-induced turbulence
                      !
                      fac = ( 0.3 + psed(7) ) * rtur(nm,kmax,1)
                      !
                      ! include upward sediment flux (erosion) by means of pickup function
                      !
                      if ( fac > psed(5) ) rhsc(nm,kmax) = rhsc(nm,kmax) + psed(9) * ( ( fac - psed(5) ) / psed(5) )**1.5
                      !
                      ! include downward sediment flux (deposition) by means of fall velocity
                      !
                      amatc(nm,kmax,1) = amatc(nm,kmax,1) + psed(1)
                      !
                   endif
                   !
                enddo
             enddo
             !
          endif
          !
          if ( psed(10) > 0. ) then
             !
             ! cohesive sediment mud
             ! (for meaning of different sediment parameters psed(), see module SwashCommdata3)
             !
             do n = nfu, nend
                do m = mfu, mend
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( hs(nm) > epsdry ) then
                      !
                      ! bed shear stress computed from turbulent kinetic energy at bed
                      !
                      fac = 0.3 * rhow * rtur(nm,kmax,1)
                      !
                      ! include upward sediment flux (erosion)
                      !
                      if( fac > psed(10) ) rhsc(nm,kmax) = rhsc(nm,kmax) + psed(12) * ( fac/psed(10) - 1. )
                      !
                      ! include downward sediment flux (deposition)
                      !
                      if( fac < psed(11) ) amatc(nm,kmax,1) = amatc(nm,kmax,1) + psed(1) * ( 1. - fac/psed(11) )
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
       ! compute implicit part of vertical terms
       !
       propsc = nint(pnums(51))
       kappa  = pnums(52)
       mbound = pnums(53)
       phieby = pnums(54)
       !
       do k = 2, kmax
          !
          kdd = max(k-2,1   )
          ku  = min(k+1,kmax)
          !
          kwu = 1.
          if ( k == kmax ) kwu = 0.
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > epsdry ) then
                   !
                   w   = wom(nm,k-1)
                   fac = hks(nm,k-1) + hks(nm,k)
                   !
                   ! advection term
                   !
                   if ( propsc == 3 .and. kappa == 1. ) then
                      !
                      ! central differences
                      !
                      fac1 = theta * w * hks(nm,k-1) / fac
                      fac2 = theta * w * hks(nm,k  ) / fac
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
                      amatc(nm,k  ,1) = amatc(nm,k  ,1) + fac1
                      amatc(nm,k-1,1) = amatc(nm,k-1,1) - fac2
                      amatc(nm,k  ,2) =                   fac2
                      amatc(nm,k-1,3) =                 - fac1
                      !
                   else
                      !
                      ! first order upwind scheme
                      !
                      fac1 = theta * w
                      !
                      amatc(nm,k  ,1) = amatc(nm,k  ,1) + max(fac1,0.)
                      amatc(nm,k-1,1) = amatc(nm,k-1,1) - min(fac1,0.)
                      amatc(nm,k  ,2) =                   min(fac1,0.)
                      amatc(nm,k-1,3) =                 - max(fac1,0.)
                      !
                      ! add higher order (flux-limited) correction, if appropriate
                      !
                      if ( propsc /= 1 ) then
                         !
                         if ( w > 0. ) then
                            !
                            grad1 = rpo(nm,k-1) - rpo(nm,k )
                            grad2 = rpo(nm,k  ) - rpo(nm,ku)
                            !
                            fac2 = 0.5 * w * fluxlim(grad1,grad2)
                            !
                            rhsc(nm,k  ) = rhsc(nm,k  ) - fac2
                            rhsc(nm,k-1) = rhsc(nm,k-1) + fac2
                            !
                         else if ( w < 0. ) then
                            !
                            grad1 = rpo(nm,k  -1) - rpo(nm,k  )
                            grad2 = rpo(nm,kdd  ) - rpo(nm,k-1)
                            !
                            fac2 = 0.5 * w * fluxlim(grad1,grad2)
                            !
                            rhsc(nm,k  ) = rhsc(nm,k  ) + fac2
                            rhsc(nm,k-1) = rhsc(nm,k-1) - fac2
                            !
                         endif
                         !
                      endif
                      !
                   endif
                   !
                   ! diffusivity term
                   !
                   fac1 = 2. * theta * vnu3d(nm,k-1) / psm / fac
                   !
                   amatc(nm,k  ,1) = amatc(nm,k  ,1) + fac1
                   amatc(nm,k-1,1) = amatc(nm,k-1,1) + fac1
                   amatc(nm,k  ,2) = amatc(nm,k  ,2) - fac1
                   amatc(nm,k-1,3) = amatc(nm,k-1,3) - fac1
                   !
                   ! include fall velocity in case of sediment transport
                   !
                   if ( l == lsed .and. psed(1) > 0. ) then
                      !
                      amatc(nm,k-1,1) = amatc(nm,k-1,1) +     psed(1)
                      amatc(nm,k  ,2) = amatc(nm,k  ,2) - kwu*psed(1)
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
       ! compute explicit part of vertical terms, if appropriate
       !
       if ( theta /= 1. ) then
          !
          do k = 2, kmax
             !
             do n = nfu, nend
                do m = mfu, mend
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( hs(nm) > epsdry ) then
                      !
                      w   = wom(nm,k-1)
                      fac = hkso(nm,k-1) + hkso(nm,k)
                      !
                      ! advection term
                      !
                      if ( propsc == 3 .and. kappa == 1. ) then
                         !
                         ! central differences
                         !
                         fac1 = (1. - theta) * w * rpo(nm,k  ) * hkso(nm,k-1) / fac
                         fac2 = (1. - theta) * w * rpo(nm,k-1) * hkso(nm,k  ) / fac
                         !
                         if ( k == 2 .or. k == kmax ) then
                            !
                            if ( w > 0. ) then
                               !
                               fac1 = (1. - theta) * w * rpo(nm,k)
                               fac2 = 0.
                               !
                            else
                               !
                               fac1 = 0.
                               fac2 = (1. - theta) * w * rpo(nm,k-1)
                               !
                            endif
                            !
                         endif
                         !
                         rhsc(nm,k  ) = rhsc(nm,k  ) - fac1  - fac2
                         rhsc(nm,k-1) = rhsc(nm,k-1) + fac1  + fac2
                         !
                      else
                         !
                         ! first order upwind scheme
                         !
                         fac1 = (1. - theta) * max(w,0.) * rpo(nm,k  )
                         fac2 = (1. - theta) * min(w,0.) * rpo(nm,k-1)
                         !
                         rhsc(nm,k  ) = rhsc(nm,k  ) - fac1 - fac2
                         rhsc(nm,k-1) = rhsc(nm,k-1) + fac1 + fac2
                         !
                      endif
                      !
                      ! diffusivity term
                      !
                      fac1 = 2. * (1. - theta) * ( rpo(nm,k-1) - rpo(nm,k) ) * vnu3d(nm,k-1) / psm / fac
                      !
                      rhsc(nm,k  ) = rhsc(nm,k  ) + fac1
                      rhsc(nm,k-1) = rhsc(nm,k-1) - fac1
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
       ! add anti-creep terms to the transport equation, if appropriate
       !
       if ( icreep == 1 ) call SwashAntiCreep2DH ( amatc, rhsc, rpo, psm )
       !
       ! solve the transport equation
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             bi = 1./amatc(nm,1,1)
             !
             amatc(nm,1,1) = bi
             amatc(nm,1,3) = amatc(nm,1,3)*bi
             rhsc (nm,1  ) = rhsc (nm,1  )*bi
             !
             do k = 2, kmax
                !
                bi = 1./(amatc(nm,k,1) - amatc(nm,k,2)*amatc(nm,k-1,3))
                amatc(nm,k,1) = bi
                amatc(nm,k,3) = amatc(nm,k,3)*bi
                rhsc (nm,k  ) = (rhsc(nm,k) - amatc(nm,k,2)*rhsc(nm,k-1))*bi
                !
             enddo
             !
             rp(nm,kmax,l) = rhsc(nm,kmax)
             do k = kmax-1, 1, -1
                rp(nm,k,l) = rhsc(nm,k) - amatc(nm,k,3)*rp(nm,k+1,l)
             enddo
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       rp(1,:,l) = 0.
       !
       ! exchange concentrations with neighbouring subdomains
       !
       call SWEXCHG ( rp(1,1,l), kgrpnt, 1, kmax )
       if (STPNOW()) return
       !
       ! synchronize concentrations at appropriate boundaries in case of repeating grid
       !
       call periodic ( rp(1,1,l), kgrpnt, 1, kmax )
       !
       ! compute higher order correction to advective flux (if appropriate)
       !
       propsc = nint(pnums(46))
       !
       if ( propsc /= 1 ) then
          !
          kappa  = pnums(47)
          mbound = pnums(48)
          phieby = pnums(49)
          !
          if ( LMXF ) then
             msta = mf + 1         ! first internal u-point
          else
             msta = mf             ! first internal u-point at subdomain interface
          endif
          if ( lreptx ) then
             mend = ml             ! last  internal u-point in case of repeating grid
          else
             mend = ml - 1         ! last  internal u-point
          endif
          !
          do k = 1, kmax
             !
             do n = nfu, nl
                !
                nmf = kgrpnt(mf,n)
                nml = kgrpnt(ml,n)
                !
                flux(nmf,k,1) = 0.      ! no correction on left boundary
                flux(nml,k,1) = 0.      ! no correction on right boundary
                !
                do m = msta, mend
                   !
                   md  = m  - 1
                   mu  = m  + 1
                   muu = mu + 1
                   !
                   if ( lreptx .and. LMXL .and. muu > mlu ) muu = mlu
                   !
                   nm   = kgrpnt(m  ,n)
                   nmd  = kgrpnt(md ,n)
                   nmu  = kgrpnt(mu ,n)
                   nmuu = kgrpnt(muu,n)
                   !
                   if ( qx(nm,k) > 0. ) then
                      !
                      grad1 = rp(nmu,k,l) - rpo(nm ,k)
                      grad2 = rp(nm ,k,l) - rpo(nmd,k)
                      !
                      dc = 0.5 * fluxlim(grad1,grad2)
                      !
                   else
                      !
                      grad1 = rp(nm ,k,l) - rpo(nmu ,k)
                      grad2 = rp(nmu,k,l) - rpo(nmuu,k)
                      !
                      dc = 0.5 * fluxlim(grad1,grad2)
                      !
                   endif
                   !
                   flux(nm,k,1) = qx(nm,k) * dc
                   !
                enddo
                !
             enddo
             !
          enddo
          !
          ! synchronize corrections at appropriate boundaries in case of repeating grid
          !
          call periodic ( flux(1,1,1), kgrpnt, 1, kmax )
          !
       endif
       !
       ! compute higher order correction to advective flux in y-direction (if appropriate)
       !
       if ( propsc /= 1 ) then
          !
          if ( LMYF ) then
             nsta = nf + 1         ! first internal v-point
          else
             nsta = nf             ! first internal v-point at subdomain interface
          endif
          if ( lrepty ) then
             nend = nl             ! last  internal v-point in case of repeating grid
          else
             nend = nl - 1         ! last  internal v-point
          endif
          !
          do k = 1, kmax
             !
             do m = mfu, ml
                !
                nfm = kgrpnt(m,nf)
                nlm = kgrpnt(m,nl)
                !
                flux(nfm,k,2) = 0.      ! no correction on lower boundary
                flux(nlm,k,2) = 0.      ! no correction on upper boundary
                !
                do n = nsta, nend
                   !
                   nd  = n  - 1
                   nu  = n  + 1
                   nuu = nu + 1
                   !
                   if ( lrepty .and. LMYL .and. nuu > nlu ) nuu = nlu
                   !
                   nm   = kgrpnt(m,n  )
                   ndm  = kgrpnt(m,nd )
                   num  = kgrpnt(m,nu )
                   nuum = kgrpnt(m,nuu)
                   !
                   if ( qy(nm,k) > 0. ) then
                      !
                      grad1 = rp(num,k,l) - rpo(nm ,k)
                      grad2 = rp(nm ,k,l) - rpo(ndm,k)
                      !
                      dc = 0.5 * fluxlim(grad1,grad2)
                      !
                   else
                      !
                      grad1 = rp(nm ,k,l) - rpo(num ,k)
                      grad2 = rp(num,k,l) - rpo(nuum,k)
                      !
                      dc = 0.5 * fluxlim(grad1,grad2)
                      !
                   endif
                   !
                   flux(nm,k,2) = qy(nm,k) * dc
                   !
                enddo
                !
             enddo
             !
          enddo
          !
          ! synchronize corrections at appropriate boundaries in case of repeating grid
          !
          call periodic ( flux(1,1,2), kgrpnt, 1, kmax )
          !
       endif
       !
       ! correct concentration associated with higher order approximation of advection terms (if appropriate)
       !
       ! not computed for end point at subdomain interface since, this end point is owned by the neighbouring subdomain
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       nend = nl - 1
       if ( LMYL ) nend = nl
       !
       if ( propsc /= 1 ) then
          !
          do k = 1, kmax
             !
             do n = nfu, nend
                do m = mfu, mend
                   !
                   md = m - 1
                   nd = n - 1
                   !
                   nm  = kgrpnt(m ,n )
                   nmd = kgrpnt(md,n )
                   ndm = kgrpnt(m ,nd)
                   !
                   if ( hs(nm) > epsdry ) then
                      !
                      ! for permanently dry neighbours, corresponding values will be mirrored
                      !
                      if ( nmd == 1 ) nmd = nm
                      if ( ndm == 1 ) ndm = nm
                      !
                      ! update total concentration based on higher order correction flux at left, right, bottom and top cell-faces
                      !
                      rp(nm,k,l) = hks(nm,k) * rp(nm,k,l) - dt * ( flux(nm,k,1) - flux(nmd,k,1) + flux(nm,k,2) - flux(ndm,k,2) ) / gsqs(nm)
                      !
                      ! compute layer-averaged concentration
                      !
                      rp(nm,k,l) = rp(nm,k,l) / hks(nm,k)
                      !
                   endif
                   !
                enddo
             enddo
             !
          enddo
          !
          ! set to zero for permanently dry points
          !
          rp(1,:,l) = 0.
          !
       endif
       !
       ! exchange concentrations with neighbouring subdomains (if appropriate)
       !
       if ( propsc /= 1 ) then
          !
          call SWEXCHG ( rp(1,1,l), kgrpnt, 1, kmax )
          if (STPNOW()) return
          !
       endif
       !
       ! calculate total mass
       !
       if ( ITEST >= 30 ) then
          !
          mass = 0.
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                do k = 1, kmax
                   !
                   mass = mass + hks(nm,k) * gsqs(nm) * rp(nm,k,l)
                   !
                enddo
                !
             enddo
          enddo
          !
          ! accumulate mass over all subdomains
          !
          call SWREDUCE ( mass, 1, SWREAL, SWSUM )
          !
          if ( l == lsal  ) write(PRINTF,101) mass
          if ( l == ltemp ) write(PRINTF,102) mass
          if ( l == lsed  ) write(PRINTF,103) mass*rhos
          !
       endif
       !
       ! give warning for instable points
       !
       if ( icistb > 0 ) then
          !
          rproc = 100.*real(icistb)/real(mcgrd)
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
 101 format (2x,'the total mass associated with saline water is ',e14.8e2)
 102 format (2x,'the total mass associated with heat is ',e14.8e2)
 103 format (2x,'the total mass associated with suspended sediment ',e14.8e2)
    !
end subroutine SwashExpLay2DHtrans
