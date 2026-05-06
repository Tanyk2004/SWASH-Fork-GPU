subroutine SwashExpDep2DHtrans ( rp, rpo, u1, qx, v1, qy, flux )
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
!
!   Purpose
!
!   Performs the time integration for the depth-averaged 2D transport equations
!
!   Method
!
!   The time integration is fully explicit.
!
!   The time integration with respect to the advective term is based on the predictor-corrector scheme
!   of MacCormack, while that for the diffusive term is based on Euler explicit.
!
!   The space discretization is based on the finite volume approach. The discretization is such that
!   it is consistent with the discretization of the global continuity equation.
!
!   The advective term is approximated by either first order upwind or higher order (flux-limited) scheme
!   (CDS, Fromm, BDF, QUICK, MUSCL, Koren, etc.). The higher order scheme is treated by defect correction
!   consistent with the MacCormack scheme.
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
    use SwashFlowdata, fluxtmp => flux, &
                       rptmp   => rp  , &
                       rpotmp  => rpo , &
                       u1tmp   => u1  , &
                       qxtmp   => qx  , &
                       v1tmp   => v1  , &
                       qytmp   => qy
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,2)     , intent(out)   :: flux ! total flux at cell-faces
    real, dimension(mcgrd)       , intent(inout) :: qx   ! mass flux in x-direction
    real, dimension(mcgrd)       , intent(inout) :: qy   ! mass flux in y-direction
    real, dimension(mcgrd,ltrans), intent(inout) :: rp   ! concentration at current time level
    real, dimension(mcgrd)       , intent(out)   :: rpo  ! concentration at previous time level
    real, dimension(mcgrd)       , intent(in)    :: u1   ! u-velocity at current time level
    real, dimension(mcgrd)       , intent(in)    :: v1   ! v-velocity at current time level
!
!   Local variables
!
    integer        :: icistb   ! counter for number of instable points
    integer, save  :: ient = 0 ! number of entries in this subroutine
    integer        :: l        ! loop counter over constituents
    integer        :: m        ! loop counter
    integer        :: md       ! index of point m-1
    integer        :: mend     ! end index of loop over wl-points in x-direction
    integer        :: msta     ! start index of loop over wl-points in x-direction
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
    real           :: dc       ! higher order correction to advective flux
    real           :: dif2d    ! horizontal eddy diffusivity coefficient in velocity point
    real           :: dxl      ! local mesh size in x-direction
    real           :: dyl      ! local mesh size in y-direction
    real           :: fac      ! a factor
    real           :: fluxlim  ! flux limiter
    real           :: grad1    ! solution gradient
    real           :: grad2    ! another solution gradient
    real           :: mass     ! total mass
    real           :: psm      ! Prandtl-Schmidt number
    real           :: rproc    ! auxiliary variable with percentage of instable points
    real           :: stabmx   ! auxiliary variable with maximum diffusivity based stability criterion
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
    if (ltrace) call strace (ient,'SwashExpDep2DHtrans')
    !
    ! if momentum equations have been skipped, compute the mass flux in x- and y-direction
    !
    if ( momskip ) then
       qx = hu * guu * u1
       qy = hv * gvv * v1
    endif
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
       ! compute advective flux and concentration at all non-periodic boundaries using boundary conditions
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
                   if ( qx(nmf) > 0. ) then
                      !
                      ! inflow
                      !
                      flux(nmf,1) = qx(nmf) * rpo(nmf)
                      !
                      if ( l == lsal ) then
                         fac = max(icretl(n,1),0.) / max(tcret,dt)
                         rp(nmf,l) = coutl(n,1) + 0.5 * ( cbndl(n,1,l) - coutl(n,1) ) * ( 1. + cos(fac*pi) )
                         if ( .not. icretl(n,1) < 0. ) icretl(n,1) = icretl(n,1) - dt
                      else
                         rp(nmf,l) = cbndl(n,1,l)
                      endif
                      !
                   else
                      !
                      ! outflow
                      !
                      flux(nmf,1) = qx(nmf) * rpo(nmfu)
                      !
                      rp(nmf,l) = rpo(nmf) - dt * guu(nmf) * u1(nmf) * ( rpo(nmfu) - rpo(nmf) ) / gsqs(nmfu)
                      !
                      if ( l == lsal ) then
                         coutl (n,1) = rp(nmf,l)
                         icretl(n,1) = tcret
                      endif
                      !
                   endif
                   !
                else
                   !
                   ! ... closed
                   !
                   flux(nmf,1) = 0.          ! no advective and diffusive fluxes at closed boundary
                   rp  (nmf,l) = rp(nmfu,l)
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
                   if ( qx(nml) > 0. ) then
                      !
                      ! outflow
                      !
                      flux(nml,1) = qx(nml) * rpo(nml)
                      !
                      rp(nmlu,l) = rpo(nmlu) - dt * guu(nml) * u1(nml) * ( rpo(nmlu) - rpo(nml) ) / gsqs(nml)
                      !
                      if ( l == lsal ) then
                         coutr (n,1) = rp(nmlu,l)
                         icretr(n,1) = tcret
                      endif
                      !
                   else
                      !
                      ! inflow
                      !
                      flux(nml,1) = qx(nml) * rpo(nmlu)
                      !
                      if ( l == lsal ) then
                         fac = max(icretr(n,1),0.) / max(tcret,dt)
                         rp(nmlu,l) = coutr(n,1) + 0.5 * ( cbndr(n,1,l) - coutr(n,1) ) * ( 1. + cos(fac*pi) )
                         if ( .not. icretr(n,1) < 0. ) icretr(n,1) = icretr(n,1) - dt
                      else
                         rp(nmlu,l) = cbndr(n,1,l)
                      endif
                      !
                   endif
                   !
                else
                   !
                   ! ... closed
                   !
                   flux(nml ,1) = 0.         ! no advective and diffusive fluxes at closed boundary
                   rp  (nmlu,l) = rp(nml,l)
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
                   if ( qy(nfm) > 0. ) then
                      !
                      ! inflow
                      !
                      flux(nfm,2) = qy(nfm) * rpo(nfm)
                      !
                      if ( l == lsal ) then
                         fac = max(icretb(m,1),0.) / max(tcret,dt)
                         rp(nfm,l) = coutb(m,1) + 0.5 * ( cbndb(m,1,l) - coutb(m,1) ) * ( 1. + cos(fac*pi) )
                         if ( .not. icretb(m,1) < 0. ) icretb(m,1) = icretb(m,1) - dt
                      else
                         rp(nfm,l) = cbndb(m,1,l)
                      endif
                      !
                   else
                      !
                      ! outflow
                      !
                      flux(nfm,2) = qy(nfm) * rpo(nfum)
                      !
                      rp(nfm,l) = rpo(nfm) - dt * gvv(nfm) * v1(nfm) * ( rpo(nfum) - rpo(nfm) ) / gsqs(nfum)
                      !
                      if ( l == lsal ) then
                         coutb (m,1) = rp(nfm,l)
                         icretb(m,1) = tcret
                      endif
                      !
                   endif
                   !
                else
                   !
                   ! ... closed
                   !
                   flux(nfm,2) = 0.          ! no advective and diffusive fluxes at closed boundary
                   rp  (nfm,l) = rp(nfum,l)
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
                   if ( qy(nlm) > 0. ) then
                      !
                      ! outflow
                      !
                      flux(nlm,2) = qy(nlm) * rpo(nlm)
                      !
                      rp(nlum,l) = rpo(nlum) - dt * gvv(nlm) * v1(nlm) * ( rpo(nlum) - rpo(nlm) ) / gsqs(nlm)
                      !
                      if ( l == lsal ) then
                         coutt (m,1) = rp(nlum,l)
                         icrett(m,1) = tcret
                      endif
                      !
                   else
                      !
                      ! inflow
                      !
                      flux(nlm,2) = qy(nlm) * rpo(nlum)
                      !
                      if ( l == lsal ) then
                         fac = max(icrett(m,1),0.) / max(tcret,dt)
                         rp(nlum,l) = coutt(m,1) + 0.5 * ( cbndt(m,1,l) - coutt(m,1) ) * ( 1. + cos(fac*pi) )
                         if ( .not. icrett(m,1) < 0. ) icrett(m,1) = icrett(m,1) - dt
                      else
                         rp(nlum,l) = cbndt(m,1,l)
                      endif
                      !
                   endif
                   !
                else
                   !
                   ! ... closed
                   !
                   flux(nlm ,2) = 0.         ! no advective and diffusive fluxes at closed boundary
                   rp  (nlum,l) = rp(nlm,l)
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
             nm   = kgrpnt(m ,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             !
             ! compute the advective flux based on first order upwind
             !
             if ( qx(nm) > 0. ) then
                !
                flux(nm,1) = qx(nm) * rpo(nm)
                !
             else
                !
                flux(nm,1) = qx(nm) * rpo(nmu)
                !
             endif
             !
             ! compute effective horizontal diffusivity coefficient in u-point
             !
             if ( hdiff > 0. ) then
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
             ! compute the diffusive flux and update total flux
             !
             flux(nm,1) = flux(nm,1) - 0.5 * ( hso(nm) + hso(nmu) ) * dif2d * guu(nm) * ( rpo(nmu) - rpo(nm) ) / gvu(nm)
             !
          enddo
          !
       enddo
       !
       ! synchronize fluxes at appropriate boundaries in case of repeating grid
       !
       call periodic ( flux(1,1), kgrpnt, 1, 1 )
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
             nm   = kgrpnt(m ,n )
             nmd  = kgrpnt(md,n )
             num  = kgrpnt(m ,nu)
             !
             ! compute the advective flux based on first order upwind
             !
             if ( qy(nm) > 0. ) then
                !
                flux(nm,2) = qy(nm) * rpo(nm)
                !
             else
                !
                flux(nm,2) = qy(nm) * rpo(num)
                !
             endif
             !
             ! compute effective horizontal diffusivity coefficient in v-point
             !
             if ( hdiff > 0. ) then
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
             ! compute the diffusive flux and update total flux
             !
             flux(nm,2) = flux(nm,2) - 0.5 * ( hso(nm) + hso(num) ) * dif2d * gvv(nm) * ( rpo(num) - rpo(nm) ) / guv(nm)
             !
          enddo
          !
       enddo
       !
       ! synchronize fluxes at appropriate boundaries in case of repeating grid
       !
       call periodic ( flux(1,2), kgrpnt, 1, 1 )
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
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             if ( hs(nm) > epsdry ) then
                !
                ! update total concentration based on total flux at right cell-face
                !
                rp(nm,l) = hso(nm) * rpo(nm) - dt * flux(nm,1) / gsqs(nm)
                !
             endif
             !
          enddo
       enddo
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,n)
             nmd = kgrpnt(md,n)
             !
             if ( hs(nm) > epsdry ) then
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd == 1 ) nmd = nm
                !
                ! update total concentration based on total flux at left cell-face
                !
                rp(nm,l) = rp(nm,l) + dt * flux(nmd,1) / gsqs(nm)
                !
             endif
             !
          enddo
       enddo
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             if ( hs(nm) > epsdry ) then
                !
                ! update total concentration based on total flux at top cell-face
                !
                rp(nm,l) = rp(nm,l) - dt * flux(nm,2) / gsqs(nm)
                !
             endif
             !
          enddo
       enddo
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nd = n - 1
             !
             nm  = kgrpnt(m,n )
             ndm = kgrpnt(m,nd)
             !
             if ( hs(nm) > epsdry ) then
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( ndm == 1 ) ndm = nm
                !
                ! update total concentration based on total flux at bottom cell-face
                !
                rp(nm,l) = rp(nm,l) + dt * flux(ndm,2) / gsqs(nm)
                !
                ! compute depth-averaged concentration
                !
                rp(nm,l) = rp(nm,l) / hs(nm)
                !
             endif
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       rp(1,l) = 0.
       !
       ! exchange concentrations with neighbouring subdomains
       !
       call SWEXCHG ( rp(1,l), kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize concentrations at appropriate boundaries in case of repeating grid
       !
       call periodic ( rp(1,l), kgrpnt, 1, 1 )
       !
       ! compute higher order correction to advective flux in x-direction (if appropriate)
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
          do n = nfu, nl
             !
             nmf = kgrpnt(mf,n)
             nml = kgrpnt(ml,n)
             !
             flux(nmf,1) = 0.      ! no correction on left boundary
             flux(nml,1) = 0.      ! no correction on right boundary
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
                if ( qx(nm) > 0. ) then
                   !
                   grad1 = rp(nmu,l) - rpo(nm )
                   grad2 = rp(nm ,l) - rpo(nmd)
                   !
                   dc = 0.5 * fluxlim(grad1,grad2)
                   !
                else
                   !
                   grad1 = rp(nm ,l) - rpo(nmu )
                   grad2 = rp(nmu,l) - rpo(nmuu)
                   !
                   dc = 0.5 * fluxlim(grad1,grad2)
                   !
                endif
                !
                flux(nm,1) = qx(nm) * dc
                !
             enddo
             !
          enddo
          !
          ! synchronize corrections at appropriate boundaries in case of repeating grid
          !
          call periodic ( flux(1,1), kgrpnt, 1, 1 )
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
          do m = mfu, ml
             !
             nfm = kgrpnt(m,nf)
             nlm = kgrpnt(m,nl)
             !
             flux(nfm,2) = 0.      ! no correction on lower boundary
             flux(nlm,2) = 0.      ! no correction on upper boundary
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
                if ( qy(nm) > 0. ) then
                   !
                   grad1 = rp(num,l) - rpo(nm )
                   grad2 = rp(nm ,l) - rpo(ndm)
                   !
                   dc = 0.5 * fluxlim(grad1,grad2)
                   !
                else
                   !
                   grad1 = rp(nm ,l) - rpo(num )
                   grad2 = rp(num,l) - rpo(nuum)
                   !
                   dc = 0.5 * fluxlim(grad1,grad2)
                   !
                endif
                !
                flux(nm,2) = qy(nm) * dc
                !
             enddo
             !
          enddo
          !
          ! synchronize corrections at appropriate boundaries in case of repeating grid
          !
          call periodic ( flux(1,2), kgrpnt, 1, 1 )
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
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > epsdry ) then
                   !
                   ! update total concentration based on higher order correction flux at right cell-face
                   !
                   rp(nm,l) = hs(nm) * rp(nm,l) - dt * flux(nm,1) / gsqs(nm)
                   !
                endif
                !
             enddo
          enddo
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                md = m - 1
                !
                nm  = kgrpnt(m ,n)
                nmd = kgrpnt(md,n)
                !
                if ( hs(nm) > epsdry ) then
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd == 1 ) nmd = nm
                   !
                   ! update total concentration based on higher order correction flux at left cell-face
                   !
                   rp(nm,l) = rp(nm,l) + dt * flux(nmd,1) / gsqs(nm)
                   !
                endif
                !
             enddo
          enddo
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                if ( hs(nm) > epsdry ) then
                   !
                   ! update total concentration based on higher order correction flux at top cell-face
                   !
                   rp(nm,l) = rp(nm,l) - dt * flux(nm,2) / gsqs(nm)
                   !
                endif
                !
             enddo
          enddo
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nd = n - 1
                !
                nm  = kgrpnt(m,n )
                ndm = kgrpnt(m,nd)
                !
                if ( hs(nm) > epsdry ) then
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( ndm == 1 ) ndm = nm
                   !
                   ! update total concentration based on higher order correction flux at bottom cell-face
                   !
                   rp(nm,l) = rp(nm,l) + dt * flux(ndm,2) / gsqs(nm)
                   !
                   ! compute depth-averaged concentration
                   !
                   rp(nm,l) = rp(nm,l) / hs(nm)
                   !
                endif
                !
             enddo
          enddo
          !
          ! set to zero for permanently dry points
          !
          rp(1,l) = 0.
          !
       endif
       !
       ! exchange concentrations with neighbouring subdomains (if appropriate)
       !
       if ( propsc /= 1 ) then
          !
          call SWEXCHG ( rp(1,l), kgrpnt, 1, 1 )
          if (STPNOW()) return
          !
       endif
       !
       ! synchronize concentrations at appropriate boundaries in case of repeating grid (if appropriate)
       !
       if ( propsc /= 1 ) call periodic ( rp(1,l), kgrpnt, 1, 1 )
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
                mass = mass + hs(nm) * gsqs(nm) * rp(nm,l)
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
end subroutine SwashExpDep2DHtrans
