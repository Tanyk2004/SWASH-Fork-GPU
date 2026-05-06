subroutine SwashExpDep1DHtrans ( rp, rpo, u1, qx, flux, ibl, ibr, kgrpnt )
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
!    1.00, July 2013: New subroutine
!
!   Purpose
!
!   Performs the time integration for the depth-averaged 1D transport equations
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
    use m_parall
    use SwashFlowdata, fluxtmp => flux, &
                       ibltmp  => ibl , &
                       ibrtmp  => ibr , &
                       rptmp   => rp  , &
                       rpotmp  => rpo , &
                       u1tmp   => u1  , &
                       qxtmp   => qx
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                             :: ibl    ! boundary condition type at left boundary
    integer, intent(in)                             :: ibr    ! boundary condition type at right boundary
    !
    integer, dimension(mxc)         , intent(in)    :: kgrpnt ! index table containing the address of each (active) grid point
                                                              ! =1: not active grid point
                                                              ! >1: active grid point
    !
    real   , dimension(mcgrd)       , intent(out)   :: flux   ! total flux at cell-faces
    real   , dimension(mcgrd)       , intent(inout) :: qx     ! mass flux in x-direction
    real   , dimension(mcgrd,ltrans), intent(inout) :: rp     ! concentration at current time level
    real   , dimension(mcgrd)       , intent(out)   :: rpo    ! concentration at previous time level
    real   , dimension(mcgrd)       , intent(in)    :: u1     ! u-velocity at current time level
!
!   Local variables
!
    integer        :: icistb   ! counter for number of instable points
    integer, save  :: ient = 0 ! number of entries in this subroutine
    integer        :: l        ! loop counter over constituents
    integer        :: m        ! loop counter
    integer        :: md       ! index of point m-1
    integer        :: mend     ! end index of loop over wl-points
    integer        :: msta     ! start index of loop over wl-points
    integer        :: mu       ! index of point m+1
    integer        :: muu      ! index of point m+2
    integer        :: nm       ! pointer to m
    integer        :: nmd      ! pointer to m-1
    integer        :: nmf      ! pointer to mf
    integer        :: nmfu     ! pointer to mfu
    integer        :: nml      ! pointer to ml
    integer        :: nmlu     ! pointer to mlu
    integer        :: nmu      ! pointer to m+1
    integer        :: nmuu     ! pointer to m+2
    !
    real           :: dc       ! higher order correction to advective flux
    real           :: dif2d    ! horizontal eddy diffusivity coefficient in velocity point
    real           :: fac      ! a factor
    real           :: fluxlim  ! flux limiter
    real           :: grad1    ! solution gradient
    real           :: grad2    ! another solution gradient
    real           :: mass     ! total mass
    real           :: psm      ! Prandtl-Schmidt number
    real           :: rdx      ! reciprocal of mesh size
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
    if (ltrace) call strace (ient,'SwashExpDep1DHtrans')
    !
    nmf  = kgrpnt(mf )
    nmfu = kgrpnt(mfu)
    nml  = kgrpnt(ml )
    nmlu = kgrpnt(mlu)
    !
    rdx = 1./dx
    !
    stabmx = 0.5 * dx * dx / dt
    !
    ! if momentum equation has been skipped, compute the mass flux
    !
    if ( momskip ) qx = hu * u1
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
       ! compute advective flux and concentration at left and right boundaries using boundary conditions
       !
       if ( LMXF ) then
          !
          ! left boundary is ...
          !
          if ( ibl /= 1 ) then
             !
             ! ... open
             !
             if ( qx(nmf) > 0. ) then
                !
                ! inflow
                !
                flux(nmf) = qx(nmf) * rpo(nmf)
                !
                if ( l == lsal ) then
                   fac = max(icretl(1,1),0.) / max(tcret,dt)
                   rp(nmf,l) = coutl(1,1) + 0.5 * ( cbndl(1,1,l) - coutl(1,1) ) * ( 1. + cos(fac*pi) )
                   if ( .not. icretl(1,1) < 0. ) icretl(1,1) = icretl(1,1) - dt
                else
                   rp(nmf,l) = cbndl(1,1,l)
                endif
                !
             else
                !
                ! outflow
                !
                flux(nmf) = qx(nmf) * rpo(nmfu)
                !
                rp(nmf,l) = rpo(nmf) - dt * rdx * u1(nmf) * ( rpo(nmfu) - rpo(nmf) )
                !
                if ( l == lsal ) then
                   coutl (1,1) = rp(nmf,l)
                   icretl(1,1) = tcret
                endif
                !
             endif
             !
          else
             !
             ! ... closed
             !
             flux(nmf  ) = 0.          ! no advective and diffusive fluxes at closed boundary
             rp  (nmf,l) = rp(nmfu,l)
             !
          endif
          !
       endif
       !
       if ( LMXL ) then
          !
          ! right boundary is ...
          !
          if ( ibr /= 1 ) then
             !
             ! ... open
             !
             if ( qx(nml) > 0. ) then
                !
                ! outflow
                !
                flux(nml) = qx(nml) * rpo(nml)
                !
                rp(nmlu,l) = rpo(nmlu) - dt * rdx * u1(nml) * ( rpo(nmlu) - rpo(nml) )
                !
                if ( l == lsal ) then
                   coutr (1,1) = rp(nmlu,l)
                   icretr(1,1) = tcret
                endif
                !
             else
                !
                ! inflow
                !
                flux(nml) = qx(nml) * rpo(nmlu)
                !
                if ( l == lsal ) then
                   fac = max(icretr(1,1),0.) / max(tcret,dt)
                   rp(nmlu,l) = coutr(1,1) + 0.5 * ( cbndr(1,1,l) - coutr(1,1) ) * ( 1. + cos(fac*pi) )
                   if ( .not. icretr(1,1) < 0. ) icretr(1,1) = icretr(1,1) - dt
                else
                   rp(nmlu,l) = cbndr(1,1,l)
                endif
                !
             endif
             !
          else
             !
             ! ... closed
             !
             flux(nml   ) = 0.         ! no advective and diffusive fluxes at closed boundary
             rp  (nmlu,l) = rp(nml,l)
             !
          endif
          !
       endif
       !
       ! compute advective and diffusive flux at internal cell-faces
       !
       if ( LMXF ) then
          msta = mf + 1   ! first internal cell-face
       else
          msta = mf       ! first internal cell-face at subdomain interface
       endif
       mend = ml - 1      ! last  internal cell-face
       !
       do m = msta, mend
          !
          mu = m + 1
          !
          nm  = kgrpnt(m )
          nmu = kgrpnt(mu)
          !
          ! compute the advective flux based on first order upwind
          !
          if ( qx(nm) > 0. ) then
             !
             flux(nm) = qx(nm) * rpo(nm)
             !
          else
             !
             flux(nm) = qx(nm) * rpo(nmu)
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
                dif2d = 0.5 * ( vnu2d(nm) + vnu2d(nmu) ) / psm
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
          if ( .not. dif2d < stabmx ) then
             dif2d  = stabmx
             icistb = icistb + 1
          endif
          !
          ! compute the diffusive flux and update total flux
          !
          flux(nm) = flux(nm) - 0.5 * ( hso(nm) + hso(nmu) ) * dif2d * rdx * ( rpo(nmu) - rpo(nm) )
          !
       enddo
       !
       ! compute concentration in internal wl-points (based on finite volume approach)
       !
       ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       do m = mfu, mend
          !
          nm = kgrpnt(m)
          !
          if ( hs(nm) > epsdry ) then
             !
             ! update total concentration based on total flux at right cell-face
             !
             rp(nm,l) = hso(nm) * rpo(nm) - dt * rdx * flux(nm)
             !
          endif
          !
       enddo
       !
       do m = mfu, mend
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          if ( hs(nm) > epsdry ) then
             !
             ! update total concentration based on total flux at left cell-face
             !
             rp(nm,l) = rp(nm,l) + dt * rdx * flux(nmd)
             !
             ! compute depth-averaged concentration
             !
             rp(nm,l) = rp(nm,l) / hs(nm)
             !
          endif
          !
       enddo
       !
       ! exchange concentrations with neighbouring subdomains
       !
       call SWEXCHG ( rp(1,l), kgrpnt, 1, 1 )
       if (STPNOW()) return
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
          flux(nmf) = 0.     ! no correction on left boundary
          flux(nml) = 0.     ! no correction on right boundary
          !
          if ( LMXF ) then
             msta = mf + 1   ! first internal cell-face
          else
             msta = mf       ! first internal cell-face at subdomain interface
          endif
          mend = ml - 1      ! last  internal cell-face
          !
          do m = msta, mend
             !
             md  = m  - 1
             mu  = m  + 1
             muu = mu + 1
             !
             nm   = kgrpnt(m  )
             nmd  = kgrpnt(md )
             nmu  = kgrpnt(mu )
             nmuu = kgrpnt(muu)
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
             flux(nm) = qx(nm) * dc
             !
          enddo
          !
       endif
       !
       ! correct concentration associated with higher order approximation of advection term (if appropriate)
       !
       ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       if ( propsc /= 1 ) then
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m)
             !
             if ( hs(nm) > epsdry ) then
                !
                ! update total concentration based on higher order correction flux at right cell-face
                !
                rp(nm,l) = hs(nm) * rp(nm,l) - dt * rdx * flux(nm)
                !
             endif
             !
          enddo
          !
          do m = mfu, mend
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             if ( hs(nm) > epsdry ) then
                !
                ! update total concentration based on higher order correction flux at left cell-face
                !
                rp(nm,l) = rp(nm,l) + dt * rdx * flux(nmd)
                !
                ! compute depth-averaged concentration
                !
                rp(nm,l) = rp(nm,l) / hs(nm)
                !
             endif
             !
          enddo
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
       ! calculate total mass
       !
       if ( ITEST >= 30 ) then
          !
          mass = 0.
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m)
             !
             mass = mass + hs(nm) * dx * rp(nm,l)
             !
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
end subroutine SwashExpDep1DHtrans
