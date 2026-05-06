subroutine SwashExpLay1DHtrans ( ibl, ibr, kgrpnt )
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
!   Performs the time integration for the layer-averaged 1D transport equations
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
    use m_parall
    use SwashFlowdata, ibltmp => ibl, &
                       ibrtmp => ibr
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                 :: ibl    ! boundary condition type at left boundary
    integer, intent(in)                 :: ibr    ! boundary condition type at right boundary
    !
    integer, dimension(mxc), intent(in) :: kgrpnt ! index table containing the address of each (active) grid point
                                                  ! =1: not active grid point
                                                  ! >1: active grid point
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
    integer        :: nm       ! pointer to m
    integer        :: nmd      ! pointer to m-1
    integer        :: nmf      ! pointer to mf
    integer        :: nmfu     ! pointer to mfu
    integer        :: nml      ! pointer to ml
    integer        :: nmlu     ! pointer to mlu
    integer        :: nmu      ! pointer to m+1
    integer        :: nmuu     ! pointer to m+2
    !
    real           :: bi       ! inverse of main diagonal of the matrix
    real           :: dc       ! higher order correction to advective flux
    real           :: dif2d    ! horizontal eddy diffusivity coefficient in velocity point
    real           :: fac      ! a factor
    real           :: fac1     ! another factor
    real           :: fac2     ! some other factor
    real           :: fluxlim  ! flux limiter
    real           :: grad1    ! solution gradient
    real           :: grad2    ! another solution gradient
    real           :: kwu      ! =1. if layer k+1 exists otherwise 0.
    real           :: mass     ! total mass
    real           :: psm      ! Prandtl-Schmidt number
    real           :: rdx      ! reciprocal of mesh size
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
    if (ltrace) call strace (ient,'SwashExpLay1DHtrans')
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
       ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       do m = mfu, mend
          !
          nm = kgrpnt(m)
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
       !
       ! compute the time derivative
       !
       do k = 1, kmax
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m)
             !
             if ( hs(nm) > epsdry ) then
                !
                amatc(nm,k,1) = hks (nm,k) / dt
                rhsc (nm,k  ) = hkso(nm,k) * rpo(nm,k) / dt
                !
             endif
             !
          enddo
          !
       enddo
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
             do k = 1, kmax
                !
                if ( qx(nmf,k) > 0. ) then
                   !
                   ! inflow
                   !
                   flux(nmf,k,1) = qx(nmf,k) * rpo(nmf,k)
                   !
                   if ( l == lsal ) then
                      fac = max(icretl(1,k),0.) / max(tcret,dt)
                      rp(nmf,k,l) = coutl(1,k) + 0.5 * ( cbndl(1,k,l) - coutl(1,k) ) * ( 1. + cos(fac*pi) )
                      if ( .not. icretl(1,k) < 0. ) icretl(1,k) = icretl(1,k) - dt
                   else
                      rp(nmf,k,l) = cbndl(1,k,l)
                   endif
                   !
                else
                   !
                   ! outflow
                   !
                   flux(nmf,k,1) = qx(nmf,k) * rpo(nmfu,k)
                   !
                   rp(nmf,k,l) = rpo(nmf,k) - dt * rdx * u1(nmf,k) * ( rpo(nmfu,k) - rpo(nmf,k) )
                   !
                   if ( l == lsal ) then
                      coutl (1,k) = rp(nmf,k,l)
                      icretl(1,k) = tcret
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
             do k = 1, kmax
                !
                if ( qx(nml,k) > 0. ) then
                   !
                   ! outflow
                   !
                   flux(nml,k,1) = qx(nml,k) * rpo(nml,k)
                   !
                   rp(nmlu,k,l) = rpo(nmlu,k) - dt * rdx * u1(nml,k) * ( rpo(nmlu,k) - rpo(nml,k) )
                   !
                   if ( l == lsal ) then
                      coutr (1,k) = rp(nmlu,k,l)
                      icretr(1,k) = tcret
                   endif
                   !
                else
                   !
                   ! inflow
                   !
                   flux(nml,k,1) = qx(nml,k) * rpo(nmlu,k)
                   !
                   if ( l == lsal ) then
                      fac = max(icretr(1,k),0.) / max(tcret,dt)
                      rp(nmlu,k,l) = coutr(1,k) + 0.5 * ( cbndr(1,k,l) - coutr(1,k) ) * ( 1. + cos(fac*pi) )
                      if ( .not. icretr(1,k) < 0. ) icretr(1,k) = icretr(1,k) - dt
                   else
                      rp(nmlu,k,l) = cbndr(1,k,l)
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
          ! compute the diffusive flux along sigma planes and update total flux
          !
          do k = 1, kmax
             !
             flux(nm,k,1) = flux(nm,k,1) - 0.5 * ( hkso(nm,k) + hkso(nmu,k) ) * dif2d * rdx * ( rpo(nmu,k) - rpo(nm,k) )
             !
          enddo
          !
       enddo
       !
       ! compute horizontal terms (based on finite volume approach)
       !
       ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       do k = 1, kmax
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
                rhsc(nm,k) = rhsc(nm,k) - rdx * ( flux(nm,k,1) - flux(nmd,k,1) )
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! compute the diffusive flux along strictly horizontal planes, if appropriate
       !
       if ( icreep == 2 ) call SwashHDiffZplane1DH ( rhsc, rpo, kgrpnt, psm, icistb )
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
             do m = mfu, mend
                !
                nm = kgrpnt(m)
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
             !
          endif
          !
          if ( psed(10) > 0. ) then
             !
             ! cohesive sediment mud
             ! (for meaning of different sediment parameters psed(), see module SwashCommdata3)
             !
             do m = mfu, mend
                !
                nm = kgrpnt(m)
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
          do m = mfu, mend
             !
             nm = kgrpnt(m)
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
          !
       enddo
       !
       ! compute explicit part of vertical terms, if appropriate
       !
       if ( theta /= 1. ) then
          !
          do k = 2, kmax
             !
             do m = mfu, mend
                !
                nm = kgrpnt(m)
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
             !
          enddo
          !
       endif
       !
       ! add anti-creep terms to the transport equation, if appropriate
       !
       if ( icreep == 1 ) call SwashAntiCreep1DH ( amatc, rhsc, rpo, kgrpnt, psm )
       !
       ! solve the transport equation
       !
       do m = mfu, mend
          !
          nm = kgrpnt(m)
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
       !
       ! exchange concentrations with neighbouring subdomains
       !
       call SWEXCHG ( rp(1,1,l), kgrpnt, 1, kmax )
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
          flux(nmf,:,1) = 0.   ! no correction on left boundary
          flux(nml,:,1) = 0.   ! no correction on right boundary
          !
          if ( LMXF ) then
             msta = mf + 1   ! first internal cell-face
          else
             msta = mf       ! first internal cell-face at subdomain interface
          endif
          mend = ml - 1      ! last  internal cell-face
          !
          do k = 1, kmax
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
          do k = 1, kmax
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
                   ! update total concentration based on higher order correction flux at left and right cell-faces
                   !
                   rp(nm,k,l) = hks(nm,k) * rp(nm,k,l) - dt * rdx * ( flux(nm,k,1) - flux(nmd,k,1) )
                   !
                   ! compute layer-averaged concentration
                   !
                   rp(nm,k,l) = rp(nm,k,l) / hks(nm,k)
                   !
                endif
                !
             enddo
             !
          enddo
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
          do m = mfu, mend
             !
             nm = kgrpnt(m)
             !
             do k = 1, kmax
                !
                mass = mass + hks(nm,k) * dx * rp(nm,k,l)
                !
             enddo
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
end subroutine SwashExpLay1DHtrans
