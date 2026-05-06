subroutine SwashWindStress
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
!    1.00, March 2010: New subroutine
!
!   Purpose
!
!   Calculates wind stresses
!
!   Method
!
!   The so-called Charnock drag coefficient formulation is proposed by Charnock (1955).
!   He assumed a logarithmic wind velocity profile in the turbulent layer above the free surface:
!
!     w10(z)     1        z
!     ------ = ----- ln (---)
!       u*     kappa      z0
!
!   where u* is the friction velocity, kappa is the Von Karman constant, z is the vertical height
!   above the free surface and z0 is the roughness height:
!
!     z0 = b * u*^2/g
!
!   with b the dimensionless Charnock coefficient and g the gravity acceleration.
!
!   Since, Cd = u*^2/w10^2 we have an implicit relation between drag coefficient Cd and wind speed w10.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use SwashTimecomm, only: dt
    use m_genarr
    use SwashFlowdata
!
    implicit none
!
!   Parameter variables
!
    integer, parameter :: maxnit = 100  ! maximum number of iterations
    !
    real   , parameter :: eps  =  0.01  ! convergence criterion
    real   , parameter :: pp   =  0.55  ! first coefficient of 2nd order polynomial fit
    real   , parameter :: qq   =  2.97  ! second coefficient of 2nd order polynomial fit
    real   , parameter :: rr   = -1.49  ! third coefficient of 2nd order polynomial fit
    real   , parameter :: uref = 31.5   ! reference wind speed
!
!   Local variables
!
    integer            :: i            ! loop counter
    integer, save      :: ient = 0     ! number of entries in this subroutine
    integer            :: m            ! loop counter
    integer            :: n            ! loop counter
    integer            :: nit          ! number of iterations
    integer            :: md           ! index of point m-1
    integer            :: mu           ! index of point m+1
    integer            :: nd           ! index of point n-1
    integer            :: ndm          ! pointer to m,n-1
    integer            :: ndmu         ! pointer to m+1,n-1
    integer            :: nm           ! pointer to m,n
    integer            :: nmd          ! pointer to m-1,n
    integer            :: nmu          ! pointer to m+1,n
    integer            :: nu           ! index of point n+1
    integer            :: num          ! pointer to m,n+1
    integer            :: numd         ! pointer to m-1,n+1
    integer            :: numu         ! pointer to m+1,n+1
    !
    real               :: alpha        ! correction factor for wind at 10 metres height to wind at surface
    real               :: angle        ! local wave direction
    real               :: cd           ! drag coefficient
    real               :: dsdt         ! vertical velocity of the free surface
    real               :: dsdx         ! slope of water level in x-direction
    real               :: dsdy         ! slope of water level in y-direction
    real               :: fac          ! factor
    real               :: rwx          ! relative wind u-velocity
    real               :: rwy          ! relative wind v-velocity
    real               :: s            ! magnitude u/ustar
    real               :: slope        ! magnitude of slope of water level
    real               :: sold         ! ratio u/ustar at previous iteration
    real               :: utl          ! ratio u10/uref
    real               :: w10          ! magnitude of local wind velocity
    real               :: wcel         ! local wave celerity
    real               :: wcrstp       ! percentage of height to apply wind stress on wave crest
    real               :: wfac         ! wind stress multiplication factor
    real               :: wspeed       ! local wind speed
    real               :: wxc          ! constant wind velocity component in x-direction
    real               :: wyc          ! constant wind velocity component in y-direction
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashWindStress')
    !
    alpha  = pwnd(4)
    wcrstp = pwnd(10)
    !
    ! compute wind stress multiplication factor
    !
    wfac = rhoa / ( rhow * alpha * alpha )
    !
    ! compute constant wind velocity components
    !
    wxc = u10 * cos(wdic)
    wyc = u10 * sin(wdic)
    !
    ! determine wind u-velocity in u-points
    !
    if ( oned ) then
       !
       ! determine wind u-velocity in u-points in 1D computational grid
       !
       if ( varwi ) then
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             windu(nm) = wxf(nm,2)
             !
          enddo
          !
       else
          !
          windu = wxc
          !
       endif
       !
       ! adjust wind velocity by wave celerity, if appropriate
       !
       if ( relwav ) then
          !
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             dsdx = ( s0(nmu) - s0(nm) ) / dx
             fac  = max( abs(dsdx), 1.e-6 )
             !
             dsdt = 0.5 * ( s1(nmu) + s1(nm) - so(nmu) - so(nm) ) / dt
             !
             wcel = min( abs(dsdt) / fac, sqrt( grav * hum(nm) ) )
             if ( dsdx < 0. ) wcel = -wcel
             !
             windu(nm) = windu(nm) - wcel
             !
          enddo
          !
       endif
       !
    else
       !
       ! determine wind u-velocity in u-points in 2D computational grid
       !
       if ( lstag(5) ) then
          !
          do n = nfu, nl
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,n)
                !
                windu(nm) = wxf(nm,2)
                !
             enddo
             !
          enddo
          !
       else if ( varwi ) then
          !
          do n = nfu, nl
             !
             do m = mf, ml
                !
                nd = n - 1
                !
                nm  = kgrpnt(m,n )
                ndm = kgrpnt(m,nd)
                !
                windu(nm) = 0.5 * ( wxf(nm,2) + wxf(ndm,2) ) * ( ycgrid(m,n) - ycgrid(m,n-1) ) - 0.5 * ( wyf(nm,2) + wyf(ndm,2) ) * ( xcgrid(m,n) - xcgrid(m,n-1) )
                windu(nm) = windu(nm) / guu(nm)
                !
             enddo
             !
          enddo
          !
       else
          !
          do n = nfu, nl
             !
             do m = mf, ml
                !
                nd = n - 1
                !
                nm  = kgrpnt(m,n )
                ndm = kgrpnt(m,nd)
                !
                windu(nm) = wxc * ( ycgrid(m,n) - ycgrid(m,n-1) ) - wyc * ( xcgrid(m,n) - xcgrid(m,n-1) )
                windu(nm) = windu(nm) / guu(nm)
                !
             enddo
             !
          enddo
          !
       endif
       !
       ! determine wind v-velocity in v-points in 2D computational grid
       !
       if ( lstag(6) ) then
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                windv(nm) = wyf(nm,2)
                !
             enddo
             !
          enddo
          !
       else if ( varwi ) then
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m ,n)
                nmd = kgrpnt(md,n)
                !
                windv(nm) = 0.5 * ( wyf(nm,2) + wyf(nmd,2) ) * ( xcgrid(m,n) - xcgrid(m-1,n) ) - 0.5 * ( wxf(nm,2) + wxf(nmd,2) ) * ( ycgrid(m,n) - ycgrid(m-1,n) )
                windv(nm) = windv(nm) / gvv(nm)
                !
             enddo
             !
          enddo
          !
       else
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m ,n)
                nmd = kgrpnt(md,n)
                !
                windv(nm) = wyc * ( xcgrid(m,n) - xcgrid(m-1,n) ) - wxc * ( ycgrid(m,n) - ycgrid(m-1,n) )
                windv(nm) = windv(nm) / gvv(nm)
                !
             enddo
             !
          enddo
          !
       endif
       !
       ! adjust wind velocity by wave celerity, if appropriate
       !
       if ( relwav ) then
          !
          do n = nfu, nl
             !
             do m = mf, ml
                !
                mu = m + 1
                nd = n - 1
                nu = n + 1
                !
                nm   = kgrpnt(m ,n )
                nmu  = kgrpnt(mu,n )
                ndm  = kgrpnt(m ,nd)
                ndmu = kgrpnt(mu,nd)
                num  = kgrpnt(m ,nu)
                numu = kgrpnt(mu,nu)
                !
                dsdx = ( s0(nmu) - s0(nm) ) / gvu(nm)
                dsdy = 0.25*( s0(num) + s0(numu) - s0(ndm) - s0(ndmu) ) / guu(nm)
                !
                slope = sqrt ( dsdx**2 + dsdy**2 )
                fac   = max( slope, 1.e-6 )
                !
                dsdt = 0.5 * ( s1(nmu) + s1(nm) - so(nmu) - so(nm) ) / dt
                !
                wcel  = min( abs(dsdt) / fac, sqrt( grav * hum(nm) ) )
                angle = atan2( dsdy, dsdx )
                !
                windu(nm) = windu(nm) - wcel * cos(angle)
                !
             enddo
             !
          enddo
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                mu = m + 1
                nu = n + 1
                !
                nm   = kgrpnt(m ,n )
                num  = kgrpnt(m ,nu)
                nmd  = kgrpnt(md,n )
                numd = kgrpnt(md,nu)
                nmu  = kgrpnt(mu,n )
                numu = kgrpnt(mu,nu)
                !
                dsdx = 0.25*( s0(nmu) + s0(numu) - s0(nmd) - s0(numd) ) / gvv(nm)
                dsdy = ( s0(num) - s0(nm) ) / guv(nm)
                !
                slope = sqrt ( dsdx**2 + dsdy**2 )
                fac   = max( slope, 1.e-6 )
                !
                dsdt = 0.5 * ( s1(num) + s1(nm) - so(num) - so(nm) ) / dt
                !
                wcel  = min( abs(dsdt) / fac, sqrt( grav * hvm(nm) ) )
                angle = atan2( dsdy, dsdx )
                !
                windv(nm) = windv(nm) - wcel * sin(angle)
                !
             enddo
             !
          enddo
          !
       endif
       !
       ! synchronize wind u- and v-velocities at appropriate boundaries in case of repeating grid
       !
       call periodic ( windu, kgrpnt, 1, 1 )
       call periodic ( windv, kgrpnt, 1, 1 )
       !
    endif
    !
    ! compute wind stress coefficient depending on drag coefficient and (relative) wind velocity
    !
    if ( oned ) then
       !
       if ( iwind == 1 ) then
          !
          ! constant drag coefficient
          !
          cd = pwnd(1)
          cd = min( cdcap, cd )
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             if ( relwnd ) then
                !
                w10 = abs(alpha*windu(nm) - u0(nm,1))
                !
             else
                !
                w10 = abs(windu(nm))
                !
             endif
             !
             cwndu(nm) = wfac * cd * w10
             !
          enddo
          !
       else if ( iwind == 2 ) then
          !
          ! drag coefficient based on Charnock drag formula
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             if ( relwnd ) then
                !
                w10 = abs(alpha*windu(nm) - u0(nm,1))
                !
             else
                !
                w10 = abs(windu(nm))
                !
             endif
             !
             nit = 0
             !
             ! initial value for s
             !
             s    = 22.4
             sold = 0.
             !
             ! Newton-Raphson iteration
             !
             do
                if ( abs(sold-s) < (eps*s) ) exit
                !
                nit  = nit + 1
                sold = s
                s    = sold*(log(pwnd(3)*grav*sold*sold / (max(0.001,pwnd(2)*w10*w10)))-2.) / (vonkar*sold-2.)
                !
                if ( .not. nit < maxnit ) then
                   call msgerr (1, 'no convergence in Charnock drag computation')
                   s = 22.4
                   exit
                endif
                !
             enddo
             !
             cd = 1./(s*s)
             cd = min( cdcap, cd )
             !
             cwndu(nm) = wfac * cd * w10
             !
          enddo
          !
       else if ( iwind == 3 ) then
          !
          ! drag coefficient linear dependent on wind speed
          !
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             if ( relwnd ) then
                !
                w10 = abs(alpha*windu(nm) - u0(nm,1))
                !
             else
                !
                w10 = abs(windu(nm))
                !
             endif
             !
             dsdx   = ( s0(nmu) - s0(nm) ) / dx
             wspeed = max(pwnd(8), min(w10,pwnd(9)))
             !
             cd = 0.001*( pwnd(5) + pwnd(6)*abs(dsdx) + pwnd(7)*wspeed )
             cd = min( cdcap, cd )
             !
             cwndu(nm) = wfac * cd * w10
             !
          enddo
          !
       else if ( iwind == 4 ) then
          !
          ! drag coefficient based on 2nd order polynomial fit
          !
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             if ( relwnd ) then
                !
                w10 = abs(alpha*windu(nm) - u0(nm,1))
                !
             else
                !
                w10 = abs(windu(nm))
                !
             endif
             !
             utl = w10/uref
             !
             cd = 0.001*( pp + qq*utl + rr*utl*utl )
             cd = min( cdcap, cd )
             !
             cwndu(nm) = wfac * cd * w10
             !
          enddo
          !
       else
          !
          call msgerr ( 4, 'unknown drag formula for wind stress' )
          return
          !
       endif
       !
    else
       !
       if ( iwind == 1 ) then
          !
          ! constant drag coefficient
          !
          cd = pwnd(1)
          cd = min( cdcap, cd )
          !
          do n = nfu, nl
             !
             do m = mf, ml
                !
                mu = m + 1
                nd = n - 1
                !
                nm   = kgrpnt(m ,n )
                nmu  = kgrpnt(mu,n )
                ndm  = kgrpnt(m ,nd)
                ndmu = kgrpnt(mu,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmu  == 1 ) nmu  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmu == 1 ) ndmu = ndm
                !
                if ( relwnd ) then
                   !
                   rwx = alpha*windu(nm) - u0(nm,1)
                   rwy = 0.25*alpha*( windv(nm) + windv(ndm) + windv(nmu) + windv(ndmu) ) - 0.25*( v0(nm,1) + v0(ndm,1) + v0(nmu,1) + v0(ndmu,1) )
                   !
                else
                   !
                   rwx = windu(nm)
                   rwy = 0.25*( windv(nm) + windv(ndm) + windv(nmu) + windv(ndmu) )
                   !
                endif
                !
                w10 = sqrt( rwx**2 + rwy**2 )
                !
                cwndu(nm) = wfac * cd * w10
                !
             enddo
             !
          enddo
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                nu = n + 1
                !
                nm   = kgrpnt(m ,n )
                num  = kgrpnt(m ,nu)
                nmd  = kgrpnt(md,n )
                numd = kgrpnt(md,nu)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( num  == 1 ) num  = nm
                if ( nmd  == 1 ) nmd  = nm
                if ( numd == 1 ) numd = nmd
                !
                if ( relwnd ) then
                   !
                   rwx = 0.25*alpha*( windu(nm) + windu(nmd) + windu(num) + windu(numd) ) - 0.25*( u0(nm,1) + u0(nmd,1) + u0(num,1) + u0(numd,1) )
                   rwy = alpha*windv(nm) - v0(nm,1)
                   !
                else
                   !
                   rwx = 0.25*( windu(nm) + windu(nmd) + windu(num) + windu(numd) )
                   rwy = windv(nm)
                   !
                endif
                !
                w10 = sqrt( rwx**2 + rwy**2 )
                !
                cwndv(nm) = wfac * cd * w10
                !
             enddo
             !
          enddo
          !
       else if ( iwind == 2 ) then
          !
          ! drag coefficient based on Charnock drag formula
          !
          do n = nfu, nl
             !
             do m = mf, ml
                !
                mu = m + 1
                nd = n - 1
                !
                nm   = kgrpnt(m ,n )
                nmu  = kgrpnt(mu,n )
                ndm  = kgrpnt(m ,nd)
                ndmu = kgrpnt(mu,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmu  == 1 ) nmu  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmu == 1 ) ndmu = ndm
                !
                if ( relwnd ) then
                   !
                   rwx = alpha*windu(nm) - u0(nm,1)
                   rwy = 0.25*alpha*( windv(nm) + windv(ndm) + windv(nmu) + windv(ndmu) ) - 0.25*( v0(nm,1) + v0(ndm,1) + v0(nmu,1) + v0(ndmu,1) )
                   !
                else
                   !
                   rwx = windu(nm)
                   rwy = 0.25*( windv(nm) + windv(ndm) + windv(nmu) + windv(ndmu) )
                   !
                endif
                !
                w10 = sqrt( rwx**2 + rwy**2 )
                !
                nit = 0
                !
                ! initial value for s
                !
                s    = 22.4
                sold = 0.
                !
                ! Newton-Raphson iteration
                !
                do
                   if ( abs(sold-s) < (eps*s) ) exit
                   !
                   nit  = nit + 1
                   sold = s
                   s    = sold*(log(pwnd(3)*grav*sold*sold / (max(0.001,pwnd(2)*w10*w10)))-2.) / (vonkar*sold-2.)
                   !
                   if ( .not. nit < maxnit ) then
                      call msgerr (1, 'no convergence in Charnock drag computation')
                      s = 22.4
                      exit
                   endif
                   !
                enddo
                !
                cd = 1./(s*s)
                cd = min( cdcap, cd )
                !
                cwndu(nm) = wfac * cd * w10
                !
             enddo
             !
          enddo
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                nu = n + 1
                !
                nm   = kgrpnt(m ,n )
                num  = kgrpnt(m ,nu)
                nmd  = kgrpnt(md,n )
                numd = kgrpnt(md,nu)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( num  == 1 ) num  = nm
                if ( nmd  == 1 ) nmd  = nm
                if ( numd == 1 ) numd = nmd
                !
                if ( relwnd ) then
                   !
                   rwx = 0.25*alpha*( windu(nm) + windu(nmd) + windu(num) + windu(numd) ) - 0.25*( u0(nm,1) + u0(nmd,1) + u0(num,1) + u0(numd,1) )
                   rwy = alpha*windv(nm) - v0(nm,1)
                   !
                else
                   !
                   rwx = 0.25*( windu(nm) + windu(nmd) + windu(num) + windu(numd) )
                   rwy = windv(nm)
                   !
                endif
                !
                w10 = sqrt( rwx**2 + rwy**2 )
                !
                nit = 0
                !
                ! initial value for s
                !
                s    = 22.4
                sold = 0.
                !
                ! Newton-Raphson iteration
                !
                do
                   if ( abs(sold-s) < (eps*s) ) exit
                   !
                   nit  = nit + 1
                   sold = s
                   s    = sold*(log(pwnd(3)*grav*sold*sold / (max(0.001,pwnd(2)*w10*w10)))-2.) / (vonkar*sold-2.)
                   !
                   if ( .not. nit < maxnit ) then
                      call msgerr (1, 'no convergence in Charnock drag computation')
                      s = 22.4
                      exit
                   endif
                   !
                enddo
                !
                cd = 1./(s*s)
                cd = min( cdcap, cd )
                !
                cwndv(nm) = wfac * cd * w10
                !
             enddo
             !
          enddo
          !
       else if ( iwind == 3 ) then
          !
          ! drag coefficient linear dependent on wind speed
          !
          do n = nfu, nl
             !
             do m = mf, ml
                !
                mu = m + 1
                nd = n - 1
                nu = n + 1
                !
                nm   = kgrpnt(m ,n )
                nmu  = kgrpnt(mu,n )
                ndm  = kgrpnt(m ,nd)
                ndmu = kgrpnt(mu,nd)
                num  = kgrpnt(m ,nu)
                numu = kgrpnt(mu,nu)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmu  == 1 ) nmu  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmu == 1 ) ndmu = ndm
                if ( num  == 1 ) num  = nm
                if ( numu == 1 ) numu = num
                !
                if ( relwnd ) then
                   !
                   rwx = alpha*windu(nm) - u0(nm,1)
                   rwy = 0.25*alpha*( windv(nm) + windv(ndm) + windv(nmu) + windv(ndmu) ) - 0.25*( v0(nm,1) + v0(ndm,1) + v0(nmu,1) + v0(ndmu,1) )
                   !
                else
                   !
                   rwx = windu(nm)
                   rwy = 0.25*( windv(nm) + windv(ndm) + windv(nmu) + windv(ndmu) )
                   !
                endif
                !
                w10 = sqrt( rwx**2 + rwy**2 )
                !
                dsdx = ( s0(nmu) - s0(nm) ) / gvu(nm)
                dsdy = 0.25*( s0(num) + s0(numu) - s0(ndm) - s0(ndmu) ) / guu(nm)
                !
                slope  = sqrt ( dsdx**2 + dsdy**2 )
                wspeed = max(pwnd(8), min(w10,pwnd(9)))
                !
                cd = 0.001*( pwnd(5) + pwnd(6)*abs(slope) + pwnd(7)*wspeed )
                cd = min( cdcap, cd )
                !
                cwndu(nm) = wfac * cd * w10
                !
             enddo
             !
          enddo
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                mu = m + 1
                nu = n + 1
                !
                nm   = kgrpnt(m ,n )
                num  = kgrpnt(m ,nu)
                nmd  = kgrpnt(md,n )
                numd = kgrpnt(md,nu)
                nmu  = kgrpnt(mu,n )
                numu = kgrpnt(mu,nu)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( num  == 1 ) num  = nm
                if ( nmd  == 1 ) nmd  = nm
                if ( numd == 1 ) numd = nmd
                if ( nmu  == 1 ) nmu  = nm
                if ( numu == 1 ) numu = nmu
                !
                if ( relwnd ) then
                   !
                   rwx = 0.25*alpha*( windu(nm) + windu(nmd) + windu(num) + windu(numd) ) - 0.25*( u0(nm,1) + u0(nmd,1) + u0(num,1) + u0(numd,1) )
                   rwy = alpha*windv(nm) - v0(nm,1)
                   !
                else
                   !
                   rwx = 0.25*( windu(nm) + windu(nmd) + windu(num) + windu(numd) )
                   rwy = windv(nm)
                   !
                endif
                !
                w10 = sqrt( rwx**2 + rwy**2 )
                !
                dsdx = 0.25*( s0(nmu) + s0(numu) - s0(nmd) - s0(numd) ) / gvv(nm)
                dsdy = ( s0(num) - s0(nm) ) / guv(nm)
                !
                slope  = sqrt ( dsdx**2 + dsdy**2 )
                wspeed = max(pwnd(8), min(w10,pwnd(9)))
                !
                cd = 0.001*( pwnd(5) + pwnd(6)*abs(slope) + pwnd(7)*wspeed )
                cd = min( cdcap, cd )
                !
                cwndv(nm) = wfac * cd * w10
                !
             enddo
             !
          enddo
          !
       else if ( iwind == 4 ) then
          !
          ! drag coefficient based on 2nd order polynomial fit
          !
          do n = nfu, nl
             !
             do m = mf, ml
                !
                mu = m + 1
                nd = n - 1
                !
                nm   = kgrpnt(m ,n )
                nmu  = kgrpnt(mu,n )
                ndm  = kgrpnt(m ,nd)
                ndmu = kgrpnt(mu,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmu  == 1 ) nmu  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmu == 1 ) ndmu = ndm
                !
                if ( relwnd ) then
                   !
                   rwx = alpha*windu(nm) - u0(nm,1)
                   rwy = 0.25*alpha*( windv(nm) + windv(ndm) + windv(nmu) + windv(ndmu) ) - 0.25*( v0(nm,1) + v0(ndm,1) + v0(nmu,1) + v0(ndmu,1) )
                   !
                else
                   !
                   rwx = windu(nm)
                   rwy = 0.25*( windv(nm) + windv(ndm) + windv(nmu) + windv(ndmu) )
                   !
                endif
                !
                w10 = sqrt( rwx**2 + rwy**2 )
                !
                utl = w10/uref
                !
                cd = 0.001*( pp + qq*utl + rr*utl*utl )
                cd = min( cdcap, cd )
                !
                cwndu(nm) = wfac * cd * w10
                !
             enddo
             !
          enddo
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                nu = n + 1
                !
                nm   = kgrpnt(m ,n )
                num  = kgrpnt(m ,nu)
                nmd  = kgrpnt(md,n )
                numd = kgrpnt(md,nu)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( num  == 1 ) num  = nm
                if ( nmd  == 1 ) nmd  = nm
                if ( numd == 1 ) numd = nmd
                !
                if ( relwnd ) then
                   !
                   rwx = 0.25*alpha*( windu(nm) + windu(nmd) + windu(num) + windu(numd) ) - 0.25*( u0(nm,1) + u0(nmd,1) + u0(num,1) + u0(numd,1) )
                   rwy = alpha*windv(nm) - v0(nm,1)
                   !
                else
                   !
                   rwx = 0.25*( windu(nm) + windu(nmd) + windu(num) + windu(numd) )
                   rwy = windv(nm)
                   !
                endif
                !
                w10 = sqrt( rwx**2 + rwy**2 )
                !
                utl = w10/uref
                !
                cd = 0.001*( pp + qq*utl + rr*utl*utl )
                cd = min( cdcap, cd )
                !
                cwndv(nm) = wfac * cd * w10
                !
             enddo
             !
          enddo
          !
       else
          !
          call msgerr ( 4, 'unknown drag formula for wind stress' )
          return
          !
       endif
       !
       ! synchronize wind stresses at appropriate boundaries in case of repeating grid
       !
       call periodic ( cwndu, kgrpnt, 1, 1 )
       call periodic ( cwndv, kgrpnt, 1, 1 )
       !
    endif
    !
    ! set to zero for permanently dry points
    !
    cwndu(1) = 0.
    if ( .not.oned ) cwndv(1) = 0.
    !
    ! compute wind stress term
    !
    do i = 1, mcgrd
       !
       windu(i) = alpha*cwndu(i)*windu(i)
       !
    enddo
    !
    if ( .not.oned ) then
       !
       do i = 1, mcgrd
          !
          windv(i) = alpha*cwndv(i)*windv(i)
          !
       enddo
       !
    endif
    !
    ! apply wind stress on wave crest only, if appropriate
    !
    if ( relwav ) then
       !
       do i = 1, mcgrd
          !
          if ( s0(i) > smax(i) ) smax(i) = s0(i)
          !
          if ( s0(i) < (1.-wcrstp) * smax(i) ) then
             !
             windu(i) = 0.
             if ( .not.oned ) windv(i) = 0.
             !
          endif
          !
       enddo
       !
    endif
    !
end subroutine SwashWindStress
