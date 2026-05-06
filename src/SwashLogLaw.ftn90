subroutine SwashLogLaw
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
!    3.01: Marcel Zijlema
!    6.02: Tom Bogaard
!
!   Updates
!
!    1.00, January 2012: New subroutine
!    3.01, January 2016: effect suspended sediment on stable stratified BBL included
!    6.02,    July 2019: extension 3D (non)linear k-eps model
!
!   Purpose
!
!   Determines friction at beds/walls based on the logarithmic or linear law of the wall
!
!   Method
!
!   The following parameters will be used
!
!   z  = normal distance to the bottom
!   hr = bottom roughness height
!
!   z+ = cmu^(1/4) sqrt(tke) z / nu
!   roughness Reynolds number = cmu^(1/4) sqrt(tke) hr / nu
!
!   It is assumed that there is no transition layer between viscous and turbulent sublayers and so,
!   the edge of viscous sublayer is determined by the following implicit relation
!
!   ev = ln(E*ev) / kappa
!
!   from which it follows that ev = 11.6
!
!   when roughness Reynolds number < 11.6, wall is considered to be smooth, i.e. hr = 0
!
!   if wall is rough (hr /= 0)
!
!      if the location of first computational point is between the roughness elements then z = hr/30
!
!      tau_w = cmu^1/4 kappa sqrt(tke) u / ln(E_rough*z/hr)
!      u+    = ln( E_rough*z/hr )/kappa
!
!   if wall is smooth (hr = 0)
!
!      if the computational point lies in the viscous sublayer then
!
!         tau_w = nu * u / z
!         u+    = z+
!
!      else if the point lies in the logarithmic layer then
!
!         tau_w = cmu^1/4 kappa sqrt(tke) u / ln(E z+)
!         u+    = ln( E z+ )/kappa
!
!   B.E. Launder and D.B. Spalding, "The numerical computation of turbulent
!   flows", Comput. Meth. Appl. Mech. Engng., vol. 3, p. 269-289, 1974
!
!   In very shallow water (depth less than 5 cm) we assume a viscous sublayer
!
!   Various studies have shown that the sediment stratification effects in the bottom boundary layer (BBL)
!   can lead to a change in BBL dynamics. The velocity profile for a stable stratified logarithmic
!   boundary layer is given by Adams and Weatherly (1981). This actually amounts to the reduction of
!   the Von Karman constant, such that the velocity near the bed is enhanced.
!
!   C.E. Adams and G.L. Weatherly, "Some effects of suspended sediment
!   stratification on an oceanic bottom boundary layer", JGR, vol. 86,
!   no. C5, p. 4161-4172, 1981
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use m_genarr
    use m_parall
    use SwashFlowdata
!
    implicit none
!
!   Parameter variables
!
    real,    parameter :: cmu    = 0.09  ! closure constant for standard k-eps model
    real,    parameter :: erough = 33.0  ! empirical constant for logarithmic log-law in case of rough beds
    real,    parameter :: esmoot = 9.0   ! empirical constant for logarithmic log-law in case of smooth beds
    real,    parameter :: ev     = 11.6  ! edge of viscous sublayer
    real,    parameter :: hmin   = 0.05  ! minimal water depth in meters
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter over vertical layers
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mlast    ! index of last internal u-point
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nfm      ! pointer to m,nf
    integer       :: nfum     ! pointer to m,nfu
    integer       :: nfumu    ! pointer to mu,nfu
    integer       :: nlast    ! index of last internal v-point
    integer       :: nlm      ! pointer to m,nl
    integer       :: nlmu     ! pointer to mu,nl
    integer       :: nlum     ! pointer to m,nlu
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmf      ! pointer to mf,n
    integer       :: nmfu     ! pointer to mfu,n
    integer       :: nml      ! pointer to ml,n
    integer       :: nmlu     ! pointer to mlu,n
    integer       :: nu       ! index of point n+1
    integer       :: numfu    ! pointer to mfu,nu
    integer       :: numl     ! pointer to ml,nu
    !
    real          :: drhodz   ! vertical gradient of density
    real          :: dudz     ! vertical gradient of u-velocity
    real          :: dvdz     ! vertical gradient of v-velocity
    real          :: dz       ! local layer thickness
    real          :: fac      ! stability factor
    real          :: hr       ! local roughness height
    real          :: r        ! roughness Reynolds number
    real          :: richn    ! gradient Richardson number
    real          :: shear    ! magnitude of shear squared
    real          :: tke      ! turbulent kinetic energy at bed/wall
    real          :: uplus    ! dimensionless flow velocity
    real          :: ustar    ! friction velocity
    real          :: yplus    ! dimensionless normal distance to wall
    real          :: z0       ! roughness parameter
    real          :: zplus    ! dimensionless normal distance to bed
    real          :: zs       ! distance of half layer thickness to bed
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashLogLaw')
    !
    if ( oned ) then
       !
       if ( varfr ) then
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m,1)
             !
             zs = 0.5 * ( zks(nm,kmax-1) - zks(nm,kmax) )
             !
             if ( hs(nm) > hmin ) then
                !
                tke = rtur(nm,kmax,1)
                !
                ustar = ( cmu**.25 ) * sqrt(tke)
                !
                zplus = ustar * zs          / kinvis
                r     = ustar * fricf(nm,2) / kinvis
                !
                if ( r > ev ) then
                   hr = fricf(nm,2)
                else
                   hr = 0.
                endif
                !
                z0 = hr / erough
                !
                uplus = zplus
                !
                if ( z0 /= 0. ) then
                   !
                   if ( .not. zs > z0 ) zs = zs + z0   ! the bottom is shifted to z = z0
                   !
                   uplus = log( zs/z0 ) / vonkar
                   !
                else
                   !
                   if ( zplus > ev ) uplus = log( esmoot*zplus ) / vonkar
                   !
                endif
                !
             else
                !
                ustar = 1.                   ! dummy value
                zplus = ustar * zs / kinvis
                uplus = zplus
                !
             endif
             !
             logfrc(nm,1) = ustar / uplus
             logfrc(nm,2) = uplus
             !
          enddo
          !
       else
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m,1)
             !
             zs = 0.5 * ( zks(nm,kmax-1) - zks(nm,kmax) )
             !
             if ( hs(nm) > hmin ) then
                !
                tke = rtur(nm,kmax,1)
                !
                ustar = ( cmu**.25 ) * sqrt(tke)
                !
                zplus = ustar * zs      / kinvis
                r     = ustar * pbot(2) / kinvis
                !
                if ( r > ev ) then
                   hr = pbot(2)
                else
                   hr = 0.
                endif
                !
                z0 = hr / erough
                !
                uplus = zplus
                !
                if ( z0 /= 0. ) then
                   !
                   if ( .not. zs > z0 ) zs = zs + z0   ! the bottom is shifted to z = z0
                   !
                   uplus = log( zs/z0 ) / vonkar
                   !
                else
                   !
                   if ( zplus > ev ) uplus = log( esmoot*zplus ) / vonkar
                   !
                endif
                !
             else
                !
                ustar = 1.                   ! dummy value
                zplus = ustar * zs / kinvis
                uplus = zplus
                !
             endif
             !
             logfrc(nm,1) = ustar / uplus
             logfrc(nm,2) = uplus
             !
          enddo
          !
       endif
       !
       if ( lsed > 0 .and. psed(13) > 0. ) then
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,1)
             nmd = kgrpnt(md,1)
             !
             if ( hs(nm) > hmin ) then
                !
                dz = 0.5 * ( hks(nm,kmax-1) + hks(nm,kmax) )
                !
                ! compute vertical shear squared near the bed
                !
                dudz = 0.5 * ( u0(nm,kmax-1) + u0(nmd,kmax-1) - u0(nm,kmax) - u0(nmd,kmax) ) / dz
                !
                shear = max ( dudz*dudz, 1.e-8 )
                !
                ! compute vertical gradient of density near the bed
                !
                drhodz = ( rho(nm,kmax-1) - rho(nm,kmax) ) / dz
                !
                ! compute gradient Richardson number near the bed
                !
                richn = -grav * drhodz / ( shear*rhow )
                !
                ! compute stability factor to reduce the bed shear stress due to the BBL stratification induced by resuspended load
                !
                fac = psed(13) * max( 0., min( 0.2, richn ) )
                !
                logfrc(nm,1) = logfrc(nm,1) / ( 1. + fac )
                logfrc(nm,2) = logfrc(nm,2) * ( 1. + fac )
                !
             endif
             !
          enddo
          !
       endif
       !
       ! exchange coefficients with neighbouring subdomains
       !
       call SWEXCHG ( logfrc(1,1), kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( logfrc(1,2), kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! at real boundaries, copy them from internal points
       !
       nmf  = kgrpnt(mf ,1)
       nmfu = kgrpnt(mfu,1)
       nml  = kgrpnt(ml ,1)
       nmlu = kgrpnt(mlu,1)
       !
       if ( LMXF ) then
          logfrc(nmf,1) = logfrc(nmfu,1)
          logfrc(nmf,2) = logfrc(nmfu,2)
       endif
       if ( LMXL ) then
          logfrc(nmlu,1) = logfrc(nml,1)
          logfrc(nmlu,2) = logfrc(nml,2)
       endif
       !
    else
       !
       if ( varfr ) then
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                zs = 0.5 * ( zks(nm,kmax-1) - zks(nm,kmax) )
                !
                if ( hs(nm) > hmin ) then
                   !
                   tke = rtur(nm,kmax,1)
                   !
                   ustar = ( cmu**.25 ) * sqrt(tke)
                   !
                   zplus = ustar * zs          / kinvis
                   r     = ustar * fricf(nm,2) / kinvis
                   !
                   if ( r > ev ) then
                      hr = fricf(nm,2)
                   else
                      hr = 0.
                   endif
                   !
                   z0 = hr / erough
                   !
                   uplus = zplus
                   !
                   if ( z0 /= 0. ) then
                      !
                      if ( .not. zs > z0 ) zs = zs + z0   ! the bottom is shifted to z = z0
                      !
                      uplus = log( zs/z0 ) / vonkar
                      !
                   else
                      !
                      if ( zplus > ev ) uplus = log( esmoot*zplus ) / vonkar
                      !
                   endif
                   !
                else
                   !
                   ustar = 1.                   ! dummy value
                   zplus = ustar * zs / kinvis
                   uplus = zplus
                   !
                endif
                !
                logfrc(nm,1) = ustar / uplus
                logfrc(nm,2) = uplus
                !
             enddo
          enddo
          !
       else
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                zs = 0.5 * ( zks(nm,kmax-1) - zks(nm,kmax) )
                !
                if ( hs(nm) > hmin ) then
                   !
                   tke = rtur(nm,kmax,1)
                   !
                   ustar = ( cmu**.25 ) * sqrt(tke)
                   !
                   zplus = ustar * zs      / kinvis
                   r     = ustar * pbot(2) / kinvis
                   !
                   if ( r > ev ) then
                      hr = pbot(2)
                   else
                      hr = 0.
                   endif
                   !
                   z0 = hr / erough
                   !
                   uplus = zplus
                   !
                   if ( z0 /= 0. ) then
                      !
                      if ( .not. zs > z0 ) zs = zs + z0   ! the bottom is shifted to z = z0
                      !
                      uplus = log( zs/z0 ) / vonkar
                      !
                   else
                      !
                      if ( zplus > ev ) uplus = log( esmoot*zplus ) / vonkar
                      !
                   endif
                   !
                else
                   !
                   ustar = 1.                   ! dummy value
                   zplus = ustar * zs / kinvis
                   uplus = zplus
                   !
                endif
                !
                logfrc(nm,1) = ustar / uplus
                logfrc(nm,2) = uplus
                !
             enddo
          enddo
       endif
       !
       if ( lsed > 0 .and. psed(13) > 0. ) then
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
                   dz = 0.5 * ( hks(nm,kmax-1) + hks(nm,kmax) )
                   !
                   ! compute vertical shear squared near the bed
                   !
                   dudz = 0.5 * ( u0(nm,kmax-1) + u0(nmd,kmax-1) - u0(nm,kmax) - u0(nmd,kmax) ) / dz
                   dvdz = 0.5 * ( v0(nm,kmax-1) + v0(ndm,kmax-1) - v0(nm,kmax) - v0(ndm,kmax) ) / dz
                   !
                   shear = max ( dudz*dudz + dvdz*dvdz, 1.e-8 )
                   !
                   ! compute vertical gradient of density near the bed
                   !
                   drhodz = ( rho(nm,kmax-1) - rho(nm,kmax) ) / dz
                   !
                   ! compute gradient Richardson number near the bed
                   !
                   richn = -grav * drhodz / ( shear*rhow )
                   !
                   ! compute stability factor to reduce the bed shear stress due to the BBL stratification induced by resuspended load
                   !
                   fac = psed(13) * max( 0., min( 0.2, richn ) )
                   !
                   logfrc(nm,1) = logfrc(nm,1) / ( 1. + fac )
                   logfrc(nm,2) = logfrc(nm,2) * ( 1. + fac )
                   !
                endif
                !
             enddo
          enddo
          !
       endif
       !
       ! set to zero for permanently dry points
       !
       logfrc(1,1) = 0.
       logfrc(1,2) = 0.
       !
       ! exchange coefficients with neighbouring subdomains
       !
       call SWEXCHG ( logfrc(1,1), kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( logfrc(1,2), kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! at real boundaries, copy them from internal points or apply periodicity
       !
       if ( .not.lreptx ) then
          !
          do n = nfu, nl
             !
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             nml  = kgrpnt(ml ,n)
             nmlu = kgrpnt(mlu,n)
             !
             if ( LMXF ) then
                !
                logfrc(nmf,1) = logfrc(nmfu,1)
                logfrc(nmf,2) = logfrc(nmfu,2)
                !
             endif
             !
             if ( LMXL ) then
                !
                logfrc(nmlu,1) = logfrc(nml,1)
                logfrc(nmlu,2) = logfrc(nml,2)
                !
             endif
             !
          enddo
          !
       else
          !
          call periodic ( logfrc(1,1), kgrpnt, 1, 1 )
          call periodic ( logfrc(1,2), kgrpnt, 1, 1 )
          !
       endif
       !
       if ( .not.lrepty ) then
          !
          do m = mfu, ml
             !
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             nlm  = kgrpnt(m,nl )
             nlum = kgrpnt(m,nlu)
             !
             if ( LMYF ) then
                !
                logfrc(nfm,1) = logfrc(nfum,1)
                logfrc(nfm,2) = logfrc(nfum,2)
                !
             endif
             !
             if ( LMYL ) then
                !
                logfrc(nlum,1) = logfrc(nlm,1)
                logfrc(nlum,2) = logfrc(nlm,2)
                !
             endif
             !
          enddo
          !
       else
          !
          call periodic ( logfrc(1,1), kgrpnt, 1, 1 )
          call periodic ( logfrc(1,2), kgrpnt, 1, 1 )
          !
       endif
       !
    endif
    !
    ! compute wall shear stresses at closed boundaries due to full 3D turbulence
    ! note: walls are smooth
    !
    if ( .not.oned .and. iturb > 1 ) then
       !
       if ( lreptx ) then
          mlast = ml               ! last  internal u-point in case of repeating grid
       else
          mlast = ml - 1           ! last  internal u-point
       endif
       !
       if ( lrepty ) then
          nlast = nl               ! last  internal v-point in case of repeating grid
       else
          nlast = nl - 1           ! last  internal v-point
       endif
       !
       if ( .not.lrepty ) then
          !
          do m = mf+1, mlast
             !
             mu = m + 1
             !
             nfum  = kgrpnt(m ,nfu)
             nfumu = kgrpnt(mu,nfu)
             nlm   = kgrpnt(m ,nl )
             nlmu  = kgrpnt(mu,nl )
             !
             if ( ibb(m) == 1 .and. LMYF ) then
                !
                do k = 1, kmax
                   !
                   tke = 0.25 * ( rtur(nfum,k-1,1) + rtur(nfumu,k-1,1) + rtur(nfum,k,1) + rtur(nfumu,k,1) )
                   !
                   ustar = ( cmu**.25 ) * sqrt(tke)
                   !
                   yplus = 0.5 * guu(nfum) * ustar / kinvis
                   !
                   uplus = log( esmoot*yplus ) / vonkar
                   !
                   logwlb(m,k) = ustar / uplus
                   !
                enddo
                !
             endif
             !
             if ( ibt(m) == 1 .and. LMYL ) then
                !
                do k = 1, kmax
                   !
                   tke = 0.25 * ( rtur(nlm,k-1,1) + rtur(nlmu,k-1,1) + rtur(nlm,k,1) + rtur(nlmu,k,1) )
                   !
                   ustar = ( cmu**.25 ) * sqrt(tke)
                   !
                   yplus = 0.5 * guu(nlm) * ustar / kinvis
                   !
                   uplus = log( esmoot*yplus ) / vonkar
                   !
                   logwlt(m,k) = ustar / uplus
                   !
                enddo
                !
             endif
             !
          enddo
          !
       endif
       !
       if ( .not.lreptx ) then
          !
          do n = nf+1, nlast
             !
             nu = n + 1
             !
             nmfu  = kgrpnt(mfu,n )
             numfu = kgrpnt(mfu,nu)
             nml   = kgrpnt(ml ,n )
             numl  = kgrpnt(ml ,nu)
             !
             if ( ibl(n) == 1 .and. LMXF ) then
                !
                do k = 1, kmax
                   !
                   tke = 0.25 * ( rtur(nmfu,k-1,1) + rtur(numfu,k-1,1) + rtur(nmfu,k,1) + rtur(numfu,k,1) )
                   !
                   ustar = ( cmu**.25 ) * sqrt(tke)
                   !
                   yplus = 0.5 * gvv(nmfu) * ustar / kinvis
                   !
                   uplus = log( esmoot*yplus ) / vonkar
                   !
                   logwll(n,k) = ustar / uplus
                   !
                enddo
                !
             endif
             !
             if ( ibr(n) == 1 .and. LMXL ) then
                !
                do k = 1, kmax
                   !
                   tke = 0.25 * ( rtur(nml,k-1,1) + rtur(numl,k-1,1) + rtur(nml,k,1) + rtur(numl,k,1) )
                   !
                   ustar = ( cmu**.25 ) * sqrt(tke)
                   !
                   yplus = 0.5 * gvv(nml) * ustar / kinvis
                   !
                   uplus = log( esmoot*yplus ) / vonkar
                   !
                   logwlr(n,k) = ustar / uplus
                   !
                enddo
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
end subroutine SwashLogLaw
