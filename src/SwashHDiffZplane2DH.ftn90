subroutine SwashHDiffZplane2DH ( rhsc, rpo, psm, icistb )
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
!    1.00, July 2016: New subroutine
!
!   Purpose
!
!   Computes horizontal diffusion along strictly horizontal planes
!
!   Method
!
!   Horizontal diffusion is computed explicitly along strictly horizontal planes
!   using the MUSCL limiter to avoid artificial vertical upwelling based on
!   the method described in
!
!   On the approximation of horizontal gradients in sigma co-ordinates for
!   bathymetry with steep bottom slopes
!   G.S. Stelling and J.A.Th.M. van Kester, IJMF, vol. 18, 915-935, 1994
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
    use m_genarr, only: kgrpnt, guu, guv, gvu, gvv, gsqs
    use m_parall
    use SwashFlowdata, only: mfu, ml, nfu, nl, hs, vnu2d, zkso
!
    implicit none
!
!   Argument variables
!
    integer                    , intent(inout) :: icistb ! counter for number of instable points
    !
    real                       , intent(in)    :: psm    ! Prandtl-Schmidt number
    real, dimension(mcgrd,kmax), intent(in)    :: rpo    ! constituent in wl-point at previous time level
    real, dimension(mcgrd,kmax), intent(inout) :: rhsc   ! right-hand side of tri-diagonal system of transport equation
!
!   Parameter variables
!
    integer, parameter         :: km = 500 !  dimension integer
!
!   Local variables
!
    integer, save              :: ient = 0 ! number of entries in this subroutine
    integer                    :: k        ! loop counter over vertical layers
    integer                    :: k1       ! auxiliary integer
    integer                    :: k2       ! auxiliary integer
    integer                    :: kd       ! k-index downward of considered flux point
    integer                    :: kf       ! loop counter
    integer                    :: kflux    ! k-index of flux
    integer, dimension(2*km+1) :: kicod    ! k-index of constituent point downward of considered flux
    integer, dimension(2*km+1) :: kicou    ! k-index of constituent point upward of considered flux
    integer                    :: ku       ! k-index upward of considered flux point
    integer                    :: m        ! loop counter
    integer                    :: md       ! index of point m-1
    integer                    :: mend     ! end index of loop over wl-points
    integer                    :: mu       ! index of point m+1
    integer                    :: n        ! loop counter
    integer                    :: nd       ! index of point n-1
    integer                    :: ndm      ! pointer to m,n-1
    integer                    :: nend     ! end index of loop over wl-points in y-direction
    integer                    :: nm       ! pointer to m,n
    integer                    :: nmd      ! pointer to m-1,n
    integer                    :: nmu      ! pointer to m+1,n
    integer                    :: nsta     ! start index of loop over wl-points in y-direction
    integer                    :: nu       ! index of point n+1
    integer                    :: num      ! pointer to m,n+1
    !
    real                       :: cd       ! interpolated constituent downward of considered flux
    real                       :: cu       ! interpolated constituent upward of considered flux
    real                       :: dif2d    ! effective horizontal eddy diffusivity coefficient
    real                       :: dxl      ! local mesh size in x-direction
    real                       :: dyl      ! local mesh size in y-direction
    real                       :: flux     ! total diffusion flux in x/y-direction
    real                       :: grad     ! strictly horizontal gradient of constituent
    real                       :: grad1    ! first strictly horizontal gradient of constituent
    real                       :: grad2    ! second strictly horizontal gradient of constituent
    real                       :: grmax    ! maximum of grad1 and grad2
    real                       :: grmin    ! minimum of grad1 and grad2
    real                       :: grsum    ! sum of grad1 and grad2
    real, dimension(km)        :: pocod    ! z-coordinate in downward wl-point
    real, dimension(km)        :: pocou    ! z-coordinate in upward wl-point
    real, dimension(2*km+1)    :: poflu    ! z-coordinate of gradient flux
    real, dimension(0:2*km+1)  :: point    ! merge point
    real                       :: stabmx   ! auxiliary variable with maximum diffusivity based stability criterion
    !
    logical                    :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashHDiffZplane2DH')
    !
    ! not computed for end point at subdomain interface since, this end point is owned by the neighbouring subdomain
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    ! compute diffusive flux in x-direction at internal cell-faces
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          mu = m + 1
          nd = n - 1
          !
          nm  = kgrpnt(m ,n )
          nmu = kgrpnt(mu,n )
          ndm = kgrpnt(m ,nd)
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
          if ( hs(nm) > epsdry ) then
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmu == 1 ) nmu = nm
             !
             do k = 1, kmax
                pocod(k) = 0.5 * ( zkso( nm,k-1) + zkso( nm,k) )
                pocou(k) = 0.5 * ( zkso(nmu,k-1) + zkso(nmu,k) )
             enddo
             !
             ! merge adjacent layer-interfaces
             !
             kd = 0
             ku = 0
             do k = 0, 2*kmax+1
                if ( zkso(nm,kd) > zkso(nmu,ku) ) then
                   point(k) = zkso(nm,kd)
                   kd       = kd + 1
                   if ( kd > kmax ) then
                      kflux = k + 1
                      point(kflux) = zkso(nmu,ku)
                      goto 10
                   endif
                else
                   point(k) = zkso(nmu,ku)
                   ku       = ku + 1
                   if ( ku > kmax ) then
                      kflux = k + 1
                      point(kflux) = zkso(nm,kd)
                      goto 10
                   endif
                endif
             enddo
             kflux = 2*kmax+1
     10      continue
             !
             ! position flux points
             !
             do k = 1, kflux
                poflu(k) = 0.5 * ( point(k) + point(k-1) )
             enddo
             !
             ! k-index concentration points left and right of flux point
             !
             kd = 1
             ku = 1
             do kf = 1, kflux
                kicod(kf) = 0
                kicou(kf) = 0
                do k = kd, kmax
                   if ( (.not. poflu(kf) > zkso(nm,k-1)) .and. (.not. poflu(kf) < zkso(nm,k)) ) then
                      kicod(kf) = k
                      kd        = k
                      exit
                   endif
                enddo
                do k = ku, kmax
                   if ( (.not. poflu(kf) > zkso(nmu,k-1)) .and. (.not. poflu(kf) < zkso(nmu,k)) ) then
                      kicou(kf) = k
                      ku        = k
                      exit
                   endif
                enddo
             enddo
             !
             ! compute diffusive flux using limiter
             !
             do kf = 1, kflux
                !
                kd = kicod(kf)
                ku = kicou(kf)
                !
                if ( kd /= 0 .and. ku /= 0 ) then
                   !
                   ! interpolation left point
                   !
                   k2 = kd
                   if ( pocou(ku) > pocod(kd) ) then
                      k1 = k2 - 1
                      if ( k1 < 1 ) then
                         cd = rpo(nm,k2)
                      else
                         cd = ((pocou(ku)-pocod(k2))/(pocod(k1)-pocod(k2))) * rpo(nm,k1) + ((pocou(ku)-pocod(k1))/(pocod(k2)-pocod(k1))) * rpo(nm,k2)
                      endif
                   else
                      k1 = k2 + 1
                      if ( k1 > kmax ) then
                         cd = rpo(nm,k2)
                      else
                         cd = ((pocou(ku)-pocod(k2))/(pocod(k1)-pocod(k2))) * rpo(nm,k1) + ((pocou(ku)-pocod(k1))/(pocod(k2)-pocod(k1))) * rpo(nm,k2)
                      endif
                   endif
                   !
                   ! interpolation right point
                   !
                   k2 = ku
                   if ( pocod(kd) > pocou(ku) ) then
                      k1 = k2 - 1
                      if ( k1 < 1 ) then
                         cu = rpo(nmu,k2)
                      else
                         cu = ((pocod(kd)-pocou(k2))/(pocou(k1)-pocou(k2))) * rpo(nmu,k1) + ((pocod(kd)-pocou(k1))/(pocou(k2)-pocou(k1))) * rpo(nmu,k2)
                      endif
                   else
                      k1 = k2 + 1
                      if ( k1 > kmax ) then
                         cu = rpo(nmu,k2)
                      else
                         cu = ((pocod(kd)-pocou(k2))/(pocou(k1)-pocou(k2))) * rpo(nmu,k1) + ((pocod(kd)-pocou(k1))/(pocou(k2)-pocou(k1))) * rpo(nmu,k2)
                      endif
                   endif
                   !
                   grad1 = rpo(nmu,ku) - cd
                   grad2 = cu - rpo(nm,kd)
                   grsum = grad1 + grad2
                   grmax = max( grad1, grad2 )
                   grmin = min( grad1, grad2 )
                   !
                   grad = 0.
                   !
                   ! Van Leer limiter
                   !if ( grmax < 0. .or. grmin > 0. ) grad = 2. * grad1 * grad2 / grsum
                   !
                   ! MUSCL limiter
                   if ( grmax < 0. .or. grmin > 0. ) then
                      grad = min( grmin, 0.25 * grsum )
                      if (grad1 < 0.) grad = max( grmax, 0.25 * grsum )
                   endif
                   grad = 2. * grad
                   !
                   ! diffusion flux in x-direction
                   !
                   flux = dif2d * ( point(kf-1) - point(kf) ) * grad * guu(nm) / gvu(nm)
                   if ( m == ml ) flux = 0.
                   !
                   rhsc( nm,kd) = rhsc( nm,kd) + flux / gsqs(nm )
                   rhsc(nmu,ku) = rhsc(nmu,ku) - flux / gsqs(nmu)
                   !
                endif
                !
             enddo
             !
          endif
          !
       enddo
    enddo
    !
    ! compute diffusive flux in y-direction at internal cell-faces
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          md = m - 1
          nu = n + 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          num = kgrpnt(m ,nu)
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
          if ( hs(nm) > epsdry ) then
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( num == 1 ) num = nm
             !
             do k = 1, kmax
                pocod(k) = 0.5 * ( zkso( nm,k-1) + zkso( nm,k) )
                pocou(k) = 0.5 * ( zkso(num,k-1) + zkso(num,k) )
             enddo
             !
             ! merge adjacent layer-interfaces
             !
             kd = 0
             ku = 0
             do k = 0, 2*kmax+1
                if ( zkso(nm,kd) > zkso(num,ku) ) then
                   point(k) = zkso(nm,kd)
                   kd       = kd + 1
                   if ( kd > kmax ) then
                      kflux = k + 1
                      point(kflux) = zkso(num,ku)
                      goto 20
                   endif
                else
                   point(k) = zkso(num,ku)
                   ku       = ku + 1
                   if ( ku > kmax ) then
                      kflux = k + 1
                      point(kflux) = zkso(nm,kd)
                      goto 20
                   endif
                endif
             enddo
             kflux = 2*kmax+1
     20      continue
             !
             ! position flux points
             !
             do k = 1, kflux
                poflu(k) = 0.5 * ( point(k) + point(k-1) )
             enddo
             !
             ! k-index concentration points bottom and top of flux point
             !
             kd = 1
             ku = 1
             do kf = 1, kflux
                kicod(kf) = 0
                kicou(kf) = 0
                do k = kd, kmax
                   if ( (.not. poflu(kf) > zkso(nm,k-1)) .and. (.not. poflu(kf) < zkso(nm,k)) ) then
                      kicod(kf) = k
                      kd        = k
                      exit
                   endif
                enddo
                do k = ku, kmax
                   if ( (.not. poflu(kf) > zkso(num,k-1)) .and. (.not. poflu(kf) < zkso(num,k)) ) then
                      kicou(kf) = k
                      ku        = k
                      exit
                   endif
                enddo
             enddo
             !
             ! compute diffusive flux using limiter
             !
             do kf = 1, kflux
                !
                kd = kicod(kf)
                ku = kicou(kf)
                !
                if ( kd /= 0 .and. ku /= 0 ) then
                   !
                   ! interpolation bottom point
                   !
                   k2 = kd
                   if ( pocou(ku) > pocod(kd) ) then
                      k1 = k2 - 1
                      if ( k1 < 1 ) then
                         cd = rpo(nm,k2)
                      else
                         cd = ((pocou(ku)-pocod(k2))/(pocod(k1)-pocod(k2))) * rpo(nm,k1) + ((pocou(ku)-pocod(k1))/(pocod(k2)-pocod(k1))) * rpo(nm,k2)
                      endif
                   else
                      k1 = k2 + 1
                      if ( k1 > kmax ) then
                         cd = rpo(nm,k2)
                      else
                         cd = ((pocou(ku)-pocod(k2))/(pocod(k1)-pocod(k2))) * rpo(nm,k1) + ((pocou(ku)-pocod(k1))/(pocod(k2)-pocod(k1))) * rpo(nm,k2)
                      endif
                   endif
                   !
                   ! interpolation top point
                   !
                   k2 = ku
                   if ( pocod(kd) > pocou(ku) ) then
                      k1 = k2 - 1
                      if ( k1 < 1 ) then
                         cu = rpo(num,k2)
                      else
                         cu = ((pocod(kd)-pocou(k2))/(pocou(k1)-pocou(k2))) * rpo(num,k1) + ((pocod(kd)-pocou(k1))/(pocou(k2)-pocou(k1))) * rpo(num,k2)
                      endif
                   else
                      k1 = k2 + 1
                      if ( k1 > kmax ) then
                         cu = rpo(num,k2)
                      else
                         cu = ((pocod(kd)-pocou(k2))/(pocou(k1)-pocou(k2))) * rpo(num,k1) + ((pocod(kd)-pocou(k1))/(pocou(k2)-pocou(k1))) * rpo(num,k2)
                      endif
                   endif
                   !
                   grad1 = rpo(num,ku) - cd
                   grad2 = cu - rpo(nm,kd)
                   grsum = grad1 + grad2
                   grmax = max( grad1, grad2 )
                   grmin = min( grad1, grad2 )
                   !
                   grad = 0.
                   !
                   ! Van Leer limiter
                   !if ( grmax < 0. .or. grmin > 0. ) grad = 2. * grad1 * grad2 / grsum
                   !
                   ! MUSCL limiter
                   if ( grmax < 0. .or. grmin > 0. ) then
                      grad = min( grmin, 0.25 * grsum )
                      if (grad1 < 0.) grad = max( grmax, 0.25 * grsum )
                   endif
                   grad = 2. * grad
                   !
                   ! diffusion flux in y-direction
                   !
                   flux = dif2d * ( point(kf-1) - point(kf) ) * grad * gvv(nm) / guv(nm)
                   if ( n == nl ) flux = 0.
                   !
                   rhsc( nm,kd) = rhsc( nm,kd) + flux / gsqs(nm )
                   rhsc(num,ku) = rhsc(num,ku) - flux / gsqs(num)
                   !
                endif
                !
             enddo
             !
          endif
          !
       enddo
    enddo
    !
end subroutine SwashHDiffZplane2DH
