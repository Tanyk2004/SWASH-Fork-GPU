subroutine SwashHDiffZplane1DH ( rhsc, rpo, kgrpnt, psm, icistb )
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
    use m_parall
    use SwashFlowdata, only: mfu, ml, hs, vnu2d, zkso
!
    implicit none
!
!   Argument variables
!
    integer                    , intent(inout) :: icistb ! counter for number of instable points
    integer, dimension(mxc)    , intent(in)    :: kgrpnt ! index table containing the address of each (active) grid point
                                                         ! =1: not active grid point
                                                         ! >1: active grid point
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
    integer                    :: kf       ! loop counter
    integer                    :: kflux    ! k-index of flux
    integer, dimension(2*km+1) :: kicol    ! k-index of constituent point left of considered flux
    integer, dimension(2*km+1) :: kicor    ! k-index of constituent point right of considered flux
    integer                    :: kl       ! k-index left of considered point
    integer                    :: kr       ! k-index right of considered point
    integer                    :: m        ! loop counter
    integer                    :: mend     ! end index of loop over wl-points
    integer                    :: mu       ! index of point m+1
    integer                    :: nm       ! pointer to m
    integer                    :: nmu      ! pointer to m+1
    !
    real                       :: cl       ! interpolated constituent left of considered point
    real                       :: cr       ! interpolated constituent right of considered point
    real                       :: dif2d    ! effective horizontal eddy diffusivity coefficient
    real                       :: flux     ! total diffusion flux in x-direction
    real                       :: grad     ! strictly horizontal gradient of constituent
    real                       :: grad1    ! first strictly horizontal gradient of constituent
    real                       :: grad2    ! second strictly horizontal gradient of constituent
    real                       :: grmax    ! maximum of grad1 and grad2
    real                       :: grmin    ! minimum of grad1 and grad2
    real                       :: grsum    ! sum of grad1 and grad2
    real, dimension(km)        :: pocol    ! z-coordinate in constituent point nm,k
    real, dimension(km)        :: pocor    ! z-coordinate in constituent point nmu,k
    real, dimension(2*km+1)    :: poflu    ! z-coordinate of gradient flux
    real, dimension(0:2*km+1)  :: point    ! merge point
    real                       :: rdx      ! reciprocal of mesh size
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
    if (ltrace) call strace (ient,'SwashHDiffZplane1DH')
    !
    rdx = 1./dx
    !
    stabmx = 0.5 * dx * dx / dt
    !
    ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    do m = mfu, mend
       !
       mu = m + 1
       !
       nm  = kgrpnt(m )
       nmu = kgrpnt(mu)
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
       if ( hs(nm) > epsdry ) then
          !
          do k = 1, kmax
             pocol(k) = 0.5 * ( zkso( nm,k-1) + zkso( nm,k) )
             pocor(k) = 0.5 * ( zkso(nmu,k-1) + zkso(nmu,k) )
          enddo
          !
          ! merge adjacent layer-interfaces
          !
          kl = 0
          kr = 0
          do k = 0, 2*kmax+1
             if ( zkso(nm,kl) > zkso(nmu,kr) ) then
                point(k) = zkso(nm,kl)
                kl       = kl + 1
                if ( kl > kmax ) then
                   kflux = k + 1
                   point(kflux) = zkso(nmu,kr)
                   goto 10
                endif
             else
                point(k) = zkso(nmu,kr)
                kr       = kr + 1
                if ( kr > kmax ) then
                   kflux = k + 1
                   point(kflux) = zkso(nm,kl)
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
          kl = 1
          kr = 1
          do kf = 1, kflux
             kicol(kf) = 0
             kicor(kf) = 0
             do k = kl, kmax
                if ( (.not. poflu(kf) > zkso(nm,k-1)) .and. (.not. poflu(kf) < zkso(nm,k)) ) then
                   kicol(kf) = k
                   kl        = k
                   exit
                endif
             enddo
             do k = kr, kmax
                if ( (.not. poflu(kf) > zkso(nmu,k-1)) .and. (.not. poflu(kf) < zkso(nmu,k)) ) then
                   kicor(kf) = k
                   kr        = k
                   exit
                endif
             enddo
          enddo
          !
          ! compute diffusive flux using limiter
          !
          do kf = 1, kflux
             !
             kl = kicol(kf)
             kr = kicor(kf)
             !
             if ( kl /= 0 .and. kr /= 0 ) then
                !
                ! interpolation left point
                !
                k2 = kl
                if ( pocor(kr) > pocol(kl) ) then
                   k1 = k2 - 1
                   if ( k1 < 1 ) then
                      cl = rpo(nm,k2)
                   else
                      cl = ((pocor(kr)-pocol(k2))/(pocol(k1)-pocol(k2))) * rpo(nm,k1) + ((pocor(kr)-pocol(k1))/(pocol(k2)-pocol(k1))) * rpo(nm,k2)
                   endif
                else
                   k1 = k2 + 1
                   if ( k1 > kmax ) then
                      cl = rpo(nm,k2)
                   else
                      cl = ((pocor(kr)-pocol(k2))/(pocol(k1)-pocol(k2))) * rpo(nm,k1) + ((pocor(kr)-pocol(k1))/(pocol(k2)-pocol(k1))) * rpo(nm,k2)
                   endif
                endif
                !
                ! interpolation right point
                !
                k2 = kr
                if ( pocol(kl) > pocor(kr) ) then
                   k1 = k2 - 1
                   if ( k1 < 1 ) then
                      cr = rpo(nmu,k2)
                   else
                      cr = ((pocol(kl)-pocor(k2))/(pocor(k1)-pocor(k2))) * rpo(nmu,k1) + ((pocol(kl)-pocor(k1))/(pocor(k2)-pocor(k1))) * rpo(nmu,k2)
                   endif
                else
                   k1 = k2 + 1
                   if ( k1 > kmax ) then
                      cr = rpo(nmu,k2)
                   else
                      cr = ((pocol(kl)-pocor(k2))/(pocor(k1)-pocor(k2))) * rpo(nmu,k1) + ((pocol(kl)-pocor(k1))/(pocor(k2)-pocor(k1))) * rpo(nmu,k2)
                   endif
                endif
                !
                grad1 = rpo(nmu,kr) - cl
                grad2 = cr - rpo(nm,kl)
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
                flux = dif2d * rdx * rdx * ( point(kf-1) - point(kf) ) * grad
                if ( m == ml ) flux = 0.
                !
                rhsc( nm,kl) = rhsc( nm,kl) + flux
                rhsc(nmu,kr) = rhsc(nmu,kr) - flux
                !
             endif
             !
          enddo
          !
       endif
       !
    enddo
    !
end subroutine SwashHDiffZplane1DH
