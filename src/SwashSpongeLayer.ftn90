subroutine SwashSpongeLayer ( sponxl, sponxr, sponyb, sponyt )
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
!   Determines damping function for flow variables due to sponge layers
!
!   Method
!
!   The sponge layer formula as described in Mayer et al. (1998), Eqs (43) and (44), is employed
!
!   S. Mayer, A. Garapon and L.S. Sorensen
!   A fractional step method for unsteady free-surface flow with applications to nonlinear wave dynamics
!   IJNMF, vol. 28, 293-315, 1998
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashCommdata4
    use m_parall
!
    implicit none
!
!   Argument variables
!
    real, dimension(MCGRDGL), intent(inout) :: sponxl ! left sponge layer in x-direction in global domain
    real, dimension(MCGRDGL), intent(inout) :: sponxr ! right sponge layer in x-direction in global domain
    real, dimension(MCGRDGL), intent(inout) :: sponyb ! lower sponge layer in y-direction in global domain
    real, dimension(MCGRDGL), intent(inout) :: sponyt ! upper sponge layer in y-direction in global domain
!
!   Parameter variables
!
    real, parameter :: grt = 0.5  ! growth rate of relaxation
!
!   Local variables
!
    integer, save   :: ient = 0   ! number of entries in this subroutine
    integer         :: indx       ! point index
    integer         :: isl        ! index of starting point of sponge layer w.r.t. left boundary in 1D
    integer         :: isr        ! index of starting point of sponge layer w.r.t. right boundary in 1D
    integer         :: ix         ! index of point in x-direction
    integer         :: iy         ! index of point in y-direction
    !
    real            :: coslat     ! cosine of latitude
    real            :: dist       ! distance from starting point of sponge layer
    real            :: dxl        ! difference between two consecutive grid points in x-direction
    real            :: dyl        ! difference between two consecutive grid points in y-direction
    real            :: xl         ! normalised distance in x-direction
    real            :: yl         ! normalised distance in y-direction
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashSpongeLayer')
    !
    if ( spwidl > 0. .or. spwidr > 0. ) then
       !
       if ( oned ) then
          !
          if ( spwidl > 0. ) then
             !
             isl = 1 + nint(spwidl/dx)
             !
             ! loop over u-points in 1D (global) computational grid
             !
             do ix = 1, MXCGL-1
                !
                indx = KGRPGL(ix,1)
                !
                xl = 0.
                !
                if ( .not. ix > isl ) then
                   !
                   dist = ( isl - ix ) * dx
                   xl   = dist / spwidl
                   !
                   if ( xl > 0. ) then
                      !
                      sponxl(indx) = grt * xl**3 + (1.-grt) * xl**6
                      !
                   else
                      !
                      sponxl(indx) = 0.
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          endif
          !
          if ( spwidr > 0. ) then
             !
             isr = MXCGL - ( 1 + nint(spwidr/dx) )
             !
             ! loop over u-points in 1D (global) computational grid
             !
             do ix = 1, MXCGL-1
                !
                indx = KGRPGL(ix,1)
                !
                xl = 0.
                !
                if ( .not. ix < isr ) then
                   !
                   dist = ( ix - isr ) * dx
                   xl   = dist / spwidr
                   !
                   if ( xl > 0. ) then
                      !
                      sponxr(indx) = grt * xl**3 + (1.-grt) * xl**6
                      !
                   else
                      !
                      sponxr(indx) = 0.
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          endif
          !
       else
          !
          ! loop over u-points in 2D (global) computational grid
          !
          do iy = 2, MYCGL-1
             !
             if ( spwidl > 0. ) then
                !
                dist = 0.
                !
                do ix = 1, MXCGL-1
                   !
                   indx = KGRPGL(ix,iy)
                   !
                   if ( .not. dist > spwidl ) then
                      !
                      xl = 1. - dist / spwidl
                      !
                      if ( xl > 0. ) then
                         !
                         sponxl(indx) = grt * xl**3 + (1.-grt) * xl**6
                         !
                      else
                         !
                         sponxl(indx) = 0.
                         !
                      endif
                      !
                   endif
                   !
                   dxl = 0.5 * ( XGRDGL(ix+1,iy) + XGRDGL(ix+1,iy-1) - XGRDGL(ix,iy) - XGRDGL(ix,iy-1) )
                   dyl = 0.5 * ( YGRDGL(ix+1,iy) + YGRDGL(ix+1,iy-1) - YGRDGL(ix,iy) - YGRDGL(ix,iy-1) )
                   !
                   if ( kspher > 0 ) then
                      !
                      coslat = cos(degrad*(0.25*(YGRDGL(ix+1,iy)+YGRDGL(ix+1,iy-1)+YGRDGL(ix,iy)+YGRDGL(ix,iy-1)) + yoffs))
                      !
                      dxl = lendeg * coslat * dxl
                      dyl = lendeg * dyl
                      !
                   endif
                   !
                   dist = dist + sqrt( dxl*dxl + dyl*dyl )
                   !
                enddo
                !
             endif
             !
             if ( spwidr > 0. ) then
                !
                dist = 0.
                !
                do ix = MXCGL-1, 2, -1
                   !
                   indx = KGRPGL(ix,iy)
                   !
                   if ( .not. dist > spwidr ) then
                      !
                      xl = 1. - dist / spwidr
                      !
                      if ( xl > 0. ) then
                         !
                         sponxr(indx) = grt * xl**3 + (1.-grt) * xl**6
                         !
                      else
                         !
                         sponxr(indx) = 0.
                         !
                      endif
                      !
                   endif
                   !
                   dxl = 0.5 * ( XGRDGL(ix,iy) + XGRDGL(ix,iy-1) - XGRDGL(ix-1,iy) - XGRDGL(ix-1,iy-1) )
                   dyl = 0.5 * ( YGRDGL(ix,iy) + YGRDGL(ix,iy-1) - YGRDGL(ix-1,iy) - YGRDGL(ix-1,iy-1) )
                   !
                   if ( kspher > 0 ) then
                      !
                      coslat = cos(degrad*(0.25*(YGRDGL(ix,iy)+YGRDGL(ix,iy-1)+YGRDGL(ix-1,iy)+YGRDGL(ix-1,iy-1)) + yoffs))
                      !
                      dxl = lendeg * coslat * dxl
                      dyl = lendeg * dyl
                      !
                   endif
                   !
                   dist = dist + sqrt( dxl*dxl + dyl*dyl )
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
    if ( spwidb > 0. .or. spwidt > 0. ) then
       !
       ! loop over v-points in 2D (global) computational grid
       !
       do ix = 2, MXCGL-1
          !
          if ( spwidb > 0. ) then
             !
             dist = 0.
             !
             do iy = 1, MYCGL-1
                !
                indx = KGRPGL(ix,iy)
                !
                if ( .not. dist > spwidb ) then
                   !
                   yl = 1. - dist / spwidb
                   !
                   if ( yl > 0. ) then
                      !
                      sponyb(indx) = grt * yl**3 + (1.-grt) * yl**6
                      !
                   else
                      !
                      sponyb(indx) = 0.
                      !
                   endif
                   !
                endif
                !
                dxl = 0.5 * ( XGRDGL(ix,iy+1) + XGRDGL(ix-1,iy+1) - XGRDGL(ix,iy) - XGRDGL(ix-1,iy) )
                dyl = 0.5 * ( YGRDGL(ix,iy+1) + YGRDGL(ix-1,iy+1) - YGRDGL(ix,iy) - YGRDGL(ix-1,iy) )
                !
                if ( kspher > 0 ) then
                   !
                   coslat = cos(degrad*(0.25*(YGRDGL(ix,iy+1)+YGRDGL(ix,iy)+YGRDGL(ix-1,iy+1)+YGRDGL(ix-1,iy)) + yoffs))
                   !
                   dxl = lendeg * coslat * dxl
                   dyl = lendeg * dyl
                   !
                endif
                !
                dist = dist + sqrt( dxl*dxl + dyl*dyl )
                !
             enddo
             !
          endif
          !
          if ( spwidt > 0. ) then
             !
             dist = 0.
             !
             do iy = MYCGL-1, 2, -1
                !
                indx = KGRPGL(ix,iy)
                !
                if ( .not. dist > spwidt ) then
                   !
                   yl = 1. - dist / spwidt
                   !
                   if ( yl > 0. ) then
                      !
                      sponyt(indx) = grt * yl**3 + (1.-grt) * yl**6
                      !
                   else
                      !
                      sponyt(indx) = 0.
                      !
                   endif
                   !
                endif
                !
                dxl = 0.5 * ( XGRDGL(ix,iy) + XGRDGL(ix-1,iy) - XGRDGL(ix,iy-1) - XGRDGL(ix-1,iy-1) )
                dyl = 0.5 * ( YGRDGL(ix,iy) + YGRDGL(ix-1,iy) - YGRDGL(ix,iy-1) - YGRDGL(ix-1,iy-1) )
                !
                if ( kspher > 0 ) then
                   !
                   coslat = cos(degrad*(0.25*(YGRDGL(ix,iy)+YGRDGL(ix,iy-1)+YGRDGL(ix-1,iy)+YGRDGL(ix-1,iy-1)) + yoffs))
                   !
                   dxl = lendeg * coslat * dxl
                   dyl = lendeg * dyl
                   !
                endif
                !
                dist = dist + sqrt( dxl*dxl + dyl*dyl )
                !
             enddo
             !
          endif
          !
       enddo
       !
    endif
    !
end subroutine SwashSpongeLayer
