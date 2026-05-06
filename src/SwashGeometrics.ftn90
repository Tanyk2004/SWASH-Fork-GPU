subroutine SwashGeometrics
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
!   Computes geometric quantities based on given grid
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashCommdata4
    use m_genarr
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: indx     ! point index
    integer       :: indxb    ! index of point below from point of consideration
    integer       :: indxl    ! index of point left from point of consideration
    integer       :: indxr    ! index of point right from point of consideration
    integer       :: indxu    ! index of point up from point of consideration
    integer       :: ix       ! loop counter
    integer       :: iy       ! loop counter
    !
    real          :: coslat   ! cosine of latitude
    real          :: dx1      ! first component of covariant base vector a_(1)
    real          :: dx2      ! first component of covariant base vector a_(2)
    real          :: dy1      ! second component of covariant base vector a_(1)
    real          :: dy2      ! second component of covariant base vector a_(2)
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashGeometrics')
    !
    if ( optg == 1 ) then
       !
       dx1 = dx * cospc
       dy1 = dx * sinpc
       !
       dx2 = -dy * sinpc
       dy2 =  dy * cospc
       !
       if ( kspher > 0 ) then
          !
          dx1 = lendeg * dx1
          dy1 = lendeg * dy1
          dx2 = lendeg * dx2
          dy2 = lendeg * dy2
          !
       endif
       !
       ! compute metric tensors g_11, g_22 and Jacobian in u-points
       !
       do ix = 1, mxc-1
          do iy = 2, myc-1
             indx = kgrpnt(ix,iy)
             if ( indx > 1 ) then
                !
                if ( kspher > 0 ) then
                   !
                   ! in case of spherical coordinates compute cos of latitude
                   ! note: latitude is in degrees
                   !
                   coslat = cos(degrad*(0.5*(ycgrid(ix,iy)+ycgrid(ix,iy-1)) + yoffs))
                   !
                else
                   !
                   coslat = 1.
                   !
                endif
                !
                gvu  (indx) = sqrt( coslat*coslat*dx1*dx1 + dy1*dy1 )
                guu  (indx) = sqrt( coslat*coslat*dx2*dx2 + dy2*dy2 )
                gsqsu(indx) = coslat * ( dy2*dx1 - dy1*dx2 )
                !
             endif
          enddo
          !
          indx  = kgrpnt(ix,1)
          indxu = kgrpnt(ix,2)
          gvu  (indx) = gvu  (indxu)
          guu  (indx) = guu  (indxu)
          gsqsu(indx) = gsqsu(indxu)
          !
          indx  = kgrpnt(ix,myc  )
          indxb = kgrpnt(ix,myc-1)
          gvu  (indx) = gvu  (indxb)
          guu  (indx) = guu  (indxb)
          gsqsu(indx) = gsqsu(indxb)
       enddo
       !
       ! compute metric tensors g_11, g_22 and Jacobian in v-points
       !
       do iy = 1, myc-1
          do ix = 2, mxc-1
             indx = kgrpnt(ix,iy)
             if ( indx > 1 ) then
                !
                if ( kspher > 0 ) then
                   !
                   coslat = cos(degrad*(0.5*(ycgrid(ix,iy)+ycgrid(ix-1,iy)) + yoffs))
                   !
                else
                   !
                   coslat = 1.
                   !
                endif
                !
                gvv  (indx) = sqrt( coslat*coslat*dx1*dx1 + dy1*dy1 )
                guv  (indx) = sqrt( coslat*coslat*dx2*dx2 + dy2*dy2 )
                gsqsv(indx) = coslat * ( dy2*dx1 - dy1*dx2 )
                !
             endif
          enddo
          !
          indx  = kgrpnt(1,iy)
          indxr = kgrpnt(2,iy)
          gvv  (indx) = gvv  (indxr)
          guv  (indx) = guv  (indxr)
          gsqsv(indx) = gsqsv(indxr)
          !
          indx  = kgrpnt(mxc  ,iy)
          indxl = kgrpnt(mxc-1,iy)
          gvv  (indx) = gvv  (indxl)
          guv  (indx) = guv  (indxl)
          gsqsv(indx) = gsqsv(indxl)
       enddo
       !
    else
       !
       ! compute metric tensors g_11, g_22 and Jacobian in u-points
       !
       do iy = 2, myc-1
          indx = kgrpnt(1,iy)
          if ( indx > 1 ) then
             !
             dx1 = 0.5*( xcgrid(2,iy) + xcgrid(2,iy-1) - xcgrid(1,iy) - xcgrid(1,iy-1) )
             dy1 = 0.5*( ycgrid(2,iy) + ycgrid(2,iy-1) - ycgrid(1,iy) - ycgrid(1,iy-1) )
             !
             dx2 = xcgrid(1,iy) - xcgrid(1,iy-1)
             dy2 = ycgrid(1,iy) - ycgrid(1,iy-1)
             !
             if ( kspher > 0 ) then
                !
                coslat = cos(degrad*(0.5*(ycgrid(1,iy)+ycgrid(1,iy-1)) + yoffs))
                !
                dx1 = lendeg * coslat * dx1
                dy1 = lendeg * dy1
                dx2 = lendeg * coslat * dx2
                dy2 = lendeg * dy2
                !
             endif
             !
             gvu  (indx) = sqrt( dx1*dx1 + dy1*dy1 )
             guu  (indx) = sqrt( dx2*dx2 + dy2*dy2 )
             gsqsu(indx) = dy2*dx1 - dy1*dx2
             !
          endif
          do ix = 2, mxc-1
             indx = kgrpnt(ix,iy)
             if ( indx > 1 ) then
                !
                dx1 = 0.25*( xcgrid(ix+1,iy) + xcgrid(ix+1,iy-1) - xcgrid(ix-1,iy) - xcgrid(ix-1,iy-1) )
                dy1 = 0.25*( ycgrid(ix+1,iy) + ycgrid(ix+1,iy-1) - ycgrid(ix-1,iy) - ycgrid(ix-1,iy-1) )
                !
                dx2 = xcgrid(ix,iy) - xcgrid(ix,iy-1)
                dy2 = ycgrid(ix,iy) - ycgrid(ix,iy-1)
                !
                if ( kspher > 0 ) then
                   !
                   coslat = cos(degrad*(0.5*(ycgrid(ix,iy)+ycgrid(ix,iy-1)) + yoffs))
                   !
                   dx1 = lendeg * coslat * dx1
                   dy1 = lendeg * dy1
                   dx2 = lendeg * coslat * dx2
                   dy2 = lendeg * dy2
                   !
                endif
                !
                gvu  (indx) = sqrt( dx1*dx1 + dy1*dy1 )
                guu  (indx) = sqrt( dx2*dx2 + dy2*dy2 )
                gsqsu(indx) = dy2*dx1 - dy1*dx2
                !
             endif
          enddo
       enddo
       !
       do ix = 1, mxc-1
          indx  = kgrpnt(ix,1)
          indxu = kgrpnt(ix,2)
          gvu  (indx) = gvu  (indxu)
          guu  (indx) = guu  (indxu)
          gsqsu(indx) = gsqsu(indxu)
          !
          indx  = kgrpnt(ix,myc  )
          indxb = kgrpnt(ix,myc-1)
          gvu  (indx) = gvu  (indxb)
          guu  (indx) = guu  (indxb)
          gsqsu(indx) = gsqsu(indxb)
       enddo
       !
       ! compute metric tensors g_11, g_22 and Jacobian in v-points
       !
       do ix = 2, mxc-1
          indx = kgrpnt(ix,1)
          if ( indx > 1 ) then
             !
             dx1 = xcgrid(ix,1) - xcgrid(ix-1,1)
             dy1 = ycgrid(ix,1) - ycgrid(ix-1,1)
             !
             dx2 = 0.5*( xcgrid(ix,2) + xcgrid(ix-1,2) - xcgrid(ix,1) - xcgrid(ix-1,1) )
             dy2 = 0.5*( ycgrid(ix,2) + ycgrid(ix-1,2) - ycgrid(ix,1) - ycgrid(ix-1,1) )
             !
             if ( kspher > 0 ) then
                !
                coslat = cos(degrad*(0.5*(ycgrid(ix,1)+ycgrid(ix-1,1)) + yoffs))
                !
                dx1 = lendeg * coslat * dx1
                dy1 = lendeg * dy1
                dx2 = lendeg * coslat * dx2
                dy2 = lendeg * dy2
                !
             endif
             !
             gvv  (indx) = sqrt( dx1*dx1 + dy1*dy1 )
             guv  (indx) = sqrt( dx2*dx2 + dy2*dy2 )
             gsqsv(indx) = dy2*dx1 - dy1*dx2
             !
          endif
          do iy = 2, myc-1
             indx = kgrpnt(ix,iy)
             if ( indx > 1 ) then
                !
                dx1 = xcgrid(ix,iy) - xcgrid(ix-1,iy)
                dy1 = ycgrid(ix,iy) - ycgrid(ix-1,iy)
                !
                dx2 = 0.25*( xcgrid(ix,iy+1) + xcgrid(ix-1,iy+1) - xcgrid(ix,iy-1) - xcgrid(ix-1,iy-1) )
                dy2 = 0.25*( ycgrid(ix,iy+1) + ycgrid(ix-1,iy+1) - ycgrid(ix,iy-1) - ycgrid(ix-1,iy-1) )
                !
                if ( kspher > 0 ) then
                   !
                   coslat = cos(degrad*(0.5*(ycgrid(ix,iy)+ycgrid(ix-1,iy)) + yoffs))
                   !
                   dx1 = lendeg * coslat * dx1
                   dy1 = lendeg * dy1
                   dx2 = lendeg * coslat * dx2
                   dy2 = lendeg * dy2
                   !
                endif
                !
                gvv  (indx) = sqrt( dx1*dx1 + dy1*dy1 )
                guv  (indx) = sqrt( dx2*dx2 + dy2*dy2 )
                gsqsv(indx) = dy2*dx1 - dy1*dx2
                !
             endif
          enddo
       enddo
       !
       do iy = 1, myc-1
          indx  = kgrpnt(1,iy)
          indxr = kgrpnt(2,iy)
          gvv  (indx) = gvv  (indxr)
          guv  (indx) = guv  (indxr)
          gsqsv(indx) = gsqsv(indxr)
          !
          indx  = kgrpnt(mxc  ,iy)
          indxl = kgrpnt(mxc-1,iy)
          gvv  (indx) = gvv  (indxl)
          guv  (indx) = guv  (indxl)
          gsqsv(indx) = gsqsv(indxl)
       enddo
       !
    endif
    !
    ! compute Jacobian in wl-points
    !
    do iy = 2, myc-1
       do ix = 2, mxc-1
          indx  = kgrpnt(ix  ,iy  )
          indxb = kgrpnt(ix  ,iy-1)
          indxl = kgrpnt(ix-1,iy  )
          if ( indx > 1 ) then
             !
             dx1 = 0.5 * ( gvv(indx) + gvv(indxb) )
             dy1 = 0.5 * ( guu(indx) + guu(indxl) )
             !
             gsqs(indx) = dx1 * dy1
             !
          endif
       enddo
    enddo
    !
end subroutine SwashGeometrics
