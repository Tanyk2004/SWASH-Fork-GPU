subroutine SwashCoorOutp ( pstype, opi, opr, x, y, mip, xp, yp, xc, yc, kgrpnt, xcgrid, ycgrid )
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
!   Calculates grid coordinates of the output points
!
!   Modules used
!
    use ocpcomm3
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata3
    use SwashCommdata4
    use m_parall
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc,myc), intent(in) :: kgrpnt   ! index table containing the address of each (active) grid point
                                                        ! =1; not active grid point
                                                        ! >1; active grid point
    integer, intent(in)                     :: mip      ! number of output points
    integer, dimension(2), intent(in)       :: opi      ! some integers
    !
    real, dimension(5), intent(in)          :: opr      ! some real coefficients
    real, dimension(mip), intent(in)        :: x        ! user coordinates of output point (curve, vertices, etc.) in x-direction
    real, dimension(mip), intent(out)       :: xc       ! computational coordinates (or broken coordinates) of output point in x-direction
    real, dimension(mxc,myc)  , intent(in)  :: xcgrid   ! coordinates of computational grid in x-direction
    real, dimension(mip), intent(out)       :: xp       ! user coordinates of output point in x-direction
    real, dimension(mip), intent(in)        :: y        ! user coordinates of output point (curve, vertices, etc.) in y-direction
    real, dimension(mip), intent(out)       :: yc       ! computational coordinates (or broken coordinates) of output point in y-direction
    real, dimension(mxc,myc)  , intent(in)  :: ycgrid   ! coordinates of computational grid in y-direction
    real, dimension(mip), intent(out)       :: yp       ! user coordinates of output point in y-direction
    !
    character(1), intent(in)                :: pstype   ! type of output point set
!
!   Local variables
!
    integer                                 :: i        ! loop counter
    integer, save                           :: ient = 0 ! number of entries in this subroutine
    integer                                 :: itmp1    ! auxiliary integer
    integer                                 :: itmp2    ! auxiliary integer
    integer                                 :: itmp3    ! auxiliary integer
    integer                                 :: itmp4    ! auxiliary integer
    integer                                 :: itmp5    ! auxiliary integer
    integer                                 :: j        ! loop counter
    integer                                 :: k        ! counter
    !
    real                                    :: dcxq     ! increment in x-direction
    real                                    :: dcyq     ! increment in y-direction
    real                                    :: rtmp1    ! auxiliary real
    real                                    :: rtmp2    ! auxiliary real
    real                                    :: rtmp3    ! auxiliary real
    real                                    :: rtmp4    ! auxiliary real
    real                                    :: xca      ! local computational x-coordinate
    real                                    :: xcmax    ! maximum user x-coordinate
    real                                    :: xcmin    ! minimum user x-coordinate
    real                                    :: xpa      ! local user x-coordinate
    real                                    :: xpmax    ! maximum x-coordinate of frame
    real                                    :: xpmin    ! minimum x-coordinate of frame
    real                                    :: xx       ! local x-coordinate
    real                                    :: yca      ! local computational y-coordinate
    real                                    :: ycmax    ! maximum user y-coordinate
    real                                    :: ycmin    ! minimum user y-coordinate
    real                                    :: ypa      ! local user y-coordinate
    real                                    :: ypmax    ! maximum y-coordinate of frame
    real                                    :: ypmin    ! minimum y-coordinate of frame
    real                                    :: yy       ! local y-coordinate
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashCoorOutp')
    !
    if ( pstype == 'F' ) then
       !
       mxq = opi(1)
       myq = opi(2)
       !
       alpq  = opr(5)
       cospq = cos(alpq)
       sinpq = sin(alpq)
       !
       xpq   = opr(1)
       ypq   = opr(2)
       xqlen = opr(3)
       yqlen = opr(4)
       !
       xqp = -xpq * cospq - ypq * sinpq
       yqp =  xpq * sinpq - ypq * cospq
       !
       if ( mxq > 1 ) then
          dxq = xqlen / float(mxq-1)
       else
          dxq = 0.01
       endif
       if ( myq > 1 ) then
          dyq = yqlen / float(myq-1)
       else
          dyq = 0.01
       endif
       !
       k = 0
       do j = 1, myq
           yy  = ( j - 1 ) * dyq
           xpa = xpq - yy * sinpq
           ypa = ypq + yy * cospq
           do i = 1, mxq
              xx = ( i - 1 ) * dxq
              k = k + 1
              xp(k) = xpa + xx * cospq
              yp(k) = ypa + xx * sinpq
           enddo
       enddo
       !
    else if ( pstype == 'H' ) then
       !
       mxq = opi(1)
       myq = opi(2)
       !
       alpq  = opr(5)
       cospq = cos(alpq)
       sinpq = sin(alpq)
       !
       xcmax = opr(1)
       ycmax = opr(2)
       xcmin = opr(3)
       ycmin = opr(4)
       !
       xpmin =  1.e09
       ypmin =  1.e09
       xpmax = -1.e09
       ypmax = -1.e09
       !
       xpq = xpmin
       ypq = ypmin
       xqp = 0.
       yqp = 0.
       !
       xqlen = xpmax - xpmin
       yqlen = ypmax - ypmin
       !
       if ( mxq > 1 ) then
          dxq = xqlen / float(mxq-1)
       else
          dxq = 0.
       endif
       if ( myq > 1 ) then
          dyq = yqlen / float(myq-1)
       else
          dyq = 0.
       endif
       !
       xc(1) = xcmin
       yc(1) = ycmin
       !
       if ( mxq > 1 ) then
          dcxq = ( xcmax - xcmin ) / float(mxq-1)
       else
          dcxq = 0.
       endif
       if ( myq > 1 ) then
          dcyq = ( ycmax - ycmin ) / float(myq-1)
       else
          dcyq = 0.
       endif
       !
       if ( ITEST >= 120 ) then
          write (PRINTF, 201) xqlen, yqlen, mxq, myq, dcxq, dcyq, xc(1), yc(1), dxq, dyq
          write (PRINTF, 202) xcmin, xcmax, ycmin, ycmax, xpmin, xpmax, ypmin, ypmax
       endif
       !
       yy = ycmin - dcyq
       k = 0
       do j = 1 ,myq
          !
          xx = xcmin - dcxq
          yy = yy    + dcyq
          !
          do i = 1 ,mxq
             !
             k = k + 1
             xx = xx + dcxq
             !
             xc(k) = xx - real(MXF) + 1.
             yc(k) = yy - real(MYF) + 1.
             !
             if ( xc(k) >= -0.01 .and. xc(k) <= real(mxc-2)+0.01 .and. yc(k) >= -0.01 .and. yc(k) <= real(myc-2)+0.01 ) then
                !
                if ( kgrpnt(nint(xc(k))+1,nint(yc(k))+1) > 1 ) then
                   !
                   call EVALF ( xc(k)+1., yc(k)+1., xpa, ypa, xcgrid, ycgrid )
                   xp(k) = xpa
                   yp(k) = ypa
                   !
                else
                   !
                   xp(k) = ovexcv(1)
                   yp(k) = ovexcv(2)
                   !
                endif
                !
             else
                !
                xp(k) = ovexcv(1)
                yp(k) = ovexcv(2)
                !
             endif
             !
             if ( ITEST >= 200 ) write (PRTEST, 204) xp(k), yp(k), xc(k), yc(k)
             !
          enddo
       enddo
       goto 100
       !
    else if ( pstype == 'C' .or. pstype == 'P' .or. pstype == 'N' .or. pstype == 'U' ) then
      !
      xp = x
      yp = y
      !
    endif
    !
    ! transform to computational grid, i.e. find broken coordinates
    !
    if ( ITEST >= 150 .and. optg == 1 ) write (PRTEST, 203) xcp, ycp, cospc, sinpc, dx, dy
    !
    ! the transformation to computational grid depends on the grid type: rectilinear or curvilinear
    ! (Note: unstructured grid does not employed broken coordinates)
    !
    if ( optg == 1 ) then
       !
       do i = 1, mip
          !
          xc(i) = ( xcp + xp(i)*cospc + yp(i)*sinpc ) / dx
          xc(i) = xc(i) - real(MXF) + 1.
          !
          if ( oned ) then
             yc(i) = 0
          else
             yc(i) = ( ycp - xp(i)*sinpc + yp(i)*cospc ) / dy
             yc(i) = yc(i) - real(MYF) + 1.
          endif
          !
          if ( ITEST >= 250 ) write (PRTEST, 204) xp(i), yp(i), xc(i), yc(i)
          !
       enddo
       !
    else if ( optg == 3 ) then
       !
       do i = 1, mip
          !
          xpa    = xp(i)
          ypa    = yp(i)
          itmp1  = mxc
          itmp2  = myc
          itmp3  = mcgrd
          itmp4  = MXF
          itmp5  = MYF
          rtmp1  = XCLMIN
          rtmp2  = XCLMAX
          rtmp3  = YCLMIN
          rtmp4  = YCLMAX
          mxc    = MXCGL
          myc    = MYCGL
          mcgrd  = MCGRDGL
          MXF    = 1
          MYF    = 1
          XCLMIN = xcgmin
          XCLMAX = xcgmax
          YCLMIN = ycgmin
          YCLMAX = ycgmax
          !
          call CVMESH ( xpa, ypa, xca, yca, KGRPGL, XGRDGL ,YGRDGL )
          !
          mxc    = itmp1
          myc    = itmp2
          mcgrd  = itmp3
          MXF    = itmp4
          MYF    = itmp5
          XCLMIN = rtmp1
          XCLMAX = rtmp2
          YCLMIN = rtmp3
          YCLMAX = rtmp4
          !
          xc(i) = xca - real(MXF) + 1.
          yc(i) = yca - real(MYF) + 1.
          !
          if ( ITEST >= 250 ) write (PRTEST, 204) xp(i), yp(i), xc(i), yc(i)
          !
       enddo
       !
    endif
    !
 100 continue
    !
 201 format (' FRAME DATA :',/,' XQLEN      ,YQLEN      ',',MXQ ,MYQ , DCXQ    ,DCYQ     ,XC(1)    ,YC(1)','    ,DXQ       ,DYQ',/,1x,2(1x,e9.3),2(1x,i4),2x,6(1x,e9.3))
 202 format (' XCMIN        ,XCMAX  ,YCMIN     ,YCMAX    ,XPMIN','   ,XPMAX     ,YPMIN    ,YPMAX   ',/,8(1x,e9.3),/)
 203 format (' grid coefficients ', 6(1x,e12.4))
 204 format (' PROBLEM  COORD:', 2(1x,f12.4),/,'         COMPUT   COORD:', 2(1x,f12.4))
    !
end subroutine SwashCoorOutp
