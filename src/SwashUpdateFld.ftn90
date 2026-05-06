subroutine SwashUpdateFld ( igr1, igr2, arrx, arry, cosfc, sinfc, arrfx, arrfy, xcgrid, ycgrid, kgrpnt )
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
!    1.00, February 2010: New subroutine
!
!   Purpose
!
!   Updates user-defined input fields, maps onto computational grid and interpolates in time
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use SwashTimecomm
    use m_parall
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                              :: igr1     ! grid number for which a scalar or first component of a vector is read
    integer, intent(in)                              :: igr2     ! grid number for which second component of a vector is read
    integer, dimension(mxc,myc), intent(in)          :: kgrpnt   ! index table containing the address of each (active) grid point
                                                                 ! =1; not active grid point
                                                                 ! >1; active grid point
    !
    real, dimension(mcgrd,3), intent(out)            :: arrfx    ! input field of scalar or x-component at three time levels mapped onto computational grid
    real, dimension(mcgrd,3), intent(out)            :: arrfy    ! input field of y-component at three time levels mapped onto computational grid
    real, dimension(mxg(igr1)*myg(igr1)), intent(in) :: arrx     ! user-defined input field of scalar or x-component
    real, dimension(*), intent(in)                   :: arry     ! user-defined input field of y-component
    real, intent(in)                                 :: cosfc    ! cosine of angle between input grid and computational grid
    real, intent(in)                                 :: sinfc    ! sine of angle between input grid and computational grid
    real, dimension(mxc,myc)  , intent(in)           :: xcgrid   ! coordinates of computational grid in x-direction
    real, dimension(mxc,myc)  , intent(in)           :: ycgrid   ! coordinates of computational grid in y-direction
!
!   Local variables
!
    integer                                          :: i        ! loop counter
    integer, save                                    :: ient = 0 ! number of entries in this subroutine
    integer                                          :: igl      ! counter in x-direction of global domain
    integer                                          :: indx     ! index point
    integer                                          :: ix       ! loop counter in x-direction of subdomain
    integer                                          :: iy       ! loop counter in y-direction of subdomain
    integer                                          :: j        ! loop counter
    integer                                          :: jgl      ! counter in y-direction of global domain
    !
    real*8                                           :: dval     ! double precision value
    real                                             :: fac      ! weighted length of vector
    real                                             :: fac1     ! length of vector at current time
    real                                             :: fac3     ! length of vector at time last read
    real                                             :: outval   ! a value taken outside input grid
    real                                             :: SVALQI   ! interpolated value of an input array for point given in problem coordinates
    real*8                                           :: timlr    ! time of one before last reading input field
    real                                             :: uu       ! local u-velocity component or scalar
    real                                             :: vtot     ! length of vector
    real                                             :: vv       ! local v-velocity component
    real                                             :: wf1      ! weighting factor for interpolation in time
    real                                             :: wf3      ! weighting factor for interpolation in time
    real                                             :: xp       ! x-coordinate of grid point
    real                                             :: yp       ! y-coordinate of grid point
    !
    logical                                          :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdateFld')
    !
    arrfx(:,1) = arrfx(:,2)
    if ( igr2 > 0 ) arrfy(:,1) = arrfy(:,2)
    !
    if ( igr1 == 2 .or. igr1 == 5 ) then
       outval = 0.
    else if ( igr1 == 4 .or. igr1 == 7 .or. igr1 == 13 ) then
       outval = NEAREST
    endif
    !
    timlr = timco - dt
    !
 100 if ( .not. timco > ifltim(igr1) ) goto 110
    !
    timlr = ifltim(igr1)
    ifltim(igr1) = ifltim(igr1) + iflint(igr1)
    if ( ifltim(igr1) > iflend(igr1) ) then
       ifltim(igr1) = 1.e10
       if ( igr2 > 0 ) ifltim(igr2) = ifltim(igr1)
       goto 110
    endif
    if ( iflnds(igr1) > 0 ) then
       if ( INODE == MASTER ) then
          call INAR2D( arrx, mxg(igr1), myg(igr1),               &
                       iflndf(igr1),                             &
                       iflnds(igr1), iflifm(igr1), iflfrm(igr1), &
                       iflidl(igr1), iflfac(igr1),               &
                       iflnhd(igr1), iflnhf(igr1) )
          if (STPNOW()) return
       endif
       call SWBROADC( iflidl(igr1), 1, SWINT )
       if ( iflidl(igr1) < 0 ) then
          ! end of file was encountered
          ifltim(igr1) = 1.e10
          if ( igr2 > 0 ) ifltim(igr2) = ifltim(igr1)
          goto 110
       else
          call SWBROADC( arrx, mxg(igr1)*myg(igr1), SWREAL )
       endif
    else
       if ( ITEST >= 20 ) then
          call msgerr (1, 'no read of input field because of undefined reference unit number')
          write (PRINTF, 201) igr1
       endif
    endif
    !
    if ( igr2 > 0 ) then
       ifltim(igr2) = ifltim(igr1)
       if ( iflnds(igr2) > 0 ) then
          if ( INODE == MASTER ) then
             call INAR2D( arry, mxg(igr2), myg(igr2), iflndf(igr2),   &
                          iflnds(igr2), iflifm(igr2), iflfrm(igr2),   &
                          iflidl(igr2), iflfac(igr2), iflnhd(igr2),   &
                          0 )
             if (STPNOW()) return
          endif
          call SWBROADC( arry, mxg(igr2)*myg(igr2), SWREAL )
       endif
    endif
    !
    ! mapping onto computational grid
    !
    if ( optg /= 5 ) then
       !
       ! structured grid
       !
       do j = 1, myc - 1
          do i = 1, mxc - 1
             indx = kgrpnt(i,j)
             if ( indx > 1 ) then
                igl = i + MXF - 1
                jgl = j + MYF - 1
                xp = xcgrid(i,j)
                yp = ycgrid(i,j)
                if ( lflgrd(igr1) ) then
                   uu = arrx((jgl-1)*mxg(igr1)+igl)
                else
                   uu = SVALQI (xp, yp, igr1, arrx, outval)
                endif
                if ( igr2 == 0 ) then
                   arrfx(indx,3) = uu
                else
                   if ( lflgrd(igr2) ) then
                      vv = arry((jgl-1)*mxg(igr2)+igl)
                   else
                      vv = SVALQI (xp, yp, igr2, arry, outval)
                   endif
                   arrfx(indx,3) =  uu*cosfc + vv*sinfc
                   arrfy(indx,3) = -uu*sinfc + vv*cosfc
                endif
             endif
          enddo
       enddo
       !
    else
       !
       ! unstructured grid
       !
       do i = 1, nverts
          xp = xcugrd(i)
          yp = ycugrd(i)
          if ( lflgrd(igr1) ) then
             uu = arrx(i)
          else
             uu = SVALQI (xp, yp, igr1, arrx, outval)
          endif
          if ( igr2 == 0 ) then
             arrfx(i,3) = uu
          else
             if ( lflgrd(igr2) ) then
                vv = arry(i)
             else
                vv = SVALQI (xp, yp, igr2, arry, outval)
             endif
             arrfx(i,3) =  uu*cosfc + vv*sinfc
             arrfy(i,3) = -uu*sinfc + vv*cosfc
          endif
       enddo
       !
    endif
    !
    goto 100
    !
    ! interpolation in time
    !
 110 dval = ( timco - timlr ) / ( ifltim(igr1) - timlr )
    wf3 = real(dval)
    wf1 = 1. - wf3
    !
    if ( ITEST >= 60 ) write (PRTEST, 202) igr1, timco, ifltim(igr1), wf1, wf3
    !
    do i = 2, mcgrd
       uu = wf1 * arrfx(i,1) + wf3 * arrfx(i,3)
       if ( igr2 == 0 ) then
          arrfx(i,2) = uu
       else
          vv = wf1 * arrfy(i,1) + wf3 * arrfy(i,3)
          vtot = sqrt (uu*uu + vv*vv)
          !
          ! procedure to prevent loss of magnitude due to interpolation
          !
          if ( vtot > 0. ) then
             !
             fac1 = sqrt (arrfx(i,1)**2 + arrfy(i,1)**2)
             fac3 = sqrt (arrfx(i,3)**2 + arrfy(i,3)**2)
             fac  = wf1 * fac1 + wf3 * fac3
             !
             arrfx(i,2) = fac * uu / vtot
             arrfy(i,2) = fac * vv / vtot
             !
          else
             !
             arrfx(i,2) = uu
             arrfy(i,2) = vv
             !
          endif
       endif
    enddo
    !
 201 format (' field nr.', i2)
 202 format (' input field', i2, ' interp at ', 2f9.0, 2f8.3)
    !
end subroutine SwashUpdateFld
