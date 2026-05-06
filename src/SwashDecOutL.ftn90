subroutine SwashDecOutL ( opi, opr, psname, pstype, mip, mxk, myk )
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
!   Decodes output point set
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata4
    use outp_data
!
    implicit none
!
!   Argument variables
!
    integer, intent(inout)            :: mip      ! number of output points
    integer, intent(out)              :: mxk      ! number of output points in x-direction
    integer, intent(out)              :: myk      ! number of output points in y-direction
    integer, dimension(2), intent(in) :: opi      ! some integers
    real, dimension(5), intent(in)    :: opr      ! some real coefficients
    character(8), intent(in)          :: psname   ! name of output point set
    character(1), intent(in)          :: pstype   ! type of output point set
!
!   Local variables
!
    integer, save                     :: ient = 0 ! number of entries in this subroutine
    !
    real                              :: xcmax    ! maximum user x-coordinate
    real                              :: xcmin    ! minimum user x-coordinate
    real                              :: xpmax    ! maximum x-coordinate of frame
    real                              :: xpmin    ! minimum x-coordinate of frame
    real                              :: ycmax    ! maximum user y-coordinate
    real                              :: ycmin    ! minimum user y-coordinate
    real                              :: ypmax    ! maximum y-coordinate of frame
    real                              :: ypmin    ! minimum y-coordinate of frame
    !
    logical                           :: EQREAL   ! compares two reals
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashDecOutL')
    !
    alpq  = 0.
    cospq = 1.
    sinpq = 0.
    !
    alcq  = 0.
    coscq = cos(alcq)
    sincq = sin(alcq)
    !
    if ( pstype == 'F' ) then
       !
       mxk = opi(1)
       myk = opi(2)
       mip = mxk * myk
       !
       xpq   = opr(1)
       ypq   = opr(2)
       alpq  = opr(5)
       cospq = cos(alpq)
       sinpq = sin(alpq)
       !
       xqp = -xpq * cospq - ypq * sinpq
       yqp =  xpq * sinpq - ypq * cospq
       !
       if ( EQREAL(outpar(2),1.) ) then
          ! directions will be w.r.t. frame coordinate system
          alcq = -alpq
       else
          ! directions will be w.r.t. user coordinate system (default)
          alcq = 0.
       endif
       !
       coscq = cos(alcq)
       sincq = sin(alcq)
       !
       if ( mxk > 1 ) then
          dxk = opr(3) / float(mxk-1)
       else
          dxk = 0.
       endif
       if ( myk > 1 ) then
          dyk = opr(4) / float(myk-1)
       else
          dyk = 0.
       endif
       !
    else if ( pstype == 'H' ) then
       !
       mxk  = opi(1)
       myk  = opi(2)
       mip  = mxk * myk
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
       !
       xqp = -xpq
       yqp = -ypq
       !
       alcq  = 0.
       coscq = cos(alcq)
       sincq = sin(alcq)
       !
       if ( mxk > 1 ) then
          dxk = ( xpmax - xpmin ) / float(mxk-1)
       else
          dxk = 0.
       endif
       if ( myk > 1 ) then
          dyk = ( ypmax - ypmin ) / float(myk-1)
       else
          dyk = 0.
       endif
       !
       if ( IOUTES >= 20 ) then
          write (PRINTF, 201) mxk, myk, xqp, yqp, dxk, dyk
          if ( pstype == 'H' ) write (PRINTF, 202) xcmin, xcmax, ycmin, ycmax, xpmax, xpmin
       endif
       !
    else if (pstype == 'C' .or. pstype == 'P' ) then
       !
       mxk = 0
       myk = 0
       !
    else if ( pstype == 'N' ) then
       !
       mxk = 0
       myk = 0
       !
    else if ( pstype == 'U' ) then
       !
       mxk = mip
       myk = 1
       !
    else
       !
       write (PRTEST, '(a)') ' error SwashDecOutL: no pstype defined'
       !
    endif
    !
    if ( psname == 'COMPGRID' ) then
       lcompgrd = .true.
    else
       lcompgrd = .false.
    endif
    !
    if ( ITEST >= 100 .or. IOUTES >= 30 ) then
       write (PRTEST, 203) psname, pstype, mip
       if ( pstype == 'F' .or. pstype == 'H' ) write (PRTEST, 204) mxk, myk, alpq, dxk, dyk
    endif
    !
 201 format (' SwashDecOutL:',/,'     MXK ,MYK , XQP     ,YQP       ,DXK','        ,DYK' ,/,2(1x,i5), 4(1x,e9.3))
 202 format (' XCMIN   ,XCMAX    ,YCMIN    ,YCMAX    ,', 'XPMAX    ,XPMIN :',/,6(1x,e9.3))
 203 format (' exit SwashDecOutL  ', a16, 2x, a1, 3i6)
 204 format (' SwashDecOutL: MXK, MYK,   ALPQ,   DXK,     DYK',/, 6x, 2i5, 4(1x,e9.3))
    !
end subroutine SwashDecOutL
