subroutine SwashInputGrid ( igrid1, igrid2, psname )
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
!    1.10,    March 2012: extension porous media
!    1.21,     July 2013: extension transport constituents
!    1.22,   August 2013: extension vegetation
!
!   Purpose
!
!   Reads parameters of an input grid
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use outp_data
    use m_parall, only: mxcgl, mycgl
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    integer,      intent(in) :: igrid1 ! grid number for which a scalar or first component of a vector is read
    integer,      intent(in) :: igrid2 ! grid number for which second component of a vector is read
    character(8), intent(in) :: psname ! name of point set
!
!   Local variables
!
    integer                  :: i        ! loop counter
    integer, save            :: ient = 0 ! number of entries in this subroutine
    real                     :: xqlen    ! length of input grid in x-direction
    real                     :: yqlen    ! length of input grid in y-direction
    logical                  :: KEYWIS   ! indicates whether keyword in user manual is found or not
    !
    type(opsdat), pointer    :: opstmp   ! list containing parameters for grid
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    call strace (ient,'SwashInputGrid')
    !
    ! ==========================================================================
    !
    !            | BOTtom   |
    !            |          |
    !            | WLEVel   |
    !            |          |
    !            | CURrent  |
    !            |          |
    !            | VX       |
    !            |          |
    !            | VY       |
    !            |          |
    !            | MWL      |
    !            |          |
    !            | ACURrent |
    !            |          |
    !            | AVX      |
    !            |          |
    !            | AVY      |
    !            |          |
    !            | FRiction |
    !            |          |
    !            | WInd     |
    !            |          |
    !            | WX       |
    ! INPgrid  (<            >) &
    !            | WY       |
    !            |          |
    !            | PRessure |
    !            |          |
    !            | CORIolis |
    !            |          |
    !            | POROsity |
    !            |          |
    !            | PSIZe    |
    !            |          |
    !            | HSTRUCtu |
    !            |          |
    !            | SALinity |
    !            |          |
    !            | TEMPerat |
    !            |          |
    !            | SEDiment |
    !            |          |
    !            | NPLAnts  |
    !            |          |
    !            | DRAFt    |
    !            |          |
    !            | LABel    |
    !
    !
    !    | REGular [xpinp] [ypinp] [alpinp] [mxinp] [myinp] [dxinp] [dyinp] |
    !    |                                                                  |
    !   <  CURVilinear STAGgered                                             > &
    !    |                                                                  |
    !    | UNSTRUCtured                                                     |
    !
    !    (EXCeption  [excval])                                                 &
    !
    !                                        | -> SEC  |
    !    (NONSTATionary [tbeginp] [deltinp] <     MIN   >  [tendinp])          &
    !                                        |    HR   |
    !                                        |    DAY  |
    !
    !    (NONUNIForm  [kmax])
    !
    ! ==========================================================================
    !
    call INKEYW ('STA', 'REG')
    if ( KEYWIS('CURV') ) then
       !
       if ( oned ) then
          call msgerr (4, '1D simulation cannot be done with curvilinear grid')
          return
       endif
       !
       if ( optg == 1 ) then
          call msgerr (4, 'curvilinear input grid cannot be used with rectangular computational grid')
          return
       endif
       !
       if ( optg == 5 ) then
          call msgerr (4, 'curvilinear input grid cannot be used with unstructured grid')
          return
       endif
       !
       igtype(igrid1) = 2
       !
       call INKEYW ('STA', ' ')
       if ( KEYWIS ('STAG') ) then
          lstag(igrid1) = .true.
          if ( igrid2 > 0 ) lstag(igrid2) = .true.
       endif
       !
       if ( mxcgl == 0 .and. mycgl == 0 ) then
          mxg(igrid1) = mxc - 1
          myg(igrid1) = myc - 1
       else
          mxg(igrid1) = mxcgl - 1
          myg(igrid1) = mycgl - 1
       endif
       !
    else if ( KEYWIS('UNSTRUC') ) then
       !
       if ( oned ) then
          call msgerr (4, '1D simulation cannot be done with unstructured grid')
          return
       endif
       !
       if ( optg == 1 ) then
          call msgerr (4, 'unstructured input grid cannot be used with rectangular computational grid')
          return
       endif
       !
       if ( optg == 3 ) then
          call msgerr (4, 'unstructured input grid cannot be used with curvilinear computational grid')
          return
       endif
       !
       igtype(igrid1) = 3
       !
       mxg(igrid1)    = nverts
       myg(igrid1)    = 1
       if ( igrid2 > 0 ) then
          igtype(igrid2) = 3
          mxg(igrid2)    = nverts
          myg(igrid2)    = 1
       endif
       !
    else
       call IGNORE('REG')
       !
       igtype(igrid1) = 1
       !
       if ( kspher == 0 ) then
          call READXY ('XPINP', 'YPINP', xpg(igrid1), ypg(igrid1), 'UNC', 0., 0.)
          call INREAL ('ALPINP',alpg(igrid1),'UNC',0.)
       else
          call READXY ('XPINP', 'YPINP', xpg(igrid1), ypg(igrid1), 'REQ', 0., 0.)
          call INREAL ('ALPINP',alpg(igrid1),'STA',0.)
       endif
       !
       ! alpg(igrid1) is made to be between -pi and pi
       !
       alpg(igrid1)  = pi2 * (alpg(igrid1)/360. - nint(alpg(igrid1)/360.))
       !
       cospg(igrid1) = cos(alpg(igrid1))
       sinpg(igrid1) = sin(alpg(igrid1))
       !
       ! user gives number of cells for each direction
       ! however, this program uses the number of grid points for each direction
       !
       mxg(igrid1) = mxg(igrid1) - 1
       myg(igrid1) = myg(igrid1) - 1
       call ININTG ('MXINP', mxg(igrid1),'RQI',-1)
       if ( oned ) then
          call ININTG ('MYINP', myg(igrid1),'STA',0)
          if ( myg(igrid1) /= 0 ) then
             call msgerr (1, '1D simulation: [myinp] set to zero !')
             myg(igrid1) = 0
          endif
       else
          call ININTG ('MYINP', myg(igrid1),'RQI',-1)
       endif
       mxg(igrid1) = mxg(igrid1) + 1
       myg(igrid1) = myg(igrid1) + 1
       !
       call INREAL ('DXINP',dxg(igrid1),'RQI',0.)
       call INREAL ('DYINP',dyg(igrid1),'STA',dxg(igrid1))
       !
    endif
    !
    call INKEYW ('STA', ' ')
    if ( KEYWIS ('EXC') ) then
       if ( oned ) call msgerr (3, 'keyword EXC not allowed in 1D mode')
       call INREAL ('EXCVAL', excfld(igrid1), 'REQ', 0.)
       if ( igrid2 > 0 ) excfld(igrid2) = excfld(igrid1)
    endif
    !
    leds(igrid1) = 1
    if ( igrid2 > 0 ) leds(igrid2) = 1
    !
    call INKEYW ('STA', ' ')
    if ( KEYWIS ('NONSTAT') ) then
       if ( nstatm == 0 ) call msgerr (3, 'keyword NONSTAT not allowed in stationary mode')
       nstatm = 1
       if ( igrid1 == 1  .or. igrid1 == 8  .or.                   &
            igrid1 == 10 .or. igrid1 == 11 .or. igrid1 == 12 .or. &
            igrid1 == 14 .or.                                     &
            igrid1 == 15 .or. igrid1 == 16 .or. igrid1 == 17 .or. &
            igrid1 == 19 .or.                                     &
            igrid1 == 20 .or. igrid1 == 21 .or. igrid1 == 22 .or. &
            igrid1 == 23 ) call msgerr (2, 'nonstationary input field not allowed in this case')
       call INCTIM (ITMOPT,'TBEGINP',iflbeg(igrid1),'REQ',0d0)
       call INITVD ('DELTINP', iflint(igrid1), 'REQ', 0d0)
       call INCTIM (ITMOPT,'TENDINP',iflend(igrid1),'STA',1.d20)
       ifldyn(igrid1) = 1
       ifltim(igrid1) = iflbeg(igrid1)
       if ( igrid2 > 0 ) then
          iflbeg(igrid2) = iflbeg(igrid1)
          iflint(igrid2) = iflint(igrid1)
          iflend(igrid2) = iflend(igrid1)
          ifldyn(igrid2) = ifldyn(igrid1)
          ifltim(igrid2) = ifltim(igrid1)
       endif
       if ( iflend(igrid1) < 0.9e20 ) then
          if ( mod(iflend(igrid1)-iflbeg(igrid1),iflint(igrid1)) > 0.01*iflint(igrid1) .and. &
               mod(iflend(igrid1)-iflbeg(igrid1),iflint(igrid1)) < 0.99*iflint(igrid1) )     &
             call msgerr (1, '[deltinp] is not a fraction of the period')
       endif
    endif
    !
    call INKEYW ('STA', ' ')
    if ( KEYWIS ('NONUNIF') .or. KEYWIS ('VERT') ) then
       if ( kmax == 1 ) call msgerr (3, 'keyword NONUNIF not allowed in depth-averaged mode')
       if ( igrid1 /= 15 .and. igrid1 /= 16 .and. igrid1 /= 17 .and. igrid1 /= 20 .and. igrid1 /= 21 ) then
          call msgerr (2, 'this input parameter cannot be nonuniform in vertical')
       endif
       call ININTG ('KMAX', ifllay(igrid1), 'RQI', 0)
       if ( ifllay(igrid1) /= kmax .and. ifllay(igrid1) /= 1 ) then
          if ( igrid1 /= 20 .and. igrid1 /= 21 ) then
             call msgerr (1, 'incorrect number of layers, input field now regarded as uniform')
             ifllay(igrid1) = 1
          else
             if (ifllay(igrid1) > kmax ) then
                call msgerr (2, 'more than maximum number of layers encountered with ambient current input')
                ifllay(igrid1) = kmax
             endif
          endif
       endif
       if ( igrid2 > 0 ) ifllay(igrid2) = ifllay(igrid1)
    endif
    !
    if ( mxg(igrid1) == 0 .or. myg(igrid1) == 0 ) return
    if ( igtype(igrid1) == 3 ) return
    !
    ! input grid is included in output data
    !
    allocate(opstmp)
    opstmp%psname = psname
    allocate(opstmp%xp(0))
    allocate(opstmp%yp(0))
    !
    if ( igtype(igrid1) == 1 ) then
       !
       opstmp%pstype = 'F'
       xqlen         = ( mxg(igrid1) - 1 ) * dxg(igrid1)
       opstmp%opr(3) = xqlen
       !
       if ( oned ) then
          yqlen = xqlen
       else
          yqlen = ( myg(igrid1) - 1 ) * dyg(igrid1)
       endif
       !
       opstmp%opr(4) = yqlen
       opstmp%opr(1) = xpg(igrid1)
       opstmp%opr(2) = ypg(igrid1)
       opstmp%opr(5) = alpg(igrid1)
       opstmp%opr(6) = 0.
       opstmp%opi(1) = mxg(igrid1)
       opstmp%opi(2) = myg(igrid1)
    else
       opstmp%pstype = 'H'
       opstmp%opr(1) = float(mxg(igrid1)-1)
       opstmp%opr(2) = float(myg(igrid1)-1)
       opstmp%opr(3) = 0.
       opstmp%opr(4) = 0.
       opstmp%opr(5) = 0.
       opstmp%opr(6) = 0.
       opstmp%opi(1) = mxg(igrid1)
       opstmp%opi(2) = myg(igrid1)
       !
    endif
    !
    if ( ITEST >= 50 .or. intes >= 5 ) then
       !
       if ( opstmp%pstype == 'F' ) then
          write (PRINTF,201) igrid1,'F',xqlen,yqlen,xpg(igrid1),ypg(igrid1),alpg(igrid1),mxg(igrid1),myg(igrid1)
       else
          write (PRINTF,202) igrid1,opstmp%pstype,mxg(igrid1)-1,myg(igrid1)-1,0,0,0.,mxg(igrid1),myg(igrid1)
       endif
       !
    endif
    !
    nullify(opstmp%nextops)
    if ( .not.lops ) then
       fops = opstmp
       cops => fops
       lops = .true.
    else
       cops%nextops => opstmp
       cops => opstmp
    endif
    !
    ! check whether rectilinear input grid equals computational grid
    !
    if ( igtype(igrid1) == 1 ) then
       !
       if ( optg /= 1 ) then
          !
          lflgrd(igrid1) = .false.
          !
       else
          !
          ! computational grid is rectilinear
          !
          if ( xpg(igrid1) /= xpc .or. ypg(igrid1) /= ypc .or. alpg(igrid1) /= alpc .or. mxg(igrid1) /= mxc-1 .or. myg(igrid1) /= myc-1 ) lflgrd(igrid1) = .false.
          !
       endif
       !
    endif
    !
    if ( igrid2 > 0 ) then
       xpg   (igrid2) = xpg   (igrid1)
       ypg   (igrid2) = ypg   (igrid1)
       alpg  (igrid2) = alpg  (igrid1)
       cospg (igrid2) = cospg (igrid1)
       sinpg (igrid2) = sinpg (igrid1)
       dxg   (igrid2) = dxg   (igrid1)
       dyg   (igrid2) = dyg   (igrid1)
       mxg   (igrid2) = mxg   (igrid1)
       myg   (igrid2) = myg   (igrid1)
       igtype(igrid2) = igtype(igrid1)
       lflgrd(igrid2) = lflgrd(igrid1)
    endif
    !
 201 format (' input grid parameters: ',/, &
             'igrid  frame type  xlenfr   ylenfr   xpfr       ypfr       alpfr   mxfr myfr',/, &
             i3,9x,a,4x, 2(1x,e8.3), 2(1x,e10.3), f7.3, 2(1x,i4))
 202 format (' input grid parameters: ',/, &
             'igrid  frame type  xmaxfr ymaxfr xminfr yminfr  alpfr    mxfr myfr',/, &
             i3,9x,a,1x,4(3x,i4),5x,e8.3,2(1x,i4))
    !
end subroutine SwashInputGrid
