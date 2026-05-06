subroutine SwashCompUTrans
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
!    1.00, January 2023: New subroutine
!
!   Purpose
!
!   Computes constituents by means of solving the transport equations on triangular mesh
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use outp_data, only: ltraoutp
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashCompUTrans')
    !
    ! calculate constituents
    !
    if ( kmax == 1 ) then
       !
       ! solve depth-averaged transport equations
       !
       call SwashExpDepUtrans ( rp(1,1,1), rpo(1,1), rp1(1,1), rp0(1,1), rpi(1,1,0), u1(1,1), qn(1,1), fluxt(1,1) )
       if (STPNOW()) return
       !
    else
       !
       ! solve layer-averaged transport equations with a fix number of layers
       !
       call SwashExpLayUtrans
       if (STPNOW()) return
       !
    endif
    !
    ! update/compute mean concentrations for output purposes
    !
    if ( ltraoutp ) call SwashAverOutp ( 2 )
    !
end subroutine SwashCompUTrans
