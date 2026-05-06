subroutine SwashCompUnstruc
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
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
!    1.00, January 2017: New subroutine
!
!   Purpose
!
!   Performs one full simulation step with unstructured grid
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashCompUnstruc')
    !
    ! compute flow
    !
!TIMG    call SWTSTA(51)
    if ( .not.momskip ) call SwashCompUFlow
!TIMG    call SWTSTO(51)
    if (STPNOW()) return
    !
    ! update water depths
    !
!TIMG    call SWTSTA(52)
    if ( kmax == 1 ) then
       call SwashUpdateUDepths ( u1 )
    else
       call SwashUpdateUDepths ( udep )
    endif
!TIMG    call SWTSTO(52)
    if (STPNOW()) return
    !
    ! update layer interfaces
    !
!TIMG    call SWTSTA(53)
    if ( kmax > 1 ) call SwashLayUIntfaces
!TIMG    call SWTSTO(53)
    !
    ! update mask arrays for wetting and drying
    !
!TIMG    call SWTSTA(54)
    call SwashUDryWet
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! update mask array for wave breaking
    !
!TIMG    call SWTSTA(54)
    if ( isurf /= 0 ) call SwashUBreakPoint
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! compute transport
    !
!TIMG    call SWTSTA(101)
    if ( itrans /= 0 ) call SwashCompUTrans
!TIMG    call SWTSTO(101)
    if (STPNOW()) return
    !
    ! compute 3D turbulence
    !
!TIMG    call SWTSTA(102)
    if ( iturb /= 0 ) call SwashCompUTurb
!TIMG    call SWTSTO(102)
    if (STPNOW()) return
    !
end subroutine SwashCompUnstruc
