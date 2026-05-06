subroutine SwashComputTurb
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
!    1.00, March 2011: New subroutine
!
!   Purpose
!
!   Computes turbulence quantities by means of solving a 3D turbulence model
!
!   Method
!
!   The well-known standard k-epsilon closure of Launder and Spalding (1974)
!   is employed for vertical mixing with which 3D turbulence is mainly generated
!   by vertical velocity gradients and bottom friction. Typically, its length
!   scale is smaller than the water depth.
!
!   B.E. Launder and D.B. Spalding, "The numerical computation of turbulent
!   flows", Comput. Meth. Appl. Mech. Engng., vol. 3, p. 269-289, 1974
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr
    use outp_data, only: lturoutp
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
    if (ltrace) call strace (ient,'SwashComputTurb')
    !
    ! calculate turbulence quantities
    !
    if ( oned ) then
       !
       ! 1D mode
       !
       if ( iturb /= 0 ) then
          !
          ! standard k-epsilon model
          !
          call SwashKepsMod1DH ( kgrpnt(1,1) )
          if (STPNOW()) return
          !
       endif
       !
    else
       !
       ! 2D mode
       !
       if ( iturb /= 0 ) then
          !
          ! standard k-epsilon model
          !
          call SwashKepsMod2DH
          if (STPNOW()) return
          !
       endif
       !
    endif
    !
    ! update/compute mean turbulence quantities for output purposes
    !
    if ( lturoutp ) call SwashAverOutp ( 3 )
    !
end subroutine SwashComputTurb
