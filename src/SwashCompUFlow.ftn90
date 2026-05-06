subroutine SwashCompUFlow
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
!   Computes water level and flow velocities by means of solving the shallow water equations on triangular mesh
!
!   Method
!
!   Time integration is based on a semi-implicit approach that is unconditionally stable
!   with respect to gravity waves
!   Alternatively, water level gradients and the depth-integrated continuity equation are
!   treated explicitly using the leap-frog technique
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use outp_data, only: lwavoutp, lcuroutp
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashCompUFlow')
    !
    if ( kmax == 1 ) then
       !
       ! solve (non-hydrostatic) depth-integrated shallow water equations
       !
       if ( mtimei == 1 ) then
          !
          ! explicit (leap-frog) approach
          !
          call SwashExpDepUflow ( u1(1,1), u0(1,1), uvc(1,1,1), qn(1,1), quf(1,1,1), q(1,1), dq(1,1), gmatu(1,1,1), rho(1,1) )
          !
       else if ( mtimei == 2 ) then
          !
          ! semi-implicit approach
          !
          call SwashImpDepUflow ( u1(1,1), u0(1,1), uvc(1,1,1), qn(1,1), quf(1,1,1), q(1,1), dq(1,1), gmatu(1,1,1), rho(1,1) )
          !
       endif
       !
    else
       !
       ! solve (non-hydrostatic) layer-averaged shallow water equations with a fix number of layers
       !
       if ( mtimei == 1 ) then
          !
          ! explicit (leap-frog) approach
          !
          call SwashExpLayUflow
          !
       else if ( mtimei == 2 ) then
          !
          ! semi-implicit approach
          !
          call SwashImpLayUflow
          !
       endif
       !
    endif
    !
    ! update/compute setup, wave height and mean current for output purposes
    !
    if ( lwavoutp .or. lcuroutp ) call SwashAverOutp ( 1 )
    !
end subroutine SwashCompUFlow
