subroutine SwashComputFlow
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
!   Computes water level and flow velocities by means of solving the shallow water equations
!
!   Method
!
!   Time integration is based on a semi-implicit approach that is unconditionally stable
!   with respect to gravity waves, bottom shear stresses and vertical eddy viscosity terms
!   Alternatively, water level gradients and the depth-integrated continuity equation are
!   treated explicitly using the leap-frog technique
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr
    use outp_data, only: lwavoutp, lcuroutp
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
    if (ltrace) call strace (ient,'SwashComputFlow')
    !
    ! calculate water level and flow velocities
    !
    if ( oned ) then
       !
       ! 1D mode
       !
       if ( kmax == 1 ) then
          !
          ! solve (non-hydrostatic) depth-averaged shallow water equations
          !
          if ( mtimei == 1 ) then
             !
             ! explicit (leap-frog) approach
             !
             call SwashExpDep1DHflow ( u1(1,1), u0(1,1), ua(1,1), up(1,1), qx(1,1), qm(1,1), q(1,1), dq(1,1), gmatu(1,1,1), rho(1,1), ibl(1), ibr(1), kgrpnt(1,1) )
             if (STPNOW()) return
             !
          else if ( mtimei == 2 ) then
             !
             ! semi-implicit approach
             !
             if ( .not.mimetic ) then
                call SwashImpDep1DHflow ( u1(1,1), u0(1,1), ua(1,1), up(1,1), qx(1,1), qm(1,1), q(1,1), dq(1,1), gmatu(1,1,1), rho(1,1), ui(1,1), dqgrd(1,1), ibl(1), ibr(1), kgrpnt(1,1) )
             else
                call SwashImpDepM1DHflow ( u1(1,1), u0(1,1), qx(1,1), qm(1,1), q(1,1), dq(1,1), gmatu(1,1,1), rho(1,1), ui(1,1), dqgrd(1,1), ibl(1), ibr(1), kgrpnt(1,1) )
             endif
             if (STPNOW()) return
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
             if ( .not.lsubg ) then
                call SwashExpLay1DHflow  ( ibl(1), ibr(1), kgrpnt(1,1) )
             else
                call SwashExpLayP1DHflow ( ibl(1), ibr(1), kgrpnt(1,1) )
             endif
             if (STPNOW()) return
             !
          else if ( mtimei == 2 ) then
             !
             ! semi-implicit approach
             !
             if ( mimetic ) then
                call SwashImpLayM1DHflow ( ibl(1), ibr(1), kgrpnt(1,1) )
             else if ( lsubg ) then
                call SwashImpLayP1DHflow ( ibl(1), ibr(1), kgrpnt(1,1) )
             else
                call SwashImpLay1DHflow ( ibl(1), ibr(1), kgrpnt(1,1) )
             endif
             if (STPNOW()) return
             !
          endif
          !
       endif
       !
    else
       !
       ! 2D mode
       !
       if ( kmax == 1 ) then
          !
          ! solve (non-hydrostatic) depth-averaged shallow water equations
          !
          if ( mtimei == 1 ) then
             !
             ! explicit (leap-frog) approach
             !
             call SwashExpDep2DHflow ( u1(1,1), u0(1,1), ua(1,1), ua2(1,1), up(1,1), qx(1,1), qxm(1,1), v1(1,1), v0(1,1), qy(1,1), qym(1,1), q(1,1), dq(1,1), gmatu(1,1,1), gmatv(1,1,1), rho(1,1) )
             if (STPNOW()) return
             !
          else if ( mtimei == 2 ) then
             !
             ! semi-implicit approach
             !
             if ( .not.mimetic ) then
                call SwashImpDep2DHflow ( u1(1,1), u0(1,1), ua(1,1), ua2(1,1), up(1,1), qx(1,1), qxm(1,1), v1(1,1), v0(1,1), qy(1,1), qym(1,1), q(1,1), dq(1,1), gmatu(1,1,1), gmatv(1,1,1), rho(1,1), ui(1,1), dqgrdu(1,1), vi(1,1), dqgrdv(1,1) )
             else
                call SwashImpDepM2DHflow ( u1(1,1), u0(1,1), qx(1,1), qxm(1,1), v1(1,1), v0(1,1), qy(1,1), qym(1,1), q(1,1), dq(1,1), gmatu(1,1,1), gmatv(1,1,1), rho(1,1), ui(1,1), dqgrdu(1,1), vi(1,1), dqgrdv(1,1) )
             endif
             if (STPNOW()) return
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
             if ( .not.lsubg ) then
                call SwashExpLay2DHflow
             else
                call SwashExpLayP2DHflow
             endif
             if (STPNOW()) return
             !
          else if ( mtimei == 2 ) then
             !
             ! semi-implicit approach
             !
             if ( .not.lsubg ) then
                call SwashImpLay2DHflow
             else
                call SwashImpLayP2DHflow
             endif
             if (STPNOW()) return
             !
          endif
          !
       endif
       !
    endif
    !
    ! update/compute setup, wave height and mean current for output purposes
    !
    if ( lwavoutp .or. lcuroutp ) call SwashAverOutp ( 1 )
    !
end subroutine SwashComputFlow
