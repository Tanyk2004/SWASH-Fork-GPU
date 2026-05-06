subroutine SwashUBotFrict ( u )
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
!    1.00, February 2020: New subroutine
!
!   Purpose
!
!   Calculates bottom friction coefficient in case of flexible mesh
!
!   Method
!
!   Based on 6 roughness methods:
!
!   1) a dimensionless constant,
!   2) Chezy formulation,
!   3) Manning formulation,
!   4) Colebrook-White formulation,
!   5) Nikuradse roughness height (logarithmic velocity profile assumed), or
!   6) linear bottom friction
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use m_genarr
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    real, dimension(nfaces), intent(in) :: u ! depth-averaged flow velocity at faces
!
!   Parameter variables
!
    integer, parameter :: maxnit = 100    ! maximum number of iterations
!
    real   , parameter :: cfix   = 1.0129 ! minimum value for argument of log10 function in Colebrook-White
    real   , parameter :: eps    = 0.01   ! convergence criterion
    real   , parameter :: erough = 33.0   ! empirical constant for logarithmic log-law in case of rough beds
    real   , parameter :: esmoot = 9.0    ! empirical constant for logarithmic log-law in case of smooth beds
    real   , parameter :: ev     = 11.6   ! edge of viscous sublayer
!
!   Local variables
!
    integer               :: icell    ! loop counter over cells
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer               :: nit      ! number of iterations
    integer, dimension(3) :: v        ! vertices of present cell
!
    real                  :: cf       ! local friction coefficient
    real                  :: cz       ! Chezy value
    real                  :: r        ! Reynolds number
    real                  :: s        ! magnitude u/ustar
    real                  :: sold     ! ratio u/ustar at previous iteration
    real                  :: utot     ! velocity magnitude
    !
    type(celltype), dimension(:), pointer :: cell ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUBotFrict')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    if ( irough == 1 .or. irough == 11 ) then
       !
       ! dimensionless constant or linear bottom friction (dimension is m/s)
       !
       if ( varfr ) then
          !
          do icell = 1, ncells
             !
             v(1) = cell(icell)%atti(CELLV1)
             v(2) = cell(icell)%atti(CELLV2)
             v(3) = cell(icell)%atti(CELLV3)
             !
             cf = ( fricf(v(1),2) + fricf(v(2),2) + fricf(v(3),2) )/ 3.
             !
             cfricu(icell) = cf
             !
          enddo
          !
       else
          !
          do icell = 1, ncells
             !
             cfricu(icell) = pbot(1)
             !
          enddo
          !
       endif
       !
    else if ( irough == 2 ) then
       !
       ! Chezy formulation
       !
       if ( varfr ) then
          !
          do icell = 1, ncells
             !
             v(1) = cell(icell)%atti(CELLV1)
             v(2) = cell(icell)%atti(CELLV2)
             v(3) = cell(icell)%atti(CELLV3)
             !
             cf = ( fricf(v(1),2) + fricf(v(2),2) + fricf(v(3),2) )/ 3.
             !
             if ( cf /= 0. ) then
                !
                cfricu(icell) = grav / ( cf * cf )
                !
             endif
             !
          enddo
          !
       else
          !
          do icell = 1, ncells
             !
             if ( pbot(1) /= 0. ) then
                !
                cfricu(icell) = grav / ( pbot(1) * pbot(1) )
                !
             endif
             !
          enddo
          !
       endif
       !
    else if ( irough == 3 ) then
       !
       ! Manning formulation
       !
       if ( varfr ) then
          !
          do icell = 1, ncells
             !
             v(1) = cell(icell)%atti(CELLV1)
             v(2) = cell(icell)%atti(CELLV2)
             v(3) = cell(icell)%atti(CELLV3)
             !
             cf = ( fricf(v(1),2) + fricf(v(2),2) + fricf(v(3),2) )/ 3.
             !
             if ( hs(icell) > 0. ) then
                !
                cfricu(icell) = grav * cf * cf / hs(icell)**(1./3.)
                !
             endif
             !
          enddo
          !
       else
          !
          do icell = 1, ncells
             !
             if ( hs(icell) > 0. ) then
                !
                cfricu(icell) = grav * pbot(1) * pbot(1) / hs(icell)**(1./3.)
                !
             endif
             !
          enddo
          !
       endif
       !
    else if ( irough == 4 ) then
       !
       ! Nikuradse roughness height
       !
       call perot ( u, 1, 1 )
       !
       if ( varfr ) then
          !
          do icell = 1, ncells
             !
             v(1) = cell(icell)%atti(CELLV1)
             v(2) = cell(icell)%atti(CELLV2)
             v(3) = cell(icell)%atti(CELLV3)
             !
             cf = ( fricf(v(1),2) + fricf(v(2),2) + fricf(v(3),2) )/ 3.
             !
             utot = sqrt( uvc(icell,1,1)*uvc(icell,1,1) + uvc(icell,1,2)*uvc(icell,1,2) )
             !
             if ( hs(icell) > 0. ) then
                !
                if ( cf /= 0. ) then
                   !
                   cfricu(icell) = ( vonkar / log( erough*hs(icell)/(exp(1.)*cf) ) )**2.
                   !
                else
                   !
                   r = utot*hs(icell)/(exp(1.)*kinvis)
                   if ( r < 0.001 ) r = 0.001
                   !
                   if ( r > ev**2 ) then
                      !
                      nit = 0
                      !
                      ! initial value for s
                      !
                      s    = sqrt(r)
                      sold = 0.
                      !
                      ! Newton-Raphson iteration
                      !
                      do
                         if ( abs(sold-s) < (eps*s) ) exit
                         !
                         nit  = nit + 1
                         sold = s
                         s    = sold*(1.+log(esmoot*r/sold))/(1.+vonkar*sold)
                         !
                         if ( .not. nit < maxnit ) then
                            call msgerr (1, 'no convergence in bottom friction computation')
                            s = sqrt(r)
                            exit
                         endif
                         !
                      enddo
                      !
                   else
                      !
                      s = sqrt(r)
                      !
                   endif
                   !
                   cfricu(icell) = 1./(s*s)
                   !
                endif
                !
             endif
             !
          enddo
          !
       else
          !
          do icell = 1, ncells
             !
             utot = sqrt( uvc(icell,1,1)*uvc(icell,1,1) + uvc(icell,1,2)*uvc(icell,1,2) )
             !
             if ( hs(icell) > 0. ) then
                !
                if ( pbot(2) /= 0. ) then
                   !
                   cfricu(icell) = ( vonkar / log( erough*hs(icell)/(exp(1.)*pbot(2)) ) )**2.
                   !
                else
                   !
                   r = utot*hs(icell)/(exp(1.)*kinvis)
                   if ( r < 0.001 ) r = 0.001
                   !
                   if ( r > ev**2 ) then
                      !
                      nit = 0
                      !
                      ! initial value for s
                      !
                      s    = sqrt(r)
                      sold = 0.
                      !
                      ! Newton-Raphson iteration
                      !
                      do
                         if ( abs(sold-s) < (eps*s) ) exit
                         !
                         nit  = nit + 1
                         sold = s
                         s    = sold*(1.+log(esmoot*r/sold))/(1.+vonkar*sold)
                         !
                         if ( .not. nit < maxnit ) then
                            call msgerr (1, 'no convergence in bottom friction computation')
                            s = sqrt(r)
                            exit
                         endif
                         !
                      enddo
                      !
                   else
                      !
                      s = sqrt(r)
                      !
                   endif
                   !
                   cfricu(icell) = 1./(s*s)
                   !
                endif
                !
             endif
             !
          enddo
          !
       endif
       !
    else if ( irough == 5 ) then
       !
       ! Colebrook-White formulation
       !
       if ( varfr ) then
          !
          do icell = 1, ncells
             !
             v(1) = cell(icell)%atti(CELLV1)
             v(2) = cell(icell)%atti(CELLV2)
             v(3) = cell(icell)%atti(CELLV3)
             !
             cf = ( fricf(v(1),2) + fricf(v(2),2) + fricf(v(3),2) )/ 3.
             !
             if ( hs(icell) > 0. .and. cf /= 0. ) then
                !
                cz = 18. * log10 ( max(12.*hs(icell) / cf, cfix) )
                cfricu(icell) = grav / cz**2
                !
             endif
             !
          enddo
          !
       else
          !
          do icell = 1, ncells
             !
             if ( hs(icell) > 0. .and. pbot(1) /= 0. ) then
                !
                cz = 18. * log10 ( max(12.*hs(icell) / pbot(1), cfix) )
                cfricu(icell) = grav / cz**2
                !
             endif
             !
          enddo
          !
       endif
       !
    else
       !
       call msgerr ( 4, 'unknown roughness method for bottom friction' )
       return
       !
    endif
    !
end subroutine SwashUBotFrict
