subroutine SwashBotFrict ( u, v )
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
!    1.00, January 2010: New subroutine
!    1.05, January 2012: extension logarithmic wall-law
!
!   Purpose
!
!   Calculates bottom friction coefficient
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
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd), intent(in) :: u ! depth-averaged u-velocity
    real, dimension(mcgrd), intent(in) :: v ! depth-averaged v-velocity
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
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    integer       :: nit      ! number of iterations
    integer       :: nm       ! pointer to m,n
!
    real          :: cz       ! Chezy value
    real          :: r        ! Reynolds number
    real          :: s        ! magnitude u/ustar
    real          :: sold     ! ratio u/ustar at previous iteration
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashBotFrict')
    !
    if ( oned ) then
       !
       if ( irough == 1 .or. irough == 11 ) then
          !
          ! dimensionless constant or linear bottom friction (dimension is m/s)
          !
          if ( varfr ) then
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 ) then
                   !
                   cfricu(nm) = fricf(nm,2)
                   !
                endif
                !
             enddo
             !
          else
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 ) then
                   !
                   cfricu(nm) = pbot(1)
                   !
                endif
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
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 .and. fricf(nm,2) /= 0. ) then
                   !
                   cfricu(nm) = grav / ( fricf(nm,2) * fricf(nm,2) )
                   !
                endif
                !
             enddo
             !
          else
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 .and. pbot(1) /= 0. ) then
                   !
                   cfricu(nm) = grav / ( pbot(1) * pbot(1) )
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
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                   !
                   cfricu(nm) = grav * fricf(nm,2) * fricf(nm,2) / hum(nm)**(1./3.)
                   !
                endif
                !
             enddo
             !
          else
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                   !
                   cfricu(nm) = grav * pbot(1) * pbot(1) / hum(nm)**(1./3.)
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
          if ( varfr ) then
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                   !
                   if ( fricf(nm,2) /= 0. ) then
                      !
                      cfricu(nm) = ( vonkar / log( erough*hum(nm)/(exp(1.)*fricf(nm,2)) ) )**2.
                      !
                   else
                      !
                      r = abs(u(nm))*hum(nm)/(exp(1.)*kinvis)
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
                      cfricu(nm) = 1./(s*s)
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          else
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                   !
                   if ( pbot(2) /= 0. ) then
                      !
                      cfricu(nm) = ( vonkar / log( erough*hum(nm)/(exp(1.)*pbot(2)) ) )**2.
                      !
                   else
                      !
                      r = abs(u(nm))*hum(nm)/(exp(1.)*kinvis)
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
                      cfricu(nm) = 1./(s*s)
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
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 .and. hum(nm) > 0. .and. fricf(nm,2) /= 0. ) then
                   !
                   cz = 18. * log10 ( max(12.*hum(nm) / fricf(nm,2), cfix) )
                   cfricu(nm) = grav / cz**2
                   !
                endif
                !
             enddo
             !
          else
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( wetu(nm) == 1 .and. hum(nm) > 0. .and. pbot(1) /= 0. ) then
                   !
                   cz = 18. * log10 ( max(12.*hum(nm) / pbot(1), cfix) )
                   cfricu(nm) = grav / cz**2
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
    else
       !
       if ( irough == 1 .or. irough == 11 ) then
          !
          ! dimensionless constant or linear bottom friction (dimension is m/s)
          !
          if ( varfr ) then
             !
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 ) then
                      !
                      cfricu(nm) = fricf(nm,2)
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 ) then
                      !
                      cfricv(nm) = fricf(nm,2)
                      !
                   endif
                   !
                enddo
             enddo
             !
          else
             !
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 ) then
                      !
                      cfricu(nm) = pbot(1)
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 ) then
                      !
                      cfricv(nm) = pbot(1)
                      !
                   endif
                   !
                enddo
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
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 .and. fricf(nm,2) /= 0. ) then
                      !
                      cfricu(nm) = grav / ( fricf(nm,2) * fricf(nm,2) )
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 .and. fricf(nm,2) /= 0. ) then
                      !
                      cfricv(nm) = grav / ( fricf(nm,2) * fricf(nm,2) )
                      !
                   endif
                   !
                enddo
             enddo
             !
          else
             !
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 .and. pbot(1) /= 0. ) then
                      !
                      cfricu(nm) = grav / ( pbot(1) * pbot(1) )
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 .and. pbot(1) /= 0. ) then
                      !
                      cfricv(nm) = grav / ( pbot(1) * pbot(1) )
                      !
                   endif
                   !
                enddo
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
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                      !
                      cfricu(nm) = grav * fricf(nm,2) * fricf(nm,2) / hum(nm)**(1./3.)
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 .and. hvm(nm) > 0. ) then
                      !
                      cfricv(nm) = grav * fricf(nm,2) * fricf(nm,2) / hvm(nm)**(1./3.)
                      !
                   endif
                   !
                enddo
             enddo
             !
          else
             !
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                      !
                      cfricu(nm) = grav * pbot(1) * pbot(1) / hum(nm)**(1./3.)
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 .and. hvm(nm) > 0. ) then
                      !
                      cfricv(nm) = grav * pbot(1) * pbot(1) / hvm(nm)**(1./3.)
                      !
                   endif
                   !
                enddo
             enddo
             !
          endif
          !
       else if ( irough == 4 ) then
          !
          ! Nikuradse roughness height
          !
          if ( varfr ) then
             !
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                      !
                      if ( fricf(nm,2) /= 0. ) then
                         !
                         cfricu(nm) = ( vonkar / log( erough*hum(nm)/(exp(1.)*fricf(nm,2)) ) )**2.
                         !
                      else
                         !
                         r = abs(u(nm))*hum(nm)/(exp(1.)*kinvis)
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
                         cfricu(nm) = 1./(s*s)
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 .and. hvm(nm) > 0. ) then
                      !
                      if ( fricf(nm,2) /= 0. ) then
                         !
                         cfricv(nm) = ( vonkar / log( erough*hvm(nm)/(exp(1.)*fricf(nm,2)) ) )**2.
                         !
                      else
                         !
                         r = abs(v(nm))*hvm(nm)/(exp(1.)*kinvis)
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
                         cfricv(nm) = 1./(s*s)
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
          else
             !
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                      !
                      if ( pbot(2) /= 0. ) then
                         !
                         cfricu(nm) = ( vonkar / log( erough*hum(nm)/(exp(1.)*pbot(2)) ) )**2.
                         !
                      else
                         !
                         r = abs(u(nm))*hum(nm)/(exp(1.)*kinvis)
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
                         cfricu(nm) = 1./(s*s)
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 .and. hvm(nm) > 0. ) then
                      !
                      if ( pbot(2) /= 0. ) then
                         !
                         cfricv(nm) = ( vonkar / log( erough*hvm(nm)/(exp(1.)*pbot(2)) ) )**2.
                         !
                      else
                         !
                         r = abs(v(nm))*hvm(nm)/(exp(1.)*kinvis)
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
                         cfricv(nm) = 1./(s*s)
                         !
                      endif
                      !
                   endif
                   !
                enddo
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
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 .and. hum(nm) > 0. .and. fricf(nm,2) /= 0. ) then
                      !
                      cz = 18. * log10 ( max(12.*hum(nm) / fricf(nm,2), cfix) )
                      cfricu(nm) = grav / cz**2
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 .and. hvm(nm) > 0. .and. fricf(nm,2) /= 0. ) then
                      !
                      cz = 18. * log10 ( max(12.*hvm(nm) / fricf(nm,2), cfix) )
                      cfricv(nm) = grav / cz**2
                      !
                   endif
                   !
                enddo
             enddo
             !
          else
             !
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetu(nm) == 1 .and. hum(nm) > 0. .and. pbot(1) /= 0. ) then
                      !
                      cz = 18. * log10 ( max(12.*hum(nm) / pbot(1), cfix) )
                      cfricu(nm) = grav / cz**2
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( wetv(nm) == 1 .and. hvm(nm) > 0. .and. pbot(1) /= 0. ) then
                      !
                      cz = 18. * log10 ( max(12.*hvm(nm) / pbot(1), cfix) )
                      cfricv(nm) = grav / cz**2
                      !
                   endif
                   !
                enddo
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
       ! synchronize bottom friction coefficients at appropriate boundaries in case of repeating grid
       !
       call periodic ( cfricu, kgrpnt, 1, 1 )
       call periodic ( cfricv, kgrpnt, 1, 1 )
       !
    endif
    !
end subroutine SwashBotFrict
