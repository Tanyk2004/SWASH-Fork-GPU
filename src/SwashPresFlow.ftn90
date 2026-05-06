subroutine SwashPresFlow
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
!    1.00: Dirk Rijnsdorp
!
!   Updates
!
!    1.00, August 2014: New subroutine
!
!   Purpose
!
!   Updates mask arrays for pressurized flow in water level and velocity points
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: m        ! loop counter
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nm       ! pointer to m,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashPresFlow')
    !
    if ( oned ) then
       !
       ! loop over water level points
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m,1)
          !
          if ( s1(nm) < -flos(nm) ) then
             !
             presp(nm) = 0
             !
          else
             !
             presp(nm) = 1
             !
          endif
          !
       enddo
       !
       ! exchange mask values for water level with neighbouring subdomains
       !
       call SWEXCHGI ( presp, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! loop over u-points
       !
       do m = mf, ml
          !
          mu = m + 1
          !
          nm  = kgrpnt(m ,1)
          nmu = kgrpnt(mu,1)
          !
          if ( presp(nm) == 0 .and. presp(nmu) == 0 ) then
             !
             presu(nm) = 0
             !
          else
             !
             presu(nm) = 1
             !
          endif
          !
       enddo
       !
       ! exchange mask values for u-velocity with neighbouring subdomains
       !
       call SWEXCHGI ( presu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    else
       !
       ! loop over water level points
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             if ( s1(nm) < -flos(nm) ) then
                !
                presp(nm) = 0
                !
             else
                !
                presp(nm) = 1
                !
             endif
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       presp(1) = 0
       !
       ! exchange mask values for water level with neighbouring subdomains
       !
       call SWEXCHGI ( presp, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize mask values for water level at appropriate boundaries in case of repeating grid
       !
       call periodici ( presp, kgrpnt, 1, 1 )
       !
       ! loop over u-points
       !
       do n = nfu, nl
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( presp(nm) == 0 .and. presp(nmu) == 0 ) then
                !
                presu(nm) = 0
                !
             else
                !
                presu(nm) = 1
                !
             endif
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       presu(1) = 0
       !
       ! exchange mask values for u-velocity with neighbouring subdomains
       !
       call SWEXCHGI ( presu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize mask values for u-velocity at appropriate boundaries in case of repeating grid
       !
       call periodici ( presu, kgrpnt, 1, 1 )
       !
       ! loop over v-points
       !
       do m = mfu, ml
          do n = nf, nl
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( presp(nm) == 0 .and. presp(num) == 0 ) then
                !
                presv(nm) = 0
                !
             else
                !
                presv(nm) = 1
                !
             endif
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       presv(1) = 0
       !
       ! exchange mask values for v-velocity with neighbouring subdomains
       !
       call SWEXCHGI ( presv, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize mask values for v-velocity at appropriate boundaries in case of repeating grid
       !
       call periodici ( presv, kgrpnt, 1, 1 )
       !
    endif
    !
end subroutine SwashPresFlow
