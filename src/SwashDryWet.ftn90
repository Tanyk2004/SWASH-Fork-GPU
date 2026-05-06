subroutine SwashDryWet
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
!    1.10, March 2012: extension porous media
!
!   Purpose
!
!   Performs wetting and drying checks in water level and velocity points
!
!   Method
!
!   Update mask arrays
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
    integer       :: md       ! index of point m-1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashDryWet')
    !
    if ( oned ) then
       !
       ! loop over u-points
       !
       do m = mf, ml
          !
          nm = kgrpnt(m,1)
          !
          if ( hu(nm) > epshu ) then
             !
             wetu(nm) = 1
             !
          else
             !
             wetu(nm) = 0
             !
          endif
          !
       enddo
       !
       ! exchange mask values for u-velocity with neighbouring subdomains
       !
       call SWEXCHGI ( wetu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! loop over water level points
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          !
          if ( ( 1 - wetu(nm) ) * ( 1 - wetu(nmd) ) == 0 ) then
             !
             wets(nm) = 1
             !
          else
             !
             wets(nm) = 0
             !
          endif
          !
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
             if ( hu(nm) > epshu ) then
                !
                wetu(nm) = 1
                !
             else
                !
                wetu(nm) = 0
                !
             endif
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       wetu(1) = 0
       !
       ! exchange mask values for u-velocity with neighbouring subdomains
       !
       call SWEXCHGI ( wetu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize mask values for u-velocity at appropriate boundaries in case of repeating grid
       !
       call periodici ( wetu, kgrpnt, 1, 1 )
       !
       ! loop over v-points
       !
       do m = mfu, ml
          do n = nf, nl
             !
             nm = kgrpnt(m,n)
             !
             if ( hv(nm) > epshu ) then
                !
                wetv(nm) = 1
                !
             else
                !
                wetv(nm) = 0
                !
             endif
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       wetv(1) = 0
       !
       ! exchange mask values for v-velocity with neighbouring subdomains
       !
       call SWEXCHGI ( wetv, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize mask values for v-velocity at appropriate boundaries in case of repeating grid
       !
       call periodici ( wetv, kgrpnt, 1, 1 )
       !
       ! loop over water level points
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             if ( ( 1 - wetu(nm) ) * ( 1 - wetu(nmd) ) * ( 1 - wetv(nm) ) * ( 1 - wetv(ndm) ) == 0 ) then
                !
                wets(nm) = 1
                !
             else
                !
                wets(nm) = 0
                !
             endif
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       wets(1) = 0
       !
       ! synchronize mask values for water level at appropriate boundaries in case of repeating grid
       !
       call periodici ( wets, kgrpnt, 1, 1 )
       !
    endif
    !
    ! exchange mask values for water level with neighbouring subdomains
    !
    call SWEXCHGI ( wets, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
end subroutine SwashDryWet
