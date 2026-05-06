subroutine SwashUpdFlowFlds
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
!   Initializes / updates flow variables based on space varying input fields
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
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: ndmd     ! pointer to m-1,n-1
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
    if (ltrace) call strace (ient,'SwashUpdFlowFlds')
    !
    if ( initsf .or. ifldyn(7) == 1 ) then
       !
       ! determine water level in wl-points
       !
       if ( oned ) then
          !
          ! loop over wl-points in 1D computational grid
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,1)
             nmd = kgrpnt(md,1)
             !
             s1(nm) = 0.5 * ( wlevf(nm,2) + wlevf(nmd,2) )
             !
          enddo
          !
       else
          !
          ! loop over wl-points in 2D computational grid
          !
          do n = nfu, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                nd = n - 1
                !
                nm   = kgrpnt(m ,n )
                nmd  = kgrpnt(md,n )
                ndm  = kgrpnt(m ,nd)
                ndmd = kgrpnt(md,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd  == 1 ) nmd  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmd == 1 ) ndmd = ndm
                !
                s1(nm) = 0.25 * ( wlevf(nm,2) + wlevf(nmd,2) + wlevf(ndm,2) + wlevf(ndmd,2) )
                !
             enddo
             !
          enddo
          !
          ! set to zero for permanently dry points
          !
          s1(1) = 0.
          !
       endif
       !
       ! exchange water levels with neighbouring subdomains
       !
       call SWEXCHG ( s1, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize water levels at appropriate boundaries in case of repeating grid
       !
       call periodic ( s1, kgrpnt, 1, 1 )
       !
       initsf = .false.
       !
    endif
    !
    if ( inituf .or. ifldyn(2) == 1 ) then
       !
       ! determine u-velocity in u-points
       !
       if ( oned ) then
          !
          ! loop over u-points in 1D computational grid
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             u1(nm,1) = uxf(nm,2)
             !
          enddo
          !
       else
          !
          ! loop over u-points in 2D computational grid
          !
          if ( lstag(2) ) then
             !
             do n = nfu, nl
                !
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   u1(nm,1) = uxf(nm,2)
                   !
                enddo
                !
             enddo
             !
          else
             !
             do n = nfu, nl
                !
                do m = mf, ml
                   !
                   nd = n - 1
                   !
                   nm  = kgrpnt(m,n )
                   ndm = kgrpnt(m,nd)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( ndm == 1 ) ndm = nm
                   !
                   u1(nm,1) = 0.5 * ( uxf(nm,2) + uxf(ndm,2) ) * ( ycgrid(m,n) - ycgrid(m,n-1) ) - 0.5 * ( uyf(nm,2) + uyf(ndm,2) ) * ( xcgrid(m,n) - xcgrid(m,n-1) )
                   u1(nm,1) = u1(nm,1) / guu(nm)
                   if ( .not.cvleft ) u1(nm,1) = -u1(nm,1)
                   !
                enddo
                !
             enddo
             !
          endif
          !
          ! set to zero for permanently dry points
          !
          u1(1,1) = 0.
          !
       endif
       !
       ! exchange u-velocities with neighbouring subdomains
       !
       call SWEXCHG ( u1(1,1), kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize u-velocities at appropriate boundaries in case of repeating grid
       !
       call periodic ( u1(1,1), kgrpnt, 1, 1 )
       !
       if ( kmax > 1 ) then
          !
          udep(:) = u1(:,1)
          !
          ! expand in vertical direction under the assumption of uniform velocity in vertical
          !
          do k = 2, kmax
             !
             u1(:,k) = udep(:)
             !
          enddo
          !
       endif
       !
       inituf = .false.
       !
    endif
    !
    if ( .not.oned .and. (initvf .or. ifldyn(3) == 1) ) then
       !
       ! determine v-velocity in v-points in 2D computational grid
       !
       if ( lstag(3) ) then
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                v1(nm,1) = uyf(nm,2)
                !
             enddo
             !
          enddo
          !
       else
          !
          do n = nf, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m ,n)
                nmd = kgrpnt(md,n)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd == 1 ) nmd = nm
                !
                v1(nm,1) = 0.5 * ( uyf(nm,2) + uyf(nmd,2) ) * ( xcgrid(m,n) - xcgrid(m-1,n) ) - 0.5 * ( uxf(nm,2) + uxf(nmd,2) ) * ( ycgrid(m,n) - ycgrid(m-1,n) )
                v1(nm,1) = v1(nm,1) / gvv(nm)
                if ( .not.cvleft ) v1(nm,1) = -v1(nm,1)
                !
             enddo
             !
          enddo
          !
       endif
       !
       ! set to zero for permanently dry points
       !
       v1(1,1) = 0.
       !
       ! exchange v-velocities with neighbouring subdomains
       !
       call SWEXCHG ( v1(1,1), kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize v-velocities at appropriate boundaries in case of repeating grid
       !
       call periodic ( v1(1,1), kgrpnt, 1, 1 )
       !
       if ( kmax > 1 ) then
          !
          vdep(:) = v1(:,1)
          !
          ! expand in vertical direction under the assumption of uniform velocity in vertical
          !
          do k = 2, kmax
             !
             v1(:,k) = vdep(:)
             !
          enddo
          !
       endif
       !
       initvf = .false.
       !
    endif
    !
end subroutine SwashUpdFlowFlds
