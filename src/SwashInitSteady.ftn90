subroutine SwashInitSteady
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
!   Initializes velocities based on Chezy formula for uniform flow
!
!   Method
!
!   The Chezy formula used
!
!      u = C * sqrt(h * i)
!
!   where
!      u = velocity
!      C = Chezy friction coefficient
!      h = water depth
!      i = slope of water level
!
!   The initialized velocities are limited by the Froude criterion so that no super-critical flow will occur
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
    integer       :: k        ! loop counter
    integer       :: m        ! loop counter
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nu       ! index of point n+1
    integer       :: nm       ! pointer to m
    integer       :: nmu      ! pointer to m+1
    integer       :: num      ! pointer to n+1
    !
    real          :: cz       ! Chezy value
    real          :: froude   ! Froude number
    real          :: slope    ! slope of water level
    real          :: vel      ! local velocity
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    call strace (ient,'SwashInitSteady')
    !
    if ( oned ) then
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
          slope = ( s1(nm) - s1(nmu) ) / dx
          !
          if ( cfricu(nm) > 0. ) then
             cz = sqrt( grav / cfricu(nm) )
          else
             cz = 0.
          endif
          !
          vel = sign(1.,slope) * cz * sqrt( hu(nm) * abs(slope) )
          !
          if ( hu(nm) > 0. ) then
             !
             froude = abs(vel)/sqrt( grav * hu(nm) )
             !
             if ( froude > 0.9999 ) vel = ( 0.9999 / froude ) * vel
             !
          endif
          !
          u1(nm,1) = vel
          !
       enddo
       !
       ! exchange u-velocities with neighbouring subdomains
       !
       call SWEXCHG ( u1(1,1), kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    else
       !
       ! loop over u-points
       !
       do m = mf, ml
          !
          do n = nfu, nl
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             if ( nmu == 1 ) nmu = nm
             !
             slope = ( s1(nm) - s1(nmu) ) / gvu(nm)
             !
             if ( cfricu(nm) > 0. ) then
                cz = sqrt( grav / cfricu(nm) )
             else
                cz = 0.
             endif
             !
             vel = sign(1.,slope) * cz * sqrt( hu(nm) * abs(slope) )
             !
             if ( hu(nm) > 0. ) then
                !
                froude = abs(vel)/sqrt( grav * hu(nm) )
                !
                if ( froude > 0.9999 ) vel = ( 0.9999 / froude ) * vel
                !
             endif
             !
             u1(nm,1) = vel
             !
          enddo
          !
       enddo
       !
       ! loop over v-points
       !
       do n = nf, nl
          !
          do m = mfu, ml
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             if ( num == 1 ) num = nm
             !
             slope = ( s1(nm) - s1(num) ) / guv(nm)
             !
             if ( cfricv(nm) > 0. ) then
                cz = sqrt( grav / cfricv(nm) )
             else
                cz = 0.
             endif
             !
             vel = sign(1.,slope) * cz * sqrt( hv(nm) * abs(slope) )
             !
             if ( hv(nm) > 0. ) then
                !
                froude = abs(vel)/sqrt( grav * hv(nm) )
                !
                if ( froude > 0.9999 ) vel = ( 0.9999 / froude ) * vel
                !
             endif
             !
             v1(nm,1) = vel
             !
          enddo
          !
       enddo
       !
       ! set to zero for permanently dry points
       !
       u1(1,1) = 0.
       v1(1,1) = 0.
       !
       ! exchange u- and v-velocities with neighbouring subdomains
       !
       call SWEXCHG ( u1(1,1), kgrpnt, 1, 1 )
       call SWEXCHG ( v1(1,1), kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize u- and v-velocities at appropriate boundaries in case of repeating grid
       !
       call periodic ( u1(1,1), kgrpnt, 1, 1 )
       call periodic ( v1(1,1), kgrpnt, 1, 1 )
       !
    endif
    !
    if ( kmax > 1 ) then
       !
       udep(:) = u1(:,1)
       if ( .not.oned ) vdep(:) = v1(:,1)
       !
       ! expand in vertical direction under the assumption of uniform velocity in vertical
       !
       do k = 2, kmax
          !
          u1(:,k) = udep(:)
          if ( .not.oned ) v1(:,k) = vdep(:)
          !
       enddo
       !
    endif
    !
end subroutine SwashInitSteady
