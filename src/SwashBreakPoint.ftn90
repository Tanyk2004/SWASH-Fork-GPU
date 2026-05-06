subroutine SwashBreakPoint
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
!    1.00: Pieter Smit
!    6.01: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00,     May 2011: New subroutine
!    1.05, January 2012: extension to correct amount of energy dissipation
!    6.01,    June 2019: extension internal wave generation
!
!   Purpose
!
!   Determines grid points where steep bore-like wave front occurs for wave breaking
!
!   Method
!
!   Update mask array at those points where the vertical speed of the free surface
!   exceeds a fraction of the wave phase speed. At those points, hydrostatic pressure
!   is assumed and remains so at the front face of the breaking wave.
!
!   This approach combined with a proper momentum conservation leads to a correct
!   amount of energy dissipation on the front face of the breaking wave. Moreover,
!   nonlinear wave properties such as asymmetry and skewness are preserved as well.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
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
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
    !
    real          :: alpha    ! threshold parameter corresponding to onset of wave breaking
    real          :: beta     ! threshold parameter corresponding to re-initiation of breaking in post-breaking area
    real          :: dsdt     ! vertical velocity of the free surface
    real          :: rootgh   ! shallow water celerity
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashBreakPoint')
    !
    alpha = psurf(1)
    beta  = psurf(2)
    !
    ! store array brks temporarily
    !
    iwrk(:,1) = brks(:)
    !
    if ( oned ) then
       !
       do m = mfu, ml
          !
          md = m - 1
          mu = m + 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          nmu = kgrpnt(mu,1)
          !
          if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
             !
             rootgh = sqrt( grav * hs(nm) )
             !
             dsdt = ( s1(nm) - s0(nm) ) / dt
             !
             if ( dsdt > alpha * rootgh ) then
                !
                brks(nm) = 1
                !
                q(nm,:) = 0.
                !
             else if ( dsdt > beta * rootgh .and. (iwrk(nmd,1) == 1 .or. iwrk(nmu,1) == 1) ) then
                !
                brks(nm) = 1
                !
                q(nm,:) = 0.
                !
             else if ( dsdt > 0. .and. brks(nm) == 1 ) then
                !
                q(nm,:) = 0.
                !
             else
                !
                brks(nm) = 0
                !
             endif
             !
          else
             !
             brks(nm) = 0
             !
          endif
          !
       enddo
       !
    else
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             nmu = kgrpnt(mu,n )
             ndm = kgrpnt(m ,nd)
             num = kgrpnt(m ,nu)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmd == 1 ) nmd = nm
             if ( nmu == 1 ) nmu = nm
             if ( ndm == 1 ) ndm = nm
             if ( num == 1 ) num = nm
             !
             if ( wets(nm) == 1 .and. presp(nm) == 0 ) then
                !
                rootgh = sqrt( grav * hs(nm) )
                !
                dsdt = ( s1(nm) - s0(nm) ) / dt
                !
                if ( dsdt > alpha * rootgh ) then
                   !
                   brks(nm) = 1
                   !
                   q(nm,:) = 0.
                   !
                else if ( dsdt > beta * rootgh .and. (iwrk(nmd,1) == 1 .or. iwrk(nmu,1) == 1 .or. iwrk(ndm,1) == 1 .or. iwrk(num,1) == 1) ) then
                   !
                   brks(nm) = 1
                   !
                   q(nm,:) = 0.
                   !
                else if ( dsdt > 0. .and. brks(nm) == 1 ) then
                   !
                   q(nm,:) = 0.
                   !
                else
                   !
                   brks(nm) = 0
                   !
                endif
                !
             else
                !
                brks(nm) = 0
                !
             endif
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       brks(1) = 0
       !
    endif
    !
    ! breaking deactivated inside the source region
    !
    if ( iwvgen /= 0 ) then
       !
       do nm = 2, mcgrd
          !
          if ( srcm(nm) /= 0. ) brks(nm) = 0
          !
       enddo
       !
    endif
    !
    ! exchange mask values with neighbouring subdomains
    !
    call SWEXCHGI ( brks, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize mask values at appropriate boundaries in case of repeating grid
    !
    call periodici ( brks, kgrpnt, 1, 1 )
    !
    ! re-update mask array for wetting and drying at water level points by taking into account the breaking points
    !
    do nm = 1, mcgrd
       !
       if ( wets(nm) == 1 .and. brks(nm) == 1 ) wets(nm) = 0
       !
    enddo
    !
end subroutine SwashBreakPoint
