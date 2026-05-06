subroutine SwashFlowDP
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
!   Determines bottom values in water level and velocity points
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
    integer       :: n        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mu       ! index of point m+1
    integer       :: nd       ! index of point n-1
    integer       :: nu       ! index of point n+1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: ndmd     ! pointer to m-1,n-1
    integer       :: nfm      ! pointer to m,nf
    integer       :: nfum     ! pointer to m,nfu
    integer       :: nlm      ! pointer to m,nl
    integer       :: nlum     ! pointer to m,nlu
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmf      ! pointer to mf,n
    integer       :: nmfu     ! pointer to mfu,n
    integer       :: nml      ! pointer to ml,n
    integer       :: nmlu     ! pointer to mlu,n
    integer       :: nmu      ! pointer to m+1,n
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
    if (ltrace) call strace (ient,'SwashFlowDP')
    !
    ! first determine bottom values in cells by means of user-defined input
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
          if ( dpsopt == 1 ) then
             dps(nm) = min( depf(nmd), depf(nm) )
          else if (dpsopt == 2 ) then
             dps(nm) = 0.5 * ( depf(nm) + depf(nmd) )
          else if ( dpsopt == 3 ) then
             dps(nm) = max( depf(nmd), depf(nm) )
          else if ( dpsopt == 4 ) then
             dps(nm) = depf(nm)
          endif
          !
          if ( s1(nm) < epsdry - dps(nm) ) then
             !
             s1(nm) = 0.99*epsdry - dps(nm)
             !
          endif
          !
       enddo
       !
       nmf  = kgrpnt(mf ,1)
       nmfu = kgrpnt(mfu,1)
       nml  = kgrpnt(ml ,1)
       nmlu = kgrpnt(mlu,1)
       !
       ! set bottom values and water levels in virtual wl-points
       !
       dps(nmf ) = dps(nmfu)
       dps(nmlu) = dps(nml )
       !
       s1(nmf ) = s1(nmfu)
       s1(nmlu) = s1(nml )
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
             if ( dpsopt == 1 ) then
                dps(nm) = min( depf(nmd), depf(nm), depf(ndmd), depf(ndm) )
             else if ( dpsopt == 2 ) then
                dps(nm) = 0.25 * ( depf(nm) + depf(nmd) + depf(ndm) + depf(ndmd) )
             else if ( dpsopt == 3 ) then
                dps(nm) = max( depf(nmd), depf(nm), depf(ndmd), depf(ndm) )
             else if ( dpsopt == 4 ) then
                dps(nm) = depf(nm)
             endif
             !
             if ( s1(nm) < epsdry - dps(nm) ) then
                !
                s1(nm) = 0.99*epsdry - dps(nm)
                !
             endif
             !
          enddo
          !
          nmf  = kgrpnt(mf ,n)
          nmfu = kgrpnt(mfu,n)
          nml  = kgrpnt(ml ,n)
          nmlu = kgrpnt(mlu,n)
          !
          ! set bottom values and water levels in virtual wl-points along left and right boundaries
          !
          if ( .not.lreptx ) then
             !
             dps(nmf ) = dps(nmfu)
             dps(nmlu) = dps(nml )
             !
             s1(nmf ) = s1(nmfu)
             s1(nmlu) = s1(nml )
             !
          endif
          !
       enddo
       !
       do m = mf, mlu
          !
          nfm  = kgrpnt(m,nf )
          nfum = kgrpnt(m,nfu)
          nlm  = kgrpnt(m,nl )
          nlum = kgrpnt(m,nlu)
          !
          ! set bottom values and water levels in virtual wl-points along lower and upper boundaries (including corner points)
          !
          if ( .not.lrepty ) then
             !
             dps(nfm ) = dps(nfum)
             dps(nlum) = dps(nlm )
             !
             s1(nfm ) = s1(nfum)
             s1(nlum) = s1(nlm )
             !
          endif
          !
       enddo
       !
       ! set to zero for permanently dry points
       !
       dps(1) = 0.
       s1 (1) = 0.
       !
    endif
    !
    ! exchange bottom values and water levels with neighbouring subdomains
    !
    call SWEXCHG ( dps, kgrpnt, 1, 1 )
    call SWEXCHG (  s1, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize bottom values and water levels at appropriate boundaries in case of repeating grid
    !
    call periodic ( dps, kgrpnt, 1, 1 )
    call periodic (  s1, kgrpnt, 1, 1 )
    !
    ! next determine bottom values in velocity points using the tiled approach
    !
    if ( oned ) then
       !
       ! loop over u-points in 1D computational grid
       !
       do m = mf, ml
          !
          mu = m + 1
          !
          nm  = kgrpnt(m ,1)
          nmu = kgrpnt(mu,1)
          !
          dpu(nm) = min( dps(nm), dps(nmu) )
          !
       enddo
       !
       ! exchange values with neighbouring subdomains
       !
       call SWEXCHG ( dpu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    else
       !
       ! loop over u-points in 2D computational grid
       !
       do n = nf, nlu
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,n)
             nmu = kgrpnt(mu,n)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmu == 1 ) nmu = nm
             !
             dpu(nm) = min( dps(nm), dps(nmu) )
             !
          enddo
       enddo
       !
       ! loop over v-points in 2D computational grid
       !
       do m = mf, mlu
          do n = nf, nl
             !
             nu = n + 1
             !
             nm  = kgrpnt(m,n )
             num = kgrpnt(m,nu)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( num == 1 ) num = nm
             !
             dpv(nm) = min( dps(nm), dps(num) )
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       dpu(1) = 0.
       dpv(1) = 0.
       !
       ! exchange values with neighbouring subdomains
       !
       call SWEXCHG ( dpu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( dpv, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize values at appropriate boundaries in case of repeating grid
       !
       call periodic ( dpu, kgrpnt, 1, 1 )
       call periodic ( dpv, kgrpnt, 1, 1 )
       !
    endif
    !
end subroutine SwashFlowDP
