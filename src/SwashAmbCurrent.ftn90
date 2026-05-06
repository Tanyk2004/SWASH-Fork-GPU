subroutine SwashAmbCurrent
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
!    1.00, July 2023: New subroutine
!
!   Purpose
!
!   Determines ambient current in water level or velocity points
!
!   Modules used
!
    use ocpcomm4
    use m_genarr
    use SwashCommdata2
    use SwashCommdata3
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
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
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
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashAmbCurrent')
    !
    if ( iamb == 1 ) then
       !
       ! determine ambient current in water level points
       !
       if ( oned ) then
          !
          ! determine ambient u-velocity in wl-points in 1D computational grid
          !
          if ( varva ) then
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nm  = kgrpnt(m ,1)
                nmd = kgrpnt(md,1)
                !
                Ubck(nm,:) = 0.5 * ( avxf(nm,:) + avxf(nmd,:) )
                !
             enddo
             !
             nmf  = kgrpnt(mf ,1)
             nmfu = kgrpnt(mfu,1)
             nml  = kgrpnt(ml ,1)
             nmlu = kgrpnt(mlu,1)
             !
             ! set ambient u-velocity in virtual wl-points
             !
             Ubck(nmf ,:) = Ubck(nmfu,:)
             Ubck(nmlu,:) = Ubck(nml ,:)
             !
          else
             !
             Ubck = uamb
             !
          endif
          !
          ! exchange ambient u-velocity with neighbouring subdomains
          !
          call SWEXCHG ( Ubck, kgrpnt, 1, kpmax )
          if (STPNOW()) return
          !
       else
          !
          ! determine ambient u-velocity in wl-points in 2D computational grid
          !
          if ( lstag(20) ) then
             !
             do n = nfu, nl
                !
                do m = mfu, ml
                   !
                   md = m - 1
                   !
                   nm  = kgrpnt(m ,n)
                   nmd = kgrpnt(md,n)
                   !
                   Ubck(nm,:) = 0.5 * ( avxf(nm,:) + avxf(nmd,:) )
                   !
                enddo
                !
                nmf  = kgrpnt(mf ,n)
                nmfu = kgrpnt(mfu,n)
                nml  = kgrpnt(ml ,n)
                nmlu = kgrpnt(mlu,n)
                !
                ! set ambient u-velocity in virtual wl-points along left and right boundaries
                !
                if ( .not.lreptx ) then
                   !
                   Ubck(nmf ,:) = Ubck(nmfu,:)
                   Ubck(nmlu,:) = Ubck(nml ,:)
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
                ! set ambient u-velocity in virtual wl-points along lower and upper boundaries (including corner points)
                !
                if ( .not.lrepty ) then
                   !
                   Ubck(nfm ,:) = Ubck(nfum,:)
                   Ubck(nlum,:) = Ubck(nlm ,:)
                   !
                endif
                !
             enddo
             !
          else if ( varva ) then
             !
             ! first, transform (Cartesian to contravariant) and interpolate to u-points
             !
             do n = nfu, nl
                !
                nd = n - 1
                !
                do m = mf, ml
                   !
                   nm  = kgrpnt(m,n )
                   ndm = kgrpnt(m,nd)
                   !
                   wrk(nm,:) = 0.5 * ( avxf(nm,:) + avxf(ndm,:) ) * ( ycgrid(m,n) - ycgrid(m,n-1) ) - 0.5 * ( avyf(nm,:) + avyf(ndm,:) ) * ( xcgrid(m,n) - xcgrid(m,n-1) )
                   wrk(nm,:) = wrk(nm,:) / guu(nm)
                   !
                enddo
                !
             enddo
             !
             ! next, interpolate to cell centers
             !
             do n = nfu, nl
                !
                do m = mfu, ml
                   !
                   md = m - 1
                   !
                   nm  = kgrpnt(m ,n)
                   nmd = kgrpnt(md,n)
                   !
                   Ubck(nm,:) = 0.5 * ( wrk(nm,:) + wrk(nmd,:) )
                   !
                enddo
                !
                nmf  = kgrpnt(mf ,n)
                nmfu = kgrpnt(mfu,n)
                nml  = kgrpnt(ml ,n)
                nmlu = kgrpnt(mlu,n)
                !
                ! set ambient u-velocity in virtual wl-points along left and right boundaries
                !
                if ( .not.lreptx ) then
                   !
                   Ubck(nmf ,:) = Ubck(nmfu,:)
                   Ubck(nmlu,:) = Ubck(nml ,:)
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
                ! set ambient u-velocity in virtual wl-points along lower and upper boundaries (including corner points)
                !
                if ( .not.lrepty ) then
                   !
                   Ubck(nfm ,:) = Ubck(nfum,:)
                   Ubck(nlum,:) = Ubck(nlm ,:)
                   !
                endif
                !
             enddo
             !
          else
             !
             ! first, transform (Cartesian to contravariant) and interpolate to u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   work(nm,1) = uamb * ( ycgrid(m,n) - ycgrid(m,n-1) ) - vamb * ( xcgrid(m,n) - xcgrid(m,n-1) )
                   work(nm,1) = work(nm,1) / guu(nm)
                   !
                enddo
             enddo
             !
             ! next, interpolate to cell centers
             !
             do n = nfu, nl
                !
                do m = mfu, ml
                   !
                   md = m - 1
                   !
                   nm  = kgrpnt(m ,n)
                   nmd = kgrpnt(md,n)
                   !
                   Ubck(nm,:) = 0.5 * ( work(nm,1) + work(nmd,1) )
                   !
                enddo
                !
                nmf  = kgrpnt(mf ,n)
                nmfu = kgrpnt(mfu,n)
                nml  = kgrpnt(ml ,n)
                nmlu = kgrpnt(mlu,n)
                !
                ! set ambient u-velocity in virtual wl-points along left and right boundaries
                !
                if ( .not.lreptx ) then
                   !
                   Ubck(nmf ,:) = Ubck(nmfu,:)
                   Ubck(nmlu,:) = Ubck(nml ,:)
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
                ! set ambient u-velocity in virtual wl-points along lower and upper boundaries (including corner points)
                !
                if ( .not.lrepty ) then
                   !
                   Ubck(nfm ,:) = Ubck(nfum,:)
                   Ubck(nlum,:) = Ubck(nlm ,:)
                   !
                endif
                !
             enddo
             !
          endif
          !
          ! determine ambient v-velocity in wl-points in 2D computational grid
          !
          if ( lstag(21) ) then
             !
             do n = nfu, nl
                !
                nd = n - 1
                !
                do m = mfu, ml
                   !
                   nm  = kgrpnt(m,n )
                   ndm = kgrpnt(m,nd)
                   !
                   Vbck(nm,:) = 0.5 * ( avyf(nm,:) + avyf(ndm,:) )
                   !
                enddo
                !
                nmf  = kgrpnt(mf ,n)
                nmfu = kgrpnt(mfu,n)
                nml  = kgrpnt(ml ,n)
                nmlu = kgrpnt(mlu,n)
                !
                ! set ambient v-velocity in virtual wl-points along left and right boundaries
                !
                if ( .not.lreptx ) then
                   !
                   Vbck(nmf ,:) = Vbck(nmfu,:)
                   Vbck(nmlu,:) = Vbck(nml ,:)
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
                ! set ambient v-velocity in virtual wl-points along lower and upper boundaries (including corner points)
                !
                if ( .not.lrepty ) then
                   !
                   Vbck(nfm ,:) = Vbck(nfum,:)
                   Vbck(nlum,:) = Vbck(nlm ,:)
                   !
                endif
                !
             enddo
             !
          else if ( varva ) then
             !
             ! first, transform (Cartesian to contravariant) and interpolate to v-points
             !
             do n = nf, nl
                do m = mfu, ml
                   !
                   md = m - 1
                   !
                   nm  = kgrpnt(m ,n)
                   nmd = kgrpnt(md,n)
                   !
                   wrk(nm,:) = 0.5 * ( avyf(nm,:) + avyf(nmd,:) ) * ( xcgrid(m,n) - xcgrid(m-1,n) ) - 0.5 * ( avxf(nm,:) + avxf(nmd,:) ) * ( ycgrid(m,n) - ycgrid(m-1,n) )
                   wrk(nm,:) = wrk(nm,:) / gvv(nm)
                   !
                enddo
             enddo
             !
             ! next, interpolate to cell centers
             !
             do n = nfu, nl
                !
                nd = n - 1
                !
                do m = mfu, ml
                   !
                   nm  = kgrpnt(m,n )
                   ndm = kgrpnt(m,nd)
                   !
                   Vbck(nm,:) = 0.5 * ( wrk(nm,:) + wrk(ndm,:) )
                   !
                enddo
                !
                nmf  = kgrpnt(mf ,n)
                nmfu = kgrpnt(mfu,n)
                nml  = kgrpnt(ml ,n)
                nmlu = kgrpnt(mlu,n)
                !
                ! set ambient v-velocity in virtual wl-points along left and right boundaries
                !
                if ( .not.lreptx ) then
                   !
                   Vbck(nmf ,:) = Vbck(nmfu,:)
                   Vbck(nmlu,:) = Vbck(nml ,:)
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
                ! set ambient v-velocity in virtual wl-points along lower and upper boundaries (including corner points)
                !
                if ( .not.lrepty ) then
                   !
                   Vbck(nfm ,:) = Vbck(nfum,:)
                   Vbck(nlum,:) = Vbck(nlm ,:)
                   !
                endif
                !
             enddo
             !
          else
             !
             ! first, transform (Cartesian to contravariant) and interpolate to v-points
             !
             do n = nf, nl
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   work(nm,1) = vamb * ( xcgrid(m,n) - xcgrid(m-1,n) ) - uamb * ( ycgrid(m,n) - ycgrid(m-1,n) )
                   work(nm,1) = work(nm,1) / gvv(nm)
                   !
                enddo
             enddo
             !
             ! next, interpolate to cell centers
             !
             do n = nfu, nl
                !
                nd = n - 1
                !
                do m = mfu, ml
                   !
                   nm  = kgrpnt(m,n )
                   ndm = kgrpnt(m,nd)
                   !
                   Vbck(nm,:) = 0.5 * ( work(nm,1) + work(ndm,1) )
                   !
                enddo
                !
                nmf  = kgrpnt(mf ,n)
                nmfu = kgrpnt(mfu,n)
                nml  = kgrpnt(ml ,n)
                nmlu = kgrpnt(mlu,n)
                !
                ! set ambient v-velocity in virtual wl-points along left and right boundaries
                !
                if ( .not.lreptx ) then
                   !
                   Vbck(nmf ,:) = Vbck(nmfu,:)
                   Vbck(nmlu,:) = Vbck(nml ,:)
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
                ! set ambient v-velocity in virtual wl-points along lower and upper boundaries (including corner points)
                !
                if ( .not.lrepty ) then
                   !
                   Vbck(nfm ,:) = Vbck(nfum,:)
                   Vbck(nlum,:) = Vbck(nlm ,:)
                   !
                endif
                !
             enddo
             !
          endif
          !
          ! set to zero for permanently dry points
          !
          Ubck(1,:) = 0.
          Vbck(1,:) = 0.
          !
          ! exchange ambient u- and v-velocities with neighbouring subdomains
          !
          call SWEXCHG ( Ubck, kgrpnt, 1, kpmax )
          if (STPNOW()) return
          call SWEXCHG ( Vbck, kgrpnt, 1, kpmax )
          if (STPNOW()) return
          !
          ! synchronize ambient u- and v-velocities at appropriate boundaries in case of repeating grid
          !
          call periodic ( Ubck, kgrpnt, 1, kpmax )
          call periodic ( Vbck, kgrpnt, 1, kpmax )
          !
       endif
       !
    else if ( iamb == 2 ) then
       !
       ! determine ambient current in velocity points
       !
       if ( oned ) then
          !
          ! determine ambient u-velocity in u-points in 1D computational grid
          !
          if ( varva ) then
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                Ubck(nm,:) = avxf(nm,:)
                !
             enddo
             !
          else
             !
             Ubck = uamb
             !
          endif
          !
          ! exchange ambient u-velocity with neighbouring subdomains
          !
          call SWEXCHG ( Ubck, kgrpnt, 1, kpmax )
          if (STPNOW()) return
          !
       else
          !
          ! determine ambient u-velocity in u-points in 2D computational grid
          !
          if ( lstag(20) ) then
             !
             do n = nfu, nl
                !
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   Ubck(nm,:) = avxf(nm,:)
                   !
                enddo
                !
             enddo
             !
          else if ( varva ) then
             !
             do n = nfu, nl
                !
                nd = n - 1
                !
                do m = mf, ml
                   !
                   nm  = kgrpnt(m,n )
                   ndm = kgrpnt(m,nd)
                   !
                   Ubck(nm,:) = 0.5 * ( avxf(nm,:) + avxf(ndm,:) ) * ( ycgrid(m,n) - ycgrid(m,n-1) ) - 0.5 * ( avyf(nm,:) + avyf(ndm,:) ) * ( xcgrid(m,n) - xcgrid(m,n-1) )
                   Ubck(nm,:) = Ubck(nm,:) / guu(nm)
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
                   nm = kgrpnt(m,n)
                   !
                   Ubck(nm,:) = uamb * ( ycgrid(m,n) - ycgrid(m,n-1) ) - vamb * ( xcgrid(m,n) - xcgrid(m,n-1) )
                   Ubck(nm,:) = Ubck(nm,:) / guu(nm)
                   !
                enddo
                !
             enddo
             !
          endif
          !
          ! determine ambient v-velocity in v-points in 2D computational grid
          !
          if ( lstag(21) ) then
             !
             do n = nf, nl
                !
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   Vbck(nm,:) = avyf(nm,:)
                   !
                enddo
                !
             enddo
             !
          else if ( varva ) then
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
                   Vbck(nm,:) = 0.5 * ( avyf(nm,:) + avyf(nmd,:) ) * ( xcgrid(m,n) - xcgrid(m-1,n) ) - 0.5 * ( avxf(nm,:) + avxf(nmd,:) ) * ( ycgrid(m,n) - ycgrid(m-1,n) )
                   Vbck(nm,:) = Vbck(nm,:) / gvv(nm)
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
                   nm = kgrpnt(m,n)
                   !
                   Vbck(nm,:) = vamb * ( xcgrid(m,n) - xcgrid(m-1,n) ) - uamb * ( ycgrid(m,n) - ycgrid(m-1,n) )
                   Vbck(nm,:) = Vbck(nm,:) / gvv(nm)
                   !
                enddo
                !
             enddo
             !
          endif
          !
          ! set to zero for permanently dry points
          !
          Ubck(1,:) = 0.
          Vbck(1,:) = 0.
          !
          ! exchange ambient u- and v-velocities with neighbouring subdomains
          !
          call SWEXCHG ( Ubck, kgrpnt, 1, kpmax )
          if (STPNOW()) return
          call SWEXCHG ( Vbck, kgrpnt, 1, kpmax )
          if (STPNOW()) return
          !
          ! synchronize ambient u- and v-velocities at appropriate boundaries in case of repeating grid
          !
          call periodic ( Ubck, kgrpnt, 1, kpmax )
          call periodic ( Vbck, kgrpnt, 1, kpmax )
          !
       endif
       !
    endif
    !
end subroutine SwashAmbCurrent
