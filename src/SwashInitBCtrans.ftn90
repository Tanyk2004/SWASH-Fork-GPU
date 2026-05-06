subroutine SwashInitBCtrans
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
!    1.00, July 2013: New subroutine
!
!   Purpose
!
!   Initializes transport constituents based on input fields and stores boundary values as Dirichlet type condition
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr
    use m_parall
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
    integer       :: ndmd     ! pointer to m-1,n-1
    integer       :: ndmf     ! pointer to mf,n-1
    integer       :: ndml     ! pointer to ml,n-1
    integer       :: nfm      ! pointer to m,nf
    integer       :: nfmd     ! pointer to m-1,nf
    integer       :: nfum     ! pointer to m,nfu
    integer       :: nlm      ! pointer to m,nl
    integer       :: nlmd     ! pointer to m-1,nl
    integer       :: nlum     ! pointer to m,nlu
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmf      ! pointer to mf,n
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
    if (ltrace) call strace (ient,'SwashInitBCtrans')
    !
    if ( lsal > 0 ) then
       !
       ! specify salinity in wl-points based on input field
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
             rp(nm,:,lsal) = 0.5 * ( salf(nm,:) + salf(nmd,:) )
             !
          enddo
          !
          ! exchange salinity with neighbouring subdomains
          !
          call SWEXCHG ( rp(1,1,lsal), kgrpnt, 1, kmax )
          if (STPNOW()) return
          !
          ! store boundary values at real boundaries
          !
          nmf  = kgrpnt(mf ,1)
          nml  = kgrpnt(ml ,1)
          nmlu = kgrpnt(mlu,1)
          !
          if ( LMXF ) then
             cbndl(1  ,:,lsal) = salf(nmf,:)
             rp   (nmf,:,lsal) = cbndl(1,:,lsal)
          endif
          !
          if ( LMXL ) then
             cbndr(1   ,:,lsal) = salf(nml,:)
             rp   (nmlu,:,lsal) = cbndr(1,:,lsal)
          endif
          !
       else
          !
          ! loop over wl-points in 2D computational grid
          !
          do n = nfu, nl
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
                rp(nm,:,lsal) = 0.25 * ( salf(nm,:) + salf(nmd,:) + salf(ndm,:) + salf(ndmd,:) )
                !
             enddo
          enddo
          !
          ! set to zero for permanently dry points
          !
          rp(1,:,lsal) = 0.
          !
          ! exchange salinity with neighbouring subdomains
          !
          call SWEXCHG ( rp(1,1,lsal), kgrpnt, 1, kmax )
          if (STPNOW()) return
          !
          ! store boundary values at real non-periodic boundaries
          !
          if ( .not.lreptx ) then
             !
             do n = nfu, nl
                !
                nd = n - 1
                !
                nmf  = kgrpnt(mf ,n )
                ndmf = kgrpnt(mf ,nd)
                nml  = kgrpnt(ml ,n )
                ndml = kgrpnt(ml ,nd)
                nmlu = kgrpnt(mlu,n )
                !
                if ( LMXF ) then
                   !
                   cbndl(n  ,:,lsal) = 0.5 * ( salf(nmf,:) + salf(ndmf,:) )
                   rp   (nmf,:,lsal) = cbndl(n,:,lsal)
                   !
                endif
                !
                if ( LMXL ) then
                   !
                   cbndr(n   ,:,lsal) = 0.5 * ( salf(nml,:) + salf(ndml,:) )
                   rp   (nmlu,:,lsal) = cbndr(n,:,lsal)
                   !
                endif
                !
             enddo
             !
          else
             !
             call periodic ( rp(1,1,lsal), kgrpnt, 1, kmax )
             !
          endif
          !
          if ( .not.lrepty ) then
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nfm  = kgrpnt(m ,nf )
                nfmd = kgrpnt(md,nf )
                nlm  = kgrpnt(m ,nl )
                nlmd = kgrpnt(md,nl )
                nlum = kgrpnt(m ,nlu)
                !
                if ( LMYF ) then
                   !
                   cbndb(m  ,:,lsal) = 0.5 * ( salf(nfm,:) + salf(nfmd,:) )
                   rp   (nfm,:,lsal) = cbndb(m,:,lsal)
                   !
                endif
                !
                if ( LMYL ) then
                   !
                   cbndt(m   ,:,lsal) = 0.5 * ( salf(nlm,:) + salf(nlmd,:) )
                   rp   (nlum,:,lsal) = cbndt(m,:,lsal)
                   !
                endif
                !
             enddo
             !
          else
             !
             call periodic ( rp(1,1,lsal), kgrpnt, 1, kmax )
             !
          endif
          !
       endif
       !
    endif
    !
    if ( ltemp > 0 ) then
       !
       ! specify temperature in wl-points based on input field
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
             rp(nm,:,ltemp) = 0.5 * ( tempf(nm,:) + tempf(nmd,:) )
             !
          enddo
          !
          ! exchange temperature with neighbouring subdomains
          !
          call SWEXCHG ( rp(1,1,ltemp), kgrpnt, 1, kmax )
          if (STPNOW()) return
          !
          ! store boundary values at real boundaries
          !
          nmf  = kgrpnt(mf ,1)
          nml  = kgrpnt(ml ,1)
          nmlu = kgrpnt(mlu,1)
          !
          if ( LMXF ) then
             cbndl(1  ,:,ltemp) = tempf(nmf,:)
             rp   (nmf,:,ltemp) = cbndl(1,:,ltemp)
          endif
          !
          if ( LMXL ) then
             cbndr(1   ,:,ltemp) = tempf(nml,:)
             rp   (nmlu,:,ltemp) = cbndr(1,:,ltemp)
          endif
          !
       else
          !
          ! loop over wl-points in 2D computational grid
          !
          do n = nfu, nl
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
                rp(nm,:,ltemp) = 0.25 * ( tempf(nm,:) + tempf(nmd,:) + tempf(ndm,:) + tempf(ndmd,:) )
                !
             enddo
          enddo
          !
          ! set to zero for permanently dry points
          !
          rp(1,:,ltemp) = 0.
          !
          ! exchange temperature with neighbouring subdomains
          !
          call SWEXCHG ( rp(1,1,ltemp), kgrpnt, 1, kmax )
          if (STPNOW()) return
          !
          ! store boundary values at real non-periodic boundaries
          !
          if ( .not.lreptx ) then
             !
             do n = nfu, nl
                !
                nd = n - 1
                !
                nmf  = kgrpnt(mf ,n )
                ndmf = kgrpnt(mf ,nd)
                nml  = kgrpnt(ml ,n )
                ndml = kgrpnt(ml ,nd)
                nmlu = kgrpnt(mlu,n )
                !
                if ( LMXF ) then
                   !
                   cbndl(n  ,:,ltemp) = 0.5 * ( tempf(nmf,:) + tempf(ndmf,:) )
                   rp   (nmf,:,ltemp) = cbndl(n,:,ltemp)
                   !
                endif
                !
                if ( LMXL ) then
                   !
                   cbndr(n   ,:,ltemp) = 0.5 * ( tempf(nml,:) + tempf(ndml,:) )
                   rp   (nmlu,:,ltemp) = cbndr(n,:,ltemp)
                   !
                endif
                !
             enddo
             !
          else
             !
             call periodic ( rp(1,1,ltemp), kgrpnt, 1, kmax )
             !
          endif
          !
          if ( .not.lrepty ) then
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nfm  = kgrpnt(m ,nf )
                nfmd = kgrpnt(md,nf )
                nlm  = kgrpnt(m ,nl )
                nlmd = kgrpnt(md,nl )
                nlum = kgrpnt(m ,nlu)
                !
                if ( LMYF ) then
                   !
                   cbndb(m  ,:,ltemp) = 0.5 * ( tempf(nfm,:) + tempf(nfmd,:) )
                   rp   (nfm,:,ltemp) = cbndb(m,:,ltemp)
                   !
                endif
                !
                if ( LMYL ) then
                   !
                   cbndt(m   ,:,ltemp) = 0.5 * ( tempf(nlm,:) + tempf(nlmd,:) )
                   rp   (nlum,:,ltemp) = cbndt(m,:,ltemp)
                   !
                endif
                !
             enddo
             !
          else
             !
             call periodic ( rp(1,1,ltemp), kgrpnt, 1, kmax )
             !
          endif
          !
       endif
       !
    endif
    !
    if ( lsed > 0 ) then
       !
       ! specify volumetric sediment in wl-points based on input field
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
             rp(nm,:,lsed) = 0.5 * ( sedf(nm,:) + sedf(nmd,:) ) / rhos
             !
          enddo
          !
          ! exchange volumetric sediment with neighbouring subdomains
          !
          call SWEXCHG ( rp(1,1,lsed), kgrpnt, 1, kmax )
          if (STPNOW()) return
          !
          ! store boundary values at real boundaries
          !
          nmf  = kgrpnt(mf ,1)
          nml  = kgrpnt(ml ,1)
          nmlu = kgrpnt(mlu,1)
          !
          if ( LMXF ) then
             cbndl(1  ,:,lsed) = sedf(nmf,:) / rhos
             rp   (nmf,:,lsed) = cbndl(1,:,lsed)
          endif
          !
          if ( LMXL ) then
             cbndr(1   ,:,lsed) = sedf(nml,:) / rhos
             rp   (nmlu,:,lsed) = cbndr(1,:,lsed)
          endif
          !
       else
          !
          ! loop over wl-points in 2D computational grid
          !
          do n = nfu, nl
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
                rp(nm,:,lsed) = 0.25 * ( sedf(nm,:) + sedf(nmd,:) + sedf(ndm,:) + sedf(ndmd,:) ) / rhos
                !
             enddo
          enddo
          !
          ! set to zero for permanently dry points
          !
          rp(1,:,lsed) = 0.
          !
          ! exchange volumetric sediment with neighbouring subdomains
          !
          call SWEXCHG ( rp(1,1,lsed), kgrpnt, 1, kmax )
          if (STPNOW()) return
          !
          ! store boundary values at real non-periodic boundaries
          !
          if ( .not.lreptx ) then
             !
             do n = nfu, nl
                !
                nd = n - 1
                !
                nmf  = kgrpnt(mf ,n )
                ndmf = kgrpnt(mf ,nd)
                nml  = kgrpnt(ml ,n )
                ndml = kgrpnt(ml ,nd)
                nmlu = kgrpnt(mlu,n )
                !
                if ( LMXF ) then
                   !
                   cbndl(n  ,:,lsed) = 0.5 * ( sedf(nmf,:) + sedf(ndmf,:) ) / rhos
                   rp   (nmf,:,lsed) = cbndl(n,:,lsed)
                   !
                endif
                !
                if ( LMXL ) then
                   !
                   cbndr(n   ,:,lsed) = 0.5 * ( sedf(nml,:) + sedf(ndml,:) ) / rhos
                   rp   (nmlu,:,lsed) = cbndr(n,:,lsed)
                   !
                endif
                !
             enddo
             !
          else
             !
             call periodic ( rp(1,1,lsed), kgrpnt, 1, kmax )
             !
          endif
          !
          if ( .not.lrepty ) then
             !
             do m = mfu, ml
                !
                md = m - 1
                !
                nfm  = kgrpnt(m ,nf )
                nfmd = kgrpnt(md,nf )
                nlm  = kgrpnt(m ,nl )
                nlmd = kgrpnt(md,nl )
                nlum = kgrpnt(m ,nlu)
                !
                if ( LMYF ) then
                   !
                   cbndb(m  ,:,lsed) = 0.5 * ( sedf(nfm,:) + sedf(nfmd,:) ) / rhos
                   rp   (nfm,:,lsed) = cbndb(m,:,lsed)
                   !
                endif
                !
                if ( LMYL ) then
                   !
                   cbndt(m   ,:,lsed) = 0.5 * ( sedf(nlm,:) + sedf(nlmd,:) ) / rhos
                   rp   (nlum,:,lsed) = cbndt(m,:,lsed)
                   !
                endif
                !
             enddo
             !
          else
             !
             call periodic ( rp(1,1,lsed), kgrpnt, 1, kmax )
             !
          endif
          !
       endif
       !
    endif
    !
end subroutine SwashInitBCtrans
