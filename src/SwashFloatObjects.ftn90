subroutine SwashFloatObjects
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
!    1.00,  August 2014: New subroutine
!    9.01, October 2022: extension moving rigid bodies
!
!   Purpose
!
!   Determines some parameters for floating objects
!
!   Method
!
!   The following parameters are computed:
!
!   1) drafts of floating objects in water level and velocity points
!   2) labels of moving rigid bodies in water level points
!   3) areas, volumes, second area moments and centers of buoyancy
!   4) pretensions of mooring lines
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr, only: kgrpnt, guu, gvv, xcgrid, ycgrid, flobjf, lbodf, iwrk
    use m_parall
    use SwashFlowdata
    use SwashRigBoddata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: l        ! label index
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: mu       ! index of point m+1
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: ndmd     ! pointer to m-1,n-1
    integer       :: nend     ! end index of loop over wl-points in y-direction
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
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
    !
    real          :: area     ! water plane area of rigid body
    real          :: dxl      ! local mesh size in x-direction
    real          :: dyl      ! local mesh size in y-direction
    real          :: fac      ! auxiliary factor
    real          :: gamma    ! inclination of mooring line
    real          :: rx       ! perpendicular distance from the y-axis of rigid body
    real          :: ry       ! perpendicular distance from the x-axis of rigid body
    real          :: vol      ! volume of displacement of rigid body
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashFloatObjects')
    !
    ! first determine drafts in cells by means of user-defined input
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
             flos(nm) = min( flobjf(nmd), flobjf(nm) )
          else if (dpsopt == 2 ) then
             flos(nm) = 0.5 * ( flobjf(nm) + flobjf(nmd) )
          else if ( dpsopt == 3 ) then
             flos(nm) = max( flobjf(nmd), flobjf(nm) )
          else if ( dpsopt == 4 ) then
             flos(nm) = flobjf(nm)
          endif
          !
       enddo
       !
       nmf  = kgrpnt(mf ,1)
       nmfu = kgrpnt(mfu,1)
       nml  = kgrpnt(ml ,1)
       nmlu = kgrpnt(mlu,1)
       !
       ! set drafts in virtual wl-points
       !
       flos(nmf ) = flos(nmfu)
       flos(nmlu) = flos(nml )
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
                flos(nm) = min( flobjf(nmd), flobjf(nm), flobjf(ndmd), flobjf(ndm) )
             else if ( dpsopt == 2 ) then
                flos(nm) = 0.25 * ( flobjf(nm) + flobjf(nmd) + flobjf(ndm) + flobjf(ndmd) )
             else if ( dpsopt == 3 ) then
                flos(nm) = max( flobjf(nmd), flobjf(nm), flobjf(ndmd), flobjf(ndm) )
             else if ( dpsopt == 4 ) then
                flos(nm) = flobjf(nm)
             endif
             !
          enddo
          !
          nmf  = kgrpnt(mf ,n)
          nmfu = kgrpnt(mfu,n)
          nml  = kgrpnt(ml ,n)
          nmlu = kgrpnt(mlu,n)
          !
          ! set drafts in virtual wl-points along left and right boundaries
          !
          if ( .not.lreptx ) then
             !
             flos(nmf ) = flos(nmfu)
             flos(nmlu) = flos(nml )
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
          ! set drafts in virtual wl-points along lower and upper boundaries (including corner points)
          !
          if ( .not.lrepty ) then
             !
             flos(nfm ) = flos(nfum)
             flos(nlum) = flos(nlm )
             !
          endif
          !
       enddo
       !
       ! set to infinity for permanently dry points
       !
       flos(1) = -9999.
       !
    endif
    !
    ! exchange drafts with neighbouring subdomains
    !
    call SWEXCHG ( flos, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize drafts at appropriate boundaries in case of repeating grid
    !
    call periodic ( flos, kgrpnt, 1, 1 )
    !
    ! next determine drafts in velocity points using the tiled approach
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
          flou(nm) = max( flos(nm), flos(nmu) )
          !
       enddo
       !
       ! exchange values with neighbouring subdomains
       !
       call SWEXCHG ( flou, kgrpnt, 1, 1 )
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
             flou(nm) = max( flos(nm), flos(nmu) )
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
             flov(nm) = max( flos(nm), flos(num) )
             !
          enddo
       enddo
       !
       ! set to infinity for permanently dry points
       !
       flou(1) = -9999.
       flov(1) = -9999.
       !
       ! exchange values with neighbouring subdomains
       !
       call SWEXCHG ( flou, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( flov, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize values at appropriate boundaries in case of repeating grid
       !
       call periodic ( flou, kgrpnt, 1, 1 )
       call periodic ( flov, kgrpnt, 1, 1 )
       !
    endif
    !
    ! check depth beneath floating object
    !
    if ( any( dps < flos ) ) call msgerr (3, 'depth beneath floating object is negative; check input grid of FLOAT')
    !
    ! determine labels of floating bodies
    !
    if ( ifloat == 1 ) then
       !
       lfbs = 1    ! is superfluous, but used to compute hydrodynamic loads as output
       !
    else if ( ifloat == 2 ) then
       !
       lfbs = int(lbodf)
       !
    endif
    !
    ! exchange labels with neighbouring subdomains
    !
    call SWEXCHG ( lfbs, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize labels at appropriate boundaries in case of repeating grid
    !
    call periodic ( lfbs, kgrpnt, 1, 1 )
    !
    ! compute some properties of floating objects for the sake of computing hydrostatic restoring forces
    !
    if ( ifloat == 2 ) then
       !
       barea = 0.
       bvol  = 0.
       bmoax = 0.
       bmoay = 0.
       bcob  = 0.
       !
       if ( oned ) then
          !
          ! loop over wl-points in 1D computational grid
          !
          ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
          !
          mend = ml - 1
          if ( LMXL ) mend = ml
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m,1)
             !
             l = lfbs(nm)
             !
             if ( l > 0 ) then
                !
                ! area and volume of body
                area = dx
                vol  = flos(nm) * area
                !
                barea(l) = barea(l) + area
                bvol (l) = bvol (l) + vol
                !
                ! distance from the y-axis
                rx = 0.5 * ( xcgrid(m,1) + xcgrid(m-1,1) ) - bcog(l,1)
                !
                ! second moment of area with respect to y-axis
                bmoay(l) = bmoay(l) + rx*rx * area
                !
                ! center of buoyancy relative to z = 0
                bcob(l) = bcob(l) - 0.5 * flos(nm) * vol
                !
             endif
             !
          enddo
          !
       else
          !
          ! loop over wl-points in 2D computational grid
          !
          ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
          !
          mend = ml - 1
          if ( LMXL ) mend = ml
          !
          nend = nl - 1
          if ( LMYL ) nend = nl
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd == 1 ) nmd = nm
                if ( ndm == 1 ) ndm = nm
                !
                dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
                dyl = 0.5 * ( guu(nm) + guu(nmd) )
                !
                l = lfbs(nm)
                !
                if ( l > 0 ) then
                   !
                   ! area and volume of body
                   area = dxl * dyl
                   vol  = flos(nm) * area
                   !
                   barea(l) = barea(l) + area
                   bvol (l) = bvol (l) + vol
                   !
                   ! distance from the y- and x-axis, respectively
                   rx = 0.5 * ( xcgrid(m,n) + xcgrid(md,n ) ) - bcog(l,1)
                   ry = 0.5 * ( ycgrid(m,n) + ycgrid(m ,nd) ) - bcog(l,2)
                   !
                   ! second moment of area with respect to x- and y-axis, respectively
                   bmoax(l) = bmoax(l) + ry*ry * area
                   bmoay(l) = bmoay(l) + rx*rx * area
                   !
                   ! center of buoyancy relative to z = 0
                   bcob(l) = bcob(l) - 0.5 * flos(nm) * vol
                   !
                endif
                !
             enddo
          enddo
          !
       endif
       !
       ! accumulate body properties over all subdomains
       !
       call SWREDUCE ( barea, mbod, SWREAL, SWSUM )
       call SWREDUCE ( bvol , mbod, SWREAL, SWSUM )
       call SWREDUCE ( bmoax, mbod, SWREAL, SWSUM )
       call SWREDUCE ( bmoay, mbod, SWREAL, SWSUM )
       call SWREDUCE ( bcob , mbod, SWREAL, SWSUM )
       if (STPNOW()) return
       !
       ! finalize computation of center of buoyancy
       bcob = bcob / bvol
       !
    endif
    !
    ! include pretension to counteract the buoyant force, if appropriate
    !
    if ( ifloat == 2 ) then
       !
       if ( mxmli > 0 ) then
          !
          do m = 1, mbod
             !
             if ( cpto(m) ) then
                !
                ! determine initial inclination of mooring line relative to vertical axis
                ! so that vertical component of pretension balances the buoyant force
                !
                fac = ( rhow * bvol(m) - bmass(m) ) * grav / float(nmlb(m))
                !
                do n = 1, nmlb(m)
                   !
                   gamma = atan2( bmli(m,n,8) - bmli(m,n,5), sqrt( (bmli(m,n,6)-bmli(m,n,3))**2 + (bmli(m,n,7)-bmli(m,n,4))**2 ) )
                   !
                   if ( gamma /= 0. ) then
                      bmli(m,n,10) = fac / sin(gamma)
                   else
                      ! horizintal line? no pretension
                      bmli(m,n,10) = 0.
                   endif
                   !
                enddo
                !
             else
                !
                bmli(m,:,10) = 0.
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
end subroutine SwashFloatObjects
