subroutine SwashHydroLoads ( nflob, forx, fory, forz, momx, momy, momz )
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
!    1.00, February 2018: New subroutine
!    9.01,  October 2022: revised
!
!   Purpose
!
!   Calculates hydrodynamic loads acting on floating body
!
!   Method
!
!   Integration of the total pressure over the wet surface of the body
!   to obtain the hydrodynamic (Froude-Krylov) forces
!
!   The moments of total pressure around the center of mass acting on the
!   body surface is integrated in the same manner
!
!   For the integration we assume a linear distribution of non-hydrostatic
!   pressure over each vertical layer
!
!   The contribution to the integration consists of two parts: a rectangle
!   and a triangle, of which the center of mass is 1/2 and 2/3 of the height,
!   respectively, and is thus the point of action
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use outp_data, only: alpobj
    use m_genarr, only: kgrpnt, guu, gvv, xcgrid, ycgrid
    use m_parall
    use SwashFlowdata
    use SwashRigBoddata
!
    implicit none
!
!   Argument variables
!
    integer,                intent(in ) :: nflob ! number of floating objects
    !
    real, dimension(nflob), intent(out) :: forx  ! hydrodynamic force in x-direction
    real, dimension(nflob), intent(out) :: fory  ! hydrodynamic force in y-direction
    real, dimension(nflob), intent(out) :: forz  ! hydrodynamic force in z-direction
    real, dimension(nflob), intent(out) :: momx  ! hydrodynamic moment in x-direction
    real, dimension(nflob), intent(out) :: momy  ! hydrodynamic moment in y-direction
    real, dimension(nflob), intent(out) :: momz  ! hydrodynamic moment in z-direction
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter of vertical layers
    integer       :: l        ! label of body
    integer       :: m        ! loop counter
    integer       :: mend     ! end index of loop over u-points
    integer       :: md       ! index of point m-1
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nend     ! end index of loop over v-points
    integer       :: nm       ! pointer to m,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
    !
    real          :: beta     ! beta = 180 deg - alpobj
    real          :: calpo    ! cosine of alpobj
    real          :: cbeta    ! cosine of beta
    real          :: dxl      ! local mesh size in x-direction
    real          :: dyl      ! local mesh size in y-direction
    real          :: dz       ! height of part of hull
    real          :: frc1     ! force contribution - part 1 (rectangle)
    real          :: frc2     ! force contribution - part 2 (triangle)
    real          :: frc1_rx  ! x-component of force contribution relative to rotated body - part 1 (rectangle)
    real          :: frc2_rx  ! x-component of force contribution relative to rotated body - part 2 (triangle)
    real          :: frc1_ry  ! y-component of force contribution relative to rotated body - part 1 (rectangle)
    real          :: frc2_ry  ! y-component of force contribution relative to rotated body - part 2 (triangle)
    real          :: qd       ! =q(:,k-1) if layer k-1 exists otherwise 0.
    real          :: qh       ! non-hydrostatic pressure at hull
    real          :: rx       ! x-component of the moment arm vector
    real          :: rxr      ! x-component of the moment arm vector relative to coordinate system of floating object
    real          :: ry       ! y-component of the moment arm vector
    real          :: ryr      ! y-component of the moment arm vector relative to coordinate system of floating object
    real          :: rz1      ! z-component of the moment arm vector - part 1 (rectangle)
    real          :: rz2      ! z-component of the moment arm vector - part 2 (triangle)
    real          :: salpo    ! sine of alpobj
    real          :: sbeta    ! sine of beta
    real          :: sloc     ! local layer interface k-1
    real          :: zloc     ! local layer interface k
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashHydroLoads')
    !
    ! initialize forces and moments
    !
    forx = 0.
    fory = 0.
    forz = 0.
    !
    momx = 0.
    momy = 0.
    momz = 0.
    !
    ! compute hydrodynamic loads acting on (moving) rigid body
    !
    if ( oned ) then
       !
       ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       do m = mfu, mend
          !
          mu = m + 1
          md = m - 1
          !
          nm  = kgrpnt(m ,1)
          nmu = kgrpnt(mu,1)
          nmd = kgrpnt(md,1)
          !
          ! vertical force and contribution to moment around y-axis
          !
  10      if ( presp(nm) == 1 ) then
             !
             l = lfbs(nm)
             if ( l == 0 .or. ( .not.bdof(l,3,1) .and. .not.bdof(l,2,2) ) ) goto 20
             !
             ! hydrostatic force
             !
             frc1 = rhow * grav * ( s1(nm) + flos(nm) ) * dx
             !
             ! point of action
             !
             rx = 0.5 * ( xcgrid(m,1) + xcgrid(md,1) ) - bcog(l,1)
             !
             if ( bdof(l,3,1) ) forz(l) = forz(l) +      frc1
             if ( bdof(l,2,2) ) momy(l) = momy(l) - rx * frc1
             !
             ! non-hydrostatic force
             !
             if ( ihydro /= 0 ) then
                !
                ! note: q is cell-centered
                frc1 = rhow * q(nm,1) * dx
                !
                if ( bdof(l,3,1) ) forz(l) = forz(l) +      frc1
                if ( bdof(l,2,2) ) momy(l) = momy(l) - rx * frc1
                !
             endif
             !
          endif
          !
          ! horizontal force and contribution to moment around y-axis
          !
  20      if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
             !
             l = lfbs(nmu)
             if ( l == 0 .or. ( .not.bdof(l,1,1) .and. .not.bdof(l,2,2) ) ) goto 30
             !
             dz = s1(nm) + flos(nmu)
             !
             ! hydrostatic force on the port side adjacent to body
             !
             frc2 = + rhow * grav * 0.5 * dz*dz
             !
             ! point of action
             !
             rz2 = s1(nm) - 2./3.*dz - bcog(l,3)
             !
             if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc2
             if ( bdof(l,2,2) ) momy(l) = momy(l) + rz2 * frc2
             !
             ! non-hydrostatic force on the port side adjacent to body
             !
             if ( ihydro /= 0 ) then
                !
                ! for each layer
                !
                do k = 1, kmax
                   !
                   ! note: q is defined at layer interface
                   if ( k == 1 ) then
                      qd = 0.
                   else
                      qd = q(nm,k-1)
                   endif
                   !
                   ! total depth or layer
                   !
                   if ( kmax == 1 ) then
                      zloc = -dps(nm)
                      sloc =   s1(nm)
                   else
                      zloc = zks(nm,k  )
                      sloc = zks(nm,k-1)
                   endif
                   !
                   if ( zloc > -flos(nmu) ) then
                      !
                      dz   = sloc - zloc
                      !
                      frc1 = + rhow * qd * dz
                      frc2 = + rhow * 0.5*( qd + q(nm,k) ) * dz
                      !
                      rz1  = sloc - 1./2.*dz - bcog(l,3)
                      rz2  = sloc - 2./3.*dz - bcog(l,3)
                      !
                      if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1 +       frc2
                      if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1 + rz2 * frc2
                      !
                   else
                      !
                      ! remaining part of the hull
                      !
                      dz = sloc + flos(nmu)
                      !
                      ! compute non-hydrostatic pressure at hull
                      qh = q(nm,k) - ( qd - q(nm,k) ) * ( zloc + flos(nmu) ) / ( sloc - zloc )
                      !
                      frc1 = + rhow * qd * dz
                      frc2 = + rhow * 0.5*( qd + qh ) * dz
                      !
                      rz1  = sloc - 1./2.*dz - bcog(l,3)
                      rz2  = sloc - 2./3.*dz - bcog(l,3)
                      !
                      if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1 +       frc2
                      if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1 + rz2 * frc2
                      !
                      exit
                      !
                   endif
                   !
                enddo
                !
             endif
             !
          endif
          !
  30      if ( presp(nm) == 0 .and. presp(nmd) == 1 ) then
             !
             l = lfbs(nmd)
             if ( l == 0 .or. ( .not.bdof(l,1,1) .and. .not.bdof(l,2,2) ) ) goto 40
             !
             dz = s1(nm) + flos(nmd)
             !
             ! hydrostatic force on the starboard side adjacent to body
             !
             frc2 = - rhow * grav * 0.5 * dz*dz
             !
             ! point of action
             !
             rz2 = s1(nm) - 2./3.*dz - bcog(l,3)
             !
             if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc2
             if ( bdof(l,2,2) ) momy(l) = momy(l) + rz2 * frc2
             !
             ! non-hydrostatic force on the starboard side adjacent to body
             !
             if ( ihydro /= 0 ) then
                !
                ! for each layer
                !
                do k = 1, kmax
                   !
                   ! note: q is defined at layer interface
                   if ( k == 1 ) then
                      qd = 0.
                   else
                      qd = q(nm,k-1)
                   endif
                   !
                   ! total depth or layer
                   !
                   if ( kmax == 1 ) then
                      zloc = -dps(nm)
                      sloc =   s1(nm)
                   else
                      zloc = zks(nm,k  )
                      sloc = zks(nm,k-1)
                   endif
                   !
                   if ( zloc > -flos(nmd) ) then
                      !
                      dz   = sloc - zloc
                      !
                      frc1 = - rhow * qd * dz
                      frc2 = - rhow * 0.5*( qd + q(nm,k) ) * dz
                      !
                      rz1  = sloc - 1./2.*dz - bcog(l,3)
                      rz2  = sloc - 2./3.*dz - bcog(l,3)
                      !
                      if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1 +       frc2
                      if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1 + rz2 * frc2
                      !
                   else
                      !
                      ! remaining part of the hull
                      !
                      dz = sloc + flos(nmd)
                      !
                      ! compute non-hydrostatic pressure at hull
                      qh = q(nm,k) - ( qd - q(nm,k) ) * ( zloc + flos(nmd) ) / ( sloc - zloc )
                      !
                      frc1 = - rhow * qd * dz
                      frc2 = - rhow * 0.5*( qd + qh ) * dz
                      !
                      rz1  = sloc - 1./2.*dz - bcog(l,3)
                      rz2  = sloc - 2./3.*dz - bcog(l,3)
                      !
                      if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1 +       frc2
                      if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1 + rz2 * frc2
                      !
                      exit
                      !
                   endif
                   !
                enddo
                !
             endif
             !
          endif
          !
  40      if ( presp(nm) == 1 ) then
             !
             l = lfbs(nm)
             if ( l == 0 .or. ( .not.bdof(l,1,1) .and. .not.bdof(l,2,2) ) ) goto 50
             !
             ! underneath floating object
             !
             if ( flos(nm) < flos(nmu) ) then
                !
                ! hydrostatic force on the port side
                !
                dz   = flos(nmu) - flos(nm)
                !
                frc1 = + rhow * grav * ( s1(nm)+flos(nm) ) * dz
                frc2 = + rhow * grav * 0.5 * dz*dz
                !
                rz1  = -flos(nm) - 1./2.*dz - bcog(l,3)
                rz2  = -flos(nm) - 2./3.*dz - bcog(l,3)
                !
                if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1 +       frc2
                if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1 + rz2 * frc2
                !
                ! non-hydrostatic force on the port side
                !
                if ( ihydro /= 0 ) then
                   !
                   do k = 1, kmax
                      !
                      if ( kmax == 1 ) then
                         zloc = -dps (nm)
                         sloc = -flos(nm)
                      else
                         zloc = zks(nm,k  )
                         sloc = zks(nm,k-1)
                      endif
                      !
                      if ( zloc > -flos(nmu) ) then
                         !
                         dz   = sloc - zloc
                         !
                         ! note: q is cell-centered
                         frc1 = + rhow * q(nm,k) * dz
                         !
                         rz1  = sloc - 1./2.*dz - bcog(l,3)
                         !
                         if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1
                         if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1
                         !
                      else
                         !
                         dz   = sloc + flos(nmu)
                         !
                         frc1 = + rhow * q(nm,k) * dz
                         !
                         rz1  = sloc - 1./2.*dz - bcog(l,3)
                         !
                         if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1
                         if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1
                         !
                         exit
                         !
                      endif
                      !
                   enddo
                   !
                endif
                !
             endif
             !
             if ( flos(nm) < flos(nmd) ) then
                !
                ! hydrostatic force on the starboard side
                !
                dz   = flos(nmd) - flos(nm)
                !
                frc1 = - rhow * grav * ( s1(nm)+flos(nm) ) * dz
                frc2 = - rhow * grav * 0.5 * dz*dz
                !
                rz1  = -flos(nm) - 1./2.*dz - bcog(l,3)
                rz2  = -flos(nm) - 2./3.*dz - bcog(l,3)
                !
                if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1 +       frc2
                if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1 + rz2 * frc2
                !
                ! non-hydrostatic force on the starboard side
                !
                if ( ihydro /= 0 ) then
                   !
                   do k = 1, kmax
                      !
                      if ( kmax == 1 ) then
                         zloc = -dps (nm)
                         sloc = -flos(nm)
                      else
                         zloc = zks(nm,k  )
                         sloc = zks(nm,k-1)
                      endif
                      !
                      if ( zloc > -flos(nmd) ) then
                         !
                         dz   = sloc - zloc
                         !
                         ! note: q is cell-centered
                         frc1 = - rhow * q(nm,k) * dz
                         !
                         rz1  = sloc - 1./2.*dz - bcog(l,3)
                         !
                         if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1
                         if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1
                         !
                      else
                         !
                         dz   = sloc + flos(nmd)
                         !
                         frc1 = - rhow * q(nm,k) * dz
                         !
                         rz1  = sloc - 1./2.*dz - bcog(l,3)
                         !
                         if ( bdof(l,1,1) ) forx(l) = forx(l) +       frc1
                         if ( bdof(l,2,2) ) momy(l) = momy(l) + rz1 * frc1
                         !
                         exit
                         !
                      endif
                      !
                   enddo
                   !
                endif
                !
             endif
             !
          endif
          !
       enddo
       !
    else
       !
       ! rotate grid of floating body
       !
       if ( ifloat == 1 ) then
          !
          ! for output purposes
          !
          if ( alpobj < 0. ) then
             beta = -alpobj - PI
          else
             beta = PI - alpobj
          endif
          !
          calpo = cos(alpobj)
          salpo = sin(alpobj)
          cbeta = cos(beta)
          sbeta = sin(beta)
          !
       else if ( ifloat == 2 ) then
          !
          ! in case of moving bodies, do not rotate!
          !
          calpo =  1.
          salpo =  0.
          cbeta = -1.
          sbeta =  0.
          !
       endif
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
             mu = m + 1
             md = m - 1
             nu = n + 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmu = kgrpnt(mu,n )
             nmd = kgrpnt(md,n )
             num = kgrpnt(m ,nu)
             ndm = kgrpnt(m ,nd)
             !
             if ( nmd == 1 ) nmd = nm
             if ( nmu == 1 ) nmu = nm
             if ( ndm == 1 ) ndm = nm
             if ( num == 1 ) num = nm
             !
             dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
             dyl = 0.5 * ( guu(nm) + guu(nmd) )
             !
             ! vertical force and contribution to moments around x- and y-axis
             !
  50         if ( presp(nm) == 1 ) then
                !
                l = lfbs(nm)
                if ( l == 0 .or. ( .not.bdof(l,3,1) .and. .not.bdof(l,1,2) .and. .not.bdof(l,2,2) ) ) goto 60
                !
                ! hydrostatic force
                !
                frc1 = rhow * grav * ( s1(nm) + flos(nm) ) * dxl * dyl
                !
                ! points of action
                !
                rx = 0.5 * ( xcgrid(m,n) + xcgrid(md,n ) ) - bcog(l,1)
                ry = 0.5 * ( ycgrid(m,n) + ycgrid(m ,nd) ) - bcog(l,2)
                !
                ! rotate these points to coordinate system of floating object
                rxr =  rx*calpo + ry*salpo
                ryr = -rx*salpo + ry*calpo
                !
                if ( bdof(l,3,1) ) forz(l) = forz(l) +       frc1
                if ( bdof(l,1,2) ) momx(l) = momx(l) + ryr * frc1
                if ( bdof(l,2,2) ) momy(l) = momy(l) - rxr * frc1
                !
                ! non-hydrostatic force
                !
                if ( ihydro /= 0 ) then
                   !
                   frc1 = rhow * q(nm,1) * dxl * dyl
                   !
                   if ( bdof(l,3,1) ) forz(l) = forz(l) +       frc1
                   if ( bdof(l,1,2) ) momx(l) = momx(l) + ryr * frc1
                   if ( bdof(l,2,2) ) momy(l) = momy(l) - rxr * frc1
                   !
                endif
                !
             endif
             !
             ! horizontal force in x-direction and contribution to moments around y- and z-axis
             !
  60         if ( presp(nm) == 0 .and. presp(nmu) == 1 ) then
                !
                l = lfbs(nmu)
                if ( l == 0 .or. ( .not.bdof(l,1,1) .and. .not.bdof(l,2,2) .and. .not.bdof(l,3,2) ) ) goto 70
                !
                dz = s1(nm) + flos(nmu)
                !
                ! hydrostatic force on the port side adjacent to body
                !
                frc2 = + rhow * grav * 0.5 * dz*dz * dyl
                !
                ! rotate forces to grid of floating object
                frc2_rx = frc2 * cbeta
                frc2_ry = frc2 * sbeta
                !
                ! points of action
                !
                rz2 = s1(nm) - 2./3.*dz - bcog(l,3)
                !
                ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                !
                if ( bdof(l,1,1) ) forx(l) = forx(l) - frc2_rx
                if ( bdof(l,1,1) ) fory(l) = fory(l) - frc2_ry
                !
                if ( bdof(l,2,2) ) momy(l) = momy(l) - rz2 * frc2_rx
                if ( bdof(l,2,2) ) momx(l) = momx(l) - rz2 * frc2_ry
                if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc2
                !
                ! non-hydrostatic force on the port side adjacent to body
                !
                if ( ihydro /= 0 ) then
                   !
                   ! for each layer
                   !
                   do k = 1, kmax
                      !
                      if ( k == 1 ) then
                         qd = 0.
                      else
                         qd = q(nm,k-1)
                      endif
                      !
                      ! total depth or layer
                      !
                      if ( kmax == 1 ) then
                         zloc = -dps(nm)
                         sloc =   s1(nm)
                      else
                         zloc = zks(nm,k  )
                         sloc = zks(nm,k-1)
                      endif
                      !
                      if ( zloc > -flos(nmu) ) then
                         !
                         dz   = sloc - zloc
                         !
                         frc1 = + rhow * qd * dz * dyl
                         frc2 = + rhow * 0.5*( qd + q(nm,k) ) * dz * dyl
                         !
                         ! rotate forces to grid of floating object
                         frc1_rx = frc1 * cbeta
                         frc1_ry = frc1 * sbeta
                         frc2_rx = frc2 * cbeta
                         frc2_ry = frc2 * sbeta
                         !
                         rz1 = sloc - 1./2.*dz - bcog(l,3)
                         rz2 = sloc - 2./3.*dz - bcog(l,3)
                         !
                         ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                         !
                         if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx - frc2_rx
                         if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                         !
                         if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx - rz2 * frc2_rx
                         if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry - rz2 * frc2_ry
                         if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1    -  ry * frc2
                         !
                      else
                         !
                         ! remaining part of the hull
                         !
                         dz = sloc + flos(nmu)
                         !
                         ! compute non-hydrostatic pressure at hull
                         qh = q(nm,k) - ( qd - q(nm,k) ) * ( zloc + flos(nmu) ) / ( sloc - zloc )
                         !
                         frc1 = + rhow * qd * dz * dyl
                         frc2 = + rhow * 0.5*( qd + qh ) * dz * dyl
                         !
                         ! rotate forces to grid of floating object
                         frc1_rx = frc1 * cbeta
                         frc1_ry = frc1 * sbeta
                         frc2_rx = frc2 * cbeta
                         frc2_ry = frc2 * sbeta
                         !
                         rz1 = sloc - 1./2.*dz - bcog(l,3)
                         rz2 = sloc - 2./3.*dz - bcog(l,3)
                         !
                         ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                         !
                         if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx - frc2_rx
                         if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                         !
                         if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx - rz2 * frc2_rx
                         if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry - rz2 * frc2_ry
                         if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1    -  ry * frc2
                         !
                         exit
                         !
                      endif
                      !
                   enddo
                   !
                endif
                !
             endif
             !
  70         if ( presp(nm) == 0 .and. presp(nmd) == 1 ) then
                !
                l = lfbs(nmd)
                if ( l == 0 .or. ( .not.bdof(l,1,1) .and. .not.bdof(l,2,2) .and. .not.bdof(l,3,2) ) ) goto 80
                !
                dz = s1(nm) + flos(nmd)
                !
                ! hydrostatic force on the starboard side adjacent to body
                !
                frc2 = - rhow * grav * 0.5 * dz*dz * dyl
                !
                ! rotate forces to grid of floating object
                frc2_rx = frc2 * cbeta
                frc2_ry = frc2 * sbeta
                !
                ! points of action
                !
                rz2 = s1(nm) - 2./3.*dz - bcog(l,3)
                !
                ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                !
                if ( bdof(l,1,1) ) forx(l) = forx(l) - frc2_rx
                if ( bdof(l,1,1) ) fory(l) = fory(l) - frc2_ry
                !
                if ( bdof(l,2,2) ) momy(l) = momy(l) - rz2 * frc2_rx
                if ( bdof(l,2,2) ) momx(l) = momx(l) - rz2 * frc2_ry
                if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc2
                !
                ! non-hydrostatic force on the starboard side adjacent to body
                !
                if ( ihydro /= 0 ) then
                   !
                   ! for each layer
                   !
                   do k = 1, kmax
                      !
                      if ( k == 1 ) then
                         qd = 0.
                      else
                         qd = q(nm,k-1)
                      endif
                      !
                      ! total depth or layer
                      !
                      if ( kmax == 1 ) then
                         zloc = -dps(nm)
                         sloc =   s1(nm)
                      else
                         zloc = zks(nm,k  )
                         sloc = zks(nm,k-1)
                      endif
                      !
                      if ( zloc > -flos(nmd) ) then
                         !
                         dz   = sloc - zloc
                         !
                         frc1 = - rhow * qd * dz * dyl
                         frc2 = - rhow * 0.5*( qd + q(nm,k) ) * dz * dyl
                         !
                         ! rotate forces to grid of floating object
                         frc1_rx = frc1 * cbeta
                         frc1_ry = frc1 * sbeta
                         frc2_rx = frc2 * cbeta
                         frc2_ry = frc2 * sbeta
                         !
                         rz1 = sloc - 1./2.*dz - bcog(l,3)
                         rz2 = sloc - 2./3.*dz - bcog(l,3)
                         !
                         ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                         !
                         if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx - frc2_rx
                         if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                         !
                         if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx - rz2 * frc2_rx
                         if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry - rz2 * frc2_ry
                         if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1    -  ry * frc2
                         !
                      else
                         !
                         ! remaining part of the hull
                         !
                         dz = sloc + flos(nmd)
                         !
                         ! compute non-hydrostatic pressure at hull
                         qh = q(nm,k) - ( qd - q(nm,k) ) * ( zloc + flos(nmd) ) / ( sloc - zloc )
                         !
                         frc1 = - rhow * qd * dz * dyl
                         frc2 = - rhow * 0.5*( qd + qh ) * dz * dyl
                         !
                         ! rotate forces to grid of floating object
                         frc1_rx = frc1 * cbeta
                         frc1_ry = frc1 * sbeta
                         frc2_rx = frc2 * cbeta
                         frc2_ry = frc2 * sbeta
                         !
                         rz1 = sloc - 1./2.*dz - bcog(l,3)
                         rz2 = sloc - 2./3.*dz - bcog(l,3)
                         !
                         ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                         !
                         if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx - frc2_rx
                         if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                         !
                         if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx - rz2 * frc2_rx
                         if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry - rz2 * frc2_ry
                         if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1    -  ry * frc2
                         !
                         exit
                         !
                      endif
                      !
                   enddo
                   !
                endif
                !
             endif
             !
  80         if ( presp(nm) == 1 ) then
                !
                l = lfbs(nm)
                if ( l == 0 .or. ( .not.bdof(l,1,1) .and. .not.bdof(l,2,2) .and. .not.bdof(l,3,2) ) ) goto 90
                !
                ! underneath floating object
                !
                if ( flos(nm) < flos(nmu) ) then
                   !
                   ! hydrostatic force on the port side
                   !
                   dz   = flos(nmu) - flos(nm)
                   !
                   frc1 = + rhow * grav * ( s1(nm)+flos(nm) ) * dz * dyl
                   frc2 = + rhow * grav * 0.5 * dz*dz * dyl
                   !
                   ! rotate forces to grid of floating object
                   frc1_rx = frc1 * cbeta
                   frc1_ry = frc1 * sbeta
                   frc2_rx = frc2 * cbeta
                   frc2_ry = frc2 * sbeta
                   !
                   rz1 = -flos(nm) - 1./2.*dz - bcog(l,3)
                   rz2 = -flos(nm) - 2./3.*dz - bcog(l,3)
                   !
                   ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                   !
                   if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx - frc2_rx
                   if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                   !
                   if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx - rz2 * frc2_rx
                   if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry - rz2 * frc2_ry
                   if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1    -  ry * frc2
                   !
                   ! non-hydrostatic force on the port side
                   !
                   if ( ihydro /= 0 ) then
                      !
                      do k = 1, kmax
                         !
                         if ( kmax == 1 ) then
                            zloc = -dps (nm)
                            sloc = -flos(nm)
                         else
                            zloc = zks(nm,k  )
                            sloc = zks(nm,k-1)
                         endif
                         !
                         if ( zloc > -flos(nmu) ) then
                            !
                            dz   = sloc - zloc
                            !
                            frc1 = + rhow * q(nm,k) * dz * dyl
                            !
                            ! rotate force to grid of floating object
                            frc1_rx = frc1 * cbeta
                            frc1_ry = frc1 * sbeta
                            !
                            rz1 = sloc - 1./2.*dz - bcog(l,3)
                            !
                            ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                            !
                            if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx
                            if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry
                            !
                            if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx
                            if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry
                            if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1
                            !
                         else
                            !
                            dz   = sloc + flos(nmu)
                            !
                            frc1 = + rhow * q(nm,k) * dz * dyl
                            !
                            ! rotate force to grid of floating object
                            frc1_rx = frc1 * cbeta
                            frc1_ry = frc1 * sbeta
                            !
                            rz1 = sloc - 1./2.*dz - bcog(l,3)
                            !
                            ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                            !
                            if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx
                            if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry
                            !
                            if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx
                            if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry
                            if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1
                            !
                            exit
                            !
                         endif
                         !
                      enddo
                      !
                   endif
                   !
                endif
                !
                if ( flos(nm) < flos(nmd) ) then
                   !
                   ! hydrostatic force on the starboard side
                   !
                   dz   = flos(nmd) - flos(nm)
                   !
                   frc1 = - rhow * grav * ( s1(nm)+flos(nm) ) * dz * dyl
                   frc2 = - rhow * grav * 0.5 * dz*dz * dyl
                   !
                   ! rotate forces to grid of floating object
                   frc1_rx = frc1 * cbeta
                   frc1_ry = frc1 * sbeta
                   frc2_rx = frc2 * cbeta
                   frc2_ry = frc2 * sbeta
                   !
                   rz1 = -flos(nm) - 1./2.*dz - bcog(l,3)
                   rz2 = -flos(nm) - 2./3.*dz - bcog(l,3)
                   !
                   ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                   !
                   if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx - frc2_rx
                   if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                   !
                   if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx - rz2 * frc2_rx
                   if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry - rz2 * frc2_ry
                   if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1    -  ry * frc2
                   !
                   ! non-hydrostatic force on the starboard side
                   !
                   if ( ihydro /= 0 ) then
                      !
                      do k = 1, kmax
                         !
                         if ( kmax == 1 ) then
                            zloc = -dps (nm)
                            sloc = -flos(nm)
                         else
                            zloc = zks(nm,k  )
                            sloc = zks(nm,k-1)
                         endif
                         !
                         if ( zloc > -flos(nmd) ) then
                            !
                            dz   = sloc - zloc
                            !
                            frc1 = - rhow * q(nm,k) * dz * dyl
                            !
                            ! rotate force to grid of floating object
                            frc1_rx = frc1 * cbeta
                            frc1_ry = frc1 * sbeta
                            !
                            rz1 = sloc - 1./2.*dz - bcog(l,3)
                            !
                            ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                            !
                            if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx
                            if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry
                            !
                            if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx
                            if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry
                            if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1
                            !
                         else
                            !
                            dz   = sloc + flos(nmd)
                            !
                            frc1 = - rhow * q(nm,k) * dz * dyl
                            !
                            ! rotate force to grid of floating object
                            frc1_rx = frc1 * cbeta
                            frc1_ry = frc1 * sbeta
                            !
                            rz1 = sloc - 1./2.*dz - bcog(l,3)
                            !
                            ry  = 0.5 * ( ycgrid(m,n) + ycgrid(m,nd) ) - bcog(l,2)
                            !
                            if ( bdof(l,1,1) ) forx(l) = forx(l) - frc1_rx
                            if ( bdof(l,1,1) ) fory(l) = fory(l) - frc1_ry
                            !
                            if ( bdof(l,2,2) ) momy(l) = momy(l) - rz1 * frc1_rx
                            if ( bdof(l,2,2) ) momx(l) = momx(l) - rz1 * frc1_ry
                            if ( bdof(l,3,2) ) momz(l) = momz(l) -  ry * frc1
                            !
                            exit
                            !
                         endif
                         !
                      enddo
                      !
                   endif
                   !
                endif
                !
             endif
             !
             ! horizontal force in y-direction and contribution to moments around x- and z-axis
             !
  90         if ( presp(nm) == 0 .and. presp(num) == 1 ) then
                !
                l = lfbs(num)
                if ( l == 0 .or. ( .not.bdof(l,2,1) .and. .not.bdof(l,1,2) .and. .not.bdof(l,3,2) ) ) goto 100
                !
                dz = s1(nm) + flos(num)
                !
                ! hydrostatic force on the port side adjacent to body
                !
                frc2 = + rhow * grav * 0.5 * dz*dz * dxl
                !
                ! rotate forces to grid of floating object
                frc2_rx = frc2 * sbeta
                frc2_ry = frc2 * cbeta
                !
                ! points of action
                !
                rz2 = s1(nm) - 2./3.*dz - bcog(l,3)
                !
                rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                !
                if ( bdof(l,2,1) ) forx(l) = forx(l) + frc2_rx
                if ( bdof(l,2,1) ) fory(l) = fory(l) - frc2_ry
                !
                if ( bdof(l,1,2) ) momy(l) = momy(l) + rz2 * frc2_rx
                if ( bdof(l,1,2) ) momx(l) = momx(l) + rz2 * frc2_ry
                if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc2
                !
                ! non-hydrostatic force on the port side adjacent to body
                !
                if ( ihydro /= 0 ) then
                   !
                   ! for each layer
                   !
                   do k = 1, kmax
                      !
                      if ( k == 1 ) then
                         qd = 0.
                      else
                         qd = q(nm,k-1)
                      endif
                      !
                      ! total depth or layer
                      !
                      if ( kmax == 1 ) then
                         zloc = -dps(nm)
                         sloc =   s1(nm)
                      else
                         zloc = zks(nm,k  )
                         sloc = zks(nm,k-1)
                      endif
                      !
                      if ( zloc > -flos(num) ) then
                         !
                         dz   = sloc - zloc
                         !
                         frc1 = + rhow * qd * dz * dxl
                         frc2 = + rhow * 0.5*( qd + q(nm,k) ) * dz * dxl
                         !
                         ! rotate forces to grid of floating object
                         frc1_rx = frc1 * sbeta
                         frc1_ry = frc1 * cbeta
                         frc2_rx = frc2 * sbeta
                         frc2_ry = frc2 * cbeta
                         !
                         rz1 = sloc - 1./2.*dz - bcog(l,3)
                         rz2 = sloc - 2./3.*dz - bcog(l,3)
                         !
                         rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                         !
                         if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx + frc2_rx
                         if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                         !
                         if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx + rz2 * frc2_rx
                         if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry + rz2 * frc2_ry
                         if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1    +  rx * frc2
                         !
                      else
                         !
                         ! remaining part of the hull
                         !
                         dz = sloc + flos(num)
                         !
                         ! compute non-hydrostatic pressure at hull
                         qh = q(nm,k) - ( qd - q(nm,k) ) * ( zloc + flos(num) ) / ( sloc - zloc )
                         !
                         frc1 = + rhow * qd * dz * dxl
                         frc2 = + rhow * 0.5*( qd + qh ) * dz * dxl
                         !
                         ! rotate forces to grid of floating object
                         frc1_rx = frc1 * sbeta
                         frc1_ry = frc1 * cbeta
                         frc2_rx = frc2 * sbeta
                         frc2_ry = frc2 * cbeta
                         !
                         rz1 = sloc - 1./2.*dz - bcog(l,3)
                         rz2 = sloc - 2./3.*dz - bcog(l,3)
                         !
                         rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                         !
                         if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx + frc2_rx
                         if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                         !
                         if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx + rz2 * frc2_rx
                         if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry + rz2 * frc2_ry
                         if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1    +  rx * frc2
                         !
                         exit
                         !
                      endif
                      !
                   enddo
                   !
                endif
                !
             endif
             !
 100         if ( presp(nm) == 0 .and. presp(ndm) == 1 ) then
                !
                l = lfbs(ndm)
                if ( l == 0 .or. ( .not.bdof(l,2,1) .and. .not.bdof(l,1,2) .and. .not.bdof(l,3,2) ) ) goto 110
                !
                dz = s1(nm) + flos(ndm)
                !
                ! hydrostatic force on the starboard side adjacent to body
                !
                frc2 = - rhow * grav * 0.5 * dz*dz * dxl
                !
                ! rotate forces to grid of floating object
                frc2_rx = frc2 * sbeta
                frc2_ry = frc2 * cbeta
                !
                ! points of action
                !
                rz2 = s1(nm) - 2./3.*dz - bcog(l,3)
                !
                rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                !
                if ( bdof(l,2,1) ) forx(l) = forx(l) + frc2_rx
                if ( bdof(l,2,1) ) fory(l) = fory(l) - frc2_ry
                !
                if ( bdof(l,1,2) ) momy(l) = momy(l) + rz2 * frc2_rx
                if ( bdof(l,1,2) ) momx(l) = momx(l) + rz2 * frc2_ry
                if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc2
                !
                ! non-hydrostatic force on the starboard side adjacent to body
                !
                if ( ihydro /= 0 ) then
                   !
                   ! for each layer
                   !
                   do k = 1, kmax
                      !
                      if ( k == 1 ) then
                         qd = 0.
                      else
                         qd = q(nm,k-1)
                      endif
                      !
                      ! total depth or layer
                      !
                      if ( kmax == 1 ) then
                         zloc = -dps(nm)
                         sloc =   s1(nm)
                      else
                         zloc = zks(nm,k  )
                         sloc = zks(nm,k-1)
                      endif
                      !
                      if ( zloc > -flos(ndm) ) then
                         !
                         dz   = sloc - zloc
                         !
                         frc1 = - rhow * qd * dz * dxl
                         frc2 = - rhow * 0.5*( qd + q(nm,k) ) * dz * dxl
                         !
                         ! rotate forces to grid of floating object
                         frc1_rx = frc1 * sbeta
                         frc1_ry = frc1 * cbeta
                         frc2_rx = frc2 * sbeta
                         frc2_ry = frc2 * cbeta
                         !
                         rz1 = sloc - 1./2.*dz - bcog(l,3)
                         rz2 = sloc - 2./3.*dz - bcog(l,3)
                         !
                         rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                         !
                         if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx + frc2_rx
                         if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                         !
                         if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx + rz2 * frc2_rx
                         if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry + rz2 * frc2_ry
                         if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1    +  rx * frc2
                         !
                      else
                         !
                         ! remaining part of the hull
                         !
                         dz = sloc + flos(ndm)
                         !
                         ! compute non-hydrostatic pressure at hull
                         qh = q(nm,k) - ( qd - q(nm,k) ) * ( zloc + flos(ndm) ) / ( sloc - zloc )
                         !
                         frc1 = - rhow * qd * dz * dxl
                         frc2 = - rhow * 0.5*( qd + qh ) * dz * dxl
                         !
                         ! rotate forces to grid of floating object
                         frc1_rx = frc1 * sbeta
                         frc1_ry = frc1 * cbeta
                         frc2_rx = frc2 * sbeta
                         frc2_ry = frc2 * cbeta
                         !
                         rz1 = sloc - 1./2.*dz - bcog(l,3)
                         rz2 = sloc - 2./3.*dz - bcog(l,3)
                         !
                         rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                         !
                         if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx + frc2_rx
                         if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                         !
                         if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx + rz2 * frc2_rx
                         if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry + rz2 * frc2_ry
                         if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1    +  rx * frc2
                         !
                         exit
                         !
                      endif
                      !
                   enddo
                   !
                endif
                !
             endif
             !
 110         if ( presp(nm) == 1 ) then
                !
                l = lfbs(nm)
                if ( l == 0 .or. ( .not.bdof(l,2,1) .and. .not.bdof(l,1,2) .and. .not.bdof(l,3,2) ) ) goto 120
                !
                ! underneath floating object
                !
                if ( flos(nm) < flos(num) ) then
                   !
                   ! hydrostatic force on the port side
                   !
                   dz   = flos(num) - flos(nm)
                   !
                   frc1 = + rhow * grav * ( s1(nm)+flos(nm) ) * dz * dxl
                   frc2 = + rhow * grav * 0.5 * dz*dz * dxl
                   !
                   ! rotate forces to grid of floating object
                   frc1_rx = frc1 * sbeta
                   frc1_ry = frc1 * cbeta
                   frc2_rx = frc2 * sbeta
                   frc2_ry = frc2 * cbeta
                   !
                   rz1 = -flos(nm) - 1./2.*dz - bcog(l,3)
                   rz2 = -flos(nm) - 2./3.*dz - bcog(l,3)
                   !
                   rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                   !
                   if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx + frc2_rx
                   if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                   !
                   if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx + rz2 * frc2_rx
                   if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry + rz2 * frc2_ry
                   if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1    +  rx * frc2
                   !
                   ! non-hydrostatic force on the port side
                   !
                   if ( ihydro /= 0 ) then
                      !
                      do k = 1, kmax
                         !
                         if ( kmax == 1 ) then
                            zloc = -dps (nm)
                            sloc = -flos(nm)
                         else
                            zloc = zks(nm,k  )
                            sloc = zks(nm,k-1)
                         endif
                         !
                         if ( zloc > -flos(num) ) then
                            !
                            dz   = sloc - zloc
                            !
                            frc1 = + rhow * q(nm,k) * dz * dxl
                            !
                            ! rotate force to grid of floating object
                            frc1_rx = frc1 * sbeta
                            frc1_ry = frc1 * cbeta
                            !
                            rz1 = sloc - 1./2.*dz - bcog(l,3)
                            !
                            rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                            !
                            if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx
                            if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry
                            !
                            if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx
                            if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry
                            if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1
                            !
                         else
                            !
                            dz   = sloc + flos(num)
                            !
                            frc1 = + rhow * q(nm,k) * dz * dxl
                            !
                            ! rotate force to grid of floating object
                            frc1_rx = frc1 * sbeta
                            frc1_ry = frc1 * cbeta
                            !
                            rz1 = sloc - 1./2.*dz - bcog(l,3)
                            !
                            rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                            !
                            if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx
                            if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry
                            !
                            if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx
                            if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry
                            if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1
                            !
                            exit
                            !
                         endif
                         !
                      enddo
                      !
                   endif
                   !
                endif
                !
                if ( flos(nm) < flos(ndm) ) then
                   !
                   ! hydrostatic force on the starboard side
                   !
                   dz   = flos(ndm) - flos(nm)
                   !
                   frc1 = - rhow * grav * ( s1(nm)+flos(nm) ) * dz * dxl
                   frc2 = - rhow * grav * 0.5 * dz*dz * dxl
                   !
                   ! rotate forces to grid of floating object
                   frc1_rx = frc1 * sbeta
                   frc1_ry = frc1 * cbeta
                   frc2_rx = frc2 * sbeta
                   frc2_ry = frc2 * cbeta
                   !
                   rz1 = -flos(nm) - 1./2.*dz - bcog(l,3)
                   rz2 = -flos(nm) - 2./3.*dz - bcog(l,3)
                   !
                   rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                   !
                   if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx + frc2_rx
                   if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry - frc2_ry
                   !
                   if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx + rz2 * frc2_rx
                   if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry + rz2 * frc2_ry
                   if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1    +  rx * frc2
                   !
                   ! non-hydrostatic force on the starboard side
                   !
                   if ( ihydro /= 0 ) then
                      !
                      do k = 1, kmax
                         !
                         if ( kmax == 1 ) then
                            zloc = -dps (nm)
                            sloc = -flos(nm)
                         else
                            zloc = zks(nm,k  )
                            sloc = zks(nm,k-1)
                         endif
                         !
                         if ( zloc > -flos(ndm) ) then
                            !
                            dz   = sloc - zloc
                            !
                            frc1 = - rhow * q(nm,k) * dz * dxl
                            !
                            ! rotate force to grid of floating object
                            frc1_rx = frc1 * sbeta
                            frc1_ry = frc1 * cbeta
                            !
                            rz1 = sloc - 1./2.*dz - bcog(l,3)
                            !
                            rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                            !
                            if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx
                            if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry
                            !
                            if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx
                            if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry
                            if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1
                            !
                         else
                            !
                            dz   = sloc + flos(ndm)
                            !
                            frc1 = - rhow * q(nm,k) * dz * dxl
                            !
                            ! rotate force to grid of floating object
                            frc1_rx = frc1 * sbeta
                            frc1_ry = frc1 * cbeta
                            !
                            rz1 = sloc - 1./2.*dz - bcog(l,3)
                            !
                            rx  = 0.5 * ( xcgrid(m,n) + xcgrid(md,n) ) - bcog(l,1)
                            !
                            if ( bdof(l,2,1) ) forx(l) = forx(l) + frc1_rx
                            if ( bdof(l,2,1) ) fory(l) = fory(l) - frc1_ry
                            !
                            if ( bdof(l,1,2) ) momy(l) = momy(l) + rz1 * frc1_rx
                            if ( bdof(l,1,2) ) momx(l) = momx(l) + rz1 * frc1_ry
                            if ( bdof(l,3,2) ) momz(l) = momz(l) +  rx * frc1
                            !
                            exit
                            !
                         endif
                         !
                      enddo
                      !
                   endif
                   !
                endif
                !
             endif
             !
 120         continue
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! accumulate hydrodynamic loads over all subdomains
    !
    call SWREDUCE ( forx , nflob, SWREAL, SWSUM )
    call SWREDUCE ( forz , nflob, SWREAL, SWSUM )
    call SWREDUCE ( momy , nflob, SWREAL, SWSUM )
    !
    if ( .not.oned ) then
       !
       call SWREDUCE ( fory , nflob, SWREAL, SWSUM )
       call SWREDUCE ( momx , nflob, SWREAL, SWSUM )
       call SWREDUCE ( momz , nflob, SWREAL, SWSUM )
       !
    endif
    !
end subroutine SwashHydroLoads
