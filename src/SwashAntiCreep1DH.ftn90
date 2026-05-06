subroutine SwashAntiCreep1DH ( amatc, rhsc, rpo, kgrpnt, psm )
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
!    1.00: Koen Hilgersom
!
!   Updates
!
!    1.00, January 2015: New subroutine
!
!   Purpose
!
!   Adds anti-creep terms to the 1DH transport equation
!
!   Method
!
!   In case of steep bottom slopes and strong stable stratification, the
!   horizontal diffusion term must be transformed, generating some
!   curvature terms. When these so-called anti-creep terms are neglected,
!   this may lead to artificial mixing of constituent in vertical direction.
!
!   Details on the derivation and discretization of the anti-creep terms
!   can be found in
!
!   Technical documentation of WAQUA / TRIWAQ
!   SIMONA report nr. 99-01
!   M. Zijlema
!   pg. 31 - 33 (derivation) and pg. 90 (discretization)
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
    use m_genarr, only: wrk
    use m_parall
    use SwashFlowdata, only: mf, mfu, ml, mlu, flux, hkso, hs, vnu2d, zkso
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc)      , intent(in)    :: kgrpnt ! index table containing the address of each (active) grid point
                                                           ! =1: not active grid point
                                                           ! >1: active grid point
    !
    real, dimension(mcgrd,kmax,3), intent(inout) :: amatc  ! tri-diagonal matrix containing vertical terms in transport equation
    real                         , intent(in)    :: psm    ! Prandtl-Schmidt number
    real, dimension(mcgrd,kmax  ), intent(in)    :: rpo    ! constituent in wl-point at previous time level
    real, dimension(mcgrd,kmax  ), intent(inout) :: rhsc   ! right-hand side of tri-diagonal system of transport equation
!
!   Local variables
!
    integer, save  :: ient = 0 ! number of entries in this subroutine
    integer        :: k        ! loop counter over vertical layers
    integer        :: kd       ! index of layer k-1
    integer        :: ku       ! index of layer k+1
    integer        :: m        ! loop counter
    integer        :: md       ! index of point m-1
    integer        :: mend     ! end index of loop over wl-points
    integer        :: msta     ! start index of loop over wl-points
    integer        :: mu       ! index of point m+1
    integer        :: nm       ! pointer to m
    integer        :: nmd      ! pointer to m-1
    integer        :: nmf      ! pointer to mf
    integer        :: nml      ! pointer to ml
    integer        :: nmu      ! pointer to m+1
    !
    real           :: ctrm1    ! first part of contribution of anti-creep terms in x-direction
    real           :: ctrm2    ! second part of contribution of anti-creep terms in x-direction
    real           :: ctrm3    ! third part of contribution of anti-creep terms in x-direction
    real           :: ctrm4    ! fourth part of contribution of anti-creep terms in x-direction
    real           :: dif2d    ! effective horizontal eddy diffusivity coefficient
    real           :: rdx      ! reciprocal of mesh size
    real           :: stabmx   ! auxiliary variable with maximum diffusivity based stability criterion
    !
    logical        :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashAntiCreep1DH')
    !
    nmf = kgrpnt(mf)
    nml = kgrpnt(ml)
    !
    rdx = 1./dx
    !
    stabmx = 0.5 * dx * dx / dt
    !
    flux(nmf,:,1) = 0.
    flux(nml,:,1) = 0.
    !
    ! compute d( D h(k) dz(k)/dx d c(k-1/2)/dz ) / dx at internal cell-faces
    ! (Eq. 4.83)
    !
    if ( LMXF ) then
       msta = mf + 1   ! first internal cell-face
    else
       msta = mf       ! first internal cell-face at subdomain interface
    endif
    mend = ml - 1      ! last  internal cell-face
    !
    do m = msta, mend
       !
       mu = m + 1
       !
       nm  = kgrpnt(m )
       nmu = kgrpnt(mu)
       !
       ! compute effective horizontal diffusivity coefficient in u-point
       !
       if ( hdiff > 0. ) then
          !
          dif2d = hdiff
          !
       else
          !
          if ( ihvisc == 2 .or. ihvisc == 3 ) then
             !
             dif2d = 0.5 * ( vnu2d(nm) + vnu2d(nmu) ) / psm
             !
          else
             !
             dif2d = 0.
             !
          endif
          !
       endif
       !
       ! clip diffusivity coefficient, if necessary
       !
       if ( .not. dif2d < stabmx ) dif2d = stabmx
       !
       do k = 1, kmax
          !
          kd = max(   1,k-1)
          ku = min(kmax,k+1)
          !
          flux(nm,k,1) = 0.5 * ( hkso(nm,k) + hkso(nmu,k) ) * dif2d * rdx *                    &

!         --- compute d z(k) / dx

                         0.5 * ( zkso(nmu,k-1) - zkso(nm,k-1) + zkso(nmu,k) - zkso(nm,k) ) *   &

!         --- compute d c(k-1/2) / dz

                         ( rpo(nmu,kd) + rpo(nm,kd) - rpo(nmu,ku) - rpo(nm,ku) ) * 2. /        &
                         ( (ku-k) * (zkso(nmu,k -1) - zkso(nmu,ku)) +                          &
                           (k-kd) * (zkso(nmu,kd-1) - zkso(nmu,k )) +                          &
                           (ku-k) * (zkso(nm ,k -1) - zkso(nm ,ku)) +                          &
                           (k-kd) * (zkso(nm ,kd-1) - zkso(nm ,k )) )
          !
       enddo
       !
    enddo
    !
    ! update right-hand side with flux parts of anti-creep terms
    !
    ! not computed for end point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    do k = 1, kmax
       !
       do m = mfu, mend
          !
          md = m - 1
          !
          nm  = kgrpnt(m )
          nmd = kgrpnt(md)
          !
          if ( hs(nm) > epsdry ) then
             !
             rhsc(nm,k) = rhsc(nm,k) - rdx * ( flux(nm,k,1) - flux(nmd,k,1) )
             !
          endif
          !
       enddo
       !
    enddo
    !
    ! compute constituent at interface between layers k and k+1 (stored in wrk)
    ! (Eq. 4.87)
    !
    do m = mf, mlu
       !
       nm = kgrpnt(m)
       !
       do k = 1, kmax-1
          wrk(nm,k) = ( hkso(nm,k) * rpo(nm,k+1) + hkso(nm,k+1) * rpo(nm,k) ) / ( hkso(nm,k) + hkso(nm,k+1) )
       enddo
       !
    enddo
    !
    ! update right-hand side with remaining (B and C) parts of anti-creep terms
    ! (Eq. 4.82)
    !
    do m = mfu, mend
       !
       md = m - 1
       mu = m + 1
       !
       nm   = kgrpnt(m )
       nmd  = kgrpnt(md)
       nmu  = kgrpnt(mu)
       !
       ! compute effective horizontal diffusivity coefficient in wl-point
       !
       if ( hdiff > 0. ) then
          !
          dif2d = hdiff
          !
       else
          !
          if ( ihvisc == 2 .or. ihvisc == 3 ) then
             !
             dif2d = vnu2d(nm) / psm
             !
          else
             !
             dif2d = 0.
             !
          endif
          !
       endif
       !
       ! clip diffusivity coefficient, if necessary
       !
       if ( .not. dif2d < stabmx ) dif2d = stabmx
       !
       do k = 1, kmax-1
          !
          ! compute d z(k) / dx at interface between layers k and k+1
          !
          ctrm1 = 0.5 * rdx * ( zkso(nmu,k) - zkso(nmd,k) )
          !
          ! compute dc / dx  at interface between layers k and k+1
          !
          ctrm2 = 0.5 * rdx * ( wrk(nmu,k) - wrk(nmd,k) )
          !
          ! compute second part of anti-creep in x-direction at interface between layers k and k+1
          ! (Eq. 4.85)
          !
          ctrm3 = dif2d * ctrm1 * ctrm2
          !
          ! compute third part of anti-creep in x-direction at interface between layers k and k+1
          ! (Eq. 4.88)
          !
          ctrm4 = 2. * dif2d * ctrm1*ctrm1 / ( hkso(nm,k) + hkso(nm,k+1) )
          !
          if ( hs(nm) > epsdry ) then
             !
             ! second part is taken explicitly
             !
             rhsc(nm,k  ) = rhsc(nm,k  ) + ctrm3            ! at interface between layers k and k+1
             rhsc(nm,k+1) = rhsc(nm,k+1) - ctrm3            ! at interface between layers k-1 and k
             !
             ! third part is taken implicitly
             !
             amatc(nm,k  ,1) = amatc(nm,k  ,1) + ctrm4      ! at interface between layers k and k+1
             amatc(nm,k  ,3) = amatc(nm,k  ,3) - ctrm4
             amatc(nm,k+1,1) = amatc(nm,k+1,1) + ctrm4      ! at interface between layers k-1 and k
             amatc(nm,k+1,2) = amatc(nm,k+1,2) - ctrm4
             !
          endif
          !
       enddo
       !
    enddo
    !
end subroutine SwashAntiCreep1DH
