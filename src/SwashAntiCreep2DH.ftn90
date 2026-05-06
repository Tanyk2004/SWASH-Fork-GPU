subroutine SwashAntiCreep2DH ( amatc, rhsc, rpo, psm )
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
!   Adds anti-creep terms to the 2DH transport equation
!
!   Method
!
!   In case of steep bottom slopes and strong stable stratification, the
!   horizontal diffusion terms must be transformed, generating some
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
    use m_genarr, only: kgrpnt, guu, guv, gvu, gvv, gsqs, wrk
    use m_parall
    use SwashFlowdata, only: mf, mfu, ml, mlu, nf, nfu, nl, nlu, flux, hkso, hs, vnu2d, zkso
!
    implicit none
!
!   Argument variables
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
    integer        :: n        ! loop counter
    integer        :: nd       ! index of point n-1
    integer        :: ndm      ! pointer to m,n-1
    integer        :: ndmd     ! pointer to m-1,n-1
    integer        :: nend     ! end index of loop over wl-points in y-direction
    integer        :: nfm      ! pointer to m,nf
    integer        :: nlm      ! pointer to m,nl
    integer        :: nm       ! pointer to m,n
    integer        :: nmd      ! pointer to m-1,n
    integer        :: nmf      ! pointer to mf,n
    integer        :: nml      ! pointer to ml,n
    integer        :: nmu      ! pointer to m+1,n
    integer        :: nsta     ! start index of loop over wl-points in y-direction
    integer        :: nu       ! index of point n+1
    integer        :: num      ! pointer to m,n+1
    !
    real           :: ctrm1    ! first part of contribution of anti-creep terms in x-direction
    real           :: ctrm2    ! second part of contribution of anti-creep terms in x-direction
    real           :: ctrm3    ! third part of contribution of anti-creep terms in x-direction
    real           :: ctrm4    ! fourth part of contribution of anti-creep terms in x-direction
    real           :: ctrm5    ! fifth part of contribution of anti-creep terms in x-direction
    real           :: ctrn1    ! first part of contribution of anti-creep terms in y-direction
    real           :: ctrn2    ! second part of contribution of anti-creep terms in y-direction
    real           :: ctrn3    ! third part of contribution of anti-creep terms in y-direction
    real           :: ctrn4    ! fourth part of contribution of anti-creep terms in y-direction
    real           :: ctrn5    ! fifth part of contribution of anti-creep terms in y-direction
    real           :: dif2d    ! effective horizontal eddy diffusivity coefficient
    real           :: dxl      ! local mesh size in x-direction
    real           :: dyl      ! local mesh size in y-direction
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
    if (ltrace) call strace (ient,'SwashAntiCreep2DH')
    !
    ! compute d( D h(k) dz(k)/dx d c(k-1/2)/dz ) / dx at internal cell-faces
    ! (Eq. 4.83)
    !
    if ( LMXF ) then
       msta = mf + 1         ! first internal u-point
    else
       msta = mf             ! first internal u-point at subdomain interface
    endif
    if ( lreptx ) then
       mend = ml             ! last  internal u-point in case of repeating grid
    else
       mend = ml - 1         ! last  internal u-point
    endif
    !
    do n = nfu, nl
       !
       nmf = kgrpnt(mf,n)
       nml = kgrpnt(ml,n)
       !
       flux(nmf,:,1) = 0.
       flux(nml,:,1) = 0.
       !
       do m = msta, mend
          !
          mu = m + 1
          nd = n - 1
          !
          nm  = kgrpnt(m ,n )
          nmu = kgrpnt(mu,n )
          ndm = kgrpnt(m ,nd)
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
                dif2d = 0.5 * ( vnu2d(nm) + vnu2d(ndm) ) / psm
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
          dxl = gvu(nm)
          dyl = guu(nm)
          !
          stabmx = 1./(2.*dt*(1./(dxl*dxl) + 1./(dyl*dyl)))
          !
          if ( .not. dif2d < stabmx ) dif2d = stabmx
          !
          do k = 1, kmax
             !
             kd = max(   1,k-1)
             ku = min(kmax,k+1)
             !
             flux(nm,k,1) = 0.5 * ( hkso(nm,k) + hkso(nmu,k) ) * dif2d * guu(nm) / gvu(nm)    *   &

!            --- compute d z(k) / dx

                            0.5 * ( zkso(nmu,k-1) - zkso(nm,k-1) + zkso(nmu,k) - zkso(nm,k) ) *   &

!            --- compute d c(k-1/2) / dz

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
    enddo
    !
    ! synchronize fluxes at appropriate boundaries in case of repeating grid
    !
    call periodic ( flux(1,1,1), kgrpnt, 1, kmax )
    !
    ! compute d( D h(k) dz(k)/dy d c(k-1/2)/dz ) / dy at internal cell-faces
    ! (Eq. 4.84)
    !
    if ( LMYF ) then
       nsta = nf + 1         ! first internal v-point
    else
       nsta = nf             ! first internal v-point at subdomain interface
    endif
    if ( lrepty ) then
       nend = nl             ! last  internal v-point in case of repeating grid
    else
       nend = nl - 1         ! last  internal v-point
    endif
    !
    do m = mfu, ml
       !
       nfm = kgrpnt(m,nf)
       nlm = kgrpnt(m,nl)
       !
       flux(nfm,:,2) = 0.
       flux(nlm,:,2) = 0.
       !
       do n = nsta, nend
          !
          md = m - 1
          nu = n + 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          num = kgrpnt(m ,nu)
          !
          ! compute effective horizontal diffusivity coefficient in v-point
          !
          if ( hdiff > 0. ) then
             !
             dif2d = hdiff
             !
          else
             !
             if ( ihvisc == 2 .or. ihvisc == 3 ) then
                !
                dif2d = 0.5 * ( vnu2d(nm) + vnu2d(nmd) ) / psm
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
          dxl = gvv(nm)
          dyl = guv(nm)
          !
          stabmx = 1./(2.*dt*(1./(dxl*dxl) + 1./(dyl*dyl)))
          !
          if ( .not. dif2d < stabmx ) dif2d = stabmx
          !
          do k = 1, kmax
             !
             kd = max(   1,k-1)
             ku = min(kmax,k+1)
             !
             flux(nm,k,2) = 0.5 * ( hkso(nm,k) + hkso(num,k) ) * dif2d * gvv(nm) / guv(nm)    *   &

!            --- compute d z(k) / dy

                            0.5 * ( zkso(num,k-1) - zkso(nm,k-1) + zkso(num,k) - zkso(nm,k) ) *   &

!            --- compute d c(k-1/2) / dz

                            ( rpo(num,kd) + rpo(nm,kd) - rpo(num,ku) - rpo(nm,ku) ) * 2. /        &
                            ( (ku-k) * (zkso(num,k -1) - zkso(num,ku)) +                          &
                              (k-kd) * (zkso(num,kd-1) - zkso(num,k )) +                          &
                              (ku-k) * (zkso(nm ,k -1) - zkso(nm ,ku)) +                          &
                              (k-kd) * (zkso(nm ,kd-1) - zkso(nm ,k )) )
             !
          enddo
          !
       enddo
       !
    enddo
    !
    ! synchronize fluxes at appropriate boundaries in case of repeating grid
    !
    call periodic ( flux(1,1,2), kgrpnt, 1, kmax )
    !
    ! update right-hand side with flux parts of anti-creep terms
    !
    ! not computed for end point at subdomain interface since, this end point is owned by the neighbouring subdomain
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    do k = 1, kmax
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
             if ( hs(nm) > epsdry ) then
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd == 1 ) nmd = nm
                if ( ndm == 1 ) ndm = nm
                !
                rhsc(nm,k) = rhsc(nm,k) - ( flux(nm,k,1) - flux(nmd,k,1) + flux(nm,k,2) - flux(ndm,k,2) ) / gsqs(nm)
                !
             endif
             !
          enddo
       enddo
       !
    enddo
    !
    ! compute constituent at interface between layers k and k+1 (stored in wrk)
    ! (Eq. 4.87)
    !
    do n = nf, nlu
       do m = mf, mlu
          !
          nm = kgrpnt(m,n)
          !
          do k = 1, kmax-1
             wrk(nm,k) = ( hkso(nm,k) * rpo(nm,k+1) + hkso(nm,k+1) * rpo(nm,k) ) / ( hkso(nm,k) + hkso(nm,k+1) )
          enddo
          !
       enddo
    enddo
    !
    ! update right-hand side with remaining (B and C) parts of anti-creep terms
    ! (Eq. 4.82)
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          md = m - 1
          mu = m + 1
          nd = n - 1
          nu = n + 1
          !
          nm   = kgrpnt(m ,n )
          nmd  = kgrpnt(md,n )
          nmu  = kgrpnt(mu,n )
          ndm  = kgrpnt(m ,nd)
          num  = kgrpnt(m ,nu)
          ndmd = kgrpnt(md,nd)
          !
          ! for permanently dry neighbours, corresponding values will be mirrored
          !
          if ( nmd  == 1 ) nmd  = nm
          if ( nmu  == 1 ) nmu  = nm
          if ( ndm  == 1 ) ndm  = nm
          if ( num  == 1 ) num  = nm
          if ( ndmd == 1 ) ndmd = nmd
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
                dif2d = 0.25 * ( vnu2d(nm) + vnu2d(ndm) + vnu2d(nmd) + vnu2d(ndmd) ) / psm
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
          dxl = 0.5 * ( gvv(nm) + gvv(ndm) )
          dyl = 0.5 * ( guu(nm) + guu(nmd) )
          !
          stabmx = 1./(2.*dt*(1./(dxl*dxl) + 1./(dyl*dyl)))
          !
          if ( .not. dif2d < stabmx ) dif2d = stabmx
          !
          do k = 1, kmax-1
             !
             ctrm1 = dif2d / dxl**2
             !
             ctrn1 = dif2d / dyl**2
             !
             ! compute d z(k) / dx and d z(k) / dy at interface between layers k and k+1
             !
             ctrm2 = 0.5 * ( zkso(nmu,k) - zkso(nmd,k) )
             !
             ctrn2 = 0.5 * ( zkso(num,k) - zkso(ndm,k) )
             !
             ! compute dc / dx  and dc / dy at interface between layers k and k+1
             !
             ctrm3 = 0.5 * ( wrk(nmu,k) - wrk(nmd,k) )
             !
             ctrn3 = 0.5 * ( wrk(num,k) - wrk(ndm,k) )
             !
             ! compute second part of anti-creep in x- and y-direction at interface between layers k and k+1
             ! (Eqs. 4.85 and 4.86)
             !
             ctrm4 = ctrm1 * ctrm2 * ctrm3
             !
             ctrn4 = ctrn1 * ctrn2 * ctrn3
             !
             ! compute third part of anti-creep in x- and y-direction at interface between layers k and k+1
             ! (Eq. 4.88)
             !
             ctrm5 = 2. * ctrm1 * ctrm2*ctrm2 / ( hkso(nm,k) + hkso(nm,k+1) )
             !
             ctrn5 = 2. * ctrn1 * ctrn2*ctrn2 / ( hkso(nm,k) + hkso(nm,k+1) )
             !
             if ( hs(nm) > epsdry ) then
                !
                ! second part is taken explicitly
                !
                rhsc(nm,k  ) = rhsc(nm,k  ) + ctrm4 + ctrn4                    ! at interface between layers k and k+1
                rhsc(nm,k+1) = rhsc(nm,k+1) - ctrm4 - ctrn4                    ! at interface between layers k-1 and k
                !
                ! third part is taken implicitly
                !
                amatc(nm,k  ,1) = amatc(nm,k  ,1) + ctrm5 + ctrn5              ! at interface between layers k and k+1
                amatc(nm,k  ,3) = amatc(nm,k  ,3) - ctrm5 - ctrn5
                amatc(nm,k+1,1) = amatc(nm,k+1,1) + ctrm5 + ctrn5              ! at interface between layers k-1 and k
                amatc(nm,k+1,2) = amatc(nm,k+1,2) - ctrm5 - ctrn5
                !
             endif
             !
          enddo
          !
       enddo
    enddo
    !
end subroutine SwashAntiCreep2DH
