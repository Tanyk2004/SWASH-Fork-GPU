subroutine SwashReynoldsStress
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
!    1.00: Tom Bogaard
!
!   Updates
!
!    1.00, January 2012: New subroutine
!
!   Purpose
!
!   Calculates the (non)linear Reynolds stress tensor
!
!   Method
!
!   Using series expansion arguments a general relationship between
!   Reynolds stresses and main strains and vorticities can be found
!
!   A linear variant of this constitutive relation results in isotropic
!   eddy viscosity and a nonlinear one yields anisotropic eddy viscosity
!
!   The employed nonlinear variant is quadratic in the mean velocity
!   gradients and is thus the lowest order. For details see
!
!   C.G. Speziale, "Analytical methods for the development of Reynolds-stress
!   closures in turbulence", Ann. Rev. Fluid Mech., vol. 23, p. 107-157, 1991
!
!   The resulting constitutive relation must be combined with a 3D k-epsilon model
!
!   The first option is the linear k-epsilon model of Launder and Spalding (1974)
!
!   The second option is the nonlinear k-epsilon model of Speziale (1987) which
!   satisfies the frame-indifference in the limit of two-dimensional turbulence
!
!   C.G. Speziale, "On nonlinear k-l and k-epsilon models of turbulence",
!   J. Fluid Mech., vol. 178, p. 459-475, 1987
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
!   Parameter variables
!
    real, parameter :: c1  =  0.4536 ! first closure constant for nonlinear Reynolds stress tensor
    real, parameter :: c2  =  0.3024 ! second closure constant for nonlinear Reynolds stress tensor
    real, parameter :: c3  = -0.1512 ! third closure constant for nonlinear Reynolds stress tensor
    real, parameter :: cmu =  0.09   ! closure constant for standard k-eps model
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter over vertical layers
    integer       :: kd       ! index of layer k-1
    integer       :: ku       ! index of layer k+1
    integer       :: m        ! loop counter in x-direction
    integer       :: md       ! index of point m-1
    integer       :: mlast    ! index of last internal u-point
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter in y-direction
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: ndmd     ! pointer to m-1,n-1
    integer       :: ndmu     ! pointer to m+1,n-1
    integer       :: nlast    ! index of last internal v-point
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
    integer       :: numd     ! pointer to m-1,n+1
    integer       :: numu     ! pointer to m+1,n+1
    !
    real          :: fac      ! a factor
    real          :: ltke     ! local turbulent kinetic energy
    real          :: lvisc    ! local eddy viscosity
    real          :: dz       ! local layer thickness
    real          :: dudx     ! gradient of u-velocity in x-direction
    real          :: dudy     ! gradient of u-velocity in y-direction
    real          :: dudz     ! gradient of u-velocity in vertical direction
    real          :: dvdx     ! gradient of v-velocity in x-direction
    real          :: dvdy     ! gradient of v-velocity in y-direction
    real          :: dvdz     ! gradient of v-velocity in vertical direction
    real          :: dwdx     ! gradient of w-velocity in x-direction
    real          :: dwdy     ! gradient of w-velocity in y-direction
    real          :: dwdz     ! gradient of w-velocity in vertical direction
    real          :: rdx      ! reciprocal of mesh size
    real          :: s1xx     ! first kind of mean strain tensor (xx comp.)
    real          :: s1xy     ! first kind of mean strain tensor (xy comp.)
    real          :: s1xz     ! first kind of mean strain tensor (xz comp.)
    real          :: s1yx     ! first kind of mean strain tensor (yx comp.)
    real          :: s1yy     ! first kind of mean strain tensor (yy comp.)
    real          :: s1yz     ! first kind of mean strain tensor (yz comp.)
    real          :: s1zx     ! first kind of mean strain tensor (zx comp.)
    real          :: s1zy     ! first kind of mean strain tensor (zy comp.)
    real          :: s1zz     ! first kind of mean strain tensor (zz comp.)
    real          :: s2xx     ! second kind of mean strain tensor (xx comp.)
    real          :: s2xy     ! second kind of mean strain tensor (xy comp.)
    real          :: s2xz     ! second kind of mean strain tensor (xz comp.)
    real          :: s2yx     ! second kind of mean strain tensor (yx comp.)
    real          :: s2yy     ! second kind of mean strain tensor (yy comp.)
    real          :: s2yz     ! second kind of mean strain tensor (yz comp.)
    real          :: s2zx     ! second kind of mean strain tensor (zx comp.)
    real          :: s2zy     ! second kind of mean strain tensor (zy comp.)
    real          :: s2zz     ! second kind of mean strain tensor (zz comp.)
    real          :: s3xx     ! third kind of mean strain tensor (xx comp.)
    real          :: s3xy     ! third kind of mean strain tensor (xy comp.)
    real          :: s3xz     ! third kind of mean strain tensor (xz comp.)
    real          :: s3yx     ! third kind of mean strain tensor (yx comp.)
    real          :: s3yy     ! third kind of mean strain tensor (yy comp.)
    real          :: s3yz     ! third kind of mean strain tensor (yz comp.)
    real          :: s3zx     ! third kind of mean strain tensor (zx comp.)
    real          :: s3zy     ! third kind of mean strain tensor (zy comp.)
    real          :: s3zz     ! third kind of mean strain tensor (zz comp.)
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashReynoldsStress')

    if ( oned ) then
       !
       rdx = 1./dx
       !
       ! first, compute linear part of the Reynolds stress tensor
       !
       ! for u-momentum equation
       !
       ! -u'u'
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,1)
             nmd = kgrpnt(md,1)
             !
             rsuu(nm,k) = ( vnu3d(nm,k-1) + vnu3d(nm,k) ) * ( u0(nm,k) - u0(nmd,k) ) * rdx
             !
          enddo
          !
       enddo
       !
       ! -u'w'
       !
       do k = 0, kmax
          !
          kd = max(k  ,1   )
          ku = min(k+1,kmax)
          !
          do m = mf+1, ml-1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             dz   =  0.5 * ( hkum(nm,kd) + hkum(nm,ku) )
             dudz = ( u0(nm ,kd) - u0(nm,ku) ) / dz
             dwdx = ( w0(nmu,k ) - w0(nm,k ) ) * rdx
             !
             rsuw(nm,k) = 0.5 * ( vnu3d(nm,k) + vnu3d(nmu,k) ) * ( dudz + dwdx )
             !
          enddo
          !
       enddo
       !
       ! for w-momentum equation
       !
       ! -w'u'
       !
       do k = 1, kmax-1
          !
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             dz   =  0.5 * ( hkum(nm,k) + hkum(nm,k+1) )
             dwdx = ( w0(nmu,k) - w0(nm,k  ) ) * rdx
             dudz = ( u0(nm ,k) - u0(nm,k+1) ) / dz
             !
             rswu(nm,k) = 0.5 * ( vnu3d(nm,k) + vnu3d(nmu,k) ) * ( dwdx + dudz )
             !
          enddo
          !
       enddo
       !
       ! -w'w'
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m,1)
             !
             rsww(nm,k) = ( vnu3d(nm,k-1) + vnu3d(nm,k) ) * ( w0(nm,k-1) - w0(nm,k) ) / hks(nm,k)
             !
          enddo
          !
       enddo
       !
       ! next, compute nonlinear parts, if appropriate
       !
       if ( iturb == 3 ) then
          !
          ! for u-momentum equation
          !
          ! -u'u'
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
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
                lvisc = 0.5 * ( vnu3d(nm,k-1  ) + vnu3d(nm,k  ) )
                ltke  = 0.5 * ( rtur (nm,k-1,1) + rtur (nm,k,1) )
                !
                fac  = lvisc*lvisc / ( cmu*ltke )
                !
                dudx = ( u0(nm,k) - u0(nmd,k) ) * rdx
                dudz = 0.25 * ( u0(nm,kd) + u0(nmd,kd) - u0(nm,ku) - u0(nmd,ku) ) / hks(nm,k)
                !
                dwdx = 0.25 * ( w0(nmu,k-1) - w0(nmd,k-1) + w0(nmu,k) - w0(nmd,k) ) * rdx
                dwdz = ( w0(nm,k-1) - w0(nm,k) ) / hks(nm,k)
                !
                s1xx = dudx * dudx + dudz * dudz
                s1zz = dwdx * dwdx + dwdz * dwdz
                !
                s2xx = dudx * dudx + dudz * dwdx
                s2zz = dwdx * dudz + dwdz * dwdz
                !
                s3xx = dudx * dudx + dwdx * dwdx
                s3zz = dudz * dudz + dwdz * dwdz
                !
                rsuu(nm,k) = rsuu(nm,k) - fac * ( c1 * ( s1xx - ( s1xx + s1zz ) / 3. ) + &
                                                  c2 * ( s2xx - ( s2xx + s2zz ) / 3. ) + &
                                                  c3 * ( s3xx - ( s3xx + s3zz ) / 3. ) )
                !
             enddo
             !
          enddo
          !
          ! -u'w'
          !
          do k = 0, kmax
             !
             kd = max(k  ,1   )
             ku = min(k+1,kmax)
             !
             do m = mf+1, ml-1
                !
                md = m - 1
                mu = m + 1
                !
                nm  = kgrpnt(m ,1)
                nmd = kgrpnt(md,1)
                nmu = kgrpnt(mu,1)
                !
                lvisc = 0.5 * ( vnu3d(nm,k  ) + vnu3d(nmu,k  ) )
                ltke  = 0.5 * ( rtur (nm,k,1) + rtur (nmu,k,1) )
                !
                fac  = lvisc*lvisc / ( cmu*ltke )
                !
                dz   =  0.5 * ( hkum(nm,kd) + hkum(nm,ku) )
                !
                dudx = 0.25 * ( u0(nmu,kd) - u0(nmd,kd) + u0(nmu,ku) - u0(nmd,ku) ) * rdx
                dudz = ( u0(nm,kd) - u0(nm,ku) ) / dz
                !
                dwdx = ( w0(nmu,k) - w0(nm,k) ) * rdx
                dwdz = 0.25 * ( w0(nm,kd-1) + w0(nmu,kd-1) - w0(nm,ku) - w0(nmu,ku) ) / dz
                !
                s1xz = dudx * dwdx + dudz * dwdz
                s2xz = 0.5 * ( dudx * dudz + dwdx * dudx + dudz * dwdz + dwdz * dwdx )
                s3xz = dudx * dudz + dwdx * dwdz
                !
                rsuw(nm,k) = rsuw(nm,k) - fac * ( c1 * s1xz + c2 * s2xz + c3 * s3xz )
                !
             enddo
             !
          enddo
          !
          ! for w-momentum equation
          !
          ! -w'u'
          !
          do k = 1, kmax-1
             !
             do m = mf, ml
                !
                md = m - 1
                mu = m + 1
                !
                if ( LMXF .and. md < mf ) md = mf
                !
                nm  = kgrpnt(m ,1)
                nmd = kgrpnt(md,1)
                nmu = kgrpnt(mu,1)
                !
                lvisc = 0.5 * ( vnu3d(nm,k  ) + vnu3d(nmu,k  ) )
                ltke  = 0.5 * ( rtur (nm,k,1) + rtur (nmu,k,1) )
                !
                fac  = lvisc*lvisc / ( cmu*ltke )
                !
                dz   =  0.5 * ( hkum(nm,k) + hkum(nm,k+1) )
                !
                dudx = 0.25 * ( u0(nmu,k) - u0(nmd,k) + u0(nmu,k+1) - u0(nmd,k+1) ) * rdx
                dudz = ( u0(nm,k) - u0(nm,k+1) ) / dz
                !
                dwdx = ( w0(nmu,k) - w0(nm,k) ) * rdx
                dwdz = 0.25 * ( w0(nm,k-1) + w0(nmu,k-1) - w0(nm,k+1) - w0(nmu,k+1) ) / dz
                !
                s1zx = dwdx * dudx + dwdz * dudz
                s2zx = 0.5 * ( dudx * dudz + dwdx * dudx + dudz * dwdz + dwdz * dwdx )
                s3zx = dudz * dudx + dwdz * dwdx
                !
                rswu(nm,k) = rswu(nm,k) - fac * ( c1 * s1zx + c2 * s2zx + c3 * s3zx )
                !
             enddo
             !
          enddo
          !
          ! -w'w'
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
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
                lvisc = 0.5 * ( vnu3d(nm,k-1  ) + vnu3d(nm,k  ) )
                ltke  = 0.5 * ( rtur (nm,k-1,1) + rtur (nm,k,1) )
                !
                fac  = lvisc*lvisc / ( cmu*ltke )
                !
                dudx = ( u0(nm ,k) - u0(nmd ,k) ) * rdx
                dudz = 0.25 * ( u0(nm,kd) + u0(nmd,kd) - u0(nm,ku) - u0(nmd,ku) ) / hks(nm,k)
                !
                dwdx = 0.25 * ( w0(nmu,k-1) - w0(nmd,k-1) + w0(nmu,k) - w0(nmd,k) ) * rdx
                dwdz = ( w0(nm,k-1) - w0(nm,k) ) / hks(nm,k)
                !
                s1xx = dudx * dudx + dudz * dudz
                s1zz = dwdx * dwdx + dwdz * dwdz
                !
                s2xx = dudx * dudx + dudz * dwdx
                s2zz = dwdx * dudz + dwdz * dwdz
                !
                s3xx = dudx * dudx + dwdx * dwdx
                s3zz = dudz * dudz + dwdz * dwdz
                !
                rsww(nm,k) = rsww(nm,k) - fac * ( c1 * ( s1zz - ( s1xx + s1zz ) / 3. ) + &
                                                  c2 * ( s2zz - ( s2xx + s2zz ) / 3. ) + &
                                                  c3 * ( s3zz - ( s3xx + s3zz ) / 3. ) )
                !
             enddo
             !
          enddo
          !
       endif
       !
    else
       !
       ! first, compute linear part of the Reynolds stress tensor
       !
       ! for u-momentum equation
       !
       if ( lreptx ) then
          mlast = ml               ! last  internal u-point in case of repeating grid
       else
          mlast = ml - 1           ! last  internal u-point
       endif
       !
       ! -u'u'
       !
       do k = 1, kmax
          !
          do n = nfu, nl
             do m = mfu, mlast+1
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
                !
                rsuu(nm,k) = 2. * ( vnu3d(nm,k-1) + vnu3d(nm,k) ) * ( u0(nm,k) - u0(nmd,k) ) / ( gvv(nm) + gvv(ndm) )
                !
             enddo
          enddo
          !
       enddo
       !
       ! -u'v'
       !
       do k = 1, kmax
          !
          do n = nf, nl
             do m = mf+1, mlast
                !
                md = m - 1
                mu = m + 1
                nu = n + 1
                !
                nm   = kgrpnt(m ,n )
                nmd  = kgrpnt(md,n )
                nmu  = kgrpnt(mu,n )
                num  = kgrpnt(m ,nu)
                numu = kgrpnt(mu,nu)
                !
                lvisc = 0.125 * ( vnu3d(nm,k-1) + vnu3d(nmu,k-1) + vnu3d(num,k-1) + vnu3d(numu,k-1) + vnu3d(nm,k) + vnu3d(nmu,k) + vnu3d(num,k) + vnu3d(numu,k) )
                !
                dudy = 2. * ( u0(num,k) - u0(nm,k) ) / ( guu(nm) + guu(num) )
                dvdx = 2. * ( v0(nmu,k) - v0(nm,k) ) / ( gvv(nm) + gvv(nmu) )
                !
                rsuv(nm,k) = lvisc * ( dudy + dvdx )
                !
             enddo
          enddo
          !
       enddo
       !
       ! -u'w'
       !
       do k = 0, kmax
          !
          kd = max(k  ,1   )
          ku = min(k+1,kmax)
          !
          do n = nfu, nl
             do m = mf+1, mlast
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                dz   =  0.5 * ( hkum(nm,kd) + hkum(nm,ku) )
                dudz = ( u0(nm ,kd) - u0(nm,ku) ) / dz
                dwdx = ( w0(nmu,k ) - w0(nm,k ) ) / gvu(nm)
                !
                rsuw(nm,k) = 0.5 * ( vnu3d(nm,k) + vnu3d(nmu,k) ) * ( dudz + dwdx )
                !
             enddo
          enddo
          !
       enddo
       !
       ! for v-momentum equation
       !
       if ( lrepty ) then
          nlast = nl               ! last  internal v-point in case of repeating grid
       else
          nlast = nl - 1           ! last  internal v-point
       endif
       !
       ! -v'u'
       !
       do k = 1, kmax
          !
          do m = mf, ml
             do n = nf+1, nlast
                !
                nd = n - 1
                nu = n + 1
                mu = m + 1
                !
                nm   = kgrpnt(m ,n )
                ndm  = kgrpnt(m ,nd)
                num  = kgrpnt(m ,nu)
                nmu  = kgrpnt(mu,n )
                numu = kgrpnt(mu,nu)
                !
                lvisc = 0.125 * ( vnu3d(nm,k-1) + vnu3d(num,k-1) + vnu3d(nmu,k-1) + vnu3d(numu,k-1) + vnu3d(nm,k) + vnu3d(num,k) + vnu3d(nmu,k) + vnu3d(numu,k) )
                !
                dvdx = 2. * ( v0(nmu,k) - v0(nm,k) ) / ( gvv(nm) + gvv(nmu) )
                dudy = 2. * ( u0(num,k) - u0(nm,k) ) / ( guu(nm) + guu(num) )
                !
                rsvu(nm,k) = lvisc * ( dvdx + dudy )
                !
             enddo
          enddo
          !
       enddo
       !
       ! -v'v'
       !
       do k = 1, kmax
          !
          do m = mfu, ml
             do n = nfu, nlast+1
                !
                nd = n - 1
                md = m - 1
                !
                nm  = kgrpnt(m ,n )
                ndm = kgrpnt(m ,nd)
                nmd = kgrpnt(md,n )
                !
                rsvv(nm,k) = 2. * ( vnu3d(nm,k-1) + vnu3d(nm,k) ) * ( v0(nm,k) - v0(ndm,k) ) / ( guu(nm) + guu(nmd) )
                !
             enddo
          enddo
          !
       enddo
       !
       ! -v'w'
       !
       do k = 0, kmax
          !
          kd = max(k  ,1   )
          ku = min(k+1,kmax)
          !
          do m = mfu, ml
             do n = nf+1, nlast
                !
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                dz   =  0.5 * ( hkvm(nm,kd) + hkvm(nm,ku) )
                dvdz = ( v0(nm ,kd) - v0(nm,ku) ) / dz
                dwdy = ( w0(num,k ) - w0(nm,k ) ) / guv(nm)
                !
                rsvw(nm,k) = 0.5 * ( vnu3d(nm,k) + vnu3d(num,k) ) * ( dvdz + dwdy )
                !
             enddo
          enddo
          !
       enddo
       !
       ! for w-momentum equation
       !
       ! -w'u'
       !
       do k = 1, kmax-1
          !
          do n = nfu, nl
             do m = mf, ml
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                dz   =  0.5 * ( hkum(nm,k) + hkum(nm,k+1) )
                dwdx = ( w0(nmu,k) - w0(nm,k  ) ) / gvu(nm)
                dudz = ( u0(nm ,k) - u0(nm,k+1) ) / dz
                !
                rswu(nm,k) = 0.5 * ( vnu3d(nm,k) + vnu3d(nmu,k) ) * ( dwdx + dudz )
                !
             enddo
          enddo
          !
       enddo
       !
       ! -w'v'
       !
       do k = 1, kmax-1
          !
          do n = nf, nl
             do m = mfu, ml
                !
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                dz   =  0.5 * ( hkvm(nm,k) + hkvm(nm,k+1) )
                dwdy = ( w0(num,k) - w0(nm,k  ) ) / guv(nm)
                dvdz = ( v0(nm ,k) - v0(nm,k+1) ) / dz
                !
                rswv(nm,k) = 0.5 * ( vnu3d(nm,k) + vnu3d(num,k) ) * ( dwdy + dvdz )
                !
             enddo
          enddo
          !
       enddo
       !
       ! -w'w'
       !
       do k = 1, kmax
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                rsww(nm,k) = ( vnu3d(nm,k-1) + vnu3d(nm,k) ) * ( w0(nm,k-1) - w0(nm,k) ) / hks(nm,k)
                !
             enddo
          enddo
          !
       enddo
       !
       ! next, compute nonlinear parts, if appropriate
       !
       if ( iturb == 3 ) then
          !
          ! for u-momentum equation
          !
          ! -u'u'
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
             !
             do n = nfu, nl
                do m = mfu, mlast+1
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
                   ndmd = kgrpnt(md,nd)
                   ndmu = kgrpnt(mu,nd)
                   num  = kgrpnt(m ,nu)
                   numd = kgrpnt(md,nu)
                   !
                   lvisc = 0.5 * ( vnu3d(nm,k-1  ) + vnu3d(nm,k  ) )
                   ltke  = 0.5 * ( rtur (nm,k-1,1) + rtur (nm,k,1) )
                   !
                   fac  = lvisc*lvisc / ( cmu*ltke )
                   !
                   dudx = 2.   * ( u0(nm ,k) - u0(nmd ,k) ) / ( gvv(nm) + gvv(ndm) )
                   dudy = 0.5  * ( u0(num,k) + u0(numd,k) - u0(ndm,k) - u0(ndmd,k) ) / ( guu(nm) + guu(nmd) )
                   dudz = 0.25 * ( u0(nm,kd) + u0(nmd,kd) - u0(nm,ku) - u0(nmd,ku) ) / hks(nm,k)
                   !
                   dvdx = 0.5  * ( v0(nmu,k) + v0(ndmu,k) - v0(nmd,k) - v0(ndmd,k) ) / ( gvv(nm) + gvv(ndm) )
                   dvdy = 2.   * ( v0(nm ,k) - v0(ndm ,k) ) / ( guu(nm) + guu(nmd) )
                   dvdz = 0.25 * ( v0(nm,kd) + v0(ndm,kd) - v0(nm,ku) - v0(ndm,ku) ) / hks(nm,k)
                   !
                   dwdx = 0.5 * ( w0(nmu,k-1) - w0(nmd,k-1) + w0(nmu,k) - w0(nmd,k) ) / ( gvv(nm) + gvv(ndm) )
                   dwdy = 0.5 * ( w0(num,k-1) - w0(ndm,k-1) + w0(num,k) - w0(ndm,k) ) / ( guu(nm) + guu(nmd) )
                   dwdz = ( w0(nm,k-1) - w0(nm,k) ) / hks(nm,k)
                   !
                   s1xx = dudx * dudx + dudy * dudy + dudz * dudz
                   s1yy = dvdx * dvdx + dvdy * dvdy + dvdz * dvdz
                   s1zz = dwdx * dwdx + dwdy * dwdy + dwdz * dwdz
                   !
                   s2xx = dudx * dudx + dudy * dvdx + dudz * dwdx
                   s2yy = dvdx * dudy + dvdy * dvdy + dvdz * dwdy
                   s2zz = dwdx * dudz + dwdy * dvdz + dwdz * dwdz
                   !
                   s3xx = dudx * dudx + dvdx * dvdx + dwdx * dwdx
                   s3yy = dudy * dudy + dvdy * dvdy + dwdy * dwdy
                   s3zz = dudz * dudz + dvdz * dvdz + dwdz * dwdz
                   !
                   rsuu(nm,k) = rsuu(nm,k) - fac * ( c1 * ( s1xx - ( s1xx + s1yy + s1zz ) / 3. ) + &
                                                     c2 * ( s2xx - ( s2xx + s2yy + s2zz ) / 3. ) + &
                                                     c3 * ( s3xx - ( s3xx + s3yy + s3zz ) / 3. ) )
                   !
                enddo
             enddo
             !
          enddo
          !
          ! -u'v'
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
             !
             do n = nf, nl
                do m = mf+1, mlast
                   !
                   md = m - 1
                   mu = m + 1
                   nd = n - 1
                   nu = n + 1
                   !
                   if ( LMYF .and. nd < nf ) nd = nf
                   !
                   nm   = kgrpnt(m ,n )
                   nmd  = kgrpnt(md,n )
                   nmu  = kgrpnt(mu,n )
                   ndm  = kgrpnt(m ,nd)
                   ndmu = kgrpnt(mu,nd)
                   num  = kgrpnt(m ,nu)
                   numd = kgrpnt(md,nu)
                   numu = kgrpnt(mu,nu)
                   !
                   lvisc = 0.125 * ( vnu3d(nm,k-1  ) + vnu3d(nmu,k-1  ) + vnu3d(num,k-1  ) + vnu3d(numu,k-1  ) + vnu3d(nm,k  ) + vnu3d(nmu,k  ) + vnu3d(num,k  ) + vnu3d(numu,k  ) )
                   ltke  = 0.125 * ( rtur (nm,k-1,1) + rtur (nmu,k-1,1) + rtur (num,k-1,1) + rtur (numu,k-1,1) + rtur (nm,k,1) + rtur (nmu,k,1) + rtur (num,k,1) + rtur (numu,k,1) )
                   !
                   fac  = lvisc*lvisc / ( cmu*ltke )
                   !
                   dz   =  0.5 * ( hkum(nm,k) + hkum(num,k) )
                   !
                   dudx = 0.5  * ( u0(nmu,k) + u0(numu,k) - u0(nmd,k) - u0(numd,k) ) / ( gvv(nm) + gvv(nmu) )
                   dudy = 2.   * ( u0(num,k) - u0(nm,k) ) / ( guu(nm) + guu(num) )
                   dudz = 0.25 * ( u0(num,kd) + u0(nm,kd) - u0(num,ku) - u0(nm,ku) ) / dz
                   !
                   dz   =  0.5 * ( hkvm(nm,k) + hkvm(nmu,k) )
                   !
                   dvdx = 2.  * ( v0(nmu,k) - v0(nm,k) ) / ( gvv(nm) + gvv(nmu) )
                   dvdy = 0.5 * ( v0(num,k) + v0(numu,k) - v0(ndm,k) - v0(ndmu,k) ) / ( guu(nm) + guu(num) )
                   dvdz = 0.25 * ( v0(nmu,kd) + v0(nm,kd) - v0(nmu,ku) - v0(nm,ku) ) / dz
                   !
                   dwdx = 0.5 * ( w0(nmu,k-1) + w0(numu,k-1) + w0(nmu,k) + w0(numu,k) - w0(nm,k-1) - w0(num,k-1) - w0(nm,k) - w0(num,k) ) / ( gvv(nm) + gvv(nmu) )
                   dwdy = 0.5 * ( w0(num,k-1) + w0(numu,k-1) + w0(num,k) + w0(numu,k) - w0(nm,k-1) - w0(nmu,k-1) - w0(nm,k) - w0(nmu,k) ) / ( guu(nm) + guu(num) )
                   !
                   s1xy = dudx * dvdx + dudy * dvdy + dudz * dvdz
                   s2xy = 0.5 * ( dudx * dudy + dvdx * dudx + dudy * dvdy + dvdy * dvdx + dudz * dwdy + dvdz * dwdx )
                   s3xy = dudx * dudy + dvdx * dvdy + dwdx * dwdy
                   !
                   rsuv(nm,k) = rsuv(nm,k) - fac * ( c1 * s1xy + c2 * s2xy + c3 * s3xy )
                   !
                enddo
             enddo
             !
          enddo
          !
          ! -u'w'
          !
          do k = 0, kmax
             !
             kd = max(k  ,1   )
             ku = min(k+1,kmax)
             !
             do n = nfu, nl
                do m = mf+1, mlast
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
                   ndmu = kgrpnt(mu,nd)
                   num  = kgrpnt(m ,nu)
                   numu = kgrpnt(mu,nu)
                   !
                   lvisc = 0.5 * ( vnu3d(nm,k  ) + vnu3d(nmu,k  ) )
                   ltke  = 0.5 * ( rtur (nm,k,1) + rtur (nmu,k,1) )
                   !
                   fac  = lvisc*lvisc / ( cmu*ltke )
                   !
                   dz   =  0.5 * ( hkum(nm,kd) + hkum(nm,ku) )
                   !
                   dudx = 0.25 * ( u0(nmu,kd) - u0(nmd,kd) + u0(nmu,ku) - u0(nmd,ku) ) / gvu(nm)
                   dudy = 0.25 * ( u0(num,kd) - u0(ndm,kd) + u0(num,ku) - u0(ndm,ku) ) / guu(nm)
                   dudz = ( u0(nm ,kd) - u0(nm,ku) ) / dz
                   !
                   dvdx = 0.25 * ( v0(nmu,kd) - v0(nm,kd) + v0(ndmu,kd) - v0(ndm,kd) + v0(nmu,ku) - v0(nm,ku) + v0(ndmu,ku) - v0(ndm,ku) ) / gvu(nm)
                   dvdz = 0.25 * ( v0(nmu,kd) + v0(nm,kd) + v0(ndmu,kd) + v0(ndm,kd) - v0(nmu,ku) - v0(nm,ku) - v0(ndmu,ku) - v0(ndm,ku) ) / dz
                   !
                   dwdx = ( w0(nmu,k) - w0(nm,k) ) / gvu(nm)
                   dwdy = 0.25 * ( w0(num,k   ) + w0(numu,k   ) - w0(ndm,k ) - w0(ndmu,k ) ) / guu(nm)
                   dwdz = 0.25 * ( w0(nm ,kd-1) + w0(nmu ,kd-1) - w0(nm ,ku) - w0(nmu ,ku) ) / dz
                   !
                   s1xz = dudx * dwdx + dudy * dwdy + dudz * dwdz
                   s2xz = 0.5 * ( dudx * dudz + dwdx * dudx + dudy * dvdz + dwdy * dvdx + dudz * dwdz + dwdz * dwdx )
                   s3xz = dudx * dudz + dvdx * dvdz + dwdx * dwdz
                   !
                   rsuw(nm,k) = rsuw(nm,k) - fac * ( c1 * s1xz + c2 * s2xz + c3 * s3xz )
                   !
                enddo
             enddo
             !
          enddo
          !
          ! for v-momentum equation
          !
          ! -v'u'
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
             !
             do m = mf, ml
                do n = nf+1, nlast
                   !
                   nd = n - 1
                   nu = n + 1
                   md = m - 1
                   mu = m + 1
                   !
                   if ( LMXF .and. md < mf ) md = mf
                   !
                   nm   = kgrpnt(m ,n )
                   ndm  = kgrpnt(m ,nd)
                   num  = kgrpnt(m ,nu)
                   nmu  = kgrpnt(mu,n )
                   numu = kgrpnt(mu,nu)
                   !
                   lvisc = 0.125 * ( vnu3d(nm,k-1  ) + vnu3d(num,k-1  ) + vnu3d(nmu,k-1  ) + vnu3d(numu,k-1  ) + vnu3d(nm,k  ) + vnu3d(num,k  ) + vnu3d(nmu,k  ) + vnu3d(numu,k  ) )
                   ltke  = 0.125 * ( rtur (nm,k-1,1) + rtur (num,k-1,1) + rtur (nmu,k-1,1) + rtur (numu,k-1,1) + rtur (nm,k,1) + rtur (nmu,k,1) + rtur (nmu,k,1) + rtur (numu,k,1) )
                   !
                   fac  = lvisc*lvisc / ( cmu*ltke )
                   !
                   dz   =  0.5 * ( hkum(nm,k) + hkum(num,k) )
                   !
                   dudx = 0.5  * ( u0(nmu,k) + u0(numu,k) - u0(nmd,k) - u0(numd,k) ) / ( gvv(nm) + gvv(nmu) )
                   dudy = 2.   * ( u0(num,k) - u0(nm,k) ) / ( guu(nm) + guu(num) )
                   dudz = 0.25 * ( u0(num,kd) + u0(nm,kd) - u0(num,ku) - u0(nm,ku) ) / dz
                   !
                   dz   =  0.5 * ( hkvm(nm,k) + hkvm(nmu,k) )
                   !
                   dvdx = 2.  * ( v0(nmu,k) - v0(nm,k) ) / ( gvv(nm) + gvv(nmu) )
                   dvdy = 0.5 * ( v0(num,k) + v0(numu,k) - v0(ndm,k) - v0(ndmu,k) ) / ( guu(nm) + guu(num) )
                   dvdz = 0.25 * ( v0(nmu,kd) + v0(nm,kd) - v0(nmu,ku) - v0(nm,ku) ) / dz
                   !
                   dwdx = 0.5 * ( w0(nmu,k-1) + w0(numu,k-1) + w0(nmu,k) + w0(numu,k) - w0(nm,k-1) - w0(num,k-1) - w0(nm,k) - w0(num,k) ) / ( gvv(nm) + gvv(nmu) )
                   dwdy = 0.5 * ( w0(num,k-1) + w0(numu,k-1) + w0(num,k) + w0(numu,k) - w0(nm,k-1) - w0(nmu,k-1) - w0(nm,k) - w0(nmu,k) ) / ( guu(nm) + guu(num) )
                   !
                   s1yx = dvdx * dudx + dvdy * dudy + dvdz * dudz
                   s2yx = 0.5 * ( dvdx * dudx + dudx * dudy + dvdy * dvdx + dudy * dvdy + dvdz * dwdx + dudz * dwdy )
                   s3yx = dudy * dudx + dvdy * dvdx + dwdy * dwdx
                   !
                   rsvu(nm,k) = rsvu(nm,k) - fac * ( c1 * s1yx + c2 * s2yx + c3 * s3yx )
                   !
                enddo
             enddo
             !
          enddo
          !
          ! -v'v'
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
             !
             do m = mfu, ml
                do n = nfu, nlast+1
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
                   ndmd = kgrpnt(md,nd)
                   ndmu = kgrpnt(mu,nd)
                   num  = kgrpnt(m ,nu)
                   numd = kgrpnt(md,nu)
                   !
                   lvisc = 0.5 * ( vnu3d(nm,k-1  ) + vnu3d(nm,k  ) )
                   ltke  = 0.5 * ( rtur (nm,k-1,1) + rtur (nm,k,1) )
                   !
                   fac  = lvisc*lvisc / ( cmu*ltke )
                   !
                   dudx = 2.   * ( u0(nm ,k) - u0(nmd ,k) ) / ( gvv(nm) + gvv(ndm) )
                   dudy = 0.5  * ( u0(num,k) + u0(numd,k) - u0(ndm,k) - u0(ndmd,k) ) / ( guu(nm) + guu(nmd) )
                   dudz = 0.25 * ( u0(nm,kd) + u0(nmd,kd) - u0(nm,ku) - u0(nmd,ku) ) / hks(nm,k)
                   !
                   dvdx = 0.5  * ( v0(nmu,k) + v0(ndmu,k) - v0(nmd,k) - v0(ndmd,k) ) / ( gvv(nm) + gvv(ndm) )
                   dvdy = 2.   * ( v0(nm ,k) - v0(ndm ,k) ) / ( guu(nm) + guu(nmd) )
                   dvdz = 0.25 * ( v0(nm,kd) + v0(ndm,kd) - v0(nm,ku) - v0(ndm,ku) ) / hks(nm,k)
                   !
                   dwdx = 0.5 * ( w0(nmu,k-1) - w0(nmd,k-1) + w0(nmu,k) - w0(nmd,k) ) / ( gvv(nm) + gvv(ndm) )
                   dwdy = 0.5 * ( w0(num,k-1) - w0(ndm,k-1) + w0(num,k) - w0(ndm,k) ) / ( guu(nm) + guu(nmd) )
                   dwdz = ( w0(nm,k-1) - w0(nm,k) ) / hks(nm,k)
                   !
                   s1xx = dudx * dudx + dudy * dudy + dudz * dudz
                   s1yy = dvdx * dvdx + dvdy * dvdy + dvdz * dvdz
                   s1zz = dwdx * dwdx + dwdy * dwdy + dwdz * dwdz
                   !
                   s2xx = dudx * dudx + dudy * dvdx + dudz * dwdx
                   s2yy = dvdx * dudy + dvdy * dvdy + dvdz * dwdy
                   s2zz = dwdx * dudz + dwdy * dvdz + dwdz * dwdz
                   !
                   s3xx = dudx * dudx + dvdx * dvdx + dwdx * dwdx
                   s3yy = dudy * dudy + dvdy * dvdy + dwdy * dwdy
                   s3zz = dudz * dudz + dvdz * dvdz + dwdz * dwdz
                   !
                   rsvv(nm,k) = rsvv(nm,k) - fac * ( c1 * ( s1yy - ( s1xx + s1yy + s1zz ) / 3. ) + &
                                                     c2 * ( s2yy - ( s2xx + s2yy + s2zz ) / 3. ) + &
                                                     c3 * ( s3yy - ( s3xx + s3yy + s3zz ) / 3. ) )
                   !
                enddo
             enddo
             !
          enddo
          !
          ! -v'w'
          !
          do k = 0, kmax
             !
             kd = max(k  ,1   )
             ku = min(k+1,kmax)
             !
             do m = mfu, ml
                do n = nf+1, nlast
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
                   numd = kgrpnt(md,nu)
                   numu = kgrpnt(mu,nu)
                   !
                   lvisc = 0.5 * ( vnu3d(nm,k  ) + vnu3d(num,k  ) )
                   ltke  = 0.5 * ( rtur (nm,k,1) + rtur (num,k,1) )
                   !
                   fac  = lvisc*lvisc / ( cmu*ltke )
                   !
                   dz   =  0.5 * ( hkvm(nm,kd) + hkvm(nm,ku) )
                   !
                   dudy = 0.25 * ( u0(num,kd) - u0(nm,kd) + u0(numd,kd) - u0(nmd,kd) + u0(num,ku) - u0(nm,ku) + u0(numd,ku) - u0(nmd,ku) ) / guv(nm)
                   dudz = 0.25 * ( u0(num,kd) + u0(nm,kd) + u0(numd,kd) + u0(nmd,kd) - u0(num,ku) - u0(nm,ku) - u0(numd,ku) - u0(nmd,ku) ) / dz
                   !
                   dvdx = 0.25 * ( v0(nmu,kd) - v0(nmd,kd) + v0(nmu,ku) - v0(nmd,ku) ) / gvv(nm)
                   dvdy = 0.25 * ( v0(num,kd) - v0(ndm,kd) + v0(num,ku) - v0(ndm,ku) ) / guv(nm)
                   dvdz = ( v0(nm ,kd) - v0(nm,ku) ) / dz
                   !
                   dwdx = 0.25 * ( w0(nmu,k   ) + w0(numu,k   ) - w0(nmd,k ) - w0(numd,k ) ) / gvv(nm)
                   dwdy = ( w0(num,k) - w0(nm,k) ) / guv(nm)
                   dwdz = 0.25 * ( w0(nm ,kd-1) + w0(num ,kd-1) - w0(nm ,ku) - w0(num ,ku) ) / dz
                   !
                   s1yz = dvdx * dwdx + dvdy * dwdy + dvdz * dwdz
                   s2yz = 0.5 * ( dvdx * dudz + dwdx * dudy + dvdy * dvdz + dwdy * dvdy + dvdz * dwdz + dwdz * dwdy )
                   s3yz = dudy * dudz + dvdy * dvdz + dwdy * dwdz
                   !
                   rsvw(nm,k) = rsvw(nm,k) - fac * ( c1 * s1yz + c2 * s2yz + c3 * s3yz )
                   !
                enddo
             enddo
             !
          enddo
          !
          ! for w-momentum equation
          !
          ! -w'u'
          !
          do k = 1, kmax-1
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   md = m - 1
                   mu = m + 1
                   nd = n - 1
                   nu = n + 1
                   !
                   if ( LMXF .and. md < mf ) md = mf
                   !
                   nm   = kgrpnt(m ,n )
                   nmd  = kgrpnt(md,n )
                   nmu  = kgrpnt(mu,n )
                   ndm  = kgrpnt(m ,nd)
                   ndmu = kgrpnt(mu,nd)
                   num  = kgrpnt(m ,nu)
                   numu = kgrpnt(mu,nu)
                   !
                   lvisc = 0.5 * ( vnu3d(nm,k  ) + vnu3d(nmu,k  ) )
                   ltke  = 0.5 * ( rtur (nm,k,1) + rtur (nmu,k,1) )
                   !
                   fac  = lvisc*lvisc / ( cmu*ltke )
                   !
                   dz   =  0.5 * ( hkum(nm,k) + hkum(nm,k+1) )
                   !
                   dudx = 0.25 * ( u0(nmu,k) - u0(nmd,k) + u0(nmu,k+1) - u0(nmd,k+1) ) / gvu(nm)
                   dudy = 0.25 * ( u0(num,k) - u0(ndm,k) + u0(num,k+1) - u0(ndm,k+1) ) / guu(nm)
                   dudz = ( u0(nm ,k) - u0(nm,k+1) ) / dz
                   !
                   dvdx = 0.25 * ( v0(nmu,k) - v0(nm,k) + v0(ndmu,k) - v0(ndm,k) + v0(nmu,k+1) - v0(nm,k+1) + v0(ndmu,k+1) - v0(ndm,k+1) ) / gvu(nm)
                   dvdz = 0.25 * ( v0(nmu,k) + v0(nm,k) + v0(ndmu,k) + v0(ndm,k) - v0(nmu,k+1) - v0(nm,k+1) - v0(ndmu,k+1) - v0(ndm,k+1) ) / dz
                   !
                   dwdx = ( w0(nmu,k) - w0(nm,k) ) / gvu(nm)
                   dwdy = 0.25 * ( w0(num,k  ) + w0(numu,k  ) - w0(ndm,k  ) - w0(ndmu,k  ) ) / guu(nm)
                   dwdz = 0.25 * ( w0(nm ,k-1) + w0(nmu ,k-1) - w0(nm ,k+1) - w0(nmu ,k+1) ) / dz
                   !
                   s1zx = dwdx * dudx + dwdy * dudy + dwdz * dudz
                   s2zx = 0.5 * ( dudx * dudz + dwdx * dudx + dudy * dvdz + dwdy * dvdx + dudz * dwdz + dwdz * dwdx )
                   s3zx = dudz * dudx + dvdz * dvdx + dwdz * dwdx
                   !
                   rswu(nm,k) = rswu(nm,k) - fac * ( c1 * s1zx + c2 * s2zx + c3 * s3zx )
                   !
                enddo
             enddo
             !
          enddo
          !
          ! -w'v'
          !
          do k = 1, kmax-1
             !
             do n = nf, nl
                do m = mfu, ml
                   !
                   md = m - 1
                   mu = m + 1
                   nd = n - 1
                   nu = n + 1
                   !
                   if ( LMYF .and. nd < nf ) nd = nf
                   !
                   nm   = kgrpnt(m ,n )
                   nmd  = kgrpnt(md,n )
                   nmu  = kgrpnt(mu,n )
                   ndm  = kgrpnt(m ,nd)
                   num  = kgrpnt(m ,nu)
                   numd = kgrpnt(md,nu)
                   numu = kgrpnt(mu,nu)
                   !
                   lvisc = 0.5 * ( vnu3d(nm,k  ) + vnu3d(num,k  ) )
                   ltke  = 0.5 * ( rtur (nm,k,1) + rtur (num,k,1) )
                   !
                   fac  = lvisc*lvisc / ( cmu*ltke )
                   !
                   dz   =  0.5 * ( hkvm(nm,k) + hkvm(nm,k+1) )
                   !
                   dudy = 0.25 * ( u0(num,k) - u0(nm,k) + u0(numd,k) - u0(nmd,k) + u0(num,k+1) - u0(nm,k+1) + u0(numd,k+1) - u0(nmd,k+1) ) / guv(nm)
                   dudz = 0.25 * ( u0(num,k) + u0(nm,k) + u0(numd,k) + u0(nmd,k) - u0(num,k+1) - u0(nm,k+1) - u0(numd,k+1) - u0(nmd,k+1) ) / dz
                   !
                   dvdx = 0.25 * ( v0(nmu,k) - v0(nmd,k) + v0(nmu,k+1) - v0(nmd,k+1) ) / gvv(nm)
                   dvdy = 0.25 * ( v0(num,k) - v0(ndm,k) + v0(num,k+1) - v0(ndm,k+1) ) / guv(nm)
                   dvdz = ( v0(nm ,k) - v0(nm,k+1) ) / dz
                   !
                   dwdx = 0.25 * ( w0(nmu,k  ) + w0(numu,k  ) - w0(nmd,k  ) - w0(numd,k  ) ) / gvv(nm)
                   dwdy = ( w0(num,k) - w0(nm,k) ) / guv(nm)
                   dwdz = 0.25 * ( w0(nm ,k-1) + w0(num ,k-1) - w0(nm ,k+1) - w0(num ,k+1) ) / dz
                   !
                   s1zy = dwdx * dvdx + dwdy * dvdy + dwdz * dvdz
                   s2zy = 0.5 * ( dvdx * dudz + dwdx * dudy + dvdy * dvdz + dwdy * dvdy + dvdz * dwdz + dwdz * dwdy )
                   s3zy = dudz * dudy + dvdz * dvdy + dwdz * dwdy
                   !
                   rswv(nm,k) = rswv(nm,k) - fac * ( c1 * s1zy + c2 * s2zy + c3 * s3zy )
                   !
                enddo
             enddo
             !
          enddo
          !
          ! -w'w'
          !
          do k = 1, kmax
             !
             kd = max(k-1,1   )
             ku = min(k+1,kmax)
             !
             do n = nfu, nl
                do m = mfu, ml
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
                   ndmd = kgrpnt(md,nd)
                   ndmu = kgrpnt(mu,nd)
                   num  = kgrpnt(m ,nu)
                   numd = kgrpnt(md,nu)
                   !
                   lvisc = 0.5 * ( vnu3d(nm,k-1  ) + vnu3d(nm,k  ) )
                   ltke  = 0.5 * ( rtur (nm,k-1,1) + rtur (nm,k,1) )
                   !
                   fac  = lvisc*lvisc / ( cmu*ltke )
                   !
                   dudx = 2.   * ( u0(nm ,k) - u0(nmd ,k) ) / ( gvv(nm) + gvv(ndm) )
                   dudy = 0.5  * ( u0(num,k) + u0(numd,k) - u0(ndm,k) - u0(ndmd,k) ) / ( guu(nm) + guu(nmd) )
                   dudz = 0.25 * ( u0(nm,kd) + u0(nmd,kd) - u0(nm,ku) - u0(nmd,ku) ) / hks(nm,k)
                   !
                   dvdx = 0.5  * ( v0(nmu,k) + v0(ndmu,k) - v0(nmd,k) - v0(ndmd,k) ) / ( gvv(nm) + gvv(ndm) )
                   dvdy = 2.   * ( v0(nm ,k) - v0(ndm ,k) ) / ( guu(nm) + guu(nmd) )
                   dvdz = 0.25 * ( v0(nm,kd) + v0(ndm,kd) - v0(nm,ku) - v0(ndm,ku) ) / hks(nm,k)
                   !
                   dwdx = 0.5 * ( w0(nmu,k-1) - w0(nmd,k-1) + w0(nmu,k) - w0(nmd,k) ) / ( gvv(nm) + gvv(ndm) )
                   dwdy = 0.5 * ( w0(num,k-1) - w0(ndm,k-1) + w0(num,k) - w0(ndm,k) ) / ( guu(nm) + guu(nmd) )
                   dwdz = ( w0(nm,k-1) - w0(nm,k) ) / hks(nm,k)
                   !
                   s1xx = dudx * dudx + dudy * dudy + dudz * dudz
                   s1yy = dvdx * dvdx + dvdy * dvdy + dvdz * dvdz
                   s1zz = dwdx * dwdx + dwdy * dwdy + dwdz * dwdz
                   !
                   s2xx = dudx * dudx + dudy * dvdx + dudz * dwdx
                   s2yy = dvdx * dudy + dvdy * dvdy + dvdz * dwdy
                   s2zz = dwdx * dudz + dwdy * dvdz + dwdz * dwdz
                   !
                   s3xx = dudx * dudx + dvdx * dvdx + dwdx * dwdx
                   s3yy = dudy * dudy + dvdy * dvdy + dwdy * dwdy
                   s3zz = dudz * dudz + dvdz * dvdz + dwdz * dwdz
                   !
                   rsww(nm,k) = rsww(nm,k) - fac * ( c1 * ( s1zz - ( s1xx + s1yy + s1zz ) / 3. ) + &
                                                     c2 * ( s2zz - ( s2xx + s2yy + s2zz ) / 3. ) + &
                                                     c3 * ( s3zz - ( s3xx + s3yy + s3zz ) / 3. ) )
                   !
                enddo
             enddo
             !
          enddo
          !
       endif
       !
    endif
    !
end subroutine SwashReynoldsStress
