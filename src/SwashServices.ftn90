!
!  SWASH service routines
!
!  Contents of this file
!
!     disprel
!     jonswap
!     fluxlim
!     sigmacoor
!     periodic
!     periodici
!
subroutine disprel ( htot, omega, kwav, cg, n )
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
!    1.00, January 2010: New subroutine
!
!   Purpose
!
!   Computes the wave number, group velocity and ratio of group and phase velocity
!
!   Method
!
!   Based on the dispersion relation for linear waves approximated with the formula of Guo (2002)
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
!
    implicit none
!
!   Argument variables
!
    real, intent(out) :: cg    ! group velocity
    real, intent(in)  :: htot  ! local water depth
    real, intent(out) :: kwav  ! wave number
    real, intent(in)  :: omega ! angular frequency
    real, intent(out) :: n     ! ratio of group and phase velocity
!
!   Parameter variables
!
    integer, parameter :: maxnit = 100 ! maximum number of iterations
    !
    real   , parameter :: eps = 0.01   ! convergence criterion
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: nit      ! number of iterations
    !
    real          :: fac      ! auxiliary factor
    real          :: fac1     ! auxiliary factor
    real          :: fac2     ! auxiliary factor
    real          :: fac3     ! auxiliary factor
    real          :: kh2      ! square of kwav*htot
    real          :: kh2o     ! square of kwav*htot at previous iteration
    real          :: rootdg   ! square root of htot/grav
    real          :: snd      ! dimensionless angular frequency
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'disprel')
    !
    rootdg = sqrt(htot/grav)
    snd    = omega * rootdg
    !
    if ( snd > 2.5 ) then
       !
       ! deep water
       !
       kwav = omega * omega / grav
       cg   = 0.5 * grav / omega
       n    = 0.5
       !
    else if ( snd < 1.e-6 ) then
       !
       ! very shallow water
       !
       kwav = snd / htot
       cg   = rootdg * grav
       n    = 1.
       !
    else
       !
       ! intermediate water depth
       !
       if ( .not.numdisp ) then
          !
          ! exact relation (Guo, 2002)
          !
          fac1 = ( omega * rootdg )**2.5
          fac2 = ( 1.-exp(-fac1) )**(-0.4)
          kwav = fac2 * omega * omega / grav
          !
       else
          !
          ! approximate relation based on the box scheme
          !
          if ( kpmax == 1 ) then
             !
             if ( snd < 1.8 ) then
                kwav = omega / sqrt( grav*htot - 0.25*omega*omega*htot*htot )
             else
                kwav = omega * omega / grav
             endif
             !
          else if ( kpmax == 2 ) then
             !
             fac  = snd * snd
             fac1 = (6.*fac-16.) / (fac/16.-1.)
             fac2 =      16.*fac / (fac/16.-1.)
             kh2  = -0.5*fac1 + 0.5*sqrt( fac1*fac1 -4.*fac2 )
             if ( .not. kh2 > 0. ) kh2 = -0.5*fac1 - 0.5*sqrt( fac1*fac1 - 4.*fac2 )
             kwav = sqrt(kh2) / htot
             !
          else if ( kpmax == 3 ) then
             !
             fac = snd * snd
             !
             fac1 =    fac/46656. - 1./1296.
             fac2 = 5.*fac/432.   - 5./54.
             fac3 = 5.*fac/12.    - 1.
             !
             nit = 0
             !
             ! initial value for (kh)^2
             !
             kh2  = fac*fac
             kh2o = 0.
             !
             ! Newton-Raphson iteration
             !
             do
                if ( abs(kh2o-kh2) < (eps*kh2) ) exit
                !
                nit  = nit + 1
                kh2o = kh2
                kh2  = kh2o - (fac1*kh2**3 + fac2*kh2**2 + fac3*kh2 + fac)/(3.*fac1*kh2**2 + 2.*fac2*kh2 + fac3)
                !
                if ( .not. nit < maxnit ) then
                   call msgerr (1, 'no convergence in dispersion relation')
                   kh2 = fac*fac
                   exit
                endif
                !
             enddo
             !
             kwav = sqrt(kh2) / htot
             !
          endif
          !
       endif
       !
       n  = 0.5 + kwav * htot / sinh( 2.*kwav*htot )
       cg = n * omega / kwav
       !
    endif
    !
end subroutine disprel
!
subroutine jonswap ( spec, fmin, df, nfreq, fp, gamma )
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
!    1.00, April 2010: New subroutine
!
!   Purpose
!
!   Calculates energy density spectrum of a Jonswap-type
!
!   Modules used
!
    use ocpcomm4
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                 :: nfreq ! number of frequencies
    !
    real, intent(in)                    :: df    ! increment in frequency space
    real, dimension(nfreq), intent(out) :: spec  ! nondimensional relative spectral density, equal to one at the peak
    real, intent(in)                    :: fmin  ! minimum frequency
    real, intent(in)                    :: fp    ! peak frequency
    real, intent(in)                    :: gamma ! peak-enhancement factor
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: n        ! loop counter
    !
    real          :: f        ! normalised frequency
    real          :: fac1     ! auxiliary factor
    real          :: fac2     ! auxiliary factor
    real          :: fac3     ! auxiliary factor
    real          :: sigma    ! peak-width parameter
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'jonswap')
    !
    f = fmin / fp
    !
    do n = 1, nfreq
       !
       if ( f < 1. ) then
          sigma = 0.07
       else
          sigma = 0.09
       endif
       !
       fac1 = f**(-5.)
       fac2 = exp(-1.25*(f**(-4.)))
       !
       if ( gamma /= 1. ) then
          !
          fac3 = gamma**(exp(-0.5*((f-1.)/sigma)**2.))
          !
       else
          !
          fac3 = 1.
          !
       endif
       !
       spec(n) = fac1 * fac2 * fac3
       !
       f = f + df/fp
       !
    enddo
    !
    spec = spec / maxval(spec)
    !
end subroutine jonswap
!
real function fluxlim ( a, b )
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
!    1.00, January 2010: New subroutine
!
!   Purpose
!
!   Computes the flux limiter for higher order correction
!
!   Method
!
!   The following classes of limiters can be considered:
!
!   1) linear higher order schemes (kappa-scheme) like BDF, CDS, QUICK, Fromm, etc.
!   2) Sweby Phi-limiter like minmod, superbee, etc.
!   3) R-kappa limiter like Van Leer
!   4) PL-kappa limiter like MUSCL, Koren, SMART, etc.
!
!   For details, see Ph.D thesis of M. Zijlema
!
!   Note that the flux limiter is multiplied with a solution gradient (=denominator)!
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
!
    implicit none
!
!   Argument variables
!
    real, intent(in) :: a        ! solution gradient as numerator
    real, intent(in) :: b        ! solution gradient as denominator
!
!   Local variables
!
    integer, save    :: ient = 0 ! number of entries in this subroutine
    real             :: fac      ! auxiliary factor
    real             :: r        ! ratio of consecutive solution gradients
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'fluxlim')
    !
    if ( abs(b) > 1.e-10 ) then
       !
       r = a / b
       !
    else
       !
       r = 0.
       !
    endif
    !
    if ( propsc == 3 ) then
       !
       ! linear higher order scheme (kappa-scheme)
       !
       fluxlim = 0.5 * ( 1. + kappa ) * a + 0.5 * ( 1. - kappa ) * b
       !
    else if ( propsc == 4 ) then
       !
       ! Sweby Phi-limiter
       !
       fluxlim = max( max( 0., min(phieby*r, 1.) ), min(r, phieby) ) * b
       !
    else if ( propsc == 5 ) then
       !
       ! R-kappa limiter
       !
       if ( r < 0. ) then
          !
          fluxlim = 0.
          !
       else if ( kappa < 0. ) then
          !
          if ( r <= 1. ) then
             !
             fluxlim = 2.*r * ( -r*r + (3.+kappa)*r - kappa )/((1.+r)**2)
             !
          else
             !
             fluxlim = ( (2.+kappa)*r - kappa )/( 1.+r )
             !
          endif
          !
       else
          !
          fluxlim = 2.*r * ( (1.+kappa)*r + 1. - kappa )/( (1.+r)**2 )
          !
       endif
       !
       fluxlim = fluxlim * b
       !
    else if ( propsc == 6 ) then
       !
       ! PL-kappa limiter
       !
       fac = min( mbound, 0.5 * ( 1. + kappa ) * r + 0.5 * (1. - kappa ) )
       !
       fluxlim = max( 0., min( mbound*r, fac ) ) * b
       !
    else
       !
       call msgerr ( 4, 'unknown higher order correction' )
       return
       !
    endif
    !
end function fluxlim
!
subroutine sigmacoor ( zk, h, ntot )
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
!    1.00, July 2010: New subroutine
!
!   Purpose
!
!   Computes layer interfaces
!
!   Method
!
!   Apply sigma transformation of Phillips (1957)
!
!   Note: zk(:,0) and zk(:,kmax) must be filled!
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: kmax, fixlay, depfix
    use m_genarr, only: indlay, hlay, hlaysh
!
    implicit none
!
!   Argument variables
!
    integer,               intent(in)    :: ntot ! number of grid elements
    !
    real   , dimension(*), intent(in)    :: h    ! water depth
    real   , dimension(*), intent(inout) :: zk   ! layer interfaces
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! index
    integer       :: k        ! loop counter over layers
    integer       :: n        ! loop counter over grid elements
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'sigmacoor')
    !
    if ( .not.fixlay ) then
       !
       ! sigma layers
       !
       do n = 1, ntot
          do k = 1, kmax-1
             !
             j = n + (k-1)*ntot
             !
             zk(j+ntot) = zk(j) - h(n) * hlay(k)
             !
          enddo
       enddo
       !
    else
       !
       ! fixed layers
       !
       do n = 1, ntot
          !
          if ( h(n) > depfix ) then
             !
             ! use sigma transformation
             !
             do k = 1, kmax-1
                !
                if ( indlay(k) == 1 ) then
                   !
                   ! sigma layer
                   !
                   j = n + (k-1)*ntot
                   !
                   zk(j+ntot) = zk(j) - (h(n) - depfix) * hlay(k)
                   !
                else
                   !
                   ! fixed layer
                   !
                   j = n + (k-1)*ntot
                   !
                   zk(j+ntot) = zk(j) - hlay(k)
                   !
                endif
                !
             enddo
             !
          else
             !
             ! use shadow sigma transformation
             !
             do k = 1, kmax-1
                !
                j = n + (k-1)*ntot
                !
                zk(j+ntot) = zk(j) - h(n) * hlaysh(k)
                !
             enddo
             !
          endif
          !
       enddo
       !
    endif
    !
end subroutine sigmacoor
!
subroutine periodic ( u, kgrpnt, ks, ke )
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
!    1.00, October 2012: New subroutine
!
!   Purpose
!
!   Fills virtual unknowns using periodicity
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata
    use m_parall
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc,myc),     intent(in)    :: kgrpnt ! index table containing the address of each (active) grid point
                                                             ! =1; not active grid point
                                                             ! >1; active grid point
    integer,                         intent(in)    :: ke     ! last index in vertical direction
    integer,                         intent(in)    :: ks     ! first index in vertical direction
    real   , dimension(mcgrd,ks:ke), intent(inout) :: u      ! unknown (water level, velocity, pressure, etc.)
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter in vertical direction
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    integer       :: nfm      ! pointer to m,nf
    integer       :: nfum     ! pointer to m,nfu
    integer       :: nlm      ! pointer to m,nl
    integer       :: nlum     ! pointer to m,nlu
    integer       :: nmf      ! pointer to mf,n
    integer       :: nmfu     ! pointer to mfu,n
    integer       :: nml      ! pointer to ml,n
    integer       :: nmlu     ! pointer to mlu,n
    integer       :: npnts    ! number of points to be send/receive
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'periodic')
    !
    ! if no periodicity, return
    !
    if ( .not.lreptx .and. .not.lrepty ) return
    !
    if ( lreptx ) then
       !
       if ( kpart /= 4 ) then
          !
          do n = nf, nlu
             !
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             nml  = kgrpnt(ml ,n)
             nmlu = kgrpnt(mlu,n)
             !
             do k = ks, ke
                u(nmf ,k) = u(nml ,k)
                u(nmlu,k) = u(nmfu,k)
             enddo
             !
          enddo
          !
       else
          !
          ! for all layers do
          !
          do k = ks, ke
             !
             npnts = nlu - nf + 1
             !
             if ( INODE == 1 ) then
                !
                ! store data to be sent in array WRKEXG
                !
                do n = nf, nlu
                   nmfu = kgrpnt(mfu,n)
                   WRKEXG(n-nf+1) = u(nmfu,k)
                enddo
                !
                ! send array WRKEXG
                !
                call SWSENDNB ( WRKEXG, npnts, SWREAL, NPROC, 2 )
                if (STPNOW()) return
                !
                ! receive next array and store in WRKEXG
                !
                call SWRECVNB ( WRKEXG, npnts, SWREAL, NPROC, 2 )
                if (STPNOW()) return
                !
                ! store the received data
                !
                do n = nf, nlu
                   nmf = kgrpnt(mf,n)
                   u(nmf,k) = WRKEXG(n-nf+1)
                enddo
                !
             else if ( INODE == NPROC ) then
                !
                ! store data to be sent in array WRKEXG
                !
                do n = nf, nlu
                   nml = kgrpnt(ml,n)
                   WRKEXG(n-nf+1) = u(nml,k)
                enddo
                !
                ! send array WRKEXG
                !
                call SWSENDNB ( WRKEXG, npnts, SWREAL, 1, 2 )
                if (STPNOW()) return
                !
                ! receive next array and store in WRKEXG
                !
                call SWRECVNB ( WRKEXG, npnts, SWREAL, 1, 2 )
                if (STPNOW()) return
                !
                ! store the received data
                !
                do n = nf, nlu
                   nmlu = kgrpnt(mlu,n)
                   u(nmlu,k) = WRKEXG(n-nf+1)
                enddo
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
    if ( lrepty ) then
       !
       if ( kpart /= 3 ) then
          !
          do m = mf, mlu
             !
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             nlm  = kgrpnt(m,nl )
             nlum = kgrpnt(m,nlu)
             !
             do k = ks, ke
                u(nfm ,k) = u(nlm ,k)
                u(nlum,k) = u(nfum,k)
             enddo
             !
          enddo
          !
       else
          !
          ! for all layers do
          !
          do k = ks, ke
             !
             npnts = mlu - mf + 1
             !
             if ( INODE == 1 ) then
                !
                ! store data to be sent in array WRKEXG
                !
                do m = mf, mlu
                   nfum = kgrpnt(m,nfu)
                   WRKEXG(m-mf+1) = u(nfum,k)
                enddo
                !
                ! send array WRKEXG
                !
                call SWSENDNB ( WRKEXG, npnts, SWREAL, NPROC, 2 )
                if (STPNOW()) return
                !
                ! receive next array and store in WRKEXG
                !
                call SWRECVNB ( WRKEXG, npnts, SWREAL, NPROC, 2 )
                if (STPNOW()) return
                !
                ! store the received data
                !
                do m = mf, mlu
                   nfm = kgrpnt(m,nf)
                   u(nfm,k) = WRKEXG(m-mf+1)
                enddo
                !
             else if ( INODE == NPROC ) then
                !
                ! store data to be sent in array WRKEXG
                !
                do m = mf, mlu
                   nlm = kgrpnt(m,nl)
                   WRKEXG(m-mf+1) = u(nlm,k)
                enddo
                !
                ! send array WRKEXG
                !
                call SWSENDNB ( WRKEXG, npnts, SWREAL, 1, 2 )
                if (STPNOW()) return
                !
                ! receive next array and store in WRKEXG
                !
                call SWRECVNB ( WRKEXG, npnts, SWREAL, 1, 2 )
                if (STPNOW()) return
                !
                ! store the received data
                !
                do m = mf, mlu
                   nlum = kgrpnt(m,nlu)
                   u(nlum,k) = WRKEXG(m-mf+1)
                enddo
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
end subroutine periodic
!
subroutine periodici ( u, kgrpnt, ks, ke )
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
!    1.00, August 2014: New subroutine
!
!   Purpose
!
!   Fills virtual unknowns of type integer using periodicity
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata
    use m_parall
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc,myc),     intent(in)    :: kgrpnt ! index table containing the address of each (active) grid point
                                                             ! =1; not active grid point
                                                             ! >1; active grid point
    integer,                         intent(in)    :: ke     ! last index in vertical direction
    integer,                         intent(in)    :: ks     ! first index in vertical direction
    integer, dimension(mcgrd,ks:ke), intent(inout) :: u      ! unknown of type integer
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter in vertical direction
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    integer       :: nfm      ! pointer to m,nf
    integer       :: nfum     ! pointer to m,nfu
    integer       :: nlm      ! pointer to m,nl
    integer       :: nlum     ! pointer to m,nlu
    integer       :: nmf      ! pointer to mf,n
    integer       :: nmfu     ! pointer to mfu,n
    integer       :: nml      ! pointer to ml,n
    integer       :: nmlu     ! pointer to mlu,n
    integer       :: npnts    ! number of points to be send/receive
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'periodici')
    !
    ! if no periodicity, return
    !
    if ( .not.lreptx .and. .not.lrepty ) return
    !
    if ( lreptx ) then
       !
       if ( kpart /= 4 ) then
          !
          do n = nf, nlu
             !
             nmf  = kgrpnt(mf ,n)
             nmfu = kgrpnt(mfu,n)
             nml  = kgrpnt(ml ,n)
             nmlu = kgrpnt(mlu,n)
             !
             do k = ks, ke
                u(nmf ,k) = u(nml ,k)
                u(nmlu,k) = u(nmfu,k)
             enddo
             !
          enddo
          !
       else
          !
          ! for all layers do
          !
          do k = ks, ke
             !
             npnts = nlu - nf + 1
             !
             if ( INODE == 1 ) then
                !
                ! store data to be sent in array IWRKEX
                !
                do n = nf, nlu
                   nmfu = kgrpnt(mfu,n)
                   IWRKEX(n-nf+1) = u(nmfu,k)
                enddo
                !
                ! send array IWRKEX
                !
                call SWSENDNB ( IWRKEX, npnts, SWINT, NPROC, 2 )
                if (STPNOW()) return
                !
                ! receive next array and store in IWRKEX
                !
                call SWRECVNB ( IWRKEX, npnts, SWINT, NPROC, 2 )
                if (STPNOW()) return
                !
                ! store the received data
                !
                do n = nf, nlu
                   nmf = kgrpnt(mf,n)
                   u(nmf,k) = IWRKEX(n-nf+1)
                enddo
                !
             else if ( INODE == NPROC ) then
                !
                ! store data to be sent in array IWRKEX
                !
                do n = nf, nlu
                   nml = kgrpnt(ml,n)
                   IWRKEX(n-nf+1) = u(nml,k)
                enddo
                !
                ! send array IWRKEX
                !
                call SWSENDNB ( IWRKEX, npnts, SWINT, 1, 2 )
                if (STPNOW()) return
                !
                ! receive next array and store in IWRKEX
                !
                call SWRECVNB ( IWRKEX, npnts, SWINT, 1, 2 )
                if (STPNOW()) return
                !
                ! store the received data
                !
                do n = nf, nlu
                   nmlu = kgrpnt(mlu,n)
                   u(nmlu,k) = IWRKEX(n-nf+1)
                enddo
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
    if ( lrepty ) then
       !
       if ( kpart /= 3 ) then
          !
          do m = mf, mlu
             !
             nfm  = kgrpnt(m,nf )
             nfum = kgrpnt(m,nfu)
             nlm  = kgrpnt(m,nl )
             nlum = kgrpnt(m,nlu)
             !
             do k = ks, ke
                u(nfm ,k) = u(nlm ,k)
                u(nlum,k) = u(nfum,k)
             enddo
             !
          enddo
          !
       else
          !
          ! for all layers do
          !
          do k = ks, ke
             !
             npnts = mlu - mf + 1
             !
             if ( INODE == 1 ) then
                !
                ! store data to be sent in array IWRKEX
                !
                do m = mf, mlu
                   nfum = kgrpnt(m,nfu)
                   IWRKEX(m-mf+1) = u(nfum,k)
                enddo
                !
                ! send array IWRKEX
                !
                call SWSENDNB ( IWRKEX, npnts, SWINT, NPROC, 2 )
                if (STPNOW()) return
                !
                ! receive next array and store in IWRKEX
                !
                call SWRECVNB ( IWRKEX, npnts, SWINT, NPROC, 2 )
                if (STPNOW()) return
                !
                ! store the received data
                !
                do m = mf, mlu
                   nfm = kgrpnt(m,nf)
                   u(nfm,k) = IWRKEX(m-mf+1)
                enddo
                !
             else if ( INODE == NPROC ) then
                !
                ! store data to be sent in array IWRKEX
                !
                do m = mf, mlu
                   nlm = kgrpnt(m,nl)
                   IWRKEX(m-mf+1) = u(nlm,k)
                enddo
                !
                ! send array IWRKEX
                !
                call SWSENDNB ( IWRKEX, npnts, SWINT, 1, 2 )
                if (STPNOW()) return
                !
                ! receive next array and store in IWRKEX
                !
                call SWRECVNB ( IWRKEX, npnts, SWINT, 1, 2 )
                if (STPNOW()) return
                !
                ! store the received data
                !
                do m = mf, mlu
                   nlum = kgrpnt(m,nlu)
                   u(nlum,k) = IWRKEX(m-mf+1)
                enddo
                !
             endif
             !
          enddo
          !
       endif
       !
    endif
    !
end subroutine periodici
