subroutine SwashBCspectrum ( bcfour, spparm )
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
!    1.00, September 2010: New subroutine
!
!   Purpose
!
!   Computes Fourier components based on energy density spectrum for synthesizing time series along open boundaries
!
!   Method
!
!   In case of multi-directional wave field it is assumed that each wave component has a unique direction of propagation.
!   Based on Miles (1989), wave directions are selected at random from the cumulative distribution function (cdf) of
!   the directional spreading function and are assigned to each frequency component.
!
!   M.D. Miles
!   A note on directional random wave synthesis by the single-summation method
!   Proceedings XXIII IAHR Congress, Ottawa, Canada, C243-C250, 1989
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
    use m_bndspec
    use m_parall
!
    implicit none
!
!   Argument variables
!
    real, dimension(5), intent(in)  :: spparm   ! parameters used for computation of incident spectrum and subsequent Fourier series
                                                ! =1; wave height (significant or rms)
                                                ! =2; wave period (peak or mean)
                                                ! =3; incident or peak wave direction w.r.t. normal on the boundary
                                                ! =4; directional distribution coefficient
                                                ! =5; cyclic period of time series
    !
    type(bfsdat), intent(inout)     :: bcfour   ! list containing parameters for Fourier series
!
!   Parameter variables
!
    integer, parameter              :: maxnit = 10 ! maximum number of iterations
    integer, parameter              :: ndir = 51   ! number of directions
    !
    real   , parameter              :: eps = 0.01  ! convergence criterion
!
!   Local variables
!
    integer                         :: i        ! loop counter
    integer, save                   :: ient = 0 ! number of entries in this subroutine
    integer                         :: j        ! loop counter
    integer                         :: nfmax    ! maximum number of frequencies possible
    integer                         :: nfp      ! nearest integer of peak frequency
    integer                         :: nfreq    ! number of frequencies
    integer                         :: nit      ! number of iterations
    integer                         :: shape    ! spectral shape
                                                ! = 1; Pierson Moskowitz
                                                ! = 2; Jonswap
                                                ! = 3; TMA
    !
    real, dimension(ndir)           :: cdf      ! cumulative distribution function of directional spreading
    real                            :: cosdir   ! cosine of directions
    real                            :: df       ! increment in frequency space
    real, dimension(ndir)           :: dspr     ! directional spreading function
    real                            :: dwidth   ! directional width
    real                            :: f        ! frequency
    real                            :: fac      ! a factor
    real                            :: fmax     ! maximum frequency
    real                            :: fmin     ! minimum frequency
    real                            :: fnyq     ! Nyquist frequency
    real                            :: fp       ! peak frequency
    real                            :: hs       ! significant wave height
    real                            :: m0       ! moment of zeroth order
    real                            :: m1       ! moment of first order
    real                            :: ms       ! spreading coefficient for cosine-power
    real                            :: rval     ! auxiliary real
    real                            :: sector   ! directional sector to be considered
    real                            :: tcycl    ! cyclic period of time record
    real, dimension(ndir)           :: theta    ! spectral directions
    real                            :: tm       ! mean wave period
    real                            :: tp       ! peak wave period
    !
    character(80)                   :: msgstr   ! string to pass message
    !
    logical                         :: peak     ! indicates whether peak (.true.) or mean (.false.) period is given
    logical                         :: STPNOW   ! indicates that program must stop
    !
    real, dimension(:), allocatable :: spec     ! nondimensional relative spectral density, equal to one at the peak
    real, dimension(:), allocatable :: tarr     ! temporary array to store random values
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashBCspectrum')
    !
    if ( spshape(2) < 0 ) then
       shape = -spshape(2)
       peak  = .false.
    else
       shape = spshape(2)
       peak  = .true.
    endif
    !
    nit = 0
    !
    ! determime significant wave height
    !
    if ( spshape(1) == 1 ) then
       hs = sqrt(2.) * spparm(1)
    else if ( spshape(1) == 2 ) then
       hs = spparm(1)
    endif
    !
    ! determine cyclic period of time series
    !
    if ( spparm(5) < 9.e9 ) then
       tcycl = spparm(5)
    else
       tcycl = tfinc - tinic + dt   ! add dt to make it cyclical for the whole simulation
    endif
    !
    ! determine frequency step
    !
    df = 1./tcycl
    !
    fnyq  = 0.5/dt
    nfmax = nint(fnyq/df) + 1
    !
    allocate(spec(nfmax))
    spec = 0.
    !
    ! determine peak frequency
    !
    tp  = spparm(2)
 10 fp  = 1./tp
    nfp = nint(fp/df)
    fp  = nfp*df
    !
    ! determine frequency range and number of frequencies (uniformly distributed)
    !
    fmin  = (nfp/2)*df
    fmax  = min(fnyq,3.*fp)
    nfreq = nint((fmax-fmin)/df)
    !
    ! compute energy spectrum
    !
    if ( shape <= 3 ) then
       call jonswap ( spec, fmin, df, nfreq, fp, gamma )
    else
       call msgerr ( 4, 'unknown frequency shape for incident spectrum' )
       return
    endif
    !
    ! calculate moment of zeroth order
    !
    m0 = df*sum(spec(1:nfreq))
    !
    ! if mean frequency is given recalculate peak period and restart
    !
    if ( .not.peak .and. nit < maxnit ) then
       !
       nit = nit + 1
       !
       m1 = 0.
       f = fmin
       do j = 1, nfreq
          m1 = m1 + f * spec(j)
          f = f + df
       enddo
       m1 = df * m1
       !
       ! calculate mean period
       !
       if ( m1 > 0. ) then
          tm = m0 / m1
       else
          call msgerr (3, 'not possible to compute mean period for incident spectrum')
       endif
       if ( ITEST >= 80 ) write (PRTEST,202) nit, spparm(2), tm, tp
       !
       if ( abs(tm-spparm(2)) > eps*spparm(2) ) then
          tp = ( spparm(2)/tm ) * tp
          goto 10
       endif
       !
    endif
    !
    if ( nit >= maxnit ) then
       call msgerr (1, 'no convergence in computing peak period for incident spectrum')
       write(PRINTF,203) tp
    else if ( .not.peak ) then
       write (msgstr, '(a,f7.2,a,f7.2,a)') ' mean wave period is ',tm,' sec, while peak wave period is ',tp,' sec'
       call msgerr( 1, trim(msgstr) )
    endif
    !
    fac  = hs*hs/16./m0
    spec = fac*spec
    !
    if ( ITEST >= 30 ) then
       m0 = df*sum(spec(1:nfreq))
       rval = 4.*sqrt(m0)
       if ( abs(rval - hs) > 0.1*hs ) write (PRINTF, 201) hs, rval
    endif
    !
    bcfour%nfreq = nfreq
    !
    if ( nfreq > 0 ) then
       !
       allocate(bcfour%ampl (nfreq))
       allocate(bcfour%omega(nfreq))
       allocate(bcfour%phase(nfreq))
       allocate(bcfour%theta(nfreq))
       !
       ! generate random values for phase
       allocate(tarr(nfreq))
       if ( INODE == MASTER ) then
          call random_seed(put=seedf)
          do j = 1, nfreq
             call random_number(rval)
             tarr(j) = rval
          enddo
       endif
       !
       ! scatter array tarr to all nodes
       !
       call SWBROADC ( tarr, nfreq, SWREAL )
       if (STPNOW()) return
       !
       f = fmin
       do j = 1, nfreq
          bcfour%ampl (j) = sqrt(2.*df*spec(j))
          bcfour%omega(j) = 2.*pi*f
          bcfour%phase(j) = 2.*pi*tarr(j)
          bcfour%theta(j) = 0.
          f = f + df
       enddo
       !
       deallocate(tarr)
       !
    endif
    !
    deallocate(spec)
    !
    ! include distribution over directions, if appropriate
    !
    if ( spparm(4) /= 0. ) then
       !
       ! calculate spreading coefficient
       !
       if ( spshape(3) == 1 ) then
          dwidth = spparm(4)
          if ( dwidth > 23. ) then
             fac = 1.2
          elseif ( dwidth > 17. ) then
             fac = 1.096
          else
             fac = 1.01
          endif
          ms = max(fac*(pi*dwidth/180.)**(-2) - 2., 1.)
       else
          ms = spparm(4)
       endif
       !
       ! determine sector of directional space depending on spreading
       !
       if ( ms <= 6. ) then
          sector = pi
       else if ( ms <= 30. ) then
          sector = 0.5*pi
       else if ( ms <= 100. ) then
          sector = 0.25*pi
       else
          sector = 0.125*pi
       endif
       !
       ! determine directional space
       !
       do i = 1, ndir
          !
          theta(i) = sector*(real(i-1)/real(ndir-1) - 0.5)
          !
       enddo
       !
       ! compute directional spreading (based on cos m-model) and its distribution function
       !
       do i = 1, ndir
          cosdir = cos(theta(i))
          if ( cosdir > 0. ) then
             dspr(i) = max(cosdir**ms, 1.e-10)
          else
             dspr(i) = 0.
          endif
          cdf(i) = sum(dspr(1:i))
       enddo
       cdf = cdf/sum(dspr)
       !
       ! pick a random number in distribution function and assign corresponding (interpolated) direction to each frequency
       !
       allocate(tarr(nfreq))
       if ( INODE == MASTER ) then
          call random_seed(put=seedt)
          do j = 1, nfreq
             call random_number(rval)
             tarr(j) = rval
          enddo
       endif
       call SWBROADC ( tarr, nfreq, SWREAL )
       if (STPNOW()) return
       !
       do j = 1, nfreq
          rval = tarr(j)
          do i = 2, ndir
             if ( rval > cdf(i-1) .and. rval <= cdf(i) ) then
                bcfour%theta(j) = theta(i-1) + (theta(i)-theta(i-1))*(rval-cdf(i-1))/(cdf(i)-cdf(i-1))
             endif
          enddo
       enddo
       !
       deallocate(tarr)
       !
    endif
    !
 201 format (' deviation in Hs of incident spectrum, should be ', f8.3, ', calculated ', f8.3)
 202 format (' computation incident spectrum: iter=', i2, '  period values:', 3f7.2)
 203 format (' peak wave period is now ',f7.2,' sec')
    !
end subroutine SwashBCspectrum
