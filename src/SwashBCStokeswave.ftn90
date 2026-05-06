subroutine SwashBCStokeswave ( bcfour, nfreq, ibgrpt, swd, istok )
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
!    1.00: Panagiotis Vasarmidis
!
!   Updates
!
!    1.00, October 2023: New subroutine
!
!   Purpose
!
!   Computes first and second order Stokes wave components for synthesizing time series along open boundaries
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_bndspec
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)         :: ibgrpt ! actual boundary grid point
    integer, intent(in)         :: istok  ! indicates the order of Stokes wave theory
                                          ! = 0; use hyperbolic cosine distribution for velocity (Airy wave theory)
                                          ! = 1; first order Stokes wave (Airy wave theory)
                                          ! = 2; second order Stokes wave
                                          ! = 3; second order sub- and super-harmonic transfer functions
    integer, intent(in)         :: nfreq  ! number of frequencies
    !
    real, intent(in)            :: swd    ! still water depth
    !
    type(bfsdat), intent(inout) :: bcfour ! list containing parameters for Fourier series
!
!   Local variables
!
    integer, save  :: ient  = 0           ! number of entries in this subroutine
    integer        :: j                   ! loop counter
    !
    real           :: ampl                ! amplitude of a Fourier component
    real           :: cph                 ! cosine of phase
    real           :: c2ph                ! cosine of two times phase
    real           :: dn                  ! denominator
    real           :: fac                 ! multiplication factor
    real           :: kd                  ! dimensionless still water depth
    real           :: kwav                ! wave number of a Fourier component
    real           :: omeg                ! angular frequency of a Fourier component
    real           :: sph                 ! sine of phase
    real           :: s2ph                ! sine of two times phase
    real           :: wn                  ! last angular frequency
    !
    logical        :: spectrum            ! indicates whether wave spectrum is dealt with or not
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashBCStokeswave')
    !
    ! is Fourier series a spectrum?
    !
    if ( istok == 3 .and. nfreq > 2 ) then
       spectrum = .true.
    else
       spectrum = .false.
    endif
    !
    wn = bcfour%omega(nfreq)
    !
    if ( kmax == 1 ) then
       !
       do j = 1, nfreq
          !
          ampl = bcfour%ampl(j)
          kwav = kwave(ibgrpt,j)
          !
          if ( .not. ampl /= 0. ) cycle
          !
          kd = kwav * swd
          !
          cph = comp1(ibgrpt,j) / ampl
          sph = comp2(ibgrpt,j) / ampl
          !
          ! first order solution surface
          !
          stkz1c(ibgrpt,j) = ampl * cph
          stkz1s(ibgrpt,j) = ampl * sph
          !
          ! first order solution velocity
          !
          fac = 2. * ampl * grav / sqrt( (4.+kd*kd)*grav*swd )
          !
          stku1c(ibgrpt,j,1) = fac * cph
          stku1s(ibgrpt,j,1) = fac * sph
          !
       enddo
       !
       if ( istok > 1 ) then
          !
          do j = 1, nfreq
             !
             ampl = bcfour%ampl (j)
             omeg = bcfour%omega(j)
             kwav = kwave(ibgrpt,j)
             !
             if ( .not. ampl /= 0. ) cycle
             if ( spectrum .and. 2.*omeg > wn ) exit
             !
             kd = kwav * swd
             !
             cph = comp1(ibgrpt,j) / ampl
             sph = comp2(ibgrpt,j) / ampl
             !
             c2ph = 2.*cph*cph - 1.
             s2ph = 2.*sph*cph
             !
             ! second order solution surface
             !
             fac = ampl * ampl * ( 4. + kd*kd ) / ( 4.*kd*kd*swd )
             !
             stkz2c(ibgrpt,j) = fac * c2ph
             stkz2s(ibgrpt,j) = fac * s2ph
             !
             ! second order solution velocity
             !
             fac = -( kd*kd - 4. ) * ampl * ampl * grav * swd*swd / ( 2.*kd*kd * sqrt( (4.+kd*kd)*grav*swd**7) )
             !
             stku2c(ibgrpt,j,1) = fac * c2ph
             stku2s(ibgrpt,j,1) = fac * s2ph
             !
          enddo
          !
       endif
       !
    else if ( kmax == 2 ) then
       !
       do j = 1, nfreq
          !
          ampl = bcfour%ampl(j)
          kwav = kwave(ibgrpt,j)
          !
          if ( .not. ampl /= 0. ) cycle
          !
          kd = kwav * swd
          !
          cph = comp1(ibgrpt,j) / ampl
          sph = comp2(ibgrpt,j) / ampl
          !
          ! first order solution surface
          !
          stkz1c(ibgrpt,j) = ampl * cph
          stkz1s(ibgrpt,j) = ampl * sph
          !
          ! first order solution velocity
          !
          dn = sqrt(256. + 96.*kd*kd + kd**4) * sqrt( (kd*kd*(16. + kd*kd)*grav)/swd ) * swd
          !
          fac = 4. * kd * ( 16. + 3.*kd*kd) * ampl * grav / dn
          !
          stku1c(ibgrpt,j,1) = fac * cph
          stku1s(ibgrpt,j,1) = fac * sph
          !
          fac = -4. * kd * (-16. + kd*kd) * ampl * grav / dn
          !
          stku1c(ibgrpt,j,2) = fac * cph
          stku1s(ibgrpt,j,2) = fac * sph
          !
       enddo
       !
       if ( istok > 1 ) then
          !
          do j = 1, nfreq
             !
             ampl = bcfour%ampl (j)
             omeg = bcfour%omega(j)
             kwav = kwave(ibgrpt,j)
             !
             if ( .not. ampl /= 0. ) cycle
             if ( spectrum .and. 2.*omeg > wn ) exit
             !
             kd = kwav * swd
             !
             cph = comp1(ibgrpt,j) / ampl
             sph = comp2(ibgrpt,j) / ampl
             !
             c2ph = 2.*cph*cph - 1.
             s2ph = 2.*sph*cph
             !
             ! second order solution surface
             !
             fac = ampl * ampl * (49152. + 33792.*kd*kd + 9792.*kd**4 + 892.*kd**6 + 7.*kd**8) / ( 12.*kd*kd * (16. + kd*kd) * (320. + 20.*kd*kd + kd**4) * swd )
             !
             stkz2c(ibgrpt,j) = fac * c2ph
             stkz2s(ibgrpt,j) = fac * s2ph
             !
             ! second order solution velocity
             !
             dn = 6.*kd * (320. + 20.*kd*kd + kd**4) * sqrt(256. + 96.*kd*kd + kd**4) * sqrt( (kd*kd*(16. + kd*kd)*grav)/swd ) *swd*swd
             !
             fac = ( 98304. + 55296.*kd*kd + 18048.*kd**4 + 3608.*kd**6 - 33.*kd**8 ) * ampl * ampl * grav / dn
             !
             stku2c(ibgrpt,j,1) = fac * c2ph
             stku2s(ibgrpt,j,1) = fac * s2ph
             !
             fac = ( 98304. - 43008.*kd*kd +  5760.*kd**4 -  904.*kd**6 + 37.*kd**8 ) * ampl * ampl * grav / dn
             !
             stku2c(ibgrpt,j,2) = fac * c2ph
             stku2s(ibgrpt,j,2) = fac * s2ph
             !
          enddo
          !
       endif
       !
    else if ( kmax == 3 ) then
       !
       do j = 1, nfreq
          !
          ampl = bcfour%ampl(j)
          kwav = kwave(ibgrpt,j)
          !
          if ( .not. ampl /= 0. ) cycle
          !
          kd = kwav * swd
          !
          cph = comp1(ibgrpt,j) / ampl
          sph = comp2(ibgrpt,j) / ampl
          !
          ! first order solution surface
          !
          stkz1c(ibgrpt,j) = ampl * cph
          stkz1s(ibgrpt,j) = ampl * sph
          !
          ! first order solution velocity
          !
          dn = sqrt( (36. + kd*kd) * (1296. + 504.*kd*kd + kd**4) ) * sqrt( (kd*kd*(12. + kd*kd)*(108. + kd*kd)*grav) / swd ) * swd
          !
          fac = 6. * kd * ( 1296. + 360.*kd*kd + 5.*kd**4 ) * ampl * grav / dn
          !
          stku1c(ibgrpt,j,1) = fac * cph
          stku1s(ibgrpt,j,1) = fac * sph
          !
          fac = -18. * ( -6. + kd ) * kd * ( 6. + kd ) * ( 12. + kd*kd ) * ampl * grav / dn
          !
          stku1c(ibgrpt,j,2) = fac * cph
          stku1s(ibgrpt,j,2) = fac * sph
          !
          fac = 6. * kd * ( -36. + kd*kd )**2 * ampl * grav / dn
          !
          stku1c(ibgrpt,j,3) = fac * cph
          stku1s(ibgrpt,j,3) = fac * sph
          !
       enddo
       !
       if ( istok > 1 ) then
          !
          do j = 1, nfreq
             !
             ampl = bcfour%ampl (j)
             omeg = bcfour%omega(j)
             kwav = kwave(ibgrpt,j)
             !
             if ( .not. ampl /= 0. ) cycle
             if ( spectrum .and. 2.*omeg > wn ) exit
             !
             kd = kwav * swd
             !
             cph = comp1(ibgrpt,j) / ampl
             sph = comp2(ibgrpt,j) / ampl
             !
             c2ph = 2.*cph*cph - 1.
             s2ph = 2.*sph*cph
             !
             ! second order solution surface
             !
             fac = ampl * ampl * ( 44079842304. + 38773935360.*kd*kd + 14751227520.*kd**4 + 2589058080.*kd**6 + 151553268.*kd**8 + 2932065.*kd**10 + 16482.*kd**12 + 25.*kd**14 ) / ( 36.*kd*kd* (1587237120. + 411505920.*kd*kd + 42620256.*kd**4 + 1963440.*kd**6 + 32337.*kd**8 + 270.*kd**10 + kd**12) * swd )
             !
             stkz2c(ibgrpt,j) = fac * c2ph
             stkz2s(ibgrpt,j) = fac * s2ph
             !
             ! second order solution velocity
             !
             dn = 6.*sqrt((36. + kd*kd)*(1296. + 504.*kd*kd + kd**4)) * (1224720. + 204120.*kd*kd + 13041.*kd**4 + 150.*kd**6 + kd**8) * ((kd*kd*(12. + kd*kd)*(108. + kd*kd)*grav)/swd)**1.5 * swd**3
             !
             fac = kd * (57127475625984. + 60830182379520.*kd*kd + 25362235192320.*kd**4 + 7842357211392.*kd**6 + 1009720893888.*kd**8 + 54844816176.*kd**10 + 1259413596.*kd**12 + 10131075.*kd**14 + 5829.*kd**16 - 124.*kd**18) * ampl * ampl * grav * grav / dn
             !
             stku2c(ibgrpt,j,1) = fac * c2ph
             stku2s(ibgrpt,j,1) = fac * s2ph
             !
             fac = 3.*kd * (19042491875328. + 3350068015104.*kd*kd + 1871488613376.*kd**4 + 744504908544.*kd**6 + 46933929792.*kd**8 - 381541104.*kd**10 - 35910540.*kd**12 + 214209.*kd**14 + 13971.*kd**16 + 68.*kd**18) * ampl * ampl * grav * grav / dn
             !
             stku2c(ibgrpt,j,2) = fac * c2ph
             stku2s(ibgrpt,j,2) = fac * s2ph
             !
             fac = kd * (57127475625984. - 15339785121792.*kd*kd + 1382800978944.*kd**4 - 317855570688.*kd**6 - 20778949440.*kd**8 + 2738368944.*kd**10 + 88138044.*kd**12 + 487431.*kd**14 - 10356.*kd**16 - 59.*kd**18) * ampl * ampl * grav * grav / dn
             !
             stku2c(ibgrpt,j,3) = fac * c2ph
             stku2s(ibgrpt,j,3) = fac * s2ph
             !
          enddo
          !
       endif
       !
    else if ( kmax == 4 ) then
       !
       do j = 1, nfreq
          !
          ampl = bcfour%ampl(j)
          kwav = kwave(ibgrpt,j)
          !
          if ( .not. ampl /= 0. ) cycle
          !
          kd = kwav * swd
          !
          cph = comp1(ibgrpt,j) / ampl
          sph = comp2(ibgrpt,j) / ampl
          !
          ! first order solution surface
          !
          stkz1c(ibgrpt,j) = ampl * cph
          stkz1s(ibgrpt,j) = ampl * sph
          !
          ! first order solution velocity
          !
          dn = sqrt( 16777216. + 7340032.*kd*kd + 286720.*kd**4 + 1792.*kd**6 + kd**8 ) * sqrt( (kd*kd*(64. + kd*kd)*(4096. + 384.*kd*kd + kd**4)*grav)/swd ) * swd
          !
          fac = 8. * kd * (262144. + 86016.*kd*kd + 2240.*kd**4 + 7.*kd**6) * ampl * grav / dn
          !
          stku1c(ibgrpt,j,1) = fac * cph
          stku1s(ibgrpt,j,1) = fac * sph
          !
          fac = -8. * (-8. + kd) * kd * (8. + kd) * (4096. + 640.*kd*kd + 5.*kd**4) * ampl * grav / dn
          !
          stku1c(ibgrpt,j,2) = fac * cph
          stku1s(ibgrpt,j,2) = fac * sph
          !
          fac = 8. * kd * (-64. + kd*kd)**2 * (64. + 3.*kd*kd) * ampl * grav / dn
          !
          stku1c(ibgrpt,j,3) = fac * cph
          stku1s(ibgrpt,j,3) = fac * sph
          !
          fac = -8. * kd * (-64. + kd*kd)**3 * ampl * grav / dn
          !
          stku1c(ibgrpt,j,4) = fac * cph
          stku1s(ibgrpt,j,4) = fac * sph
          !
       enddo
       !
       if ( istok > 1 ) then
          !
          do j = 1, nfreq
             !
             ampl = bcfour%ampl (j)
             omeg = bcfour%omega(j)
             kwav = kwave(ibgrpt,j)
             !
             if ( .not. ampl /= 0. ) cycle
             if ( spectrum .and. 2.*omeg > wn ) exit
             !
             kd = kwav * swd
             !
             cph = comp1(ibgrpt,j) / ampl
             sph = comp2(ibgrpt,j) / ampl
             !
             c2ph = 2.*cph*cph - 1.
             s2ph = 2.*sph*cph
             !
             ! second order solution surface
             !
             fac = ampl * ampl * ( 54043195528445952. + 53198770598313984.*kd*kd + 22935812555407360.*kd**4 + 5028960027017216.*kd**6 + 466739464765440.*kd**8 + 20187671691264.*kd**10 + 420926980096.*kd**12 + 3866914816.*kd**14 + 12686592.*kd**16 + 15504.*kd**18 + 7.*kd**20 ) / ( 12.*kd*kd * (5910974510923776. + 2031897488130048.*kd*kd + 281406257233920.*kd**4 + 19189913878528.*kd**6 + 658254069760.*kd**8 + 11703681024.*kd**10 + 105091072.*kd**12 + 446208.*kd**14 + 1008.*kd**16 + kd**18) * swd )
             !
             stkz2c(ibgrpt,j) = fac * c2ph
             stkz2s(ibgrpt,j) = fac * s2ph
             !
             ! second order solution velocity
             !
             dn = 3.*sqrt(16777216. + 7340032.*kd*kd + 286720.*kd**4 + 1792.*kd**6 + kd**8) * (22548578304. + 5284823040.*kd*kd + 456916992.*kd**4 + 14110720.*kd**6 + 166656.*kd**8 + 560.*kd**10 + kd**12) * ( (kd*kd*(64. + kd*kd)*(4096. + 384.*kd*kd + kd**4)*grav) / swd )**1.5 * swd**3
             !
             fac = kd * (64. + kd*kd) * (442721857769029238784. + 574154909294209794048.*kd*kd + 271224783958760751104.*kd**4 + 91730443110189105152.*kd**6 + 15671717462691282944.*kd**8 + 1315611292366536704.*kd**10 + 56798533012946944.*kd**12 + 1249202982092800.*kd**14 + 12676835573760.*kd**16 + 46869209088.*kd**18 + 31214848.*kd**20 - 139648.*kd**22 - 187.*kd**24) * ampl * ampl * grav * grav / dn
             !
             stku2c(ibgrpt,j,1) = fac * c2ph
             stku2s(ibgrpt,j,1) = fac * s2ph
             !
             fac = kd * (64. + kd*kd) * (442721857769029238784. + 242113515967437864960.*kd*kd + 98286558267733704704.*kd**4 + 41733731846935543808.*kd**6 + 5539005329200644096.*kd**8 + 294519582662590464.*kd**10 + 6540981788737536.*kd**12 + 28983311728640.*kd**14 - 362364272640.*kd**16 + 9686278144.*kd**18 + 107953920.*kd**20 + 331968.*kd**22 + 311.*kd**24) * ampl * ampl * grav * grav / dn
             !
             stku2c(ibgrpt,j,2) = fac * c2ph
             stku2s(ibgrpt,j,2) = fac * s2ph
             !
             fac = kd * (64. + kd*kd) * (442721857769029238784. + 20752587082923245568.*kd*kd + 29111267991322886144.*kd**4 + 12437816270890467328.*kd**6 + 1071082655128223744.*kd**8 + 24421252764532736.*kd**10 + 128063039864832.*kd**12 + 37424264642560.*kd**14 + 1296101277696.*kd**16 + 2277367808.*kd**18 - 41250048.*kd**20 - 170560.*kd**22 - 189.*kd**24) * ampl * ampl * grav * grav / dn
             !
             stku2c(ibgrpt,j,3) = fac * c2ph
             stku2s(ibgrpt,j,3) = fac * s2ph
             !
             fac = kd * (64. + kd*kd) * (442721857769029238784. - 89927877359334064128.*kd*kd + 8358680908399640576.*kd**4 - 3291005427700989952.*kd**6 - 412536762741555200.*kd**8 + 4378805057617920.*kd**10 + 720708397170688.*kd**12 - 2290626854912.*kd**14 - 534823567360.*kd**16 - 1883348992.*kd**18 + 11237632.*kd**20 + 56960.*kd**22 + 73.*kd**24) * ampl * ampl * grav * grav / dn
             !
             stku2c(ibgrpt,j,4) = fac * c2ph
             stku2s(ibgrpt,j,4) = fac * s2ph
             !
          enddo
          !
       endif
       !
    endif
    !
end subroutine SwashBCStokeswave
