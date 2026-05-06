subroutine SwashPrintSettings
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
!    1.00, March 2010: New subroutine
!
!   Purpose
!
!   Prints all the settings used in SWASH run
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use SwashRigBoddata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    character(18) :: DTTIWR   ! to write a time string
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashPrintSettings')
    !
    write (PRINTF, 201) 'SWASH'
    !
    if ( oned ) then
       write (PRINTF, *) 'One-dimensional mode of SWASH is activated'
    endif
    !
    if ( linswe ) then
       write (PRINTF, *) 'The linearized shallow water equations are employed'
    endif
    !
    if ( .not.momskip ) then
       !
       if ( ihydro /= 0 ) then
          write (PRINTF, *) 'Non-hydrostatic mode of SWASH is activated'
       endif
       !
       if ( kmax == 1 ) then
          write (PRINTF, *) 'Depth-averaged mode of SWASH is activated'
       else
          if ( .not.lsubg ) then
             write (PRINTF,'(a,i3,a)') ' Water depth is divided into ',kmax,' layers'
          else
             write (PRINTF,'(a,i3,a)') ' Water depth is divided into ',kmax ,' velocity layers'
             write (PRINTF,'(a,i3,a)') ' Water depth is divided into ',kpmax,' pressure layers'
          endif
       endif
       !
       if ( qlay /= 0 ) then
          if ( qmax == 1 ) then
             write (PRINTF,'(a)') ' Reduced pressure equation method is employed: 1 pressure layer'
          else if ( qmax > 1 ) then
             write (PRINTF,'(a,i3,a)') ' Reduced pressure equation method is employed: ',qmax,' pressure layers'
          endif
       endif
       !
       if ( numdisp ) then
          write (PRINTF,*) 'Approximate dispersion relation is employed'
       endif
       !
       if ( .not.mimetic ) then
          if ( strictmom ) then
             write (PRINTF,*) 'Advection approximation is strictly momentum conservative'
          endif
          if ( stricthead ) then
             write (PRINTF,*) 'Advection approximation is strictly head energy conservative'
          endif
       else
          write (PRINTF,*) 'Mimetic discretizations are employed'
       endif
       !
    else
       !
       write (PRINTF,*) 'Momentum and continuity equations are skipped'
       !
    endif
    !
    if ( bnaut ) then
       write (PRINTF, 202) 'nautical '
    else
       write (PRINTF, 202) 'Cartesian'
    endif
    !
    if ( kspher == 0 ) then
       write (PRINTF, *) 'Cartesian coordinate system is used'
    else
       write (PRINTF, *) 'Spherical coordinate system is used'
    endif
    !
    if ( optg /= 5 ) then
       if ( optg == 1 ) then
          write (PRINTF, *) 'Computational grid is rectilinear'
       else
          write (PRINTF, *) 'Computational grid is curvilinear'
       endif
       if ( lreptx ) then
          write (PRINTF, *) 'This (partitioned) grid repeats itself in x-direction'
       else if ( lrepty ) then
          write (PRINTF, *) 'This (partitioned) grid repeats itself in y-direction'
       endif
       write (PRINTF, 203) mxc, myc
       write (PRINTF, 204) mcgrd
       if ( optg == 1 ) then
          if ( oned ) then
             write (PRINTF, 205) dx, 0.
          else
             write (PRINTF, 205) dx, dy
          endif
       endif
    else
       write (PRINTF, *) 'Computational mesh is unstructured'
    endif
    !
    if ( mtimei == 2 ) write (PRINTF, 206) nstatc, mtc
    !
    write (PRINTF, 207) epsdry, dpsopt
    !
    if ( .not.momskip ) then
       !
       write (PRINTF, 208) grav, rhow
       write (PRINTF, 209) rhoa, dynvis
       !
       write (PRINTF, 210) u10, wdic/degrad
       !
       if ( iwind == 1 ) then
          write (PRINTF, 211) pwnd(1)
       else if ( iwind == 2 ) then
          write (PRINTF, 212) pwnd(2), pwnd(3)
       else if ( iwind == 3 ) then
          write (PRINTF, 213) pwnd(5), pwnd(7)
          write (PRINTF, 214) pwnd(8), pwnd(9)
       else if ( iwind == 4 ) then
          write (PRINTF, *) 'Winddrag based on 2nd order polynomial fit'
       else
          write (PRINTF, *) 'Wind stress is off'
       endif
       !
       if ( svwp ) write (PRINTF, *) 'Space varying wind and pressure is included'
       !
       if ( iamb /= 0 ) write (PRINTF, *) 'Ambient current is included'
       !
       if ( irough == 1 ) then
          write (PRINTF, 215) pbot(1)
       else if ( irough == 2 ) then
          write (PRINTF, 216) pbot(1)
       else if ( irough == 3 ) then
          write (PRINTF, 217) pbot(1)
       else if ( irough == 5 ) then
          write (PRINTF, 218) pbot(1)
       else if ( irough == 4 .and. pbot(2) /= 0. ) then
          write (PRINTF, 219) pbot(2)
       else if ( irough == 4 ) then
          write (PRINTF, *) 'Logarithmic wall-law for smooth bed'
       else if ( irough == 11 ) then
          write (PRINTF, 220) pbot(1)
       else
          write (PRINTF, *) 'Bottom friction is off'
       endif
       !
       if ( iveg /= 0 ) then
          write (PRINTF, *) 'Wave damping due to vegetation is included'
       else
          write (PRINTF, *) 'Vegetation is off'
       endif
       !
       if ( isurf /= 0 ) then
          write (PRINTF, 221) psurf(1), psurf(2)
       else
          write (PRINTF, *) 'No wave breaking control'
       endif
       !
       if ( idens /= 0 ) then
          write (PRINTF, 222) salw, tempw
       else
          write (PRINTF, *) 'Baroclinic forcing is off'
       endif
       !
       if ( coriolis ) then
          if ( leds(23) /= 2 ) then
             if ( pcor(1) /= -999. ) then
                write (PRINTF, 254) pcor(1)
             else
                if ( kspher == 0 ) then
                   write (PRINTF, 254) 4.*pi * sin(degrad*ylatc) / ( 24. * 3600. )
                else
                   write (PRINTF, *) 'Coriolis forcing is included'
                endif
             endif
          else
             write (PRINTF, *) 'Coriolis forcing is included'
          endif
       else
          write (PRINTF, *) 'Coriolis forcing is off'
       endif
       !
       if ( itrans /= 0 ) then
          write (PRINTF, 223) hdiff
          if ( lsal > 0 ) then
             write (PRINTF, *) 'Transport of salinity is included'
             write (PRINTF, 224) tcret
          endif
          if ( ltemp > 0 ) write (PRINTF, *) 'Transport of temperature is included'
          if ( lsed > 0 ) then
             if ( lmixt ) then
                write (PRINTF, *) 'Transport of suspended sediment is included'
                if ( psed(10) > 0. ) then
                   write (PRINTF, 225) psed(1), psed(12)
                else if ( psed(2) > 0. ) then
                   write (PRINTF, 226) psed(1), psed(2)
                else
                   write (PRINTF, *) 'No mass exchange between bed and flow'
                endif
             else
                write (PRINTF, *) 'Transport of tracer is included'
             endif
          endif
       else
          write (PRINTF, *) 'No transport of constituent'
       endif
       !
       if ( ihvisc == 1 ) then
          write (PRINTF, 227) hvisc
       else if ( ihvisc == 2 ) then
          write (PRINTF, 228) csmag
       else if ( ihvisc == 3 ) then
          write (PRINTF, 229) lmix
       else if ( iturb < 2 ) then
          write (PRINTF, *) 'Horizontal viscosity is off'
       endif
       !
       if ( iturb == 1 ) then
          write (PRINTF, *) 'Standard k-eps model'
       else if ( iturb == 2 ) then
          write (PRINTF, *) '3D viscosity with linear k-eps model'
       else if ( iturb == 3 ) then
          write (PRINTF, *) '3D viscosity with nonlinear k-eps model'
       else if ( bvisc > 0. ) then
         write (PRINTF, 230) bvisc
       else
          write (PRINTF, *) 'Vertical viscosity is off'
       endif
       !
       if ( iporos == 1 ) then
          write (PRINTF, 231) ppor(3), ppor(4)
       else
          write (PRINTF, *) 'Porosity is off'
       endif
       !
       if ( mtimei == 1 ) then
          write (PRINTF, 232) pnums(2), pnums(3)
       else if ( mtimei == 2 .and. .not.mimetic ) then
          write (PRINTF, 233) pnums(1), pnums(4)
       endif
       !
       if ( inewt /= 0 ) write (PRINTF,*) 'Newton-type iteration is employed'
       !
       if ( ihydro /= 0 ) then
          write (PRINTF, 234) pnums(5)
       endif
       !
       if ( kmax > 1 .and. .not.mimetic ) then
          write (PRINTF, 235) pnums(31)
          if ( ihydro == 1 .or. ihydro == 2 ) write (PRINTF, 236) pnums(32)
          if ( itrans /= 0 ) write (PRINTF, 237) pnums(33)
       endif
       !
       if ( ifloat /= 0 ) then
          write (PRINTF, 238) pship(2)
       endif
       !
       if ( .not.mimetic ) then
          write (PRINTF, 239) nint(pnums(6)), pnums(7)
          write (PRINTF, 246) pnums(8), pnums(9)
       endif
       if ( kmax > 1 .and. .not.mimetic ) then
          write (PRINTF, 240) nint(pnums(36)), pnums(37)
          write (PRINTF, 246) pnums(38), pnums(39)
       endif
       if ( horwinc .and. ihydro /= 0 ) then
          write (PRINTF, 241) nint(pnums(16)), pnums(17)
          write (PRINTF, 246) pnums(18), pnums(19)
       endif
       if ( kmax > 1 .and. (ihydro == 1 .or. ihydro == 2) .and. .not.mimetic ) then
          if ( verwinc ) then
             write (PRINTF, 242) nint(pnums(41)), pnums(42)
             write (PRINTF, 246) pnums(43), pnums(44)
          endif
       endif
       if ( .not.mimetic ) then
          write (PRINTF, 243) nint(pnums(11)), pnums(12)
          write (PRINTF, 246) pnums(13), pnums(14)
       endif
       !
    else
       !
       write (PRINTF, 223) hdiff
       write (PRINTF, *) 'Transport of tracer is included'
       !
    endif
    !
    if ( itrans /= 0 ) then
       write (PRINTF, 244) nint(pnums(46)), pnums(47)
       write (PRINTF, 246) pnums(48), pnums(49)
       if ( kmax > 1 ) then
          write (PRINTF, 245) nint(pnums(51)), pnums(52)
          write (PRINTF, 246) pnums(53), pnums(54)
       endif
    endif
    !
    if ( ihydro /= 0 ) then
       if ( iproj == 1 ) then
          write (PRINTF,*) 'Pressure correction method is employed'
       else if ( iproj == 2 ) then
          write (PRINTF,*) 'Pressure projection method is employed'
          if ( lpproj ) then
             write (PRINTF, 247) pnums(58), nint(pnums(59))
          endif
       endif
    endif
    !
    if ( ifloat == 2 ) then
       if ( mfoti == 1 ) then
          write (PRINTF,*) 'Newmark scheme is employed'
       else if ( mfoti == 2 ) then
          write (PRINTF,*) 'Chung-Hulbert scheme is employed'
       else if ( mfoti == 3 ) then
          write (PRINTF,*) 'Hilber-Hughes-Taylor scheme is employed'
       else if ( mfoti == 4 ) then
          write (PRINTF,*) 'Wood-Bossak-Zienkiewicz scheme is employed'
       endif
       write (PRINTF, 248) alfam   , alfaf
       write (PRINTF, 249) pship(3), pship(4)
       write (PRINTF, 250) pship(8), nint(pship(9))
    endif
    !
    if ( optg /= 5 ) then
       write (PRINTF, 251) spwidl, spwidr
       write (PRINTF, 252) spwidb, spwidt
    endif
    !
    if ( ITMOPT /= 7 ) write (PRINTF, 253) DTTIWR(ITMOPT, 0d0)
    !
 201 format (/,                                                          &
     '----------------------------------------------------------------'  &
     ,/,                                                                 &
     '                  COMPUTATIONAL PART OF ', a                       &
     ,/,                                                                 &
     '----------------------------------------------------------------'  &
     ,/)
 202 format (' The ',a9, ' convention for velocity direction is used')
 203 format (' Gridresolution       : MXC    ',i12  ,' MYC   ',i12)
 204 format ('                      : MCGRD  ',i12)
 205 format (' Mesh sizes           : DX     ',e12.4,' DY    ',e12.4)
 206 format ('                      : NSTATC ',i12  ,' MTC   ',i12)
 207 format (' Drying/flooding      : DEPMIN ',e12.4,' DPSOPT',i12)
 208 format (' Physical constants   : GRAV   ',e12.4,' RHOW  ',e12.4)
 209 format ('                      : RHOA   ',e12.4,' DYNVIS',e12.4)
 210 format (' Wind input           : WSPEED ',e12.4,' DIR   ',e12.4)
 211 format (' Constant wind stress : CD     ',e12.4)
 212 format (' Charnock formulation : BETA   ',e12.4,' HEIGHT',e12.4)
 213 format (' Winddrag formulation : A1     ',e12.4,' B     ',e12.4)
 214 format ('                      : WLOW   ',e12.4,' WHIGH ',e12.4)
 215 format (' Const bottom friction: CF     ',e12.4)
 216 format (' Chezy formulation    : CF     ',e12.4)
 217 format (' Manning formulation  : CF     ',e12.4)
 218 format (' Colebrook-White form : HEIGHT ',e12.4)
 219 format (' Log-law for rough bed: HEIGHT ',e12.4)
 220 format (' Linear bottom frict  : K      ',e12.4)
 221 format (' Wave breaking control: ALPHA  ',e12.4,' BETA  ',e12.4)
 222 format (' Nonuniform density   : SALW   ',e12.4,' TEMPW ',e12.4)
 223 format (' Horz diffusivity     : DIFF   ',e12.4)
 224 format (' Return time          : TCRET  ',e12.4)
 225 format (' Cohesive sediment    : FALLV  ',e12.4,' ERATE ',e12.4)
 226 format (' Noncohesive sediment : FALLV  ',e12.4,' SIZE  ',e12.4)
 227 format (' Const horz viscosity : VISC   ',e12.4)
 228 format (' Smagorinsky model    : CS     ',e12.4)
 229 format (' Mixing length model  : LM     ',e12.4)
 230 format (' Const vert viscosity : VVISC  ',e12.4)
 231 format (' Porosity             : ALPHA0 ',e12.4,' BETA0 ',e12.4)
 232 format (' Leap-frog scheme     : CFLLOW ',e12.4,' CFLHIG',e12.4)
 233 format (' Semi-implicit method : THETAC ',e12.4,' THETAS',e12.4)
 234 format (' Non-hydrostatic pres : THETA  ',e12.4)
 235 format (' Vertical terms u-mom : THETA  ',e12.4)
 236 format (' Vertical terms w-mom : THETA  ',e12.4)
 237 format (' Vertical terms trans : THETA  ',e12.4)
 238 format (' Beneath floating obj : THETA  ',e12.4)
 239 format (' Horz advection u-mom : PROPSC ',i12  ,' KAPPA ',e12.4)
 240 format (' Vert advection u-mom : PROPSC ',i12  ,' KAPPA ',e12.4)
 241 format (' Horz advection w-mom : PROPSC ',i12  ,' KAPPA ',e12.4)
 242 format (' Vert advection w-mom : PROPSC ',i12  ,' KAPPA ',e12.4)
 243 format (' Correction water dep : PROPSC ',i12  ,' KAPPA ',e12.4)
 244 format (' Horz advection trans : PROPSC ',i12  ,' KAPPA ',e12.4)
 245 format (' Vert advection trans : PROPSC ',i12  ,' KAPPA ',e12.4)
 246 format ('                      : M      ',e12.4,' PHI   ',e12.4)
 247 format (' Iterative process    : ACCUR  ',e12.4,' MAXIT ',i12)
 248 format (' Rigid body solver    : ALFAM  ',e12.4,' ALFAF ',e12.4)
 249 format ('                      : BETA   ',e12.4,' GAMMA ',e12.4)
 250 format ('                      : ACCUR  ',e12.4,' MAXIT ',i12)
 251 format (' Sponge layer widths  : LEFT   ',e12.4,' RIGHT ',e12.4)
 252 format ('                      : LOWER  ',e12.4,' UPPER ',e12.4)
 253 format (' Reference date and time         ',a15)
 254 format (' Coriolis forcing     : FACTOR ',e12.4)
    !
end subroutine SwashPrintSettings
