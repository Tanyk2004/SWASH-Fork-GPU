subroutine SwashInputField ( logcom )
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
!    1.00, February 2010: New subroutine
!    1.10,    March 2012: extension porous media
!    1.21,     July 2013: extension transport constituents
!    1.22,   August 2013: extension vegetation
!    1.50, December 2013: extension to read nonuniform input
!    4.05,  January 2017: extension floating objects
!
!   Purpose
!
!   Reads and stores actual input field
!
!   Modules used
!
    use ocpcomm2, only: LENFNM
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use m_genarr
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    logical, dimension(6), intent(inout) :: logcom   ! indicates which commands have been given to know if all the
                                                     ! information for a certain command is available. Meaning:
                                                     ! (1) no meaning
                                                     ! (2) command CGRID has been carried out
                                                     ! (3) command READINP BOTTOM has been carried out
                                                     ! (4) command READ COOR has been carried out
                                                     ! (5) command READ UNSTRUC has been carried out
                                                     ! (6) arrays s1, u1 and v1 have been allocated
!
!   Local variables
!
    integer                              :: i        ! loop counter
    integer, save                        :: ient = 0 ! number of entries in this subroutine
    integer                              :: igr1     ! grid number for which a scalar or first component of a vector is read
    integer                              :: igr2     ! grid number for which second component of a vector is read
    integer                              :: iostat   ! I/O status in call FOR
    integer                              :: j        ! loop counter
    integer                              :: k        ! loop counter over layers
    integer                              :: nhedc    ! number of heading lines in the file before each array or vector component
    !
    real, dimension(:), allocatable      :: tarr     ! temporary array
    !
    logical                              :: KEYWIS   ! indicates whether keyword in user manual is found or not
    logical                              :: STPNOW   ! indicates that program must stop
    logical                              :: vector   ! indicates whether parameter is a vector or not
    !
    character(LEN=LENFNM)                :: filenm   ! file name
    character(80)                        :: hedlin   ! heading line
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    call strace (ient,'SwashInputField')
    !
    ! ==========================================================================
    !
    ! READinp  BOTtom/WLevel/CURrent/FRiction/WInd/PRes/COOR/CORI/PORO/PSIZ/   &
    !          HSTRUC/NPLAnts/DRAFt/LABel/SALinity/TEMPerature/SEDiment/MWL/   &
    !          ACURrent                                                        &
    !
    !              | 'fname1'        |
    !      [fac]  <  SERIes 'fname2'  >  [idla] [nhedf] ([nhedt]) (nhedvec])   &
    !              | LAYers 'fname3' |
    !
    !      FREE / FORMAT 'form' / [idfm] / UNFORMATTED
    !
    ! ==========================================================================
    !
    call INKEYW ('REQ',' ')
    igr2 = 0
    if ( KEYWIS ('BOT') ) then
       igr1      = 1
       logcom(3) = .true.
    else if ( KEYWIS ('CUR') .or. KEYWIS ('VEL') ) then
       igr1 = 2
       igr2 = 3
    else if ( KEYWIS ('FR') ) then
       igr1  = 4
       varfr = .true.
    else if ( KEYWIS ('WI') ) then
       igr1  = 5
       igr2  = 6
       varwi = .true.
    else if ( KEYWIS ('WL') ) then
       igr1 = 7
    else if ( KEYWIS ('COOR') ) then
       igr1      = 8
       igr2      = 9
       logcom(4) = .true.
       mxg(igr1) = mxc - 1
       myg(igr1) = myc - 1
       mxg(igr2) = mxc - 1
       myg(igr2) = myc - 1
    else if ( KEYWIS ('PORO') ) then
       igr1   = 10
       iporos = 1
    else if ( KEYWIS ('PSIZ') ) then
       igr1  = 11
       vargs = .true.
    else if ( KEYWIS ('HSTRUC') ) then
       igr1  = 12
       varsh = .true.
    else if ( KEYWIS ('PR') ) then
       igr1  = 13
       svwp  = .true.
    else if ( KEYWIS ('NPLA') ) then
       igr1   = 14
       varnpl = .true.
    else if ( KEYWIS ('SAL') ) then
       igr1   = 15
       itrans = 1
       ltrans = ltrans + 1
       lsal   = ltrans
    else if ( KEYWIS ('TEMP') ) then
       igr1   = 16
       itrans = 1
       ltrans = ltrans + 1
       ltemp  = ltrans
    else if ( KEYWIS ('SED') ) then
       igr1   = 17
       itrans = 1
       ltrans = ltrans + 1
       lsed   = ltrans
    else if ( KEYWIS ('SHIP') .or. KEYWIS ('VESS') .or. KEYWIS ('WEC') .or. KEYWIS ('FLOAT') .or. KEYWIS ('DRAF') .or. KEYWIS ('PONT') .or. KEYWIS ('PFLOW') ) then
       igr1   = 18
       ifloat = 1
    else if ( KEYWIS ('LAB') .or. KEYWIS ('BODY') ) then
       igr1   = 19
       ifloat = 2
    else if ( KEYWIS ('ACUR') .or. KEYWIS ('BCUR') .or. KEYWIS ('CCUR') ) then
       igr1 = 20
       igr2 = 21
       varva = .true.
    else if ( KEYWIS ('MWL') ) then
       igr1 = 22
       varea = .true.
    else if ( KEYWIS ('CORI') ) then
       igr1     = 23
       coriolis = .true.
    else
       call WRNKEY
    endif
    !
    ! read multiplication factor
    !
    call INREAL ('FAC', iflfac(igr1), 'STA', 1.)
    if ( igr1 == 19 .and. iflfac(19) /= 1. ) then
       call msgerr (1, 'for labelling rigid bodies, [fac] must be set to 1.0')
       iflfac(19) = 1.
    endif
    !
    if ( igr2 > 0 ) then
       vector       = .true.
       iflfac(igr2) = iflfac(igr1)
    else
       vector       = .false.
    endif
    !
    call REPARM (iflndf(igr1), iflnds(igr1), iflidl(igr1) , &
                 iflifm(igr1), iflfrm(igr1), iflnhf(igr1) , &
                 ifldyn(igr1), iflnhd(igr1), vector, nhedc)
    if (STPNOW()) return
    !
    if ( ITEST >= 60 ) write (PRTEST, 201) igr1, iflndf(igr1),   &
       iflnds(igr1), iflidl(igr1), iflifm(igr1), iflfrm(igr1),   &
       iflnhf(igr1), ifldyn(igr1), iflnhd(igr1), vector, nhedc
    !
    iflnhd(igr1) = iflnhd(igr1) + nhedc
    if ( igr2 > 0 ) then
       iflndf(igr2) = iflndf(igr1)
       iflnds(igr2) = iflnds(igr1)
       iflidl(igr2) = iflidl(igr1)
       iflifm(igr2) = iflifm(igr1)
       iflfrm(igr2) = iflfrm(igr1)
       iflnhf(igr2) = iflnhf(igr1)
       iflnhd(igr2) = nhedc
    endif
    !
    if ( leds(igr1) == 0 .and. igr1 /= 8 ) then
       call msgerr (4, 'input grid is not given')
       return
    endif
    !
    if ( igr1 == 1 .and. grid_generator == meth_adcirc ) then
       call msgerr (1, 'depth will be taken from grid file fort.14 ')
       igtype(1) = 3
       leds(1)   = 2
       return
    endif
    !
    do k = 1, ifllay(igr1)
       !
       if (.not.allocated(tarr) ) then
          allocate( tarr(mxg(igr1)*myg(igr1)) )
          tarr = 0.
       endif
       !
       call INAR2D( tarr      , mxg(igr1), myg(igr1), iflndf(igr1), &
                    iflnds(igr1), iflifm(igr1), iflfrm(igr1),       &
                    iflidl(igr1), iflfac(igr1),                     &
                    iflnhd(igr1), iflnhf(igr1))
       if (STPNOW()) return
       !
       if ( igr1 == 1 ) then
          !
          if (.not.allocated(depth)) allocate(depth(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, depth, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 2 ) then
          !
          if (.not.allocated(uxb)) allocate(uxb(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, uxb, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 4 ) then
          !
          if (.not.allocated(fric)) allocate(fric(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, fric, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 5 ) then
          !
          if (.not.allocated(wxi)) allocate(wxi(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, wxi, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 7 ) then
          !
          if (.not.allocated(wlevl)) allocate(wlevl(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, wlevl, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 8 ) then
          !
          do j = 1, myg(igr1)
             do i = 1, mxg(igr1)
                xcgrid(i,j) = tarr((j-1)*mxg(igr1)+i)
             enddo
             xcgrid(mxc,j) = xcgrid(mxc-1,j)
          enddo
          xcgrid(:,myc) = xcgrid(:,myc-1)
          !
       else if ( igr1 == 10 ) then
          !
          if (.not.allocated(npor)) allocate(npor(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, npor, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 11 ) then
          !
          if (.not.allocated(gsiz)) allocate(gsiz(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, gsiz, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 12 ) then
          !
          if (.not.allocated(hpor)) allocate(hpor(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, hpor, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 13 ) then
          !
          if (.not.allocated(pres)) allocate(pres(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, pres, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 14 ) then
          !
          if (.not.allocated(npla)) allocate(npla(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, npla, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 15 ) then
          !
          if (.not.allocated(sal)) allocate(sal(mxg(igr1)*myg(igr1)*ifllay(igr1)))
          call SWCOPR( tarr, sal(1+(k-1)*mxg(igr1)*myg(igr1):k*mxg(igr1)*myg(igr1)), mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 16 ) then
          !
          if (.not.allocated(temp)) allocate(temp(mxg(igr1)*myg(igr1)*ifllay(igr1)))
          call SWCOPR( tarr, temp(1+(k-1)*mxg(igr1)*myg(igr1):k*mxg(igr1)*myg(igr1)), mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 17 ) then
          !
          if (.not.allocated(sed)) allocate(sed(mxg(igr1)*myg(igr1)*ifllay(igr1)))
          call SWCOPR( tarr, sed(1+(k-1)*mxg(igr1)*myg(igr1):k*mxg(igr1)*myg(igr1)), mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 18 ) then
          !
          if (.not.allocated(flobj)) allocate(flobj(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, flobj, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 19 ) then
          !
          if (.not.allocated(lbod)) allocate(lbod(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, lbod, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 20 ) then
          !
          if (.not.allocated(avxi)) allocate(avxi(mxg(igr1)*myg(igr1)*ifllay(igr1)))
          call SWCOPR( tarr, avxi(1+(k-1)*mxg(igr1)*myg(igr1):k*mxg(igr1)*myg(igr1)), mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 22 ) then
          !
          if (.not.allocated(mwl)) allocate(mwl(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, mwl, mxg(igr1)*myg(igr1) )
          !
       else if ( igr1 == 23 ) then
          !
          if (.not.allocated(corp)) allocate(corp(mxg(igr1)*myg(igr1)))
          call SWCOPR( tarr, corp, mxg(igr1)*myg(igr1) )
          !
       endif
       !
       deallocate(tarr)
       !
       if ( igr2 > 0 ) then
          !
          if (.not.allocated(tarr) ) then
             allocate( tarr(mxg(igr2)*myg(igr2)) )
             tarr = 0.
          endif
          !
          call INAR2D( tarr      , mxg(igr2), myg(igr2), iflndf(igr2), &
                       iflnds(igr2), iflifm(igr2), iflfrm(igr2),       &
                       iflidl(igr2), iflfac(igr2),                     &
                       iflnhd(igr2), 0)
          if (STPNOW()) return
          !
          if ( igr2 == 3 ) then
             !
             if (.not.allocated(uyb)) allocate(uyb(mxg(igr2)*myg(igr2)))
             call SWCOPR( tarr, uyb, mxg(igr2)*myg(igr2) )
             !
          else if ( igr2 == 6 ) then
             !
             if (.not.allocated(wyi)) allocate(wyi(mxg(igr2)*myg(igr2)))
             call SWCOPR( tarr, wyi, mxg(igr2)*myg(igr2) )
             !
          else if ( igr2 == 9 ) then
             !
             do j = 1, myg(igr2)
                do i = 1, mxg(igr2)
                   ycgrid(i,j) = tarr((j-1)*mxg(igr2)+i)
                enddo
                ycgrid(mxc,j) = ycgrid(mxc-1,j)
             enddo
             ycgrid(:,myc) = ycgrid(:,myc-1)
             !
          else if ( igr2 == 21 ) then
             !
             if (.not.allocated(avyi)) allocate(avyi(mxg(igr2)*myg(igr2)*ifllay(igr2)))
             call SWCOPR( tarr, avyi(1+(k-1)*mxg(igr2)*myg(igr2):k*mxg(igr2)*myg(igr2)), mxg(igr2)*myg(igr2) )
             !
          endif
          !
          deallocate(tarr)
          !
       endif
       !
       if ( k < ifllay(igr1) ) then
          close(iflnds(igr1))
          read (iflndf(igr1), '(a)', end=901) filenm
          if ( iflifm(igr1) /= -1 ) then
             iostat = 0
             call FOR ( iflnds(igr1), filenm, 'OF', iostat )
             if (STPNOW()) return
          else
             iostat = 0
             call FOR ( iflnds(igr1), filenm, 'OU', iostat )
             if (STPNOW()) return
          endif
          !
          ! read heading lines
          !
          if ( iflnhf(igr1) > 0 ) then
             if ( iflifm(igr1) < 0 ) then
                if ( ITEST >= 30 ) write (PRINTF, '(i3,a,a)') iflnhf(igr1),' heading lines at begin of file ', trim(filenm)
                do j = 1, iflnhf(igr1)
                   read (iflnds(igr1))
                enddo
             else
                write (PRINTF, '(a,a,a)') ' **  Heading lines file ', trim(filenm), ' **'
                do j = 1, iflnhf(igr1)
                   read (iflnds(igr1), '(a80)') hedlin
                   write (PRINTF, '(a4,a80)') ' -> ', hedlin
                enddo
             endif
          endif
       endif
       !
    enddo
    !
    ! set time of reading
    !
    if ( ifldyn(igr1) ==  1 ) then
       if ( nstatm == 0 ) call msgerr (2, 'nonstationary input field requires MODE NONSTAT')
       ifltim(igr1) = iflbeg(igr1)
       if ( igr2 > 0 ) ifltim(igr2) = ifltim(igr1)
    endif
    !
    leds(igr1) = 2
    if ( vector ) leds(igr2) = 2
    !
    return
    !
 201 format (' Reading parameters: ', 5i4, a, /, 12x,i5,i5,i5,l2,i5)
    !
 901 filenm = 'UNKNOWN_FILE'
     inquire ( unit=iflndf(igr1), name=filenm )
     call msgerr (2, 'layers of input files ended in '//trim(filenm))
    !
end subroutine SwashInputField
