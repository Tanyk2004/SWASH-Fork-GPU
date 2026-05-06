subroutine SwashReqOutQ ( found )
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
!   41.95: Marcel Zijlema
!
!   Updates
!
!    1.00, February 2010: New subroutine
!   41.95,     July 2022: extension to write block output to VTK files
!
!   Purpose
!
!   Reading and processing of the output quantities
!
!   Modules used
!
    use ocpcomm2
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata3
    use SwashCommdata4
    use outp_data
    use m_parall
!
    implicit none
!
!   Argument variables
!
    logical, intent(inout)      :: found          ! keyword found
!
!   Local variables
!
    integer                     :: i              ! loop counter
    integer                     :: idlao          ! lay-out indicator
    integer, save               :: ient = 0       ! number of entries in this subroutine
    integer                     :: ierr           ! error indicator: ierr=0: no error, otherwise error
    integer                     :: ilpos          ! actual length of error message filename
    integer                     :: iproc          ! loop counter over processors
    integer                     :: iostat         ! I/O status in call FOR
    integer                     :: ivar           ! index of output variable
    integer                     :: ivtype         ! type of output quantity
    integer                     :: mip            ! number of output points
    integer                     :: nref           ! reference number
    integer                     :: nvar           ! actual number of output quantities
    !
    real                        :: dfac           ! multiplication factor of block output
    !
    character(80)               :: msgstr         ! string to pass message
    character(len=lenfnm)       :: outdir         ! output directory for time-varying VTK files
    character(4)                :: pnum           ! string containing processor number
    character(16)               :: psname         ! name of point set
    character(1)                :: pstype         ! type of point set ('F', 'H', 'U', etc.)
    character(4)                :: rtype          ! type of output request ('BLKD', 'BLKP', 'TABP', 'TABS', etc.)
    !
    logical                     :: lc             ! indicates low characters in string
    logical                     :: KEYWIS         ! indicates whether keyword in user manual is found or not
    logical                     :: STPNOW         ! indicates whether program must be terminated or not
    logical                     :: vtk            ! indicates whether VTK files are used
    !
    type(orqdat), pointer       :: orqtmp         ! list containing parameters for request output
    type(orqdat), save, pointer :: corq           ! current item in list of request output
    !
    type auxpt                                    ! auxiliary linked list
       integer              :: i
       real                 :: r
       type(auxpt), pointer :: nexti
    end type auxpt
    type(auxpt), target     :: frst
    type(auxpt), pointer    :: curr, tmp
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashReqOutQ')
    !
    ! ==========================================================================
    !
    ! BLOCK   'sname'  HEAD / NOHEAD  'fname' (LAY-OUT [idla])                 &
    !       < TSEC|XP|YP|DEP|BOTL|WATL|DRAF|VMAG|VDIR|VEL|VKSI|VETA|           &
    !         PRESS|NHPRES|QMAG|QDIR|DISCH|QKSI|QETA|VORT|WMAG|WDIR|WIND|      &
    !         FRC|USTAR|UFRIC|ZK|HK|VMAGK|VDIRK|VELK|VKSIK|VETAK|              &
    !         VZ|VOMEGA|QMAGK|QDIRK|DISCHK|QKSIK|QETAK|PRESSK|NHPRSK|          &
    !         TKE|EPS|VISC|HS|HRMS|SETUP|MVMAG|MVDIR|MVEL|MVKSI|MVETA|         &
    !         MVMAGK|MVDIRK|MVELK|MVKSIK|MVETAK|MTKE|MEPS|MVISC|               &
    !         SAL|TEMP|SED|MSAL|MTEMP|MSED|SALK|TEMPK|SEDK|MSALK|MTEMPK|       &
    !         MSEDK [unit] >                                                   &
    !         (OUTPUT [tbegblk] [deltblk] SEC/MIN/HR/DAY)
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('BLO') ) then
       !
       call SWNMPS ( psname, pstype, mip, ierr )
       if ( ierr /= 0 ) goto 800
       if ( pstype /= 'F' .and. pstype /= 'H' .and. pstype /= 'U' ) then
          call msgerr (2, 'set of output locations is not correct type')
          goto 800
       endif
       if ( psname(1:6) == 'NOGRID' ) then
          call msgerr (2, 'empty set is not allowed')
          goto 800
       endif
       !
       ! output frame exists
       !
       allocate(orqtmp)
       nreoq = nreoq + 1
       if ( nreoq > max_outp_req ) call msgerr (2, 'too many output requests')
       !
       idlao = 1
       !
       call INKEYW ('REQ',' ')
       if ( KEYWIS('NOHEAD') .or. KEYWIS ('FIL') ) then
          call INCSTR ('FNAME', FILENM, 'REQ', ' ')
          dfac  = 1.
          rtype = 'BLKD'
       else
          call IGNORE ('HEAD')
          call INCSTR ('FNAME', FILENM, 'STA', ' ')
          dfac  = -1.
          rtype = 'BLKP'
          if ( index ( FILENM, '.MAT' ) /= 0 .or. index ( FILENM, '.mat' ) /= 0 ) then
             call msgerr (4, 'no header allowed for Matlab files')
             return
          endif
       endif
       vtk = index ( FILENM, '.VT' ) /= 0 .or. index (FILENM, '.vt' ) /= 0
       if ( vtk ) then
          dfac  = 1.
          rtype = 'BLKV'
       endif
       if ( FILENM == ' ' ) then
          nref = PRINTF
       else
          nref = 0
          !
          ! append node number to FILENM in case of parallel computing
          !
          if ( parll .and. .not.vtk ) then
             ilpos = index ( FILENM, ' ' )-1
             write (FILENM(ilpos+1:ilpos+4),201) inode
          endif
       endif
       call INKEYW ('STA', ' ')
       if ( KEYWIS('LAY') ) then
          call ININTG ('IDLA', idlao, 'REQ', 0)
          call INKEYW ('REQ', ' ')
          if ( idlao /= 1 .and. idlao /= 3 .and. idlao /= 4 ) call msgerr (2, 'illegal value for [idla]')
       endif
       orqtmp%oqr(1) = -1.
       orqtmp%oqr(2) = -1.
       orqtmp%rqtype = rtype
       orqtmp%psname = psname
       orqtmp%oqi(1) = nref
       orqtmp%oqi(2) = nreoq
       nvar = 0
       orqtmp%oqi(3) = nvar
       orqtmp%oqi(4) = idlao
       outp_files(nreoq) = FILENM
       !
       ! read types of variables to be printed in block
       !
       frst%i = 0
       frst%r = 0.
       nullify(frst%nexti)
       curr => frst
 100   call svartp (ivtype)
       if ( ivtype == 98 ) goto 110
       if ( ivtype /= 999 ) then
          call INREAL ('UNIT', dfac, 'STA', -1.)
          if ( (ivtype == 95 .or. ivtype == 96) .and. optg /= 5 ) then
             call msgerr (1, 'horizontal divergence operators will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( ivtype > 50 .and. ivtype < 100 .and. ivtype /= 95 .and. ivtype /= 96 .and. kmax == 1 ) then
             call msgerr (1, 'this layer-dependent quantity will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( (ivtype == 10 .or. ivtype == 11 .or. ivtype == 75 .or. ivtype == 76) .and. optg == 5 ) then
             call msgerr (1, 'this grid-oriented velocity component will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( (ivtype == 36 .or. ivtype == 37 .or. ivtype == 87 .or. ivtype == 88) .and. optg == 5 ) then
             call msgerr (1, 'this mean grid-oriented velocity component will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( (ivtype == 17 .or. ivtype == 18 .or. ivtype == 80 .or. ivtype == 81) .and. optg == 5 ) then
             call msgerr (1, 'this grid-oriented discharge component will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( ivtype > 0 ) then
             nvar = nvar + 1
             allocate(tmp)
             tmp%i = ivtype
             tmp%r = dfac
             nullify(tmp%nexti)
             curr%nexti => tmp
             curr => tmp
             if ( ivtype == 22 .or. ivtype == 23 .or. ivtype == 24 ) lwavoutp = .true.
             if ( ivtype >= 33 .and. ivtype <= 37 ) lcuroutp = .true.
             if ( ivtype == 42 .or. ivtype == 43 .or. ivtype == 44 ) ltraoutp = .true.
             if ( ivtype == 57 .or. ivtype == 58 .or. ivtype == 59 ) lturoutp = .true.
             if ( ivtype >= 84 .and. ivtype <= 88 ) lcuroutp = .true.
             if ( ivtype == 92 .or. ivtype == 93 .or. ivtype == 94 ) ltraoutp = .true.
          endif
          goto 100
       endif
       !
 110   if ( nvar > 0 ) then
          allocate(orqtmp%ivtyp(nvar))
          allocate(orqtmp%fac(nvar))
          curr => frst%nexti
          do i = 1, nvar
             orqtmp%ivtyp(i) = curr%i
             orqtmp%fac  (i) = curr%r
             curr => curr%nexti
          enddo
          deallocate(tmp)
       endif
       !
       if ( ivtype == 98 ) then
          if ( nstatm == 0 ) call msgerr (3, 'time information not allowed in stationary mode')
          nstatm = 1
          call INCTIM (ITMOPT, 'TBEG', orqtmp%oqr(1), 'REQ', 0d0)
          call INITVD ('DELT', orqtmp%oqr(2), 'REQ', 0d0)
          !
          do i = 1, nvar
             ivar = orqtmp%ivtyp(i)
             if ( ivar == 22 .or.  ivar == 23  .or. ivar == 24 .or.  &
                  ivar == 42 .or.  ivar == 43  .or. ivar == 44 .or.  &
                  ivar == 57 .or.  ivar == 58  .or. ivar == 59 .or.  &
                  ivar == 92 .or.  ivar == 93  .or. ivar == 94 .or.  &
                 (ivar >= 33 .and. ivar <= 37) .or.                  &
                 (ivar >= 84 .and. ivar <= 88) ) then
                call msgerr (1, 'command OUTput ignored for this stationary quantity')
                write (PRINTF, *) ' -> ', ovkeyw(ivar)
                orqtmp%oqr(1) = -1.
                orqtmp%oqr(2) = -1.
                goto 115
             endif
          enddo
          !
          ! a PVD file is created to collect time-varying output
          !
          if ( vtk .and. INODE == MASTER ) then
             lc = .false.
             ilpos = index( FILENM, '.VT' )
             if ( ilpos == 0 ) then
                lc = .true.
                ilpos = index( FILENM, '.vt' )
             endif
             if ( lc ) then
                write (FILENM(ilpos+1:ilpos+3),202) 'pvd'
             else
                write (FILENM(ilpos+3:ilpos+3),202) 'PVD'
             endif
             nref   =  0
             iostat = -1
             call FOR (nref, FILENM, 'UF', iostat)
             if (STPNOW()) return
             upvdf(nreoq) = nref
             !
             ! write the header lines
             !
             write (nref,'(a)') trim(xmllin1)
             write (nref,'(a)') trim(xmllin2)
             vtkline = '  This file was generated by SWASH version '//  &
                       trim(VERTXT)//'; project: '//trim(PROJID)//      &
                       '; run number: '//trim(PROJNR)
             write (nref,'(a)') trim(vtkline)
             write (nref,'(a)') trim(xmllin3)
             write (nref,'(a)') trim(pvdlin1)
             write (nref,'(a)') trim(pvdlin2)
             !
             ! make output directory
             !
             if ( lc ) then
                outdir = FILENM(1:ilpos-1)//'_output'
             else
                outdir = FILENM(1:ilpos-1)//'_OUTPUT'
             endif
             call MKPATH ( outdir, ierr )
             if ( ierr /= 0 ) outdir = '.'
             vtkdir(nreoq) = outdir
             !
             ! create subdirectories to store each piece of output data
             !
             if ( PARLL ) then
                do iproc = 1, NPROC
                   if ( lc ) then
                      write (pnum(1:4),203) 'p',iproc
                   else
                      write (pnum(1:4),203) 'P',iproc
                   endif
                   call MKPATH ( trim(outdir)//DIRCH2//pnum, ierr )
                enddo
                if ( ierr /= 0 ) then
                   write (msgstr, '(a,i5)') 'Error while creating folders '//  &
                                            '- status error =',ierr
                   call msgerr( 3, trim(msgstr) )
                endif
             endif
          endif
       endif
       !
 115   continue
       orqtmp%oqi(3) = nvar
       if ( nvar == 0 ) then
          allocate(orqtmp%ivtyp(0))
          allocate(orqtmp%fac(0))
       endif
       nullify(orqtmp%nextorq)
       if ( .not.lorq ) then
           forq = orqtmp
           corq => forq
           lorq = .true.
       else
          corq%nextorq => orqtmp
          corq => orqtmp
       endif
       goto 800
       !
    endif
    !
    ! ==========================================================================
    !
    ! TABLE   'sname'  HEAD | NOHEAD | STAB | SWASH | IND  'fname'             &
    !       < TIME|TSEC|XP|YP|DIST|DEP|BOTL|WATL|DRAF|VMAG|VDIR|VEL|VKSI|VETA| &
    !         PRESS|NHPRES|QMAG|QDIR|DISCH|QKSI|QETA|VORT|WMAG|WDIR|WIND|FRC|  &
    !         USTAR|UFRIC|ZK|HK|VMAGK|VDIRK|VELK|VKSIK|VETAK|VZ|VOMEGA|        &
    !         QMAGK|QDIRK|DISCHK|QKSIK|QETAK|PRESSK|NHPRSK|TKE|EPS|VISC|       &
    !         HS|HRMS|SETUP|MVMAG|MVDIR|MVEL|MVKSI|MVETA|MVMAGK|MVDIRK|        &
    !         MVELK|MVKSIK|MVETAK|MTKE|MEPS|MVISC|SAL|TEMP|SED|MSAL|MTEMP|     &
    !         MSED|SALK|TEMPK|SEDK|MSALK|MTEMPK|MSEDK|FORCEX|FORCEY|FORCEZ|    &
    !         MOMX|MOMY|MOMZ|TRAX|TRAY|TRAZ|ROTX|ROTY|ROTZ|PTOP|RUNUP >        &
    !         (OUTPUT [tbegtbl] [delttbl] SEC/MIN/HR/DAY)
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('TAB') ) then
       !
       call SWNMPS ( psname, pstype, mip, ierr )
       if ( ierr /= 0 ) goto 800
       !
       ! output points exist
       !
       allocate(orqtmp)
       nreoq = nreoq + 1
       if ( nreoq > max_outp_req ) call msgerr (2, 'too many output requests')
       !
       call INKEYW ('STA','HEAD')
       if ( KEYWIS('NOHEAD') .or. KEYWIS ('FIL') ) then
          rtype = 'TABD'
       else if ( KEYWIS ('IND') ) then
          rtype = 'TABI'
       else if ( KEYWIS ('SWASH') ) then
          rtype = 'TABS'
       else if ( KEYWIS ('STAB') ) then
          rtype = 'TABT'
       else
          call IGNORE ('HEAD')
          rtype = 'TABP'
       endif
       orqtmp%oqr(1) = -1.
       orqtmp%oqr(2) = -1.
       orqtmp%rqtype = rtype
       call INCSTR ('FNAME', FILENM, 'STA', ' ')
       if ( FILENM /= '    ' ) then
          nref = 0
          !
          ! append node number to FILENM in case of parallel computing
          !
          if ( parll ) then
             ilpos = index ( FILENM, ' ' )-1
             write (FILENM(ilpos+1:ilpos+4),201) inode
          endif
       else
          nref = PRINTF
       endif
       orqtmp%psname = psname
       orqtmp%oqi(1) = nref
       orqtmp%oqi(2) = nreoq
       outp_files(nreoq) = FILENM
       !
       nvar = 0
       orqtmp%oqi(3) = nvar
       !
       ! read types of variables to be printed in table
       !
       frst%i = 0
       nullify(frst%nexti)
       curr => frst
 120   call svartp (ivtype)
       if ( ivtype == 98 ) goto 130
       if ( ivtype /= 999 ) then
          if ( ivtype < 100 .and. ivtype /= 41 .and. psname(1:6) == 'NOGRID' ) then
             call msgerr (2, 'empty set not allowed for this quantity')
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          endif
          if ( ivtype == 115 .and. .not.oned ) then
             call msgerr (1, 'wave runup height will not be computed in 2D mode')
          endif
          if ( (ivtype == 95 .or. ivtype == 96) .and. optg /= 5 ) then
             call msgerr (1, 'horizontal divergence operators will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( ivtype > 50 .and. ivtype < 100 .and. ivtype /= 95 .and. ivtype /= 96 .and. kmax == 1 ) then
             call msgerr (1, 'this layer-dependent quantity will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( (ivtype == 10 .or. ivtype == 11 .or. ivtype == 75 .or. ivtype == 76) .and. optg == 5 ) then
             call msgerr (1, 'this grid-oriented velocity component will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( (ivtype == 36 .or. ivtype == 37 .or. ivtype == 87 .or. ivtype == 88) .and. optg == 5 ) then
             call msgerr (1, 'this mean grid-oriented velocity component will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( (ivtype == 17 .or. ivtype == 18 .or. ivtype == 80 .or. ivtype == 81) .and. optg == 5 ) then
             call msgerr (1, 'this grid-oriented discharge component will not be outputted' )
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( ivtype > 0 ) then
             nvar = nvar + 1
             allocate(tmp)
             tmp%i = ivtype
             nullify(tmp%nexti)
             curr%nexti => tmp
             curr => tmp
             if ( ivtype == 22 .or. ivtype == 23 .or. ivtype == 24 ) lwavoutp = .true.
             if ( ivtype >= 33 .and. ivtype <= 37 ) lcuroutp = .true.
             if ( ivtype == 42 .or. ivtype == 43 .or. ivtype == 44 ) ltraoutp = .true.
             if ( ivtype == 57 .or. ivtype == 58 .or. ivtype == 59 ) lturoutp = .true.
             if ( ivtype >= 84 .and. ivtype <= 88 ) lcuroutp = .true.
             if ( ivtype == 92 .or. ivtype == 93 .or. ivtype == 94 ) ltraoutp = .true.
             call INKEYW ('STA', ' ')
             if ( KEYWIS('UNIT') ) then
                call msgerr (1, '[unit] is ignored in this version')
             endif
          endif
          goto 120
       endif
       !
 130   if ( nvar > 0 ) then
          allocate(orqtmp%ivtyp(nvar))
          curr => frst%nexti
          do i = 1, nvar
             orqtmp%ivtyp(i) = curr%i
             curr => curr%nexti
          enddo
          deallocate(tmp)
       endif
       allocate(orqtmp%fac(0))
       !
       if ( ivtype == 98 ) then
          if ( nstatm == 0 ) call msgerr (3, 'time information not allowed in stationary mode')
          nstatm = 1
          call INCTIM (ITMOPT, 'TBEG', orqtmp%oqr(1), 'REQ', 0d0)
          call INITVD ('DELT', orqtmp%oqr(2), 'REQ', 0d0)
          do i = 1, nvar
             ivar = orqtmp%ivtyp(i)
             if ( ivar == 22 .or.  ivar == 23  .or. ivar == 24 .or.  &
                  ivar == 42 .or.  ivar == 43  .or. ivar == 44 .or.  &
                  ivar == 57 .or.  ivar == 58  .or. ivar == 59 .or.  &
                  ivar == 92 .or.  ivar == 93  .or. ivar == 94 .or.  &
                 (ivar >= 33 .and. ivar <= 37) .or.                  &
                 (ivar >= 84 .and. ivar <= 88) ) then
                call msgerr (1, 'command OUTput ignored for this stationary quantity')
                write (PRINTF, *) ' -> ', ovkeyw(ivar)
                orqtmp%oqr(1) = -1.
                orqtmp%oqr(2) = -1.
             endif
          enddo
       endif
       orqtmp%oqi(3) = nvar
       if ( nvar == 0 ) then
          allocate(orqtmp%ivtyp(0))
          allocate(orqtmp%fac(0))
       endif
       nullify(orqtmp%nextorq)
       if ( .not.lorq ) then
          forq = orqtmp
          corq => forq
          lorq = .true.
       else
          corq%nextorq => orqtmp
          corq => orqtmp
       endif
       goto 800
    endif
    !
    ! command not found
    !
    return
    !
 800 found = .true.
    return
    !
 201 format('-',i3.3)
 202 format(a3)
 203 format(a1,i3.3)
    !
end subroutine SwashReqOutQ
