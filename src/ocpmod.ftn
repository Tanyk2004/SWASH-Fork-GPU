!                 COMMON VARIABLES RELATED MODULES FOR OCEAN PACK
!
!     Contents of this file
!
!     OCPCOMM1           contains common variables for Ocean Pack
!     OCPCOMM2           contains common variables for Ocean Pack
!     OCPCOMM3           contains common variables for Ocean Pack
!     OCPCOMM4           contains common variables for Ocean Pack
!
      MODULE OCPCOMM1
!
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
!  0. Authors
!
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file OCPCOMM1.INC
!
!  2. Purpose
!
!     Common variables used by the Ocean Pack Service Routines and in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     LINELN [ 180]  max length of input lines (in command file)
!
      INTEGER LINELN
      PARAMETER  (LINELN=180)
!
!  7. Local variables
!
!     *** character data used by the command reading system ***
!
! BLANK  ['    '] blank string
! COMID  [   '$'] character which distinguishes comments in the command input
! ELTEXT [      ] contents of the last string read by reading system
! ELTYPE [      ] type of the element last read by reading system
!                 ='CHAR'; last read data element is the string in ELTEXT
!                 ='EMPT;  empty data field
!                 ='EOF';  end of file has been reached
!                 ='EOR';  end of repeat has been reached
!                 ='ERR';  incorrect data field was encountered
!                 ='INT';  last read data element is the integer in ELINT
!                 ='KEY';  last read data element is the keyword in KEYWRD
!                 ='OTHR'; other
!                 ='REAL'; last read data element is the real in ELREAL
!                 ='USED'; last read data element is processed, new can be read
! KAART  [      ] contents of the input line last read by reading system
! KAR    [      ] character last read by reading system
!                 =COMID; begin or end of comment
!                 =TABC; data element separation mark
!                 =' '; data element separation mark
!                 ='&'; continuation mark
!                 ='('; begin of data group
!                 ='*'; repetition mark
!                 =','; data element seperation mark (only for numbers)
!                 ='/'; end of repetition mark
!                 =':'; assignment mark
!                 =';'; end of record or input line
!                 ='='; assignment mark
!                 ='@'; end of file mark
!                 ='_'; continuation mark
!                 other: letter or digit to be processed by reading system
! KEYWRD [      ] contents of the last keyword read by reading system
! TABC   [CALCUL] =CHAR(9); tabular character
!
      CHARACTER*4        BLANK
      CHARACTER          COMID
      CHARACTER*(LINELN) ELTEXT
      CHARACTER*4        ELTYPE
      CHARACTER*(LINELN) KAART
      CHARACTER          KAR
      CHARACTER*8        KEYWRD
      CHARACTER          TABC
!
!     *** numerical data used by the command reading system ***
!
! CHGVAL [      ] whether last read value is different from a given value for
!                 subroutines INREAL, ININTG, INCSTR, INCTIM
! ELINT  [      ] last element read from user command, when integer
! KARNR  [      ] position on the input line of character last processed
!                 by the reading system,
!                 =0; no characters read yet
!                 =81; next input line has to be read to the common KAART first
! LENCST [      ] length of the string stored in ELTEXT
! ELREAL [      ] last element read from user command, when real or double
!
      INTEGER          ELINT, KARNR, LENCST
      DOUBLE PRECISION ELREAL
      LOGICAL          CHGVAL
!
!     *** origin for day and time ***
!
! REFDAY    [    ] Day number of the reference day. The first day entered is used
!                  as reference day, the reference time is 0:00 of the reference day.
!
      INTEGER REFDAY
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ---
!
! 13. Source text
!
      END MODULE OCPCOMM1

      MODULE OCPCOMM2
!
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
!  0. Authors
!
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file OCPCOMM2.INC
!
!  2. Purpose
!
!     Common variables used by the Ocean Pack Service Routines and in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     LENFNM [ 140]  length of file names (including path)
!
      INTEGER LENFNM
      PARAMETER (LENFNM=140)
!
!  7. Local variables
!
!     *** names and other character strings ***
!
! DIRCH1 [\     ] directory separation character as appears in input file
! DIRCH2 [\     ] directory separation character replacing DIRCH1
! FILEA  [      ] not used
! FILEB  [      ] not used
! FILENM [      ] file name of the file currently used for I/O
! INST   ['Delft University of Technology'] name of the institute
!                 can be changed in the file swashinit
! PROJID ['SWASH'] acronym of the project for which the computation is taking place
!                 ='NAME'; set by command PROJ 'NAME' ...
! PROJNR [CALCUL] =BLANK; run number for the computation
!                 ='NR'; set by command PROJ ... 'NR' ...
! PROJT1 [CALCUL] =BLANK; 1st line of the project title
!                 ='title1'; set by command PROJ ... 'title1' ...
! PROJT2 [CALCUL] =BLANK; 2nd line of the project title
!                 ='title2'; set by command PROJ ... 'title2' ...
! PROJT3 [CALCUL] =BLANK; 3rd line of the project title
!                 ='title3'; set by command PROJ ... 'title3'
! PTITLE [      ] not used
! VERTXT [calcul] program version, character representation
!
      CHARACTER (LEN=1)      :: DIRCH1, DIRCH2
      CHARACTER (LEN=LENFNM) :: FILEA , FILEB , FILENM
      CHARACTER (LEN=40)     :: INST
      CHARACTER (LEN=16)     :: PROJID
      CHARACTER (LEN=4)      :: PROJNR
      CHARACTER (LEN=72)     :: PROJT1, PROJT2, PROJT3
      CHARACTER (LEN=36)     :: PTITLE
      CHARACTER (LEN=20)     :: VERTXT
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ---
!
! 13. Source text
!
      END MODULE OCPCOMM2

      MODULE OCPCOMM3
!
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
!  0. Authors
!
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file OCPCOMM3.INC
!
!  2. Purpose
!
!     Common variables used by the Ocean Pack Service Routines and in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     ---
!
!  7. Local variables
!
!     *** data for output, mainly plotting ***
!
! DXQ    [      ] mesh size of the output frame in X-direction
!                 =0.01; if MXQ=1
!                 =XQLEN/(MXQ-1); if MXQ>1
! DYQ    [      ] mesh size of the output frame in Y-direction
!                 =0.01; if MYQ=1
!                 =YQLEN/(MYQ-1); if MYQ>1
! MXQ    [CALCUL] number of grid points of the output frame in X-direction
! MYQ    [CALCUL] number of grid points of the output frame in Y-direction
! VERNUM [  1.00] version number of SWASH
!
      INTEGER MXQ, MYQ
      REAL    DXQ, DYQ, VERNUM
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ---
!
! 13. Source text
!
      END MODULE OCPCOMM3

      MODULE OCPCOMM4
!
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
!  0. Authors
!
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Oct. 04: taken from the include file OCPCOMM4.INC
!
!  2. Purpose
!
!     Common variables used by the Ocean Pack Service Routines and in SWAN
!
!  3. Method
!
!     MODULE construct
!
!  4. Modules used
!
!     ---
!
      IMPLICIT NONE
!
!  5. Argument variables
!
!     ---
!
!  6. Parameter variables
!
!     INAN   [   1]  integer representing not a number
!     RNAN   [   1]  real representing not a number
!
      INTEGER INAN
      REAL    RNAN
      PARAMETER (INAN=-1073750760,
     &           RNAN=-1.07374515E+09)
!
!  7. Local variables
!
!     *** file unit reference numbers ***
!
! EXPORT [     ] not used
! FUNHI  [  CAL] highest free unit number
!                =IUNMAX
! FUNLO  [   21] lowest free unit number
! HIOPEN [     ] highest unit number of an open file
! IMPORT [     ] not used
! INPUTF [    3] unit number for the file with command input ('INPUT')
!                set by file swashinit
! ITMOPT [    1] time coding option
!                =1; yyyymmdd.hhmmss (ISO notation)
!                =7; hhmmss.msc (SWASH notation)
!                set by file swashinit
! IUNMIN [    0] minimum unit number
! IUNMAX [99999] maximum unit number
!                set by file swashinit
! PRINTF [    4] unit number for the file with standard output ('PRINT')
!                set by file swashinit
! PRTEST [  CAL] unit number for the print file containing test output
!                =PRINTF
!                set by file swashinit
! SCREEN [    6] unit number for the screen
!                (is for batch-oriented systems equal to PRINTF)
!
      INTEGER EXPORT, FUNHI , FUNLO , HIOPEN
      INTEGER IMPORT, INPUTF, ITMOPT, IUNMAX
      INTEGER IUNMIN, PRINTF, PRTEST, SCREEN

!     *** test parameters ***
!
! ITEST  [  0] indicates the amount of test output requested
!              =30; for command TEST
!              =itest; set by command TEST [itest] [itrace]
! ITRACE [  0] a message is printed up to ITRACE times
!              =itrace; set by command TEST [itest] [itrace]
! LEVERR [  0] severity of the errors encountered.
! LTRACE [.F.] indicates whether to call STRACE
!              =.T.; when ITRACE>0
! MAXERR [  1] maximum severity of errors allowed, if larger no computation
!              =1; warnings
!              =2; errors
!              =3; severe errors
!              =4; terminating errors
!              =maxerr; set by command SET ... [maxerr] ...
!
      INTEGER ITEST, ITRACE, LEVERR, MAXERR
      LOGICAL LTRACE
!
!  8. Subroutines and functions used
!
!     ---
!
!  9. Subroutines and functions calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ---
!
! 13. Source text
!
      END MODULE OCPCOMM4
