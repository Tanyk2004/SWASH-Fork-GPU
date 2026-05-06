      MODULE M_PARALL
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
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
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     Dec. 03: New Module
!     Jul. 04: introduction logicals
!
!  2. Purpose
!
!     Contains data with respect to parallel process
!     based on distributed-memory apprach using MPI
!
!  3. Method
!
!     ---
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
!     IHALOX  : width of halo area in x-direction
!     IHALOY  : width of halo area in y-direction
!     MASTER  : rank of master process
!
      INTEGER MASTER
      INTEGER IHALOX, IHALOY
      PARAMETER (MASTER = 1,
     &           IHALOX = 3, IHALOY = 3)
!
!  7. Local variables
!
!     *** variables for parallel process with MPI:
!
!     INODE   : rank of present node
!     NPROC   : number of nodes
!     PARLL   : flag to denote run as parallel (.TRUE.) or not (.FALSE.)
!     SWCHAR  : MPI datatype for characters
!     SWDBLE  : MPI datatype for double precision
!     SWINT   : MPI datatype for integers
!     SWMAX   : MPI collective maximum operation
!     SWMIN   : MPI collective minimum operation
!     SWREAL  : MPI datatype for reals
!     SWSUM   : MPI collective summation
!
      INTEGER INODE, NPROC
      INTEGER SWCHAR, SWINT, SWREAL, SWDBLE
      INTEGER SWMAX, SWMIN, SWSUM
      LOGICAL PARLL
!
!     *** information related to global domain and subdomains
!
!     DWRKEX  : double precision work array to store data to be sent to or received from neighbour
!     IBLKAD  : administration array for subdomain interfaces
!               contents:
!               pos. 1                     number of neighbouring subdomains
!                                          =m
!               pos. 3*i-1                 number of i-th neighbour
!               pos. 3*i                   position of i-th neighbour with
!                                          respect to present subdomain
!               pos. 3*i+1                 pointer of i-th neighbour in
!                                          last part of this array
!               pos. 3*m+2                 number of overlapping unknowns
!                                          on subdomain interface
!               pos. 3*m+3 ... 3*m+2+n     position of unknown in array
!                                          to be sent to neighbour
!               pos. 3*m+3+n ... 3*m+2*n+2 position of unknown in array
!                                          to be received from neighbour
!     IMRGE   : indicator merging output files
!               = 0; do not merge the individual output files
!               = 1; merge the individual output files into a single file
!                    and then delete all the output files
!     INCORN  : logical indicating whether a grid point belongs to a corner
!               of subdomain (=.TRUE.) or not (=.FALSE.)
!     IWEIG   : weights to determine load per part
!     IWRKEX  : integer work array to store data to be sent to or received from neighbour
!     KGRPGL  : index table containing the address of each (active) grid point in global domain
!               =1; not active grid point
!               >1; active grid point
!     KPART   : grid partitioning method
!               =-1; automatically setting based on computation mode
!               = 1; stripwise manner
!               = 2; orthogonal recursive bisection
!               = 3; stripwise partitioning along x-axis
!               = 4; stripwise partitioning along y-axis
!     LMXF    : logical indicating whether first x-point of subdomain equals
!               first x-point of global domain (=.TRUE.) or not (=.FALSE.)
!     LMXL    : logical indicating whether last x-point of subdomain equals
!               last x-point of global domain (=.TRUE.) or not (=.FALSE.)
!     LMYF    : logical indicating whether first y-point of subdomain equals
!               first y-point of global domain (=.TRUE.) or not (=.FALSE.)
!     LMYL    : logical indicating whether last y-point of subdomain equals
!               last y-point of global domain (=.TRUE.) or not (=.FALSE.)
!     LORB    : logical indicating which partition method will be carried out:
!               true, in case of ORB
!               false, in case of stripwise manner
!     MCGRDGL : number of wet grid points in global computational grid
!     MXCGL   : number of grid points in x-direction in global
!               computational grid
!     MXF     : first index w.r.t. global grid in x-direction
!     MXL     : last index w.r.t. global grid in x-direction
!     MYCGL   : number of grid points in y-direction in global
!               computational grid
!     MYF     : first index w.r.t. global grid in y-direction
!     MYL     : last index w.r.t. global grid in y-direction
!     NBGGL   : number of boundary grid points in global domain
!     WRKEXG  : work array to store data to be sent to or received from neighbour
!     XCLMAX  : maximum x-coordinate in subdomain
!     XCLMIN  : minimum x-coordinate in subdomain
!     YCLMAX  : maximum y-coordinate in subdomain
!     YCLMIN  : minimum y-coordinate in subdomain
!     XGRDGL  : x-coordinate of computational grid in global domain
!     YGRDGL  : y-coordinate of computational grid in global domain
!
      INTEGER IMRGE, KPART
      INTEGER MCGRDGL, MXCGL, MYCGL
      INTEGER MXF, MXL, MYF, MYL
      INTEGER NBGGL
      REAL    XCLMAX, XCLMIN, YCLMAX, YCLMIN

      LOGICAL LMXF, LMXL, LMYF, LMYL
      LOGICAL LORB

      INTEGER, SAVE, ALLOCATABLE :: IBLKAD(:)
      INTEGER, SAVE, ALLOCATABLE :: IWEIG(:)
      INTEGER, SAVE, ALLOCATABLE :: IWRKEX(:)
      INTEGER, SAVE, ALLOCATABLE :: KGRPGL(:,:)
      REAL   , SAVE, ALLOCATABLE :: XGRDGL(:,:), YGRDGL(:,:)
      REAL   , SAVE, ALLOCATABLE :: WRKEXG(:)
      REAL*8 , SAVE, ALLOCATABLE :: DWRKEX(:)
      LOGICAL, SAVE, ALLOCATABLE :: INCORN(:,:)
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
      END MODULE M_PARALL
!MPI!/impi
!MPI!/impi      MODULE MPI
!MPI!/impi      INCLUDE 'mpif.h'
!MPI!/impi      END MODULE MPI
