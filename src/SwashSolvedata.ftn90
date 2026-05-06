module SwashSolvedata
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
!    1.00, April 2010: New Module
!
!   Purpose
!
!   Module containing data for the linear solvers
!
!   Method
!
!   Data with respect to solution of system of equations
!
!   Modules used
!
!   none
!
    implicit none
!
!   Module parameters
!
!   ---
!
!   Module variables
!
    integer                                     :: iamout ! control parameter indicating the amount of output required
                                                          ! 0: no output
                                                          ! 1: only fatal errors will be printed
                                                          ! 2: gives output concerning the iteration process
                                                          ! 3: additional information about the iteration is printed
    integer                                     :: icond  ! control parameter indicating the choice of preconditioner for the BiCGSTAB solver
                                                          ! 1: ILUD as left-right or split preconditioner
                                                          ! 2: ILUD as right preconditioner
                                                          ! 3: ILU  as right preconditioner
    !
    integer, dimension(:)   , save, allocatable :: di     ! array containing the column indices of the diagonal elements in each row (CSR format used in bicgstab)
    integer, dimension(:)   , save, allocatable :: dis    ! array containing the column indices of the diagonal elements in each row (CSR format used in pcg)
    integer, dimension(:)   , save, allocatable :: ia     ! integer array containing the pointers to the beginning of each row in arrays ax and ja (CSR format used in bicgstab)
    integer, dimension(:)   , save, allocatable :: ias    ! integer array containing the pointers to the beginning of each row in arrays ax and ja (CSR format used in pcg)
    integer, dimension(:)   , save, allocatable :: iplu   ! integer array containing the pointers to the entries in L and U in each row (used in pcgu/bicgstabu)
    integer, dimension(:)   , save, allocatable :: iplu3  ! integer array containing the pointers to the entries in L and U in each row (used in bicgstab3)
    integer, dimension(:,:) , save, allocatable :: irw    ! integer array containing the order of matrix columns per row (used in bicgstab)
    integer, dimension(:,:) , save, allocatable :: irws   ! integer array containing the order of matrix columns per row (used in pcg)
    integer, dimension(:)   , save, allocatable :: ja     ! integer array containing the column indices of matrix elements stored in ax (CSR format used in bicgstab)
    integer, dimension(:)   , save, allocatable :: jas    ! integer array containing the column indices of matrix elements stored in ax (CSR format used in pcg)
    !
    real*8, dimension(:,:)  , save, allocatable :: amata  ! adapted coefficient matrix used in nested Newton iteration method
    real  , dimension(:)    , save, allocatable :: ax2    ! matrix stored in CSR format (used in bicgstab)
    real*8, dimension(:)    , save, allocatable :: axs    ! matrix stored in CSR format (used in pcg)
    real  , dimension(:)    , save, allocatable :: axs2   ! matrix stored in CSR format (used in pcg)
    real*8, dimension(:)    , save, allocatable :: ba     ! adapted main diagonal of tri-diagonal matrix used in nested Newton iteration method
    real  , dimension(:,:)  , save, allocatable :: cmat   ! the matrix containing an ILU factorization for the SIP solver
    real  , dimension(:)    , save, allocatable :: diag   ! the diagonal elements for scaling the system in the CG method (single precision)
    real*8, dimension(:)    , save, allocatable :: diag2  ! the diagonal elements for scaling the system in the CG method (double precision)
    real  , dimension(:,:)  , save, allocatable :: diagb  ! the diagonal elements for scaling the system in the BiCGSTAB solver
    real*8, dimension(:)    , save, allocatable :: da     ! adapted right-hand side of system of equations used in nested Newton iteration method
    real  , dimension(:,:)  , save, allocatable :: p      ! direction vector in BiCGSTAB
    real*8, dimension(:)    , save, allocatable :: p0     ! diagonal matrix p at previous iteration in nested Newton iteration method
    real*8, dimension(:)    , save, allocatable :: p1     ! diagonal matrix p at current iteration in nested Newton iteration method
    real  , dimension(:,:,:), save, allocatable :: prec   ! the matrix containing a RILU factorization for the BiCGSTAB solver
    real  , dimension(:)    , save, allocatable :: prec2  ! the matrix containing an ILU factorization for the BiCGSTAB solver (used in pcgu/bicgstabu)
    real*8, dimension(:)    , save, allocatable :: prec3  ! the matrix containing an ILU factorization for the BiCGSTAB solver (used in pcgu2)
    real  , dimension(:)    , save, allocatable :: prec4  ! the matrix containing an ILU factorization for the BiCGSTAB solver (used in bicgstab3)
    real*8, dimension(:)    , save, allocatable :: q0     ! diagonal matrix q at previous iteration in nested Newton iteration method
    real*8, dimension(:)    , save, allocatable :: q1     ! diagonal matrix q at current iteration in nested Newton iteration method
    real  , dimension(:)    , save, allocatable :: r      ! auxiliary vector containing product of matrix and conjugate vector in the CG method (single precision)
    real*8, dimension(:)    , save, allocatable :: r2     ! auxiliary vector containing product of matrix and conjugate vector in the CG method (double precision)
    real  , dimension(:,:)  , save, allocatable :: res    ! the residual vector
    real  , dimension(:,:)  , save, allocatable :: res0   ! the quasi-residual vector
    real  , dimension(:)    , save, allocatable :: resd   ! the residual vector (2DH)
    real*8, dimension(:)    , save, allocatable :: resd2  ! the residual vector (double precision)
    real*8, dimension(:)    , save, allocatable :: rhsa   ! adapted right-hand side used in nested Newton iteration method
    real  , dimension(:)    , save, allocatable :: s      ! the conjugate vector acting as search direction in the CG method (single precision)
    real*8, dimension(:)    , save, allocatable :: s2     ! the conjugate vector acting as search direction in the CG method (double precision)
    real*8, dimension(:)    , save, allocatable :: sol0   ! solution vector from previous iteration used in Newton-type iteration method
    real*8, dimension(:)    , save, allocatable :: sol1   ! solution vector from previous iteration used in Newton-type iteration method
    real  , dimension(:)    , save, allocatable :: sol3   ! solution vector (used in bicgstab3)
    real*8, dimension(:)    , save, allocatable :: sola   ! solution vector used in nested Newton iteration method
    real  , dimension(:,:)  , save, allocatable :: t      ! auxiliary vector in BiCGSTAB
    real  , dimension(:,:)  , save, allocatable :: u      ! auxiliary vector in BiCGSTAB
    real  , dimension(:,:)  , save, allocatable :: v      ! auxiliary vector containing product of matrix and direction vector in BiCGSTAB
    real*8, dimension(:)    , save, allocatable :: vt     ! auxiliary variable to store a vector temporarily
    real  , dimension(:,:)  , save, allocatable :: vt2    ! another auxiliary variable to store a vector temporarily
    real  , dimension(:,:)  , save, allocatable :: w      ! auxiliary variable to store a vector temporarily
    real  , dimension(:)    , save, allocatable :: w2     ! auxiliary variable to store a vector temporarily (used in bicgstab)
    real  , dimension(:)    , save, allocatable :: z      ! the preconditioned residual in the CG method (single precision)
    real*8, dimension(:)    , save, allocatable :: z3     ! the preconditioned residual in the CG method (double precision)
    real  , dimension(:)    , save, allocatable :: z1     ! auxiliary vector in DAC solver (single precision)
    real  , dimension(:)    , save, allocatable :: z2     ! another auxiliary vector in DAC solver (single precision)
    real*8, dimension(:)    , save, allocatable :: z12    ! auxiliary vector in DAC solver (double precision)
    real*8, dimension(:)    , save, allocatable :: z22    ! another auxiliary vector in DAC solver (double precision)
!
!   Source text
!
end module SwashSolvedata
