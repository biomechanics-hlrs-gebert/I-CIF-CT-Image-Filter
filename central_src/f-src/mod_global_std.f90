!------------------------------------------------------------------------------
! MODULE: global standards
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @description: 
!> Module containing all recurring definitions of kinds and numbers.
!------------------------------------------------------------------------------
MODULE global_std

IMPLICIT NONE

! Debugging
CHARACTER(LEN=*), PARAMETER :: out_amount  = "DEBUG" ! "PRODUCTION" 

! General constants
INTEGER, PARAMETER :: sik = 2   ! INTEGER Kind
INTEGER, PARAMETER :: ik  = 8   ! INTEGER Kind
INTEGER, PARAMETER :: pik = 4   ! PETSc INTEGER Kind; Compile with corresponding mpi.
INTEGER, PARAMETER :: mik = 4   ! MPI INTEGER Kind; Compile with corresponding mpi.
INTEGER, PARAMETER :: rk  = 8   ! Real Kind
INTEGER, PARAMETER :: mcl = 512 ! Max   character length
INTEGER, PARAMETER :: hcl = 256 ! Half  character length
INTEGER, PARAMETER :: scl = 64  ! Short character length

!-- File handles, debug_lvl and suffix
INTEGER(KIND=ik), PARAMETER :: timer_level = 3 ! 1 ! 2
INTEGER(KIND=ik), PARAMETER :: dbg_lvl = 1

CHARACTER(LEN=mcl) :: mssg = ''

! Dynamically assigned, must fit to program
INTEGER(KIND=ik)            :: std_out
INTEGER(KIND=ik), PARAMETER :: std_in  = 5
INTEGER(KIND=ik), PARAMETER :: std_err = 0

! Unsusual, as you may want to see the error messages asap in your root directory.
INTEGER(KIND=ik), PARAMETER :: fhsterr  = 10

!-- StdOut Characters
CHARACTER(len=5) :: creturn = achar(13)

! Provide versioning information for transparent data tracking
INCLUDE 'include_f90/revision_meta.f90'

END MODULE global_std