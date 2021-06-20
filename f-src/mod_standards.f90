!---------------------------------------------------------------------------------------------------
! mod_standards.f90
! Module with standard precision definitions
!
! Author Johannes Gebert
! Date 04.01.2021

MODULE standards

IMPLICIT NONE

INTEGER           , PARAMETER  :: sik        = 2         ! INTEGER Kind
INTEGER           , PARAMETER  :: ik         = 4         ! INTEGER Kind
INTEGER           , PARAMETER  :: rk         = 8         ! Real Kind
INTEGER           , PARAMETER  :: mcl        = 512       ! Maximal character length
INTEGER           , PARAMETER  :: scl        = 64        ! Short character length

!-- StdOut Characters
CHARACTER(len=5)               :: creturn    = achar(13)

!-- MPI-specific kinds
INTEGER           , PARAMETER  :: mik        = 4         ! Kind MPI was compiled with

!-- Math standards
REAL(KIND=rk)     , PARAMETER  :: is_zero    = 1.E-9_rk
REAL(KIND=rk)     , PARAMETER  :: sq2        = sqrt(2._rk)
REAL(KIND=rk)     , PARAMETER  :: pi         = 4.D0*DATAN(1.D0)
REAL(KIND=rk)     , PARAMETER  :: inv180     = 1._rk/180._rk
REAL(KIND=rk)     , PARAMETER  :: pi_div_180 = acos(-1._rk)/180._rk

!-- Higher dimensional numbers
TYPE Quaternion
   REAL (KIND=rk)              :: w,x,y,z
END TYPE Quaternion

!-- Format Files
CHARACTER(len=100)             :: std_lnbrk  = &
"----------------------------------------------------------------------------------------------------"

CONTAINS

!-- Mechanical Standard Object - maybe to put into another module sooner or later

FUNCTION tensor_iso(E, v) RESULT (t_iso)

   REAL (KIND=rk)                 :: E, v
   REAL (KIND=rk), DIMENSION(6,6) :: t_iso

   t_iso(:,1)=(/ 1._rk/E,  -v/E  , -v/E,          .0_rk, .0_rk, .0_rk /)
   t_iso(:,2)=(/ -v/E , 1._rk/E  , -v/E,          .0_rk, .0_rk, .0_rk /)
   t_iso(:,3)=(/ -v/E ,  -v/E,  1._rk/E,          .0_rk, .0_rk, .0_rk /)
   t_iso(:,4)=(/ .0_rk, .0_rk, .0_rk, 2._rk*(1._rk+v)/E, .0_rk, .0_rk /)
   t_iso(:,5)=(/ .0_rk, .0_rk, .0_rk, .0_rk, 2._rk*(1._rk+v)/E, .0_rk /)
   t_iso(:,6)=(/ .0_rk, .0_rk, .0_rk, .0_rk, .0_rk, 2._rk*(1._rk+v)/E /)

END FUNCTION tensor_iso

END MODULE standards