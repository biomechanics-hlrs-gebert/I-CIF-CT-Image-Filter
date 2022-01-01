!------------------------------------------------------------------------------
! MODULE: math
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @description: 
!> Module containing recurring math options
!------------------------------------------------------------------------------
MODULE math

USE global_std
USE strings

IMPLICIT NONE

REAL(KIND=rk), PARAMETER :: is_zero    = 1.E-9_rk
REAL(KIND=rk), PARAMETER :: sq2        = sqrt(2._rk)
REAL(KIND=rk), PARAMETER :: pi         = 4.D0*DATAN(1.D0) !acos(-1._rk)
REAL(KIND=rk), PARAMETER :: inv180     = 1._rk/180._rk
REAL(KIND=rk), PARAMETER :: pi_div_180 = acos(-1._rk)/180._rk

!-- Higher dimensional numbers
TYPE Quaternion
   REAL (KIND=rk)            :: w,x,y,z
END TYPE Quaternion

Interface zero_thres
   Module Procedure zerothres_num
   Module Procedure zerothres_OnD
   Module Procedure zerothres_TwD
   Module Procedure zerothres_ThD
End Interface zero_thres

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: check_sym
!------------------------------------------------------------------------------  
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Check the symmetry of an arbitrarily sized square matrix.
!
!> @Description
!> Evaluates the symmetry of the matrix dependent of the quotient of 
!> the L2-Norm of the subtracted minor diagonals devided by the 
!> L2-norm of the input matrix
!
!> @param[in]  fh File handle to write to
!> @param[in]  mi Input Matrix
!> @param[in]  name Name of the matrix to evaluate
!> @param[out] mo Output Matrix
!> @param[out] sym_out Sum of differences of all ii,jj entries
!------------------------------------------------------------------------------  
SUBROUTINE check_sym(fh, mi, name, mo, sym_out)

INTEGER(KIND=ik),                 INTENT(IN) :: fh
REAL(KIND=rk),    DIMENSION(:,:), INTENT(IN) :: mi
CHARACTER(LEN=*),                 INTENT(IN) , OPTIONAL :: name
REAL(KIND=rk)   , DIMENSION(:,:), INTENT(OUT), OPTIONAL :: mo
REAL(KIND=rk)                   , INTENT(OUT), OPTIONAL :: sym_out

REAL(KIND=rk), DIMENSION(:,:), ALLOCATABLE :: mat 
REAL(KIND=rk) :: sym 
INTEGER, DIMENSION(2) :: lb, ub
INTEGER(KIND=ik) :: ii, jj, q
CHARACTER(LEN=scl) :: name_u, sym_out_str

name_u = ''

lb = LBOUND(mi)
ub = UBOUND(mi)

ALLOCATE(mat(lb(1), ub(2)))
mat = 0._rk

!------------------------------------------------------------------------------
! Calculate the differences to get the information of symmetry
!------------------------------------------------------------------------------
sym = 0._rk

Do jj = lb(2), ub(2)
    DO ii = lb(1), ub(1)
        sym = sym + ABS(mi(ii,jj) -  mi(jj,ii))
    End DO
End Do

IF(PRESENT(sym_out)) sym_out = sym

!------------------------------------------------------------------------------
! Write matrix out with zeros to show check_sym
!------------------------------------------------------------------------------
mat = mi
IF(sym <= 10E-06) THEN
    q = 1
    DO jj=lb(2), ub(2)-1 ! columns
        q = q + 1 ! begins with 2
        DO ii=lb(1), ub(2) ! rows
            mat (ii,jj) = 0._rk
        END DO
    END DO
END IF

IF(PRESENT(mo)) mo = mat
DEALLOCATE(mat)

IF(PRESENT(name)) name_u = ' '//TRIM(ADJUSTL(name))

!------------------------------------------------------------------------------
! Format output
!------------------------------------------------------------------------------
WRITE(sym_out_str, '(F40.20)') sym

CALL trimzero(sym_out_str)

WRITE(fh, '(4A)')"-- check_sym",TRIM(name_u),": ", TRIM(ADJUSTL(sym_out_str))

END SUBROUTINE check_sym

!------------------------------------------------------------------------------
! SUBROUTINE: zerothres_num
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Sets a scalar=0 in case it is less than 10^(-11) by default
!
!> @param[in] num Scalar input
!> @param[in] threshold to change to 0
!------------------------------------------------------------------------------ 
SUBROUTINE zerothres_num(num, thres)

REAL(KIND=rk), INTENT(INOUT) :: num 
REAL(KIND=rk), INTENT(IN), OPTIONAL :: thres

REAL(KIND=rk) :: thres_u

thres_u = 1E-11
IF(PRESENT(thres)) thres_u = thres

IF (num >= 0._rk) THEN
    IF (num <=  thres_u) num = 0._rk
ELSE
    IF (num >= -thres_u) num = 0._rk
END IF

END SUBROUTINE zerothres_num


!------------------------------------------------------------------------------
! SUBROUTINE: zerothres_OnD
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Sets a vector=0 element wise in case it is less than 10^(-11) by default
!
!> @param[in] oneD Scalar input
!> @param[in] threshold to change to 0
!------------------------------------------------------------------------------ 
SUBROUTINE zerothres_OnD(oneD, thres)

REAL(KIND=rk), DIMENSION(:), INTENT(INOUT) :: oneD 
REAL(KIND=rk), INTENT(IN), OPTIONAL :: thres

REAL(KIND=rk) :: thres_u
INTEGER(KIND=ik) :: ii

thres_u = 1E-11
IF(PRESENT(thres)) thres_u = thres

DO ii=1, SIZE(oneD)
    IF (oneD(ii) >= 0._rk) THEN
        IF (oneD(ii) <=  thres_u) oneD(ii) = 0._rk
    ELSE
        IF (oneD(ii) >= -thres_u) oneD(ii) = 0._rk
    END IF
END DO
END SUBROUTINE zerothres_OnD


!------------------------------------------------------------------------------
! SUBROUTINE: zerothres_TwD
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Sets an array=0 element wise in case it is less than 10^(-11) by default
!
!> @param[in] oneD Scalar input
!> @param[in] threshold to change to 0
!------------------------------------------------------------------------------ 
SUBROUTINE zerothres_TwD(TwD, thres)

REAL(KIND=rk), DIMENSION(:,:), INTENT(INOUT) :: TwD 
REAL(KIND=rk), INTENT(IN), OPTIONAL :: thres

REAL(KIND=rk) :: thres_u
INTEGER(KIND=ik) :: ii, jj

thres_u = 1E-11
IF(PRESENT(thres)) thres_u = thres

DO jj=1, SIZE(TwD, 2)
DO ii=1, SIZE(TwD, 1)
    IF (TwD(ii, jj) >= 0._rk) THEN
        IF (TwD(ii, jj) <=  thres_u) TwD(ii, jj) = 0._rk
    ELSE
        IF (TwD(ii, jj) >= -thres_u) TwD(ii, jj) = 0._rk
    END IF
END DO
END DO
END SUBROUTINE zerothres_TwD


!------------------------------------------------------------------------------
! SUBROUTINE: zerothres_ThD
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Sets an array=0 element wise in case it is less than 10^(-11) by default
!
!> @param[in] oneD Scalar input
!> @param[in] threshold to change to 0
!------------------------------------------------------------------------------ 
SUBROUTINE zerothres_ThD(ThD, thres)

REAL(KIND=rk), DIMENSION(:, :, :), INTENT(INOUT) :: ThD 
REAL(KIND=rk), INTENT(IN), OPTIONAL :: thres

REAL(KIND=rk) :: thres_u
INTEGER(KIND=ik) :: ii, jj, kk

thres_u = 1E-11
IF(PRESENT(thres)) thres_u = thres

DO kk=1, SIZE(ThD, 3)
DO jj=1, SIZE(ThD, 2)
DO ii=1, SIZE(ThD, 1)
    IF (ThD(ii, jj, kk) >= 0._rk) THEN
        IF (ThD(ii, jj, kk) <=  thres_u) ThD(ii, jj, kk) = 0._rk
    ELSE
        IF (ThD(ii, jj, kk) >= -thres_u) ThD(ii, jj, kk) = 0._rk
    END IF
END DO
END DO
END DO
END SUBROUTINE zerothres_ThD


END MODULE math