!------------------------------------------------------------------------------
! MODULE: math
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @description: 
!> Module containing recurring math options
!------------------------------------------------------------------------------
MODULE math

USE ISO_FORTRAN_ENV
USE global_std
USE strings

IMPLICIT NONE
REAL(KIND=rk), PARAMETER :: num_zero   = 1.E-9
REAL(KIND=rk), PARAMETER :: sq2        = sqrt(2._rk)
REAL(KIND=rk), PARAMETER :: pi         = 4.D0*DATAN(1.D0) !acos(-1._rk)
REAL(KIND=rk), PARAMETER :: inv180     = 1._rk/180._rk
REAL(KIND=rk), PARAMETER :: pi_div_180 = acos(-1._rk)/180._rk

!-- Higher dimensional numbers
TYPE Quaternion
   REAL (KIND=rk) :: w,x,y,z
END TYPE Quaternion

INTERFACE zero_thres
   MODULE PROCEDURE zerothres_num
   MODULE PROCEDURE zerothres_OnD
   MODULE PROCEDURE zerothres_TwD
   MODULE PROCEDURE zerothres_ThD
END INTERFACE zero_thres

CONTAINS

!------------------------------------------------------------------------------
! FUNCTION: rot_x
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @param[in] alpha Angle
!> @return[out] aa Output transformation matrix
!------------------------------------------------------------------------------  
Function rot_x(alpha) Result(aa)

    Real(kind=rk), intent(in) :: alpha
    Real(kind=rk), Dimension(3,3) :: aa

    aa(1,:) = [ 1._rk ,   0._rk    ,   0._rk     ]
    aa(2,:) = [ 0._rk , cos(alpha) , -sin(alpha) ]
    aa(3,:) = [ 0._rk , sin(alpha) ,  cos(alpha) ]

End Function rot_x

!------------------------------------------------------------------------------
! FUNCTION: rot_y
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @param[in] alpha Angle
!> @return[out] aa Output transformation matrix
!------------------------------------------------------------------------------  
Function rot_y(alpha) Result(aa)

    Real(kind=rk), intent(in) :: alpha
    Real(kind=rk), Dimension(3,3) :: aa

    aa(1,:) = [ cos(alpha), 0._rk,  sin(alpha) ]
    aa(2,:) = [   0._rk   , 1._rk,   0._rk     ]
    aa(3,:) = [-sin(alpha), 0._rk,  cos(alpha) ]

End Function rot_y

!------------------------------------------------------------------------------
! FUNCTION: rot_z
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @param[in] alpha Angle
!> @return[out] aa Output transformation matrix
!------------------------------------------------------------------------------  
Function rot_z(alpha) Result(aa)

    Real(kind=rk), intent(in) :: alpha
    Real(kind=rk), Dimension(3,3) :: aa

    aa(1,:) = [ cos(alpha), -sin(alpha), 0._rk ]
    aa(2,:) = [ sin(alpha),  cos(alpha), 0._rk ]
    aa(3,:) = [   0._rk   ,   0._rk    , 1._rk ]

End Function rot_z

!------------------------------------------------------------------------------
! FUNCTION: rot_alg
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> last edited : on  05.03.2009
!
!> @param[in]  axis Axis
!> @param[in]  angle Angle
!> @return[out] rr Output transformation matrix
!------------------------------------------------------------------------------  
Function rot_alg(axis, angle) Result(rr)

    real(kind=rk), dimension(3), Intent(In) :: axis
    real(kind=rk)              , Intent(In) :: angle
    real(kind=rk), dimension(3, 3) :: rr

    !** normalize rotation axis ***********************************************
    !axis = axis / sqrt(sum(axis*axis))

    !** Setup transformation matrix *******************************************
    rr(1,1) = cos(angle) + axis(1)*axis(1)* (1._8 - cos(angle))
    rr(1,2) = axis(1)*axis(2)* (1._8 - cos(angle)) - axis(3) * sin(angle)
    rr(1,3) = axis(1)*axis(3)* (1._8 - cos(angle)) + axis(2) * sin(angle)

    rr(2,1) = axis(2)*axis(1)* (1._8 - cos(angle)) + axis(3) * sin(angle)
    rr(2,2) = cos(angle) + axis(2)*axis(2)* (1._8 - cos(angle))
    rr(2,3) = axis(2)*axis(3)* (1._8 - cos(angle)) - axis(1) * sin(angle)

    rr(3,1) = axis(3)*axis(1)* (1._8 - cos(angle)) - axis(2) * sin(angle)
    rr(3,2) = axis(3)*axis(2)* (1._8 - cos(angle)) + axis(1) * sin(angle)
    rr(3,3) = cos(angle) + axis(3)*axis(3)* (1._8 - cos(angle))

End Function rot_alg

!------------------------------------------------------------------------------
! FUNCTION: tra_R6
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @brief
!> Transformation matrix for R6x6
!
!> @param[in]  aa Input matrix
!> @return[out] BB Output matrix
!------------------------------------------------------------------------------  
Function tra_R6(aa) Result(BB)

    Real(kind=rk), Dimension(3,3), intent(in) :: aa
    Real(kind=rk), Dimension(6,6) :: BB

    BB(1,:) = [ aa(1,1)**2 , aa(1,2)**2 , aa(1,3)**2 , sq2*aa(1,1)*aa(1,2) , sq2*aa(1,1)*aa(1,3), sq2*aa(1,2)*aa(1,3) ]
    BB(2,:) = [ aa(2,1)**2 , aa(2,2)**2 , aa(2,3)**2 , sq2*aa(2,1)*aa(2,2) , sq2*aa(2,1)*aa(2,3), sq2*aa(2,2)*aa(2,3) ]
    BB(3,:) = [ aa(3,1)**2 , aa(3,2)**2 , aa(3,3)**2 , sq2*aa(3,1)*aa(3,2) , sq2*aa(3,1)*aa(3,3), sq2*aa(3,2)*aa(3,3) ]

    BB(4,:) = [ sq2*aa(2,1)*aa(1,1) , sq2*aa(2,2)*aa(1,2) , sq2*aa(2,3)*aa(1,3) , &
        aa(2,1)*aa(1,2)+aa(2,2)*aa(1,1) , aa(2,1)*aa(1,3)+aa(2,3)*aa(1,1) , aa(2,2)*aa(1,3)+aa(2,3)*aa(1,2) ]
    BB(5,:) = [ sq2*aa(1,1)*aa(3,1) , sq2*aa(1,2)*aa(3,2) , sq2*aa(1,3)*aa(3,3) ,  &
        aa(1,1)*aa(3,2)+aa(1,2)*aa(3,1) , aa(1,1)*aa(3,3)+aa(1,3)*aa(3,1) , aa(1,2)*aa(3,3)+aa(1,3)*aa(3,2) ]
    BB(6,:) = [ sq2*aa(2,1)*aa(3,1) , sq2*aa(2,2)*aa(3,2) , sq2*aa(2,3)*aa(3,3) ,  &
        aa(2,1)*aa(3,2)+aa(2,2)*aa(3,1) , aa(2,1)*aa(3,3)+aa(2,3)*aa(3,1) , aa(2,2)*aa(3,3)+aa(2,3)*aa(3,2) ]

End Function tra_R6


!------------------------------------------------------------------------------
! SUBROUTINE: check_sym
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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

INTEGER(KIND=ik),              INTENT(IN) :: fh
REAL(KIND=rk), DIMENSION(:,:), INTENT(IN) :: mi
CHARACTER(LEN=*),              INTENT(IN) , OPTIONAL :: name
REAL(KIND=rk), DIMENSION(:,:), INTENT(OUT), OPTIONAL :: mo
REAL(KIND=rk)                , INTENT(OUT), OPTIONAL :: sym_out

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
! q counts the rows/columns to suppress a half of the matrix
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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

thres_u = num_zero
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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

thres_u = num_zero
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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

thres_u = num_zero
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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

thres_u = num_zero
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