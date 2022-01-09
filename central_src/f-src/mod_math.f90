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
REAL(KIND=rk), PARAMETER :: pi         = 4.D0*DATAN(1.D0)       !acos(-1._rk)
REAL(KIND=rk), PARAMETER :: pihalf     = 4.D0*DATAN(1.D0)/2._rk !acos(-1._rk)
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
! SUBROUTINE: transpose_mat
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Spatially transpose a R6x6 matrix (2nd rank tensor).
!
!> @param[in] tensor_in Input tensor
!> @param[in] pos_in Requested combination of euler angles.
!> @param[out] tensor_out Output tensor
!------------------------------------------------------------------------------
SUBROUTINE transpose_mat (tensor_in, pos_in, tensor_out)

     REAL(KIND=rk), DIMENSION(6,6), INTENT(IN)  :: tensor_in
     REAL(KIND=rk), DIMENSION(3)  , INTENT(IN)  :: pos_in
     REAL(KIND=rk), DIMENSION(6,6), INTENT(OUT) :: tensor_out

     REAL(KIND=rk)                 :: alpha, phi, eta
     REAL(KIND=rk), DIMENSION(3)   :: n
     REAL(KIND=rk), DIMENSION(3,3) :: aa
     REAL(KIND=rk), DIMENSION(6,6) :: BB

     !------------------------------------------------------------------------------
     !  Degrees as input, radian as output to sin/cos
     !------------------------------------------------------------------------------
     alpha = REAL(pos_in(1)) * pi / 180._rk
     phi   = REAL(pos_in(2)) * pi / 180._rk
     eta   = REAL(pos_in(3)) * pi / 180._rk

     n = [ COS(phi)*SIN(eta), SIN(phi)*SIN(eta), COS(eta) ]

     n = n / SQRT(SUM(n*n))

     aa = rot_alg(n, alpha)

     BB = tra_R6(aa)

     tensor_out = MATMUL(MATMUL(TRANSPOSE(BB), tensor_in), BB)

END SUBROUTINE transpose_mat

!------------------------------------------------------------------------------
! FUNCTION: rot_x
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @param[in] angle Angle
!> @return[out] aa Output transformation matrix
!------------------------------------------------------------------------------  
Function rot_x(angle) Result(aa)

    Real(kind=rk), intent(in) :: angle
    Real(kind=rk), Dimension(3,3) :: aa

    aa(1,:) = [ 1._rk ,   0._rk    ,   0._rk     ]
    aa(2,:) = [ 0._rk , cos(angle) , -sin(angle) ]
    aa(3,:) = [ 0._rk , sin(angle) ,  cos(angle) ]

End Function rot_x

!------------------------------------------------------------------------------
! FUNCTION: rot_y
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @param[in] angle Angle
!> @return[out] aa Output transformation matrix
!------------------------------------------------------------------------------  
Function rot_y(angle) Result(aa)

    Real(kind=rk), intent(in) :: angle
    Real(kind=rk), Dimension(3,3) :: aa

    aa(1,:) = [ cos(angle), 0._rk,  sin(angle) ]
    aa(2,:) = [   0._rk   , 1._rk,   0._rk     ]
    aa(3,:) = [-sin(angle), 0._rk,  cos(angle) ]

End Function rot_y

!------------------------------------------------------------------------------
! FUNCTION: rot_z
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @param[in] angle Angle
!> @return[out] aa Output transformation matrix
!------------------------------------------------------------------------------  
Function rot_z(angle) Result(aa)

    Real(kind=rk), intent(in) :: angle
    Real(kind=rk), Dimension(3,3) :: aa

    aa(1,:) = [ cos(angle), -sin(angle), 0._rk ]
    aa(2,:) = [ sin(angle),  cos(angle), 0._rk ]
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
!> Check the symmetry of an arbitrarily sized matrix.
!
!> @Description
!> Average percentage of deviation of minor diagonals.
!
!> @param[in]  matin Input Matrix
!> @param[out] matout Output Matrix
!> @param[out] sym Sum of differences of all ii,jj entries
!------------------------------------------------------------------------------  
SUBROUTINE check_sym(matin, sym, matout)

REAL(KIND=rk), DIMENSION(:,:), INTENT(IN)  :: matin
REAL(KIND=rk)                , INTENT(OUT) :: sym
REAL(KIND=rk), DIMENSION(:,:), INTENT(OUT), OPTIONAL :: matout

! INTEGER, DIMENSION(2) :: lb, ub
INTEGER(KIND=ik) :: ii, jj
REAL(KIND=rk) :: cummu, entry_counter

! lb = LBOUND(matin)
! ub = UBOUND(matin)

! ALLOCATE(matout(lb(1), ub(2)))
! matout = 0._rk

!------------------------------------------------------------------------------
! Calculate the differences to get the information of symmetry
! Earlier version...
!------------------------------------------------------------------------------
! sym = 0._rk

! DO jj = lb(2), ub(2)
!     DO ii = lb(1), ub(1)
!         sym = sym + ABS(matin(ii,jj) -  matin(jj,ii))
    ! END DO
! END DO

cummu = 0._rk
ii=1_ik
entry_counter = 0._rk
DO WHILE (ii < SIZE(matin, DIM=1))
 
    jj=2_ik
    DO WHILE (jj <= SIZE(matin, DIM=2))
        cummu = cummu + (matin(ii,jj) /  matin(jj,ii))  

        !------------------------------------------------------------------------------
        ! How many entries are averaged?
        !------------------------------------------------------------------------------
        entry_counter = entry_counter + 1._rk     

        jj = jj + 1_ik
    END DO

    ii = ii + 1_ik
END DO

!------------------------------------------------------------------------------
! 1 - sym quotient to compare to 0
!------------------------------------------------------------------------------
sym = 1._rk - (cummu / entry_counter)

!------------------------------------------------------------------------------
! Write matout
!------------------------------------------------------------------------------
IF(PRESENT(matout)) matout = matin 

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