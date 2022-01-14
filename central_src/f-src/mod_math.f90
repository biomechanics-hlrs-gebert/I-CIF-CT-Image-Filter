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
REAL(KIND=rk), PARAMETER :: pi_div_180 = acos(-1._rk)/180._rk   ! * deg = rad

!-- Higher dimensional numbers
TYPE quaternion
   REAL (KIND=rk) :: w,x,y,z
END TYPE quaternion

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
!------------------------------------------------------------------------------  
SUBROUTINE check_sym(matin, sym)

    REAL(KIND=rk), DIMENSION(:,:), INTENT(IN)  :: matin
    REAL(KIND=rk)                , INTENT(OUT) :: sym

    ! INTEGER, DIMENSION(2) :: lb, ub
    INTEGER(KIND=ik) :: ii, jj
    REAL(KIND=rk) :: cummu, entry_counter

    !------------------------------------------------------------------------------
    ! Calculate the differences to get the information of symmetry
    ! Earlier version...
    !------------------------------------------------------------------------------

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

!------------------------------------------------------------------------------
! FUNCTION: eul_2_un_quat
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Returns a unit Quaternion based on 3 Euler Angles [radians]
!
!> @description
!> euler_radians(1) -> x-axis -> "alpha"
!> euler_radians(2) -> y-axis -> "beta"
!> euler_radians(3) -> z-axis -> "gamma"
!
!> @param[in] euler_radians Euler radians
!> @return eul_2_un_quat Unit Quaternion
!------------------------------------------------------------------------------  
FUNCTION eul_2_un_quat (euler_radians)

    REAL (KIND=rk) , DIMENSION(3) :: euler_radians
    TYPE(Quaternion) :: eul_2_un_quat

    REAL (KIND=rk) :: ca, cb, cg, sa, sb, sg

    euler_radians = euler_radians * 0.5_rk ! angle/2 to get a quaternion!

    ca = COS(euler_radians(1))
    sa = SIN(euler_radians(1))
    cb = COS(euler_radians(2))
    sb = SIN(euler_radians(2))
    cg = COS(euler_radians(3))
    sg = SIN(euler_radians(3))

    eul_2_un_quat%w = ca * cb * cg + sa * sb * sg
    eul_2_un_quat%x = sa * cb * cg - ca * sb * sg
    eul_2_un_quat%y = ca * sb * cg + sa * cb * sg
    eul_2_un_quat%z = ca * cb * sg - sa * sb * cg

END FUNCTION eul_2_un_quat

!------------------------------------------------------------------------------
! FUNCTION: check_un_quat
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Check whether it is a uni quaternion
!
!> @param[in] quat Quaternion
!> @return check_un_quat True or false
!------------------------------------------------------------------------------  
FUNCTION check_un_quat (quat)

    TYPE(Quaternion) :: quat
    LOGICAL :: check_un_quat

    REAL(KIND=rk) :: rslt

    rslt = quat%w**2 + quat%x**2 + quat%y**2 + quat%z**2

    IF (rslt-1_rk <= num_zero) THEN
        check_un_quat = .TRUE.
    ELSE
        check_un_quat = .FALSE.
    END IF

END FUNCTION check_un_quat

!------------------------------------------------------------------------------
! FUNCTION: quat_norm
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion norm
!
!> @param[in] quat Quaternion
!> @return quat_norm Norm of the quaternion
!------------------------------------------------------------------------------  
FUNCTION quat_norm (quat)

    TYPE(Quaternion) :: quat
    REAL(KIND=rk) :: quat_norm

    quat_norm = quat%w**2 + quat%x**2 + quat%y**2 + quat%z**2

END FUNCTION quat_norm

!------------------------------------------------------------------------------
! FUNCTION: q_add
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Add a scalar to a quaternion. Also works with negative scalars.
!
!> @param[in] q1 Quaternion
!> @param[in] scalar Scalar
!> @return q_add Result of the addition
!------------------------------------------------------------------------------  
FUNCTION q_add(q1, scalar)

    TYPE(Quaternion) :: q1, q_add
    REAL(KIND=rk) :: scalar

    q_add%w = q1%w + scalar
    q_add%x = q1%x + scalar
    q_add%y = q1%y + scalar
    q_add%z = q1%z + scalar

END FUNCTION q_add

!------------------------------------------------------------------------------
! FUNCTION: q_mult
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion multiplication. Element wise, not to be confused with quat_prod!
!
!> @param[in] q1 Quaternion
!> @param[in] scalar Scalar
!> @return q_add Result of the multiplication
!------------------------------------------------------------------------------  
FUNCTION q_mult (q1, scalar)

    TYPE(Quaternion) :: q1, q_mult
    REAL(KIND=rk) :: scalar

    q_mult%w = q1%w * scalar
    q_mult%x = q1%x * scalar
    q_mult%y = q1%y * scalar
    q_mult%z = q1%z * scalar

END FUNCTION q_mult

!------------------------------------------------------------------------------
! FUNCTION: q_dvd
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion division. Element wise, not to be confused with quat_prod!
!
!> @param[in] q1 Quaternion
!> @param[in] scalar Scalar
!> @return q_add Result of the division
!------------------------------------------------------------------------------ 
FUNCTION q_dvd (q1, scalar)

    TYPE(Quaternion) :: q1, q_dvd
    REAL(KIND=rk) :: scalar
    
    q_dvd%w = q1%w / scalar
    q_dvd%x = q1%x / scalar
    q_dvd%y = q1%y / scalar
    q_dvd%z = q1%z / scalar

END FUNCTION q_dvd

!------------------------------------------------------------------------------
! FUNCTION: quat_prod
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion division. Element wise, not to be confused with quat_prod!
!
!> @param[in] q1 Quaternion
!> @param[in] q2 Quaternion
!> @return quat_prod Result of the quaternion multiplication
!------------------------------------------------------------------------------ 
FUNCTION quat_prod (q1, q2)

    TYPE(Quaternion) :: q1, q2, quat_prod

    quat_prod%w = q1%w*q2%w-q1%x*q2%x-q1%y*q2%y-q1%z*q2%z
    quat_prod%x = q1%w*q2%x+q1%x*q2%w+q1%y*q2%z-q1%z*q2%y
    quat_prod%y = q1%w*q2%y-q1%x*q2%z+q1%y*q2%w+q1%z*q2%x
    quat_prod%z = q1%w*q2%z+q1%x*q2%y-q1%y*q2%x+q1%z*q2%w

END FUNCTION quat_prod

!------------------------------------------------------------------------------
! FUNCTION: conjugate_quat
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Conjugate quaternion.
!
!> @param[in] quat Quaternion
!> @return conjugate_quat Conjugate quaternion
!------------------------------------------------------------------------------ 
FUNCTION conjugate_quat(quat)

    TYPE(Quaternion) :: quat, conjugate_quat

    conjugate_quat%w =   quat%w
    conjugate_quat%x = - quat%x
    conjugate_quat%y = - quat%y
    conjugate_quat%z = - quat%z

END FUNCTION conjugate_quat

!------------------------------------------------------------------------------
! FUNCTION: rotate_point
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Rotate a point with quaternions. Takes a unit quaternion as input
!
!> @param[in] quat Quaternion
!> @param[in] point Point to rotate with quaternion.
!> @return rotate_point Rotated point.
!------------------------------------------------------------------------------
FUNCTION rotate_point(quat, point)

    TYPE(Quaternion) :: quat
    REAL(KIND=rk), DIMENSION(3) :: point, rotate_point

    TYPE(Quaternion) :: pnt, rtt_pnt

    pnt%w = 0._rk
    pnt%x = point(1)
    pnt%y = point(2)
    pnt%z = point(3)

    rtt_pnt = quat_prod(quat_prod(quat, pnt), conjugate_quat(quat))

    rotate_point(1) = rtt_pnt%x
    rotate_point(2) = rtt_pnt%y
    rotate_point(3) = rtt_pnt%z

END FUNCTION rotate_point


!------------------------------------------------------------------------------
! FUNCTION: quat_2_rot_mat
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion to rotation matrix.
!
!> @param[in] q Quaternion
!> @return quat_2_rot_mat Rotation matrix
!------------------------------------------------------------------------------
FUNCTION quat_2_rot_mat(q)

    TYPE(Quaternion) :: q
    REAL(KIND=rk), DIMENSION(3,3) :: quat_2_rot_mat

    quat_2_rot_mat(1,:) = [ q%w**2 + q%x**2 - q%y**2 - q%z**2  , &
        2._rk*q%x*q%y - 2._rk*q%w*q%z , 2._rk*q%x*q%z - 2._rk*q%w*q%y  ]
    quat_2_rot_mat(2,:) = [ 2._rk*q%x*q%y + 2._rk*q%w*q%z , &
        q%w**2 - q%x**2 + q%y**2 + - q%z**2 , 2._rk*q%y*q%z - 2._rk*q%w*q%x ]
    quat_2_rot_mat(3,:) = [ 2._rk*q%x*q%z - 2._rk*q%w*q%y , &
        2._rk*q%y*q%z + 2._rk*q%w*q%x , q%w**2 - q%x**2 - q%y**2 + q%z**2   ]

END FUNCTION quat_2_rot_mat


!------------------------------------------------------------------------------
! FUNCTION: crpr
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Cross product of two vectors
!
!> @param[in] a Vector
!> @param[in] b Vector
!> @return crpr Cross product
!------------------------------------------------------------------------------
FUNCTION crpr(a, b)

    REAL(KIND=rk), DIMENSION(3) :: crpr,a,b

    crpr(1) = a(2) * b(3) - a(3) * b(2)
    crpr(2) = a(3) * b(1) - a(1) * b(3)
    crpr(3) = a(1) * b(2) - a(2) * b(1)

END FUNCTION crpr

END MODULE math