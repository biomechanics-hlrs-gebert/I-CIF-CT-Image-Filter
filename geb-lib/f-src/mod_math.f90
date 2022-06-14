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

REAL(rk), PARAMETER :: num_zero   = 1.E-9
REAL(rk), PARAMETER :: sq2        = sqrt(2._rk)
REAL(rk), PARAMETER :: pi         = 4.D0*DATAN(1.D0)       !acos(-1._rk)
REAL(rk), PARAMETER :: pihalf     = 4.D0*DATAN(1.D0)/2._rk !acos(-1._rk)
REAL(rk), PARAMETER :: inv180     = 1._rk/180._rk
REAL(rk), PARAMETER :: pi_div_180 = acos(-1._rk)/180._rk   ! * deg = rad
REAL(rk), PARAMETER :: div_180_pi = 180._rk/acos(-1._rk)   ! * rad = deg

!-- Higher dimensional numbers
TYPE duaternion
   REAL (rk) :: w,x,y,z
END TYPE duaternion
TYPE suaternion
   REAL (rk) :: w,x,y,z
END TYPE suaternion
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

REAL(rk), DIMENSION(6,6), INTENT(IN)  :: tensor_in
REAL(rk), DIMENSION(3)  , INTENT(IN)  :: pos_in
REAL(rk), DIMENSION(6,6), INTENT(OUT) :: tensor_out

REAL(rk)                 :: alpha, phi, eta
REAL(rk), DIMENSION(3)   :: n
REAL(rk), DIMENSION(3,3) :: aa
REAL(rk), DIMENSION(6,6) :: BB

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

REAL(rk), DIMENSION(:,:), INTENT(IN)  :: matin
REAL(rk), INTENT(OUT) :: sym

INTEGER(ik) :: ii, jj
REAL(rk) :: cummu, entry_counter

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

REAL(rk), INTENT(INOUT) :: num 
REAL(rk), INTENT(IN), OPTIONAL :: thres

REAL(rk) :: thres_u

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

REAL(rk), DIMENSION(:), INTENT(INOUT) :: oneD 
REAL(rk), INTENT(IN), OPTIONAL :: thres

REAL(rk) :: thres_u
INTEGER(ik) :: ii

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

REAL(rk), DIMENSION(:,:), INTENT(INOUT) :: TwD 
REAL(rk), INTENT(IN), OPTIONAL :: thres

REAL(rk) :: thres_u
INTEGER(ik) :: ii, jj

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

REAL(rk), DIMENSION(:, :, :), INTENT(INOUT) :: ThD 
REAL(rk), INTENT(IN), OPTIONAL :: thres

REAL(rk) :: thres_u
INTEGER(ik) :: ii, jj, kk

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
! FUNCTION: eul_2_un_quat_rk
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
!> @return eul_2_un_quat_rk Unit Quaternion
!------------------------------------------------------------------------------  
FUNCTION eul_2_un_quat_rk (euler_radians)

REAL(rk), DIMENSION(3) :: euler_radians, radians
TYPE(duaternion) :: eul_2_un_quat_rk

REAL (rk) :: ca, cb, cg, sa, sb, sg

radians = euler_radians * 0.5_rk ! angle/2 to get a quaternion!

ca = COS(radians(1))
sa = SIN(radians(1))
cb = COS(radians(2))
sb = SIN(radians(2))
cg = COS(radians(3))
sg = SIN(radians(3))

eul_2_un_quat_rk%w = ca * cb * cg + sa * sb * sg
eul_2_un_quat_rk%x = sa * cb * cg - ca * sb * sg
eul_2_un_quat_rk%y = ca * sb * cg + sa * cb * sg
eul_2_un_quat_rk%z = ca * cb * sg - sa * sb * cg

END FUNCTION eul_2_un_quat_rk

!------------------------------------------------------------------------------
! FUNCTION: check_un_quat_rk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Check whether it is a uni quaternion
!
!> @param[in] quat Quaternion
!> @return check_un_quat_rk True or false
!------------------------------------------------------------------------------  
FUNCTION check_un_quat_rk (quat)

TYPE(duaternion) :: quat
LOGICAL :: check_un_quat_rk

REAL(rk) :: rslt

rslt = quat%w**2 + quat%x**2 + quat%y**2 + quat%z**2

IF (rslt-1_rk <= num_zero) THEN
    check_un_quat_rk = .TRUE.
ELSE
    check_un_quat_rk = .FALSE.
END IF

END FUNCTION check_un_quat_rk

!------------------------------------------------------------------------------
! FUNCTION: quat_norm_rk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion norm
!
!> @param[in] quat Quaternion
!> @return quat_norm_rk Norm of the quaternion
!------------------------------------------------------------------------------  
FUNCTION quat_norm_rk (quat)

TYPE(duaternion) :: quat
REAL(rk) :: quat_norm_rk

quat_norm_rk = quat%w**2 + quat%x**2 + quat%y**2 + quat%z**2

END FUNCTION quat_norm_rk

!------------------------------------------------------------------------------
! FUNCTION: q_add_rk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Add a scalar to a quaternion. Also works with negative scalars.
!
!> @param[in] q1 Quaternion
!> @param[in] scalar Scalar
!> @return q_add_rk Result of the addition
!------------------------------------------------------------------------------  
FUNCTION q_add_rk(q1, scalar)

TYPE(duaternion) :: q1, q_add_rk
REAL(rk) :: scalar

q_add_rk%w = q1%w + scalar
q_add_rk%x = q1%x + scalar
q_add_rk%y = q1%y + scalar
q_add_rk%z = q1%z + scalar

END FUNCTION q_add_rk

!------------------------------------------------------------------------------
! FUNCTION: q_mult_rk
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
FUNCTION q_mult_rk (q1, scalar)

TYPE(duaternion) :: q1, q_mult_rk
REAL(rk) :: scalar

q_mult_rk%w = q1%w * scalar
q_mult_rk%x = q1%x * scalar
q_mult_rk%y = q1%y * scalar
q_mult_rk%z = q1%z * scalar

END FUNCTION q_mult_rk

!------------------------------------------------------------------------------
! FUNCTION: q_dvd_rk
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
FUNCTION q_dvd_rk (q1, scalar)

TYPE(duaternion) :: q1, q_dvd_rk
REAL(rk) :: scalar

q_dvd_rk%w = q1%w / scalar
q_dvd_rk%x = q1%x / scalar
q_dvd_rk%y = q1%y / scalar
q_dvd_rk%z = q1%z / scalar

END FUNCTION q_dvd_rk

!------------------------------------------------------------------------------
! FUNCTION: quat_prod_rk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion division. Element wise, not to be confused with quat_prod_rk!
!
!> @param[in] q1 Quaternion
!> @param[in] q2 Quaternion
!> @return quat_prod_rk Result of the quaternion multiplication
!------------------------------------------------------------------------------ 
FUNCTION quat_prod_rk (q1, q2)

TYPE(duaternion) :: q1, q2, quat_prod_rk

quat_prod_rk%w = q1%w*q2%w-q1%x*q2%x-q1%y*q2%y-q1%z*q2%z
quat_prod_rk%x = q1%w*q2%x+q1%x*q2%w+q1%y*q2%z-q1%z*q2%y
quat_prod_rk%y = q1%w*q2%y-q1%x*q2%z+q1%y*q2%w+q1%z*q2%x
quat_prod_rk%z = q1%w*q2%z+q1%x*q2%y-q1%y*q2%x+q1%z*q2%w

END FUNCTION quat_prod_rk

!------------------------------------------------------------------------------
! FUNCTION: conjugate_quat_rk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Conjugate quaternion.
!
!> @param[in] quat Quaternion
!> @return conjugate_quat_rk Conjugate quaternion
!------------------------------------------------------------------------------ 
FUNCTION conjugate_quat_rk(quat)

TYPE(duaternion) :: quat, conjugate_quat_rk

conjugate_quat_rk%w =   quat%w
conjugate_quat_rk%x = - quat%x
conjugate_quat_rk%y = - quat%y
conjugate_quat_rk%z = - quat%z

END FUNCTION conjugate_quat_rk

!------------------------------------------------------------------------------
! FUNCTION: rotate_point_rk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Rotate a point with quaternions. Takes a unit quaternion as input
!
!> @param[in] quat Quaternion
!> @param[in] point Point to rotate with quaternion.
!> @return rotate_point_rk Rotated point.
!------------------------------------------------------------------------------
FUNCTION rotate_point_rk(quat, point)

TYPE(duaternion) :: quat
REAL(rk), DIMENSION(3) :: point, rotate_point_rk

TYPE(duaternion) :: pnt, rtt_pnt

pnt%w = 0._rk
pnt%x = point(1)
pnt%y = point(2)
pnt%z = point(3)

rtt_pnt = quat_prod_rk(quat_prod_rk(quat, pnt), conjugate_quat_rk(quat))

rotate_point_rk(1) = rtt_pnt%x
rotate_point_rk(2) = rtt_pnt%y
rotate_point_rk(3) = rtt_pnt%z

END FUNCTION rotate_point_rk


!------------------------------------------------------------------------------
! FUNCTION: quat_2_rot_mat_rk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion to rotation matrix.
!
!> @param[in] q Quaternion
!> @return quat_2_rot_mat_rk Rotation matrix
!------------------------------------------------------------------------------
FUNCTION quat_2_rot_mat_rk(q)

TYPE(duaternion) :: q
REAL(rk), DIMENSION(3,3) :: quat_2_rot_mat_rk

quat_2_rot_mat_rk(1,:) = [ q%w**2 + q%x**2 - q%y**2 - q%z**2  , &
    2._rk*q%x*q%y - 2._rk*q%w*q%z , 2._rk*q%x*q%z - 2._rk*q%w*q%y  ]
quat_2_rot_mat_rk(2,:) = [ 2._rk*q%x*q%y + 2._rk*q%w*q%z , &
    q%w**2 - q%x**2 + q%y**2 + - q%z**2 , 2._rk*q%y*q%z - 2._rk*q%w*q%x ]
quat_2_rot_mat_rk(3,:) = [ 2._rk*q%x*q%z - 2._rk*q%w*q%y , &
    2._rk*q%y*q%z + 2._rk*q%w*q%x , q%w**2 - q%x**2 - q%y**2 + q%z**2   ]

END FUNCTION quat_2_rot_mat_rk


!------------------------------------------------------------------------------
! FUNCTION: eul_2_un_quat_sk
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
!> @return eul_2_un_quat_sk Unit Quaternion
!------------------------------------------------------------------------------  
FUNCTION eul_2_un_quat_sk (euler_radians)

REAL(sk) , DIMENSION(3) :: euler_radians, radians
TYPE(suaternion) :: eul_2_un_quat_sk

REAL (sk) :: ca, cb, cg, sa, sb, sg

radians = euler_radians * 0.5_rk ! angle/2 to get a quaternion!

ca = COS(radians(1))
sa = SIN(radians(1))
cb = COS(radians(2))
sb = SIN(radians(2))
cg = COS(radians(3))
sg = SIN(radians(3))

eul_2_un_quat_sk%w = ca * cb * cg + sa * sb * sg
eul_2_un_quat_sk%x = sa * cb * cg - ca * sb * sg
eul_2_un_quat_sk%y = ca * sb * cg + sa * cb * sg
eul_2_un_quat_sk%z = ca * cb * sg - sa * sb * cg

END FUNCTION eul_2_un_quat_sk

!------------------------------------------------------------------------------
! FUNCTION: check_un_quat_sk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Check whether it is a uni quaternion
!
!> @param[in] quat Quaternion
!> @return check_un_quat_sk True or false
!------------------------------------------------------------------------------  
FUNCTION check_un_quat_sk (quat)

TYPE(suaternion) :: quat
LOGICAL :: check_un_quat_sk

REAL(sk) :: rslt

rslt = quat%w**2 + quat%x**2 + quat%y**2 + quat%z**2

IF (rslt-1_rk <= num_zero) THEN
    check_un_quat_sk = .TRUE.
ELSE
    check_un_quat_sk = .FALSE.
END IF

END FUNCTION check_un_quat_sk

!------------------------------------------------------------------------------
! FUNCTION: quat_norm_sk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion norm
!
!> @param[in] quat Quaternion
!> @return quat_norm_sk Norm of the quaternion
!------------------------------------------------------------------------------  
FUNCTION quat_norm_sk (quat)

TYPE(suaternion) :: quat
REAL(sk) :: quat_norm_sk

quat_norm_sk = quat%w**2 + quat%x**2 + quat%y**2 + quat%z**2

END FUNCTION quat_norm_sk

!------------------------------------------------------------------------------
! FUNCTION: q_add_sk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Add a scalar to a quaternion. Also works with negative scalars.
!
!> @param[in] q1 Quaternion
!> @param[in] scalar Scalar
!> @return q_add_sk Result of the addition
!------------------------------------------------------------------------------  
FUNCTION q_add_sk(q1, scalar)

TYPE(suaternion) :: q1, q_add_sk
REAL(sk) :: scalar

q_add_sk%w = q1%w + scalar
q_add_sk%x = q1%x + scalar
q_add_sk%y = q1%y + scalar
q_add_sk%z = q1%z + scalar

END FUNCTION q_add_sk

!------------------------------------------------------------------------------
! FUNCTION: q_mult_sk
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
FUNCTION q_mult_sk (q1, scalar)

TYPE(suaternion) :: q1, q_mult_sk
REAL(sk) :: scalar

q_mult_sk%w = q1%w * scalar
q_mult_sk%x = q1%x * scalar
q_mult_sk%y = q1%y * scalar
q_mult_sk%z = q1%z * scalar

END FUNCTION q_mult_sk

!------------------------------------------------------------------------------
! FUNCTION: q_dvd_sk
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
FUNCTION q_dvd_sk (q1, scalar)

TYPE(suaternion) :: q1, q_dvd_sk
REAL(sk) :: scalar

q_dvd_sk%w = q1%w / scalar
q_dvd_sk%x = q1%x / scalar
q_dvd_sk%y = q1%y / scalar
q_dvd_sk%z = q1%z / scalar

END FUNCTION q_dvd_sk

!------------------------------------------------------------------------------
! FUNCTION: quat_prod_sk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion division. Element wise, not to be confused with quat_prod_sk!
!
!> @param[in] q1 Quaternion
!> @param[in] q2 Quaternion
!> @return quat_prod_sk Result of the quaternion multiplication
!------------------------------------------------------------------------------ 
FUNCTION quat_prod_sk (q1, q2)

TYPE(suaternion) :: q1, q2, quat_prod_sk

quat_prod_sk%w = q1%w*q2%w-q1%x*q2%x-q1%y*q2%y-q1%z*q2%z
quat_prod_sk%x = q1%w*q2%x+q1%x*q2%w+q1%y*q2%z-q1%z*q2%y
quat_prod_sk%y = q1%w*q2%y-q1%x*q2%z+q1%y*q2%w+q1%z*q2%x
quat_prod_sk%z = q1%w*q2%z+q1%x*q2%y-q1%y*q2%x+q1%z*q2%w

END FUNCTION quat_prod_sk

!------------------------------------------------------------------------------
! FUNCTION: conjugate_quat_sk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Conjugate quaternion.
!
!> @param[in] quat Quaternion
!> @return conjugate_quat_sk Conjugate quaternion
!------------------------------------------------------------------------------ 
FUNCTION conjugate_quat_sk(quat)

TYPE(suaternion) :: quat, conjugate_quat_sk

conjugate_quat_sk%w =   quat%w
conjugate_quat_sk%x = - quat%x
conjugate_quat_sk%y = - quat%y
conjugate_quat_sk%z = - quat%z

END FUNCTION conjugate_quat_sk

!------------------------------------------------------------------------------
! FUNCTION: rotate_point_sk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Rotate a point with quaternions. Takes a unit quaternion as input
!
!> @param[in] quat Quaternion
!> @param[in] point Point to rotate with quaternion.
!> @return rotate_point_sk Rotated point.
!------------------------------------------------------------------------------
FUNCTION rotate_point_sk(quat, point)

TYPE(suaternion) :: quat
REAL(sk), DIMENSION(3) :: point, rotate_point_sk

TYPE(suaternion) :: pnt, rtt_pnt

pnt%w = 0._rk
pnt%x = point(1)
pnt%y = point(2)
pnt%z = point(3)

rtt_pnt = quat_prod_sk(quat_prod_sk(quat, pnt), conjugate_quat_sk(quat))

rotate_point_sk(1) = rtt_pnt%x
rotate_point_sk(2) = rtt_pnt%y
rotate_point_sk(3) = rtt_pnt%z

END FUNCTION rotate_point_sk


!------------------------------------------------------------------------------
! FUNCTION: quat_2_rot_mat_sk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Quaternion to rotation matrix.
!
!> @param[in] q Quaternion
!> @return quat_2_rot_mat_sk Rotation matrix
!------------------------------------------------------------------------------
FUNCTION quat_2_rot_mat_sk(q)

TYPE(suaternion) :: q
REAL(rk), DIMENSION(3,3) :: quat_2_rot_mat_sk

quat_2_rot_mat_sk(1,:) = [ q%w**2 + q%x**2 - q%y**2 - q%z**2  , &
    2._rk*q%x*q%y - 2._rk*q%w*q%z , 2._rk*q%x*q%z - 2._rk*q%w*q%y  ]
quat_2_rot_mat_sk(2,:) = [ 2._rk*q%x*q%y + 2._rk*q%w*q%z , &
    q%w**2 - q%x**2 + q%y**2 + - q%z**2 , 2._rk*q%y*q%z - 2._rk*q%w*q%x ]
quat_2_rot_mat_sk(3,:) = [ 2._rk*q%x*q%z - 2._rk*q%w*q%y , &
    2._rk*q%y*q%z + 2._rk*q%w*q%x , q%w**2 - q%x**2 - q%y**2 + q%z**2   ]

END FUNCTION quat_2_rot_mat_sk


!------------------------------------------------------------------------------
! FUNCTION: crpr_rk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Cross product of two vectors
!
!> @param[in] a Vector
!> @param[in] b Vector
!> @return crpr_rk Cross product
!------------------------------------------------------------------------------
FUNCTION crpr_rk(a, b)

REAL(rk), DIMENSION(3) :: crpr_rk,a,b

crpr_rk(1) = a(2) * b(3) - a(3) * b(2)
crpr_rk(2) = a(3) * b(1) - a(1) * b(3)
crpr_rk(3) = a(1) * b(2) - a(2) * b(1)

END FUNCTION crpr_rk

!------------------------------------------------------------------------------
! FUNCTION: crpr_sk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Cross product of two vectors
!
!> @param[in] a Vector
!> @param[in] b Vector
!> @return crpr_sk Cross product
!------------------------------------------------------------------------------
FUNCTION crpr_sk(a, b)

REAL(sk), DIMENSION(3) :: crpr_sk,a,b

crpr_sk(1) = a(2) * b(3) - a(3) * b(2)
crpr_sk(2) = a(3) * b(1) - a(1) * b(3)
crpr_sk(3) = a(1) * b(2) - a(2) * b(1)

END FUNCTION crpr_sk

!------------------------------------------------------------------------------
! FUNCTION: dist_pnt_pln_by_pnts_rk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Calculates the distance of a point to a plane, given by three points.
!> https://www.mathepower.com/ebenengleichungen.php :-)
!
!> @param[in] pp1 Point 1 of the plane
!> @param[in] pp2 Point 2 of the plane
!> @param[in] pp3 Point 3 of the plane
!> @param[in] p Point
!> @return[out] di Distance
!------------------------------------------------------------------------------  
FUNCTION dist_pnt_pln_by_pnts_rk(pp1, pp2, pp3, p) Result(di)

REAL(rk), dimension(3), Intent(IN) :: pp1, pp2, pp3, p
REAL(rk), dimension(3) :: nn, unitv 
REAL(rk):: di, magni

!------------------------------------------------------------------------------
! The normal to the plane
!------------------------------------------------------------------------------
nn = crpr_rk(pp1-pp3, pp1-pp2)

!------------------------------------------------------------------------------
! Magnitude of the normal to the plane
!------------------------------------------------------------------------------
magni = SQRT(DOT_PRODUCT(nn, nn))

!------------------------------------------------------------------------------
! Unit vector
!------------------------------------------------------------------------------
unitv = nn / magni

di = DOT_PRODUCT(unitv, [p(1)-pp1(1), p(2)-pp1(2), p(3)-pp1(3)])

END FUNCTION dist_pnt_pln_by_pnts_rk

!------------------------------------------------------------------------------
! FUNCTION: dist_pnt_pln_by_pnts_sk
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Calculates the distance of a point to a plane, given by three points.
!> https://www.mathepower.com/ebenengleichungen.php :-)
!
!> @param[in] pp1 Point 1 of the plane
!> @param[in] pp2 Point 2 of the plane
!> @param[in] pp3 Point 3 of the plane
!> @param[in] p Point
!> @return[out] di Distance
!------------------------------------------------------------------------------  
FUNCTION dist_pnt_pln_by_pnts_sk(pp1, pp2, pp3, p) Result(di)

REAL(sk), dimension(3), Intent(IN) :: pp1, pp2, pp3, p
REAL(sk), dimension(3) :: nn, unitv 
REAL(sk):: di, magni

!------------------------------------------------------------------------------
! The normal to the plane
!------------------------------------------------------------------------------
nn = crpr_sk(pp1-pp3, pp1-pp2)

!------------------------------------------------------------------------------
! Magnitude of the normal to the plane
!------------------------------------------------------------------------------
magni = SQRT(DOT_PRODUCT(nn, nn))

!------------------------------------------------------------------------------
! Unit vector
!------------------------------------------------------------------------------
unitv = nn / magni

di = DOT_PRODUCT(unitv, [p(1)-pp1(1), p(2)-pp1(2), p(3)-pp1(3)])

END FUNCTION dist_pnt_pln_by_pnts_sk


!------------------------------------------------------------------------------
! FUNCTION: diagonal3D
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Returns the diagonal of a box
!
!> @param[in] box Input box/cube
!> @return diagonal Diagonal distance
!------------------------------------------------------------------------------
FUNCTION diagonal3D(box) result (diagonal)

REAL(rk) :: box(3), diagonal

diagonal = SQRT(box(1)**2 + box(2)**2 + box(3)**2)

END FUNCTION diagonal3D

END MODULE math
