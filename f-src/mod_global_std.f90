!------------------------------------------------------------------------------
! MODULE: global standards
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
! @description: 
!> Module containing all recurring definitions of kinds and numbers.
!------------------------------------------------------------------------------
MODULE global_std

IMPLICIT NONE

! General constants
INTEGER, PARAMETER :: sik = 2    ! INTEGER Kind
INTEGER, PARAMETER :: ik  = 4    ! INTEGER Kind
INTEGER, PARAMETER :: rk  = 8    ! Real    Kind
INTEGER, PARAMETER :: mcl = 512  ! Maximal character  length
INTEGER, PARAMETER :: hcl = 256  ! Half    character  length
INTEGER, PARAMETER :: scl = 64   ! Short   character  length

!-- File handles, debug_lvl and suffix
INTEGER(KIND=ik), PARAMETER :: timer_level = 3 ! 1 ! 2
INTEGER(KIND=ik), PARAMETER :: dbg_lvl = 1
CHARACTER(LEN=mcl)          :: mssg = ''

CHARACTER(len=*), PARAMETER :: std_lnbrk = REPEAT("_", 100) ! "(100('-'))"

! Dynamically assigned, must fit to program
INTEGER(KIND=ik)            :: std_out
INTEGER(KIND=ik), PARAMETER :: std_in  = 5
INTEGER(KIND=ik), PARAMETER :: std_err = 0

! Unsusual, as you may want to see the error messages asap in your root directory.
INTEGER(KIND=ik), PARAMETER :: fhsterr  = 10

!-- Mpi-specific kinds
INTEGER, PARAMETER :: mpi_ik = 4 ! MPI INTEGER Kind; Compile with corresponding mpi!!

!-- StdOut Characters
CHARACTER(len=5) :: creturn = achar(13)

!-- Mpi-specific kinds
INTEGER,PARAMETER :: mik = 4       ! MPI INTEGER Kind; Compile with corresponding mpi!!

!-- Math standards
REAL(KIND=rk), PARAMETER :: is_zero    = 1.E-9_rk
REAL(KIND=rk), PARAMETER :: sq2        = sqrt(2._rk)
REAL(KIND=rk), PARAMETER :: pi         = 4.D0*DATAN(1.D0) !acos(-1._rk)
REAL(KIND=rk), PARAMETER :: inv180     = 1._rk/180._rk
REAL(KIND=rk), PARAMETER :: pi_div_180 = acos(-1._rk)/180._rk

!-- Higher dimensional numbers
TYPE Quaternion
REAL (KIND=rk)            :: w,x,y,z
END TYPE Quaternion

END MODULE global_std

!------------------------------------------------------------------------------
! MODULE: mechanical
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
! DESCRIPTION: 
!> Module containing all recurring definitions of kinds and vmbers.
!
! REVISION HISTORY:
! 26 09 2021 - Initial Version
!------------------------------------------------------------------------------
MODULE mechanical

USE global_std

IMPLICIT NONE

! Add other parameters if necessary.
TYPE materialcard
    REAL(KIND=rk) :: E
    REAL(KIND=rk) :: nu
    ! For use with effective nummerical stiffness calculations 
    REAL(KIND=rk), DIMENSION(3) :: phdsize ! Physical domain/ Macro element size 
END TYPE materialcard

CONTAINS

!------------------------------------------------------------------------------
! FUNCTION: lamee_lambda
!------------------------------------------------------------------------------ 
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to return the lamé constant lambda
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return lambda
!------------------------------------------------------------------------------  
FUNCTION lamee_lambda(E, v) RESULT (lambda)
 
   REAL (KIND=rk) :: E, v
   REAL (KIND=rk), DIMENSION(6,6) :: lambda

   lambda = E*v/((1._rk+v)*(1._rk-2._rk*v))

END FUNCTION lamee_lambda

!------------------------------------------------------------------------------
! FUNCTION: lamee_mu_shear
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to return the lamé constant µ / the shear modulus G
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return G
!------------------------------------------------------------------------------  
FUNCTION lamee_mu_shear(E, v) RESULT (G)

   REAL (KIND=rk) :: E, v
   REAL (KIND=rk), DIMENSION(6,6) :: G ! (shear modulus) 

   G = E / (2._rk*(1._rk+v))

   ! Just for info
   ! mu = G 

END FUNCTION lamee_mu_shear


!------------------------------------------------------------------------------
! FUNCTION: bulk_modulus
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to return the bulk modulus
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return k
!------------------------------------------------------------------------------  
FUNCTION bulk_modulus(E, v) RESULT (k)

   REAL (KIND=rk) :: E, v
   REAL (KIND=rk), DIMENSION(6,6) :: k ! (shear modulus) 

   k = E / (3._rk*(1._rk-(2._rk*v)))

END FUNCTION bulk_modulus


!------------------------------------------------------------------------------
! FUNCTION: iso_compliance_voigt
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to quickly generate an isotropic 2nd rank compliance tensor
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return t_iso_inv
!------------------------------------------------------------------------------  
FUNCTION iso_compliance_voigt(E, v) RESULT (t_iso_inv)

   REAL (KIND=rk) :: E, v
   REAL (KIND=rk), DIMENSION(6,6) :: t_iso_inv

   REAL (KIND=rk) :: fctr

   fctr = 1._rk/E

!  sum of symmetric components - therefore: eps_12+eps_21 => 2*eps_12 etc.
   t_iso_inv(:,1)=[ 1._rk,    -v,    -v,  .0_rk,           .0_rk,                  .0_rk    ]
   t_iso_inv(:,2)=[    -v, 1._rk,    -v,  .0_rk,           .0_rk,                  .0_rk    ]
   t_iso_inv(:,3)=[    -v,    -v, 1._rk,  .0_rk,           .0_rk,                  .0_rk    ]
   t_iso_inv(:,4)=[ .0_rk, .0_rk, .0_rk,  2._rk*(1._rk+v), .0_rk,                  .0_rk    ]
   t_iso_inv(:,5)=[ .0_rk, .0_rk, .0_rk,  .0_rk,           2._rk*(1._rk+v),        .0_rk    ]
   t_iso_inv(:,6)=[ .0_rk, .0_rk, .0_rk,  .0_rk,           .0_rk,           2._rk*(1._rk+v) ]

   t_iso_inv = t_iso_inv*fctr

END FUNCTION iso_compliance_voigt


!------------------------------------------------------------------------------
! FUNCTION: iso_compliance_kelvin
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to quickly generate an isotropic 2nd rank compliance tensor
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return t_iso_inv
!------------------------------------------------------------------------------  
FUNCTION iso_compliance_kelvin(E, v) RESULT (t_iso_inv)

   REAL (KIND=rk) :: E, v
   REAL (KIND=rk), DIMENSION(6,6) :: t_iso_inv

   REAL (KIND=rk) :: fctr

   fctr = 1._rk/E

   t_iso_inv(:,1)=(/ 1._rk,    -v,    -v,  .0_rk,   .0_rk,   .0_rk   /)
   t_iso_inv(:,2)=(/    -v, 1._rk,    -v,  .0_rk,   .0_rk,   .0_rk   /)
   t_iso_inv(:,3)=(/    -v,    -v, 1._rk,  .0_rk,   .0_rk,   .0_rk   /)
   t_iso_inv(:,4)=(/ .0_rk, .0_rk, .0_rk,  1._rk+v, .0_rk,   .0_rk   /)
   t_iso_inv(:,5)=(/ .0_rk, .0_rk, .0_rk,  .0_rk,   1._rk+v, .0_rk   /)
   t_iso_inv(:,6)=(/ .0_rk, .0_rk, .0_rk,  .0_rk,   .0_rk,   1._rk+v /)

   t_iso_inv = t_iso_inv*fctr

END FUNCTION iso_compliance_kelvin



!------------------------------------------------------------------------------
! FUNCTION: iso_stiffness_voigt
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to quickly generate an isotropic 2nd rank stiffness tensor
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return t_iso
!------------------------------------------------------------------------------  
FUNCTION iso_stiffness_voigt(E, v) RESULT (t_iso)

   REAL (KIND=rk) :: E, v
   REAL (KIND=rk), DIMENSION(6,6) :: t_iso

   REAL (KIND=rk) :: fctr, fctr_shear
   
   fctr       = E / ((1._rk+v)*(1._rk-(2._rk*v)))
   fctr_shear = 1._rk-(2._rk*v)
   t_iso(:,1)=[ 1._rk-v,  v,     v,       .0_rk, .0_rk, .0_rk ]
   t_iso(:,2)=[    v , 1._rk-v,  v,       .0_rk, .0_rk, .0_rk ]
   t_iso(:,3)=[    v ,     v, 1._rk-v,    .0_rk, .0_rk, .0_rk ]
   t_iso(:,4)=[ .0_rk, .0_rk, .0_rk, fctr_shear/2, .0_rk, .0_rk ]
   t_iso(:,5)=[ .0_rk, .0_rk, .0_rk, .0_rk, fctr_shear/2, .0_rk ]
   t_iso(:,6)=[ .0_rk, .0_rk, .0_rk, .0_rk, .0_rk, fctr_shear/2 ]

   t_iso = t_iso*fctr ! Elementwise operation
END FUNCTION iso_stiffness_voigt


!------------------------------------------------------------------------------
! FUNCTION: iso_stiffness_kelvin
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/vM
!
!> @brief
!> Function to quickly generate an isotropic 2nd rank stiffness tensor
!
!> @param[in] E Young moduluss     
!> @param[in] v      
!> @return t_iso
!------------------------------------------------------------------------------  
FUNCTION iso_stiffness_kelvin(E, v) RESULT (t_iso)

   REAL (KIND=rk) :: E, v
   REAL (KIND=rk), DIMENSION(6,6) :: t_iso

   REAL (KIND=rk) ::fctr, fctr_shear
   
   fctr       = E / ((1._rk+v)*(1._rk-(2._rk*v))) 
   fctr_shear = 1._rk-(2._rk*v)

   ! Factor of 2 not required as for normal components: sigma=E* eps
   ! Factor of 2 not required as for shear  components: sigma=E*2eps
   t_iso(:,1)=[ 1._rk-v,  v,     v,       .0_rk, .0_rk, .0_rk ]
   t_iso(:,2)=[    v , 1._rk-v,  v,       .0_rk, .0_rk, .0_rk ]
   t_iso(:,3)=[    v ,     v, 1._rk-v,    .0_rk, .0_rk, .0_rk ]
   t_iso(:,4)=[ .0_rk, .0_rk, .0_rk, fctr_shear, .0_rk, .0_rk ]
   t_iso(:,5)=[ .0_rk, .0_rk, .0_rk, .0_rk, fctr_shear, .0_rk ]
   t_iso(:,6)=[ .0_rk, .0_rk, .0_rk, .0_rk, .0_rk, fctr_shear ]

   t_iso = t_iso*fctr ! Elementwise operation

END FUNCTION iso_stiffness_kelvin

END MODULE mechanical