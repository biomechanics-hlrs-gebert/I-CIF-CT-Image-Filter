!------------------------------------------------------------------------------
! MODULE: mechanical
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! DESCRIPTION: 
!> Module containing all recurring definitions of kinds and vmbers.
!------------------------------------------------------------------------------
MODULE mechanical

USE global_std

IMPLICIT NONE

!----------------------------------------------------------TD00-0_tc_Dev_tensor_data.meta--------------------
! Describe a tensor and its state in respect to the position of the control
! volume. Required to fully trace the origin of a stiffness matrix (tensor).
! 
! It is expected that all other information like domain size, filter options 
! etc. are described by the meta file format!
!------------------------------------------------------------------------------
TYPE tensor_2nd_rank_R66
   INTEGER(KIND=ik) :: dmn     ! Number of the control volume
   REAL(KIND=rk) :: density    ! Percentage of monolothic young modulus
   REAL(KIND=rk) :: doa_zener  ! Degree of anisotropy
   REAL(KIND=rk) :: doa_gebert ! Degree of anisotropy
   REAL(KIND=rk) :: sym        ! Symmetry deviation (quotient)
   REAL(KIND=rk), DIMENSION(3)   :: pos = 0._rk ! Positioon (deg) of alpha, eta, phi
   REAL(KIND=rk), DIMENSION(6,6) :: mat = 0._rk
END TYPE tensor_2nd_rank_R66

!------------------------------------------------------------------------------
! Characterize a material
!------------------------------------------------------------------------------
TYPE materialcard
    REAL(KIND=rk) :: E
    REAL(KIND=rk) :: nu
    ! For use with effective nummerical stiffness calculations 
    REAL(KIND=rk), DIMENSION(3) :: phdsize ! Physical domain/ Macro element size 
    REAL(KIND=rk), DIMENSION(3) :: delta   ! Spacing of a given discretized material
END TYPE materialcard

CONTAINS

!------------------------------------------------------------------------------
! FUNCTION: doa_zener
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Function to calculate the zener ratio, an anisotropy measure for 
!> orthotropic cases. https://en.wikipedia.org/wiki/Zener_ratio
!
!> @param[in] mat Young moduluss     
!> @return doa Degree of Zener-anisotropy
!------------------------------------------------------------------------------  
FUNCTION doa_zener(mat) RESULT (doa)

   REAL(KIND=rk), DIMENSION(6,6), INTENT(IN) :: mat
   REAL(KIND=rk) :: doa

   doa = 2*mat(4,4)/(mat(1,1)-mat(1,2))

END FUNCTION doa_zener


!------------------------------------------------------------------------------
! FUNCTION: doa_gebert
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Function to calculate the gebert ratio, an anisotropy measure for 
!> anisotropic cases.
!> FIRST DRAFT
!
!> @param[in] mat Young moduluss     
!> @return doa Degree of Gebert-anisotropy
!------------------------------------------------------------------------------  
FUNCTION doa_gebert(mat) RESULT (doa)

   REAL(KIND=rk), DIMENSION(6,6), INTENT(IN) :: mat
   REAL(KIND=rk) :: doa

   doa = ((mat(4,4)+mat(5,5)+mat(6,6))/3._rk + &
         (SUM(mat(1:3, 4:6)) + SUM(mat(4, 5:6)) + mat(5,6) / 12._rk)) &
            / (mat(1,1)-mat(1,2)) 

END FUNCTION doa_gebert

!------------------------------------------------------------------------------
! FUNCTION: lamee_lambda
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
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


!------------------------------------------------------------------------------
! FUNCTION: gebert_density_voigt
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Calculates the density of the control volume by its mechanical performance
!> in respect to its monolithical stiffness.
!
!> @param[in] mat Input tensor/matrix     
!> @param[in] E Young modulus of the monolithic material
!> @param[in] v Poissions ratio of the monolithic material
!> @return density Returns the gebert-density
!------------------------------------------------------------------------------  
FUNCTION gebert_density_voigt(mat, E, v) RESULT (density)

   REAL (KIND=rk), DIMENSION(6,6) :: mat
   REAL (KIND=rk) :: E, v, density

   REAL (KIND=rk), DIMENSION(6,6) :: dmat, voigt_mat

   voigt_mat = iso_stiffness_voigt(E, v)

   !------------------------------------------------------------------------------  
   ! Minor diagonals/zero entries of the density matrix dmat are inf.  
   !------------------------------------------------------------------------------    
   dmat = mat / voigt_mat 

   density =((dmat(1,1)+dmat(2,2)+dmat(3,3))/3._rk + &     
             (dmat(4,4)+dmat(5,5)+dmat(6,6))/3._rk + &     
             (dmat(1,2)+dmat(1,3)+dmat(2,3))/3._rk) / 3._rk

   ! density = SUM(mat)/36._rk / E

END FUNCTION gebert_density_voigt

END MODULE mechanical