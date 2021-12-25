!------------------------------------------------------------------------------
! MODULE: kernels
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
! DESCRIPTION: 
!> Module containing kernels/convolutional matrices for image processing.
!------------------------------------------------------------------------------
MODULE kernels

USE global_std

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: kernel_identity_2d
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Provides the 2-dimensional kernel to process the image
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_identity_2d(kernel, sizeKernel)

    INTEGER(KIND=ik), INTENT(IN)  :: sizeKernel
    REAL(KIND=rk), DIMENSION(sizeKernel, sizeKernel), INTENT(OUT) :: kernel
    
    INTEGER(KIND = ik) :: mp
    
    mp = 1 + (sizeKernel - 1) / 2
    
    kernel = 0
    kernel(mp, mp) = 1.0
END SUBROUTINE kernel_identity_2d

!------------------------------------------------------------------------------
! SUBROUTINE: kernel_identity_3d
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Provides the 3-dimensional kernel to process the image
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_identity_3d(kernel, sizeKernel)
   INTEGER(KIND = ik), INTENT(IN)  :: sizeKernel
   REAL(KIND = rk), DIMENSION(sizeKernel, sizeKernel, sizeKernel), INTENT(OUT) :: kernel

   INTEGER(KIND = ik) :: mp

   mp = 1 + (sizeKernel - 1) / 2

   kernel = 0
   kernel(mp, mp, mp) = 1.0
END SUBROUTINE kernel_identity_3d


!------------------------------------------------------------------------------
! SUBROUTINE: kernel_gauss_2d
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Provides the 2-dimensional Gauss kernel to process the image
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_gauss_2d(kernel, sizeKernel, sigma)
   ! Return the Gaussian kernel
   ! https://en.wikipedia.org/wiki/Gaussian_filter

   INTEGER(KIND=ik), INTENT(IN)  :: sizeKernel
   REAL   (KIND=rk), INTENT(IN)  :: sigma
   REAL   (KIND=rk), DIMENSION(sizeKernel, sizeKernel), INTENT(OUT) :: kernel

   INTEGER(KIND=ik) :: i, j
   REAL   (KIND=rk) :: x0, y0, x, y
   REAL   (KIND=rk) :: sumKernel

   kernel = 0

   x0 = (sizeKernel + 1) / 2
   y0 = (sizeKernel + 1) / 2

   DO i = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
           DO j = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
                   x = i + x0
                   y = j + y0

                   kernel(INT(x), INT(y)) = (1 / (2 * pi * sigma**2)) &
                           * EXP(- (((x - x0)**2 + (y - y0)**2) / (2 * sigma**2)))
           END DO
   END DO

   sumKernel = SUM(kernel)
   kernel = kernel / sumKernel
END SUBROUTINE kernel_gauss_2d

!------------------------------------------------------------------------------
! SUBROUTINE: kernel_gauss_3d
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Provides the 3-dimensional Gauss kernel to process the image
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_gauss_3d(kernel, sizeKernel, sigma)

   INTEGER(KIND=ik), INTENT(IN) :: sizeKernel
   REAL   (KIND=rk), INTENT(IN) :: sigma
   REAL   (KIND=rk), DIMENSION(sizeKernel, sizeKernel, sizeKernel), INTENT(OUT) :: kernel

   INTEGER(KIND=ik) :: i, j, k
   REAL   (KIND=rk) :: x0, y0, z0, x, y, z
   REAL   (KIND=rk) :: sumKernel

   kernel = 0

   x0 = (sizeKernel + 1) / 2
   y0 = (sizeKernel + 1) / 2
   z0 = (sizeKernel + 1) / 2

   DO i = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
           DO j = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
                DO k = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
                        x = i + x0
                        y = j + y0
                        z = k + z0

                        kernel(INT(x), INT(y), INT(z)) = (1 / (2 * pi * sigma**2)) &
                                * EXP(- (((x - x0)**2 + (y - y0)**2 + (z - z0)**2) / (2 * sigma**2)));
                END DO
           END DO
   END DO

   sumKernel = SUM(kernel)
   kernel = kernel / sumKernel
END SUBROUTINE kernel_gauss_3d


!------------------------------------------------------------------------------
! SUBROUTINE: kernel_box_2d
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Provides the 2-dimensional Box kernel to process the image
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_box_2d(kernel, sizeKernel)
   ! Return the Box filter kernel
   ! 6.869.csail.mit.edu/fa16/lecture/lecture3linearfilters.pdf

   ! Externel variables
   INTEGER(KIND=ik), INTENT(IN) :: sizeKernel
   REAL   (KIND=rk), DIMENSION(sizeKernel, sizeKernel, sizeKernel), INTENT(OUT) :: kernel

   kernel = (1.0_rk / SIZE(kernel))

END SUBROUTINE kernel_box_2d

!------------------------------------------------------------------------------
! SUBROUTINE: kernel_box_3d
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Provides the 3-dimensional Box kernel to process the image
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_box_3d(kernel, sizeKernel)
   ! Return the Box filter kernel
   ! 6.869.csail.mit.edu/fa16/lecture/lecture3linearfilters.pdf

   INTEGER(KIND=ik), INTENT(IN) :: sizeKernel
   REAL   (KIND=rk), DIMENSION(sizeKernel, sizeKernel, sizeKernel), INTENT(OUT) :: kernel

   kernel = (1.0_rk / SIZE(kernel))

END SUBROUTINE kernel_box_3d


!------------------------------------------------------------------------------
! FUNCTION: log_2d
!------------------------------------------------------------------------------  
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Auxiliary function for kernel_log_2d
!
!> @param[in] sigma Gaussian
!> @param[in] x Value
!> @param[in] y Value
!------------------------------------------------------------------------------
FUNCTION log_2d(sigma, x, y)
     
    INTEGER(KIND=ik), INTENT(IN)  :: x, y
    REAL(KIND=rk), INTENT(IN)  :: sigma

    REAL(KIND=rk):: log_2d

    log_2d = - (1 / (pi * sigma**4)) * &
    (1 - ((x**2 + y**2) / (2 * sigma**2))) * &
    EXP(- ((x**2 + y**2) / (2 * sigma**2)))

END FUNCTION log_2d

!------------------------------------------------------------------------------
! SUBROUTINE: kernel_log_2d
!------------------------------------------------------------------------------  
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Provides the 2-dimensional laplacian of the gaussian kernel
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_log_2d(kernel, sizeKernel, sigma)
    ! Return the Laplacian of Gaussian kernel
    ! https://academic.mu.edu/phys/matthysd/web226/Lab02.htm
    
    INTEGER(KIND=ik), INTENT(IN) :: sizeKernel
    REAL   (KIND=rk), INTENT(IN) :: sigma
    REAL   (KIND=rk), DIMENSION(sizeKernel, sizeKernel), INTENT(OUT) :: kernel
    
    INTEGER(KIND=ik) :: i, j
    REAL   (KIND=rk) :: x0, y0, x, y

    kernel = 0

    x0 = (sizeKernel + 1) / 2
    y0 = (sizeKernel + 1) / 2
    
    DO i = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
        DO j = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
            x = i + x0
            y = j + y0
            kernel(INT(x), INT(y)) = log_2d(sigma, INT(x - x0, KIND = ik), INT(y - y0, KIND = ik))
        END DO
    END DO

    IF ((sigma == 1.4_rk) .AND. (sizeKernel == 9_ik)) THEN
        ! https://academic.mu.edu/phys/matthysd/web226/Lab02.htm
        kernel = ANINT(kernel * (-40 / log_2d(sigma, 0_ik, 0_ik)))
    END IF
    
    IF((SUM(kernel) > 0.1) .OR. (SUM(kernel) < -0.1) .OR. (sigma <= 0.5) .OR. ((sigma / SUM(kernel)) > 0.09)) THEN
        WRITE(*,*) 'Kernel values not good! Change sigma or kernel size'
    END IF

    !sumKernel = SUM(kernel)
    
END SUBROUTINE kernel_log_2d


!------------------------------------------------------------------------------
! FUNCTION: log_3d
!------------------------------------------------------------------------------  
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Auxiliary function for kernel_log_3d
!
!> @param[in] sigma Gaussian
!> @param[in] x Value
!> @param[in] y Value
!> @param[in] z Value
!------------------------------------------------------------------------------
FUNCTION log_3d(sigma, x, y, z)
    ! Return the Gaussian kernel
    ! https://en.wikipedia.org/wiki/Gaussian_filter

    IMPLICIT NONE
    
    ! Externel variables
    
    INTEGER(KIND = ik)                                                         , INTENT(IN)  :: x, y, z
    REAL(KIND = rk)                                                            , INTENT(IN)  :: sigma

    ! Internel variables

    REAL(KIND = rk)                                                                          :: log_3d
    
    !-------------------------------
    ! Calculation
    !-------------------------------

    log_3d = - (1 / (pi * sigma**4)) * (1 - ((x**2 + y**2 + z**2) / (2 * sigma**2))) * &
    EXP(- ((x**2 + y**2 + z**2) / (2 * sigma**2)))

END FUNCTION log_3d

!------------------------------------------------------------------------------
! SUBROUTINE: kernel_log_3d
!------------------------------------------------------------------------------  
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Provides the 3-dimensional laplacian of the gaussian kernel
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_log_3d(kernel, sizeKernel, sigma)
    ! Return the Laplacian of Gaussian kernel
    ! https://academic.mu.edu/phys/matthysd/web226/Lab02.htm
    
    INTEGER(KIND=ik), INTENT(IN) :: sizeKernel
    REAL   (KIND=rk), INTENT(IN) :: sigma
    REAL   (KIND=rk), DIMENSION(sizeKernel, sizeKernel, sizeKernel), INTENT(OUT) :: kernel
    
    INTEGER(KIND=ik) :: i, j, k
    REAL   (KIND=rk) :: x0, y0, z0, x, y, z
    
    kernel = 0

    x0 = (sizeKernel + 1) / 2
    y0 = (sizeKernel + 1) / 2
    z0 = (sizeKernel + 1) / 2
    
    DO i = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
        DO j = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
            DO k = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
                x = i + x0
                y = j + y0
                z = k + z0
                kernel(INT(x), INT(y), INT(z)) = log_3d(sigma, INT(x - x0, KIND = ik), &
                INT(y - y0, KIND = ik), INT(z - z0, KIND = ik))
            END DO
        END DO
    END DO

    IF ((sigma == 1.4_rk) .AND. (sizeKernel == 9_ik)) THEN
        ! https://academic.mu.edu/phys/matthysd/web226/Lab02.htm
        kernel = ANINT(kernel * (-40 / log_3d(sigma, 0_ik, 0_ik, 0_ik)))
    END IF
    
    IF((SUM(kernel) > 0.1) .OR. (SUM(kernel) < -0.1) .OR. (sigma <= 0.5) .OR. ((sigma / SUM(kernel)) > 0.09)) THEN
        WRITE(*,*) 'Kernel values not good! Change sigma or kernel size'
    END IF
END SUBROUTINE kernel_log_3d

END MODULE kernels
