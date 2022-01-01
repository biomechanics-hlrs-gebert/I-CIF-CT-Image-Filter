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

USE ISO_FORTRAN_ENV
USE global_std
USE math

IMPLICIT NONE

INTERFACE filter
   MODULE PROCEDURE filter_ik2
   MODULE PROCEDURE filter_ik4
END INTERFACE filter

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: filter_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Filter routine accepting 3d kernels of the shape dim1=dim2=dim3 of data type
!> integer kind 2.
!
!> @param[in] subarray Input data
!> @param[in] kernel Filter kernel
!> @param[in] srb Subarray reduced boundaries. Subarray without overlap to other control volume/domain.
!> @param[out] result_subarray Output data
!------------------------------------------------------------------------------  
SUBROUTINE filter_ik2(subarray, kernel, srb, result_subarray)

    INTEGER(KIND=INT16), DIMENSION(:,:,:), INTENT(IN) :: subarray
    REAL(KIND=rk)      , DIMENSION(:,:,:), INTENT(IN) :: kernel
    INTEGER(KIND=ik)   , DIMENSION(6),     INTENT(IN) :: srb ! subarray_reduced_bndaries
    INTEGER(KIND=INT16), DIMENSION(:,:,:), INTENT(OUT) :: result_subarray

    INTEGER(KIND=ik) :: accumulator, ii, jj, kk, ll, mm, nn, border
    INTEGER(KIND=ik), DIMENSION(3) :: kernel_size

    kernel_size = SIZE(kernel, DIM=1)

    border = (kernel_size(1)-1)/2
    
    DO kk = srb(3), srb(6)
    DO jj = srb(2), srb(5)
    DO ii = srb(1), srb(4)
        accumulator = 0
        DO nn = -border, border
        DO mm = -border, border
        DO ll = -border, border
            accumulator = accumulator + (kernel( ll+border+1_ik, mm+border+1_ik, nn+border+1_ik) * &
                    subarray(ii + ll, jj + mm, kk + nn))
        END DO
        END DO
        END DO
        result_subarray(ii - border, jj - border, kk - border) = accumulator
    END DO
    END DO
    END DO

END SUBROUTINE filter_ik2


!------------------------------------------------------------------------------
! SUBROUTINE: filter_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Filter routine accepting 3d kernels of the shape dim1=dim2=dim3 of data type
!> integer kind 4.
!
!> @param[in] subarray Input data
!> @param[in] kernel Filter kernel
!> @param[in] srb Subarray reduced boundaries. Subarray without overlap to other control volume/domain.
!> @param[out] result_subarray Output data
!------------------------------------------------------------------------------  
SUBROUTINE filter_ik4(subarray, kernel, srb, result_subarray)

    INTEGER(KIND=INT32), DIMENSION(:,:,:), INTENT(IN) :: subarray
    REAL(KIND=rk)      , DIMENSION(:,:,:), INTENT(IN) :: kernel
    INTEGER(KIND=ik), DIMENSION(6),     INTENT(IN) :: srb ! subarray_reduced_bndaries
    INTEGER(KIND=INT32), DIMENSION(:,:,:), INTENT(OUT) :: result_subarray

    INTEGER(KIND=ik) :: accumulator, ii, jj, kk, ll, mm, nn, border
    INTEGER(KIND=ik), DIMENSION(3) :: kernel_size

    kernel_size = SIZE(kernel, DIM=1)

    border = (kernel_size(1)-1)/2

    DO kk = srb(3), srb(6)
    DO jj = srb(2), srb(5)
    DO ii = srb(1), srb(4)
        accumulator = 0
        DO nn = -border, border
        DO mm = -border, border
        DO ll = -border, border
            accumulator = accumulator + (kernel( ll+border+1_ik, mm+border+1_ik, nn+border+1_ik) * &
                    subarray(ii + ll, jj + mm, kk + nn))
        END DO
        END DO
        END DO
        result_subarray(ii - border, jj - border, kk - border) = accumulator
    END DO
    END DO
    END DO

END SUBROUTINE filter_ik4


!------------------------------------------------------------------------------
! SUBROUTINE: kernel_identity_2d
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Provides the 2-dimensional kernel to process the image. 
!> Regarded "obsolete" since 2D-kernels do not work out well on randomly
!> positioned 3D images. Kept for eventual future purposes.
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
!> Regarded "obsolete" since 2D-kernels do not work out well on randomly
!> positioned 3D images. Kept for eventual future purposes.
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_gauss_2d(kernel, sizeKernel, sigma)
   ! Return the Gaussian kernel
   ! https://en.wikipedia.org/wiki/Gaussian_filter

   INTEGER(KIND=ik), INTENT(IN)  :: sizeKernel
   REAL(KIND=rk), INTENT(IN)  :: sigma
   REAL(KIND=rk), DIMENSION(:,:), INTENT(OUT) :: kernel

   INTEGER(KIND=ik) :: i, j
   REAL(KIND=rk) :: x0, y0, x, y
   REAL(KIND=rk) :: sumKernel

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
   REAL(KIND=rk), INTENT(IN) :: sigma
   REAL(KIND=rk), DIMENSION(:,:,:), INTENT(OUT) :: kernel

   INTEGER(KIND=ik) :: i, j, k
   REAL(KIND=rk) :: x0, y0, z0, x, y, z
   REAL(KIND=rk) :: sumKernel

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
!> Regarded "obsolete" since 2D-kernels do not work out well on randomly
!> positioned 3D images. Kept for eventual future purposes.
!
!> @param[in] sizeKernel [Voxel]     
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_box_2d(kernel)
   ! Return the Box filter kernel
   ! 6.869.csail.mit.edu/fa16/lecture/lecture3linearfilters.pdf

   ! Externel variables
   REAL(KIND=rk), DIMENSION(:,:,:), INTENT(OUT) :: kernel

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
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_box_3d(kernel)
   ! Return the Box filter kernel
   ! 6.869.csail.mit.edu/fa16/lecture/lecture3linearfilters.pdf

   REAL(KIND=rk), DIMENSION(:,:,:), INTENT(OUT) :: kernel

   kernel = (1.0_rk / SIZE(kernel))

END SUBROUTINE kernel_box_3d


!------------------------------------------------------------------------------
! FUNCTION: log_2d
!------------------------------------------------------------------------------  
!> @author Benjamin Schnabel - HLRS - NUM - schnabel@hlrs.de
!
!> @brief
!> Auxiliary function for kernel_log_2d
!> Regarded "obsolete" since 2D-kernels do not work out well on randomly
!> positioned 3D images. Kept for eventual future purposes.
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
    REAL(KIND=rk), INTENT(IN) :: sigma
    REAL(KIND=rk), DIMENSION(:,:), INTENT(OUT) :: kernel
    
    INTEGER(KIND=ik) :: i, j
    REAL(KIND=rk) :: x0, y0, x, y

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
    REAL(KIND=rk), INTENT(IN) :: sigma
    REAL(KIND=rk), DIMENSION(:,:,:), INTENT(OUT) :: kernel
    
    INTEGER(KIND=ik) :: i, j, k
    REAL(KIND=rk) :: x0, y0, z0, x, y, z
    
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
