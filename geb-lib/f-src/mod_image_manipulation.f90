!------------------------------------------------------------------------------
! MODULE: image_manipulation
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @description: 
!> Collection of image manipulating routines.
!------------------------------------------------------------------------------
MODULE image_manipulation

USE ISO_FORTRAN_ENV
USE global_std

IMPLICIT NONE

INTERFACE downscale
    MODULE PROCEDURE downscale_ik2
    MODULE PROCEDURE downscale_ik4
END INTERFACE downscale

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: downscale_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Reduce the resolution of a given image by a given integer factor.
!
!> @param[in] in_array Input array
!> @param[in] scale_factor Factor of high resolution voxels to get a mean ave.
!> @param[out] out_array Output array
!------------------------------------------------------------------------------  
SUBROUTINE downscale_ik2(in_array, scale_factor, out_array)

    INTEGER(INT16), DIMENSION(:,:,:), INTENT(OUT) :: out_array
    INTEGER(INT16), DIMENSION(:,:,:), INTENT(IN) :: in_array
    INTEGER(ik), DIMENSION(3), INTENT(IN) :: scale_factor

    INTEGER(ik), DIMENSION(3) :: shp_in
    INTEGER(ik) :: ii, jj, kk, ll, mm, nn
    REAL(rk) :: scale_volume

    scale_volume = REAL(PRODUCT(scale_factor), rk)
    shp_in = SHAPE(in_array)

    nn = 1_ik
    DO kk=1, SIZE(in_array, DIM=3)-scale_factor(3), scale_factor(3)
        IF (shp_in(3) < kk+scale_factor(3)) CYCLE

        mm = 1_ik
        DO jj=1, SIZE(in_array, DIM=2)-scale_factor(2), scale_factor(2)
            IF (shp_in(2) < jj+scale_factor(2)) CYCLE

            ll = 1_ik
            DO ii=1, SIZE(in_array, DIM=1)-scale_factor(1), scale_factor(1)
                IF (shp_in(1) < ii+scale_factor(1)) CYCLE

                !------------------------------------------------------------------------------  
                ! DO NOT interchange SUM and REAL as the sum will overflow without conversion.
                !------------------------------------------------------------------------------  
                out_array(ll, mm, nn) = INT(SUM(REAL(in_array(&
                    ii:ii+scale_factor(1)-1_ik, &
                    jj:jj+scale_factor(2)-1_ik, &
                    kk:kk+scale_factor(3)-1_ik), rk)) / scale_volume, INT16)
                ll = ll + 1_ik
            END DO

            mm = mm + 1_ik
        END DO

        nn = nn + 1_ik
    END DO
END SUBROUTINE downscale_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: downscale_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Reduce the resolution of a given image by a given integer factor.
!
!> @param[in] in_array Input array
!> @param[in] scale_factor Factor of high resolution voxels to get a mean ave.
!> @param[out] out_array Output array
!------------------------------------------------------------------------------  
SUBROUTINE downscale_ik4(in_array, scale_factor, out_array)

    INTEGER(INT32), DIMENSION(:,:,:), INTENT(OUT) :: out_array
    INTEGER(INT32), DIMENSION(:,:,:), INTENT(IN) :: in_array
    INTEGER(ik), DIMENSION(3), INTENT(IN) :: scale_factor

    INTEGER(ik), DIMENSION(3) :: shp_in
    INTEGER(ik) :: ii, jj, kk, ll, mm, nn
    REAL(rk) :: scale_volume

    scale_volume = REAL(PRODUCT(scale_factor), rk)
    shp_in = SHAPE(in_array)

    nn = 1_ik
    DO kk=1, SIZE(in_array, DIM=3)-scale_factor(3), scale_factor(3)
        IF (shp_in(3) < kk+scale_factor(3)) CYCLE

        mm = 1_ik
        DO jj=1, SIZE(in_array, DIM=2)-scale_factor(2), scale_factor(2)
            IF (shp_in(2) < jj+scale_factor(2)) CYCLE

            ll = 1_ik
            DO ii=1, SIZE(in_array, DIM=1)-scale_factor(1), scale_factor(1)
                IF (shp_in(1) < ii+scale_factor(1)) CYCLE

                !------------------------------------------------------------------------------  
                ! DO NOT interchange SUM and REAL as the sum will overflow without conversion.
                !------------------------------------------------------------------------------  
                out_array(ll, mm, nn) = INT(SUM(REAL(in_array(&
                    ii:ii+scale_factor(1)-1_ik, &
                    jj:jj+scale_factor(2)-1_ik, &
                    kk:kk+scale_factor(3)-1_ik), rk)) / scale_volume, INT32)
                ll = ll + 1_ik
            END DO

            mm = mm + 1_ik
        END DO

        nn = nn + 1_ik
    END DO
END SUBROUTINE downscale_ik4

END MODULE image_manipulation

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

INTEGER(INT16), DIMENSION(:,:,:), INTENT(INOUT) :: subarray, result_subarray

INCLUDE "./include_f90/filter.f90"

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

INTEGER(INT32), DIMENSION(:,:,:), INTENT(INOUT) :: subarray, result_subarray

INCLUDE "./include_f90/filter.f90"

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

    INTEGER(ik), INTENT(IN)  :: sizeKernel
    REAL(rk), DIMENSION(sizeKernel, sizeKernel), INTENT(OUT) :: kernel
    
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

INTEGER(ik), INTENT(IN)  :: sizeKernel
REAL(rk), INTENT(IN)  :: sigma
REAL(rk), DIMENSION(:,:), INTENT(OUT) :: kernel

INTEGER(ik) :: i, j
REAL(rk) :: x0, y0, x, y
REAL(rk) :: sumKernel

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

INTEGER(ik), INTENT(IN) :: sizeKernel
REAL(rk), INTENT(IN) :: sigma
REAL(rk), DIMENSION(:,:,:), INTENT(OUT) :: kernel

INTEGER(ik) :: i, j, k
REAL(rk) :: x0, y0, z0, x, y, z
REAL(rk) :: sumKernel

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
! SUBROUTINE: kernel_gauss_3d_fix_5_15
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Provides the 3-dimensional Gauss kernel of size 5 and sigma 1.5
!
!> @param[out] kernel Values to multiply with 
!------------------------------------------------------------------------------  
SUBROUTINE kernel_gauss_3d_fix_5_15(kernel)

REAL(rk), DIMENSION(:,:,:), INTENT(OUT) :: kernel

kernel(:,1,1)=[ 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213 ]
kernel(:,2,1)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,3,1)=[ 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872 ]
kernel(:,4,1)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,5,1)=[ 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213 ]

kernel(:,1,2)=[ 0.120078384243213, 0.120078384243214, 0.120078384243213, 0.120078384243214, 0.120078384243213 ]
kernel(:,2,2)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,3,2)=[ 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872 ]
kernel(:,4,2)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,5,2)=[ 0.120078384243213, 0.120078384243214, 0.120078384243213, 0.120078384243214, 0.120078384243213 ]

kernel(:,1,3)=[ 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213 ]
kernel(:,2,3)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,3,3)=[ 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872 ]
kernel(:,4,3)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,5,3)=[ 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213 ]

kernel(:,1,4)=[ 0.120078384243213, 0.120078384243214, 0.120078384243213, 0.120078384243214, 0.120078384243213 ]
kernel(:,2,4)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,3,4)=[ 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872 ]
kernel(:,4,4)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,5,4)=[ 0.120078384243213, 0.120078384243214, 0.120078384243213, 0.120078384243214, 0.120078384243213 ]

kernel(:,1,5)=[ 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213 ]
kernel(:,2,5)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,3,5)=[ 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872, 0.292081718342872 ]
kernel(:,4,5)=[ 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350, 0.233880756585350 ]
kernel(:,5,5)=[ 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213, 0.120078384243213 ]

END SUBROUTINE kernel_gauss_3d_fix_5_15

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
REAL(rk), DIMENSION(:,:,:), INTENT(OUT) :: kernel

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

REAL(rk), DIMENSION(:,:,:), INTENT(OUT) :: kernel

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
    
INTEGER(ik), INTENT(IN)  :: x, y
REAL(rk), INTENT(IN)  :: sigma

REAL(rk):: log_2d

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

INTEGER(ik), INTENT(IN) :: sizeKernel
REAL(rk), INTENT(IN) :: sigma
REAL(rk), DIMENSION(:,:), INTENT(OUT) :: kernel

INTEGER(ik) :: i, j
REAL(rk) :: x0, y0, x, y

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

INTEGER(ik), INTENT(IN) :: sizeKernel
REAL(rk), INTENT(IN) :: sigma
REAL(rk), DIMENSION(:,:,:), INTENT(OUT) :: kernel

INTEGER(ik) :: i, j, k
REAL(rk) :: x0, y0, z0, x, y, z

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
