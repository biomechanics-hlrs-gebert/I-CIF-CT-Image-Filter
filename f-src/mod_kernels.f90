!---------------------------------------------------------------------------------------------------
! Module containing Kernels/Convolutional matrices for image processing
! mod_kernels.f90
!
! Authors Benjamin Schnabel     HLRS, NUM       schnabel@hlrs.de
!         Johannes Gebert       HLRS, NUM         gebert@hlrs.de
! Date    23.04.2021

! SUBROUTINE kernel_identity_2d(kernel, sizeKernel)
! SUBROUTINE kernel_identity_3d(kernel, sizeKernel)
! SUBROUTINE kernel_gauss_2d(kernel, sizeKernel, sigma)
! SUBROUTINE kernel_gauss_3d(kernel, sizeKernel, sigma)

MODULE kernels

USE ISO_FORTRAN_ENV
USE standards

IMPLICIT NONE

CONTAINS

SUBROUTINE kernel_identity_2d(kernel, sizeKernel)
   ! Return the identity kernel
   ! https://en.wikipedia.org/wiki/Kernel_(image_processing)

   ! Parameter

   ! External variables
        INTEGER(KIND = ik)                                              , INTENT(IN)  :: sizeKernel
        REAL   (KIND = rk)           , DIMENSION(sizeKernel, sizeKernel), INTENT(OUT) :: kernel
     
        ! Internal variables
        INTEGER(KIND = ik) :: mp
     
        !-------------------------------
        ! Calculation
        !-------------------------------
        
        mp = 1 + (sizeKernel - 1) / 2
     
        kernel = 0
        kernel(mp, mp) = 1.0
END SUBROUTINE kernel_identity_2d

SUBROUTINE kernel_identity_3d(kernel, sizeKernel)
   ! Return the identity kernel
   ! https://en.wikipedia.org/wiki/Kernel_(image_processing)

   ! Parameter

   ! External variables
   INTEGER(KIND = ik)                                                          , INTENT(IN)  :: sizeKernel
   REAL   (KIND = rk)           , DIMENSION(sizeKernel, sizeKernel, sizeKernel), INTENT(OUT) :: kernel

   ! Internal variables
   INTEGER(KIND = ik) :: mp

   !-------------------------------
   ! Calculation
   !-------------------------------
   
   mp = 1 + (sizeKernel - 1) / 2

   kernel = 0
   kernel(mp, mp, mp) = 1.0
END SUBROUTINE kernel_identity_3d

SUBROUTINE kernel_gauss_2d(kernel, sizeKernel, sigma)
   ! Return the Gaussian kernel
   ! https://en.wikipedia.org/wiki/Gaussian_filter

   ! Externel variables
   INTEGER(KIND = ik)                                              , INTENT(IN)  :: sizeKernel
   REAL   (KIND = rk)                                              , INTENT(IN)  :: sigma
   REAL   (KIND = rk)           , DIMENSION(sizeKernel, sizeKernel), INTENT(OUT) :: kernel

   ! Internel variables
   INTEGER(KIND = ik)                                                            :: i, j
   REAL   (KIND = rk)                                                            :: x0, y0, x, y
   REAL   (KIND = rk)                                                            :: sumKernel

   !-------------------------------
   ! Calculation
   !-------------------------------

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

SUBROUTINE kernel_gauss_3d(kernel, sizeKernel, sigma)
   ! Return the Gaussian kernel
   ! https://en.wikipedia.org/wiki/Gaussian_filter

   ! Externel variables
   INTEGER(KIND = ik)                                                          , INTENT(IN)   :: sizeKernel
   REAL   (KIND = rk)                                                          , INTENT(IN)   :: sigma
   REAL   (KIND = rk)           , DIMENSION(sizeKernel, sizeKernel, sizeKernel), INTENT(OUT)  :: kernel

   ! Internel variables
   INTEGER(KIND = ik)                                                                         :: i, j, k
   REAL   (KIND = rk)                                                                         :: x0, y0, z0, x, y, z
   REAL   (KIND = rk)                                                                         :: sumKernel

   !-------------------------------
   ! Calculation
   !-------------------------------

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

END MODULE kernels
