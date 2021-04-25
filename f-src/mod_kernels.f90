!---------------------------------------------------------------------------------------------------
! Module containing Kernels/Convolutional matrices for image processing
! mod_kernels.f90
!
! Author Benjamin Schnabel
! Date    23.04.2021

! SUBROUTINE convolution(dims, image, sizeKernel, kernel)
! SUBROUTINE kernel_identity(kernel, sizeKernel)
! SUBROUTINE kernel_gauss(kernel, sizeKernel, sigma)

MODULE kernels

USE ISO_FORTRAN_ENV
USE standards
! USE OMP_LIB

IMPLICIT NONE

CONTAINS

SUBROUTINE convolution(dimImage, image, sizeKernel, kernel)

   ! External variables
   INTEGER(KIND = ik)           , DIMENSION(3)                        , INTENT(IN)    :: dimImage
   INTEGER(KIND = ik)                                                 , INTENT(IN)    :: sizeKernel

   REAL   (KIND = rk)           , DIMENSION(dimImage(1), dimImage(2), dimImage(3)), INTENT(INOUT) :: image
   REAL   (KIND = rk)           , DIMENSION(sizeKernel, sizeKernel)   , INTENT(IN)    :: kernel

   ! Parameter
   INTEGER(KIND = ik), PARAMETER                                                      :: nthreads = 2

   ! Internal variables
   INTEGER(KIND = ik)                                                                 :: i, j, k, l, m
   INTEGER(KIND = ik)                                                                 :: border
   INTEGER(KIND = ik)                                                                 :: borderPO
   INTEGER(KIND = ik)           , DIMENSION(3)                                        :: dimImagePadded

   REAL   (KIND = rk)                                                                 :: accumulator
   REAL   (KIND = rk)           , DIMENSION(:,:,:)                    , ALLOCATABLE   :: imagePadded
   REAL   (KIND = rk)           , DIMENSION(dimImage(1), dimImage(2), dimImage(3))    :: imageCopy
   
   !-------------------------------
   ! Calculation
   !-------------------------------

!    call omp_set_num_threads(nthreads)
    
   border = FLOOR(REAL(sizeKernel) / 2)

   ALLOCATE(imagePadded(dimImage(1) + 2 * border, dimImage(2) + 2 * border, dimImage(3) + 2 * border))

   dimImagePadded = SHAPE(imagePadded)
   imagePadded = 0
   
   borderPO = border + 1

! TODO: Setup export to log with DBG_LVL == 2 !
!    write(*,*) border
!    write(*,*) borderPO
!    write(*,*) dimImage
!    write(*,*) SHAPE(image)
!    write(*,*) SHAPE(imagePadded)

   imagePadded(borderPO : dimImage(1) + border, borderPO : dimImage(2) + border, borderPO : dimImage(3) + border) = image

   DO i = borderPO, dimImage(3)
           DO j = borderPO, dimImage(2)
                   DO k = borderPO, dimImage(1)
                           accumulator = 0
                           DO l = - border, border
                                   !$OMP PARALLEL
                                   DO m = - border, border
                                           accumulator = accumulator + (kernel(l + border + 1, m + border + 1) &
                                                   * imagePadded(k + l, j + m, i))
                                   END DO
                                   !$OMP END PARALLEL
                           END DO
                           imageCopy(k - border, j - border, i) = accumulator
                   END DO
           END DO
   END DO

   image = imageCopy
END SUBROUTINE convolution

SUBROUTINE kernel_identity(kernel, sizeKernel)
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
END SUBROUTINE kernel_identity

SUBROUTINE kernel_gauss(kernel, sizeKernel, sigma)
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

   DO i = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
           DO j = -(sizeKernel - 1) / 2, (sizeKernel - 1) / 2
                   x0 = (sizeKernel + 1) / 2
                   y0 = (sizeKernel + 1) / 2
                   x = i + x0
                   y = j + y0

                   kernel(INT(x), INT(y)) = (1 / (2 * pi * sigma**2)) &
                           * EXP(- (((x - x0)**2 + (y - y0)**2) / (2 * sigma**2)))
           END DO
   END DO

   sumKernel = sum(kernel)
   kernel = kernel / sumKernel
END SUBROUTINE kernel_gauss

END MODULE kernels
