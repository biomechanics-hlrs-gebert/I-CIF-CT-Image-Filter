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

SUBROUTINE convolution_un_padded_input(dims, image, sizeKernel, kernel, padding_switch, result_image)
        ! padding_switch .EQV. .TRUE.  - image is     padded
        ! padding_switch .EQV. .FALSE. - image is not padded

        ! External variables
        INTEGER(KIND = ik)           , DIMENSION(3)                        , INTENT(IN)    :: dims
        INTEGER(KIND = ik)                                                 , INTENT(IN)    :: sizeKernel
        REAL   (KIND = rk)           , DIMENSION(sizeKernel, sizeKernel)   , INTENT(IN)    :: kernel
        LOGICAL                                                            , INTENT(IN)    :: padding_switch
        INTEGER(KIND = ik)           , DIMENSION(dims(1), dims(2), dims(3)), INTENT(IN)    :: image
        INTEGER(KIND = ik)           , DIMENSION(dims(1), dims(2), dims(3)), INTENT(OUT)   :: result_image

        ! Parameter
        INTEGER(KIND = ik), PARAMETER                                                      :: nthreads = 2

        ! Internal variables
        INTEGER(KIND = ik)                                                                 :: ii, jj, kk, ll, mm
        INTEGER(KIND = ik)                                                                 :: border
        INTEGER(KIND = ik)                                                                 :: borderPO
        INTEGER(KIND = ik)           , DIMENSION(3)                                        :: dimsPadded

        REAL   (KIND = rk)                                                                 :: accumulator
        REAL   (KIND = rk)           , DIMENSION(:,:,:)                    , ALLOCATABLE   :: imagePadded

        !-------------------------------
        ! Calculation
        !-------------------------------

        !    call omp_set_num_threads(nthreads)

        border = FLOOR(REAL(sizeKernel) / 2)
        borderPO = border + 1

        IF (padding_switch .EQV. .FALSE.) THEN
                ! Image inserted here is not padded (with anything)
                ALLOCATE(imagePadded(dims(1) + 2 * border, dims(2) + 2 * border, dims(3) + 2 * border))

                dimsPadded = SHAPE(imagePadded)
                imagePadded(borderPO : dims(1) + border, borderPO : dims(2) + border, borderPO : dims(3) + border) = image
        ELSE
                ! Image was send with padding. However, it must fit to the Kernel provided. Otherwise, the result will be corrupted
                dimsPadded = dims
                imagePadded = image
        END IF

! TODO: Setup export to log with DBG_LVL == 2 !
!    write(*,*) border
!    write(*,*) borderPO
!    write(*,*) dims
!    write(*,*) SHAPE(image)
!    write(*,*) SHAPE(imagePadded)


        DO ii = borderPO, dims(3)
                DO jj = borderPO, dims(2)
                        DO kk = borderPO, dims(1)
                                accumulator = 0
                                DO ll = - border, border
                                        !$OMP PARALLEL
                                        DO mm = - border, border
                                                accumulator = accumulator + (kernel(ll + border + 1, mm + border + 1) &
                                                        * imagePadded(kk + ll, jj + mm, ii))
                                        END DO
                                        !$OMP END PARALLEL
                                END DO
                                result_image(kk - border, jj - border, ii - border) = INT(accumulator, KIND=ik)
                        END DO
                END DO
        END DO
END SUBROUTINE convolution_un_padded_input

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
