!-------------------------------------------------------------------------------------------------
! mod_histogram.f90
! Module to extract a histogram out a Vtk structure. The Module may inherit other subroutines as well.
!
! Author          Johannes Gebert
! Date Original   12.01.2021
! Date Modified   28.04.2021

! SUBROUTINE extract_histogram_scalar_array (array, hbnds, entries, histogram, fl_nm, ldng_rnk, debug, debug_u)

MODULE histogram

  USE ISO_FORTRAN_ENV
  USE standards

  IMPLICIT NONE

CONTAINS

   SUBROUTINE extract_histogram_scalar_array (array, hbnds, entries, histogram)

   REAL     (KIND=rk)    , DIMENSION(:,:,:)                       , INTENT(IN)        :: array
   INTEGER  (KIND=ik)    , DIMENSION(3)                           , INTENT(IN)        :: hbnds    ! histogram lower/upper bounds
   INTEGER  (KIND=ik)                                             , INTENT(IN)        :: entries
   INTEGER  (KIND=ik)    , DIMENSION(:)    , ALLOCATABLE          , INTENT(INOUT)     :: histogram

   ! Internal variables
   INTEGER  (KIND=ik)    , DIMENSION(65536)                                           :: histogram_result
   INTEGER  (KIND=ik)                                                                 :: fl_un = 50_ik
   INTEGER  (KIND=ik)                                                                 :: ii

   ! Calculate Notation of the histogram - done by slave processes to save time during communication
   ! Histogram is scaled to INT2 because 65.535 entries are sufficient to calculate and print any histogram
   IF( -3276700_ik > hbnds(1) .OR.  3276600 < hbnds(2)) THEN   ! bounds - 1 and *100 to scale
      recv_bffr_1D = recv_bffr_1D / &                          ! 32767 is slightly wrong, but neglibile
            REAL(ABS(hbnds(3)), KIND= REAL64) * 3276700._rk     ! ABS to preserve sign
   END IF

   histogram(:)       =0_mik
   histogram_result(:)=0_mik

   DO ii=1, SIZE(array)
      histogram(NINT(recv_bffr_1D(ii)/100.0_rk))=&
            histogram(NINT(recv_bffr_1D(ii)/100.0_rk))+1_mik
   END DO

   END SUBROUTINE extract_histogram_scalar_array

END MODULE histogram
