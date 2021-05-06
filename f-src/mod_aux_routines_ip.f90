!-------------------------------------------------------------------------------------------------
! mod_aux_routines_IP.f90
! Module to extract a histogram out a Vtk structure. The Module may inherit other subroutines as well.
!
! Author          Johannes Gebert
! Date Original   12.01.2021
! Date Modified   28.04.2021

! SUBROUTINE extract_histogram_scalar_array (array, hbnds, entries, histogram)
! SUBROUTINE TD_Array_Scatter (size_mpi, sections)

MODULE aux_routines_ip

USE ISO_FORTRAN_ENV
USE standards

IMPLICIT NONE

CONTAINS

SUBROUTINE extract_histogram_scalar_array (array, hbnds, histogram)
! This is an inherently unflexible subroutine. It delivers exactly this kind of histogram and nothing else...

INTEGER  (KIND=ik)    , DIMENSION(:,:,:)                                           :: array
INTEGER  (KIND=ik)    , DIMENSION(3)                           , INTENT(IN)        :: hbnds    ! histogram lower/upper bounds
INTEGER  (KIND=ik)    , DIMENSION(:)    , ALLOCATABLE          , INTENT(OUT)       :: histogram

! Internal variables
INTEGER  (KIND=ik)                                                                 :: ii, jj, kk
INTEGER  (KIND=ik)    , DIMENSION(3)                                               :: shp

! Calculate Notation of the histogram - done by slave processes to save time during communication
! Histogram is scaled to INT2 because 65.535 entries are sufficient to calculate and print any histogram
IF( -3276700_ik > hbnds(1) .OR.  3276600 < hbnds(2)) THEN   ! bounds - 1 and *100 to scale
  array = array / REAL(ABS(hbnds(3)), KIND= REAL64) * 3276700._rk     ! ABS to preserve sign
END IF

ALLOCATE(histogram(-32768_ik:32767_ik))

shp = SHAPE(array)

DO ii=1, shp(1)
  DO jj=1, shp(2)
    DO kk=1, shp(3)
      histogram(INT(array(ii, jj, kk)/100_ik)) = histogram(INT(array(ii, jj, kk)/100_ik))+1_ik
    END DO
  END DO
END DO

END SUBROUTINE extract_histogram_scalar_array


SUBROUTINE TD_Array_Scatter (size_mpi, sections)
  ! Three-Dimensional Array Scatter Routine
  ! Calculate how many sections per direction are ideal to split an array via mpi
  ! It was written to scatter a 3D Array, however it may suit to differnt purposes
  ! IT IS HIGHLY RECOMMENDED TO PROVIDE A COMMON AMOUNT OF PROCESSORS AS SOME MAY BE UNUSED due to the nature of 3D arrays.
  
  INTEGER(kind=ik)                , INTENT(IN)  :: size_mpi
  INTEGER(kind=ik), DIMENSION(3)  , INTENT(OUT) :: sections

  ! Internal Variables
  INTEGER(KIND=ik)                              :: sw=0, true_size

  ! Power of 2 is handled here, because with the algorithm of CASE DEFAULT, Greedy suboptimality kicks in!  
  ! Have a look at the corresponding Matlab/Octave testing file!
  ! In example at size_mpi = 128 Processors, where CASE DEFAULT will deliver 125 Processors!
  SELECT CASE(size_mpi)
  CASE(2);     sections = (/  2,  1,  1 /)
  CASE(3);     sections = (/  3,  1,  1 /) ! Implementing it prevents crashing :)
  CASE(4);     sections = (/  2,  2,  1 /)
  CASE(5);     sections = (/  5,  1,  1 /) ! Implementing it prevents crashing :)
  CASE(6);     sections = (/  3,  2,  1 /)
  CASE(7);     sections = (/  7,  2,  1 /) ! Implementing it prevents crashing :)
  CASE(8);     sections = (/  2,  2,  2 /)
  CASE(10);    sections = (/  5,  2,  1 /)
  CASE(12);    sections = (/  3,  2,  2 /)
  CASE(16);    sections = (/  4,  2,  2 /)
  CASE(20);    sections = (/  5,  2,  2 /)
  CASE(24);    sections = (/  4,  3,  2 /)
  CASE(36);    sections = (/  4,  3,  3 /)
  CASE(40);    sections = (/  5,  4,  2 /)
  CASE(48);    sections = (/  4,  4,  3 /)
  CASE(64);    sections = (/  4,  4,  4 /)
  CASE(72);    sections = (/  6,  4,  3 /)
  CASE(80);    sections = (/  5,  4,  4 /)
  CASE(96);    sections = (/  6,  4,  4 /)
  CASE(128);   sections = (/  8,  4,  4 /)
  CASE(144);   sections = (/  6,  6,  4 /)
  CASE(216);   sections = (/  6,  6,  6 /)
  CASE(256);   sections = (/  8,  8,  4 /)
  CASE(512);   sections = (/  8,  8,  8 /)
  CASE(1024);  sections = (/ 16,  8,  8 /)
  CASE(2048);  sections = (/ 16, 16,  8 /)
  CASE(4096);  sections = (/ 16, 16, 16 /)
  CASE(8192);  sections = (/ 32, 16, 16 /)
  CASE(16384); sections = (/ 32, 32, 16 /)
  CASE(32768); sections = (/ 32, 32, 32 /)
  CASE DEFAULT
      ! If no option fits - a power of 2 will suffice. A lot of the processors may not get utilized then(!)
      true_size  = size_mpi - MODULO(size_mpi, 2_ik)
      sections   = (/ 1, 1, 1 /)
      DO WHILE (PRODUCT(sections) < true_size)
                                              sections(1) = sections(1) + 1_ik; sw=1_ik
            IF (PRODUCT(sections) < true_size) sections(2) = sections(2) + 1_ik; sw=2_ik
            IF (PRODUCT(sections) < true_size) sections(3) = sections(3) + 1_ik; sw=3_ik
      END DO
      sections(sw) = sections(sw) - 1_ik 
  END SELECT

END SUBROUTINE TD_Array_Scatter

END MODULE aux_routines_IP
