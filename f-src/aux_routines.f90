!==============================================================================
!> \file aux_routines.f90
!> Module with standard precision definitions
!>
!> \author Johannes Gebert
!> \date 04.01.2021
!> \date 20.01.2021

MODULE aux_routines

USE kinds
USE ISO_FORTRAN_ENV
USE strings


IMPLICIT NONE

CONTAINS

 !---------------------------------------------------------------------------------------------------
 !---------------------------------------------------------------------------------------------------

 SUBROUTINE write_vtk (fl_un, fl_nm, array, spcng, dims)

   INTEGER  (KIND=INT64) , INTENT(IN)                     :: fl_un
   CHARACTER(len=*)                                       :: fl_nm
   REAL     (KIND=REAL64), INTENT(IN), DIMENSION(:,:,:)   :: array
   REAL     (KIND=REAL64), INTENT(IN), DIMENSION(3)       :: spcng
   INTEGER  (KIND=INT64) , INTENT(IN), DIMENSION(3)       :: dims
   !-- internal variables
   LOGICAL                                                :: exist=.FALSE.

   INQUIRE(FILE=fl_nm,EXIST=exist)

   IF (exist .EQV. .TRUE.) THEN
      WRITE(*,'(A)')  "-- An output file (*.vtk) already exists! Program aborted."
      WRITE(*,'(2A)') "-- ",TRIM(fl_nm)
      STOP
   END IF

   OPEN(UNIT=fl_un, FILE=TRIM(fl_nm), ACTION="write", STATUS="new")
   WRITE(fl_un,'(A)')            "# vtk DataFile Version 5.1"
   WRITE(fl_un,'(A)')            "vtk output"
   WRITE(fl_un,'(A)')            "BINARY"
   WRITE(fl_un,'(A)')            "DATASET STRUCTURED_POINTS"
   WRITE(fl_un,'(A,3(I5,A))')    "DIMENSIONS", dims(1)," ",dims(2)," ",dims(3),""
   WRITE(fl_un,'(A,3(F11.6,A))') "SPACING ", spcng(1)," ",spcng(2)," ",spcng(3),""
   WRITE(fl_un,'(A)')            "ORIGIN 0 0 0"
   WRITE(fl_un,'(A, I11)')       "POINT_DATA", dims(1)*dims(2)*dims(3)
   WRITE(fl_un,'(A)')            "SCALARS DICOMImage short"
   WRITE(fl_un,'(A)')            "LOOKUP_TABLE default"

   CLOSE(UNIT=fl_un)

   OPEN(UNIT=fl_un, FILE=TRIM(fl_nm),CONVERT='big_endian', &
        ACCESS="stream", FORM="unformatted", STATUS="old", POSITION="append")

   WRITE(UNIT=fl_un) INT( ANINT(array(:,:,:), KIND=REAL64), KIND=INT16)

   CLOSE(UNIT=fl_un)

   OPEN(UNIT=fl_un, FILE=TRIM(fl_nm), ACTION="write", STATUS="old", POSITION="append")

   WRITE(fl_un ,'(A)')          "METADATA"
   WRITE(fl_un ,'(A)')          "INFORMATION 0"
   WRITE(fl_un ,'(A)')

   CLOSE(UNIT=fl_un)

 END SUBROUTINE write_vtk

 !---------------------------------------------------------------------------------------------------
 !---------------------------------------------------------------------------------------------------

 SUBROUTINE write_raw (fl_un, fl_m_nm, fl_d_nm, array, spcng, dims)

   INTEGER  (KIND=INT64) , INTENT(IN)                     :: fl_un
   CHARACTER(LEN=*)      , INTENT(IN)                     :: fl_m_nm, fl_d_nm
   REAL     (KIND=REAL64), INTENT(IN), DIMENSION(:,:,:)   :: array
   REAL     (KIND=REAL64),             DIMENSION(3)       :: spcng
   INTEGER  (KIND=INT64) , INTENT(IN), DIMENSION(3)       :: dims
   !-- internal variables
   INTEGER  (KIND=INT64)                                  :: fl_m_un, fl_d_un
   LOGICAL                                                :: exist_m=.FALSE., exist_d=.FALSE.

   fl_m_un=fl_un
   fl_d_un=fl_un+1_ik

   INQUIRE(FILE=fl_m_nm, EXIST=exist_m)
   INQUIRE(FILE=fl_d_nm, EXIST=exist_d)

   IF (exist_m .EQV. .TRUE. .OR. exist_d .EQV. .TRUE.) THEN
      WRITE(*,'(A)')  "-- An output file (*.raw or *.meta) already exists! Program aborted."
      WRITE(*,'(2A)') "-- ", TRIM(fl_m_nm)
      WRITE(*,'(2A)') "-- ", TRIM(fl_d_nm)
      STOP
   END IF

   spcng=spcng/1000._rk

   OPEN(UNIT=fl_m_un, FILE=TRIM(fl_m_nm), ACTION="write", STATUS="new")

   WRITE(fl_m_un,'(A)')          "# raw DataFile"
   WRITE(fl_m_un,'(A)')          "raw output"
   WRITE(fl_m_un,'(A)')          "BINARY"
   WRITE(fl_m_un,'(A)')          "DATASET STRUCTURED_POINTS"
   WRITE(fl_m_un,'(A,3I5)')      "DIMENSIONS", dims(1)," ",dims(2)," ",dims(3),""
   WRITE(fl_m_un,'(A,3(F11.6))') "SPACING ", spcng(1)," ",spcng(2)," ",spcng(3),""
   WRITE(fl_m_un,'(A)')          "ORIGIN 0 0 0"
   WRITE(fl_m_un,'(A, I11)')     "POINT_DATA", dims(1)*dims(2)*dims(3)
   WRITE(fl_m_un,'(A)')          "# raw DataFile"
   WRITE(fl_m_un,'(A)')          "raw output"

   CLOSE(UNIT=fl_m_un)

   OPEN(UNIT=fl_d_un, FILE=TRIM(fl_d_nm), ACCESS="stream", FORM="unformatted", STATUS="new", POSITION="append")

   WRITE(UNIT=fl_d_un) INT( ANINT(array(:,:,:), KIND=REAL64), KIND=INT16)

   CLOSE(UNIT=fl_d_un)

END SUBROUTINE write_raw

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

SUBROUTINE read_vtk(fun, fl, array, dims, spcng, size_o, fov_o, bnds_o, log_un, status_o)
  ! log_un exists means "print log"!
  ! status  = 0 - everything is ok
  ! status /= 0 - Error
  ! status  = 1 - file does not exist
  ! status  = 2 - not a *.vtk file
  ! status  = 3 - file does not contain STRUCTURED_POINTS

  INTEGER  (KIND=INT64)                                           , INTENT(IN)     :: fun
  CHARACTER(len=*)                                                , INTENT(IN)     :: fl
  REAL     (KIND=REAL64), DIMENSION (:,:,:), ALLOCATABLE          , INTENT(OUT)    :: array
  INTEGER  (KIND=INT64) , DIMENSION(3)                            , INTENT(OUT)    :: dims
  REAL     (KIND=REAL64), DIMENSION(3)                            , INTENT(OUT)    :: spcng
  INTEGER  (KIND=INT32)                                 , OPTIONAL, INTENT(OUT)    :: size_o ! int32 mpi!
  REAL     (KIND=REAL64), DIMENSION(3)                  , OPTIONAL, INTENT(OUT)    :: fov_o
  INTEGER  (KIND=INT64) , DIMENSION(3,2)                , OPTIONAL, INTENT(OUT)    :: bnds_o
  INTEGER  (KIND=INT64)                                 , OPTIONAL, INTENT(IN)     :: log_un
  INTEGER  (KIND=INT64)                                 , OPTIONAL, INTENT(OUT)    :: status_o
  !-- Initialize variables in case they're not used
  INTEGER  (KIND=4_ik)                                                            :: sze
  REAL     (KIND=REAL64), DIMENSION(3)                                             :: fov
  INTEGER  (KIND=INT64) , DIMENSION(3,2)                                           :: bnds
  INTEGER  (KIND=INT64)                                                            :: status

  !-- Internal Variables
  INTEGER  (KIND=INT16) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_two
  INTEGER  (KIND=INT32) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_four
  REAL     (KIND=REAL32), DIMENSION (:,:,:), ALLOCATABLE                           :: array_r_four
  LOGICAL                                                                          :: exist, log_exist
  REAL     (KIND=REAL64)                                                           :: start, end
  CHARACTER(len=mcl)                                                               :: suf, line
  INTEGER  (KIND=INT64)                                                            :: ii, hdr_lngth, lui
  INTEGER  (KIND=INT64)                                                            :: file_size, kind
  INTEGER  (KIND=INT32)                                                            :: ios, ntokens
  CHARACTER(len=mcl)                                                               :: tokens(100)
  CHARACTER(len=mcl)    , DIMENSION(3)                                             :: token

  CALL CPU_TIME(start)

  !-- Check existence of optional variables
  IF (PRESENT(size_o)  ) sze   = size_o
  IF (PRESENT(fov_o)   ) fov    = fov_o
  IF (PRESENT(bnds_o)  ) bnds   = bnds_o
  IF (PRESENT(status_o)) status = status_o

  INQUIRE(FILE=TRIM(fl), EXIST=exist, SIZE=file_size)
  log_exist = PRESENT(log_un)

  IF (PRESENT(log_un) .EQV. .TRUE.) THEN
     lui = log_un     ! lui = log uni internal
  ELSE
     lui = 6_ik       ! sdtout of Fortran
  END IF

  !-- Read Meta
  IF (exist .EQV. .TRUE.) THEN

     suf=fl(LEN_TRIM(fl)-2 : LEN_TRIM(fl))

     IF(suf=="vtk") THEN

        OPEN(UNIT=fun, FILE=TRIM(fl), STATUS="OLD")

        hdr_lngth=0

        DO ii=1,10
           READ(fun,'(A)') line
           hdr_lngth=hdr_lngth+LEN(TRIM(line))+2_ik                   ! eol characters, whitechar
           CALL parse(str=line,delims=" ",args=tokens,nargs=ntokens)
           IF (ntokens > 0) THEN
              IF (tokens(1)=="DIMENSIONS") THEN
                 CALL value_di(tokens(2), dims(1), ios=ios)
                 CALL value_di(tokens(3), dims(2), ios=ios)
                 CALL value_di(tokens(4), dims(3), ios=ios)
                 bnds(:,1) = 1
                 bnds(:,2) = dims(:)
                 sze      = dims(1)*dims(2)*dims(3)
              ELSEIF (tokens(1)=="SPACING") THEN
                 CALL value_dr(tokens(2), spcng(1), ios=ios)
                 CALL value_dr(tokens(3), spcng(2), ios=ios)
                 CALL value_dr(tokens(4), spcng(3), ios=ios)
              ELSEIF (tokens(1)=="DATASET") THEN
                 IF (tokens(2) /= "STRUCTURED_POINTS") THEN
                    WRITE(lui,'(3A)') "The input file ",fl," does not contain STRUCTURED_POINTS!"
                    status = 3_ik
                 ENDIF
              ELSEIF (tokens(1)=="ORIGIN") THEN
                 IF (tokens(2)/="0" .OR. tokens(3)/="0" .OR. tokens(4)/="0") THEN
                    WRITE(lui,'(A)') "Can't deal with origin /= (0 0 0) yet. Assumes this is not an issue."
                 END IF
              ELSEIF (tokens(1)=="SCALARS") THEN
                 !-- Get data type of the vtk-file
                 token(2) = tokens(2)
                 token(3) = tokens(3)
              END IF
           END IF !ntokens <0
        END DO

        CLOSE(fun)

        fov = dims*spcng

        !-- Allocate memory an read array
        kind = INT( ANINT( REAL((file_size-hdr_lngth), KIND=REAL64) / REAL(sze, KIND=REAL64), KIND=INT64))

        OPEN(UNIT=fun, FILE=fl, CONVERT='BIG_ENDIAN', ACCESS="STREAM", FORM="UNFORMATTED", STATUS="OLD")

        IF (TRIM(token(2)) == "float" .OR. TRIM(token(3)) == "float" ) THEN

           IF(kind==4_ik) THEN
              ALLOCATE( array_r_four(dims(1),dims(2),dims(3)))
              READ(UNIT=fun, POS=hdr_lngth) array_r_four(:,:,:)
              array = REAL(array_r_four, KIND=REAL64)
              DEALLOCATE(array_r_four)
           END IF

        ELSE IF (TRIM(token(2)) == "double" .OR. TRIM(token(3)) == "double" ) THEN

           IF(kind==8_ik) THEN
              ALLOCATE(array(dims(1),dims(2),dims(3)))
              READ(UNIT=fun, POS=hdr_lngth) array(:,:,:)
           END IF

        ELSE IF (TRIM(token(2)) == "short" .OR. TRIM(token(3)) == "short") THEN

           ALLOCATE(array_i_two(dims(1),dims(2),dims(3)))
           READ(UNIT=fun, POS=hdr_lngth) array_i_two(:,:,:)
           array = REAL(array_i_two, KIND=REAL64)
           DEALLOCATE(array_i_two)

        ELSE IF (TRIM(token(2)) == "int" .OR. TRIM(token(3)) == "int") THEN
           WRITE(*,*) 'dims(1)', dims(1)
           WRITE(*,*) 'dims(2)', dims(2)
           WRITE(*,*) 'dims(3)', dims(3)
           WRITE(*,*) 'hdr_lngth', hdr_lngth
           ALLOCATE(array_i_four(dims(1),dims(2),dims(3)))
           READ(UNIT=fun, POS=hdr_lngth) array_i_four(:,:,:)
           !WRITE(*,*) 'MAXVAL:', MAXVAL(array_i_four)
           !WRITE(*,*) 'MINVAL:', MINVAL(array_i_four)
           array = REAL(array_i_four, KIND=REAL64)
           DEALLOCATE(array_i_four)

        ELSE
           status = 4
           WRITE(lui,'(3A)')    "The data type of the input file ",fl," was not identified."
           WRITE(lui,'(A)')    "Check header/input file. Program aborted."
        END IF ! if token ==  for data type, allocation and reading

        CLOSE(fun)

        IF (PRESENT(log_un)) THEN
           WRITE(lui,'(A)')
           WRITE(lui,'(A)')           "Input file"
           WRITE(lui,'(A)')           TRIM(fl)
           WRITE(lui,'(A)')           "Read vtk module assumes Big-Endian while reading array!"
           WRITE(lui,'(A)')
           WRITE(lui,'((A,I2))')      "The data type is KIND=                ", kind
           WRITE(lui,'(A,I5,A)')      "Header length                      ", hdr_lngth," Bytes"
           WRITE(lui,'(A,F8.3,A)')    "Resolution        - x               ", spcng(1)*1000._rk ," µm / Voxel"
           WRITE(lui,'(A,F8.3,A)')    "Resolution        - y               ", spcng(2)*1000._rk ," µm / Voxel"
           WRITE(lui,'(A,F8.3,A)')    "Resolution        - z               ", spcng(3)*1000._rk ," µm / Voxel"
           WRITE(lui,'((A,I5))')      "Voxels & bounds   - x              ", dims(1)
           WRITE(lui,'((A,I5))')      "Voxels & bounds   - y              ", dims(2)
           WRITE(lui,'((A,I5))')      "Voxels & bounds   - z              ", dims(3)
           WRITE(lui,'(A,F6.1,A)')    "Field of View     - x               ", fov(1) , " mm"
           WRITE(lui,'(A,F6.1,A)')    "Field of View     - y               ", fov(2) , " mm"
           WRITE(lui,'(A,F6.1,A)')    "Field of View     - z               ", fov(3) , " mm"
           WRITE(lui,'(A,I13,A)')     "Size of the internal array:",          sze, " Elements"
           CALL CPU_TIME(end)             ! deliberately put here to get most realistic impression
           WRITE(lui,'(A,F9.4)')      "Time to read file:               ", end-start
        END IF  ! print log output

     ELSE
        status = 2_ik
     END IF
  ELSE
     status = 1_ik
  END IF    ! exist == true

END SUBROUTINE read_vtk

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

SUBROUTINE read_raw(fun, fl, kind, type, dims, array, log_un, status)
  ! log_un exists means "print log"!
  ! type = "real" or "int"
  ! status  = 0 - everything is ok
  ! status /= 0 - Error
  ! status  = 1 - file does not exist
  ! status  = 2 - not a *.raw file
  ! status  = 3 - not a valid kind
  ! status  = 4 - not a valid type

  INTEGER  (KIND=INT64)                                           , INTENT(IN)     :: fun
  CHARACTER(len=*)                                                , INTENT(IN)     :: fl
  INTEGER  (KIND=INT64)                                           , INTENT(IN)     :: kind
  CHARACTER(len=*)                                                , INTENT(IN)     :: type
  INTEGER  (KIND=INT64) , DIMENSION(3)                            , INTENT(IN)     :: dims
  REAL     (KIND=REAL64), DIMENSION (:,:,:), ALLOCATABLE          , INTENT(OUT)    :: array
  INTEGER  (KIND=INT64)                                 , OPTIONAL, INTENT(IN)     :: log_un
  INTEGER  (KIND=INT64)                                 , OPTIONAL, INTENT(OUT)    :: status
  !-- Internal Variables
  INTEGER  (KIND=INT16) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_two
  INTEGER  (KIND=INT32) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_four
  INTEGER  (KIND=INT64) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_eight
  REAL     (KIND=REAL32), DIMENSION (:,:,:), ALLOCATABLE                           :: array_r_four
  LOGICAL                                                                          :: exist, log_exist
  REAL     (KIND=REAL64)                                                           :: start, end
  CHARACTER(len=mcl)                                                               :: suf
  INTEGER  (KIND=INT64)                                                            :: ii, lui
  INTEGER  (KIND=INT64)                                                            :: file_size

  CALL CPU_TIME(start)

  !-- General
  INQUIRE(FILE=TRIM(fl), EXIST=exist, SIZE=file_size)
  log_exist = PRESENT(log_un)

  IF (log_exist .EQV. .TRUE.) THEN
     lui = log_un     ! lui = log uni internal
  ELSE
     lui = 6_ik       ! sdtout of Fortran
  END IF

  !-- Check Meta
  IF (exist .EQV. .TRUE.) THEN
     suf=fl(LEN_TRIM(fl)-2 : LEN_TRIM(fl))
     IF(suf/="raw") THEN
        status = 2_ik
     END IF

     !-- Allocate memory an read array

     OPEN(UNIT=fun, FILE=fl, ACCESS="STREAM", FORM="UNFORMATTED", STATUS="OLD")

     IF (TRIM(type)=="real") THEN

        IF(kind==4_ik) THEN
           ALLOCATE( array_r_four(dims(1),dims(2),dims(3)))
           READ(fun) array_r_four(:,:,:)
           array = REAL(array_r_four, KIND=REAL64)
           DEALLOCATE(array_r_four)
        ELSE IF(kind==8_ik) THEN
           ALLOCATE(array(dims(1),dims(2),dims(3)))
           READ(fun) array(:,:,:)
        ELSE
           status = 3
        END IF

     ELSE IF (TRIM(type)=="int") THEN

        IF(kind==2_ik) THEN
           ALLOCATE(array_i_two(dims(1),dims(2),dims(3)))
           READ(UNIT=fun) array_i_two(:,:,:)
           array = REAL(array_i_two, KIND=REAL64)
           DEALLOCATE(array_i_two)
        ELSE IF(kind==4_ik) THEN
           ALLOCATE(array_i_four(dims(1),dims(2),dims(3)))
           READ(UNIT=fun) array_i_four(:,:,:)
           array = REAL(array_i_four, KIND=REAL64)
           DEALLOCATE(array_i_four)
        ELSE IF(kind==8_ik) THEN
           ALLOCATE(array_i_eight(dims(1),dims(2),dims(3)))
           READ(UNIT=fun) array_i_eight(:,:,:)
           array = REAL(array_i_eight, KIND=REAL64)
           DEALLOCATE(array_i_eight)
        ELSE
           status = 3

        END IF

     ELSE
        status = 4
     END IF

     CLOSE(fun)

     IF (lui==log_un) THEN
        ! All native meta data - including dims - shall be logged by another instance
        WRITE(lui,'(A)')
        WRITE(lui,'(A,A)')         "Input file                           ", TRIM(fl)
        WRITE(lui,'(A)')           "Read raw module assumes Little-Endian while reading array!"
        CALL CPU_TIME(end)             ! deliberately put here to get most realistic impression
        WRITE(lui,'(A,F9.4)')      "Time to read file:               ", end-start
     END IF  ! print log output

ELSE
  status = 1_ik
END IF    ! exist == true

END SUBROUTINE read_raw

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

END MODULE aux_routines
 !---------------------------------------------------------------------------------------------------
