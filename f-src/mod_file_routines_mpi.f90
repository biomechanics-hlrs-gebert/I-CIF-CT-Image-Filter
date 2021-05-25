!---------------------------------------------------------------------------------------------------
! mod_file_routines.f90
! Module with standard precision definitions
!
! \author Johannes Gebert
! \date 04.01.2021
! \date 25.04.2021

! subroutine mpi_err(ierr, mssg)
! SUBROUTINE check_file_exist(filename, must_exist, mpi)
! SUBROUTINE write_vtk (fl_un, fl_nm, array, spcng, dims, debug_in, debug_u_in)
! SUBROUTINE write_raw (fl_un, fl_m_nm, fl_d_nm, array, spcng, dims, debug, debug_u)
! SUBROUTINE read_vtk_meta(fh, filename, dims, spcng, typ, displacement, sze_o, fov_o, bnds_o, rd_o, status_o)
! SUBROUTINE read_raw(fh, fl, kind, type, dims, array, log_un, status)
! SUBROUTINE write_matrix(matrix, title, u, frmwrk)

!-- Tex routines
! SUBROUTINE tikz_std_plot (un, fl_nm, part, dims, tick)
! SUBROUTINE export_voxel_star_tex(vs, star, fl_nm, debug, debug_u)
! SUBROUTINE export_tikz_vox_star_voxel(t ,u ,v ,vs , shp, cntr, hiR, fl_nm)

MODULE file_routines_mpi

USE ISO_FORTRAN_ENV
USE MPI
USE standards
USE strings

IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------------------------------
!-- Subroutine to evaluate allocation errors
!-- Copy and Pasted from struct-process by Dr.-Ing.Ralf Schneider (HLRS - Head of NUM)
!
subroutine mpi_err(ierr, mssg)

   !-- Dummy parameters
   integer(kind=mik), intent(in)    :: ierr
   character(len=*), intent(in)      :: mssg
   
   if (ierr /= MPI_SUCCESS) THEN
      write(*, "(100('!'))")
      write(*, '(A,I0,A)') 'MPI ERROR :', ierr,";"
      write(*, '(A)') trim(mssg)
      write(*, "(100('!'))")
      write(*, *) 'Exit ...'
      stop
    end if
   
end subroutine mpi_err

 !---------------------------------------------------------------------------------------------------

 SUBROUTINE check_file_exist(filename, must_exist, mpi)

   INTEGER  (KIND=ik)    , INTENT(IN)                     :: must_exist
   CHARACTER(len=*)      , INTENT(IN)                     :: filename
   LOGICAL               , INTENT(IN)                     :: mpi
   !-- Internal Variables
   LOGICAL                                                :: exist=.FALSE.
   INTEGER  (KIND=ik)                                     :: ierr

   INQUIRE (FILE = TRIM(filename), EXIST = exist)

   IF (exist .EQV. .TRUE.) THEN
      
      IF ( must_exist .EQ. 0_ik ) THEN
         WRITE(*,'(2A)') TRIM(filename), ' already exists! Program aborted.'; WRITE(*,'(A)')

         IF ( mpi .EQV. .TRUE. ) THEN
            CALL MPI_ABORT(MPI_COMM_WORLD, 28_mik, ierr) ! ierr is de facto a dummy
         ELSE
            STOP
         END IF
      END IF

   ELSE

      IF ( must_exist .EQ. 1_ik ) THEN
         WRITE(*,'(2A)') TRIM(filename), ' does not exist! Program aborted.'; WRITE(*,'(A)')

         IF ( mpi .EQV. .TRUE. ) THEN
            CALL MPI_ABORT(MPI_COMM_WORLD, 42_mik, ierr) ! ierr is de facto a dummy
         ELSE
            STOP
         END IF
      END IF

   END IF

 END SUBROUTINE check_file_exist

 !---------------------------------------------------------------------------------------------------

 SUBROUTINE write_vtk_meta (fh, filename, type, atStart, spcng, dims)

   ! It's HIGHLY recommended to check the existence of the output file prior to CALLing this
   ! Subroutine! Otherwise the program will crash. It's not double-checkd here, because this
   ! sequence often is placed at the very end of a program, which may run some time.

   INTEGER  (KIND=ik), INTENT(IN)                             :: fh
   CHARACTER(len=*)                                           :: filename
   CHARACTER(LEN=*)  , INTENT(IN)              , OPTIONAL     :: type
   LOGICAL           , INTENT(IN)                             :: atStart
   REAL     (KIND=rk), INTENT(IN), DIMENSION(3), OPTIONAL     :: spcng
   INTEGER  (KIND=ik), INTENT(IN), DIMENSION(3), OPTIONAL     :: dims

   IF (atStart .EQV. .TRUE.) THEN
      OPEN(UNIT=fh, FILE=TRIM(filename), ACTION='WRITE', STATUS='NEW')

      WRITE(fh,'(A)')            "# vtk DataFile Version 5.1"
      WRITE(fh,'(A)')            "vtk output"
      WRITE(fh,'(A)')            "BINARY"
      WRITE(fh,'(A)')            "DATASET STRUCTURED_POINTS"
      WRITE(fh,'(A,3(I5,A))')    "DIMENSIONS", dims(1)," ",dims(2)," ",dims(3),""
      WRITE(fh,'(A,3(F11.6,A))') "SPACING ", spcng(1)," ",spcng(2)," ",spcng(3),""
      WRITE(fh,'(A)')            "ORIGIN 0 0 0"
      WRITE(fh,'(A, I11)')       "POINT_DATA", dims(1)*dims(2)*dims(3)

      IF (TRIM(type) .EQ. 'int2') WRITE(fh,'(A)') "SCALARS DICOMImage short"    
      IF (TRIM(type) .EQ. 'int4') WRITE(fh,'(A)') "SCALARS DICOMImage int"

      WRITE(fh,'(A)')            "LOOKUP_TABLE default"
   ELSE
      OPEN(UNIT=fh, FILE=TRIM(filename), ACTION='WRITE', STATUS='OLD', POSITION='APPEND')
      ! WRITE(fh ,'(A)')          ''
      WRITE(fh ,'(A)')          "METADATA"
      WRITE(fh ,'(A)')          "INFORMATION 0"
      WRITE(fh ,'(A)')
   END IF

   CLOSE(UNIT=fh)
 END SUBROUTINE write_vtk_meta

 !---------------------------------------------------------------------------------------------------
 
 SUBROUTINE write_raw_mpi (type, hdr_lngth, filename, dims, subarray_dims, subarray_origin, subarray)
! type = 'int2', 'int4'

CHARACTER(LEN=*)                        , INTENT(IN)                         :: type
INTEGER  (KIND=MPI_OFFSET_KIND)                                              :: hdr_lngth
CHARACTER(LEN=*)                        , INTENT(IN)                         :: filename
INTEGER  (KIND=ik)   , DIMENSION(3)     , INTENT(IN)                         :: dims
INTEGER  (KIND=ik)   , DIMENSION(3)     , INTENT(IN)                         :: subarray_dims
INTEGER  (KIND=ik)   , DIMENSION(3)     , INTENT(IN)                         :: subarray_origin
INTEGER  (KIND=int16), DIMENSION (:,:,:)                                     :: subarray

! Internal Variables
INTEGER  (KIND=ik)                                                           :: fh

! MPI
INTEGER  (KIND=ik)                                                           :: my_rank, size_mpi, ierr
INTEGER  (KIND=ik)                                                           :: type_subarray

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_RANK couldn't be retrieved")

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_SIZE couldn't be retrieved")

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

! IF (TRIM(type) .EQ. 'int2') THEN
   CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
   dims                                , &
   subarray_dims                       , &
   subarray_origin - 1_mik             , &
   MPI_ORDER_FORTRAN                   , &
   MPI_INTEGER2                        , &
   type_subarray                       , &
   ierr)

   CALL MPI_TYPE_COMMIT(type_subarray, ierr)

   CALL MPI_FILE_SET_VIEW( fh , &
   hdr_lngth                  , &
   MPI_INTEGER2               , &
   type_subarray              , &
   'EXTERNAL32'               , &
   MPI_INFO_NULL              , &
   ierr)

   CALL MPI_FILE_WRITE_ALL(fh, subarray, SIZE(subarray), MPI_INTEGER2, MPI_STATUS_IGNORE, ierr)

! ELSE IF (TRIM(type) .EQ. 'int4') THEN
   ! CHANGE TYPE DEFINITION FIRST!

   ! CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
   ! dims                                , &
   ! subarray_dims                       , &
   ! subarray_origin - 1_mik             , &
   ! MPI_ORDER_FORTRAN                   , &
   ! MPI_INTEGER                         , &
   ! type_subarray                       , &
   ! ierr)

   ! CALL MPI_TYPE_COMMIT(type_subarray, ierr)

   ! CALL MPI_FILE_SET_VIEW( fh , &
   ! hdr_lngth                  , &
   ! MPI_INTEGER                , &
   ! type_subarray              , &
   ! 'EXTERNAL32'               , &
   ! MPI_INFO_NULL              , &
   ! ierr)

   ! CALL MPI_FILE_WRITE_ALL(fh, subarray, SIZE(subarray), MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
! END IF

CALL MPI_TYPE_FREE(type_subarray, ierr)
CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE write_raw_mpi

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

SUBROUTINE read_vtk_meta(fh, filename, dims, spcng, typ, displacement, sze_o, fov_o, bnds_o, rd_o, status_o)
! log_un exists means "print log"!
! status  = 0 - everything is ok
! status /= 0 - Error
! status  = 1 - file does not exist
! status  = 2 - not a *.vtk file
! status  = 3 - file does not contain STRUCTURED_POINTS

INTEGER  (KIND=ik)                                                            :: fh
CHARACTER(len=*)                                                , INTENT(IN)  :: filename
INTEGER  (KIND=ik)    , DIMENSION(3)                            , INTENT(OUT) :: dims
REAL     (KIND=rk)    , DIMENSION(3)                            , INTENT(OUT) :: spcng
CHARACTER(len=*)                                                , INTENT(OUT) :: typ
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(OUT) :: displacement
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(OUT) :: sze_o ! int32 mpi!
REAL     (KIND=rk)    , DIMENSION(3)                  , OPTIONAL, INTENT(OUT) :: fov_o
INTEGER  (KIND=ik)    , DIMENSION(3,2)                , OPTIONAL, INTENT(OUT) :: bnds_o
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(IN)  :: rd_o
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(OUT) :: status_o

!-- Initialize variables in case they're not used
INTEGER  (KIND=ik)                                                            :: sze
REAL     (KIND=rk)    , DIMENSION(3)                                          :: fov
INTEGER  (KIND=ik)    , DIMENSION(3,2)                                        :: bnds
INTEGER  (KIND=ik)                                                            :: status=0, ii=0, hdr_lngth, lui=6, ios, ntokens, knd

CHARACTER(len=mcl)                                                            :: line
CHARACTER(len=mcl)                                                            :: tokens(100)
CHARACTER(len=mcl)    , DIMENSION(3)                                          :: token

!-- Check existence of optional variables
IF (PRESENT(sze_o)   ) sze    = sze_o
IF (PRESENT(fov_o)   ) fov    = fov_o
IF (PRESENT(bnds_o)  ) bnds   = bnds_o
IF (PRESENT(status_o)) status = status_o
IF (PRESENT(rd_o)    ) lui=rd_o

OPEN(UNIT=fh, FILE=TRIM(filename), STATUS="OLD")

hdr_lngth=0

DO ii=1,10
   READ(fh,'(A)') line
   hdr_lngth=hdr_lngth+LEN(TRIM(line))+2_ik                   ! eol characters, whitechar
   CALL parse(str=line,delims=" ",args=tokens,nargs=ntokens)
   IF (ntokens > 0) THEN
      IF (tokens(1) .EQ. "DIMENSIONS") THEN
         CALL value_di(tokens(2), dims(1), ios=ios)
         CALL value_di(tokens(3), dims(2), ios=ios)
         CALL value_di(tokens(4), dims(3), ios=ios)
         bnds(:,1) = 1
         bnds(:,2) = dims(:)
         sze       = dims(1)*dims(2)*dims(3)
      ELSEIF (tokens(1) .EQ. "SPACING") THEN
         CALL value_dr(tokens(2), spcng(1), ios=ios)
         CALL value_dr(tokens(3), spcng(2), ios=ios)
         CALL value_dr(tokens(4), spcng(3), ios=ios)
      ELSEIF (tokens(1) .EQ. "DATASET") THEN
         IF (tokens(2) /= "STRUCTURED_POINTS") THEN
            WRITE(lui,'(3A)') "The input file ",filename," does not contain STRUCTURED_POINTS!"
            status = 3_ik
         ENDIF
      ELSEIF (tokens(1) .EQ. "ORIGIN") THEN
         IF (tokens(2)/="0" .OR. tokens(3)/="0" .OR. tokens(4)/="0") THEN
            WRITE(lui,'(A)') "Can't deal with origin /= (0 0 0) yet. Assumes this is not an issue."
         END IF
      ELSEIF (tokens(1) .EQ. "SCALARS") THEN
         !-- Get data type of the vtk-file
         token(2) = tokens(2)
         token(3) = tokens(3)
         IF (TRIM(token(2)) .EQ. "float"  .OR. TRIM(token(3)) .EQ. "float")  typ = 'real4'
         IF (TRIM(token(2)) .EQ. "double" .OR. TRIM(token(3)) .EQ. "double") typ = 'real8'
         IF (TRIM(token(2)) .EQ. "short"  .OR. TRIM(token(3)) .EQ. "short")  typ = 'int2'
         IF (TRIM(token(2)) .EQ. "int"    .OR. TRIM(token(3)) .EQ. "int")    typ = 'int4'
      END IF
   END IF !ntokens <0
END DO

CLOSE(fh)

fov = dims*spcng

IF (PRESENT(rd_o)) THEN
   WRITE(rd_o,'(A)')           "Input file"
   WRITE(rd_o,'(A)')           TRIM(filename)
   WRITE(rd_o,'(A)')           "Read vtk module assumes Big-Endian while reading array!"
   WRITE(rd_o,'(A)')
   WRITE(rd_o,'((A,I2))')      "The data type is KIND=                ", knd
   WRITE(rd_o,'(A,I5,A)')      "Header length                      ", hdr_lngth," Bytes"
   WRITE(rd_o,'(A,F8.3,A)')    "Resolution        - x               ", spcng(1)*1000._rk ," µm / Voxel"
   WRITE(rd_o,'(A,F8.3,A)')    "Resolution        - y               ", spcng(2)*1000._rk ," µm / Voxel"
   WRITE(rd_o,'(A,F8.3,A)')    "Resolution        - z               ", spcng(3)*1000._rk ," µm / Voxel"
   WRITE(rd_o,'((A,I5))')      "Voxels & bounds   - x              ", dims(1)
   WRITE(rd_o,'((A,I5))')      "Voxels & bounds   - y              ", dims(2)
   WRITE(rd_o,'((A,I5))')      "Voxels & bounds   - z              ", dims(3)
   WRITE(rd_o,'(A,F6.1,A)')    "Field of View     - x               ", fov(1) , " mm"
   WRITE(rd_o,'(A,F6.1,A)')    "Field of View     - y               ", fov(2) , " mm"
   WRITE(rd_o,'(A,F6.1,A)')    "Field of View     - z               ", fov(3) , " mm"
   WRITE(rd_o,'(A,I13,A)')     "Size of the internal array:",          sze, " Elements"
END IF  ! print log output

!-- Check existence of optional variables
IF (PRESENT(sze_o)   ) sze_o    = sze
IF (PRESENT(fov_o)   ) fov_o    = fov
IF (PRESENT(bnds_o)  ) bnds_o   = bnds
IF (PRESENT(status_o)) status_o = status
IF (PRESENT(displacement)) displacement = hdr_lngth

END SUBROUTINE read_vtk_meta

!---------------------------------------------------------------------------------------------------

SUBROUTINE read_raw_mpi(filename, type, hdr_lngth, dims, subarray_dims, subarray_origin, subarray, displacement, log_un, status_o)
! MPI Parallel read always reads subarrays.
! log_un exists means "print log"!
! type = 'real4', 'real8, 'int2', 'int4'

CHARACTER(LEN=*)                                                , INTENT(IN)     :: filename
CHARACTER(LEN=*)                                                , INTENT(IN)     :: type
INTEGER  (KIND=MPI_OFFSET_KIND)                                                  :: hdr_lngth
INTEGER  (KIND=ik)    , DIMENSION(3)                            , INTENT(IN)     :: dims
INTEGER  (KIND=ik)    , DIMENSION(3)                            , INTENT(IN)     :: subarray_dims
INTEGER  (KIND=ik)    , DIMENSION(3)                            , INTENT(IN)     :: subarray_origin
INTEGER  (KIND=ik)    , DIMENSION (:,:,:), ALLOCATABLE          , INTENT(OUT)    :: subarray
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(IN)     :: displacement
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(IN)     :: log_un
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(OUT)    :: status_o

! Internal Variables
INTEGER  (KIND=INT16) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_two
INTEGER  (KIND=INT32) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_four
REAL     (KIND=REAL64), DIMENSION (:,:,:), ALLOCATABLE                           :: array_r_eight
REAL     (KIND=REAL32), DIMENSION (:,:,:), ALLOCATABLE                           :: array_r_four
INTEGER  (KIND=ik)                                                               :: status=0, rd_o
INTEGER  (KIND=ik)                                                               :: fh

! MPI
INTEGER  (KIND=ik)                                                               :: my_rank, size_mpi, ierr
INTEGER  (KIND=ik)                                                               :: type_subarray

IF (PRESENT(displacement)) hdr_lngth = displacement
IF (PRESENT(log_un))            rd_o = log_un

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_RANK couldn't be retrieved")

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_SIZE couldn't be retrieved")

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

IF (TRIM(type) .EQ. 'real4') THEN

   CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
   dims                                , &
   subarray_dims                       , &
   subarray_origin - 1_mik             , &
   MPI_ORDER_FORTRAN                   , &
   MPI_REAL                            , &
   type_subarray                       , &
   ierr)

   CALL MPI_TYPE_COMMIT(type_subarray, ierr)

   CALL MPI_FILE_SET_VIEW( fh , &
   hdr_lngth                  , &
   MPI_REAL                   , &
   type_subarray              , &
   'EXTERNAL32'               , &
   MPI_INFO_NULL              , &
   ierr)

   ALLOCATE( array_r_four( subarray_dims(1), subarray_dims(2), subarray_dims(3) ))
   CALL MPI_FILE_READ_ALL(fh, array_r_four, SIZE(array_r_four), MPI_REAL, MPI_STATUS_IGNORE, ierr)

   subarray = INT(FLOOR(array_r_four), KIND=INT32)
   DEALLOCATE(array_r_four)

ELSE IF (TRIM(type) .EQ. 'real8') THEN

   CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
   dims                                , &
   subarray_dims                       , &
   subarray_origin - 1_mik             , &
   MPI_ORDER_FORTRAN                   , &
   MPI_DOUBLE_PRECISION                , &
   type_subarray                       , &
   ierr)

   CALL MPI_TYPE_COMMIT(type_subarray, ierr)

   CALL MPI_FILE_SET_VIEW( fh , &
   hdr_lngth                  , &
   MPI_DOUBLE_PRECISION       , &
   type_subarray              , &
   'EXTERNAL32'               , &
   MPI_INFO_NULL              , &
   ierr)

   ALLOCATE( array_r_eight( subarray_dims(1), subarray_dims(2), subarray_dims(3) ))
   CALL MPI_FILE_READ_ALL(fh, array_r_eight, SIZE(array_r_eight), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

   subarray = INT(FLOOR(array_r_four), KIND=INT32)
   DEALLOCATE(array_r_four)

   WRITE(rd_o,'(A)') 'WARNING: Converted real 8 to integer 4 during file read. Check validity.'

ELSE IF (TRIM(type) .EQ. 'int2') THEN

   CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
   dims                                , &
   subarray_dims                       , &
   subarray_origin - 1_mik             , &
   MPI_ORDER_FORTRAN                   , &
   MPI_INTEGER2                        , &
   type_subarray                       , &
   ierr)

   CALL MPI_TYPE_COMMIT(type_subarray, ierr)

   CALL MPI_FILE_SET_VIEW( fh , &
   hdr_lngth                  , &
   MPI_INTEGER2               , &
   type_subarray              , &
   'EXTERNAL32'               , &
   MPI_INFO_NULL              , &
   ierr)

   ALLOCATE( array_i_two( subarray_dims(1), subarray_dims(2), subarray_dims(3)) )
   CALL MPI_FILE_READ_ALL(fh, array_i_two, SIZE(array_i_two), MPI_INTEGER2, MPI_STATUS_IGNORE, ierr)
   
   subarray = INT(array_i_two, KIND=INT32)

   IF (MINVAL(array_i_two) .GT. 32767_ik) THEN
      WRITE(log_un,'(A)') 'WARNING: INVALID INPUT - UNSIGNED SHORT PROBABLY COLLIDING WITH SIGNED INTEGERS. CHECK DATA.'
      status = 1_ik
   END IF

   DEALLOCATE(array_i_two)

ELSE IF (TRIM(type) .EQ. 'int4') THEN

   CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
   dims                                , &
   subarray_dims                       , &
   subarray_origin - 1_mik             , &
   MPI_ORDER_FORTRAN                   , &
   MPI_INTEGER                         , &
   type_subarray                       , &
   ierr)

   CALL MPI_TYPE_COMMIT(type_subarray, ierr)

   CALL MPI_FILE_SET_VIEW( fh, &
   hdr_lngth                  , &
   MPI_INTEGER                , &
   type_subarray              , &
   'EXTERNAL32'               , &
   MPI_INFO_NULL              , &
   ierr)

   ALLOCATE( array_i_four( subarray_dims(1), subarray_dims(2), subarray_dims(3)) )
   CALL MPI_FILE_READ_ALL(fh, array_i_four, SIZE(array_i_four), MPI_INTEGER, MPI_STATUS_IGNORE, ierr)

   subarray = INT(array_i_four, KIND=INT16)

   IF (MINVAL(array_i_four) .LT. -32768_ik .OR. MAXVAL(array_i_four) .GT. 32767_ik) THEN
      WRITE(log_un,'(A)') 'WARNING: INVALID CONVERSION FROM INT4 TO INT2. CHECK DATA.'
      status = 1_ik
   END IF

   DEALLOCATE(array_i_four)
ELSE 
   status = 1_ik
END IF

CALL MPI_TYPE_FREE(type_subarray, ierr)
CALL MPI_FILE_CLOSE(fh, ierr)

IF (PRESENT(status_o)) status_o = status
END SUBROUTINE read_raw_mpi

!---------------------------------------------------------------------------------------------------

SUBROUTINE write_matrix(matrix, title, u, frmwrk)

  !-- Prints a matrix according to its dimensions
  REAL       (KIND=rk)  , DIMENSION(:,:)                          , INTENT(IN)     :: matrix
  CHARACTER  (LEN=*)                                              , INTENT(IN)     :: title
  INTEGER    (KIND=ik)                                            , INTENT(IN)     :: u
  CHARACTER  (LEN=*)                                    , OPTIONAL, INTENT(IN)     :: frmwrk
  !-- Internal variables
  INTEGER    (KIND=ik)  , DIMENSION(2)                                             :: shp
  INTEGER    (KIND=ik)                                                             :: kk
  CHARACTER  (LEN=mcl)                                                             :: frmwrk_used, fmt_std_out

  shp = SHAPE(matrix)

  IF (PRESENT(frmwrk)) THEN
     frmwrk_used = frmwrk
  ELSE
     frmwrk_used = "std_out"
  END IF

  IF ( TRIM(frmwrk_used) .EQ. "std_out" ) THEN

     WRITE(u, "( 10('-'), A, 10('-') )") TRIM(title)
     !-- Generate format
     WRITE (fmt_std_out, '(A, I3, A)' ) '(A,', shp(1),'(F8.2,A),A )'
     DO kk=1, shp(2)
        WRITE (u, fmt_std_out) '[', TRANSPOSE(matrix), ']'
     END DO

  ELSE IF ( TRIM(frmwrk_used) .EQ. "tex" ) THEN

     WRITE(u, '(A)') "The export of nxn dimensional matrices is not implemented at the moment."

  END IF

  END SUBROUTINE write_matrix

END MODULE file_routines_mpi