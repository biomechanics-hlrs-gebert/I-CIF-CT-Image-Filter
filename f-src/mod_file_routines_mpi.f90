!---------------------------------------------------------------------------------------------------
! mod_file_routines.f90
! Module with standard precision definitions
!
! \author Johannes Gebert
! \date 04.01.2021
! \date 25.04.2021

! subroutine mpi_err(ierr, mssg)
! SUBROUTINE check_file_exist(filename, todo, status)
! SUBROUTINE write_vtk (fl_un, fl_nm, array, spcng, dims, debug, debug_u)
! SUBROUTINE write_raw (fl_un, fl_m_nm, fl_d_nm, array, spcng, dims, debug, debug_u)
! SUBROUTINE read_vtk(fun, fl, array, dims, spcng, size, fov, bnds, log_un, status)
! SUBROUTINE read_raw(fun, fl, kind, type, dims, array, log_un, status)

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
 !---------------------------------------------------------------------------------------------------

 SUBROUTINE check_file_exist(filename, must_exist, mpi, status)

   INTEGER  (KIND=ik)    , INTENT(IN)                     :: must_exist
   CHARACTER(len=*)      , INTENT(IN)                     :: filename
   LOGICAL               , INTENT(IN)                     :: mpi
   INTEGER  (KIND=ik)    , INTENT(OUT), OPTIONAL          :: status
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
 !---------------------------------------------------------------------------------------------------

 SUBROUTINE write_vtk (fl_un, fl_nm, array, spcng, dims, debug_in, debug_u_in)

   ! It's HIGHLY recommended to check the existence of the output file prior to calling this
   ! Subroutine! Otherwise the program will crash. It's not double-checkd here, because this
   ! sequence often is placed at the very end of a program, which may run some time.

   INTEGER  (KIND=ik), INTENT(IN)                                    :: fl_un
   CHARACTER(len=*)                                                  :: fl_nm
   INTEGER  (KIND=ik), INTENT(IN), DIMENSION(:,:,:)                  :: array
   REAL     (KIND=rk), INTENT(IN), DIMENSION(3)                      :: spcng
   INTEGER  (KIND=ik), INTENT(IN), DIMENSION(3)                      :: dims
   !-- Debugging variables
   INTEGER  (KIND=ik)                               , OPTIONAL       :: debug_in    ! debug request
   INTEGER  (KIND=ik)                                                :: debug       ! debug request
   INTEGER  (KIND=ik)                               , OPTIONAL       :: debug_u_in
   INTEGER  (KIND=ik)                                                :: debug_u
   REAL     (KIND=rk)                                                :: start, end

   CALL CPU_TIME(start)

   IF (PRESENT(debug_in) .EQV. .FALSE.) THEN
      debug = 0_ik
   ELSE
      debug = debug_in
   END IF

   IF (PRESENT(debug_u_in) .EQV. .FALSE.) THEN
      debug_u = 6_ik
   ELSE
      debug_u = debug_u_in
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

   WRITE(UNIT=fl_un) INT( array(:,:,:), KIND=INT16)

   CLOSE(UNIT=fl_un)

   OPEN(UNIT=fl_un, FILE=TRIM(fl_nm), ACTION="write", STATUS="old", POSITION="append")

   WRITE(fl_un ,'(A)')          "METADATA"
   WRITE(fl_un ,'(A)')          "INFORMATION 0"
   WRITE(fl_un ,'(A)')

   CLOSE(UNIT=fl_un)

   CALL CPU_TIME(end)

   IF (debug >= 1_ik) THEN
      WRITE(debug_u,'(A)')        std_lnbrk
      WRITE(debug_u,'(A,A)')      "File name: ", TRIM(fl_nm)
      WRITE(debug_u,'(A,F9.4,A)') "Time to write file:                 ", end-start,  " seconds."
   END IF

 END SUBROUTINE write_vtk

 !---------------------------------------------------------------------------------------------------
 !---------------------------------------------------------------------------------------------------

 SUBROUTINE write_raw (fl_un, fl_m_nm, fl_d_nm, array, spcng, dims, debug, debug_u)

   ! It's HIGHLY recommended to check the existence of the output file prior to calling this
   ! Subroutine! Otherwise the program will crash. It's not double-checkd here, because this
   ! sequence often is placed at the very end of a program, which may run some time.

   INTEGER  (KIND=ik)    , INTENT(IN)                     :: fl_un
   CHARACTER(LEN=*)      , INTENT(IN)                     :: fl_m_nm, fl_d_nm
   REAL     (KIND=REAL64), INTENT(IN), DIMENSION(:,:,:)   :: array
   REAL     (KIND=REAL64),             DIMENSION(3)       :: spcng
   INTEGER  (KIND=ik)    , INTENT(IN), DIMENSION(3)       :: dims
   !-- internal variables
   INTEGER  (KIND=ik)                                     :: fl_m_un, fl_d_un
   !-- Debugging variables
   INTEGER  (KIND=ik)                               , OPTIONAL       :: debug    ! debug request
   INTEGER  (KIND=ik)                               , OPTIONAL       :: debug_u
   REAL     (KIND=rk)                                                :: start, end

   ! CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_APPEND + MPI_MODE_CREATE, MPI_INFO_NULL, fun, ierr)


   CALL CPU_TIME(start)

   IF (PRESENT(debug) .EQV. .FALSE.) THEN
      IF(PRESENT(debug_u) .EQV. .FALSE.) THEN
         debug   = 0_ik
         debug_u = 6_ik
      END IF
   END IF

   fl_m_un=fl_un
   fl_d_un=fl_un+1_ik

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

   CALL CPU_TIME(end)

   IF (debug >= 1) THEN
      WRITE(debug_u,'(A,A)')      "File name of meta data file: ", fl_m_nm
      WRITE(debug_u,'(A,A)')      "File name of      data file: ", fl_d_nm
      WRITE(debug_u,'(A,F9.4,A)') "Time to write file:          ", end-start,  " seconds."
      WRITE(debug_u,'(A)')        std_lnbrk
   END IF

END SUBROUTINE write_raw

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

SUBROUTINE read_vtk_meta(fun, filename, dims, spcng, typ, displacement, sze_o, fov_o, bnds_o, rd_o, status_o)
! log_un exists means "print log"!
! status  = 0 - everything is ok
! status /= 0 - Error
! status  = 1 - file does not exist
! status  = 2 - not a *.vtk file
! status  = 3 - file does not contain STRUCTURED_POINTS

INTEGER  (KIND=ik)                                                            :: fun
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

OPEN(UNIT=fun, FILE=TRIM(filename), STATUS="OLD")

hdr_lngth=0

DO ii=1,10
   READ(fun,'(A)') line
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
         typ      = tokens(3)
      END IF
   END IF !ntokens <0
END DO

CLOSE(fun)

fov = dims*spcng

IF (PRESENT(rd_o)) THEN
   WRITE(rd_o,'(A)')           std_lnbrk
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

END SUBROUTINE read_vtk_meta

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

SUBROUTINE read_raw_mpi(fun, filename, type, dims, subarray_dims, subarray, displacement, log_un, status_o)
! MPI Parallel read always reads subarrays.
! log_un exists means "print log"!
! type = "real" or "int"
! status  = 0 - everything is oka
! status /= 0 - Error
! status  = 1 - file does not exist
! status  = 2 - not a *.raw file
! status  = 3 - not a valid kind
! status  = 4 - not a valid type

INTEGER  (KIND=ik)                                                               :: fun
CHARACTER(len=*)                                                , INTENT(IN)     :: filename
CHARACTER(len=*)                                                , INTENT(IN)     :: type
INTEGER  (KIND=ik)    , DIMENSION(3)                            , INTENT(IN)     :: dims
INTEGER  (KIND=ik)    , DIMENSION(3)                            , INTENT(IN)     :: subarray_dims
INTEGER  (KIND=ik)    , DIMENSION (:,:,:), ALLOCATABLE          , INTENT(OUT)    :: subarray
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(IN)     :: displacement
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(IN)     :: log_un
INTEGER  (KIND=ik)                                    , OPTIONAL, INTENT(OUT)    :: status_o
!-- Internal Variables
INTEGER  (KIND=INT16) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_two
INTEGER  (KIND=INT32) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_four
INTEGER  (KIND=INT64) , DIMENSION (:,:,:), ALLOCATABLE                           :: array_i_eight
REAL     (KIND=REAL32), DIMENSION (:,:,:), ALLOCATABLE                           :: array_r_four
INTEGER  (KIND=ik)                                                               :: status, rd_o
INTEGER  (KIND=ik)                                                               :: kind
INTEGER  (KIND=MPI_OFFSET_KIND)                                                  :: hdr_lngth
INTEGER  (KIND=MPI_OFFSET_KIND)                                                  :: file_size
INTEGER  (KIND=ik)    , DIMENSION(3)                                             :: subarray_origin

! MPI
INTEGER  (KIND=ik)                                                               :: my_rank, size_mpi, ierr
INTEGER  (KIND=ik)                                                               :: type_subarray

IF (PRESENT(displacement)) hdr_lngth = displacement
IF (PRESENT(log_un))            rd_o = log_un

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_RANK couldn't be retrieved")

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_SIZE couldn't be retrieved")

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
dims                                , & ! Original array as all the addresses must fit
subarray_dims                       , &
subarray_origin - 1_mik             , & ! array_of_starts indexed from 0
MPI_ORDER_FORTRAN                   , &
MPI_INTEGER                         , &
type_subarray                       , &
ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

! OPEN(UNIT=fun, FILE=filename, CONVERT='LITTLE_ENDIAN', ACCESS="STREAM", FORM="UNFORMATTED", STATUS="OLD")
CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fun, ierr)

CALL MPI_FILE_GET_SIZE(fun, file_size, ierr)

kind = NINT( REAL((file_size-hdr_lngth), KIND=rk) / REAL(dims(1)*dims(2)*dims(3), KIND=rk))

! File Read

!-- Allocate memory an read array

! OPEN(UNIT=fun, FILE=fl, ACCESS="STREAM", FORM="UNFORMATTED", STATUS="OLD")

IF (TRIM(type) .EQ. "real") THEN

      IF(kind .EQ. 4_ik) THEN

         ALLOCATE( array_r_four( subarray_dims(1), subarray_dims(2), subarray_dims(3) ))
         READ(fun) array_r_four(:,:,:)    
         subarray = REAL(array_r_four, KIND=REAL64)
         DEALLOCATE(array_r_four)

      ELSE IF(kind .EQ. 8_ik) THEN

         ALLOCATE( subarray(subarray_dims(1), subarray_dims(2), subarray_dims(3) ))
         READ(fun) subarray(:,:,:)

      ELSE
         status = 3
      END IF

   ELSE IF (TRIM(type) .EQ. "int") THEN

      IF(kind .EQ. 2_ik) THEN

         ALLOCATE( array_i_two(subarray_dims(1), subarray_dims(2), subarray_dims(3)) )
         READ(UNIT=fun) array_i_two(:,:,:)
         subarray = REAL(array_i_two, KIND=REAL64)
         DEALLOCATE(array_i_two)
                     
      ELSE IF(kind .EQ. -2_ik) THEN

         ALLOCATE( array_i_two( subarray_dims(1), subarray_dims(2), subarray_dims(3)) )

         READ(UNIT=fun, POS=hdr_lngth) array_i_two(:,:,:)

         IF (MINVAL(array_i_two) .LT. 0_ik) THEN
            WRITE(rd_o,'(A)') 'INVALID INPUT - UNSIGNED SHORT PROBABLY COLLIDING WITH SIGNED INTEGERS. CHECK DATA.'
            status_o = 1_ik
         END IF 

         ! Input Copy/Pasted. Maybe reading of unsigned int will be implemented...
         subarray = INT(array_i_two, KIND=ik)
         DEALLOCATE(array_i_two)


      ELSE IF(kind .EQ. 4_ik) THEN

         CALL MPI_FILE_SET_VIEW( fun, &
         hdr_lngth,                   &
         MPI_INTEGER,                 &
         type_subarray,               &
         'native',                    &
         MPI_INFO_NULL,               &
         ierr)

         ALLOCATE( array_i_four( subarray_dims(1), subarray_dims(2), subarray_dims(3)) )
         READ(UNIT=fun) array_i_four(:,:,:)
         subarray = REAL(array_i_four, KIND=REAL64)
         DEALLOCATE(array_i_four)

      ELSE IF(kind .EQ. 8_ik) THEN

         ALLOCATE( array_i_eight( subarray_dims(1), subarray_dims(2), subarray_dims(3) ))
         READ(UNIT=fun) array_i_eight(:,:,:)
         subarray = REAL(array_i_eight, KIND=REAL64)
         DEALLOCATE(array_i_eight)

      ELSE
         status = 3
      END IF

      CALl  MPI_FILE_READ(fun, subarray, 1_mik, type_subarray, MPI_STATUS_IGNORE, ierr)
   ELSE
      status = 4
   END IF

   CLOSE(fun)

   CALL MPI_TYPE_FREE(type_subarray, ierr)

   CALL MPI_FILE_CLOSE(fun, ierr)

   IF (PRESENT(log_un) .AND. my_rank .EQ. 0_ik) THEN
      WRITE(log_un,'(A)')
      WRITE(log_un,'(A,A)')         "Input file                           ", TRIM(filename)
      WRITE(log_un,'(A)')           "Read raw module assumes Little-Endian while reading array!"
   END IF  ! print log output

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

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

SUBROUTINE tikz_std_plot (un, fl_nm, part, dims, tick, HU)
  !-- part =  1 Begin  a tikz picture
  !-- part = 98 Finish a plot (\end{axis})
  !-- part = 99 Finish a tikz picture
  !-- part =  2 write  a 2D plot header
  !-- part =  3 write  a 3D plot header
  !-- part =  5 write  a Bar plot with HU

  !-- Tikz export not always 100% precise. Debugging/Modifications on tikz source can be done afterwards.
  INTEGER(KIND=ik)                             , INTENT(IN)      :: un, part
  CHARACTER(len=*)                             , INTENT(IN)      :: fl_nm
  !-- OPTIONAL
  REAL     (KIND=rk) , DIMENSION(:,:), OPTIONAL, INTENT(IN)      :: dims        !-- (/ dims, min/max /)
  REAL     (KIND=rk)                 , OPTIONAL, INTENT(IN)      :: tick
  INTEGER  (KIND=ik) , DIMENSION(3)  , OPTIONAL, INTENT(IN)      :: HU
  !-- Internal variables
  INTEGER  (KIND=ik)                                             :: dd=2        !-- Default: 2D. Add 3D by hand if necessary


  IF (part .EQ. 1_ik) THEN
     OPEN(UNIT=un, FILE=TRIM(fl_nm), ACTION="write", STATUS="new")

     WRITE(un,'(A)') "\documentclass[border=5mm]{report}    % slides, article, book, ..."
     WRITE(un,'(A)') "\usepackage{tikz}"
     WRITE(un,'(A)') "\usepackage{graphicx}"
     WRITE(un,'(A)') "\usepackage{pgfplots}"
     WRITE(un,'(A)') "\usepackage{caption}"
     WRITE(un,'(A)') "\pgfplotsset{compat=1.8}"
     WRITE(un,'(A)') "\definecolor{hlrsblue}{RGB}{19, 176, 243}"
     WRITE(un,'(A)') "\usetikzlibrary{arrows.meta,shapes.misc}"
     WRITE(un,'(A)') "\begin{document}"

     !-- Define a coordinate system
     WRITE(un,'(A)') "\begin{figure}"
     WRITE(un,'(A)') "\centering"
     WRITE(un,'(A)') "\tikzset{cross/.style={cross out, draw=black, fill=none, minimum &
          &size=2*(#1-\pgflinewidth), inner sep=0pt, outer sep=0pt}, cross/.default={2pt}}"
  ELSE IF (part .EQ. 2_ik) THEN

  ELSE IF (part .EQ. 3_ik) THEN
     WRITE(un,'(A)') "\begin{tikzpicture}[]"
     WRITE(un,'(A)') "\begin{axis}["
     WRITE(un,'(A)') "%axis lines=center,"                          ! default: commented out
     WRITE(un,'(A)') "width=12cm,height=12cm,"

     IF (PRESENT(dims)) THEN
        IF(SIZE(dims) .EQ. 6_ik) dd=3_ik                                ! check whether plot is 3D

        WRITE(un,'(4(A,F6.3),A)') &
             &  "xmin=",dims(1,2) ,&
             & ",xmax=",dims(1,1) ,&
             & ",ymin=",dims(2,2) ,&
             & ",ymax=",dims(2,1), ","
        IF (dd .EQ. 3_ik)  THEN
           WRITE(un,'(2(A,F6.3),A)') &
                &  "zmin=",dims(3,2) ,&
                & ",zmax=",dims(3,1) ,","
        END IF
     END IF
     IF (PRESENT(tick)) THEN
        WRITE(un,'(3(A,F6.3),A)') "xtick={",-tick,", ", -tick/2._rk,",...,",tick,"},"
        WRITE(un,'(3(A,F6.3),A)') "ytick={",-tick,", ", -tick/2._rk,",...,",tick,"},"
        IF(dd .EQ. 3_ik) THEN
           WRITE(un,'(3(A,F6.3),A)') "ztick={",-tick,", ", -tick/2._rk,",...,",tick,"},"
        END IF
     END IF
     WRITE(un,'(A)') "minor x tick num=1,"                          ! default - maybe 4 is better - depends on plot
     WRITE(un,'(A)') "xlabel=$x (mm)$,"
     WRITE(un,'(A)') "ylabel=$y (mm)$,"
     WRITE(un,'(A)') "zlabel=$z (mm)$,"
     WRITE(un,'(A)') "title={Some random title},"
     WRITE(un,'(A)') "grid=major]"
  ELSE IF (part .EQ. 5_ik) THEN
     IF (PRESENT(HU)) THEN
        !-- Draw Colorbar Hounsfield units
        WRITE(un,'(A)') "\begin{tikzpicture}[]"
        WRITE(un,'(A)') "\begin{axis}["
        WRITE(un,'(A)') "hide axis,"
        WRITE(un,'(A)') "colormap={hlrscolormap}{rgb255=(255,255,255) rgb255=(19, 176, 243)},"
        WRITE(un,'(A)') "colorbar horizontal,"
        WRITE(un,'(A,I10,A)') "point meta min=",HU(1),","
        WRITE(un,'(A,I10,A)') "point meta max=",HU(3),","
        WRITE(un,'(A)') "colorbar style={"
        WRITE(un,'(A)') "width=12cm,"
        WRITE(un,'(A)') "xlabel=Hounsfield Units $H_U$,"
        WRITE(un,'(3(A,I10),A)') "xtick={",HU(1) ,",",HU(1) +HU(2)  ,",...,",HU(3) ,"}"
        WRITE(un,'(A)') "}]"
        WRITE(un,'(A)') "\node [ below] at (axis cs: 0,-2) {$(HU)$};"
        WRITE(un,'(A)') "\end{axis}"
        WRITE(un,'(A)') "\end{tikzpicture}"
     END IF
  ELSE IF (part .EQ. 98_ik) THEN
     WRITE(un,'(A)') "\end{axis}"
     WRITE(un,'(A)') "\end{tikzpicture}"
  ELSE IF (part .EQ. 99_ik) THEN
     !-- End of Tikz Figure
     WRITE(un,'(A)') "\caption{Please adjust the settings properly.} \label{fig:Prttp}"
     WRITE(un,'(A)') "\end{figure}"

     WRITE(un,'(A)') "\end{document}"

     CLOSE(un)

  END IF

END SUBROUTINE tikz_std_plot

END MODULE file_routines_mpi
