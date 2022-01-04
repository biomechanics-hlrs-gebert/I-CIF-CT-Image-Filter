!------------------------------------------------------------------------------
! MODULE: raw_binary
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief:
!> Module for reading/writing raw binary files parallely.
!> PureDat provides an extended funcitonality and calls a binary blob a
!> stream file (for example .int4.st).
!
!> @description
!> There's one routine for reading raw blobs for each datatype. 
!> And a switch for big/little endian. From the viewpoint of clean code, it is 
!> not well done. But creating subroutines for preparing a read raw and 
!> releasing/finishing read raw results in three calls and an allocate in the 
!> main program. Otherwise, it results in an mpi_read_raw_prepare subroutine 
!> with a specific data type, that itself must call the actual 
!> MPI_FILE_READ_ALL routine. Since the preparation is the biggest part, the 
!> calls won't be shorter significantly.
!------------------------------------------------------------------------------
MODULE raw_binary

USE ISO_FORTRAN_ENV
USE MPI
USE global_std
USE user_interaction

IMPLICIT NONE
INTERFACE mpi_read_raw
   MODULE PROCEDURE mpi_read_raw_rk4
   MODULE PROCEDURE mpi_read_raw_rk8
   MODULE PROCEDURE mpi_read_raw_ik4
   MODULE PROCEDURE mpi_read_raw_ik2
END INTERFACE mpi_read_raw

INTERFACE mpi_write_raw
   MODULE PROCEDURE mpi_write_raw_ik4
   MODULE PROCEDURE mpi_write_raw_ik2
END INTERFACE mpi_write_raw

INTERFACE ser_write_raw
   MODULE PROCEDURE ser_write_raw_ik2
   MODULE PROCEDURE ser_write_raw_ik4
END INTERFACE ser_write_raw

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: get_rank_section
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Get the section address/number of a specific domain
!
!> @param[in] domain No of the control volume. 
!> @param[in] sections x/y/z mesh of domains
!> @param[out] rank_section position of domain in x/y/z mesh
!------------------------------------------------------------------------------  
SUBROUTINE get_rank_section(domain, sections, rank_section)
  
INTEGER(KIND=ik), INTENT(IN)  :: domain
INTEGER(KIND=ik), DIMENSION(3), INTENT(IN)  :: sections
INTEGER(KIND=ik), DIMENSION(3), INTENT(OUT) :: rank_section

INTEGER(KIND=ik) :: rank, yrmndr, zrmndr ! remainder

!------------------------------------------------------------------------------
! Power of 2 is handled here, because with the algorithm of CASE DEFAULT, Greedy suboptimality kicks in!  
! Have a look at the corresponding Matlab/Octave testing file!
! In example at size_mpi = 128 Processors, where CASE DEFAULT will deliver 125 Processors!
!------------------------------------------------------------------------------
IF ( domain == 0_ik ) THEN
   rank_section = [ 1, 1, 1 ]
ELSE
   rank = domain + 1 ! MPI starts at 0

   !------------------------------------------------------------------------------
   ! Calculate the rank_section out of my_rank and sections [ x, y, z ]
   ! Tested via Octave. Not fully implemented by 20210503
   !------------------------------------------------------------------------------
   zrmndr = MODULO(rank, sections(1)*sections(2))
   IF (zrmndr == 0_ik) THEN
      rank_section = [ sections(1), sections(2), (rank - zrmndr) / (sections(1)*sections(2)) ]
   ELSE
      rank_section(3) = (rank - zrmndr) / (sections(1) * sections(2)) 
      yrmndr = MODULO(zrmndr, sections(1))

      IF (yrmndr == 0_ik) THEN
         rank_section = [ sections(1), (zrmndr - yrmndr) / sections(1), rank_section(3)+1 ]
      ELSE
         rank_section = [ yrmndr, (zrmndr - yrmndr) / sections(1) + 1, rank_section(3) + 1 ]
      END IF
   END IF
END IF

END SUBROUTINE get_rank_section


!------------------------------------------------------------------------------
! SUBROUTINE: mpi_read_raw_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Read the raw int2 data of a binary blob
!
!> @param[in] filename File name
!> @param[in] disp Length of the header (bytes)
!> @param[in] dims Amount of voxels per direction
!> @param[in] subarray_dims Amount of voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the data set
!> @param[out] subarray data
!------------------------------------------------------------------------------  
SUBROUTINE mpi_read_raw_ik2(filename, disp, dims, subarray_dims, subarray_origin, subarray)

CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(KIND=ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
INTEGER(KIND=INT16), DIMENSION (:,:,:), ALLOCATABLE, INTENT(OUT) :: subarray

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(KIND=mik) :: ierr, type_subarray, my_rank, size_mpi, fh
CHARACTER(LEN=scl) :: datarep

datarep = 'EXTERNAL32'

! Required to open files
CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, KIND=mik), INT(subarray_dims, KIND=mik), &
   INT(subarray_origin, KIND=mik), MPI_ORDER_FORTRAN, MPI_INTEGER2, type_subarray, ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER2, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)

ALLOCATE(subarray(subarray_dims(1), subarray_dims(2), subarray_dims(3)))

CALL MPI_FILE_READ_ALL(fh, subarray, INT(SIZE(subarray), KIND=mik), MPI_INTEGER2, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)

CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_read_raw_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: mpi_read_raw_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Read the raw int4 data of a binary blob
!
!> @param[in] filename File name
!> @param[in] disp Length of the header (bytes)
!> @param[in] dims Amount of voxels per direction
!> @param[in] subarray_dims Amount of voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the data set
!> @param[out] subarray data
!------------------------------------------------------------------------------  
SUBROUTINE mpi_read_raw_ik4(filename, disp, dims, subarray_dims, subarray_origin, subarray)

CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(KIND=ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
INTEGER(KIND=INT32), DIMENSION (:,:,:), ALLOCATABLE, INTENT(OUT) :: subarray

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(KIND=mik) :: ierr, type_subarray, my_rank, size_mpi, fh
CHARACTER(LEN=scl) :: datarep

datarep = 'EXTERNAL32'

! Required to open files
CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, KIND=mik), INT(subarray_dims, KIND=mik), & 
   INT(subarray_origin, KIND=mik), MPI_ORDER_FORTRAN, MPI_INTEGER4, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)

ALLOCATE(subarray(subarray_dims(1), subarray_dims(2), subarray_dims(3)))

CALL MPI_FILE_READ_ALL(fh, subarray, INT(SIZE(subarray), KIND=mik), MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)

CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_read_raw_ik4


!------------------------------------------------------------------------------
! SUBROUTINE: uik2_to_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Convert unsigned int2 to int 2 data.
!
!> @description
!> Fortran does not know this shit. Therefore a workaround...
!
!> @param[inout] subarray int data
!------------------------------------------------------------------------------  
SUBROUTINE uik2_to_ik2(subarray)

INTEGER(KIND=INT16), DIMENSION (:,:,:), INTENT(INOUT) :: subarray
INTEGER(KIND=ik) :: ii, jj, kk

! Not so pretty workaround
DO kk=1, SIZE(subarray,3)
DO jj=1, SIZE(subarray,2)
DO ii=1, SIZE(subarray,1)
   IF(subarray(ii,jj,kk)<=0) subarray(ii,jj,kk) = subarray(ii,jj,kk) + 65536_ik
END DO
END DO
END DO

subarray = subarray - 32768_ik

END SUBROUTINE uik2_to_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: mpi_read_raw_rk4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Read the raw real 4 (single precision) data of a binary blob
!
!> @param[in] filename File name
!> @param[in] disp Length of the header (bytes)
!> @param[in] dims Amount of voxels per direction
!> @param[in] subarray_dims Amount of voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the data set
!> @param[out] subarray data
!------------------------------------------------------------------------------  
SUBROUTINE mpi_read_raw_rk4(filename, disp, dims, subarray_dims, subarray_origin, subarray)

CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(KIND=ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
REAL(KIND=REAL32), DIMENSION (:,:,:), ALLOCATABLE, INTENT(OUT) :: subarray

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(KIND=mik) :: ierr, type_subarray, my_rank, size_mpi, fh
CHARACTER(LEN=scl) :: datarep

datarep = 'EXTERNAL32'

! Required to open files
CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
  
CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, KIND=mik), INT(subarray_dims, KIND=mik), & 
   INT(subarray_origin, KIND=mik), MPI_ORDER_FORTRAN, MPI_REAL, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_REAL, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)

ALLOCATE(subarray(subarray_dims(1), subarray_dims(2), subarray_dims(3)))

CALL MPI_FILE_READ_ALL(fh, subarray, INT(SIZE(subarray), KIND=mik), MPI_REAL, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)

CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_read_raw_rk4

!------------------------------------------------------------------------------
! SUBROUTINE: mpi_read_raw_rk8
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Read the raw real 8 (double precision) data of a binary blob
!
!> @param[in] fh File handle
!> @param[in] filename File name
!> @param[in] disp Length of the header (bytes)
!> @param[in] dims Amount of voxels per direction
!> @param[in] subarray_dims Amount of voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the data set
!> @param[out] subarray data
!------------------------------------------------------------------------------  
SUBROUTINE mpi_read_raw_rk8(filename, disp, dims, subarray_dims, subarray_origin, subarray)

CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(KIND=ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
REAL(KIND=REAL64), DIMENSION (:,:,:), ALLOCATABLE, INTENT(OUT) :: subarray

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(KIND=mik) :: ierr, type_subarray, my_rank, size_mpi, fh
CHARACTER(LEN=scl) :: datarep

datarep = 'EXTERNAL32'

! Required to open files
CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, KIND=mik), INT(subarray_dims, KIND=mik), & 
   INT(subarray_origin, KIND=mik), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_DOUBLE_PRECISION, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)

ALLOCATE(subarray(subarray_dims(1), subarray_dims(2), subarray_dims(3)))

CALL MPI_FILE_READ_ALL(fh, subarray, INT(SIZE(subarray), KIND=mik), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)

CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_read_raw_rk8


!------------------------------------------------------------------------------
! SUBROUTINE: mpi_write_raw_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Write raw binary data
!
!> @param[in] fh File handle
!> @param[in] disp Length of the header (bytes) - position to write to
!> @param[in] filename File name
!> @param[in] dims Voxels per direction
!> @param[in] subarray_dims Voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the subarray
!> @param[in] subarray Scalar field / Image data
!------------------------------------------------------------------------------  
 SUBROUTINE mpi_write_raw_ik2 (filename, disp, dims, subarray_dims, subarray_origin, subarray)
! type = 'int2', 'int4'
! IF type = uint2 - send an int4 and let it convert into int2 (!) Have a look at the src for details

CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(KIND=ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
INTEGER(KIND=INT16), DIMENSION (:,:,:), INTENT(IN) :: subarray

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(KIND=mik)  :: fh, ierr, type_subarray
CHARACTER(LEN=scl) :: datarep = 'EXTERNAL32'

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, KIND=mik), INT(subarray_dims, KIND=mik), & 
   INT(subarray_origin, KIND=mik), MPI_ORDER_FORTRAN, MPI_INTEGER2, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER2, type_subarray, TRIM(datarep), MPI_INFO_NULL, ierr)

CALL MPI_FILE_WRITE_ALL(fh, subarray, INT(SIZE(subarray), KIND=mik), MPI_INTEGER2, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)
CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_write_raw_ik2


!------------------------------------------------------------------------------
! SUBROUTINE: mpi_write_raw_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Write raw binary data
!
!> @description
!
!> @param[in] disp Length of the header (bytes) - position to write to
!> @param[in] filename File name
!> @param[in] dims Voxels per direction
!> @param[in] subarray_dims Voxels per direction of the subarray
!> @param[in] subarray_origin Physical origin of the subarray
!> @param[in] subarray Scalar field / Image data
!------------------------------------------------------------------------------  
 SUBROUTINE mpi_write_raw_ik4 (filename, disp, dims, subarray_dims, subarray_origin, subarray)
! type = 'int2', 'int4'
! IF type = uint2 - send an int4 and let it convert into int2 (!) Have a look at the src for details

CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: disp
INTEGER(KIND=ik),DIMENSION(3), INTENT(IN) :: dims, subarray_dims, subarray_origin
INTEGER(KIND=INT32), DIMENSION (:,:,:), INTENT(IN) :: subarray

! file handle fh is provided by mpi itself and mustn't be given by the program/call/user
INTEGER(KIND=mik)  :: fh, ierr, type_subarray
CHARACTER(LEN=scl) :: datarep = 'EXTERNAL32'

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), &
   MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, INT(dims, KIND=mik), INT(subarray_dims, KIND=mik), & 
   INT(subarray_origin, KIND=mik), MPI_ORDER_FORTRAN, MPI_INTEGER4, type_subarray,ierr)

CALL MPI_TYPE_COMMIT(type_subarray, ierr)

CALL MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, type_subarray, &
   TRIM(datarep), MPI_INFO_NULL, ierr)

CALL MPI_FILE_WRITE_ALL(fh, subarray, INT(SIZE(subarray), KIND=mik), &
   MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)

CALL MPI_TYPE_FREE(type_subarray, ierr)
CALL MPI_FILE_CLOSE(fh, ierr)

END SUBROUTINE mpi_write_raw_ik4


!------------------------------------------------------------------------------
! SUBROUTINE: ser_write_raw_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!------------------------------------------------------------------------------
SUBROUTINE ser_write_raw_ik2(fh, filename, array)

INTEGER(KIND=ik), INTENT(IN) :: fh
INTEGER(KIND=INT16), DIMENSION(:,:,:), INTENT(IN) :: array
CHARACTER(len=*), INTENT(IN) :: filename

OPEN (UNIT=fh, FILE=filename, ACCESS="STREAM", FORM="UNFORMATTED", &
   CONVERT='BIG_ENDIAN', STATUS="OLD", POSITION="APPEND")                                       
WRITE(UNIT=fh) array
CLOSE(UNIT=fh)

END SUBROUTINE ser_write_raw_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: ser_write_raw_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write raw binary data serially. 
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!------------------------------------------------------------------------------
SUBROUTINE ser_write_raw_ik4(fh, filename, array)

INTEGER(KIND=ik), INTENT(IN) :: fh
INTEGER(KIND=INT32), DIMENSION(:,:,:), INTENT(IN) :: array
CHARACTER(len=*), INTENT(IN) :: filename

OPEN (UNIT=fh, FILE=filename, ACCESS="STREAM", FORM="UNFORMATTED", &
   CONVERT='BIG_ENDIAN', STATUS="OLD", POSITION="APPEND")                                       
WRITE(UNIT=fh) array
CLOSE(UNIT=fh)

END SUBROUTINE ser_write_raw_ik4

END MODULE raw_binary


!------------------------------------------------------------------------------
! MODULE: vtk_meta_data
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief:
!> Module for reading/writing vtk-structured-point files
!------------------------------------------------------------------------------
MODULE vtk_meta_data

USE ISO_FORTRAN_ENV
USE global_std
USE strings
USE user_interaction

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: write_vtk_struct_points_header
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write a *.vtk structured points header
!
!> @param[in] fh File handle
!> @param[in] filename Name of the file
!> @param[in] type Data type
!> @param[in] spcng Distance between to voxels/scalars
!> @param[in] origin physical origin (mm)
!> @param[in] dims Amount of voxels/scalars per direction
!------------------------------------------------------------------------------
SUBROUTINE write_vtk_struct_points_header (fh, filename, type, spcng, origin, dims)

! It's HIGHLY recommended to check the existence of the output file prior to CALLing this
! Subroutine! Otherwise the program will crash. It's not double-checkd here, because this
! sequence often is placed at the very end of a program, which may run some time.

INTEGER  (KIND=ik), INTENT(IN) :: fh
CHARACTER(len=*)  , INTENT(IN) :: filename
CHARACTER(LEN=*)  , INTENT(IN) :: type
REAL     (KIND=rk), INTENT(IN), DIMENSION(3) :: spcng
REAL     (KIND=rk), INTENT(IN), DIMENSION(3) :: origin
INTEGER  (KIND=ik), INTENT(IN), DIMENSION(3) :: dims

CHARACTER(LEN=scl) :: datatype=''
LOGICAL :: exist

INQUIRE(UNIT=fh, exist=exist)

IF(.NOT. exist) THEN
   OPEN(UNIT=fh, FILE=TRIM(filename), ACTION='WRITE', STATUS='OLD', POSITION='APPEND')
END IF

SELECT CASE(TRIM(ADJUSTL(type)))
   CASE('uik2'); datatype = "unsigned_short"
   CASE('ik2'); datatype = "short"
   CASE('ik4'); datatype = "int"
   CASE('rk4'); datatype = "float"
   CASE('rk8'); datatype = "double"
   CASE DEFAULT
      CALL print_err_stop(fh, "No valid datatype given.", 1)
END SELECT

OPEN(UNIT=fh, FILE=TRIM(filename), ACTION='WRITE', STATUS='NEW')

WRITE(fh,'(A)')          "# vtk DataFile Version 4.2" ! Compatibility issue
WRITE(fh,'(A)')          "vtk output"
WRITE(fh,'(A)')          "BINARY"
WRITE(fh,'(A)')          "DATASET STRUCTURED_POINTS"
WRITE(fh,'(A,3(I5))')    "DIMENSIONS", dims
WRITE(fh,'(A,3(F11.6))') "SPACING ", spcng
WRITE(fh,'(A,3(F11.6))') "ORIGIN ", origin
WRITE(fh,'(A, I0)')      "POINT_DATA ", PRODUCT(INT(dims, KIND=INT64))
WRITE(fh,'(A)')          "SCALARS DICOMImage "//TRIM(ADJUSTL(datatype))
WRITE(fh,'(A)')          "LOOKUP_TABLE default"

CLOSE(UNIT=fh)
END SUBROUTINE write_vtk_struct_points_header

!------------------------------------------------------------------------------
! SUBROUTINE: write_vtk_struct_points_footer
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write a *.vtk structured points footer
!
!> @param[in] fh File handle
!> @param[in] filename File name
!------------------------------------------------------------------------------
SUBROUTINE write_vtk_struct_points_footer (fh, filename)

INTEGER(KIND=ik), INTENT(IN) :: fh
CHARACTER(len=*) :: filename

LOGICAL :: opened

INQUIRE(UNIT=fh, opened=opened)

IF(.NOT. opened) THEN
   OPEN(UNIT=fh, FILE=TRIM(filename), ACTION='WRITE', STATUS='OLD', POSITION='APPEND')
END IF

WRITE(fh ,'(A)')''
WRITE(fh ,'(A)')"METADATA"
WRITE(fh ,'(A)')"INFORMATION 0"
WRITE(fh ,'(A)')

CLOSE(UNIT=fh)
END SUBROUTINE write_vtk_struct_points_footer


!------------------------------------------------------------------------------
! SUBROUTINE: read_vtk_meta
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Read a *.vtk structured points header (footer not relevant)
!
!> @param[in] filename File name
!> @param[out] disp Displacement, length of header (bytes)
!> @param[out] dims Voxels per direction
!> @param[out] origin Physical (mm) origin
!> @param[out] spcng Physical distance between two voxels (mm)
!> @param[out] type Data type contained in binary blob
!------------------------------------------------------------------------------
SUBROUTINE read_vtk_meta(filename, disp, dims, origin, spcng, type)

CHARACTER(len=*)  , INTENT(IN)  :: filename
INTEGER  (KIND=ik), INTENT(OUT) :: disp
INTEGER  (KIND=ik), DIMENSION(3) , INTENT(OUT) :: dims
REAL     (KIND=rk), DIMENSION(3) , INTENT(OUT) :: origin
REAL     (KIND=rk), DIMENSION(3) , INTENT(OUT) :: spcng
CHARACTER(len=*), INTENT(OUT) :: type

!-- Initialize variables in case they're not used
INTEGER  (KIND=ik) :: ii=0, ntokens, fh

CHARACTER(len=mcl) :: line
CHARACTER(len=mcl) :: tokens(100)
CHARACTER(len=mcl), DIMENSION(3) :: token

!------------------------------------------------------------------------------
! Determine a new unit
!------------------------------------------------------------------------------  
fh = give_new_unit()
OPEN(UNIT=fh, FILE=TRIM(filename), STATUS="OLD")

disp=0

DO ii=1,10
   READ(fh,'(A)') line
   disp=disp+LEN(TRIM(line))+1_ik ! eol characters, white charcter
   
   CALL parse(str=line,delims=" ",args=tokens,nargs=ntokens)

   IF (ntokens > 0) THEN
   
      SELECT CASE(tokens(1))
         CASE('DIMENSIONS'); READ(tokens(2:4),'(I15)') dims(1:3)
         CASE('SPACING'); READ(tokens(2:4),'(F15.6)') spcng(1:3)  
         CASE('ORIGIN'); READ(tokens(2:4),'(F15.6)') origin(1:3)  
         CASE('DATASET')
            IF (tokens(2) /= "STRUCTURED_POINTS") THEN
               mssg = "The input file "//TRIM(filename)//" does not contain STRUCTURED_POINTS!"
               CALL print_err_stop(std_out, mssg, 1)
            END IF

         CASE('SCALARS')
            token(3) = tokens(3)

            SELECT CASE( TRIM( token(3) ) )
               CASE('float') ; type = 'rk4'
               CASE('double'); type = 'rk8'
               CASE('int')   ; type = 'ik4'
               CASE('short') ; type = 'ik2'
               CASE('unsigned_short'); type = 'uik2'
               CASE DEFAULT
                  WRITE(*,'(A)') "No valid type given in *.vtk File." 
            END SELECT
      END SELECT
   END IF !ntokens <0
END DO

CLOSE(fh)

END SUBROUTINE read_vtk_meta

END MODULE vtk_meta_data