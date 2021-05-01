PROGRAM main
!-------------------------------
! Image processing
!
! Author: Benjamin Schnabel, M.Sc.
! Author: Johannes Gebert - HLRS - NUM
! GitHub: https://github.com/bennyschnabel/image_processing
! Date: 23.04.2021
!-------------------------------

USE ISO_FORTRAN_ENV
USE MPI_F08
USE file_routines     
USE kernels
USE strings
! histogram may be located in tis file as well. However it may extend to inherit additional routines
USE histogram
USE standards

IMPLICIT NONE

! For use with MPI (Calls), variables must be declared the same type like MPI was built.
! MPI: Kind=32 Bit / 4 Byte / ik=4 - Change in «working_directory/f-src/mod_standards.f90 

! Parameter
INTEGER  (KIND = ik), PARAMETER                                 :: fun1 = 5     ! File unit
INTEGER  (KIND = ik), PARAMETER                                 :: fun2 = 10    ! File unit

! Internal Variables
CHARACTER(LEN = mcl)                                            :: n2s, fileName, fileNameExportVtk
INTEGER  (KIND = ik)                                            :: i
INTEGER  (KIND = ik)           , DIMENSION(2)                   :: kernel_spec
CHARACTER(LEN = mcl)                                            :: selectKernel_str, sizekernel_str
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: dims, original_image_padding

REAL     (KIND = rk)                                            :: sigma, start, finish
CHARACTER(LEN = mcl)                                            :: sigma_str, version
REAL     (KIND = rk)           , DIMENSION(3)                   :: spcng
REAL     (KIND = rk)           , DIMENSION(:,:)  , ALLOCATABLE  :: kernel
INTEGER  (KIND = ik)           , DIMENSION(:,:,:), ALLOCATABLE  :: array, subarray

! Histogram Variables
CHARACTER(LEN = mcl)                                            :: histogram_filename_pre__Filter, & 
                                                                   histogram_filename_post_Filter
INTEGER  (KIND = ik), PARAMETER                                 :: fl_un_H_pre=41, fl_un_H_post=42
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: hist_boundaries
! Histogram is scaled to INT2 because 65.535 entries are sufficient to calculate and print virtually any histogram
INTEGER  (KIND = ik)           , DIMENSION(65536), ALLOCATABLE  :: histogram_pre__F       , histogram_post_F,     &
                                                                   histogram_pre__F_global, histogram_post_F_global

! MPI Variables
INTEGER  (KIND=mik)                                             :: ierr, my_rank, size_mpi, remainder
INTEGER  (KIND=mik)            , DIMENSION(:)    , ALLOCATABLE  :: send_cnt, dsplcmnts
REAL     (KIND=rk)             , DIMENSION(:)    , ALLOCATABLE  :: recv_bffr_1D
INTEGER  (KIND=mik)                                             :: recv_bffr_sz_1D
INTEGER  (KIND=mik)            , DIMENSION(3)                   :: sections, offset_per_dir, vox_per_dir_and_sec &
                                                                   dims_reduced, remainder_per_dir

! Debug Variables
INTEGER  (KIND = ik)                                            :: debug
INTEGER  (KIND = ik), PARAMETER                                 :: rd_o = 31    ! redirected StdOut
CHARACTER(LEN = mcl)                                            :: debug_str, log_file
LOGICAL                                                         :: log_exist = .FALSE.

!-------------------------------
! Initialize MPI Environment
!-------------------------------

CALL MPI_INIT(ierr)
CALL MPI_ERR(ierr,"MPI_INIT didn't succeed")

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_RANK couldn't be retrieved")

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_SIZE couldn't be retrieved")

IF (size_mpi < 2) THEN
        WRITE(rd_o,*)"At least two ranks required to execute this program."
        CLOSE(rd_o)
        CALL MPI_ABORT(MPI_COMM_WORLD, 1_mik, ierr)
END IF

! Start serialized sequence
IF (my_rank==0) THEN

        !-------------------------------
        ! Read input (from environment sh)
        !-------------------------------
        
        CALL GET_COMMAND_ARGUMENT(3, debug_str)
        READ(debug_str,'(I4)') debug

        !-- Create log-file
        CALL GET_COMMAND_ARGUMENT(2, log_file)

        INQUIRE(FILE=TRIM(log_file), EXIST=log_exist)

        IF ( log_exist .EQV. .TRUE. ) THEN
        WRITE(*,'(A)') "Log file already exists! Program aborted."
        CALL MPI_ABORT(MPI_COMM_WORLD, MPI_ERR_FILE_EXISTS, ierr)
        ELSE
        OPEN( UNIT = rd_o, file = TRIM(log_file), action="WRITE", status="new")
        ENDIF

        ! VTK file name
        CALL GET_COMMAND_ARGUMENT(1, fileName)
        ! Kernel for image processing ['0': Identity kernel, '1': Gaussian filter]
        CALL GET_COMMAND_ARGUMENT(4, selectKernel_str)
        READ(selectKernel_str,'(I4)') kernel_spec(1)
        ! Kernel size
        CALL GET_COMMAND_ARGUMENT(5, sizeKernel_str)
        READ(sizeKernel_str,'(I4)') kernel_spec(2)
        ! For gaussian filter kernel, sigma is required
        CALL GET_COMMAND_ARGUMENT(6, sigma_str)
        READ(sigma_str,'(F8.3)') sigma

        CALL GET_COMMAND_ARGUMENT(7, version)    

        !-------------------------------
        ! Calculation
        !-------------------------------

        ! Track calculation time
        CALL CPU_TIME(start)

        ! Titlescreen - Switched to Layout from standards.f90 - Can be changed without an issue
        WRITE(*,'(A)')  std_lnbrk
        WRITE(*,'(A)')  'Image Processing'
        WRITE(*,'(A)')  ''
        WRITE(*,'(A)')  'Author:  Benjamin Schnabel'
        WRITE(*,'(2A)') 'Version: ', TRIM(version)
        WRITE(*,'(A)')  std_lnbrk
        WRITE(*,'(A)')

ENDIF ! (my_rank==0)

! kernel_spec 0 (/ selectKernel, sizeKernel /) (in the first iteration of this program and for get_cmd_arg)
CALL MPI_BCAST (kernel_spec, 2_mik, MPI_INT             , 0_mik, MPI_COMM_WORLD)
CALL MPI_BCAST (sigma      , 1_mik, MPI_DOUBLE_PRECISION, 0_mik, MPI_COMM_WORLD)

ALLOCATE(kernel(kernel_spec(2), kernel_spec(2)))

SELECT CASE(kernel_spec(1))
      ! CASE(0) = identity Kernel, which also refers to CASE DEFAULT because it won't alter the image.
        CASE(1)
                IF (my_rank==0) WRITE(*,'(A)') 'Gaussian filter kernel selected'
                CALL kernel_gauss(kernel, kernel_spec(2), sigma)
        CASE DEFAULT
                IF (my_rank==0) WRITE(*,'(A)') 'Identity kernel selected'
                CALL kernel_identity(kernel, kernel_spec(2))
END SELECT

IF (my_rank==0) THEN
 
        WRITE(*,'(A)')

        ! Write kernel to terminal
        WRITE(n2s,*) sizeKernel
        WRITE(*,'(A)') 'Kernel:'
        DO i = 1, SIZE(kernel, 1)
                WRITE(*,'(' // TRIM(ADJUSTL(n2s)) //'F7.4)') kernel(i,:)
        END  DO
        WRITE(*,'(A)')

        ! Import VTK file
        WRITE(*,'(A)') 'VTK File import started'
        CALL read_vtk(fun1, fileName, array, dims, spcng)
        WRITE(*,'(A)') 'VTK File import done'


        ! Get the boundaries of the Histograms
        ! This procedure needs to be run on the full image. MPI Parallel file reading may be a solution
        hist_boundaries    = (/ NINT( MINVAL(array), KIND=INT32), NINT( MAXVAL(array), KIND=INT32) , 0_mik /)
        ! Apparent data Range
        hist_boundaries(3) = hist_boundaries(2)-hist_boundaries(1)
        
        ! Get the output Filenames of the Histograms
        histogram_filename_pre__Filter = TRIM(log_file(1:(LEN_TRIM(log_file)-4)))//"_hist_PRE__FILTER.csv"
        histogram_filename_post_Filter = TRIM(log_file(1:(LEN_TRIM(log_file)-4)))//"_hist_POST_FILTER.csv"

ENDIF ! (my_rank==0)

CALL MPI_BCAST (hist_boundaries, 3_mik, MPI_INT, 0_mik, MPI_COMM_WORLD)

! Get sections per direction
CALL TD_Array_Scatter (size_mpi, sections)

! Calculate Padding to decrease "size of array" to a corresponding size
! NOT the same as kernel_spec(2) !
border                 = FLOOR(REAL(kernel_spec(2)) / 2)
original_image_padding = border * 2

! Remainder per direction gets almost fully ignored (!) 
! It's assumed, that even if we split into 32768 processes (/ 32, 32, 32 /) sections,
! max. 31 x 31 x 31 Voxel get lost (Worst Case). Considering large input sets,
! which are imperative for utilizing large amounts of processors, losing 31 voxel at
! ~ 2000 ... 4000 Voxel per direction is not an issue.
! On the other hand, Distribution of the array gets way easier and presumably quicker.

remainder_per_dir = MODULO(array / sections)

! If remainder per direction is larger than the padding, simply shift the array to distribute to ranks
! into the center of array. Then add the border. This way, no artifical padding is required.
IF (remainder_per_dir(1) < original_image_padding(1) .OR. 
    remainder_per_dir(2) < original_image_padding(2) .OR.
    remainder_per_dir(3) < original_image_padding(3)) THEN
        remainder_per_dir = MODULO( (array - original_image_padding) / sections))
END IF

! Let's use the remainder to shift the filtered image to the center of the original one.
! This will ensure filtering either the whole image or at least filtering without an artificial 
! Padding at the outer surfaces (borders of the image).
     offset_per_dir = FLOOR(remainder_per_dir / 2_ik)

vox_per_dir_and_sec = ( (array-remainder_per_dir) / sections) + original_image_padding

dims_reduced        = dims - remainder_per_dir

DO ii = 1, sections(1) 
        DO jj = 1, sections(2) 
                DO kk = 1, sections(3) 
                        ! Converting address of subarray into rank is tested in Octave.
                        address = (kk-1)*sections(1)*sections(2) + (jj-1_ik)*sections(1) + ii
                        rank_section = sections
                        IF (address .EQ. my_rank) THEN
                                ! Add original image Data as padding
                                subarray_origin = sections * vox_per_dir_and_sec - border + offset_per_dir

                                ! MPI_TYPE_CREATE_DARRAY may fit better, however dealing with overlaps isn't clear
                                CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
                                        INT(SHAPE(array), KIND=mik) , &
                                        vox_per_dir_and_sec         , &
                                        subarray_origin - 1_mik     , & ! array_of_starts indexed from 0
                                        MPI_ORDER_FORTRAN           , &
                                        MPI_INT                     , &
                                        type_subarray               , &
                                        ierr)

                                CALL MPI_TYPE_COMMIT(type_subarray)
                        END IF
                END DO
        END DO
END DO 

ALLOCATE( subarray(vox_per_dir_and_sec) )

CALL MPI_SENDRECV (array, 1_mik, type_subarray, my_rank, my_rank, &
        subarray, 1_mik, type_subarray, 0_mik, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
     
CALL MPI_TYPE_FREE(type_subarray)

IF (my_rank == 0_ik) DEALLOCATE(array)

! subarray_reduced_boundaries
srb (1:3) = border + 1_ik
srb (4:6) = vox_per_dir_and_sec - border - 1_ik

! Prior to image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (subarray(srb(1):srb(4), srb(2):srb(5), srb(3):srb(6)), histogram_pre__F)            

! Start image processing
CALL convolution_un_padded_input (vox_per_dir_and_sec, subarray, kernel_spec(2), kernel, .TRUE., subarray)

! After image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (subarray, histogram_post_F)            

! Collect the subarrays to assemble a global vtk file.
! No Overlapping and no remainder.
vox_per_dir_and_sec = vox_per_dir_and_sec - original_image_padding
subarray_origin     = vox_per_dir_and_sec * rank_section

! MPI_TYPE_CREATE_DARRAY may fit better. Implementation not entirely clear at the moment...
CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
        INT(SHAPE(array), KIND=mik) , &
        vox_per_dir_and_sec         , &
        subarray_origin - 1_mik     , & ! array_of_starts indexed from 0
        MPI_ORDER_FORTRAN           , &
        MPI_INT                     , &
        type_subarray               , &
        ierr)

CALL MPI_TYPE_COMMIT(type_subarray)

ALLOCATE( result_array(dims_reduced) )

CALL MPI_SENDRECV (subarray, 1_mik, type_subarray, 0_mik, my_rank + size_mpi, &
        result_array, 1_mik, type_subarray, my_rank, my_rank + size_mpi, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
     
CALL MPI_TYPE_FREE(type_subarray) 

DEALLOCATE(subarray)

! Collect the data of the histogram pre filtering
CALL MPI_REDUCE (histogram_pre__F, histogram_global_pre__F, &
        65536_mik, MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

! Collect the data of the histogram post filtering
CALL MPI_REDUCE (histogram_pre__F, histogram_global_post__F, &
        65536_mik, MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

IF (my_rank==0) THEN

        ! Export Histogram of Scalar Array pre Filtering
        OPEN(UNIT = fl_un_H_pre, FILE=histogram_filename_pre__Filter, ACTION="WRITE", STATUS="new")
                WRITE(fl_un_H_pre,'(A)') "Scalar value / 100, Amount of Voxels per Scalar value"
                DO ii=-340,340
                        WRITE(fl_un_H_pre,'(I4,A,I18)') ii," , ",histogram_result_pre__F(ii)
                END DO
        CLOSE(fl_un_H_pre)
        
        ! Export Histogram of Scalar Array post Filtering
        OPEN(UNIT = fl_un_H_post, FILE=histogram_filename_psot_Filter, ACTION="WRITE", STATUS="new")
                WRITE(fl_un_H_post,'(A)') "Scalar value / 100, Amount of Voxels per Scalar value"
                DO ii=-340,340
                        WRITE(fl_un_H_post,'(I4,A,I18)') ii," , ",histogram_result_post_F(ii)
                END DO
        CLOSE(fl_un_H_post)

        ! Export VTK file (testing)
        WRITE(*,'(A)') 'VTK File export started'
        WRITE(n2s,*) selectKernel
        fileNameExportVtk = fileName(1:LEN_TRIM(fileName)-4) // '_Kernel_'// TRIM(ADJUSTL(n2s))  // '.vtk'
        CALL write_vtk(fun2, fileNameExportVtk, result_array, spcng, dims_reduced)
        WRITE(*,'(A)') 'VTK File export done'

        DEALLOCATE(result_array)
        DEALLOCATE(kernel)

        CALL CPU_TIME(finish)
        WRITE(*,*)
        WRITE(*,'(A7, F8.3, A8)') 'Time = ', (finish - start) / 60,' Minutes'
        WRITE(*,'(A)') std_lnbrk

ENDIF ! (my_rank==0)

CALL MPI_FINALIZE(ierr)
CALL MPI_ERR(ierr,"MPI_FINALIZE didn't succeed")

END PROGRAM main