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
USE aux_routines_ip
USE standards

IMPLICIT NONE

! For use with MPI (Calls), variables must be declared the same type like MPI was built.
! MPI: Kind=32 Bit / 4 Byte / ik=4 - Change in «working_directory/f-src/mod_standards.f90 

! Parameter
INTEGER  (KIND = ik), PARAMETER                                 :: debug = 2
INTEGER  (KIND = ik), PARAMETER                                 :: fun1 = 5     ! File unit
INTEGER  (KIND = ik), PARAMETER                                 :: fun2 = 10    ! File unit

! Internal Variables
CHARACTER(LEN = mcl)                                            :: n2s, fileName, fileNameExportVtk
INTEGER  (KIND = ik)                                            :: ii, jj, kk, address
INTEGER  (KIND = ik)           , DIMENSION(2)                   :: kernel_spec
CHARACTER(LEN = mcl)                                            :: selectKernel_str, sizekernel_str

INTEGER  (KIND = ik)           , DIMENSION(3)                   :: dims, original_image_padding, border, subarray_origin, clipboard
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: sections, rank_section, offset_per_dir, vox_per_dir_and_sec
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: dims_reduced, remainder_per_dir

INTEGER  (KIND = ik)           , DIMENSION(6)                   :: srb ! subarray_reduced_boundaries
REAL     (KIND = rk)                                            :: sigma, start, finish
CHARACTER(LEN = mcl)                                            :: sigma_str, version
REAL     (KIND = rk)           , DIMENSION(3)                   :: spcng
REAL     (KIND = rk)           , DIMENSION(:,:)  , ALLOCATABLE  :: kernel
INTEGER  (KIND = ik)           , DIMENSION(:,:,:), ALLOCATABLE  :: array, result_array, subarray, result_subarray 
CHARACTER(LEN =   8)                                            :: date
CHARACTER(LEN =  10)                                            :: time

! Histogram Variables
CHARACTER(LEN = mcl)                                            :: histogram_filename_pre__Filter 
CHARACTER(LEN = mcl)                                            :: histogram_filename_post_Filter
INTEGER  (KIND = ik), PARAMETER                                 :: fl_un_H_pre=41, fl_un_H_post=42
INTEGER  (KIND = ik)                                            :: histo_bound_global_lo, histo_bound_global_hi
INTEGER  (KIND = ik)                                            :: histo_bound_local_lo, histo_bound_local_hi
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: hbnds

! Histogram is scaled to INT2 because 65.535 entries are sufficient to calculate and print virtually any histogram
INTEGER  (KIND = ik)           , DIMENSION(:)    , ALLOCATABLE  :: histogram_pre__F       , histogram_post_F
INTEGER  (KIND = ik)           , DIMENSION(:)    , ALLOCATABLE  :: histogram_pre__F_global, histogram_post_F_global

! MPI Variables
INTEGER  (KIND = mik)                                           :: ierr, my_rank, size_mpi
TYPE(MPI_DATATYPE)                                              :: type_subarray, type_result_subarray

! Debug Variables
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
        CALL GET_COMMAND_ARGUMENT(3, selectKernel_str)
        READ(selectKernel_str,'(I4)') kernel_spec(1)
        ! Kernel size
        CALL GET_COMMAND_ARGUMENT(4, sizeKernel_str)
        READ(sizeKernel_str,'(I4)') kernel_spec(2)
        ! For gaussian filter kernel, sigma is required
        CALL GET_COMMAND_ARGUMENT(5, sigma_str)
        READ(sigma_str,'(F8.3)') sigma

        CALL GET_COMMAND_ARGUMENT(6, version)    

        !-------------------------------
        ! Calculation
        !-------------------------------

        ! Track calculation time
        CALL CPU_TIME(start)

        ! Titlescreen - Switched to Layout from standards.f90 - Can be changed without an issue
        WRITE(*,'(A)')  std_lnbrk
        WRITE(*,'(A)')  'Image Processing'
        WRITE(*,'(A)')  ''
        WRITE(*,'(A)')  'Authors: Benjamin Schnabel, Johannes Gebert | HLRS - NUM'
        WRITE(*,'(2A)') 'Version: ', TRIM(version)
        WRITE(*,'(A)')  std_lnbrk
        WRITE(*,'(A)')

        CALL DATE_AND_TIME(date, time)

        IF ( (debug .GE. 1_ik) .AND. (my_rank .EQ. 0_ik)) THEN
                WRITE(rd_o,'(2A)')     "This Logfile: ", TRIM(log_file)
                WRITE(rd_o,'(A)')   
                WRITE(rd_o,'(A)')      "HLRS - NUM"
                WRITE(rd_o,'(2A)')     "Image Processing ", TRIM(version)
                WRITE(rd_o,'(4A)')     date, " [ccyymmdd]    ", time, " [hhmmss.sss]"  
                WRITE(rd_o,'(A)')      std_lnbrk
                WRITE(rd_o,'(A, I7)')  "Debug Level:            ", debug
                WRITE(rd_o,'(A, I7)')  "Processors:             ", size_mpi  
                WRITE(rd_o,'(A)')      std_lnbrk
        END IF

ENDIF ! (my_rank==0)

! kernel_spec 0 (/ selectKernel, sizeKernel /) (in the first iteration of this program and for get_cmd_arg)
CALL MPI_BCAST (kernel_spec, 2_mik, MPI_INTEGER         , 0_mik, MPI_COMM_WORLD)
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
        WRITE(n2s,*) kernel_spec(2)
        WRITE(*,'(A)') 'Kernel:'
        DO ii = 1, SIZE(kernel, 1)
                WRITE(*,'(' // TRIM(ADJUSTL(n2s)) //'F7.4)') kernel(ii,:)
        END  DO
        WRITE(*,'(A)')

        ! Import VTK file
        WRITE(*,'(A)') 'VTK File import started'
        CALL read_vtk(fun1, fileName, array, dims, spcng)
        WRITE(*,'(A)') 'VTK File import done'

        ! Get the output Filenames of the Histograms
        histogram_filename_pre__Filter = TRIM(log_file(1:(LEN_TRIM(log_file)-4)))//"_hist_PRE__FILTER.csv"
        histogram_filename_post_Filter = TRIM(log_file(1:(LEN_TRIM(log_file)-4)))//"_hist_POST_FILTER.csv"

ENDIF ! (my_rank==0)

CALL MPI_BCAST (dims       , 3_mik, MPI_INTEGER         , 0_mik, MPI_COMM_WORLD)

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

remainder_per_dir = MODULO(dims, sections)

IF ( (debug .GE. 1_ik) .AND. (my_rank .EQ. 0_ik)) THEN
        WRITE(rd_o,'(A)')      "Calculation of domain sectioning:"
        WRITE(rd_o,'(A)')
        WRITE(rd_o,'(A, 3I5)') "sections:               ", sections
        WRITE(rd_o,'(A, 3I5)') "dims:                   ", dims
        WRITE(rd_o,'(A, 3I5)') "border:                 ", border
        WRITE(rd_o,'(A, 3I5)') "original_image_padding: ", original_image_padding
        WRITE(rd_o,'(A, 3I5)') "remainder_per_dir:      ", remainder_per_dir
END IF

! If remainder per direction is larger than the padding, simply shift the array to distribute to ranks
! into the center of array. Then add the border. This way, no artifical padding is required.
IF ((remainder_per_dir(1) < original_image_padding(1)) .OR. & 
    (remainder_per_dir(2) < original_image_padding(2)) .OR. &
    (remainder_per_dir(3) < original_image_padding(3))) THEN
        dims_reduced   = dims - original_image_padding
ELSE
        dims_reduced   = dims - remainder_per_dir
END IF

! Let's use the remainder to shift the filtered image to the center of the original one.
! This will ensure filtering either the whole image or at least filtering without an artificial 
! Padding at the outer surfaces (borders of the image).

vox_per_dir_and_sec = (dims_reduced / sections) + original_image_padding

IF ( (debug .GE. 1_ik) .AND. (my_rank .EQ. 0_ik)) THEN
        WRITE(rd_o,'(A, 3I5)') "dims-orig_img_padding:  ", clipboard
        WRITE(rd_o,'(A, 3I5)') "new remainder_per_dir:  ", remainder_per_dir
        WRITE(rd_o,'(A, 3I5)') "dims_reduced:           ", dims_reduced
        WRITE(rd_o,'(A, 3I5)') "vox_per_dir_and_sec:    ", vox_per_dir_and_sec
        CLOSE(rd_o)
END IF

if (my_rank==0) THEN
        DO ii = 1, sections(1) 
        DO jj = 1, sections(2) 
        DO kk = 1, sections(3) 
        ! Converting address of subarray into rank is tested in Octave.
        address = (kk-1)*sections(1)*sections(2) + (jj-1_ik)*sections(1) + ii
!            rank_section = sections

        IF ( debug .EQ. 2_ik ) WRITE(*,'(2(A, I5))') "Address: ", address, " My Rank: ",my_rank

!            IF (address .EQ. (my_rank + 1_ik)) THEN
                ! Add original image Data as padding

        subarray_origin = (sections-1_ik) * (vox_per_dir_and_sec - original_image_padding) + border
        
        ! IF ( debug .EQ. 2_ik ) WRITE(*,'(A, I5)') "Pre MPI_Barrier to distribute array across ranks. My Rank:", my_rank
        ! Why the heck is this useful?! Gets a Fatal Error without it. Need someone who has more experience.
        ! CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ! IF ( debug .EQ. 2_ik ) WRITE(*,'(A, I5)') "Post MPI_Barrier to distribute array across ranks. My Rank: ", my_rank

        ! if (my_rank==1) dims = (/ 171, 171 ,171/)
        write(*,*)  "3"                                           ! "My Rank: ", my_rank, 
        write(*,*)  (/ dims(1), dims(2), dims(3) /)               ! "My Rank: ", my_rank,                             
        write(*,*)  vox_per_dir_and_sec                           ! "My Rank: ", my_rank,                 
        write(*,*)  subarray_origin - 1_mik                       ! "My Rank: ", my_rank,                     
        write(*,*)  MPI_ORDER_FORTRAN                             ! "My Rank: ", my_rank,               
        write(*,*)  MPI_INTEGER                                   ! "My Rank: ", my_rank,         
        write(*,*)  type_subarray                                 ! "My Rank: ", my_rank,           
        write(*,*)  ierr                                          ! "My Rank: ", my_rank,  

        ! MPI_TYPE_CREATE_DARRAY may fit better, however dealing with overlaps isn't clear
        CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
        (/ dims(1), dims(2), dims(3) /)     , & ! Original array as all the addresses must fit
        (/ vox_per_dir_and_sec(1),     vox_per_dir_and_sec(2),       vox_per_dir_and_sec(3) /)   , &
        (/ subarray_origin(1) - 1_mik, subarray_origin(2) - 1_mik, subarray_origin(3) - 1_mik /) , & ! array_of_starts indexed from 0
        MPI_ORDER_FORTRAN                   , &
        MPI_INTEGER                         , &
        type_subarray                       , &
        ierr)
        CALL MPI_TYPE_COMMIT(type_subarray, ierr)

        IF (address .NE. 1_ik) THEN
                CALL MPI_SEND(array, 1_mik, type_subarray, address-1_mik, address-1_mik, MPI_COMM_WORLD, ierr )
        END IF

        IF ( debug .EQ. 2_ik ) WRITE(*,'(A, I7, A, I15)') "My Rank: ", address-1_ik, " Size of Subarray:", SIZE(vox_per_dir_and_sec)

        CALL MPI_TYPE_FREE(type_subarray)

!           END IF
        END DO
        END DO
        END DO 
ENDIF

ALLOCATE( subarray(vox_per_dir_and_sec(1), vox_per_dir_and_sec(2), vox_per_dir_and_sec(3) ) )

IF (my_rank > 0) THEN
        
        ! ! Calculate the rank_section out of my_rank and sections(/ x, y, z/)
        ! zremainder = MODULO(my_rank, sections(1)*sections(2))
        ! IF (zremainder .EQ 0_ik) THEN
        !         rank_section = (/ sections(1), sections(2), (my_rank - zremainder) / (sections(1)*sections(2)) /)
        ! ELSE
        !         rank_section(3) = (my_rank - zremainder) / (sections(1) * sections(2)) 
        ! yremainder = MODULO(zremainder, sections(1))
        
        ! IF (yremainder .EQ. 0_ik) THEN
        !         rank_section = (/ sections(1), (zremainder - yremainder) / sections(1), rank_section(3)+1 /)
        ! ELSE
        !         rank_section = (/ yremainder, (zremainder - yremainder) / sections(1) + 1_ik, rank_section(3) + 1_ik /)
        ! ENDIF
        ! ENDIF
        IF ( debug .EQ. 2_ik ) THEN
                WRITE(*,'(A)') "Initialize Recv to distribute array across ranks."
                WRITE(*,'(A, I7, A, 3I5)') "My Rank: ", my_rank, " vox_per_dir_and_sec:    ", vox_per_dir_and_sec
                WRITE(*,'(A, I7, A, I15)') "My Rank: ", my_rank, " Size of Subarray:", SIZE(subarray)
        END IF
        CALL MPI_RECV(subarray, SIZE(subarray), MPI_INTEGER, 0_mik, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
END IF


! CALL MPI_SENDRECV (array, 1_mik, type_subarray, my_rank, my_rank, &
!       subarray, SIZE(subarray), MPI_INTEGER, 0_mik, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

        
! CALL MPI_SENDRECV (array, 1_mik, type_subarray, my_rank, my_rank, &
!         subarray, SIZE(subarray), MPI_INTEGER, 0_mik, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

IF ( debug .EQ. 2_ik ) WRITE(*,'(A)') "Initializing recv to distribute array across ranks done."

! CALL MPI_TYPE_FREE(type_subarray)

IF (my_rank == 0_ik) DEALLOCATE(array)

! Prepare collecting the subarrays to assemble a global vtk file.
! subarray_reduced_boundaries
! No Overlapping and no remainder.

srb (1:3) = border + 1_ik
srb (4:6) = vox_per_dir_and_sec - border - 1_ik

vox_per_dir_and_sec = vox_per_dir_and_sec - original_image_padding
subarray_origin     = vox_per_dir_and_sec * rank_section

ALLOCATE( result_subarray (vox_per_dir_and_sec(1), vox_per_dir_and_sec(2), vox_per_dir_and_sec(3) ) )

! Get information about the data range of the Histogram globally. 
histo_bound_local_lo = MINVAL(subarray)
histo_bound_local_hi = MAXVAL(subarray)

CALL MPI_REDUCE(histo_bound_local_lo, histo_bound_global_lo, 1_mik, MPI_INTEGER, MPI_MIN, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_REDUCE(histo_bound_local_hi, histo_bound_global_hi, 1_mik, MPI_INTEGER, MPI_MIN, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (                      histo_bound_global_lo, 1_mik, MPI_INTEGER,          0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (                      histo_bound_global_hi, 1_mik, MPI_INTEGER,          0_mik, MPI_COMM_WORLD, ierr)

hbnds    = (/ histo_bound_global_lo, histo_bound_global_hi , histo_bound_global_hi - histo_bound_global_lo /)

! Prior to image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (subarray(srb(1):srb(4), srb(2):srb(5), srb(3):srb(6)), hbnds, histogram_pre__F)            

! Start image processing
CALL convolution_un_padded_input (vox_per_dir_and_sec, subarray, kernel_spec(2), kernel, .TRUE., result_subarray)

! After image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (subarray, hbnds, histogram_post_F)            

! MPI_TYPE_CREATE_DARRAY may fit better. Implementation not entirely clear at the moment...
CALL MPI_TYPE_CREATE_SUBARRAY (3_mik, &
        INT(SHAPE(array), KIND=mik) , &
        vox_per_dir_and_sec         , &
        subarray_origin - 1_mik     , & ! array_of_starts indexed from 0
        MPI_ORDER_FORTRAN           , &
        MPI_INTEGER                 , &
        type_result_subarray        , &
        ierr)

CALL MPI_TYPE_COMMIT(type_result_subarray, ierr)

ALLOCATE( result_array(dims_reduced(1), dims_reduced(2), dims_reduced(3) ) )

CALL MPI_SENDRECV (result_subarray, 1_mik, type_result_subarray, 0_mik, my_rank + size_mpi, &
        result_array, 1_mik, MPI_INTEGER, my_rank, my_rank + size_mpi, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
     
CALL MPI_TYPE_FREE(type_result_subarray) 

DEALLOCATE(subarray)

! Collect the data of the histogram pre filtering
CALL MPI_REDUCE (histogram_pre__F, histogram_pre__F_global, 65536_mik, MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

! Collect the data of the histogram post filtering
CALL MPI_REDUCE (histogram_post_F, histogram_post_F_global, 65536_mik, MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

IF (my_rank==0) THEN

        ! Export Histogram of Scalar Array pre Filtering
        OPEN(UNIT = fl_un_H_pre, FILE=histogram_filename_pre__Filter, ACTION="WRITE", STATUS="new")
                WRITE(fl_un_H_pre,'(A)') "Scalar value / 100, Amount of Voxels per Scalar value"
                DO ii=-3276, 3277
                        WRITE(fl_un_H_pre,'(I4,A,I18)') ii," , ",histogram_pre__F_global(ii)
                END DO
        CLOSE(fl_un_H_pre)
        
        ! Export Histogram of Scalar Array post Filtering
        OPEN(UNIT = fl_un_H_post, FILE=histogram_filename_post_Filter, ACTION="WRITE", STATUS="new")
                WRITE(fl_un_H_post,'(A)') "Scalar value / 100, Amount of Voxels per Scalar value"
                DO ii=-3276, 3277       
                     WRITE(fl_un_H_post,'(I4,A,I18)') ii," , ",histogram_post_F_global(ii)
                END DO
        CLOSE(fl_un_H_post)

        ! Export VTK file (testing)
        WRITE(*,'(A)') 'VTK File export started'
        WRITE(n2s,*) kernel_spec(1)
        fileNameExportVtk = fileName(1:LEN_TRIM(fileName)-4) // '_Kernel_'// TRIM(ADJUSTL(n2s))  // '.vtk'
        CALL write_vtk(fun2, fileNameExportVtk, result_array, spcng, dims_reduced)
        WRITE(*,'(A)') 'VTK File export done'

        DEALLOCATE(result_array)
        DEALLOCATE(kernel)

        CALL CPU_TIME(finish)
        WRITE(*,*)
        WRITE(*,'(A7, F8.3, A8)') 'Time = ', (finish - start) / 60,' Minutes'
        WRITE(*,'(A)') std_lnbrk

        IF (debug >= 1) CLOSE(rd_o)

ENDIF ! (my_rank==0)

CALL MPI_FINALIZE(ierr)
CALL MPI_ERR(ierr,"MPI_FINALIZE didn't succeed")

END PROGRAM main