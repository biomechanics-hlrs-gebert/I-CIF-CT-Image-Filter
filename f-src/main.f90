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
CHARACTER(LEN = mcl)                                            :: n2s
CHARACTER(LEN = mcl)                                            :: fileName
CHARACTER(LEN = mcl)                                            :: fileNameExportVtk

INTEGER  (KIND = ik)                                            :: i
INTEGER  (KIND = ik)                                            :: selectKernel    , sizekernel
CHARACTER(LEN = mcl)                                            :: selectKernel_str, sizekernel_str
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: dims

REAL     (KIND = rk)                                            :: start, finish
REAL     (KIND = rk)                                            :: sigma
CHARACTER(LEN = mcl)                                            :: sigma_str, version
REAL     (KIND = rk)           , DIMENSION(3)                   :: spcng
REAL     (KIND = rk)           , DIMENSION(:,:)  , ALLOCATABLE  :: kernel
REAL     (KIND = rk)           , DIMENSION(:,:,:), ALLOCATABLE  :: array

! Histogram Variables
CHARACTER(LEN = mcl)                                            :: histogram_filename_pre__Filter
CHARACTER(LEN = mcl)                                            :: histogram_filename_post_Filter
INTEGER  (KIND = ik), PARAMETER                                 :: fl_un_H_pre=41, fl_un_H_post=42
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: hist_boundaries
! Histogram is scaled to INT2 because 65.535 entries are sufficient to calculate and print virtually any histogram
INTEGER  (KIND = ik)           , DIMENSION(65536), ALLOCATABLE  :: histogram_pre__F       , histogram_post_F
INTEGER  (KIND = ik)           , DIMENSION(65536), ALLOCATABLE  :: histogram_pre__F_global, histogram_post_F_global

! MPI Variables
INTEGER  (KIND=mik)                                             :: ierr
INTEGER  (KIND=mik)                                             :: my_rank, size_mpi, remainder
INTEGER  (KIND=mik)            , DIMENSION(:)    , ALLOCATABLE  :: send_cnt, dsplcmnts
REAL     (KIND=rk)             , DIMENSION(:)    , ALLOCATABLE  :: recv_bffr_1D
INTEGER  (KIND=mik)                                             :: recv_bffr_sz_1D

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
        READ(selectKernel_str,'(I4)') selectKernel
        ! Kernel size
        CALL GET_COMMAND_ARGUMENT(5, sizeKernel_str)
        READ(sizeKernel_str,'(I4)') sizeKernel
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

        ALLOCATE(kernel(sizeKernel, sizeKernel))
        
        SELECT CASE(selectKernel)
                CASE(0)
                        WRITE(*,'(A)') 'Identity kernel selected'
                        CALL kernel_identity(kernel, sizeKernel)
                CASE(1)
                        WRITE(*,'(A)') 'Gaussian filter kernel selected'
                        CALL kernel_gauss(kernel, sizeKernel, sigma)
                CASE DEFAULT
                        WRITE(*,'(A)') 'Identity kernel selected'
                        CALL kernel_identity(kernel, sizeKernel)
        END SELECT
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

! Calculate and Allocate information for scattering the global array.
ALLOCATE(send_cnt(size_mpi))

remainder = MODULO(entries, size_mpi)     ! remainder set to "last rank"
recv_bffr_sz_1D = (entries - remainder) / size_mpi
send_cnt(1:size_mpi ) = recv_bffr_sz_1D
send_cnt(1:remainder) = send_cnt(1:remainder)+1_mik

ALLOCATE(recv_bffr_1D(send_cnt(my_rank+1_mik)))
ALLOCATE(dsplcmnts(size_mpi))

dsplcmnts(:)=0_mik

! Distribute the array globally
CALL MPI_SCATTERV(array, send_cnt, dsplcmnts, MPI_DOUBLE_PRECISION, recv_bffr_1D,  &
                     send_cnt(my_rank+1_mik), MPI_DOUBLE_PRECISION ,0_mik, MPI_COMM_WORLD)

! Prior to image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (array, hist_boundaries, dims, histogram_pre__F)            

! Start image processing
CALL convolution (dims, array, sizeKernel, kernel)

! After image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (array, hist_boundaries, dims, histogram_post_F)

! Collect the subarrays to assemble a global vtk file
CALL MPI_GATHERV(array, send_cnt, dsplcmnts, MPI_DOUBLE_PRECISION, recv_bffr_1D,  &
        send_cnt(my_rank+1_mik), MPI_DOUBLE_PRECISION ,0_mik, MPI_COMM_WORLD)

! Collect the data of the histogram pre filtering
CALL MPI_REDUCE (histogram_pre__F, histogram_global_pre__F, &
        65536_mik, MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

! Collect the data of the histogram post filtering
CALL MPI_REDUCE (histogram_pre__F, histogram_global_post__F, &
        65536_mik, MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

DEALLOCATE(recv_bffr_1D)
DEALLOCATE(send_cnt)
DEALLOCATE(dsplcmnts)

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
        CALL write_vtk(fun2, fileNameExportVtk, array, spcng, dims)
        WRITE(*,'(A)') 'VTK File export done'

        DEALLOCATE(kernel)

        CALL CPU_TIME(finish)
        WRITE(*,*)
        WRITE(*,'(A7, F8.3, A8)') 'Time = ', (finish - start) / 60,' Minutes'
        WRITE(*,'(A)') std_lnbrk

ENDIF ! (my_rank==0)

CALL MPI_FINALIZE(ierr)
CALL MPI_ERR(ierr,"MPI_FINALIZE didn't succeed")

END PROGRAM main