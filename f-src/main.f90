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
        USE standards

        IMPLICIT NONE

        ! For use with MPI (Calls), variables must be declared the same type like MPI was built.
        ! MPI: Kind=32 Bit / 4 Byte / ik=4 - Change in «working_directory/f-src/mod_standards.f90 

        ! Parameter
        INTEGER  (KIND = ik), PARAMETER                                 :: fun1 = 5     ! File unit
        INTEGER  (KIND = ik), PARAMETER                                 :: fun2 = 10    ! File unit

        ! Internal variables
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

        ! MPI Variables
        INTEGER  (KIND=mik)                                             :: ierr
        INTEGER  (KIND=mik)                                             :: my_rank, size_mpi

        ! Debug Variables
        INTEGER  (KIND = ik)                                            :: debug
        INTEGER  (KIND = ik), PARAMETER                                 :: rd_o = 31
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

        !!!!! DEV !!!!!
        ! Currently, everything is handled singlethreaded !
        IF (my_rank==0) THEN
        !!!!! DEV !!!!!

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
        WRITE(*,'(A)') std_lnbrk
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
        WRITE(*,*)

        ! Import VTK file
        WRITE(*,'(A)') 'VTK File import started'
        CALL read_vtk(fun1, fileName, array, dims, spcng)
        WRITE(*,'(A)') 'VTK File import done'

        ! Start image processing
         CALL convolution(dims, array, sizeKernel, kernel)

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