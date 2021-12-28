!------------------------------------------------------------------------------
! PROGRAM: Computed Tomography Image Filter
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!------------------------------------------------------------------------------
PROGRAM CTIF

USE global_std
USE messages_errors
USE meta
USE MPI
USE file_routines_mpi  
USE kernels
USE aux_routines_ip


  USE meta
  USE meta_puredat_interface
  USE chain_routines
  USE MPI
  USE decomp 
  USE sp_aux_routines
  USE PETSC
  USE petsc_opt

IMPLICIT NONE 

! Parameter
INTEGER(KIND=ik), PARAMETER :: mov_avg_width = 100   ! Choose an even integer!!

! Internal Variables
CHARACTER(LEN=mcl) :: filename, filenameExportVtk, typ
CHARACTER(LEN=mcl) :: selectKernel
REAL     (KIND=rk) :: global_start, init_finish, read_t_vtk, prep_Histo
REAL     (KIND=rk) :: calculation, extract_Histo, global_finish
INTEGER  (KIND=ik) :: ii, jj, kk, ll, mm, nn
INTEGER  (KIND=ik) :: displacement, counter, border, frcs=0
INTEGER  (KIND=ik), DIMENSION(2) :: kernel_spec

INTEGER  (KIND=ik), DIMENSION(3) :: dims, original_image_padding, subarray_origin
INTEGER  (KIND=ik), DIMENSION(3) :: sections, rank_section, subarray_dims
INTEGER  (KIND=ik), DIMENSION(3) :: dims_reduced, remainder_per_dir, subarray_dims_overlap

REAL     (KIND=rk) :: sigma, accumulator
CHARACTER(LEN=mcl) :: version, basename, inp_para
CHARACTER(LEN=  8) :: date
CHARACTER(LEN= 10) :: time
INTEGER  (KIND=ik), DIMENSION(6) :: srb ! subarray_reduced_bndaries
REAL     (KIND=rk), DIMENSION(3) :: spcng, origin
REAL     (KIND=rk), DIMENSION(:,:)  , ALLOCATABLE  :: kernel2d
REAL     (KIND=rk), DIMENSION(:,:,:), ALLOCATABLE  :: kernel3d
INTEGER  (KIND=ik), DIMENSION(:,:,:), ALLOCATABLE  :: subarray, result_subarray     ! Dealt with internally as int32

! Histogram Variables
CHARACTER(LEN=mcl) :: histogram_fn_pre__Filter
CHARACTER(LEN=mcl) :: histogram_fn_post_Filter
CHARACTER(LEN=mcl) :: histo_avg_fn_pre__Filter
CHARACTER(LEN=mcl) :: histo_avg_fn_post_Filter
CHARACTER(LEN=mcl) :: histogram_filename_tex_Filter
INTEGER  (KIND=ik) :: histo_bnd_global_lo, histo_bnd_global_hi
INTEGER  (KIND=ik) :: histo_bnd_local_lo,  histo_bnd_local_hi
INTEGER  (KIND=ik), DIMENSION(3) :: hbnds
INTEGER  (KIND=ik), DIMENSION(:), ALLOCATABLE  :: histogram_pre__F       , histogram_post_F
INTEGER  (KIND=ik), DIMENSION(:), ALLOCATABLE  :: histogram_pre__F_global, histogram_post_F_global
INTEGER(KIND=ik) :: fh_csv_prf  = 200, fh_csv_pof  = 210, fh_csv_aprf = 220, fh_csv_apof = 230, fh_csv_fihi = 240

! Read Input file
CHARACTER(len=mcl) :: line, parameterfile, prefix
INTEGER  (KIND=ik) :: io_status, ntokens
CHARACTER(len=mcl) :: tokens(100)
CHARACTER(len=mcl) :: tkns(100)

! MPI Variables
INTEGER  (KIND = mik) :: ierr, my_rank, size_mpi, status
INTEGER  (KIND = mik) :: wr_vtk_hdr_lngth


! Initialize MPI Environment
CALL MPI_INIT(ierr)
CALL MPI_ERR(ierr,"MPI_INIT didn't succeed")

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_RANK couldn't be retrieved")

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_SIZE couldn't be retrieved")

IF (size_mpi < 2) THEN
    WRITE(*,*) "my_rank: ", my_rank, " size_mpi: ", size_mpi
    WRITE(std_out,*)"At least two ranks required to execute this program."
    CLOSE(std_out)
    CALL MPI_ABORT(MPI_COMM_WORLD, 1_mik, ierr)
END IF

! Initialize program itself
IF (my_rank == 0) THEN

    ! Initialization
    CALL CPU_TIME(global_start)

        !------------------------------------------------------------------------------
        ! Parse the command arguments
        !------------------------------------------------------------------------------
        IF (command_argument_count() == 0) THEN 
            CALL usage()
            GOTO 1001
        END IF

        DO ii=0, 15 ! Read up to 15 command arguments.
         
            CALL GET_COMMAND_ARGUMENT(ii, cmd_arg)

            IF (cmd_arg == '') EXIT

            IF (ii == 0) binary = cmd_arg

            infile = TRIM(cmd_arg)
            
            cmd_arg_history = TRIM(cmd_arg_history)//' '//TRIM(cmd_arg)

            IF (cmd_arg(1:1) .EQ. '-') THEN
                SELECT CASE( cmd_arg(2:LEN_TRIM(cmd_arg)) )
                CASE('-restart', '-Restart')
                    restart_cmdarg = 'Y'
                CASE('v', '-Version', '-version')
                        CALL show_title(revision)
                        GOTO 1001   ! Jump to the end of the program.
                CASE('h', '-Help', '-help')
                    CALL usage()
                    GOTO 1001   ! Jump to the end of the program.
                END SELECT
                !
                SELECT CASE( cmd_arg(3:4) )
                CASE('NO', 'no', 'No', 'nO');  restart_cmdarg = 'N'
                END SELECT
            END IF
        END DO

        IF(TRIM(infile) == '') THEN
            CALL print_err_stop(std_out, 'No input file given via command argument.', 1)
        END IF 

        in%full = TRIM(infile)

        !------------------------------------------------------------------------------
        ! Check and open the input file; Modify the Meta-Filename / Basename
        ! Define the new application name first
        !------------------------------------------------------------------------------
        global_meta_prgrm_mstr_app = 'CTIF' 
        global_meta_program_keyword = 'CT-IMAGE_FILTER'
        CALL meta_append(m_rry)
      
        !------------------------------------------------------------------------------
        ! Parse input
        !------------------------------------------------------------------------------
        CALL meta_read (std_out, 'DEBUG_LVL'    , m_rry, debug)
        CALL meta_read (std_out, 'FILTER_DIM'   , m_rry, kernel_spec(1))
        CALL meta_read (std_out, 'FILTER_SIZE'  , m_rry, kernel_spec(2))
        CALL meta_read (std_out, 'FILTER_KERNEL', m_rry, selectKernel )
        CALL meta_read (std_out, 'FILTER_SIGMA' , m_rry, sigma)

        !------------------------------------------------------------------------------
        ! Restart handling
        ! Done after meta_io to decide based on keywords
        !------------------------------------------------------------------------------
        IF (restart_cmdarg /= 'U') THEN
            mssg = "The keyword »restart« was overwritten by the command flag --"
            IF (restart_cmdarg == 'N') THEN
                restart = restart_cmdarg
                mssg=TRIM(mssg)//"no-"
            ELSE IF (restart_cmdarg == 'Y') THEN
                restart = restart_cmdarg
            END IF

            mssg=TRIM(mssg)//"restart"
            WRITE(std_out, FMT_WRN) TRIM(mssg)
            WRITE(std_out, FMT_WRN_SEP)
        END IF

        CALL meta_handle_lock_file(restart)

        !------------------------------------------------------------------------------
        ! Spawn a log file and the csv-data for the tex environment
        !------------------------------------------------------------------------------
        CALL meta_start_ascii(fhl, log_suf)

        CALL meta_start_ascii(fh_csv_prf, "_hist_PRE__FILTER"//csv_suf)
        CALL meta_start_ascii(fh_csv_pof, "_hist_POST_FILTER"//csv_suf)
        CALL meta_start_ascii(fh_csv_aprf, "_hist_avg_PRE__FILTER"//csv_suf)
        CALL meta_start_ascii(fh_csv_apof, "_hist_avg_POST_FILTER"//csv_suf)
        CALL meta_start_ascii(fh_csv_fihi, "_Filter_Histogram"//csv_suf)

        CALL meta_start_ascii(fht, tex_suf)
        
        IF (std_out/=6) THEN
            CALL meta_start_ascii(std_out, '.std_out')

            CALL show_title(revision) 
        END IF


        CALL DATE_AND_TIME(date, time)

        IF ( (debug .GE. 0_ik) .AND. (my_rank == 0_ik)) THEN
            WRITE(std_out, FMT_TXT) "This Logfile:"//TRIM(log_file)
            WRITE(std_out, FMT_TXT) 
            WRITE(std_out, FMT_TXT) "HLRS - NUM"
            WRITE(std_out, FMT_TXT) "CT Image Filtering"//TRIM(version)
            WRITE(std_out, FMT_TXT) date//" [ccyymmdd]"//time//" [hhmmss.sss]"  
            WRITE(std_out, FMT_TXT_SEP)
            WRITE(std_out, FMT_MSG_AI0) "Debug Level:", debug
            WRITE(std_out, FMT_MSG_AI0) "Processors:", size_mpi  
            WRITE(std_out, FMT_MSG_AI0) "Filter Dimension:", kernel_spec(1)
            WRITE(std_out, FMT_MSG_AI0) "Filter Size:", kernel_spec(2)
            WRITE(std_out, FMT_MSG)     "Filter Kernel:"//TRIM(selectKernel)
            WRITE(std_out, FMT_TXT_SEP)
            FLUSH(std_out)
        END IF
        
        CALL CPU_TIME(init_finish)

!!!!!!!!!!!!!!!!!!!!!!!!! GO ON !!!!!!!!! VTK TO RAW AT FIRST? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
ENDIF ! (my_rank == 0)

! kernel_spec 0 (/ selectKernel, sizeKernel /) (in the first iteration of this program and for get_cmd_arg)
! May be packed into less BCasts
CALL MPI_BCAST (filename    , INT(mcl, KIND=mik), MPI_CHAR            , 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (typ         , INT(mcl, KIND=mik), MPI_CHAR            , 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (kernel_spec , 2_mik             , MPI_INTEGER         , 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (sigma       , 1_mik             , MPI_DOUBLE_PRECISION, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (selectKernel, INT(mcl, KIND=mik), MPI_CHAR            , 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (dims        , 3_mik             , MPI_INTEGER         , 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (displacement, 1_mik             , MPI_INTEGER         , 0_mik, MPI_COMM_WORLD, ierr)

! Get sections per direction
CALL MPI_DIMS_CREATE (size_mpi, 3_mik, sections, ierr)
CALL get_rank_section (domain=my_rank, sections=sections, rank_section=rank_section)

! Calculate Padding to decrease "size of array" to a corresponding size
! Size if Kernel always an odd number. Num-1 = Voxels on both sides of filtered Voxel
original_image_padding =  kernel_spec(2)-1_ik         ! 1D Array
border         = (kernel_spec(2)-1_ik) / 2_ik     ! 0D Array (scalar)

! Remainder per direction gets almost fully ignored (!) 
! It's assumed, that even if we split into 32768 processes (/ 32, 32, 32 /) sections,
! max. 31 x 31 x 31 Voxel get lost (Worst Case). Considering large input sets,
! which are imperative for utilizing large amounts of processors, losing 31 voxel at
! ~ 2000 ... 4000 Voxel per direction is not an issue.
! On the other hand, Distribution of the array gets way easier and presumably quicker.

remainder_per_dir = MODULO(dims, sections)

! If remainder per direction is larger than the padding, simply shift the array to distribute to ranks
! into the center of array. Then add the border. This way, no artifical padding is required.
IF ((remainder_per_dir(1) <= original_image_padding(1)) .OR. & 
    (remainder_per_dir(2) <= original_image_padding(2)) .OR. &
    (remainder_per_dir(3) <= original_image_padding(3))) THEN
    dims_reduced   = dims - original_image_padding
ELSE
    dims_reduced   = dims - remainder_per_dir
END IF

subarray_dims     = (dims_reduced / sections)
subarray_dims_overlap = subarray_dims + original_image_padding

IF ( (debug .GE. 1_ik) .AND. (my_rank == 0_ik) ) THEN
    WRITE(std_out,'(A)')      std_lnbrk
    WRITE(std_out,'(A)')      "Calculation of domain sectioning:"
    WRITE(std_out,'(A)')
    WRITE(std_out,'(A, 3I5)') "sections:           ", sections
    WRITE(std_out,'(A, 3I5)') "dims:           ", dims
    WRITE(std_out,'(A, 3I5)') "border:         ", border
    WRITE(std_out,'(A, 3I5)') "original_image_padding: ", original_image_padding
    WRITE(std_out,'(A, 3I5)') "remainder_per_dir:      ", remainder_per_dir
    WRITE(std_out,'(A, 3I5)') "dims_reduced:       ", dims_reduced
    WRITE(std_out,'(A, 3I5)') "subarray_dims:      ", subarray_dims
    WRITE(std_out,'(A)')      std_lnbrk
    FLUSH(std_out)
END IF

ALLOCATE( subarray(subarray_dims_overlap(1), subarray_dims_overlap(2), subarray_dims_overlap(3) ) )

subarray_origin = (rank_section-1_ik) * (subarray_dims) !+ 1_ik

CALL read_raw_mpi(filename=filename, type_in=TRIM(typ), type_out=typ, hdr_lngth=INT(displacement, KIND=8), &
            dims=dims, subarray_dims=subarray_dims_overlap, subarray_origin=subarray_origin, subarray=subarray, &
            displacement=displacement, log_un=std_out, status_o=status)

IF ((status == 1_ik) .AND. (my_rank == 1_ik)) THEN
    WRITE(std_out,'(A)')  'Something during MPI File read went wrong. Please check/debug.'
    CLOSE(std_out)
    CALL MPI_ABORT (MPI_COMM_WORLD, 1_mik, ierr)
END IF

! Prepare collecting the subarrays to assemble a global vtk file.
IF (my_rank == 0_ik) CALL CPU_TIME(read_t_vtk)

! subarray_reduced_bndaries           ! No Overlap
srb (1:3) = 1_ik + border
srb (4:6) = subarray_dims_overlap - border

! Define new subarray_origin
! subarray_origin = subarray_dims * (rank_section-1_ik) + 1_ik 

ALLOCATE( result_subarray (subarray_dims(1), subarray_dims(2), subarray_dims(3) ) )

! Get information about the data range of the Histogram globally. 
histo_bnd_local_lo = MINVAL(subarray)
histo_bnd_local_hi = MAXVAL(subarray)

CALL MPI_ALLREDUCE(histo_bnd_local_lo, histo_bnd_global_lo, 1_mik, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
CALL MPI_ALLREDUCE(histo_bnd_local_hi, histo_bnd_global_hi, 1_mik, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

hbnds    = (/ histo_bnd_global_lo, histo_bnd_global_hi , histo_bnd_global_hi - histo_bnd_global_lo /)

! Prior to image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (subarray(srb(1):srb(4), & 
                            srb(2):srb(5), srb(3):srb(6)), hbnds, histogram_pre__F)        

IF (my_rank == 0_ik) CALL CPU_TIME(prep_Histo)

! Start image processing
! result_image is necessary, because otherwise, filtered Voxels will be used for filtering following voxels.
! Therefore, doesn't really work in place
! Ensure image padding manually in front of this branch
IF (kernel_spec(1) == 2_ik) THEN
    ! 2D must be ordered explicitly.
    ALLOCATE( kernel2d ( kernel_spec(2), kernel_spec(2) ))

    SELECT CASE(selectKernel)
        CASE("Gaussian"); CALL kernel_gauss_2d   (kernel2d, kernel_spec(2), sigma)
        CASE DEFAULT;     CALL kernel_identity_2d(kernel2d, kernel_spec(2))
    END SELECT

    DO ii = srb(1), srb(4)
    DO jj = srb(2), srb(5)
    DO kk = srb(3), srb(6)
        accumulator = 0
        DO ll = -border, border
        DO mm = -border, border
            accumulator = accumulator + (kernel2d( ll+border+1_ik, mm+border+1_ik ) * &
                    subarray(ii + ll, jj + mm, ii))
        END DO
        END DO
        result_subarray(ii - border, jj - border, kk - border) = accumulator
    END DO
    END DO
    END DO

    DEALLOCATE(kernel2d)
ELSE    
    ! 3D is considered a default
    ALLOCATE( kernel3d(kernel_spec(2), kernel_spec(2), kernel_spec(2)))

    SELECT CASE(selectKernel)
        CASE("Gaussian"); CALL kernel_gauss_3d   (kernel3d, kernel_spec(2), sigma)
        CASE DEFAULT;     CALL kernel_identity_3d(kernel3d, kernel_spec(2))
    END SELECT
    
    DO ii = srb(1), srb(4)
    DO jj = srb(2), srb(5)
    DO kk = srb(3), srb(6)
        accumulator = 0
        DO ll = -border, border
        DO mm = -border, border
        DO nn = -border, border
            accumulator = accumulator + (kernel3d( ll+border+1_ik, mm+border+1_ik, nn+border+1_ik) * &
                    subarray(ii + ll, jj + mm, kk + nn))
        END DO
        END DO
        END DO
        result_subarray(ii - border, jj - border, kk - border) = accumulator
    END DO
    END DO
    END DO

    DEALLOCATE(kernel3d)
ENDIF

DEALLOCATE(subarray)

IF (my_rank == 0_ik) CALL CPU_TIME(calculation)

! After image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (result_subarray, hbnds, histogram_post_F)        

! Allocate memory for global histogram
IF (my_rank == 0_ik) THEN
    ALLOCATE(histogram_pre__F_global(hbnds(1):hbnds(2)))
    ALLOCATE(histogram_post_F_global(hbnds(1):hbnds(2)))
END IF

! Collect the data of the histogram post filtering       
CALL MPI_REDUCE (histogram_pre__F, histogram_pre__F_global, INT(SIZE(histogram_pre__F), KIND=mik), &
    MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_REDUCE (histogram_post_F, histogram_post_F_global, INT(SIZE(histogram_post_F), KIND=mik), &
    MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

CALL CPU_TIME(extract_Histo)

IF (my_rank == 0_ik) THEN

    ! Export Histograms
    CALL write_histo_csv (fh_pre,      histogram_fn_pre__Filter, hbnds, 0_ik       , histogram_pre__F_global)
    CALL write_histo_csv (fh_post,     histogram_fn_post_Filter, hbnds, 0_ik       , histogram_post_F_global)
    CALL write_histo_csv (fh_pre__avg, histo_avg_fn_pre__Filter, hbnds, mov_avg_width, histogram_pre__F_global)
    CALL write_histo_csv (fh_post_avg, histo_avg_fn_post_Filter, hbnds, mov_avg_width, histogram_post_F_global)

    CALL write_tex_for_histogram (fun3, &
        histogram_filename_tex_Filter, &
        histogram_fn_pre__Filter, &
        histogram_fn_post_Filter, &
        histo_avg_fn_pre__Filter, &
        histo_avg_fn_post_Filter)

    filenameExportVtk = TRIM(basename)//'_'//TRIM(inp_para)//'.vtk'

    CALL write_vtk_meta (fh=fh_data_out, &
                filename=filenameExportVtk, & 
                type=TRIM(typ), &
                atStart=.TRUE., &
                spcng=spcng, &
                origin=origin, &
                dims=sections*subarray_dims)

    INQUIRE(FILE=filenameExportVtk, SIZE=wr_vtk_hdr_lngth)

END IF ! (my_rank == 0_ik)

! BCAST used in some way of a Barrier.
CALL MPI_BCAST (filenameExportVtk, INT(mcl, KIND=mik), MPI_CHAR   , 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST (wr_vtk_hdr_lngth , 1_mik         , MPI_INTEGER, 0_mik, MPI_COMM_WORLD, ierr)

IF (TRIM(typ) == "int2") THEN
CALL write_raw_mpi (type=TRIM(typ), &
            hdr_lngth=INT(wr_vtk_hdr_lngth, KIND=8), &
            filename=filenameExportVtk, &
            dims=sections*subarray_dims, &
            subarray_dims=subarray_dims, &
            subarray_origin=subarray_origin, &
            subarray2=INT(result_subarray, KIND=INT16))
ELSE
CALL write_raw_mpi (type=TRIM(typ), &
            hdr_lngth=INT(wr_vtk_hdr_lngth, KIND=8) , &
            filename=filenameExportVtk, &
            dims=sections*subarray_dims, &
            subarray_dims=subarray_dims, &
            subarray_origin=subarray_origin, &
            subarray4=result_subarray)
END IF

DEALLOCATE(result_subarray)

IF (my_rank == 0_ik) THEN

    CALL write_vtk_meta (fh=fh_data_out, filename=filenameExportVtk, atStart=.FALSE.)

    CALL CPU_TIME(global_finish)

        WRITE(std_out,'(A         )')  std_lnbrk
        WRITE(std_out,'(A, F13.3, A)') 'Init and parsing     = ', (init_finish        - global_start),' Seconds'
        WRITE(std_out,'(A, F13.3, A)') 'Read File            = ', (read_t_vtk         - init_finish) ,' Seconds'
        WRITE(std_out,'(A, F13.3, A)') 'Prep Histograms      = ', (prep_Histo - read_t_vtk)          ,' Seconds'
        WRITE(std_out,'(A, F13.3, A)') 'Calculation          = ', (calculation        - prep_Histo)  ,' Seconds'
        WRITE(std_out,'(A, F13.3, A)') 'Extract Histograms   = ', (extract_Histo - calculation)      ,' Seconds'
        WRITE(std_out,'(A, F13.3, A)') 'Calculate Histograms = ', (prep_Histo - read_t_vtk + extract_Histo - calculation) ,' Seconds'
        WRITE(std_out,'(A, F13.3, A)') 'Write all data       = ', (global_finish      - extract_Histo)                    ,' Seconds'
        WRITE(std_out,'(A, F13.3, A)') 'Overall Time         = ', (global_finish      - global_start)                     ,' Seconds'
        WRITE(std_out,'(A, F13.3, A)') 'Overall Time         = ', (global_finish      - global_start) / 60                ,' Minutes'
        WRITE(std_out,'(A         )')  std_lnbrk
        WRITE(std_out,'(A, F13.3, A)') 'CPU time             = ', (global_finish      - global_start) / 60 / 60 * size_mpi,' Hours'

    CLOSE(std_out)     
ENDIF 

1001 Continue

IF(rank_mpi == 0) THEN
    CALL meta_signing(binary)
    CALL meta_close()

    CALL meta_stop_ascii(fhl, log_suf)

    CALL meta_stop_ascii(fh_csv_prf, "_hist_PRE__FILTER"//csv_suf)
    CALL meta_stop_ascii(fh_csv_pof, "_hist_POST_FILTER"//csv_suf)
    CALL meta_stop_ascii(fh_csv_aprf, "_hist_avg_PRE__FILTER"//csv_suf)
    CALL meta_stop_ascii(fh_csv_apof, "_hist_avg_POST_FILTER"//csv_suf)
    CALL meta_stop_ascii(fh_csv_fihi, "_Filter_Histogram"//csv_suf)

    CALL meta_stop_ascii(fht, tex_suf)

    IF (std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

END IF ! (rank_mpi == 0)

Call MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE didn't succeed", INT(ierr, KIND=ik))

END PROGRAM CTIF
