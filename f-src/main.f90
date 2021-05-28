PROGRAM main
!-------------------------------
! Image processing
!
! Author:  Benjamin Schnabel, M.Sc.
! Author:  Johannes Gebert - HLRS - NUM
! GitHub:  https://github.com/JoGebert/3D_Convolusional_Filtering
! Date:    23.04.2021
! LastMod: 25.05.2021
!-------------------------------


USE ISO_FORTRAN_ENV
USE MPI
USE file_routines_mpi  
USE kernels
USE strings
USE aux_routines_ip
USE standards

IMPLICIT NONE

! For use with MPI (Calls), variables must be declared the same type like MPI was built.
! MPI: Kind=32 Bit / 4 Byte / ik=4 - Change in «working_directory/f-src/mod_standards.f90 

! Parameter
INTEGER  (KIND = ik), PARAMETER                                 :: debug       = 1
INTEGER  (KIND = ik), PARAMETER                                 :: mov_avg_width = 100   ! Choose an even integer!!
INTEGER  (KIND = ik), PARAMETER                                 :: fun_input   = 999
INTEGER  (KIND = ik), PARAMETER                                 :: fh_data_in  = 5  
INTEGER  (KIND = ik), PARAMETER                                 :: fh_data_out = 10   ! write vtk
INTEGER  (KIND = ik), PARAMETER                                 :: fun3        = 15   ! write tex

! Internal Variables
REAL     (KIND = rk)                                            :: global_start, init_finish, read_t_vtk, prep_Histo
REAL     (KIND = rk)                                            :: calculation, extract_Histo, global_finish
CHARACTER(LEN = mcl)                                            :: filename, filenameExportVtk, typ
INTEGER  (KIND = ik)                                            :: ii, jj, kk, ll, mm, nn, displacement, counter, border
INTEGER  (KIND = ik)           , DIMENSION(2)                   :: kernel_spec
CHARACTER(LEN = mcl)                                            :: selectKernel

INTEGER  (KIND = ik)           , DIMENSION(3)                   :: dims, original_image_padding, subarray_origin
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: sections, rank_section, subarray_dims
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: dims_reduced, remainder_per_dir, subarray_dims_overlap

INTEGER  (KIND = ik)           , DIMENSION(6)                   :: srb ! subarray_reduced_bndaries
REAL     (KIND = rk)                                            :: sigma, accumulator
CHARACTER(LEN = mcl)                                            :: version, basename
REAL     (KIND = rk)           , DIMENSION(3)                   :: spcng
REAL     (KIND = rk)           , DIMENSION(:,:)  , ALLOCATABLE  :: kernel2d
REAL     (KIND = rk)           , DIMENSION(:,:,:), ALLOCATABLE  :: kernel3d
INTEGER  (KIND = ik)           , DIMENSION(:,:,:), ALLOCATABLE  :: subarray, result_subarray     ! Dealt with internally as int32
CHARACTER(LEN =   8)                                            :: date
CHARACTER(LEN =  10)                                            :: time
CHARACTER(LEN = mcl)                                            :: inp_para

! Histogram Variables
CHARACTER(LEN = mcl)                                            :: histogram_fn_pre__Filter
CHARACTER(LEN = mcl)                                            :: histogram_fn_post_Filter
CHARACTER(LEN = mcl)                                            :: histo_avg_fn_pre__Filter
CHARACTER(LEN = mcl)                                            :: histo_avg_fn_post_Filter
CHARACTER(LEN = mcl)                                            :: histogram_filename_tex_Filter
INTEGER  (KIND = ik)                                            :: histo_bnd_global_lo, histo_bnd_global_hi
INTEGER  (KIND = ik)                                            :: histo_bnd_local_lo,  histo_bnd_local_hi
INTEGER  (KIND = ik), PARAMETER                                 :: fh_pre     =41, fh_post    =42
INTEGER  (KIND = ik), PARAMETER                                 :: fh_pre__avg=51, fh_post_avg=52
INTEGER  (KIND = ik)           , DIMENSION(3)                   :: hbnds
INTEGER  (KIND = ik)           , DIMENSION(:)    , ALLOCATABLE  :: histogram_pre__F       , histogram_post_F
INTEGER  (KIND = ik)           , DIMENSION(:)    , ALLOCATABLE  :: histogram_pre__F_global, histogram_post_F_global

! Read Input file
CHARACTER(len=mcl)                                              :: line, parameterfile, prefix
INTEGER  (KIND=ik)                                              :: io_status, ntokens
CHARACTER(len=mcl)                                              :: tokens(100)
CHARACTER(len=mcl)                                              :: tkns(100)

! MPI Variables
INTEGER  (KIND = mik)                                           :: ierr, my_rank, size_mpi, status
INTEGER  (KIND = mik)                                           :: wr_vtk_hdr_lngth

! Debug Variables
INTEGER  (KIND = ik), PARAMETER                                 :: rd_o = 31    ! redirected StdOut
CHARACTER(LEN = mcl)                                            :: log_file
LOGICAL                                                         :: log_exist = .FALSE.

! Initialize MPI Environment
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

! Initialize program itself
IF (my_rank .EQ. 0) THEN

        ! Initialization
        CALL CPU_TIME(global_start)

        CALL GET_COMMAND_ARGUMENT(1, prefix)
        CALL GET_COMMAND_ARGUMENT(2, parameterfile)

        ! Check validity of the input file
        IF ( parameterfile(1:1) .EQ. '.'  ) parameterfile = parameterfile ( 2:LEN_TRIM(parameterfile) )

        parameterfile = TRIM(prefix)//'/'//TRIM(parameterfile)

        CALL check_file_exist( filename = parameterfile, must_exist=1_ik, mpi=.TRUE.)
               
        ! Open the input file
        OPEN(UNIT=fun_input, FILE=TRIM(parameterfile), STATUS="OLD")

        ! Count lines of input file
        counter = 0_ik
        DO
                READ(fun_input, '(A)', iostat=io_status) line
                if ( io_status .NE. 0_ik ) exit
                counter = counter + 1_ik
        END DO

        CLOSE(fun_input)
        OPEN(UNIT=fun_input, FILE=TRIM(parameterfile), STATUS="OLD")

        ! Parse input file
        DO ii=1, counter
                READ(fun_input,'(A)') line
                
                CALL parse(str=line,delims=" ",args=tokens,nargs=ntokens)
                
                IF (ntokens > 0) THEN
        
                        IF ( tokens(1) .NE. "#" ) THEN

                                CALL parse(str=tokens(2), delims="=", args=tkns, nargs=ntokens)

                                SELECT CASE( tkns(1) )
                                        CASE("IP_DATA_IN");     filename = TRIM(prefix)//tkns(2)
                                        CASE("IP_MODE_K");      READ(tkns(2),'(I4)') kernel_spec(1)
                                        CASE("IP_SELECT_K");    selectKernel = tkns(2)
                                        CASE("IP_SIZE_K");      READ(tkns(2),'(I4)') kernel_spec(2)
                                        CASE("IP_GS");          READ(tkns(2),'(F8.3)') sigma
                                END SELECT
                        END IF
                END IF
         END DO
       
        CLOSE(fun_input)

        ! Get basename of vtk file
        CALL parse( str=filename, delims="/", args=tokens, nargs=ntokens)

        basename = tokens(ntokens)
        basename = TRIM(filename(1:(LEN_TRIM(filename) - 4_ik )))

        CALL parse(str=parameterfile,delims=".",args=tokens,nargs=ntokens)
        inp_para = tokens(2)

        ! Log in dir of vtk - not its basename!!
        log_file  = TRIM(basename)//'_'//TRIM(inp_para)//".log"
        CALL check_file_exist( filename = log_file, must_exist=0_ik, mpi=.TRUE.)

        INQUIRE(FILE=TRIM(log_file), EXIST=log_exist)

        OPEN(UNIT = rd_o, file = TRIM(log_file), action="WRITE", status="new")

        CALL DATE_AND_TIME(date, time)

        IF ( (debug .GE. 0_ik) .AND. (my_rank .EQ. 0_ik)) THEN
                WRITE(rd_o,'(2A)')     "This Logfile: ", TRIM(log_file)
                WRITE(rd_o,'(A)')   
                WRITE(rd_o,'(A)')      "HLRS - NUM"
                WRITE(rd_o,'(2A)')     "Image Processing ", TRIM(version)
                WRITE(rd_o,'(4A)')     date, " [ccyymmdd]    ", time, " [hhmmss.sss]"  
                WRITE(rd_o,'(A)')      std_lnbrk
                WRITE(rd_o,'(A, I7)')  "Debug Level:            ", debug
                WRITE(rd_o,'(A, I7)')  "Processors:             ", size_mpi  
                WRITE(rd_o,'(A, I7)')  "Filter Dimension:       ", kernel_spec(1)
                WRITE(rd_o,'(A, I7)')  "Filter Size:            ", kernel_spec(2)
                WRITE(rd_o,'(2A)')     "Filter Kernel:          ", TRIM(selectKernel)
                WRITE(rd_o,'(A)')      std_lnbrk
        END IF
        
        CALL CPU_TIME(init_finish)

        IF( filename(LEN_TRIM(filename)-2 : LEN_TRIM(filename)) .NE. "vtk" )  THEN
                WRITE(rd_o,'(A)') 'Input file to read VTK subroutine has no vtk suffix. Therefor pressumably is none.'  
                CLOSE(rd_o)  
                CALL MPI_ABORT (MPI_COMM_WORLD, 1_mik, ierr)
        ENDIF

        CALL check_file_exist( filename=filename, must_exist=1_ik, mpi=.TRUE.)

        ! Read VTK file header
        CALL read_vtk_meta (    fh=fh_data_in                   , &
                                filename=filename               , &
                                dims=dims                       , &
                                spcng=spcng                     , &
                                typ=typ                         , &
                                displacement=displacement       , &
                                rd_o=rd_o                       , &
                                status_o=status)

        IF (status .EQ. 1_ik) THEN
                WRITE(rd_o,'(A)')  'Something went wrong while reading IP_DATA_IN header. Program Aborted.'
                CLOSE(rd_o)     
                CALL MPI_ABORT(MPI_COMM_WORLD, 1_mik, ierr)      
        ENDIF
       
        ! Get the output filenames of the Histograms
        histogram_fn_pre__Filter = TRIM(basename)//'_'//TRIM(inp_para)//'_hist_PRE__FILTER.csv'
        histogram_fn_post_Filter = TRIM(basename)//'_'//TRIM(inp_para)//'_hist_POST_FILTER.csv'
        histo_avg_fn_pre__Filter = TRIM(basename)//'_'//TRIM(inp_para)//'_hist_avg_PRE__FILTER.csv'
        histo_avg_fn_post_Filter = TRIM(basename)//'_'//TRIM(inp_para)//'_hist_avg_POST_FILTER.csv'
        histogram_filename_tex_Filter  = TRIM(basename)//'_'//TRIM(inp_para)//'_Filter_Histogram.tex'
ENDIF ! (my_rank .EQ. 0)

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
CALL r3_array_sectioning (domains=size_mpi, sections=sections, domain=my_rank, rank_section=rank_section)

! Calculate Padding to decrease "size of array" to a corresponding size
! Size if Kernel always an odd number. Num-1 = Voxels on both sides of filtered Voxel
original_image_padding =  kernel_spec(2)-1_ik             ! 1D Array
border                 = (kernel_spec(2)-1_ik) / 2_ik     ! 0D Array (scalar)

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

subarray_dims         = (dims_reduced / sections)
subarray_dims_overlap = subarray_dims + original_image_padding

IF ( (debug .GE. 1_ik) .AND. (my_rank .EQ. 0_ik) ) THEN
        WRITE(rd_o,'(A)')      std_lnbrk
        WRITE(rd_o,'(A)')      "Calculation of domain sectioning:"
        WRITE(rd_o,'(A)')
        WRITE(rd_o,'(A, 3I5)') "sections:               ", sections
        WRITE(rd_o,'(A, 3I5)') "dims:                   ", dims
        WRITE(rd_o,'(A, 3I5)') "border:                 ", border
        WRITE(rd_o,'(A, 3I5)') "original_image_padding: ", original_image_padding
        WRITE(rd_o,'(A, 3I5)') "remainder_per_dir:      ", remainder_per_dir
        WRITE(rd_o,'(A, 3I5)') "dims_reduced:           ", dims_reduced
        WRITE(rd_o,'(A, 3I5)') "subarray_dims:          ", subarray_dims
        WRITE(rd_o,'(A)')      std_lnbrk
END IF

ALLOCATE( subarray(subarray_dims_overlap(1), subarray_dims_overlap(2), subarray_dims_overlap(3) ) )

subarray_origin = (rank_section-1_ik) * (subarray_dims) !+ 1_ik

CALL read_raw_mpi(      filename=filename                       , &
                        type=TRIM(typ)                          , &
                        hdr_lngth=INT(displacement, KIND=8)     , &
                        dims=dims                               , &
                        subarray_dims=subarray_dims_overlap     , &
                        subarray_origin=subarray_origin         , &
                        subarray=subarray                       , &
                        displacement=displacement               , &
                        log_un=rd_o                             , &
                        status_o=status)

IF (status .EQ. 1_ik) THEN
        WRITE(rd_o,'(A)')  'Something during MPI File read went wrong. Please check/debug.'
        CLOSE(rd_o)
        CALL MPI_ABORT (MPI_COMM_WORLD, 1_mik, ierr)
END IF

! Prepare collecting the subarrays to assemble a global vtk file.
IF (my_rank .EQ. 0_ik) CALL CPU_TIME(read_t_vtk)

! subarray_reduced_bndaries                   ! No Overlap
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

hbnds        = (/ histo_bnd_global_lo, histo_bnd_global_hi , histo_bnd_global_hi - histo_bnd_global_lo /)

! Prior to image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (   subarray(       srb(1):srb(4)   , & 
                                                        srb(2):srb(5)   , &
                                                        srb(3):srb(6))  , &
                                                        hbnds           , &
                                                        histogram_pre__F)            

IF (my_rank .EQ. 0_ik) CALL CPU_TIME(prep_Histo)

! Start image processing
! result_image is necessary, because otherwise, filtered Voxels will be used for filtering following voxels.
! Therefore, doesn't really work in place
! Ensure image padding manually in front of this branch
IF (kernel_spec(1) .EQ. 2_ik) THEN
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

IF (my_rank .EQ. 0_ik) CALL CPU_TIME(calculation)

! After image filtering
! Get Histogram of Scalar Values
CALL extract_histogram_scalar_array (result_subarray, hbnds, histogram_post_F)            

! Allocate memory for global histogram
IF (my_rank .EQ. 0_ik) THEN
        ALLOCATE(histogram_pre__F_global(hbnds(1):hbnds(2)))
        ALLOCATE(histogram_post_F_global(hbnds(1):hbnds(2)))
END IF

! Collect the data of the histogram post filtering       
CALL MPI_REDUCE (histogram_pre__F, histogram_pre__F_global, INT(SIZE(histogram_pre__F), KIND=mik), &
        MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_REDUCE (histogram_post_F, histogram_post_F_global, INT(SIZE(histogram_post_F), KIND=mik), &
        MPI_INT, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

CALL CPU_TIME(extract_Histo)

IF (my_rank .EQ. 0_ik) THEN

        ! Export Histograms
        CALL write_histo_csv (fh_pre,      histogram_fn_pre__Filter, hbnds, 0_ik       , histogram_pre__F_global)
        CALL write_histo_csv (fh_post,     histogram_fn_post_Filter, hbnds, 0_ik       , histogram_post_F_global)
        CALL write_histo_csv (fh_pre__avg, histo_avg_fn_pre__Filter, hbnds, mov_avg_width, histogram_pre__F_global)
        CALL write_histo_csv (fh_post_avg, histo_avg_fn_post_Filter, hbnds, mov_avg_width, histogram_post_F_global)

        CALL write_tex_for_histogram (fun3      , &
                histogram_filename_tex_Filter   , &
                histogram_fn_pre__Filter        , &
                histogram_fn_post_Filter        , &
                histo_avg_fn_pre__Filter        , &
                histo_avg_fn_post_Filter        )

        filenameExportVtk = TRIM(basename)//'_'//TRIM(inp_para)//'.vtk'

        CALL write_vtk_meta (   fh=fh_data_out                          , &
                                filename=filenameExportVtk              , & 
                                type=TRIM(typ)                          , &
                                atStart=.TRUE.                          , &
                                spcng=spcng                             , &
                                dims=sections*subarray_dims)

        INQUIRE(FILE=filenameExportVtk, SIZE=wr_vtk_hdr_lngth)

END IF ! (my_rank .EQ. 0_ik)

        ! BCAST used in some way of a Barrier.
        CALL MPI_BCAST (filenameExportVtk, INT(mcl, KIND=mik), MPI_CHAR   , 0_mik, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST (wr_vtk_hdr_lngth , 1_mik             , MPI_INTEGER, 0_mik, MPI_COMM_WORLD, ierr)

        CALL write_raw_mpi (    type=TRIM(typ)                          , &
                                hdr_lngth=INT(wr_vtk_hdr_lngth, KIND=8) , &
                                filename=filenameExportVtk              , &
                                dims=sections*subarray_dims             , &
                                subarray_dims=subarray_dims             , &
                                subarray_origin=subarray_origin         , &
                                subarray=INT(result_subarray, KIND=INT16))

        DEALLOCATE(result_subarray)

IF (my_rank .EQ. 0_ik) THEN

        CALL write_vtk_meta (   fh=fh_data_out                          , &
                                filename=filenameExportVtk              , & 
                                atStart=.FALSE.)

        CALL CPU_TIME(global_finish)

        WRITE(rd_o,'(A         )')  std_lnbrk
        WRITE(rd_o,'(A, F13.3, A)') 'Init and parsing     = ', (init_finish        - global_start)                     ,' Seconds'
        WRITE(rd_o,'(A, F13.3, A)') 'Read File            = ', (read_t_vtk         - init_finish)                      ,' Seconds'
        WRITE(rd_o,'(A, F13.3, A)') 'Prep Histograms      = ', (prep_Histo - read_t_vtk)                               ,' Seconds'
        WRITE(rd_o,'(A, F13.3, A)') 'Calculation          = ', (calculation        - prep_Histo)                       ,' Seconds'
        WRITE(rd_o,'(A, F13.3, A)') 'Extract Histograms   = ', (extract_Histo - calculation)                           ,' Seconds'
        WRITE(rd_o,'(A, F13.3, A)') 'Calculate Histograms = ', (prep_Histo - read_t_vtk + extract_Histo - calculation) ,' Seconds'
        WRITE(rd_o,'(A, F13.3, A)') 'Write all data       = ', (global_finish      - extract_Histo)                    ,' Seconds'
        WRITE(rd_o,'(A, F13.3, A)') 'Overall Time         = ', (global_finish      - global_start)                     ,' Seconds'
        WRITE(rd_o,'(A, F13.3, A)') 'Overall Time         = ', (global_finish      - global_start) / 60                ,' Minutes'
        WRITE(rd_o,'(A         )')  std_lnbrk
        WRITE(rd_o,'(A, F13.3, A)') 'CPU time             = ', (global_finish      - global_start) / 60 / 60 * size_mpi,' Hours'

        CLOSE(rd_o)     
ENDIF 

CALL MPI_FINALIZE(ierr)
CALL MPI_ERR(ierr,"MPI_FINALIZE didn't succeed")

END PROGRAM main