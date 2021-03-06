!------------------------------------------------------------------------------
! PROGRAM: Computed Tomography Image Filter
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!------------------------------------------------------------------------------
PROGRAM CTIF

USE global_std
USE meta
USE raw_binary  
USE user_interaction
USE formatted_plain
USE MPI
USE kernels
USE histogram_routines

IMPLICIT NONE 

! Parameter
INTEGER(ik), PARAMETER :: debug = 2   ! Choose an even integer!!
INTEGER(ik), PARAMETER :: mov_avg_width = 100   ! Choose an even integer!!

INTEGER(mik) :: sections(3)

INTEGER(INT16), DIMENSION(:,:,:), ALLOCATABLE  :: subarray_ik2, result_subarray_ik2
INTEGER(INT32), DIMENSION(:,:,:), ALLOCATABLE  :: subarray_ik4, result_subarray_ik4

CHARACTER(mcl), DIMENSION(:), ALLOCATABLE :: m_rry
CHARACTER(scl) :: type='', selectKernel='', restart='', restart_cmd_arg='', &
    db_order='', datarep='', histog=''

    CHARACTER(  8) :: date=''
CHARACTER( 10) :: time=''

REAL(rk) :: global_start, init_finish, read_t_vtk, prep_Histo, &
    calculation, extract_Histo, global_finish, sigma, spcng(3)

REAL(rk), DIMENSION(:,:,:), ALLOCATABLE  :: kernel

CHARACTER(mcl) :: binary, cmd_arg_history='', stat=''
CHARACTER(scl) :: suf_csv_prf, suf_csv_pof, suf_csv_aprf, suf_csv_apof, suf_csv_fihi

INTEGER(ik) :: histo_bnd_global_lo, histo_bnd_global_hi, histo_bnd_local_lo, &
    histo_bnd_local_hi, fh_csv_prf, fh_csv_pof, fh_csv_aprf, fh_csv_apof, &
    border, kernel_size, srb(6) ! subarray_reduced_bndaries

INTEGER(ik), DIMENSION(3) :: hbnds, dims, in_img_padding, subarray_origin, &
    sections_ik, rank_section, srry_dims, &
    dims_reduced, rmndr_dir, srry_dims_overlap

INTEGER(ik), DIMENSION(:), ALLOCATABLE :: histogram_pre__F, histogram_post_F, &
    pre_F_global, post_F_global

LOGICAL :: abrt=.FALSE.

! MPI Variables
INTEGER(mik) :: ierr, my_rank, size_mpi

! Initialize MPI Environment
CALL MPI_INIT(ierr)
CALL MPI_ERR(ierr,"MPI_INIT didn't succeed")

CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_RANK couldn't be retrieved")

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL MPI_ERR(ierr,"MPI_COMM_SIZE couldn't be retrieved")

IF (size_mpi < 2) CALL print_err_stop(std_out, "At least two ranks required to execute this program.", 1)

! Initialize program itself
IF(my_rank == 0) THEN

    CALL CPU_TIME(global_start)

    !------------------------------------------------------------------------------
    ! Parse the command arguments
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, in%full, stat, restart_cmd_arg, cmd_arg_history)
    IF(stat /= '') GOTO 1001
    
    IF (in%full=='') THEN
        CALL usage(binary)    

        !------------------------------------------------------------------------------
        ! On std_out since file of std_out is not spawned
        !------------------------------------------------------------------------------
        CALL print_err_stop(6, "No input file given", 1)
    END IF

    !------------------------------------------------------------------------------
    ! Check and open the input file; Modify the Meta-Filename / Basename
    ! Define the new application name first
    !------------------------------------------------------------------------------
    global_meta_prgrm_mstr_app = 'cif' 
    global_meta_program_keyword = 'CT_IMAGE_FILTER'
    CALL meta_append(m_rry, size_mpi, stat)
    IF(stat/="") CALL print_err_stop(std_out, stat, 1)
    
    !------------------------------------------------------------------------------
    ! Redirect std_out into a file in case std_out is not useful by environment.
    ! Place these lines before handle_lock_file :-)
    !------------------------------------------------------------------------------
    std_out = determine_stout()

    !------------------------------------------------------------------------------
    ! Spawn standard out after(!) the basename is known
    !------------------------------------------------------------------------------
    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

    CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM) "])
 
    IF(debug >=0) WRITE(std_out, FMT_MSG) "Post mortem info probably in ./datasets/temporary.std_out"

    !------------------------------------------------------------------------------
    ! Parse input
    !------------------------------------------------------------------------------
    CALL meta_read('DATA_BYTE_ORDER', m_rry, db_order, stat); IF(stat/="") abrt=.TRUE.
    CALL meta_read('RESTART'   , m_rry, restart, stat); IF(stat/="") abrt=.TRUE.
    CALL meta_read('TYPE_RAW'  , m_rry, type, stat); IF(stat/="") abrt=.TRUE.
    CALL meta_read('SPACING'   , m_rry, spcng, stat); IF(stat/="") abrt=.TRUE.
    CALL meta_read('DIMENSIONS', m_rry, dims, stat); IF(stat/="") abrt=.TRUE.

    CALL meta_read('HISTOGRAM', m_rry, histog, stat); IF(stat/="") abrt=.TRUE.

    CALL meta_read('FILTER_SIZE'  , m_rry, kernel_size, stat); IF(stat/="") abrt=.TRUE.
    CALL meta_read('FILTER_KERNEL', m_rry, selectKernel, stat); IF(stat/="") abrt=.TRUE.
    CALL meta_read('FILTER_SIGMA' , m_rry, sigma, stat); IF(stat/="") abrt=.TRUE.

    IF (abrt) CALL print_err_stop(std_out, "A keyword failed.", 1)
    
    IF((type /= "ik2") .AND. (type /= "ik4")) THEN
        mssg = "Program only supports ik2 and ik4 for 'TYPE_RAW'"
        CALL print_err_stop(std_out, mssg, 1)
    END IF

    !------------------------------------------------------------------------------
    ! Check the endianess
    !------------------------------------------------------------------------------
    IF (db_order == "BigEndian") THEN
        datarep = "BIG_ENDIAN"
     ELSE
        datarep = "LITTLE_ENDIAN"
     END IF      

    !------------------------------------------------------------------------------
    ! Restart handling
    ! Done after meta_io to decide based on keywords
    !------------------------------------------------------------------------------
    CALL meta_handle_lock_file(restart, restart_cmd_arg)

    IF ((histog == "YES") .OR. (histog=='Y')) THEN
        !------------------------------------------------------------------------------
        ! Spawn files for the csv-data and for the tex environment
        !------------------------------------------------------------------------------
        suf_csv_prf = "_hist_pre_filter"//csv_suf
        suf_csv_pof = "_hist_post_filter"//csv_suf
        suf_csv_aprf = "_hist_avg_pre_filter"//csv_suf
        suf_csv_apof = "_hist_avg_post_filter"//csv_suf
        suf_csv_fihi = "_filter_histogram"//csv_suf

        fh_csv_prf  = give_new_unit()  
        CALL meta_start_ascii(fh_csv_prf, suf_csv_prf)
        
        fh_csv_pof  = give_new_unit()
        CALL meta_start_ascii(fh_csv_pof, suf_csv_pof)
        
        fh_csv_aprf = give_new_unit()
        CALL meta_start_ascii(fh_csv_aprf, suf_csv_aprf)
        
        fh_csv_apof = give_new_unit()
        CALL meta_start_ascii(fh_csv_apof, suf_csv_apof)

        CALL meta_start_ascii(fht, tex_suf)
    END IF 
    CALL DATE_AND_TIME(date, time)

    IF((debug >= 0) .AND. (my_rank == 0)) THEN
        WRITE(std_out, FMT_TXT) "Date: "//date//" [ccyymmdd]"
        WRITE(std_out, FMT_TXT) "Time: "//time//" [hhmmss.sss]"  
        WRITE(std_out, FMT_TXT) "Program invocation:"//TRIM(cmd_arg_history)          
        WRITE(std_out, FMT_TXT_SEP)
        WRITE(std_out, FMT_MSG_xAI0) "Debug Level:", debug
        WRITE(std_out, FMT_MSG_xAI0) "Processors:", size_mpi  
        WRITE(std_out, FMT_MSG_xAI0) "Filter Size:", kernel_size
        WRITE(std_out, FMT_MSG)     "Filter Kernel: "//TRIM(selectKernel)
        WRITE(std_out, FMT_TXT_SEP)
        FLUSH(std_out)
    END IF
    
    CALL CPU_TIME(init_finish)
      
ENDIF ! (my_rank == 0)

CALL MPI_BCAST( in%p_n_bsnm, INT(meta_mcl, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(out%p_n_bsnm, INT(meta_mcl, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(selectKernel, INT(scl, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(type        , INT(scl, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(datarep     , INT(scl, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(histog      , INT(scl, mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(sigma       , 1_mik, MPI_DOUBLE_PRECISION , 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(kernel_size , 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(dims, 3_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

! DEBUG INFORMATION
IF((debug >= 2) .AND. (my_rank == 0)) THEN
    WRITE(std_out,FMT_MSG) "Creating dimensions." 
    WRITE(std_out,FMT_MSG) "Data type: ", TRIM(type) 
    WRITE(std_out,FMT_MSG_xAF0) "Sigma selected: ", sigma 
END IF

!------------------------------------------------------------------------------
! Get sections per direction
!------------------------------------------------------------------------------
sections=0
CALL MPI_DIMS_CREATE(size_mpi, 3_mik, sections, ierr)
sections_ik = INT(sections, ik)

CALL get_rank_section(domain=INT(my_rank, ik), sections=sections_ik, rank_section=rank_section)

!------------------------------------------------------------------------------
! Calculate Padding to decrease "size of array" to a corresponding size
! Size if Kernel always an odd number. 
! Num-1 = Voxels on both sides of filtered Voxel
!------------------------------------------------------------------------------
in_img_padding = kernel_size-1_ik ! 1D Array
border = (kernel_size-1_ik) / 2_ik ! 0D Array (scalar)

!------------------------------------------------------------------------------
! Remainder per direction gets almost fully ignored (!) 
! It's assumed, that even if we split into 32768 processes [32, 32, 32]
! sections, max. 31 x 31 x 31 Voxel get lost (Worst Case). Considering large 
! input sets, which are imperative for utilizing large amounts of processors, 
! losing 31 voxel at ~ 2000 ... 4000 Voxel per direction is not an issue.
! On the other hand, Distribution of the array gets way easier and presumably 
! quicker.
!------------------------------------------------------------------------------
rmndr_dir = MODULO(dims, sections_ik)

!------------------------------------------------------------------------------
! If remainder per direction is larger than the padding, simply shift the 
! array to distribute to ranks into the center of array. Then add the border. 
! This way, no artifical padding is required.
!------------------------------------------------------------------------------
IF((rmndr_dir(1) <= in_img_padding(1)) .OR. & 
    (rmndr_dir(2) <= in_img_padding(2)) .OR. &
    (rmndr_dir(3) <= in_img_padding(3))) THEN
    dims_reduced   = dims - in_img_padding
ELSE
    dims_reduced   = dims - rmndr_dir
END IF

srry_dims     = (dims_reduced / sections)
srry_dims_overlap = srry_dims + in_img_padding
subarray_origin = ((rank_section-1_ik) * (srry_dims))

! DEBUG INFORMATION
IF((debug >= 2) .AND. (my_rank == 0)) THEN
    WRITE(std_out,FMT_MSG) "Calculation of domain sectioning:"
    WRITE(std_out,FMT_MSG)
    WRITE(std_out,FMT_MSG_AxI0) "sections: ", sections_ik
    WRITE(std_out,FMT_MSG_AxI0) "dims: ", dims
    WRITE(std_out,FMT_MSG_AxI0) "border: ", border
    WRITE(std_out,FMT_MSG_AxI0) "input_img_padding: ", in_img_padding
    WRITE(std_out,FMT_MSG_AxI0) "rmndr_dir: ", rmndr_dir
    WRITE(std_out,FMT_MSG_AxI0) "dims_reduced: ", dims_reduced
    WRITE(std_out,FMT_MSG_AxI0) "srry_dims: ", srry_dims
    WRITE(std_out,FMT_MSG_AxI0) "srry_dims_overlap: ", srry_dims_overlap
    WRITE(std_out,FMT_MSG_AxI0) "subarray_origin: ", subarray_origin
    WRITE(std_out,FMT_MSG_SEP)
    FLUSH(std_out)
END IF

!------------------------------------------------------------------------------
! Allocate memory and read data
!------------------------------------------------------------------------------
IF(TRIM(datarep) == "BIG_ENDIAN") THEN
    datarep = "external32"
ELSE
    datarep = "native"
END IF 
write(*,*) "my_rank: ", my_rank, "subarray_origin", subarray_origin
write(*,*) "my_rank: ", my_rank, "subarray_max", subarray_origin + srry_dims_overlap


SELECT CASE(type)
    CASE('ik2') 
        ALLOCATE(subarray_ik2(srry_dims_overlap(1), srry_dims_overlap(2), srry_dims_overlap(3)))
        ALLOCATE(result_subarray_ik2(srry_dims(1), srry_dims(2), srry_dims(3)))

        result_subarray_ik2 = 0

        CALL mpi_read_raw(TRIM(in%p_n_bsnm)//raw_suf, 0_8, dims, &
            srry_dims_overlap, subarray_origin, subarray_ik2, TRIM(datarep))

        !------------------------------------------------------------------------------
        ! Adjust the the data range of the Histogram to min/max of the subarray and
        ! exchange information with all other ranks to get the global min/max.
        !------------------------------------------------------------------------------
        IF ((histog == "YES") .OR. (histog=='Y')) THEN
            histo_bnd_local_lo = MINVAL(subarray_ik2)
            histo_bnd_local_hi = MAXVAL(subarray_ik2)
        END IF 

    CASE('ik4') 
        ALLOCATE(subarray_ik4(srry_dims_overlap(1), srry_dims_overlap(2), srry_dims_overlap(3)))
        ALLOCATE(result_subarray_ik4(srry_dims(1), srry_dims(2), srry_dims(3)))

        result_subarray_ik4 = 0

        CALL mpi_read_raw(TRIM(in%p_n_bsnm)//raw_suf, 0_8, dims, &
            srry_dims_overlap, subarray_origin, subarray_ik4, TRIM(datarep))

        IF ((histog == "YES") .OR. (histog=='Y')) THEN
            histo_bnd_local_lo = MINVAL(subarray_ik4)
            histo_bnd_local_hi = MAXVAL(subarray_ik4)
        END IF 
END SELECT

!------------------------------------------------------------------------------
! Prepare collecting the subarrays to assemble a global binary data.
!------------------------------------------------------------------------------
IF(my_rank == 0) CALL CPU_TIME(read_t_vtk)

srb(1:3) = 1_ik + border
srb(4:6) = srry_dims_overlap - border

CALL MPI_ALLREDUCE(histo_bnd_local_lo, histo_bnd_global_lo, 1_mik, MPI_INTEGER8, MPI_MIN, MPI_COMM_WORLD, ierr)
CALL MPI_ALLREDUCE(histo_bnd_local_hi, histo_bnd_global_hi, 1_mik, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, ierr)

hbnds=[histo_bnd_global_lo, histo_bnd_global_hi , histo_bnd_global_hi - histo_bnd_global_lo]  

IF(my_rank == 0) CALL CPU_TIME(prep_Histo)

!------------------------------------------------------------------------------
! Start image processing
! result_image is necessary, because otherwise, filtered Voxels will be used
! for filtering following voxels. Therefore, doesn't really work in place.
! Ensure image padding manually in front of this branch.
!------------------------------------------------------------------------------
ALLOCATE(kernel(kernel_size, kernel_size, kernel_size))

SELECT CASE(selectKernel)
    CASE("Gaussian"); CALL kernel_gauss_3d   (kernel, kernel_size, sigma)
    CASE DEFAULT;     CALL kernel_identity_3d(kernel, kernel_size)
END SELECT

! DEBUG INFORMAITON
IF((debug >= 2) .AND. (my_rank == 0)) THEN
    WRITE(std_out,FMT_MSG) "Kernel - slice of filtered voxel, XY plane."

    CALL write_matrix(std_out, "Kernel slice XY", mat=kernel(:,:,border+1), fmti='std', unit='-')

    WRITE(std_out,FMT_MSG) "Prior to image filtering."
    WRITE(std_out,FMT_MSG_SEP)
    WRITE(std_out,FMT_MSG) "Subarray reduced boundaries: "
    WRITE(std_out,FMT_MSG_AxI0) "srb(1), srb(4): ", srb(1),  srb(4) 
    WRITE(std_out,FMT_MSG_AxI0) "srb(2), srb(5): ", srb(2),  srb(5) 
    WRITE(std_out,FMT_MSG_AxI0) "srb(3), srb(6): ", srb(3),  srb(6) 
    WRITE(std_out,FMT_MSG_SEP)
    FLUSH(std_out)
END IF

SELECT CASE(type)
    CASE('ik2') 
        IF ((histog == "YES") .OR. (histog=='Y')) THEN

            !------------------------------------------------------------------------------
            ! Prior to image filtering
            ! Get Histogram of Scalar Values
            !------------------------------------------------------------------------------
            CALL extract_histogram_scalar_array(&
                subarray_ik2(srb(1):srb(4), srb(2):srb(5), srb(3):srb(6)), hbnds, histogram_pre__F)    
        END IF

        CALL filter(subarray_ik2, kernel, srb, result_subarray_ik2)
        DEALLOCATE(subarray_ik2)

        !------------------------------------------------------------------------------
        ! After image filtering
        ! Get Histogram of Scalar Values
        !------------------------------------------------------------------------------
        IF ((histog == "YES") .OR. (histog=='Y')) CALL extract_histogram_scalar_array(&
            result_subarray_ik2, hbnds, histogram_post_F)        

    CASE('ik4') 
        IF ((histog == "YES") .OR. (histog=='Y')) CALL extract_histogram_scalar_array(&
                subarray_ik4(srb(1):srb(4), srb(2):srb(5), srb(3):srb(6)), hbnds, histogram_pre__F)    

        CALL filter(subarray_ik4, kernel, srb, result_subarray_ik4)
        DEALLOCATE(subarray_ik4)

        IF ((histog == "YES") .OR. (histog=='Y')) CALL extract_histogram_scalar_array(&
            result_subarray_ik4, hbnds, histogram_post_F)
END SELECT

DEALLOCATE(kernel)

!------------------------------------------------------------------------------
! Allocate memory for global histogram
!------------------------------------------------------------------------------
IF(my_rank == 0) THEN
    CALL CPU_TIME(calculation)

    IF ((histog == "YES") .OR. (histog=='Y')) THEN
        IF(debug >= 2) THEN
            WRITE(std_out,FMT_MSG) "Allocating memory for global histograms now."
            WRITE(std_out,FMT_MSG_AxI0) "hbnds: ", hbnds(1), hbnds(2), hbnds(3)  
            WRITE(std_out,FMT_MSG_SEP)
            FLUSH(std_out)
        END IF

        ALLOCATE(pre_F_global(hbnds(1):hbnds(2)))
        ALLOCATE(post_F_global(hbnds(1):hbnds(2)))
    END IF     
END IF

IF ((histog == "YES") .OR. (histog=='Y')) THEN

    !------------------------------------------------------------------------------
    ! Collect the data of the histogram post filtering       
    !------------------------------------------------------------------------------
    ! DEBUG INFORMATION
    IF((debug >= 2) .AND. (my_rank == 0)) THEN
        WRITE(std_out,FMT_MSG) "Collecting the data of the histogram post filter now."
        WRITE(std_out,FMT_MSG_SEP)
        FLUSH(std_out)
    END IF

    CALL MPI_REDUCE (histogram_pre__F, pre_F_global, INT(SIZE(histogram_pre__F), mik), &
        MPI_INTEGER8, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

    CALL MPI_REDUCE (histogram_post_F, post_F_global, INT(SIZE(histogram_post_F), mik), &
        MPI_INTEGER8, MPI_SUM, 0_mik, MPI_COMM_WORLD, ierr)

    CALL CPU_TIME(extract_Histo)
END IF

!------------------------------------------------------------------------------
! Write binary data. Since no other data types are required within this 
! doctoral project, no other types are implemented yet.
!------------------------------------------------------------------------------

! DEBUG INFORMATION
IF((debug >= 2) .AND. (my_rank == 0)) THEN
    WRITE(std_out,FMT_MSG) "Calculation of write raw mpi:"
    WRITE(std_out,FMT_MSG)
    WRITE(std_out,FMT_MSG_AxI0) "sections: ", sections
    WRITE(std_out,FMT_MSG_AxI0) "dims_reduced: ", dims_reduced
    WRITE(std_out,FMT_MSG_AxI0) "dims(_reduced): ", srry_dims*sections
    WRITE(std_out,FMT_MSG_AxI0) "srry_dims: ", srry_dims
    WRITE(std_out,FMT_MSG_AxI0) "subarray_origin: ", subarray_origin
    WRITE(std_out,FMT_MSG_SEP)
    FLUSH(std_out)
END IF

SELECT CASE(type)
    CASE('ik2') 
        CALL mpi_write_raw(TRIM(out%p_n_bsnm)//raw_suf, 0_8, srry_dims*sections, srry_dims, subarray_origin, result_subarray_ik2)
        DEALLOCATE(result_subarray_ik2)
    CASE('ik4') 
        CALL mpi_write_raw(TRIM(out%p_n_bsnm)//raw_suf, 0_8, srry_dims*sections, srry_dims, subarray_origin, result_subarray_ik4)
        DEALLOCATE(result_subarray_ik4)
END SELECT

!------------------------------------------------------------------------------
! Write final information
!------------------------------------------------------------------------------
IF(my_rank == 0) THEN
    CALL CPU_TIME(global_finish)

    WRITE(std_out,FMT_TXT_xAF0) 'Init and parsing     = ', (init_finish- global_start),' Seconds'
    WRITE(std_out,FMT_TXT_xAF0) 'Read File            = ', (read_t_vtk - init_finish) ,' Seconds'
    WRITE(std_out,FMT_TXT_xAF0) 'Prep Histograms      = ', (prep_Histo - read_t_vtk)  ,' Seconds'
    WRITE(std_out,FMT_TXT_xAF0) 'Calculation          = ', (calculation - prep_Histo) ,' Seconds'
    WRITE(std_out,FMT_TXT_xAF0) 'Extract Histograms   = ', (extract_Histo - calculation) ,' Seconds'
    WRITE(std_out,FMT_TXT_xAF0) 'Calculate Histograms = ', (prep_Histo - read_t_vtk + extract_Histo - calculation) ,' Seconds'
    WRITE(std_out,FMT_TXT_xAF0) 'Write all data       = ', (global_finish - extract_Histo),' Seconds'
    WRITE(std_out,FMT_TXT_xAF0) 'Overall Time         = ', (global_finish - global_start) ,' Seconds'
    WRITE(std_out,FMT_TXT_xAF0) 'Overall Time         = ', (global_finish - global_start) / 60,' Minutes'
    WRITE(std_out,FMT_TXT_SEP)  
    WRITE(std_out,FMT_TXT_xAF0) 'CPU time             = ', (global_finish - global_start) / 60 / 60 * size_mpi,' Hours'
    WRITE(std_out,FMT_MSG_SEP)

    IF ((histog == "YES") .OR. (histog=='Y')) THEN
        !------------------------------------------------------------------------------
        ! Export Histograms
        !------------------------------------------------------------------------------
        CALL write_histo_csv(fh_csv_prf,  "scaledHU, Voxels", hbnds, 0_ik, pre_F_global)
        CALL write_histo_csv(fh_csv_pof,  "scaledHU, Voxels", hbnds, 0_ik, post_F_global)
        CALL write_histo_csv(fh_csv_aprf, "scaledHU, Voxels", hbnds, mov_avg_width, pre_F_global)
        CALL write_histo_csv(fh_csv_apof, "scaledHU, Voxels", hbnds, mov_avg_width, post_F_global)

        CALL write_tex_for_histogram(fht, suf_csv_prf, suf_csv_pof, suf_csv_aprf, suf_csv_apof)
    END IF 
    
    !------------------------------------------------------------------------------
    ! Finish the program
    !------------------------------------------------------------------------------
    CALL meta_write('FIELD_OF_VIEW', '(mm)' , srry_dims*sections*spcng)
    CALL meta_write('DIMENSIONS', '(-)',      srry_dims*sections)
    CALL meta_write('ENTRIES', '(-)', PRODUCT(srry_dims*sections))
END IF

!------------------------------------------------------------------------------
! Only used in specific cases to finish more gracefully. (grep -i "GOTO")
!------------------------------------------------------------------------------
1001 Continue

IF(my_rank == 0) THEN

    CALL meta_signing(binary)
    CALL meta_close()

    IF ((histog == "YES") .OR. (histog=='Y')) THEN
        CALL meta_stop_ascii(fh_csv_prf, suf_csv_prf)
        CALL meta_stop_ascii(fh_csv_pof, suf_csv_pof)
        CALL meta_stop_ascii(fh_csv_aprf, suf_csv_aprf)
        CALL meta_stop_ascii(fh_csv_apof, suf_csv_apof)
    
        CALL meta_stop_ascii(fht, tex_suf)

    END IF 

    WRITE(std_out,FMT_TXT) 'Program finished successfully.'
    WRITE(std_out,FMT_MSG_SEP)

    IF(std_out/=6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

END IF ! (my_rank == 0)

CALL MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE didn't succeed", INT(ierr, ik))

END PROGRAM CTIF
