!------------------------------------------------------------------------------
! MODULE: aux_routines_ip 
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! DESCRIPTION: 
!> Module containing auxiliary routines for image processing
!------------------------------------------------------------------------------

MODULE aux_routines_ip

USE global_std

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: extract_histogram_scalar_array
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Extracts a histogram from a 3-dimensional array
!
!> @param[in] array Actual image
!> @param[in] hbnds Boundary values
!> @param[out] histogram Returns the histogram
!------------------------------------------------------------------------------  
SUBROUTINE extract_histogram_scalar_array (array, hbnds, histogram)
! This is an inherently unflexible subroutine. It delivers exactly this kind of histogram and nothing else...
INTEGER(KIND=ik), DIMENSION(:,:,:) :: array
INTEGER(KIND=ik), DIMENSION(3)             , INTENT(IN)  :: hbnds    ! histogram lower/upper bounds
INTEGER(KIND=ik), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: histogram

! Internal variables
INTEGER(KIND=ik) :: ii, jj, kk
INTEGER(KIND=ik), DIMENSION(3) :: shp

ALLOCATE(histogram(hbnds(1):hbnds(2)))

histogram(:) = 0_ik

shp = SHAPE(array)

! Take care of sign of hmin!! Not that intuitive
DO ii=1, shp(1)
  DO jj=1, shp(2)
    DO kk=1, shp(3)
      histogram( array(ii, jj, kk) ) = histogram(  array(ii, jj, kk) ) + 1_ik
    END DO
  END DO
END DO

END SUBROUTINE extract_histogram_scalar_array

!------------------------------------------------------------------------------
! SUBROUTINE: get_rank_section
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Extracts a histogram from a 3-dimensional array
!
!> @description
!> !!! May be obsolete! Check MPI intrinsics
!
!> @param[in] domain Control volume to address
!> @param[in] sections Actual decomposition
!> @param[out] rank_section Returns the specific section
!------------------------------------------------------------------------------  
SUBROUTINE get_rank_section (domain, sections, rank_section)
  ! Three-Dimensional Array Scatter Routine
  ! Calculate how many sections per direction are ideal to split an array via mpi
  ! It was written to scatter a 3D Array, however it may suit to differnt purposes
  ! IT IS HIGHLY RECOMMENDED TO PROVIDE A COMMON AMOUNT OF PROCESSORS AS SOME MAY BE UNUSED due to the nature of 3D arrays.
  
  INTEGER(KIND=ik), INTENT(IN)                :: domain
  INTEGER(KIND=ik), DIMENSION(3), INTENT(IN)  :: sections
  INTEGER(KIND=ik), DIMENSION(3), INTENT(OUT) :: rank_section
              
  ! Internal Variables              
  INTEGER(KIND=ik) :: rank, yremainder, zremainder

  ! Power of 2 is handled here, because with the algorithm of CASE DEFAULT, Greedy suboptimality kicks in!  
  ! Have a look at the corresponding Matlab/Octave testing file!
  ! In example at size_mpi = 128 Processors, where CASE DEFAULT will deliver 125 Processors!

  IF ( domain .EQ. 0_ik ) THEN
    rank_section = (/ 1_ik, 1_ik, 1_ik /)
  ELSE
    rank = domain + 1_ik ! MPI start at 0
    ! Calculate the rank_section out of my_rank and sections (/ x, y, z /)
    ! Tested via Octave. Not fully implemented by 20210503
    zremainder = MODULO(rank, sections(1)*sections(2))
    IF (zremainder .EQ. 0_ik) THEN
            rank_section = (/ sections(1), sections(2), (rank - zremainder) / (sections(1)*sections(2)) /)
    ELSE
            rank_section(3) = (rank - zremainder) / (sections(1) * sections(2)) 
    yremainder = MODULO(zremainder, sections(1))
    
            IF (yremainder .EQ. 0_ik) THEN
                    rank_section = (/ sections(1), (zremainder - yremainder) / sections(1), rank_section(3)+1 /)
            ELSE
                    rank_section = (/ yremainder, (zremainder - yremainder) / sections(1) + 1_ik, rank_section(3) + 1_ik /)
            ENDIF
    ENDIF
END IF
END SUBROUTINE get_rank_section

!------------------------------------------------------------------------------
! SUBROUTINE: write_tex_for_histogram
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Writes a tex histogram to file
!
!> @param[in] fun file unit / file handle
!> @param[in] fn_tex Actual decomposition
!> @param[in] fn_pre File name pre filter histogram
!> @param[in] fn_post File name post filter histogram
!> @param[in] fn_pre_avg File name pre filter averaged histogram
!> @param[in] fn_post_avg File name post filter averaged histogram
!------------------------------------------------------------------------------  
SUBROUTINE write_tex_for_histogram (fun, fn_tex, fn_pre, fn_post, fn_pre_avg, fn_post_avg)

  INTEGER    (KIND = ik) , INTENT(IN)        :: fun
  CHARACTER  (LEN  = mcl), INTENT(IN)        :: fn_tex, fn_pre, fn_post, fn_pre_avg, fn_post_avg
  CHARACTER  (LEN  = mcl)                    :: title

  title = TRIM(fn_tex(1:(LEN_TRIM(fn_tex) - 4_ik )))

  CALL underscore_to_blank(title, title)
  CALL basepath(title           , title)
  CALL basepath(fn_pre          , fn_pre)
  CALL basepath(fn_post         , fn_post)
  CALL basepath(fn_pre_avg      , fn_pre_avg)
  CALL basepath(fn_post_avg     , fn_post_avg)

  OPEN( UNIT = fun, file = TRIM(fn_tex), action="WRITE", status="new")

  WRITE(fun, '(A)')  "\documentclass{standalone}"
  WRITE(fun, '(A)')  "\usepackage{pgfplots}"
  WRITE(fun, '(A)')  ""
  WRITE(fun, '(A)')  "\definecolor{hlrsblue1}{RGB}{40, 172, 226}" 
  WRITE(fun, '(A)')  "\definecolor{hlrsblue2}{RGB}{106, 206, 248}" 
  WRITE(fun, '(A)')  "\definecolor{hlrsblue3}{RGB}{161, 224, 251}" 
  WRITE(fun, '(A)')  "\definecolor{hlrsblue4}{RGB}{208, 240, 253}" 
  WRITE(fun, '(A)')  "\definecolor{hlrsgray1}{RGB}{128, 128, 128}" 
  WRITE(fun, '(A)')  "\definecolor{hlrsgray2}{RGB}{160, 160, 160}" 
  WRITE(fun, '(A)')  "\definecolor{hlrsgray3}{RGB}{191, 191, 191}" 
  WRITE(fun, '(A)')  "\definecolor{hlrsgray4}{RGB}{210, 210, 210}" 
  WRITE(fun, '(A)')  ""
  WRITE(fun, '(A)')  "\begin{document}"
  WRITE(fun, '(A)')  ""
  WRITE(fun, '(A)')  "\begin{tikzpicture}"
  WRITE(fun, '(A)')  "    \begin{axis}["
  WRITE(fun, '(A)')  "        % xmode=log,"
  WRITE(fun, '(A)')  "        ymode=log,"
  WRITE(fun, '(A)')  "        xlabel=$scaledHU$,"
  WRITE(fun, '(A)')  "        ylabel=$Amount\space of\space Voxels$ (-),"
  WRITE(fun, '(3A)') "        title=",ADJUSTL(TRIM(title)),","
  WRITE(fun, '(A)')  "        grid=both,"
  WRITE(fun, '(A)')  "        minor grid style={gray!15},"
  WRITE(fun, '(A)')  "        major grid style={gray!15},"
  WRITE(fun, '(A)')  "        width=0.75\linewidth,"
  WRITE(fun, '(A)')  "        legend style={at={(1.03,0.5)},anchor=west},"
  WRITE(fun, '(A)')  "        legend cell align={left},"
  WRITE(fun, '(A)')  "        no marks]"
  WRITE(fun, '(A)')  "    \addplot[line width=1pt,solid,color=hlrsgray4] %"
  WRITE(fun, '(3A)') "        table[x=scaledHU,y=Voxels,col sep=comma]{", ADJUSTL(TRIM(fn_pre)),"};"
  WRITE(fun, '(A)')  "    \addlegendentry{Raw};"
  WRITE(fun, '(A)')  "    \addplot[line width=1pt,solid,color=hlrsblue4] %"
  WRITE(fun, '(3A)') "        table[x=scaledHU,y=Voxels,col sep=comma]{", ADJUSTL(TRIM(fn_post)),"};"
  WRITE(fun, '(A)')  "    \addlegendentry{Filtered};"
  WRITE(fun, '(A)')  "    \addplot[line width=1pt,solid,color=hlrsgray1] %"
  WRITE(fun, '(3A)') "        table[x=scaledHU,y=Voxels,col sep=comma]{", ADJUSTL(TRIM(fn_pre_avg)),"};"
  WRITE(fun, '(A)')  "    \addlegendentry{Raw, averaged};"
  WRITE(fun, '(A)')  "    \addplot[line width=1pt,solid,color=hlrsblue1] %"
  WRITE(fun, '(3A)') "        table[x=scaledHU,y=Voxels,col sep=comma]{", ADJUSTL(TRIM(fn_post_avg)),"};"
  WRITE(fun, '(A)')  "    \addlegendentry{Filtered, averaged};"
  WRITE(fun, '(A)')  "    \end{axis}"
  WRITE(fun, '(A)')  "    \end{tikzpicture}"
  WRITE(fun, '(A)')  ""
  WRITE(fun, '(A)')  "\end{document}"

  CLOSE(fun)

END SUBROUTINE write_tex_for_histogram

!------------------------------------------------------------------------------
! SUBROUTINE: underscore_to_blank
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Replaces underscores with blanks
!
!> @param[in] instring Input string
!> @param[in] outstring Output string
!------------------------------------------------------------------------------  
SUBROUTINE underscore_to_blank (instring, outstring)
  ! This whole subroutine is a workaround :-)
  CHARACTER(LEN=*) :: instring
  CHARACTER(LEN=*) :: outstring
  INTEGER(KIND=ik) :: ii

  outstring=instring
  DO ii=1, LEN_TRIM(instring)
    IF (instring(ii:ii) == '_')  outstring(ii:ii) = ' '
  END DO

  outstring=ADJUSTL(TRIM(outstring))
END SUBROUTINE underscore_to_blank

!------------------------------------------------------------------------------
! SUBROUTINE: basepath
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Returns the base path
!
!> @param[in] instring Input string (of the path)
!> @param[in] outstring Output string (of the path)
!------------------------------------------------------------------------------  
SUBROUTINE basepath (instring, outstring)
  ! This whole subroutine is a workaround :-)
  CHARACTER(LEN=*) :: instring
  CHARACTER(LEN=*) :: outstring
  INTEGER(KIND=ik) :: ii, blanks

  outstring=instring
  DO ii=1, LEN_TRIM(instring)
    IF (instring(ii:ii) == '/') blanks = ii
  END DO
  outstring(1:blanks) = ' '

  outstring=ADJUSTL(TRIM(outfile))
END SUBROUTINE basepath

END MODULE aux_routines_IP
