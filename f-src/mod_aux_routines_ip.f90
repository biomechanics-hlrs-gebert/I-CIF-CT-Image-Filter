!-------------------------------------------------------------------------------------------------
! mod_aux_routines_IP.f90
! Module to extract a histogram out a Vtk structure. The Module may inherit other subroutines as well.
!
! Author          Johannes Gebert
! Date Original   12.01.2021
! Date Modified   09.05.2021

! SUBROUTINE extract_histogram_scalar_array (array, hbnds, histogram)
! SUBROUTINE r3_array_sectioning (domains, sections, domain, rank_section)
! SUBROUTINE write_tex_for_histogram (fun, fn_tex, fn_pre, fn_post)
! SUBROUTINE underscore_to_blank (infile, outfile)
! SUBROUTINE basepath (infile, outfile)

MODULE aux_routines_ip

USE ISO_FORTRAN_ENV
USE standards

IMPLICIT NONE

CONTAINS

SUBROUTINE extract_histogram_scalar_array (array, hbnds, histogram)
! This is an inherently unflexible subroutine. It delivers exactly this kind of histogram and nothing else...
INTEGER  (KIND=ik)    , DIMENSION(:,:,:)                                           :: array
INTEGER  (KIND=ik)    , DIMENSION(3)                           , INTENT(IN)        :: hbnds    ! histogram lower/upper bounds
INTEGER  (KIND=ik)    , DIMENSION(:)    , ALLOCATABLE          , INTENT(OUT)       :: histogram

! Internal variables
INTEGER  (KIND=ik)                                                                 :: ii, jj, kk
INTEGER  (KIND=ik)    , DIMENSION(3)                                               :: shp

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


SUBROUTINE r3_array_sectioning (domains, sections, domain, rank_section)
  ! Three-Dimensional Array Scatter Routine
  ! Calculate how many sections per direction are ideal to split an array via mpi
  ! It was written to scatter a 3D Array, however it may suit to differnt purposes
  ! IT IS HIGHLY RECOMMENDED TO PROVIDE A COMMON AMOUNT OF PROCESSORS AS SOME MAY BE UNUSED due to the nature of 3D arrays.
  
  INTEGER(KIND = ik)                        , INTENT(IN)  :: domains
  INTEGER(KIND = ik), DIMENSION(3)          , INTENT(OUT) :: sections
  INTEGER(KIND = ik)              , OPTIONAL, INTENT(IN)  :: domain
  INTEGER(KIND = ik), DIMENSION(3), OPTIONAL, INTENT(OUT) :: rank_section

  ! Internal Variables
  INTEGER(KIND = ik)                                      :: sw=0, true_size, rank
  INTEGER(KIND = ik)                                      :: yremainder, zremainder

  ! Power of 2 is handled here, because with the algorithm of CASE DEFAULT, Greedy suboptimality kicks in!  
  ! Have a look at the corresponding Matlab/Octave testing file!
  ! In example at size_mpi = 128 Processors, where CASE DEFAULT will deliver 125 Processors!
  SELECT CASE(domains)
  CASE(2);     sections = (/  2,  1,  1 /)
  CASE(3);     sections = (/  3,  1,  1 /) ! Implementing it prevents crashing :)
  CASE(4);     sections = (/  2,  2,  1 /)
  CASE(5);     sections = (/  5,  1,  1 /) ! Implementing it prevents crashing :)
  CASE(6);     sections = (/  3,  2,  1 /)
  CASE(7);     sections = (/  7,  2,  1 /) ! Implementing it prevents crashing :)
  CASE(8);     sections = (/  2,  2,  2 /)
  CASE(10);    sections = (/  5,  2,  1 /)
  CASE(12);    sections = (/  3,  2,  2 /)
  CASE(16);    sections = (/  4,  2,  2 /)
  CASE(20);    sections = (/  5,  2,  2 /)
  CASE(24);    sections = (/  4,  3,  2 /)
  CASE(36);    sections = (/  4,  3,  3 /)
  CASE(40);    sections = (/  5,  4,  2 /)
  CASE(48);    sections = (/  4,  4,  3 /)
  CASE(64);    sections = (/  4,  4,  4 /)
  CASE(72);    sections = (/  6,  4,  3 /)
  CASE(80);    sections = (/  5,  4,  4 /)
  CASE(96);    sections = (/  6,  4,  4 /)
  CASE(120);   sections = (/  6,  5,  4 /)
  CASE(128);   sections = (/  8,  4,  4 /)
  CASE(144);   sections = (/  6,  6,  4 /)
  CASE(160);   sections = (/  8,  5,  4 /)
  CASE(216);   sections = (/  6,  6,  6 /)
  CASE(256);   sections = (/  8,  8,  4 /)
  CASE(512);   sections = (/  8,  8,  8 /)
  CASE(1024);  sections = (/ 16,  8,  8 /)
  CASE(2048);  sections = (/ 16, 16,  8 /)
  CASE(4096);  sections = (/ 16, 16, 16 /)
  CASE(8192);  sections = (/ 32, 16, 16 /)
  CASE(16384); sections = (/ 32, 32, 16 /)
  CASE(32768); sections = (/ 32, 32, 32 /)
  CASE DEFAULT
    ! If no option fits - a power of 2 will suffice. A lot of the processors may not get utilized then(!)
    true_size  = domains - MODULO(domains, 2_ik)
    sections   = (/ 1, 1, 1 /)
    DO WHILE (PRODUCT(sections) < true_size)
                                             sections(1) = sections(1) + 1_ik; sw=1_ik
          IF (PRODUCT(sections) < true_size) sections(2) = sections(2) + 1_ik; sw=2_ik
          IF (PRODUCT(sections) < true_size) sections(3) = sections(3) + 1_ik; sw=3_ik
    END DO
    sections(sw) = sections(sw) - 1_ik 
  END SELECT

  IF (PRESENT(domain)) THEN

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
END IF
END SUBROUTINE r3_array_sectioning

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

!---------------------------------------------------------------------------------------------------

SUBROUTINE underscore_to_blank (infile, outfile)
  ! This whole subroutine is a workaround :-)
  CHARACTER  (LEN = *)       :: infile
  CHARACTER  (LEN = *)       :: outfile
  INTEGER    (KIND = ik)     :: ii, blanks

  outfile=infile
  DO ii=1, LEN_TRIM(infile)
    IF (infile(ii:ii) == '_')  outfile(ii:ii) = ' '
  END DO

  outfile=ADJUSTL(TRIM(outfile))
END SUBROUTINE underscore_to_blank

!---------------------------------------------------------------------------------------------------

SUBROUTINE basepath (infile, outfile)
  ! This whole subroutine is a workaround :-)
  CHARACTER  (LEN = *)      :: infile
  CHARACTER  (LEN = *)      :: outfile
  INTEGER    (KIND = ik)    :: ii, blanks

  outfile=infile
  DO ii=1, LEN_TRIM(infile)
    IF (infile(ii:ii) == '/')  blanks         = ii
  END DO
  outfile(1:blanks) = ' '

  outfile=ADJUSTL(TRIM(outfile))
END SUBROUTINE basepath

END MODULE aux_routines_IP
