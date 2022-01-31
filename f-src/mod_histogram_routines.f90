!------------------------------------------------------------------------------
! MODULE: histogram_routines 
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief: 
!> Auxiliary routines for calculating and writing the histograms.
!------------------------------------------------------------------------------
MODULE histogram_routines

USE ISO_FORTRAN_ENV
USE global_std
USE meta
USE formatted_plain

IMPLICIT NONE

INTERFACE extract_histogram_scalar_array
   MODULE PROCEDURE extract_histogram_scalar_array_ik2
   MODULE PROCEDURE extract_histogram_scalar_array_ik4
END INTERFACE extract_histogram_scalar_array

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: extract_histogram_scalar_array_ik2
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Extracts a histogram from a 3-dimensional array
!> !!! Scaling may need extrensive Rework !!! 
!
!> @param[in] array Actual image
!> @param[in] hbnds Boundary values
!> @param[out] histogram Returns the histogram
!------------------------------------------------------------------------------  
SUBROUTINE extract_histogram_scalar_array_ik2(array, hbnds, histogram)
! This is an inherently unflexible subroutine. It delivers exactly this kind of histogram and nothing else...
INTEGER(KIND=INT16), DIMENSION(:,:,:) :: array
INTEGER(KIND=ik), DIMENSION(3), INTENT(IN)  :: hbnds    ! histogram lower/upper bounds
INTEGER(KIND=ik), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: histogram

! Internal variables
INTEGER(KIND=ik) :: ii, jj, kk
INTEGER(KIND=ik), DIMENSION(3) :: shp

ALLOCATE(histogram(hbnds(1):hbnds(2)))
histogram(:) = 0_ik

shp = INT(SHAPE(array), KIND=ik)

! Take care of sign of hmin!! Not that intuitive
DO kk=1_ik, shp(3)
DO jj=1_ik, shp(2)
DO ii=1_ik, shp(1)
  histogram(array(ii, jj, kk)) = histogram(array(ii, jj, kk)) + 1_ik
END DO
END DO
END DO

END SUBROUTINE extract_histogram_scalar_array_ik2

!------------------------------------------------------------------------------
! SUBROUTINE: extract_histogram_scalar_array_ik4
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Extracts a histogram from a 3-dimensional array
!> !!! Scaling may need extrensive Rework !!! 
!
!> @param[in] array Actual image
!> @param[in] hbnds Boundary values
!> @param[out] histogram Returns the histogram
!------------------------------------------------------------------------------  
SUBROUTINE extract_histogram_scalar_array_ik4(array, hbnds, histogram)
! This is an inherently unflexible subroutine. It delivers exactly this kind of histogram and nothing else...
INTEGER(KIND=INT32), DIMENSION(:,:,:) :: array
INTEGER(KIND=ik), DIMENSION(3)             , INTENT(IN)  :: hbnds    ! histogram lower/upper bounds
INTEGER(KIND=ik), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: histogram

! Internal variables
INTEGER(KIND=ik) :: ii, jj, kk
INTEGER(KIND=ik), DIMENSION(3) :: shp

ALLOCATE(histogram(hbnds(1):hbnds(2)))

histogram(:) = 0_ik

shp = INT(SHAPE(array), KIND=ik)

DO kk=1, shp(3)
DO jj=1, shp(2)
DO ii=1, shp(1)
  histogram(array(ii, jj, kk)) = histogram(array(ii, jj, kk)) + 1_ik
END DO
END DO
END DO

END SUBROUTINE extract_histogram_scalar_array_ik4

!------------------------------------------------------------------------------
! SUBROUTINE: write_histo_csv
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Write the csv file of the histogram. Provide an odd number for the moving
!> average to get the best result out of it.
!
!> @param[in] fh File handle
!> @param[in] hdr_str String of the histograms header
!> @param[in] hbnds Histogram boundaries
!> @param[in] mov_avg_width Gliding average - width of smoothing
!> @param[in] histogram Actual histogram data
!------------------------------------------------------------------------------  
 SUBROUTINE write_histo_csv (fh, hdr_str, hbnds, mov_avg_width, histogram)
  ! Arg_divider acts as a parameter defining a moving average (!)
  ! It has an immediate effect like a filtered graph.
  INTEGER(KIND=ik), INTENT(IN) :: fh
  CHARACTER(len=*), INTENT(IN) :: hdr_str
  INTEGER(KIND=ik), DIMENSION(3), INTENT(IN) :: hbnds ! histogram lower/upper bounds
  INTEGER(KIND=ik)              , INTENT(IN) :: mov_avg_width
  INTEGER(KIND=ik), DIMENSION(:), INTENT(IN) :: histogram

  INTEGER(KIND=ik) :: ii, avg, span, step, huwritten
  
  !------------------------------------------------------------------------------
  ! Choose steps of loop according to the histogram boundaries.
  !------------------------------------------------------------------------------
  IF (mov_avg_width < 1) THEN
    step = 1
    span = 0
  ELSE
    step = CEILING(REAL(mov_avg_width, KIND=rk)/10._rk, KIND=ik) ! Factor ~10 or 1 
    span = (mov_avg_width-1_ik)/2_ik                             ! int division
  END IF

  !------------------------------------------------------------------------------
  ! Check histogram boundaries to adjust do loop.
  !------------------------------------------------------------------------------
  huwritten = hbnds(1)

  WRITE(fh,'(A)') TRIM(ADJUSTL(hdr_str))
  
  DO ii = 1+span, hbnds(2)-span, step
    avg = SUM(histogram(ii-span:ii+span))/(mov_avg_width + 1_ik)

    WRITE(fh,'(I0,A,I0)') huwritten," , ",avg
    
    huwritten = huwritten + step
  END DO

 END SUBROUTINE write_histo_csv

!------------------------------------------------------------------------------
! SUBROUTINE: write_tex_for_histogram
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Writes a tex histogram to file. The file was connected with the meta file
!> format.
!
!> @param[in] fun file unit / file handle
!> @param[in] suf_csv_prf Suffix for pre filter data
!> @param[in] suf_csv_pof Suffix for post filter data
!> @param[in] suf_csv_aprf Suffix for average pre filter data
!> @param[in] suf_csv_apof Suffix for average post filter data
!------------------------------------------------------------------------------  
SUBROUTINE write_tex_for_histogram (fun, suf_csv_prf, suf_csv_pof, suf_csv_aprf, suf_csv_apof)

  INTEGER(KIND=ik), INTENT(IN) :: fun
  CHARACTER(LEN=*), INTENT(IN) :: suf_csv_prf, suf_csv_pof, suf_csv_aprf, suf_csv_apof

  CHARACTER(LEN=mcl) :: title

  title = out%bsnm

  CALL underscore_to_blank(title, title)

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
  WRITE(fun, '(A)')  "        % ymode=log,"
  WRITE(fun, '(A)')  "        xlabel=$Grayscale (-)$,"
  WRITE(fun, '(A)')  "        ylabel=$Amount\space of\space Voxels$ (-),"
  WRITE(fun, '(3A)') "        title=",ADJUSTL(TRIM(title)),","
  WRITE(fun, '(A)')  "        grid=both,"
  WRITE(fun, '(A)')  "        minor grid style={gray!15},"
  WRITE(fun, '(A)')  "        major grid style={gray!15},"
  WRITE(fun, '(A)')  "        width=0.75\linewidth,"
  WRITE(fun, '(A)')  "        legend style={at={(1.03,0.5)},anchor=west},"
  WRITE(fun, '(A)')  "        legend cell align={left},"
  WRITE(fun, '(A)')  "        xtick={10000,15000,...,55000},"
  WRITE(fun, '(A)')  "        ytick={0,500000,...,2500000},"
  WRITE(fun, '(A)')  "        xmin=10000,"
  WRITE(fun, '(A)')  "        xmax=55000,"
  WRITE(fun, '(A)')  "        ymin=0,"
  WRITE(fun, '(A)')  "        ymax=2500000,"
  WRITE(fun, '(A)')  "        no marks]"
  WRITE(fun, '(A)')  "    \addplot[line width=1pt,solid,color=hlrsgray4] %"
  WRITE(fun, '(3A)') "        table[x=scaledHU,y=Voxels,col sep=comma]{", TRIM(out%bsnm)//TRIM(suf_csv_prf),"};"
  WRITE(fun, '(A)')  "    \addlegendentry{Raw};"
  WRITE(fun, '(A)')  "    \addplot[line width=1pt,solid,color=hlrsblue4] %"
  WRITE(fun, '(3A)') "        table[x=scaledHU,y=Voxels,col sep=comma]{", TRIM(out%bsnm)//TRIM(suf_csv_pof),"};"
  WRITE(fun, '(A)')  "    \addlegendentry{Filtered};"
  WRITE(fun, '(A)')  "    \addplot[line width=1pt,solid,color=hlrsgray1] %"
  WRITE(fun, '(3A)') "        table[x=scaledHU,y=Voxels,col sep=comma]{", TRIM(out%bsnm)//TRIM(suf_csv_aprf),"};"
  WRITE(fun, '(A)')  "    \addlegendentry{Raw, averaged};"
  WRITE(fun, '(A)')  "    \addplot[line width=1pt,solid,color=hlrsblue1] %"
  WRITE(fun, '(3A)') "        table[x=scaledHU,y=Voxels,col sep=comma]{", TRIM(out%bsnm)//TRIM(suf_csv_apof),"};"
  WRITE(fun, '(A)')  "    \addlegendentry{Filtered, averaged};"
  WRITE(fun, '(A)')  "    \end{axis}"
  WRITE(fun, '(A)')  "    \end{tikzpicture}"
  WRITE(fun, '(A)')  ""
  WRITE(fun, '(A)')  "\end{document}"

END SUBROUTINE write_tex_for_histogram

END MODULE histogram_routines
