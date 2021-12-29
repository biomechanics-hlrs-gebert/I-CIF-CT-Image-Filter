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

  outstring=ADJUSTL(TRIM(outstring))
END SUBROUTINE basepath

END MODULE aux_routines_IP
