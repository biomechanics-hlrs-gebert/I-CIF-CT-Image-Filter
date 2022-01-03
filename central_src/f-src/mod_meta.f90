!------------------------------------------------------------------------------
! MODULE: meta
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @Description:
!> Module containing all meta file read/write routines.
!------------------------------------------------------------------------------
MODULE meta

   USE ISO_FORTRAN_ENV
   USE strings
   USE user_interaction

IMPLICIT NONE

   INTEGER, PARAMETER :: meta_ik = 8
   INTEGER, PARAMETER :: meta_rk = 8
   INTEGER, PARAMETER :: meta_mcl = 512
   INTEGER, PARAMETER :: meta_scl = 64

   ! Character lengths
   INTEGER, PARAMETER :: kcl    = 25   ! Keyword character  length
   INTEGER, PARAMETER :: ucl    = 8    ! Unit    character  length
   INTEGER, PARAMETER :: stdspc = 39   ! Keyword standard space

   CHARACTER(LEN=kcl) :: global_meta_program_keyword
   CHARACTER(LEN=kcl) :: global_meta_prgrm_mstr_app

   ! Standard files
   INTEGER(KIND=meta_ik), PARAMETER :: fh_meta_in  = 20, fhmei  = 20
   INTEGER(KIND=meta_ik), PARAMETER :: fh_meta_put = 21, fhmeo  = 21
   INTEGER(KIND=meta_ik), PARAMETER :: fh_mon      = 22, fhmon  = 22
   INTEGER(KIND=meta_ik), PARAMETER :: fh_out      = 23, fho    = 23
   INTEGER(KIND=meta_ik), PARAMETER :: fh_log      = 24, fhl    = 24
   INTEGER(KIND=meta_ik), PARAMETER :: fh_res      = 25, fhr    = 25
   INTEGER(KIND=meta_ik), PARAMETER :: fh_csv      = 26, fhc    = 26
   INTEGER(KIND=meta_ik), PARAMETER :: fh_head     = 27, fhh    = 27
   INTEGER(KIND=meta_ik), PARAMETER :: fh_tex      = 28, fht    = 28
   INTEGER(KIND=meta_ik), PARAMETER :: fh_vtk      = 29, fhv    = 29
   INTEGER(KIND=meta_ik), PARAMETER :: fh_raw      = 30, fhra   = 30
   CHARACTER(LEN=*), PARAMETER :: log_suf  = '.log'
   CHARACTER(LEN=*), PARAMETER :: lock_suf = '.lock'
   CHARACTER(LEN=*), PARAMETER :: head_suf = '.head'
   CHARACTER(LEN=*), PARAMETER :: meta_suf = '.meta'
   CHARACTER(LEN=*), PARAMETER :: mon_suf  = '.mon'
   CHARACTER(LEN=*), PARAMETER :: res_suf  = '.result'
   CHARACTER(LEN=*), PARAMETER :: csv_suf  = '.csv'
   CHARACTER(LEN=*), PARAMETER :: tex_suf  = '.tex'
   CHARACTER(LEN=*), PARAMETER :: vtk_suf  = '.vtk'
   CHARACTER(LEN=*), PARAMETER :: raw_suf  = '.raw'

   ! Meta data basename handling
   TYPE basename
      ! For the use in filenames, a max. length of a part of a basename of kcl characters must suffice.
      ! Nomenclature: dataset_type_purpose_app_features
      CHARACTER(LEN=meta_mcl) :: full     = '' ! Including suffix and path
      CHARACTER(LEN=meta_mcl) :: path     = '' ! Only the path to the file
      CHARACTER(LEN=meta_mcl) :: p_n_bsnm = '' ! Just the path and the basename
      CHARACTER(LEN=meta_mcl) :: bsnm     = '' ! Just the basename
      CHARACTER(LEN=kcl) :: dataset  = '' ! For example FH01-1 (Femoral Head 1, Scan1)
      CHARACTER(LEN=2)   :: type     = '' ! 'cl' - clinical or 'mu' - microfocus
      CHARACTER(LEN=3)   :: purpose  = '' ! 'Dev' or 'Pro' (Development or Production)
      CHARACTER(LEN=kcl) :: app      = '' ! Application. For example "Binarization"
      CHARACTER(LEN=kcl) :: features = '' ! Features. For example the parametrization
   END TYPE basename

   ! Always provide in/out for meta driven environments
   TYPE(basename) :: in, out

   !> Interface: meta_read
   !> \author Johannes Gebert
   !> \date 10.11.2021
   INTERFACE meta_read
      MODULE PROCEDURE meta_read_C 
      MODULE PROCEDURE meta_read_I0D 
      MODULE PROCEDURE meta_read_I1D
      MODULE PROCEDURE meta_read_R0D
      MODULE PROCEDURE meta_read_R1D
   END INTERFACE meta_read

   !> Interface: meta_write
   !> \author Johannes Gebert
   !> \date 10.11.2021
   INTERFACE meta_write
      MODULE PROCEDURE meta_write_C 
      MODULE PROCEDURE meta_write_I0D 
      MODULE PROCEDURE meta_write_R0D 
      MODULE PROCEDURE meta_write_I1D
      MODULE PROCEDURE meta_write_R1D
   END INTERFACE meta_write

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: meta_handle_lock_file
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to encapsule the lock file handling
!
!> @param[in] restart Whether to restart or not to.
!> @param[in] restart_cmdarg Possible cmd argument override
!------------------------------------------------------------------------------  
SUBROUTINE meta_handle_lock_file(restart, restart_cmdarg)

CHARACTER, INTENT(INOUT) :: restart
CHARACTER, INTENT(IN), OPTIONAL :: restart_cmdarg

LOGICAL :: exist=.FALSE.
INTEGER  (KIND=meta_ik) :: ios
CHARACTER(LEN=meta_mcl) :: lockname

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

!------------------------------------------------------------------------------
! Automatically aborts if there is no input file found on the drive
!------------------------------------------------------------------------------
lockname=TRIM(in%path)//'.'//TRIM(in%bsnm)//lock_suf

INQUIRE (FILE = TRIM(lockname), EXIST = exist)

IF((restart == 'N') .AND. (exist)) THEN
   mssg='The .*.lock file is set and a restart prohibited by default or the user.'

   INQUIRE (FILE = out%full, EXIST = exist)

   ! Delete out meta if the lock file was set.
   IF (exist) CALL execute_command_line ('rm '//TRIM(out%full))

   CALL print_err_stop(std_out, TRIM(ADJUSTL(mssg)), err=1_meta_ik)
END IF



!------------------------------------------------------------------------------
! Create a new lock file.
!------------------------------------------------------------------------------
IF(((restart == 'Y') .AND. (.NOT. exist)) .OR. ((restart == 'N') .AND. (.NOT. exist))) THEN
   CALL execute_command_line ('touch '//TRIM(lockname), CMDSTAT=ios)
   CALL print_err_stop(std_out, 'The .*.lock file could not be set.', err=ios)
END IF

IF((restart == 'Y') .AND. (exist)) CONTINUE

END SUBROUTINE meta_handle_lock_file


!------------------------------------------------------------------------------
! SUBROUTINE: meta_append
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to open a meta file to append data/ keywords
!
!> @param[inout] meta_as_rry Meta data written into a character array
!------------------------------------------------------------------------------  
SUBROUTINE meta_append(meta_as_rry)

CHARACTER(LEN=meta_mcl), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: meta_as_rry      

CALL meta_invoke(meta_as_rry)
CALL meta_continue(meta_as_rry)

END SUBROUTINE meta_append


!------------------------------------------------------------------------------
! SUBROUTINE: meta_create_new
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to create a new meta file
!
!> @param[inout] basename_requested Input basename
!------------------------------------------------------------------------------  
SUBROUTINE meta_create_new(filename_with_suffix)

CHARACTER(LEN=*), INTENT(IN) :: filename_with_suffix      

INTEGER  (KIND=meta_ik) :: ntokens
CHARACTER(LEN=meta_mcl) :: tokens(30)
LOGICAL :: exist

!------------------------------------------------------------------------------
! Automatically aborts if there is no input file found on the drive
!------------------------------------------------------------------------------
INQUIRE (FILE = TRIM(filename_with_suffix), EXIST = exist)
IF (.NOT. exist) THEN
   mssg = "The file "//TRIM(filename_with_suffix)//" does not exist."
   CALL print_err_stop(std_out, TRIM(mssg), 1)
END IF

CALL parse( str=filename_with_suffix, delims=".", args=tokens, nargs=ntokens)

!------------------------------------------------------------------------------
! Accepts any input suffix
!------------------------------------------------------------------------------
CALL parse_basename(filename_with_suffix, "."//tokens(ntokens))

!------------------------------------------------------------------------------
! Create the meta input file
!------------------------------------------------------------------------------
INQUIRE (FILE = TRIM(in%p_n_bsnm)//meta_suf, EXIST = exist)
IF (exist) THEN
   mssg = "The file "//TRIM(filename_with_suffix)//" already exists."
   CALL print_err_stop(std_out, TRIM(mssg), 1)
END IF

OPEN(UNIT=fhmeo, FILE=TRIM(in%p_n_bsnm)//meta_suf, &
   ACTION='READWRITE', ACCESS='SEQUENTIAL', STATUS='NEW')

END SUBROUTINE meta_create_new

!------------------------------------------------------------------------------
! SUBROUTINE: meta_invoke
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to open and prepare a meta file for use
!
!> @param[inout] meta_as_rry Meta data written into a character array
!------------------------------------------------------------------------------  
SUBROUTINE meta_invoke(meta_as_rry)

CHARACTER(LEN=meta_mcl), DIMENSION(:), INTENT(INOUT), ALLOCATABLE :: meta_as_rry      

! Internal Variables
INTEGER  (KIND=meta_ik) :: lines, ii, ntokens
CHARACTER(LEN=meta_mcl) :: tokens(30)
LOGICAL :: exist

!------------------------------------------------------------------------------
! Automatically aborts if there is no input file found on the drive
!------------------------------------------------------------------------------
INQUIRE (FILE = TRIM(in%full), EXIST = exist)
IF (.NOT. exist) CALL print_err_stop(std_out, "The file "//TRIM(in%full)//" does not exist.", 1)

CALL parse( str=in%full, delims=".", args=tokens, nargs=ntokens)

IF ( '.'//TRIM(tokens(ntokens)) == meta_suf) THEN
   CALL parse_basename(in%full, meta_suf)
ELSE
   ! File is not a meta file
   CALL print_err_stop(std_out, "The input file is not a *"//meta_suf//" file.", 1)
END IF

!------------------------------------------------------------------------------
! Open the meta input file
!------------------------------------------------------------------------------
OPEN(UNIT=fhmei, FILE=TRIM(in%full), ACTION='READWRITE', ACCESS='SEQUENTIAL', STATUS='OLD')

lines = count_lines(fhmei)

ALLOCATE(meta_as_rry(lines))
!------------------------------------------------------------------------------
! Read all lines into the file
!------------------------------------------------------------------------------
DO ii=1, lines
   READ(fhmei,'(A)') meta_as_rry(ii)
END DO


END SUBROUTINE meta_invoke


!------------------------------------------------------------------------------
! SUBROUTINE: meta_continue
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to invoke the output meta file
!
!> @param[inout] m_in Meta data written into a character array
!------------------------------------------------------------------------------  
SUBROUTINE meta_continue(m_in)

CHARACTER(LEN=meta_mcl), DIMENSION(:), INTENT(IN) :: m_in      

! Internal Variables
INTEGER(KIND=meta_ik) :: ios

!------------------------------------------------------------------------------
! Alter the meta file name
! The variable »alter« must be given and must be true, 
! because its a dangerous operation which may lead to data loss.
!------------------------------------------------------------------------------
CALL meta_read (fhmon, 'NEW_BSNM_FEATURE', m_in, out%features)
CALL meta_read (fhmon, 'NEW_BSNM_PURPOSE', m_in, out%purpose)

IF ((out%purpose == in%purpose) .AND. (out%features == in%features)) THEN
   WRITE(std_out,FMT_WRN) 'The basename (in part) did not change.'
END IF

!------------------------------------------------------------------------------
! Build the new outfile path
!------------------------------------------------------------------------------
! Nomenclature: dataset_type_purpose_app_features
! This assignment requres the out = in assignment before
out%bsnm =     TRIM(out%dataset)//&
          '_'//TRIM(out%type)//&
          '_'//TRIM(out%purpose)//&
          '_'//TRIM(global_meta_prgrm_mstr_app)//&
          '_'//TRIM(out%features)

out%p_n_bsnm = TRIM(out%path)//&
               TRIM(out%bsnm)

out%full = TRIM(out%p_n_bsnm)//meta_suf

!------------------------------------------------------------------------------
! System call to update the file name of the meta file
!------------------------------------------------------------------------------
CALL execute_command_line ('cp '//TRIM(in%full)//' '//TRIM(out%full), CMDSTAT=ios)
CALL print_err_stop(std_out, 'The update of the meta filename went wrong.', ios)

!------------------------------------------------------------------------------
! Open the meta output file
!------------------------------------------------------------------------------
OPEN(UNIT=fhmeo, FILE=TRIM(out%full), ACTION='WRITE', ACCESS='APPEND', STATUS='OLD')

WRITE(fhmeo, '(A)')

END SUBROUTINE meta_continue




!------------------------------------------------------------------------------
! SUBROUTINE: meta_start_ascii
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to deal with the logging and renaming of the log file in context 
!> of the meta data approach. 
!
!> @Description
!> Meta log only gets called __after__ meta append or meta_close respectively. 
!> This way, a log file is optional.
!
!> In the meta file format, the restart procedure only is checked for the
!> meta file itself.
!>
!> If the variable restart is not explicitly set .TRUE., the program will not 
!> restart.
!> If the variable in/out are not set, the program will not start/stop 
!> accordingly.
!
!> The Name of the temporary logfile is hardcoded!
!> »'.temporary.'//log_suf«
!
!> @param[in] fh File handle of the input
!> @param[in] suf Suffix of the file
!> @param[in] restart Logfiles (temporary and permanent)
!------------------------------------------------------------------------------  
SUBROUTINE meta_start_ascii(fh, suf)

INTEGER(KIND=meta_ik), INTENT(IN) :: fh
CHARACTER(LEN=*), INTENT(IN) :: suf

CHARACTER(LEN=meta_mcl) :: temp_f_suf, perm_f_suf
INTEGER(KIND=meta_ik) :: ios

LOGICAL :: exist_temp, exist_perm

! The temporaray file is a hidden one.
temp_f_suf = TRIM(out%path)//'.temporary'//TRIM(suf)
perm_f_suf = TRIM(out%p_n_bsnm)//TRIM(suf)
  
!------------------------------------------------------------------------------
! Check for a temporary file
! Check for a permanent file
!------------------------------------------------------------------------------
INQUIRE (FILE = temp_f_suf, EXIST = exist_temp)
INQUIRE (FILE = out%p_n_bsnm//TRIM(suf), EXIST = exist_perm)

!------------------------------------------------------------------------------
! Check whether file needs to be deleted.
!------------------------------------------------------------------------------
IF(exist_temp) THEN
   CALL execute_command_line ('rm -r '//TRIM(temp_f_suf), CMDSTAT=ios)   
   CALL print_err_stop(std_out, '»'//TRIM(temp_f_suf)//'« not deletable.',ios)
END IF

IF(exist_perm) THEN
   CALL execute_command_line ('rm -r '//TRIM(out%p_n_bsnm)//TRIM(suf), CMDSTAT=ios)
   CALL print_err_stop(std_out, '»'//TRIM(out%full)//'« not deletable.', ios)
END IF

OPEN(UNIT=fh, FILE=TRIM(temp_f_suf), ACTION='WRITE', ACCESS='SEQUENTIAL', STATUS='NEW')

END SUBROUTINE meta_start_ascii

!------------------------------------------------------------------------------
! SUBROUTINE: meta_stop_ascii
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to stop the logging and renaming of additional ascii files.
!
!> @param[in] fh File handle of the input
!> @param[in] suf Suffix of the file
!------------------------------------------------------------------------------  
SUBROUTINE meta_stop_ascii(fh, suf)

INTEGER(KIND=meta_ik), INTENT(IN) :: fh
CHARACTER(LEN=*), INTENT(IN) :: suf

CHARACTER(LEN=meta_mcl) :: temp_f_suf, perm_f_suf
INTEGER  (KIND=meta_ik) :: ios

temp_f_suf = TRIM(out%path)//'.temporary'//TRIM(suf)
perm_f_suf = TRIM(out%p_n_bsnm)//TRIM(suf)

CLOSE (fh)

!------------------------------------------------------------------------------
! The temporary log file must be renamed to a permanent one
!------------------------------------------------------------------------------
CALL execute_command_line ('mv '//TRIM(temp_f_suf)//' '//TRIM(out%p_n_bsnm)//TRIM(suf), CMDSTAT=ios)

IF(ios /= 0_meta_ik) THEN
   mssg='Can not rename the suffix_file from »'//TRIM(temp_f_suf)//'« to the proper basename.'
   CALL print_err_stop(std_out, mssg, 0)
END IF

END SUBROUTINE meta_stop_ascii



!------------------------------------------------------------------------------
! SUBROUTINE: count_lines
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Truncate a keyword which was too long. Could do other stuff as well.
!
!> @param[in] fh File handle 
!> @param[in] keyword Keyword to check
!------------------------------------------------------------------------------  
function count_lines(un) result(no_lines)

Integer, Intent(in) :: un
Integer(kind=ik)    :: no_lines

Integer :: io_stat
Character(len=2) :: temp_char

io_stat = 0
no_lines=0

Rewind(un)

Do While (io_stat == 0)

Read(un,'(A)', End=1000, iostat=io_stat) temp_char
no_lines = no_lines + 1

End Do

1000 Continue

Rewind(un)

End function count_lines

!------------------------------------------------------------------------------
! SUBROUTINE: check_keyword
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Truncate a keyword which was too long. Could do other stuff as well.
!
!> @param[in] fh File handle 
!> @param[in] keyword Keyword to check
!------------------------------------------------------------------------------  
SUBROUTINE check_keyword(fh, keyword)

INTEGER(KIND=meta_ik) :: fh 
CHARACTER(LEN=*)   :: keyword
CHARACTER(LEN=kcl) :: kywd_lngth

kywd_lngth = ''

IF(LEN_TRIM(keyword) .GT. LEN(kywd_lngth)) THEN

   WRITE(fh, '(A)') ''
   
   WRITE(std_out,FMT_WRN) "The keyword »"//TRIM(keyword)//"« is longer"
   WRITE(std_out,FMT_WRN) "than the convention allows and therefore truncated!"
   
   kywd_lngth = keyword(1:LEN(kywd_lngth))
ELSE
   kywd_lngth = keyword
END IF

END SUBROUTINE check_keyword


!------------------------------------------------------------------------------
! SUBROUTINE: parse_basename
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Parse the basename
!
!> @param[in] filename Full name of the file
!> @param[in] suf Expected suffix
!------------------------------------------------------------------------------  
SUBROUTINE parse_basename(filename, suf)

CHARACTER(LEN=*) :: filename, suf
INTEGER  (KIND=meta_ik) :: ntokens
CHARACTER(LEN=meta_mcl) :: tokens(30)

in%full = TRIM(ADJUSTL(filename))

!------------------------------------------------------------------------------
! Parse all basename and path details.
!------------------------------------------------------------------------------
in%p_n_bsnm = in%full(1:LEN_TRIM(in%full)-LEN_TRIM(TRIM(ADJUSTL(suf)))) 

CALL parse( str=TRIM(in%p_n_bsnm), delims="/", args=tokens, nargs=ntokens)

in%path = in%p_n_bsnm(1:LEN_TRIM(in%p_n_bsnm) - LEN_TRIM(tokens(ntokens)))     
in%bsnm = TRIM(tokens(ntokens))

CALL parse( str=TRIM(in%bsnm), delims="_", args=tokens, nargs=ntokens)

IF (ntokens /= 5) THEN
   mssg = "Invalid basename given: "//TRIM(in%p_n_bsnm)
   CALL print_err_stop(std_out, mssg, 1)
END IF

in%dataset = TRIM(tokens(1))
in%type    = TRIM(tokens(2))
in%purpose = TRIM(tokens(3))
in%app      = TRIM(tokens(4))
in%features = TRIM(tokens(5))

out = in  

END SUBROUTINE parse_basename

!------------------------------------------------------------------------------
! SUBROUTINE: check_unit
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Routine to truncate a keyword which was too long. Could do other stuff 
!> as well.
!
!> @param[in] fh File handle 
!> @param[in] keyword Keyword to check
!------------------------------------------------------------------------------  
SUBROUTINE check_unit(fh, unit)

INTEGER  (KIND=meta_ik) :: fh 
CHARACTER(LEN=*)   :: unit
CHARACTER(LEN=ucl) :: unit_lngth

! Check unit length for convention and proper formatting
IF(LEN_TRIM(unit) .GT. LEN(unit_lngth)) THEN

   WRITE(fh, '(A)') ''

   WRITE(std_out,FMT_WRN) "The unit "//TRIM(unit)//" is longer than"
   WRITE(std_out,FMT_WRN) "the convention allows and therefore truncated!"

   unit_lngth = unit(1:LEN(unit_lngth))
ELSE
   unit_lngth = unit
END IF

END SUBROUTINE check_unit


!------------------------------------------------------------------------------
! SUBROUTINE: meta_extract_keyword_data
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to extract the data string of keywords. 
!
!> @Description
!> Module to parse information of keywords. 
!> An arbitrary Keyword with up to »kcl« characters may be specified.
!> The program reads keywords as long as they are before or withing the owns 
!> programs scope.
! 
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] dims Dimensions requested
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars Datatype to read in
!------------------------------------------------------------------------------
SUBROUTINE meta_extract_keyword_data (fh, keyword, dims, m_in, res_tokens, res_ntokens)
   
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
INTEGER(KIND=meta_ik), INTENT(IN) :: dims
CHARACTER(LEN=meta_mcl), DIMENSION(:), INTENT(IN) :: m_in
CHARACTER(LEN=meta_mcl) :: res_tokens(30)
INTEGER(KIND=meta_ik) :: res_ntokens

! Internal variables
INTEGER(KIND=meta_ik) :: kywd_found, ii, ntokens
CHARACTER(LEN=meta_mcl) :: tokens(30)
LOGICAL :: override

kywd_found = 0
override = .FALSE.

CALL check_keyword(fh, keyword)

!------------------------------------------------------------------------------
! Parse Data out of the input array
!------------------------------------------------------------------------------
DO ii =1, SIZE(m_in) 
   CALL parse(str=m_in(ii), delims=' ', args=tokens, nargs=ntokens)

   SELECT CASE(tokens(1))
      CASE('*', 'd', 'r', 'w')

         IF (tokens(2) == TRIM(keyword)) THEN
            kywd_found = 1

            !------------------------------------------------------------------------------
            ! Store the keywords data.
            ! Following m_in(ii) - lines -  will overwrite this information.
            !------------------------------------------------------------------------------
            res_tokens = tokens
            res_ntokens = ntokens

            !------------------------------------------------------------------------------
            ! Exit, if the keyword appears the first time in the programs scope.
            !------------------------------------------------------------------------------
            IF ((override) .AND. (tokens(1) /= 'w')) GOTO 2011

         END IF
      CASE('p')
         !------------------------------------------------------------------------------
         ! User tells, that the scope of another program begins.
         !------------------------------------------------------------------------------
         IF ((.NOT. override) .AND. (tokens(2) == TRIM(global_meta_program_keyword))) THEN
            override = .TRUE.
         ELSE IF (override) THEN
            !------------------------------------------------------------------------------
            ! Leave, if the scope of the next program begins.
            !------------------------------------------------------------------------------
            GOTO 2011
         END IF

   END SELECT
END DO

2011 CONTINUE

IF((res_ntokens < dims+2) .AND. (kywd_found /= 0)) THEN
   CALL print_err_stop(std_out, "Data of keyword '"//TRIM(ADJUSTL(keyword))//"' invalid.", 1)
END IF

IF(kywd_found == 0) THEN
   mssg = "Keyword '"//TRIM(ADJUSTL(keyword))//"' not found in the meta file!"
   CALL print_err_stop(std_out, TRIM(ADJUSTL(mssg)), 1)
END IF

END SUBROUTINE meta_extract_keyword_data


!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_C
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Wrapper to parse character Keywords.
!
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] chars Datatype to read in
!------------------------------------------------------------------------------
SUBROUTINE meta_read_C (fh, keyword, m_in, chars)
   
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=meta_mcl), DIMENSION(:), INTENT(IN)  :: m_in      
CHARACTER(LEN=*), INTENT(OUT) :: chars 

! Internal variables
CHARACTER(LEN=meta_mcl) :: tokens(30)
INTEGER(KIND=meta_ik) :: ntokens

CALL meta_extract_keyword_data (fh, keyword, 1, m_in, tokens, ntokens)

chars = TRIM(ADJUSTL(tokens(3)))

END SUBROUTINE meta_read_C

!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_I0D
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Wrapper to parse Keywords with 0D integer data.
! 
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] int_0D Datatype to read in
!------------------------------------------------------------------------------
SUBROUTINE meta_read_I0D (fh, keyword, m_in, int_0D)
     
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=meta_mcl), DIMENSION(:), INTENT(IN) :: m_in      
INTEGER(KIND=meta_ik), INTENT(OUT) :: int_0D 

! Internal variables
CHARACTER(LEN=meta_mcl) :: tokens(30)
INTEGER(KIND=meta_ik) :: ntokens

CALL meta_extract_keyword_data (fh, keyword, 1, m_in, tokens, ntokens)

READ(tokens(3), '(I12)') int_0D 

END SUBROUTINE meta_read_I0D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_R0D
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Wrapper to parse Keywords with 0D floating point data.
! 
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] real_0D Datatype to read in
!------------------------------------------------------------------------------
SUBROUTINE meta_read_R0D (fh, keyword, m_in, real_0D)
     
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=meta_mcl), DIMENSION(:), INTENT(IN) :: m_in      
REAL(KIND=meta_rk), INTENT(OUT) :: real_0D 

! Internal variables
CHARACTER(LEN=meta_mcl) :: tokens(30)
INTEGER(KIND=meta_ik) :: ntokens

CALL meta_extract_keyword_data (fh, keyword, 1, m_in, tokens, ntokens)

READ(tokens(3), '(F39.10)') real_0D 

END SUBROUTINE meta_read_R0D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_I1D
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Wrapper to parse Keywords with 1D integer data.
! 
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] int_1D Datatype to read in
!------------------------------------------------------------------------------
SUBROUTINE meta_read_I1D (fh, keyword, m_in, int_1D)

INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=meta_mcl), DIMENSION(:), INTENT(IN)  :: m_in      
INTEGER(KIND=meta_ik), DIMENSION(:), INTENT(OUT) :: int_1D 

! Internal variables
CHARACTER(LEN=meta_mcl) :: tokens(30)
INTEGER(KIND=meta_ik) :: ntokens

CALL meta_extract_keyword_data (fh, keyword, SIZE(int_1D), m_in, tokens, ntokens)

READ(tokens(3:2+SIZE(int_1D)), '(I12)') int_1D

END SUBROUTINE meta_read_I1D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_read_R1D
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Wrapper to parse Keywords with 1D integer data.
! 
!> @param[in] fh File handle to read a keyword from.
!> @param[in] keyword Keyword to read
!> @param[in] m_in Array of lines of ascii meta file
!> @param[in] real_1D Datatype to read in
!------------------------------------------------------------------------------
SUBROUTINE meta_read_R1D (fh, keyword, m_in, real_1D)

INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=meta_mcl), DIMENSION(:), INTENT(IN) :: m_in      
REAL(KIND=meta_rk), DIMENSION(:), INTENT(OUT) :: real_1D 

! Internal variables
CHARACTER(LEN=meta_mcl) :: tokens(30)
INTEGER(KIND=meta_ik) :: ntokens

CALL meta_extract_keyword_data (fh, keyword, SIZE(real_1D), m_in, tokens, ntokens)

READ(tokens(3:2+SIZE(real_1D)), '(F39.10)') real_1D

END SUBROUTINE meta_read_R1D


!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_keyword
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write finalized strings of keywords
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] stdspcfill String with data
!> @param[in] unit Unit of the value
!------------------------------------------------------------------------------
SUBROUTINE meta_write_keyword (fh, keyword, stdspcfill, unit)
   
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: stdspcfill 
CHARACTER(LEN=*), INTENT(IN) :: unit

INTEGER(KIND=ik) :: rpt
CHARACTER(LEN=meta_scl) :: fmt, str
CHARACTER(LEN=8)  :: date
CHARACTER(LEN=10) :: time
CHARACTER(LEN=5)  :: timezone

CALL check_keyword(fh, keyword)
CALL check_unit(fh, unit)

WRITE(fmt, '(A,I0,A)') "(2A, T", kcl, ")"
WRITE(fh, fmt, ADVANCE='NO') "w ", keyword

WRITE(fmt, '(A,I0,A)') "(A, T", stdspc+1, ")"
WRITE(fh, fmt, ADVANCE='NO') TRIM(ADJUSTL(stdspcfill))

!------------------------------------------------------------------------------
! Only write if stdspcfill (are of actual information/data) was not overflowed
! < instead of <= to get 1 space clearance
!------------------------------------------------------------------------------  
IF(LEN_TRIM(ADJUSTL(stdspcfill)) <  stdspc) THEN
   WRITE(fmt, '(A,I0,A)') "(A, T", ucl+1, ")"
   WRITE(fh, fmt, ADVANCE='NO') unit
   rpt=0
ELSE
   rpt = stdspc+ucl-LEN_TRIM(ADJUSTL(stdspcfill))
END IF

! Same as comment before
IF(LEN_TRIM(ADJUSTL(stdspcfill)) <  stdspc+ucl) THEN
   CALL DATE_AND_TIME(DATE=date, TIME=time, ZONE=timezone)

   str = '' ! Clear string
   str = REPEAT(' ', rpt)//date(7:8)//'.'//date(5:6)//'.'//date(1:4)
   str = TRIM(str)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)
   str = TRIM(str)//' '//timezone

   WRITE(fh, '(A)') TRIM(str)
ELSE
   WRITE(fh, '(A)') "" ! Basically a newline
END IF
END SUBROUTINE meta_write_keyword



!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_sha256sum
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write finalized strings of keywords
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] stdspcfill String with data
!> @param[in] unit Unit of the value
!------------------------------------------------------------------------------
SUBROUTINE meta_write_sha256sum(binary_name)
   
CHARACTER(LEN=*), INTENT(IN) :: binary_name

CHARACTER(LEN=kcl-1) :: keyword = ''
CHARACTER(LEN=meta_scl) :: fmt, stdspcfill
INTEGER(KIND=meta_ik), DIMENSION(5) :: stat = 0
INTEGER(KIND=meta_ik) :: ios

LOGICAL :: exist

!------------------------------------------------------------------------------
! Write "Keyword"
!------------------------------------------------------------------------------
keyword = "w SHA256SUM_OF_BINARY"

WRITE(fmt, '(A,I0,A)') "(2A, T", kcl, ")"
WRITE(fhmeo, fmt, ADVANCE='NO') keyword

!------------------------------------------------------------------------------
! Check the buffer file
!------------------------------------------------------------------------------
INQUIRE(FILE = 'temp_buffer', EXIST = exist)

IF (exist) THEN
   CALL EXECUTE_COMMAND_LINE ('rm -r temp_buffer', CMDSTAT=ios)   

   IF(ios /= 0_meta_ik) THEN
      mssg='Can not delete the temp_buffer'
      CALL print_err_stop(std_out, mssg, 0)
      stat(1) = 1      
   END IF
END IF

!------------------------------------------------------------------------------
! Check for auxiliary programs
!------------------------------------------------------------------------------
CALL EXECUTE_COMMAND_LINE("which cut > /dev/null 2> /dev/null", CMDSTAT=stat(2))
CALL EXECUTE_COMMAND_LINE("which sha256sum > /dev/null 2> /dev/null", CMDSTAT=stat(3))

!------------------------------------------------------------------------------
! Deal with the buffer file
!------------------------------------------------------------------------------
IF(SUM(stat)==0) THEN
   OPEN(UNIT=9, FILE='temp_buffer', ACTION='READWRITE', STATUS='NEW')

   CALL EXECUTE_COMMAND_LINE("sha256sum "//TRIM(ADJUSTL(binary_name))//" | cut -d ' ' -f 1 >> 'temp_buffer'", CMDSTAT=stat(4))

   READ(9, '(A)', iostat=stat(5)) stdspcfill

   CLOSE(9)
END IF

IF (SUM(stat) == 0) THEN
   WRITE(fhmeo, fmt) TRIM(ADJUSTL(stdspcfill))
ELSE
   WRITE(fhmeo, fmt) "Could not get sha256sum. One of the previous system calls failed."
END IF

INQUIRE(FILE='temp_buffer', EXIST=exist)

IF (exist) CALL EXECUTE_COMMAND_LINE ('rm -r temp_buffer', CMDSTAT=ios)   

END SUBROUTINE meta_write_sha256sum


!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_C
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords with character output. 
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] stdspcfill Characters to write
!------------------------------------------------------------------------------
SUBROUTINE meta_write_C (fh, keyword, stdspcfill)
   
INTEGER  (KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: stdspcfill 

CALL meta_write_keyword (fh, keyword, stdspcfill, '')

END SUBROUTINE meta_write_C

!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_I0D
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords of type integer dim 0.
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] unit Unit of the value
!> @param[in] int_0D Datatype to read in
!------------------------------------------------------------------------------
SUBROUTINE meta_write_I0D (fh, keyword, unit, int_0D)
   
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: unit
INTEGER(KIND=meta_ik), INTENT(IN) :: int_0D 

CHARACTER(LEN=meta_scl) :: stdspcfill

WRITE(stdspcfill, '(I0)') int_0D

CALL meta_write_keyword (fh, keyword, stdspcfill, unit)

END SUBROUTINE meta_write_I0D


!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_R0D
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords of type Real dim 0.
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] unit Unit of the value
!> @param[in] real_0D Datatype to read in
!------------------------------------------------------------------------------
SUBROUTINE meta_write_R0D (fh, keyword, unit, real_0D)
   
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: unit
REAL(KIND=meta_ik), INTENT(IN) :: real_0D 

CHARACTER(LEN=meta_scl) :: stdspcfill, fmt

WRITE(fmt, '(A,I0,A)') "(F", stdspc, ".7)"
WRITE(stdspcfill, fmt) real_0D ! '(F30.7)'

CALL trimzero(stdspcfill)

CALL meta_write_keyword (fh, keyword, stdspcfill, unit)

END SUBROUTINE meta_write_R0D


!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_I1D
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords of type integer dim 1.
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] unit Unit of the value
!> @param[in] int_0D Datatype
!------------------------------------------------------------------------------
SUBROUTINE meta_write_I1D (fh, keyword, unit, int_1D)
   
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: unit
INTEGER(KIND=meta_ik), INTENT(IN), DIMENSION(:) :: int_1D 

CHARACTER(LEN=meta_scl) :: stdspcfill, str, fmt
INTEGER  (KIND=meta_ik) :: ii

stdspcfill = ''
str = ''

DO ii=1, SIZE(int_1D)
   str = ''
   WRITE(fmt, '(A,I0,A)') "(I", (stdspc/3)-1, ")"
   WRITE(str, fmt) int_1D(ii) ! '(I0)'
   stdspcfill = TRIM(stdspcfill)//' '//TRIM(str)
END DO

CALL meta_write_keyword (fh, keyword, stdspcfill, unit)

END SUBROUTINE meta_write_I1D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_write_R1D
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Module to write keywords of type Real dim 1.
!
!> @param[in] fh File handle to write a log/mon or text to.
!> @param[in] keyword Keyword to write
!> @param[in] unit Unit of the value
!> @param[in] real_1D Datatype
!------------------------------------------------------------------------------
SUBROUTINE meta_write_R1D (fh, keyword, unit, real_1D)
   
INTEGER(KIND=meta_ik), INTENT(IN) :: fh 
CHARACTER(LEN=*), INTENT(IN) :: keyword
CHARACTER(LEN=*), INTENT(IN) :: unit
REAL(KIND=meta_rk), INTENT(IN), DIMENSION(:) :: real_1D 

CHARACTER(LEN=meta_scl) :: stdspcfill, str, fmt
INTEGER  (KIND=meta_ik) :: ii

stdspcfill = ''
str = ''

DO ii=1, SIZE(real_1D)
   str = ''
   WRITE(fmt, '(A,I0,A)') "(F", (stdspc/3)-1, ".7)"
   WRITE(str, fmt) real_1D(ii)
   ! CALL trimzero(str) ! Choose preferred formatting

   stdspcfill = TRIM(stdspcfill)//' '//TRIM(str)
END DO

CALL meta_write_keyword (fh, keyword, stdspcfill, unit)

END SUBROUTINE meta_write_R1D

!------------------------------------------------------------------------------
! SUBROUTINE: meta_signing
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to close a meta file.
!
!> @description
!> Requires a "revision.meta" or similar inclusion of verisoning info, 
!> provided by a makefile. Furhermore, it requires a global_stds file.
!> Some of the variables are retrieved from global_std module.
!
!> @param[in] binary Name of the executable
!------------------------------------------------------------------------------
SUBROUTINE meta_signing(binary)

CHARACTER(LEN=*), INTENT(IN) :: binary

WRITE(fhmeo, '(A)')
CALL meta_write(fhmeo, 'PROGRAM_VERSION', revision)
CALL meta_write(fhmeo, 'PROGRAM_GIT_HASH', hash)

CALL meta_write_sha256sum(binary)

CALL meta_write(fhmeo, 'COMPUTATION_FINISHED' , 'Succesfully')

END SUBROUTINE meta_signing


!------------------------------------------------------------------------------
! SUBROUTINE: meta_close
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to close a meta file.
!
!> @description
!> provided by a makefile. Furhermore, it requires a global_stds file.
!------------------------------------------------------------------------------
SUBROUTINE meta_close()

LOGICAL :: opened

WRITE(fhmeo, '(A)')
WRITE(fhmeo, "(100('-'))")

!------------------------------------------------------------------------------
! Check and close files - Routine: (fh, filename, abrt, stat)
!------------------------------------------------------------------------------
INQUIRE(UNIT=fhmei, OPENED=opened)
IF(opened) CLOSE (fhmei)

INQUIRE(UNIT=fhmeo, OPENED=opened)
IF(opened) CLOSE (fhmeo)

 
END SUBROUTINE meta_close


!------------------------------------------------------------------------------
! SUBROUTINE: meta_delete_empty_file
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Check filesize and delete if 0.
!
!> @param(in) filename Path and name of the file to be checked and deleted.
!------------------------------------------------------------------------------
SUBROUTINE meta_delete_empty_file(filename)

CHARACTER(LEN=*), INTENT(IN) :: filename

INTEGER(KIND=meta_ik) :: filesize = 0, ios

!------------------------------------------------------------------------------
! Check and close files - Routine: (fh, filename, abrt, stat)
!------------------------------------------------------------------------------
INQUIRE(FILE=TRIM(ADJUSTL(filename)), SIZE=filesize)

IF(filesize == 0) THEN
   CALL execute_command_line ('rm -r '//TRIM(ADJUSTL(filename)), CMDSTAT=ios)
END IF
 
END SUBROUTINE meta_delete_empty_file

END MODULE meta


!------------------------------------------------------------------------------
! MODULE: meta_puredat_interface
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
! @Description:
!> Module to convert a *.raw/*.meta data description into a PureDat description
!
! REVISION HISTORY:
! 27 11 2021 - Initial version
!------------------------------------------------------------------------------
MODULE meta_puredat_interface

   USE meta
   USE user_interaction

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: convert_meta_to_puredat
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> One-off special purpose interface to convert the meta to the PureDat format
!> The routine requires an opened meta file!
!
!> @param[in] free_file_handle File handle for use in this routine
!> @param[in] p_n_bsnm Path and basename of the *.meta and the *.raw file
!------------------------------------------------------------------------------
SUBROUTINE convert_meta_to_puredat(free_file_handle, m_rry)

INTEGER(KIND=meta_ik), INTENT(IN) :: free_file_handle
CHARACTER(LEN=meta_mcl), DIMENSION(:), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: m_rry      

CHARACTER(LEN=meta_mcl) :: suf, datatype, branch_description, field_content_desc, osagcs

REAL(KIND=meta_rk), DIMENSION(3) :: grid_spacings, origin_shift

! INTEGER(KIND=meta_ik), DIMENSION(7) :: bytesizes = [ 1, 2, 4, 8, 8, 1, 1]
INTEGER(KIND=meta_ik), DIMENSION(7,3) :: stda ! Stream data, 7 streams; no_of_data, lb, ub
INTEGER(KIND=meta_ik), DIMENSION(3) :: vox_per_dim, origin
INTEGER(KIND=meta_ik) :: stdout, rawsize, rawdata, ii, stat, my_size, my_pos

LOGICAL :: opened, fex

!------------------------------------------------------------------------------
! PureDat/Meta Module Formatters
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: PDM_arrowsA  = "('<',A,'> ', A)"
CHARACTER(Len=*), PARAMETER :: PDM_arrowsL  = "('<',A,'> ', L)"
CHARACTER(Len=*), PARAMETER :: PDM_arrowsI0 = "('<',A,'> ', I0)"
CHARACTER(Len=*), PARAMETER :: PDM_branch  = "('<==branch==>')"

!------------------------------------------------------------------------------
! Initialize variables
!------------------------------------------------------------------------------
stdout = 6_meta_ik
stda = 0_meta_ik
opened = .FALSE.
fex = .FALSE.

!------------------------------------------------------------------------------
! Check whether the meta file is opened and establish the proper status
!------------------------------------------------------------------------------
INQUIRE(UNIT=fhmei, OPENED=opened)

IF(.NOT. opened) THEN
   CALL print_err_stop(stdout, TRIM(in%full)//" not opened. Check your implementation!", 1)
END IF

!------------------------------------------------------------------------------
! Gather the size of the resulting stream
!------------------------------------------------------------------------------
INQUIRE(FILE=TRIM(in%p_n_bsnm)//raw_suf, EXIST=fex, SIZE=rawsize)
IF(.NOT. fex) CALL print_err_stop(stdout, TRIM(in%p_n_bsnm)//raw_suf//" does not exist.", 1)

!------------------------------------------------------------------------------
! Read all required keywords to write the information ino the PureDat format
! stdout will be given to redirect errors to the command line or error file
!------------------------------------------------------------------------------
CALL meta_read (stdout, 'CT_SCAN', m_rry, branch_description)

branch_description = TRIM(branch_description)//" scalar data structure"
 
CALL meta_read (stdout, 'INTERNAL_ID', m_rry, field_content_desc)
CALL meta_read (stdout, 'DIMENSIONS', m_rry, vox_per_dim)
CALL meta_read (stdout, 'SPACING', m_rry, grid_spacings)
CALL meta_read (stdout, 'ORIGIN', m_rry, origin)
CALL meta_read (stdout, 'ORIGIN_SHIFT_GLBL', m_rry, origin_shift)
CALL meta_read (stdout, 'TYPE_RAW', m_rry, datatype)

!------------------------------------------------------------------------------
! Rename the raw file and enter the stda
!------------------------------------------------------------------------------
suf=''
rawdata=0
SELECT CASE(TRIM(datatype))
   CASE('ik1')
      suf = ".int1.st"
      rawdata = 1
      stda(rawdata,:) = [rawsize, 1_meta_ik, rawsize] 
   CASE('ik2')
      suf = ".int2.st"
      rawdata = 2 
      stda(rawdata,:) = [rawsize/2, 1_meta_ik, rawsize/2] 
   CASE('ik4')
      suf = ".int4.st"
      rawdata = 3 
      stda(rawdata,:) = [rawsize/4, 1_meta_ik, rawsize/4] 
   CASE('ik8')
      suf = ".int8.st"
      rawdata = 4 
      stda(rawdata,:) = [rawsize/8, 1_meta_ik, rawsize/8] 
   CASE('rk8')
      suf = ".real8.st"
      rawdata = 5 
      stda(rawdata,:) = [rawsize/8, 1_meta_ik, rawsize/8] 
END SELECT

IF(suf /= '') THEN
   INQUIRE(FILE=TRIM(out%p_n_bsnm)//TRIM(suf), EXIST=fex)
   IF(fex) CALL print_err_stop(stdout, TRIM(out%p_n_bsnm)//TRIM(suf)//" already exists.", 1)

   CALL EXECUTE_COMMAND_LINE &
      ("mv "//TRIM(in%p_n_bsnm)//".raw "//TRIM(in%p_n_bsnm)//TRIM(suf), CMDSTAT=stat)

   IF(stat /= 0) CALL print_err_stop(stdout, &
      "Renaming "//TRIM(in%p_n_bsnm)//".raw to "//TRIM(in%p_n_bsnm)//TRIM(suf)//" failed.", 1)
ELSE
   CALL print_err_stop(stdout, 'No datatype given to convert meta/raw to PureDat.', 1)
END IF

!------------------------------------------------------------------------------
! Invoke the header file
!------------------------------------------------------------------------------
INQUIRE(FILE=TRIM(in%p_n_bsnm)//head_suf, EXIST=fex)
IF(fex) CALL print_err_stop(stdout, TRIM(in%p_n_bsnm)//head_suf//" already exists.", 1)

OPEN(UNIT=fhh, FILE=TRIM(in%p_n_bsnm)//head_suf, &
   ACTION='WRITE', ACCESS='SEQUENTIAL', STATUS='NEW')

!------------------------------------------------------------------------------
! Write the header header :-) 
! This stuff is hardcoded and not flexible yet.
!------------------------------------------------------------------------------
WRITE(fhh, PDM_branch)
WRITE(fhh, PDM_arrowsA) "description", "'"//TRIM(branch_description)//"'"
WRITE(fhh, PDM_arrowsI0) "no_of_branches", 0
WRITE(fhh, PDM_arrowsI0) "no_of_leaves", 6 
WRITE(fhh, PDM_arrowsL) "streams_allocated", .TRUE.  
WRITE(fhh, PDM_arrowsI0) "size_int1_stream", stda(1,1)
WRITE(fhh, PDM_arrowsI0) "size_int2_stream", stda(2,1)
WRITE(fhh, PDM_arrowsI0) "size_int4_stream", stda(3,1) + 6 ! origin, d)
WRITE(fhh, PDM_arrowsI0) "size_int8_stream", stda(4,1)
WRITE(fhh, PDM_arrowsI0) "size_real_stream", stda(5,1) + 6 ! origin_shift, sp)
WRITE(fhh, PDM_arrowsI0) "size_char_stream", stda(6,1) + LEN_TRIM(field_content_desc)
WRITE(fhh, PDM_arrowsI0) "size_log_stream" , stda(7,1)

!------------------------------------------------------------------------------
! Write the meta data into the stream files
! Only int4, real8 and chars do get additional data.
!------------------------------------------------------------------------------
DO ii=1, 6 
   ! Not integrated with 2nd SELECT CASE(ii) to get a proper CONTINUE and suf
   SELECT CASE(ii)
      CASE(1); suf = ".int1.st"
      CASE(2); suf = ".int2.st"
      CASE(3); suf = ".int4.st"
      CASE(4); suf = ".int8.st"
      CASE(5); suf = ".real8.st"
      CASE(6); suf = ".char.st"
   END SELECT

   IF(ii == rawdata) THEN

      INQUIRE(FILE=TRIM(in%p_n_bsnm)//TRIM(suf), SIZE=my_size)
      

      ! IF ((my_size/PRODUCT(vox_per_dim)) /= bytesizes(ii)) THEN
      !    WRITE(std_out,'(A,I0)') "Size of scalar data input: ", my_size
      !    WRITE(std_out,'(A,I0)') "Amount of points in image: ", PRODUCT(vox_per_dim)
      !    WRITE(std_out,'(A,I0)') "Expected size per entry: ", bytesizes(ii), " Byte(s)"
      !    WRITE(std_out,'(2A)') "Filename: ", TRIM(in%p_n_bsnm)//TRIM(suf)

      !    mssg = ""
      !    mssg = "Filesize of the raw data does not match the dimensions and data type.&
      !       & Please check the raw data or the fortran source code of the meta/raw converter.&
      !       & Converter assumes that only one blob of data set is containted in the raw file!&
      !       & Converter cannot deal with multiple binary images within 1 raw file at the moment."
      !    CALL print_err_stop(stdout, mssg, 1) 
      ! END IF

      ! It is the file of the original *.raw data
      OPEN(UNIT=free_file_handle, FILE=TRIM(in%p_n_bsnm)//TRIM(suf), &
         ACCESS='STREAM', STATUS='OLD')


   ELSE IF (ii == 6) THEN
      ! It is the character stream file
      INQUIRE(FILE=TRIM(in%p_n_bsnm)//TRIM(suf), EXIST=fex)
      IF(fex) CALL print_err_stop(stdout, TRIM(in%p_n_bsnm)//TRIM(suf)//" already exists.", 1)

      OPEN(UNIT=free_file_handle, FILE=TRIM(in%p_n_bsnm)//TRIM(suf), &
         ACCESS='STREAM', STATUS='NEW')
   ELSE
      ! It is another stream file
      INQUIRE(FILE=TRIM(in%p_n_bsnm)//TRIM(suf), EXIST=fex)
      IF(fex) CALL print_err_stop(stdout, TRIM(in%p_n_bsnm)//TRIM(suf)//" already exists.", 1)

      OPEN(UNIT=free_file_handle, FILE=TRIM(in%p_n_bsnm)//TRIM(suf), &
         ACCESS='STREAM', STATUS='NEW')
   END IF

   INQUIRE(UNIT=free_file_handle, SIZE=my_size)

   my_pos = my_size+1

   SELECT CASE(ii)
      CASE(3)
         CALL write_leaf_to_header(fhh, "Scalar data", ii, stda(ii,:))
 
         stda(ii,:) = [3, stda(ii,3)+1, stda(ii,3)+3]
         CALL write_leaf_to_header(fhh, "Number of voxels per direction", ii, stda(ii,:))           
         WRITE(free_file_handle, POS=my_pos) INT(vox_per_dim, KIND=4)

         ! stda relative to values before (!)
         stda(ii,:) = [3, stda(ii,3)+1, stda(ii,3)+3]
         CALL write_leaf_to_header(fhh, "Origin", ii, stda(ii,:))
         WRITE(free_file_handle, POS=my_pos+4*3) INT(origin, KIND=4)

      CASE(5)
         stda(ii,:) = [3, stda(ii,3)+1, stda(ii,3)+3]
         CALL write_leaf_to_header(fhh, "Grid spacings", ii, stda(ii,:))
         WRITE(free_file_handle, POS=my_pos) grid_spacings

         stda(ii,:) = [3, stda(ii,3)+1, stda(ii,3)+3]

         osagcs = "Origin shift against global coordinate system"
         CALL write_leaf_to_header(fhh, TRIM(osagcs), ii, stda(ii,:))
         WRITE(free_file_handle, POS=my_pos+8*3) origin_shift

      CASE(6)
         osagcs = "Field content description"
         CALL write_leaf_to_header(fhh, TRIM(osagcs), &
            ii, stda(ii,:)+[LEN_TRIM(field_content_desc), &
            stda(ii,3)+1, stda(ii,3)+LEN_TRIM(field_content_desc)])
         WRITE(free_file_handle, POS=my_pos) TRIM(field_content_desc)
   END SELECT

   CLOSE (free_file_handle)

   CALL meta_delete_empty_file(TRIM(in%p_n_bsnm)//TRIM(suf))
END DO

CLOSE(fhh)

END SUBROUTINE convert_meta_to_puredat


!------------------------------------------------------------------------------
! SUBROUTINE: write_leaf_to_header
!------------------------------------------------------------------------------
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write the leaf data to the PureDat header.
!
!> @param[in] fh File handle
!> @param[in] desc Description of the leaf
!> @param[in] type Type of data o.t.l.
!> @param[in] stda Stream data o.t.l.
!------------------------------------------------------------------------------
SUBROUTINE write_leaf_to_header(fh, desc, type, stda)

CHARACTER(LEN=*), INTENT(IN) :: desc
INTEGER(KIND=meta_ik), INTENT(IN) :: fh, type
INTEGER(KIND=meta_ik), DIMENSION(3), INTENT(IN):: stda ! no_of_data, lb, ub

WRITE(fh, "('<--leaf-->')")
WRITE(fh, "('<description> ', A)") "'"//TRIM(desc)//"'"
WRITE(fh, "('<no_of_data> ', I0)") stda(1)
WRITE(fh, "('<type_of_data> ', I0)") type
WRITE(fh, "('<lower_bound> ', I0)") stda(2)
WRITE(fh, "('<upper_bound> ', I0)") stda(3)

END SUBROUTINE write_leaf_to_header

END MODULE meta_puredat_interface