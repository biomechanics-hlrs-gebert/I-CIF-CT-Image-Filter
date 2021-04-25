  !---------------------------------------------------------------------------------------------------
  !> \mainpage HLRS - vtk STURCUTRED_POINTS data manipulation tool
  !>
  !> <hr>
  !> \section desc Description
  !>
  !> This program manipulates vtk header datasets to i.e. test other Fortran programs.
  !> All meta data need to be set by hand to prevent bloatware :-)
  !> And to urge you to take a moment of consideration :-)
  !>
  !> <hr>
  !> \section developers Developers
  !>
  !> Johannes Gebert
  !>
  !>  \section modified Last modified:
  !>  by: Johannes Gebert \n
  !>  on: 21.01.2021
  !---------------------------------------------------------------------------------------------------


PROGRAM manipulate_vtk_header

  USE iso_fortran_env

  IMPLICIT NONE

  INTEGER      (KIND=INT64) , PARAMETER                       :: ik=8, rk=8
  INTEGER      (KIND=INT32) , DIMENSION (:,:,:), ALLOCATABLE  :: array_in
  INTEGER      (KIND=INT16) , DIMENSION (:,:,:), ALLOCATABLE  :: array
  INTEGER      (KIND=INT64) , DIMENSION (3)                   :: dims, dims_in
  REAL         (KIND=REAL64), DIMENSION (3)                   :: spcng, spcng_in
  INTEGER      (KIND=INT64) , DIMENSION (3)                   :: flh
  INTEGER      (KIND=INT64)                                   :: xx, yy, zz, hdr_lngth_in, x_lb, x_ub, y_lb, y_ub, z_lb, z_ub
  INTEGER      (KIND=INT64) , PARAMETER                       :: fl_in_un=22, fl_out_un=23
  CHARACTER    (LEN=512)                                      :: fl_in, fl_out

  !-- Manual parametrization -----------------------------------------------------------------------
  fl_in        = "/home/GEB/AAB_programs/cCT_mfCT_Spatial_Registration_Tool/testing/SR_TC1_chequered_1000_2000.vtk"
  fl_out       = "/home/GEB/AAB_programs/cCT_mfCT_Spatial_Registration_Tool/testing/SR_TC1_chequered_1000_2000_low_res.vtk"
  !-- Input
  hdr_lngth_in = 224_ik+10_ik                                  ! 10 bytes newline characters!!!       Obtain with head *.vtk | wc -c 
  dims_in      = (/ 200 , 300 , 400 /)                         ! vtk header
  spcng_in     = (/ 0.0140023, 0.0140023, 0.0140023 /)         ! vtk header
  !-- Output
  spcng        = (/ 0.098, 0.098, 0.098 /)                     ! Dimensions dims are calculated.
  !-- Manual parametrization -----------------------------------------------------------------------


  !-- Programs task
  OPEN(UNIT=fl_in_un, FILE=TRIM(fl_in), CONVERT='little_endian', &
       ACCESS="stream", FORM="unformatted", STATUS="old")
  ALLOCATE(array_in(dims_in(1),dims_in(2),dims_in(3)))
  READ(UNIT=fl_in_un, POS=hdr_lngth_in) array_in(:,:,:)
  CLOSE(fl_in_un)

  !-- Write new vtk File
  OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), ACTION="write", STATUS="new")
  WRITE(fl_out_un,'(A)')            "# vtk DataFile Version 5.1"
  WRITE(fl_out_un,'(A)')            "vtk output"
  WRITE(fl_out_un,'(A)')            "BINARY"
  WRITE(fl_out_un,'(A)')            "DATASET STRUCTURED_POINTS"
  WRITE(fl_out_un,'(A,3(I5,A))')    "DIMENSIONS", dims(1)," ",dims(2)," ",dims(3),""
  WRITE(fl_out_un,'(A,3(F11.6,A))') "SPACING ", spcng(1)," ",spcng(2)," ",spcng(3),""
  WRITE(fl_out_un,'(A)')            "ORIGIN 0 0 0"
  WRITE(fl_out_un,'(A, I11)')       "POINT_DATA", dims(1)*dims(2)*dims(3)
  WRITE(fl_out_un,'(A)')            "SCALARS DICOMImage short"            ! tricky part of this stuff!
  WRITE(fl_out_un,'(A)')            "LOOKUP_TABLE default"
  CLOSE(UNIT=fl_out_un)

  OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), CONVERT='big_endian', &
       ACCESS="stream", FORM="unformatted", STATUS="old", POSITION="append")
  WRITE(UNIT=fl_out_un) array(:,:,:)
  CLOSE(UNIT=fl_out_un)

  OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), ACTION="write", STATUS="old", POSITION="append")
  WRITE(fl_out_un ,'(A)')          "METADATA"
  WRITE(fl_out_un ,'(A)')          "INFORMATION 0"
  WRITE(fl_out_un ,'(A)')
  CLOSE(UNIT=fl_out_un)

END PROGRAM manipulate_vtk_header
