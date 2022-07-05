CHARACTER(*), INTENT(IN) :: filename, type
INTEGER(ik), INTENT(IN) :: dims(3)
REAL(rk), INTENT(IN) :: origin(3), spcng(3)

INTEGER(ik) :: fh

fh = give_new_unit()

CALL write_vtk_struct_points_header(fh, TRIM(filename), TRIM(type), spcng, origin, dims)

CALL ser_write_raw(fh, TRIM(filename), array, 'BIG_ENDIAN')

CALL write_vtk_struct_points_footer(fh, TRIM(filename))