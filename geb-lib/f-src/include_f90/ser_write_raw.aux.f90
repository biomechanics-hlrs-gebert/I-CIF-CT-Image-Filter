INTEGER(ik), INTENT(IN) :: fh
CHARACTER(*), INTENT(IN) :: filename
CHARACTER(*), INTENT(IN), OPTIONAL :: representation

CHARACTER(scl) :: cnvrt

IF(PRESENT(representation)) THEN
   cnvrt = TRIM(representation)
ELSE
   cnvrt = 'LITTLE_ENDIAN'
END IF

OPEN (UNIT=fh, FILE=TRIM(filename), ACCESS="STREAM", FORM="UNFORMATTED", &
   CONVERT=TRIM(cnvrt), STATUS="UNKNOWN", POSITION="APPEND")                                       
WRITE(UNIT=fh) array
CLOSE(UNIT=fh)