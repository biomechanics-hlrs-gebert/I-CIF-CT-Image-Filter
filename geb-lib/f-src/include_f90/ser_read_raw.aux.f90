INTEGER(ik), INTENT(IN) :: fh
CHARACTER(*), INTENT(IN) :: filename
CHARACTER(*), INTENT(IN), OPTIONAL :: representation

CHARACTER(scl) :: cnvrt

cnvrt = 'NATIVE'
IF(PRESENT(representation)) cnvrt = TRIM(representation)
   
OPEN (UNIT=fh, FILE=TRIM(filename), ACCESS="STREAM", FORM="UNFORMATTED", &
   CONVERT=TRIM(cnvrt), STATUS="UNKNOWN")                                       
READ(UNIT=fh) array
CLOSE(UNIT=fh)