!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
SUBROUTINE getnframes(filenm,fileInt,nframesdcd)
USE atoms
IMPLICIT NONE
INTEGER :: fileInt
CHARACTER(LEN=60) :: filenm
INTEGER :: i, natomdcd, nframesdcd

REAL             :: dummyr
INTEGER          :: dummyi
CHARACTER(LEN=4) :: dummyc

   OPEN(UNIT=fileInt,FILE=filenm,STATUS='OLD',FORM='UNFORMATTED')
   READ(fileInt) dummyc, nframesdcd, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
   READ(fileInt) dummyi, dummyr
   READ(fileInt) natomdcd
   CLOSE(fileInt)

   IF (natom.NE.natomdcd) THEN
      WRITE(*,*) 'SOMETHING IS WRONG ...'
      WRITE(*,*) 'THE NUMBER OF ATOMS READ FROM THE'
      WRITE(*,*) 'PSF FILE DOES NOT MATCH THE DCD FILE.'
      STOP
   END IF

RETURN
END SUBROUTINE getnframes
!-----------------------------------------------------------------------------------------
SUBROUTINE dcdread(filenm,fileInt)
USE atoms
USE coords
IMPLICIT NONE
INTEGER :: fileInt
CHARACTER(LEN=60) :: filenm

REAL             :: dummyr
INTEGER          :: dummyi
CHARACTER(LEN=4) :: dummyc
DOUBLE PRECISION, DIMENSION(1:6) :: d

INTEGER :: i, j, natomdcd, nframesdcd
LOGICAL, PARAMETER :: DCDUnitCell=.TRUE.

   OPEN(UNIT=fileInt,FILE=filenm,STATUS='OLD',FORM='UNFORMATTED')
   READ(fileInt) dummyc, nframesdcd, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
   READ(fileInt) dummyi, dummyr
   READ(fileInt) natomdcd

   DO i = 1, nframes
      IF (DCDUnitCell) THEN
         READ(fileInt) (d(j), j=1, 6)
         box(i,1) = d(1); box(i,2) = d(3); box(i,3) = d(6)
      ENDIF
      READ(fileInt) (rxf(i,j),j=1,natom)
      READ(fileInt) (ryf(i,j),j=1,natom)
      READ(fileInt) (rzf(i,j),j=1,natom)
   END DO

   CLOSE(fileInt)
 
END SUBROUTINE dcdread
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
