!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
SUBROUTINE getnum(filenm,fileInt,istart,num)
IMPLICIT NONE
INTEGER :: fileInt, num, i, iwhere, istart
CHARACTER(LEN=60) :: filenm
CHARACTER(LEN=60) :: format1
CHARACTER(LEN=20) :: tmpChar1, tmpChar2
   iwhere = 0
   OPEN(UNIT=fileInt,FILE=filenm,FORM='FORMATTED')
   ! search for '!NATOM' in the psf file
   ! the next line(s) SHOULD be the atomic information
   DO i=1,1000
      READ(fileInt,*) tmpChar1, tmpChar2
      IF (TRIM(tmpChar2).EQ.'!NATOM') THEN
         READ(tmpChar1,*) num
         iwhere = i
         EXIT
      END IF
   END DO
   CLOSE(fileInt)
   istart=iwhere
RETURN
END SUBROUTINE getnum
!-----------------------------------------------------------------------------------------
SUBROUTINE psfread(filenm,fileInt,istart,num)
USE coords
USE atoms
IMPLICIT NONE
INTEGER :: i, num, istart, itmp1
INTEGER :: fileInt
CHARACTER(LEN=60) :: filenm
CHARACTER(LEN=60) :: format1
CHARACTER(LEN=20) :: tmpChar1, tmpChar2
CHARACTER(LEN=5)  :: typ1

   format1 = "(I8,1X,2(A5),A4,1X,A4,1X,A5,1X,F9.6,6X,F8.4,6X,I6)"

   OPEN(UNIT=fileInt,FILE=filenm,FORM='FORMATTED')
   DO i=1, istart
      READ(fileInt,*) tmpChar1, tmpChar2
   END DO
   DO i=1, num
      READ(fileInt,*) aindex(i), typ1, resid(i), resname(i), &
                      atomname(i), atomtype(i), charge(i), mass(i) ,itmp1
!      WRITE(*,format1) aindex(i), typ1, resid(i), resname(i), &
!                       atomname(i), atomtype(i), charge(i), mass(i) ,itmp1
   END DO
   CLOSE(fileInt)
RETURN
END SUBROUTINE psfread
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
