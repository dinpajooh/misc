! READS PDBS GENERATED FROM PSFGEN
! THE HEADER DOES NOT CONTAIN BOX INFORMATION
SUBROUTINE pdbread(filenm,fileInt,num)
USE coords
USE atoms
IMPLICIT NONE
DOUBLE PRECISION :: x, y, z, cor, bta
INTEGER :: itmp1, itmp2, itmp3
INTEGER :: i, num, atomind
INTEGER :: fileInt
CHARACTER(LEN=60)  :: filenm
CHARACTER(LEN=6)   :: typ1, typ2, typ3, typ4, typ5, typ6
CHARACTER(LEN=60)  :: format1
CHARACTER(LEN=200) :: dummyheader

   format1 = "(A4,I7,A6,A5,I4,4X,3(F8.3),2(F6.2),3X,A6)"

   OPEN(UNIT=fileInt,FILE=filenm,FORM='FORMATTED')
   !   READ(fileInt,*) dummyheader
   DO i=1, num
      READ(fileInt,format1) typ3, pdbaindex(i), pdbatomname(i), pdbresname(i), pdbresid(i), &
                            xini(i), yini(i), zini(i), coor(i), beta(i), typ4
!      WRITE(*,format1) typ3, pdbaindex(i), pdbatomname(i), pdbresname(i), pdbresid(i), &
!                       xini(i), yini(i), zini(i), coor(i), beta(i), typ4
   END DO
   CLOSE(fileInt)
RETURN
END SUBROUTINE pdbread
