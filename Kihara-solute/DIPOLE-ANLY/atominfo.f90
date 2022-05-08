!-----------------------------------------------------------------------------------------
SUBROUTINE configinfo(num,wattyp,nwater,nbeta0,nbeta1,nbeta2,nnwata)
USE coords
USE atoms
IMPLICIT NONE
INTEGER :: wattyp
INTEGER :: i, num, nwata, nnwata, nwater, nbeta0, nbeta1, nbeta2
DOUBLE PRECISION, DIMENSION(1:3) :: rbeta1

   nwata  = 0
   nnwata = 0
   nwater = 0
   nbeta0 = 0
   nbeta1 = 0
   nbeta2 = 0

   IF (wattyp.EQ.0) THEN

      DO i=1,num
         ! total waters with beta equal to 2
         !IF (TRIM(resname(i)).EQ.'TIP3' .AND. INT(beta(i)).EQ.2) THEN
         IF (TRIM(resname(i)).EQ.'TIP') THEN
            nwata = nwata + 1
         END IF
         IF (TRIM(resname(i)).NE.'TIP') THEN
            nnwata = nnwata + 1
         END IF
         IF (INT(beta(i)).EQ.0) THEN
            nbeta0 = nbeta0 + 1
         END IF
         IF (INT(beta(i)).EQ.1) THEN
            nbeta1 = nbeta1 + 1
         END IF
         IF (INT(beta(i)).EQ.2) THEN
            nbeta2 = nbeta2 + 1
         END IF
      END DO

      IF (MOD(nwata,3).NE.0) THEN
         WRITE(*,*) 'TIP3P : COUNTED ', nwata,' NUMBER OF WATER ATOMS WHICH'
         WRITE(*,*) 'IS NOT DIVISIBLE BY 3 ...'
         STOP
      END IF

      nwater = nwata/3

   ELSE IF (wattyp.EQ.1) THEN

      DO i=1,num
         ! total waters with beta equal to 2
         !IF (TRIM(resname(i)).EQ.'TIP4' .AND. INT(beta(i)).EQ.2) THEN
         IF (TRIM(resname(i)).EQ.'TIP4') THEN
            nwata = nwata + 1
         END IF
         IF (TRIM(resname(i)).NE.'TIP4') THEN
            nnwata = nnwata + 1
         END IF
         IF (INT(beta(i)).EQ.0) THEN
            nbeta0 = nbeta0 + 1
         END IF
         IF (INT(beta(i)).EQ.1) THEN
            nbeta1 = nbeta1 + 1
         END IF
         IF (INT(beta(i)).EQ.2) THEN
            nbeta2 = nbeta2 + 1
         END IF
      END DO

      IF (MOD(nwata,4).NE.0) THEN
         WRITE(*,*) 'TIP4P : COUNTED ', nwata,' NUMBER OF WATER ATOMS WHICH '
         WRITE(*,*) 'IS NOT DIVISIBLE BY 4 ...'
         STOP
      END IF

      nwater = nwata/4

   ELSE IF (wattyp.EQ.2) THEN

      DO i=1,num
         ! total waters with beta equal to 2
         !IF (TRIM(resname(i)).EQ.'SPCE' .AND. INT(beta(i)).EQ.2) THEN
         IF (TRIM(resname(i)).EQ.'SPCE') THEN
            nwata = nwata + 1
         END IF
         IF (TRIM(resname(i)).NE.'SPCE') THEN
            nnwata = nnwata + 1
         END IF
         IF (INT(beta(i)).EQ.0) THEN
            nbeta0 = nbeta0 + 1
         END IF
         IF (INT(beta(i)).EQ.1) THEN
            nbeta1 = nbeta1 + 1
         END IF
         IF (INT(beta(i)).EQ.2) THEN
            nbeta2 = nbeta2 + 1
         END IF
      END DO

      IF (MOD(nwata,3).NE.0) THEN
         WRITE(*,*) 'SPCE : COUNTED ', nwata,' NUMBER OF WATER ATOMS WHICH'
         WRITE(*,*) 'IS NOT DIVISIBLE BY 3 ...'
         STOP
      END IF

      nwater = nwata/3

   ELSE

      WRITE(*,*) 'UNRECOGNIZED WATER TYPE ...'

   END IF

RETURN
END SUBROUTINE configinfo
!-----------------------------------------------------------------------------------------
