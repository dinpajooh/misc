!-----------------------------------------------------------------------------------------
PROGRAM Dipanly
USE constants
USE water
USE files
USE coords
USE atoms
USE ewaldMod
USE cavResMod
IMPLICIT NONE
DOUBLE PRECISION :: cavRad, dielect, temp, nnn
DOUBLE PRECISION :: units, kT, rHS, rs, rCA, rCARed
INTEGER :: istart, wattype, nwater, nbeta0, nbeta1, nbeta2, nnwata, rwatCnt
INTEGER :: i, j, inconf, oldnconf
INTEGER :: noInt, nbin, pp
DOUBLE PRECISION :: hbox, rinterval
character(len=20) :: str


   READ(*,'(A)') psffile
   READ(*,'(A)') pdbfile
   READ(*,'(A)') dcdfile
   READ(*,*) wattype
   READ(*,*) rcut
   READ(*,*) isqMax
   READ(*,*) maxInd
   READ(*,*) llMax 
   READ(*,*) mmMax
   READ(*,*) nnMax
   READ(*,*) dielect
   READ(*,*) temp
   READ(*,*) cavRad       ! ASUMPTION : cavRad = rHS + sigLJ
   READ(*,*) noInt


   ALLOCATE(kvec(maxInd))
   kvec(:) = 0.0d0

   units      = r4pie0/100.0d0
   kT         = kbMol*temp
!   eVpDebye   = ke/meterpAng/meterpAng*debye_Cm
   rHS        = cavRad - sigLJ
   rs         = (sigLJ + sigTIP3H20)/2.0d0
   rCA        = rHS + sigLJ
   rCARed     = rCA/sigTIP3H20

   rcutsq     = rcut*rcut
   rcavCut    = rHS + rcut
   rcavCutsq  = rcavCut*rcavCut

   ! initialization of local variables
   oldnconf = 0

   WRITE(*,*) '----------------------------------------------------------------------'
   ! define water type
   IF (wattype.EQ.0) THEN
      watertype = 'TIP3'
      oxname    = 'OH2'
      hy1name   = 'H1'
      hy2name   = 'H2'
      omname    = 'XX'
      qow       = -0.834d0
      qhw       =  0.417d0
      WRITE(*,*) 'USING TIP3P WATER TYPE ...'
   ELSE IF (wattype.EQ.1) THEN
      watertype = 'TIP4'
      oxname    = 'OH2'
      hy1name   = 'H1'
      hy2name   = 'H2'
      omname    = 'OM'
      qow       =  0.0d0
      qhw       =  0.52d0
      qmw       = -1.04d0
      WRITE(*,*) 'USING TIP4P WATER TYPE ...'
   ELSE IF (wattype.EQ.2) THEN
      watertype = 'SPCE'
      oxname    = 'OH2'
      hy1name   = 'H1'
      hy2name   = 'H2'
      omname    = 'XX'
      qow       = -0.8476d0
      qhw       =  0.4238d0
      WRITE(*,*) 'USING SPCE WATER TYPE ...'
   ELSE
      WRITE(*,*) 'UNRECOGNIZED WATER TYPE :'
      WRITE(*,*) 'TIP3P -> wattype = 0 in input file'
      WRITE(*,*) 'TIP4P -> wattype = 1 in input file'
      WRITE(*,*) 'SPCE  -> wattype = 2 in input file'
      STOP
   END IF
   WRITE(*,*) '----------------------------------------------------------------------'

   CALL getnum(psffile,psfInt,istart,natom)
   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'THE NUMBER OF ATOMS REFERENCED IN THE PSF FILE IS ', natom
   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'ALLOCATING APPROPRIATE STRUCTUAL ARRAYS ...'
   CALL alStructArrays
   WRITE(*,*) 'DONE ALLOCATING APPROPRIATE STRUCTUAL ARRAYS ...'
   WRITE(*,*) '----------------------------------------------------------------------'

   CALL getnframes(dcdfile,dcdInt,nframes)
   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'THE NUMBER OF FRAMES REFERENCED IN THE DCD FILE IS ', nframes
   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'ALLOCATING APPROPRIATE DYNAMIC ARRAYS ...'
   CALL alCoordArrays
   WRITE(*,*) 'DONE ALLOCATING APPROPRIATE DYNAMIC ARRAYS ...'
   WRITE(*,*) '----------------------------------------------------------------------'

   ! psfread returns the charge and mass arrays
   CALL psfread(psffile,psfInt,istart,natom)

   ! pdbread returns initial pdb positions and betas
   CALL pdbread(pdbfile,pdbInt,natom)

   ! get number of water atom(s)
   ! get number of beta = 1.0 atom(s) 
   ! and number of beta = 2.0 atom(s)
   CALL configinfo(natom,wattype,nwater,nbeta0,nbeta1,nbeta2,nnwata)
   WRITE(*,*) '     # NATOM  |   WATERS  |  # BETA 0  |  # BETA 1  |  # BETA 2  |  # NON-WATERS '
   WRITE(*,*) '    ------------------------------------------------------------------'
   WRITE(*,*) natom, nwater, nbeta0, nbeta1, nbeta2, nnwata

   WRITE(*,*) '----------------------------------------------------------------------'
   WRITE(*,*) 'READING TRAJECTORY FROM DCD FILE ...'
   CALL dcdread(dcdfile,dcdInt)
   WRITE(*,*) 'DONE READING TRAJECTORY FROM DCD FILE ...'
   WRITE(*,*) '----------------------------------------------------------------------'

   WRITE(*,*) '     STEP  |  FRAME  |    Ecav_X    |    Ecav_Y    |   Ecav_Z   |   Ecav_tot   '
   WRITE(*,*) '-------------------------------------------------------------------------------'

   boxx       = box(1,1); boxy = box(1,2); boxz = box(1,3)
   boxxSq     = boxx*boxx; boxySq = boxy*boxy; boxzSq = boxz*boxz
   alphab     = alpha*boxx
   rcutb      = rcut/boxx
   rcutbsq    = rcutb*rcutb
   rcavCutb   = rcavCut/boxx
   rcavCutbsq = rcavCutb*rcavCutb

 
   hbx = boxx/2.0d0
   hby = boxy/2.0d0
   hbz = boxz/2.0d0

   !CALL ewSetup
 
   rinterval = 0.3
   hbox = 30.0 - 0.1d0   ! Ang
   nbin = int(hbox/rinterval)+1

   do pp = 2, nbin
         open(334+pp, file='m0'//trim(str(pp))//'.dat')
   enddo


   nnn = 0.0d0
   i=0
   !######################
   !## loop over frames ##
   !######################
   DO inconf=1+oldnconf,nframes+oldnconf
      i = i + 1
      nnn = nnn + 1.0d0

      rx(:) = DBLE(rxf(i,:))
      ry(:) = DBLE(ryf(i,:))
      rz(:) = DBLE(rzf(i,:))

      ! get box unit values for Ewald
      CALL Dipole(nnn)   


      IF (MOD(inconf,100).EQ.0) THEN
         WRITE(*,707) inconf, i
      END IF

   END DO
   !##########################
   !## end loop over frames ##
   !##########################


   do pp = 2, nbin
         close(334+pp)
   enddo

   CALL deCoordArrays
   CALL deStructArrays


704 FORMAT(' ELECTROSTATIC FIELD              :', 3(F20.10))
705 FORMAT(' ELECTROSTATIC FIELD (NAMD UNITS) :', 3(F20.10))
706 FORMAT(' DIPOLE MOMENT                    :', 3(F20.10))
707 FORMAT(2(I8),4(F20.10))
708 FORMAT(' CALCULATED CAVITY FIELD          :', 2(I8),4(F20.10))
915   FORMAT(I8,3(F40.20))
916   FORMAT(2(I8))
917   FORMAT(5(F20.10))
918   FORMAT(I8,6(F20.10))
919   FORMAT(I8,8(F20.10))
920   FORMAT(I8,4(F20.10))
921   FORMAT(7(F20.10))
922   FORMAT(3(F20.10))

END PROGRAM Dipanly 
!-----------------------------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   "Convert an integer to string."
      character(len=20) function str(k)
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
      end function str


