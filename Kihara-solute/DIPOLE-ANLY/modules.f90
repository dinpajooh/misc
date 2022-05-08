!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
MODULE files
   IMPLICIT NONE
   CHARACTER(LEN=60) :: psffile, pdbfile, dcdfile
   INTEGER, PARAMETER :: psfInt = 100
   INTEGER, PARAMETER :: pdbInt = 101
   INTEGER, PARAMETER :: dcdInt = 102
   INTEGER, PARAMETER :: uelInt = 142
END MODULE files
!------------------------------------------------------------------------------
MODULE coords
   IMPLICIT NONE
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xini, yini, zini
   REAL*4, ALLOCATABLE, DIMENSION(:,:) :: rxf,  ryf,  rzf
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rx,  ry,  rz
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: box
   DOUBLE PRECISION :: boxx, boxy, boxz, boxxSq, boxySq, boxzSq
   DOUBLE PRECISION :: hbx, hby, hbz
   DOUBLE PRECISION :: rcavCut, rcavCutsq, rcavCutb, rcavCutbsq
END MODULE coords
!------------------------------------------------------------------------------
MODULE water
   IMPLICIT NONE
   CHARACTER(LEN=4) :: watertype
   CHARACTER(LEN=3) :: oxname
   CHARACTER(LEN=2) :: hy1name
   CHARACTER(LEN=2) :: hy2name
   CHARACTER(LEN=2) :: omname
END MODULE water
!------------------------------------------------------------------------------
MODULE atoms
   IMPLICIT NONE
   INTEGER :: natom, nframes
   INTEGER, ALLOCATABLE, DIMENSION(:) :: aindex
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: resid, resname
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: atomname, atomtype
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: charge, mass
   INTEGER, ALLOCATABLE, DIMENSION(:) :: pdbaindex, pdbresid
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: pdbresname
   CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: pdbatomname
   REAL, ALLOCATABLE, DIMENSION(:) :: coor, beta
   DOUBLE PRECISION :: qow, qhw, qmw
END MODULE atoms
!------------------------------------------------------------------------------
MODULE ewaldMod
   IMPLICIT NONE
   INTEGER :: maxInd, isqMax
   INTEGER :: llMax, mmMax, nnMax 
   DOUBLE PRECISION :: alpha, alphab, calp
   DOUBLE PRECISION :: rcut, rcutb, rcutsq, rcutbsq
   DOUBLE PRECISION, POINTER, DIMENSION(:) :: kvec
   COMPLEX, POINTER, DIMENSION(:) :: eikx, eiky, eikz
END MODULE ewaldMod
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
MODULE cavResMod
   IMPLICIT NONE
   DOUBLE PRECISION :: resx, resy, resz,                                       &
                       resxSum, resySum, reszSum,                              &
                       resxAvg, resyAvg, reszAvg
   DOUBLE PRECISION :: dExdMxAvg, dEydMyAvg, dEzdMzAvg
   DOUBLE PRECISION :: dEdMAvg
END MODULE cavResMod
!------------------------------------------------------------------------------
MODULE constants
   IMPLICIT NONE
   DOUBLE PRECISION, PARAMETER :: half           = 0.5d0
   DOUBLE PRECISION, PARAMETER :: pi             = 3.141592653589793d0
   DOUBLE PRECISION, PARAMETER :: twopi          = 6.283185307179586d0
   DOUBLE PRECISION, PARAMETER :: sqrtpi         = 1.77245385090552d0
   DOUBLE PRECISION, PARAMETER :: r4pie0         = 138935.4835d0
   DOUBLE PRECISION, PARAMETER :: ke             = 8.9875517873681764d9
   DOUBLE PRECISION, PARAMETER :: NA             = 6.02214179d23
   DOUBLE PRECISION, PARAMETER :: Cpe            = 1.602176487d-19
   DOUBLE PRECISION, PARAMETER :: enUnit         = 1.6605402d-23
   DOUBLE PRECISION, PARAMETER :: kbMol          = 8.31447215d-3
   DOUBLE PRECISION, PARAMETER :: meterpAng      = 1.0d-10
   DOUBLE PRECISION, PARAMETER :: kcalPmol2Joule = 6.947700141d-21
   DOUBLE PRECISION, PARAMETER :: sigTIP3H20     = 3.15d0
   DOUBLE PRECISION, PARAMETER :: sigLJ          = 3.0d0
   DOUBLE PRECISION, PARAMETER :: c              = 2.99792458d8
   DOUBLE PRECISION, PARAMETER :: debye_Cm       = 3.335640952d-30
   DOUBLE PRECISION, PARAMETER :: eVpDebye       = 2.99792458001661d0
END MODULE constants
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
