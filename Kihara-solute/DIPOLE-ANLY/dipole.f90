SUBROUTINE Dipole(nnn)   
USE coords   
USE files   
USE constants
USE ewaldMod 
USE atoms 
IMPLICIT NONE

DOUBLE PRECISION :: xij,yij,zij,rijsq,rij
DOUBLE PRECISION :: exij,eyij,ezij
DOUBLE PRECISION :: erfunc
DOUBLE PRECISION :: nnn
INTEGER :: ii, l, m, n
INTEGER :: lmnSq, kCnt
DOUBLE PRECISION :: velectso, sum_velectso, uelectso
DOUBLE PRECISION :: bigvelectso, bigsum_velectso, biguelectso
DOUBLE PRECISION :: velectreal, velectrecip
DOUBLE PRECISION :: right_term, left_term
DOUBLE PRECISION, PARAMETER ::  qqfact = 167072.28


character(len=20) :: str
double precision :: rijdiel, rmindiel, rmaxdiel, rinterval
double precision :: hbox, temp, rcut2
integer :: pp,qq, jj, nbin, ibin, nmolecule
DOUBLE PRECISION, Allocatable, DIMENSION(:) :: mudot
DOUBLE PRECISION, Allocatable, DIMENSION(:) :: mux_bin, muy_bin, muz_bin
DOUBLE PRECISION, Allocatable, DIMENSION(:) :: mux, muy, muz
DOUBLE PRECISION, Allocatable, DIMENSION(:) :: rox, roy, roz


nmolecule = natom/3
rinterval = 0.3 ! Ang
hbox = 30.0-0.1
nbin = int(hbox/rinterval)+1
allocate(mudot(nbin))
allocate(mux_bin(nbin))
allocate(muy_bin(nbin))
allocate(muz_bin(nbin))
allocate(mux(nmolecule))
allocate(muy(nmolecule))
allocate(muz(nmolecule))
allocate(rox(nmolecule))
allocate(roy(nmolecule))
allocate(roz(nmolecule))
mudot(:) = 0.0d0
mux_bin(:) = 0.0d0
muy_bin(:) = 0.0d0
muz_bin(:) = 0.0d0
mux(:) = 0.0d0
muy(:) = 0.0d0
muz(:) = 0.0d0
rox(:) = 0.0d0
roy(:) = 0.0d0
roz(:) = 0.0d0


qq= 0
DO ii = 1,nmolecule
   do jj = 1,3
      qq = qq+1
      mux(ii) = mux(ii) + charge(qq)*rx(qq)
      muy(ii) = muy(ii) + charge(qq)*ry(qq)
      muz(ii) = muz(ii) + charge(qq)*rz(qq)
      if ( jj .eq. 1 ) then
         rox(ii) = rx(qq) 
         roy(ii) = ry(qq) 
         roz(ii) = rz(qq) 
      endif
   enddo
ENDDO

DO ii=1,nmolecule
   xij = rox(ii) 
   yij = roy(ii) 
   zij = roz(ii)
! *** minimum image the pair separations ***
   call mimage(xij,yij,zij,1)
   rijsq = xij*xij + yij*yij + zij*zij
   rij = sqrt(rijsq)
   ibin = int(rij/rinterval)+1


   if(rij .lt. hbox) then
     exij = mux(ii) 
     eyij = muy(ii) 
     ezij = muz(ii)
  ! accumulation of the radial dipole within r-delr and r
     mudot(ibin) = mudot(ibin) +  (exij*xij+eyij*yij+ezij*zij)/rij

     ! \hat{r_j} = xij/rij \hat{x} + yij/rij \hat{y} + zij/rij \hat{z}
     ! \vec{m}_j = exij \hat{x} + eyij \hat{y} + ezij \hat{z}
     right_term = (exij*xij+eyij*yij+ezij*zij)/rij
     ! \hat{r_alpha} = \hat{x} or \hat{y} or \hat{z}
     ! \hat{r_j} = xij/rij \hat{x} + yij/rij \hat{y} + zij/rij \hat{z}

     mux_bin(ibin) = mux_bin(ibin) +  xij/rij*right_term
     muy_bin(ibin) = muy_bin(ibin) +  yij/rij*right_term
     muz_bin(ibin) = muz_bin(ibin) +  zij/rij*right_term
   endif

ENDDO

do ibin = 2, nbin
   write(334+ibin,903) nnn, mux_bin(ibin), muy_bin(ibin), muz_bin(ibin), mudot(ibin)
enddo


903  format(f10.1,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5)


Deallocate(mudot)
Deallocate(mux)
Deallocate(muy)
Deallocate(muz)
Deallocate(mux_bin)
Deallocate(muy_bin)
Deallocate(muz_bin)
Deallocate(rox)
Deallocate(roy)
Deallocate(roz)
!write(uelInt,*) nnn, mudot(ibin) ! K/e

END SUBROUTINE 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EWALD UPDATE SUBROUTINE 
SUBROUTINE ewSetup
USE coords
USE ewaldMod 
USE constants
IMPLICIT NONE 
DOUBLE PRECISION ::  bbb
INTEGER :: l, m, n
INTEGER :: ll,mm,nn, lmnSq, kCnt
DOUBLE PRECISION :: ksq, kxx, kyy, kzz


calp = 6.4d0/boxx
llMax = dint(boxx*calp)+1
mmMax = dint(boxy*calp)+1
nnMax = dint(boxz*calp)+1
alphab = calp*boxx
isqMax =4.0*alphab*alphab
bbb = 1.0d0/4.0d0/alphab/alphab

print*, 'isqMax',  isqMax
print*, 'maxInd',  maxInd
print*, 'llMax',  llMax
print*, 'mmMax',  mmMax
print*, 'nnMax',  nnMax
print*, 'alphab', alphab

!    ** LOOP OVER K-VECTORS. NOTE KX IS NON-NEGATIVE **
kCnt = 0
DO l = -llMax, llMax
   kxx = twopi*DBLE(l)
   DO m = -mmMax, mmMax
      kyy = twopi*DBLE(m)
      DO n = -nnMax, nnMax
         kzz = twopi*DBLE(n)
         lmnSq = l*l + m*m + n*n
         IF ((lmnSq.LT.isqMax) .AND. (lmnSq.NE.0)) THEN
             kCnt = kCnt + 1
             IF (kCnt.GT.maxInd) THEN
                 print*, 'maxInd',  maxInd
                 print*, 'kCnt',  kCnt
                 STOP 'KVEC IS TOO SMALL'
             END IF
             ksq = kxx*kxx + kyy*kyy + kzz*kzz
             kvec(kCnt) = twopi*EXP(-bbb*ksq)/ksq
         ENDIF
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE ewSetup



!=========================================================================================
FUNCTION  erfunc(x) 
IMPLICIT NONE
!   *******************************************************************
!   ** approximation to the complementary error function             **
!   ** reference:                                                    **
!   ** abramowitz and stegun, handbook of mathematical functions,    **
!   **    national bureau of standards, formula 7.1.26               **
!   *******************************************************************
DOUBLE PRECISION :: erfunc, t, x, xsq, tp
DOUBLE PRECISION, PARAMETER :: a1 = 0.254829592d0
DOUBLE PRECISION, PARAMETER :: a2 = -0.284496736d0
DOUBLE PRECISION, PARAMETER :: a3 = 1.421413741d0
DOUBLE PRECISION, PARAMETER :: a4 = -1.453152027d0
DOUBLE PRECISION, PARAMETER :: a5 = 1.061405429d0
DOUBLE PRECISION, PARAMETER :: p  =  0.3275911d0
t  = 1.d0/(1.d0 + p*x)
xsq = x*x
tp = t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))
erfunc = tp*EXP(-xsq)
RETURN
END FUNCTION erfunc 
!========================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mimage (rxuij,ryuij,rzuij,ibox)
USE coords
implicit none
INTEGER ibox
DOUBLE PRECISION rxuij,ryuij,rzuij

if ( rxuij .gt. hbx ) then
     rxuij=rxuij-boxx
else
    if (rxuij.lt.-hbx) rxuij=rxuij+boxx
endif

if ( ryuij .gt. hby ) then
     ryuij=ryuij-boxy
else
     if (ryuij.lt.-hby) ryuij=ryuij+boxy
endif

if (rzuij.gt.hbz) then
    rzuij=rzuij-boxz
else
    if (rzuij.lt.-hbz) rzuij=rzuij+boxz
endif

return
end
