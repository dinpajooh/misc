!*****************************************************************************************
!*****************************************************************************************
!*****************************************************************************************
SUBROUTINE trapzd(func,a,b,s,n)
IMPLICIT NONE
INTEGER :: n
DOUBLE PRECISION :: a,b,s,func
EXTERNAL :: func
! This routine computes the nth stage of refinement of an extended trapezoidal rule.
! func is input as the name of the function to be integrated between limits a and b,
! also input. When called with n=1, the routine returns as s the crudest estimate of 
! int_a^b f(x) dx. Subsequent a calls with n=2,3,... (in that sequential order) will
! improve the accuracy of s by adding 2n-2 additional interior points. s should not 
! be modified between sequential calls.
INTEGER :: it,j
DOUBLE PRECISION :: del,sumf,tnm,x
   IF (n.EQ.1) THEN
      s = 0.5d0*(b-a)*(func(a) + func(b)) 
   ELSE
      it = 2**(n-2)
      tnm = it 
      del = (b - a)/tnm
      x = a + 0.5d0*del
      sumf = 0.0d0
      DO j=1,it 
         sumf = sumf + func(x)
         x = x + del 
      END DO
      s = 0.5d0*(s + (b - a)*sumf/tnm)
   END IF
RETURN
END SUBROUTINE trapzd
!*****************************************************************************************
!*****************************************************************************************
!*****************************************************************************************
SUBROUTINE polint(xa,ya,n,x,y,dy)
IMPLICIT NONE
INTEGER n
DOUBLE PRECISION :: dy,x,y,xa(n),ya(n)
INTEGER, PARAMETER :: NMAX=10 ! Largest anticipated value of n.
! Given arrays xa and ya, each of length n, and given a value x, this routine returns a
! value y, and an error estimate dy. If P (x) is the polynomial of degree N 1 such that
! P (xai ) = yai ; i = 1; : : : ; n, then the returned value y = P (x).
INTEGER :: i,m,ns
DOUBLE PRECISION :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
   ns = 1
   dif = ABS(x - xa(1))
   ! Here we nd the index ns of the closest table entry,
   DO i = 1,n
      dift = ABS(x - xa(i))
      IF (dift.LT.dif) THEN
         ns=i
         dif = dift
      END IF
      ! and initialize the tableau of c's and d's.
      c(i) = ya(i)
      d(i) = ya(i)
   END DO
   ! This is the initial approximation to y.
   y = ya(ns)
   ns = ns - 1
   DO m=1,n-1 ! For each column of the tableau,
      DO i=1,n-m ! we loop over the current c's and d's and update them.
         ho = xa(i) - x
         hp = xa(i + m)- x 
         w = c(i + 1) - d(i)
         den = ho - hp
         IF (REAL(den).EQ.0.) THEN
            stop 'failure in polint'
         END IF
         ! This error can occur only if two input xa's are (to within roundo ) identical.
         den = w/den
         ! Here the c's and d's are updated.
         d(i) = hp*den
         c(i) = ho*den
      END DO
      ! After each column in the tableau is completed, we decide
      ! which correction, c or d, we want to add to our accu-
      ! mulating value of y, i.e., which path to take through
      ! the tableau - forking up or down. We do this in such a
      ! way as to take the most "straight line" route through the
      ! tableau to its apex, updating ns accordingly to keep track
      ! of where we are. This route keeps the partial approxima-
      ! tions centered (insofar as possible) on the target x. The
      ! last dy added is thus the error indication.
      IF (2*ns.LT.(n - m)) THEN
         dy = c(ns + 1)
      ELSE
         dy = d(ns)
         ns = ns - 1
      END IF
      y = y + dy
   END DO
RETURN
END SUBROUTINE polint
!*****************************************************************************************
!*****************************************************************************************
!*****************************************************************************************
SUBROUTINE qromb(func,a,b,ss)
IMPLICIT NONE
DOUBLE PRECISION :: a,b,func,ss
EXTERNAL :: func
DOUBLE PRECISION, PARAMETER :: EPS=1.0d-6
INTEGER, PARAMETER :: JMAX=20
INTEGER, PARAMETER :: JMAXP=JMAX+1
INTEGER, PARAMETER :: K=5
INTEGER, PARAMETER :: KM=K-1
! USESpolint,trapzd
! Returns as ss the integral of the function func from a to b. Integration is
! performed by Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule.
! Parameters: EPS is the fractional accuracy desired, as determined by the
! extrapolation error estimate; JMAX limits the total number of steps;
! K is the number of points used in the extrapolation.
INTEGER :: j 
DOUBLE PRECISION :: dss,h(JMAXP),s(JMAXP) 
   h(1)=1. 
   DO j=1,JMAX
      CALL trapzd(func,a,b,s(j),j) 
      IF (j.GE.K) THEN
         CALL polint(h(j - KM),s(j - KM),K,0.0d0,ss,dss)
         IF (ABS(dss).LE.EPS*ABS(ss)) THEN
            RETURN
         END IF
      END IF
      s(j+1) = s(j)
      h(j+1) = 0.25d0*h(j)
   END DO
   stop 'too many steps in qromb'
END SUBROUTINE qromb
!*****************************************************************************************
!*****************************************************************************************
!*****************************************************************************************

