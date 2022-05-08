!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
SUBROUTINE alStructArrays
USE coords
USE atoms
IMPLICIT NONE

   ALLOCATE(resid(1:natom),resname(1:natom))
   ALLOCATE(aindex(1:natom),atomname(1:natom),atomtype(1:natom))
   ALLOCATE(charge(1:natom),mass(1:natom))

   ALLOCATE(pdbresid(1:natom),pdbresname(1:natom))
   ALLOCATE(pdbaindex(1:natom),pdbatomname(1:natom))
   ALLOCATE(coor(1:natom),beta(1:natom))

   ALLOCATE(xini(1:natom),yini(1:natom),zini(1:natom))

RETURN
END SUBROUTINE alStructArrays
!-------------------------------------------------------------------------------
SUBROUTINE deStructArrays
USE coords
USE atoms
IMPLICIT NONE

   DEALLOCATE(resid,resname)
   DEALLOCATE(aindex,atomname,atomtype)
   DEALLOCATE(charge,mass)

   DEALLOCATE(pdbresid,pdbresname)
   DEALLOCATE(pdbaindex,pdbatomname)
   DEALLOCATE(coor,beta)

   DEALLOCATE(xini,yini,zini)

RETURN
END SUBROUTINE deStructArrays
!-------------------------------------------------------------------------------
SUBROUTINE alCoordArrays
USE coords
USE atoms
IMPLICIT NONE
  ALLOCATE(rxf(1:nframes,1:natom),ryf(1:nframes,1:natom),rzf(1:nframes,1:natom))
  ALLOCATE(box(1:nframes,3))
  ALLOCATE(rx(1:natom),ry(1:natom),rz(1:natom))
RETURN
END SUBROUTINE alCoordArrays
!-------------------------------------------------------------------------------
SUBROUTINE deCoordArrays
USE coords
USE atoms
IMPLICIT NONE
  DEALLOCATE(rxf,ryf,rzf)
  DEALLOCATE(box)
  DEALLOCATE(rx,ry,rz)
RETURN
END SUBROUTINE deCoordArrays
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
