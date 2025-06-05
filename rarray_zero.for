      SUBROUTINE RARRAY_ZERO(n,x)
!***********************************************************************
!RARRAY_ZERO sets the elements of array x to 0.0
!  W.A.Houlberg 12/98
!Input:
!  n-number of elements to be zeroed
!  x-array to be zeroed
!Output:
!  x-zeroed array
!***********************************************************************
      use constant, only : dkind
	IMPLICIT NONE
!Declaration of input variables
      INTEGER        n
      real(kind=dkind)           x(*)
!Declaration of local variables
      INTEGER        i
      DO i=1,n
        x(i)=0.0
      ENDDO   
      RETURN
      END 
