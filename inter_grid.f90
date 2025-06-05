! Course-to-fine prolongation by bilinear interpolation.  The coarse
! grid has nc grid points and the fine grid has nf = 2*nc-1 grid
! points.  The coarse-grid solution is input as uc and the 
! fine-grid solution is returned in uf.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine interp(uc,nc,uf,nf,err)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  implicit none

  integer, intent(in) :: nc,nf
  integer, intent(inout) :: err
  real (kind=dkind), dimension(1:nc,1:nc), intent(in) :: uc ! Coarse Grid
  real (kind=dkind), dimension(1:nf,1:nf), intent(out) :: uf ! Fine Grid
  integer ic,iif,jc,jf

  if ( nf .ne. 2*nc - 1 ) then
     print *, "[interp]: Array Dimensions Incorrect"
     err = 1 ! Error
     return
  end if

  ! Do elements that are copies
  jf = 1
  do jc = 1, nc
     iif = 1
     do ic = 1, nc
        uf(iif,jf) = uc(ic,jc)

        iif = iif + 2
     end do
     jf = jf + 2
  end do

  ! Interpolate horizontally for iif odd and jf even
  do jf = 1, nf, 2
     do iif = 2, nf-1, 2
        uf(iif,jf) = 0.5*(uf(iif-1,jf) + uf(iif+1,jf))
     end do
  end do

  ! Interpolate vertically for all iif and jf odd
  do jf = 2, nf-1, 2
     do iif = 1, nf
        uf(iif,jf) = 0.5*(uf(iif,jf-1) + uf(iif,jf+1))
     end do
  end do

  err = 0 ! Success
end subroutine interp

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine interp_nonuni(uc,nc,uf,nf,err)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  implicit none

  integer, intent(in) :: nc,nf
  integer, intent(inout) :: err
  real (kind=dkind), dimension(1:nc,1:nc), intent(in) :: uc ! Coarse Grid
  real (kind=dkind), dimension(1:nf,1:nf), intent(out) :: uf ! Fine Grid
  integer ic,iif,jc,jf

  if ( nf .ne. 2*nc - 1 ) then
     print *, "[interp]: Array Dimensions Incorrect"
     err = 1 ! Error
     return
  end if

  ! Do elements that are copies
  jf = 1
  do jc = 1, nc
     iif = 1
     do ic = 1, nc
        uf(iif,jf) = uc(ic,jc)

        iif = iif + 2
     end do
     jf = jf + 2
  end do

  ! Interpolate horizontally for iif odd and jf even
  do jf = 1, nf, 2
     do iif = 2, nf-1, 2
!        uf(iif,jf) = 0.5*(uf(iif-1,jf) + uf(iif+1,jf))
		uf(iif,jf) = (dx_a(iif)*uf(iif-1,jf) + dx_a(iif-1)*uf(iif+1,jf))  &
						/(dx_a(iif-1) + dx_a(iif+1))
     end do
  end do

  ! Interpolate vertically for all iif and jf odd
  do jf = 2, nf-1, 2
     do iif = 1, nf
!        uf(iif,jf) = 0.5*(uf(iif,jf-1) + uf(iif,jf+1))
		uf(iif,jf) = (dz_a(jf)*uf(iif,jf-1) + dz_a(jf-1)*uf(iif,jf+1))  &
						/(dz_a(jf-1) + dz_a(jf+1))
     end do
  end do

  err = 0 ! Success
end subroutine interp_nonuni