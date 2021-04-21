module m_ensrank
contains
subroutine ensrank(A,ave,nx,nrens,it,sig)
   implicit none
   integer, intent(in)    :: nx
   integer, intent(in)    :: nrens
   integer, intent(in)    :: it
   real,    intent(in)    :: A(nx,nrens)
   real,    intent(in)    :: ave(nx)
   real,    intent(inout) :: sig(nrens)
   integer iens,msx,lwork,ierr
   character(len=4) tag4

   real U(1,1),VT(1,1)
   real, allocatable :: AA(:,:)
   real, allocatable :: work(:)

   logical ex

   msx=min(nrens,nx)

   allocate(AA(nx,nrens))
   do iens=1,nrens
      AA(:,iens)=A(:,iens)-ave(:)
   enddo

! Compute SVD of current ensemble
   lwork=20*MIN(nx,nrens);  print *,'lwork=',lwork
   allocate(work(lwork))
   call dgesvd('N', 'O', nx, nrens, AA, nx, sig, U, 1, VT, 1, work, lwork, ierr)
   if (ierr /= 0) print *, 'ierr',ierr

   write(tag4,'(i4.4)')it

   inquire(file='Enssigma',exist=ex)
   if ( .not.ex ) call system('mkdir Enssigma')
   open(10,file='Enssigma/enssigma_'//tag4//'.dat')
      do iens=1,msx
         write(10,'(i4,2g13.5)')iens,sig(iens),sig(iens)/sig(1)
      enddo
   close(10)


   deallocate(AA,work)

end subroutine ensrank
end module m_ensrank
