module m_ensrank
contains
subroutine ensrank(A,ave,nx,nrens,it,sig)
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: nrens
   integer, intent(in) :: it
   real, intent(in)    :: A(nx,nrens)
   real, intent(in)    :: ave(nx)
   real, intent(inout) :: sig(nrens)
   integer iens,msx,lwork,ierr
   character(len=4) tag4

   real, allocatable :: AA(:,:)
   real, allocatable :: U(:,:)
   real, allocatable :: VT(:,:)
   real, allocatable :: work(:)

   logical ex

   msx=min(nrens,nx)

   allocate(AA(nx,nrens))
   do iens=1,nrens
      AA(:,iens)=A(:,iens)-ave(:)
   enddo

   print *,'ensrank:  A'
! Compute SVD of current ensemble
   lwork=2*max(3*nrens+max(nx,nrens),5*nrens)
   lwork=5*MIN(nx,nrens)  
   lwork=2*lwork
   lwork=61500
   print *,'lwork=',lwork
   print *,'ensrank:  A'
   allocate( U(nx,msx), VT(msx,msx), work(lwork))
   print *,'ensrank:  dgesvd'
   call dgesvd('N', 'N', nx, nrens, AA, nx, sig, U, nx, VT, nrens, work, lwork, ierr)
   print *,'ensrank:  dgesvd done',work(1)
   if (ierr /= 0) print *, 'ierr',ierr

   print *,'ensrank:  A'
   write(tag4,'(i4.4)')it

   print *,'ensrank:  A'
   inquire(file='Enssigma',exist=ex)
   if ( .not.ex ) call system('mkdir Enssigma')
   open(10,file='Enssigma/enssigma_'//tag4//'.dat')
      do iens=1,msx
         write(10,'(i4,2g13.5)')iens,sig(iens),sig(iens)/sig(1)
      enddo
   close(10)

   print *,'ensrank:  A'

   deallocate(AA,U,VT,work)

end subroutine ensrank
end module m_ensrank
