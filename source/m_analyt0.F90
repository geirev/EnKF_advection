module m_analyt0
contains
subroutine analyt0(sol,nx,u)
   implicit none
   integer, intent(in)  :: nx
   real,    intent(in)  :: u
   real,    intent(out) :: sol(nx)
   
   integer i,ia,ib
   real tmp(nx)

   if (u > 0.0) then
      do i=1,nx
         ia=mod(i-2+nx,nx)+1
         tmp(i)=sol(ia)
      enddo

   else
      do i=1,nx
         ib=mod(i,nx)+1
         tmp(i)=sol(ib)
      enddo
   endif
   sol=tmp


end subroutine analyt0
end module m_analyt0
