module m_analyt
contains
subroutine analyt(sol,nx,time,dx,u,const)
   implicit none
   integer, intent(in)  :: nx
   real,    intent(in)  :: time
   real,    intent(in)  :: dx
   real,    intent(in)  :: u
   real,    intent(in)  :: const
   real,    intent(out) :: sol(nx)
   
   real, parameter :: pi=3.1415927
   real x,length
   integer i,k,nrmodes

   length=(nx-1)*dx
   nrmodes=16

   do i=1,nx
      sol(i)=const
      x=float(i-1)*dx
      do k=1,nrmodes
         sol(i)=sol(i)+sin(float(k)*2.0*pi*(x-u*time-float(k**3))/length)
      enddo 
   enddo

end subroutine analyt
end module m_analyt
