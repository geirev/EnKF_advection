module m_advect
contains 
   subroutine advect(sol,dt,dx,u,mode)
   use mod_dimensions
   implicit none
   integer, intent(in)   :: mode
   real,    intent(in)   :: dx
   real,    intent(in)   :: dt
   real,    intent(in)   :: u
   real,    intent(inout):: sol(nx)

   integer i,ia,ib
   real alpha

   integer, parameter :: lda=4
   real B(nx)

   alpha=u*dt/dx

! Upstream
   if (mode == 1) then
      do i=1,nx
         ia=mod(i-2+nx,nx)+1
         ib=mod(i,nx)+1
         if (u < 0.0) then
            b(i)= sol(i) - alpha*(sol(ib)-sol(i))
         else
            b(i)= sol(i) - alpha*(sol(i)-sol(ia))
         endif
      enddo
      sol=b
   elseif (mode == 2) then
      do i=1,nx
         ia=mod(i-2+nx,nx)+1
         ib=mod(i,nx)+1
         b(i)=sol(i)-0.5*alpha*(sol(ib)-sol(ia))+0.5*alpha**2*(sol(ib)-2.0*sol(i)+sol(ia))
      enddo
      sol=b
   else
      stop 'advect: unknown mode'
   endif

end subroutine advect
end module m_advect
