module m_comp_residual
contains
real function comp_residual(ave,ana,nx)
   implicit none
   integer, intent(in) :: nx
   real, intent(in)    :: ave(nx)
   real, intent(in)    :: ana(nx)
   integer i

   comp_residual=0.0
   do i=1,nx
      comp_residual=comp_residual+(ave(i)-ana(i))**2
   enddo

end function
end module
