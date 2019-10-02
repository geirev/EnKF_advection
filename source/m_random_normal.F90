module m_random_normal
contains

subroutine random_normal(work1,n)
!  Returns a vector of random values N(variance=1,mean=0)
   implicit none
   integer, intent(in) :: n
   real,   intent(out) :: work1(n)
   real, parameter     ::  pi=3.141592653589

   real w1,w2
   integer i

   do i=1,n
      !!do
      call random_number(w1)
      call random_number(w2)
      work1(i)= sqrt(-2.0*log(w1+tiny(1.0)))*cos(2.0*pi*w2)
   enddo

end subroutine random_normal






end module m_random_normal
