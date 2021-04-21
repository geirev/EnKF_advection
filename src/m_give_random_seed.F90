module m_give_random_seed
contains
subroutine give_random_seed(val1,val2)
! Sets a random seed based on the system and wall clock time
   implicit none 
   integer, intent(in) :: val1,val2
   integer sze
   logical, save :: first=.true.

   integer, allocatable, dimension(:):: pt

   call RANDOM_SEED(size=sze)
   allocate(pt(sze))


   if (first) then
      call RANDOM_SEED(get=pt)
      print *,'Seed dimension=',sze,pt(:)
      first=.false.
   endif


   pt(1) = val1
   pt(2) = val2
   call RANDOM_SEED(put=pt)
   deallocate(pt)
end subroutine give_random_seed
end module m_give_random_seed
