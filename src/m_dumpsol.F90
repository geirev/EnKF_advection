module m_dumpsol
contains
subroutine dumpsol(time,ana,ave,var,nx,iprt,dx,obs,obsvar,obspos,nrobs,mem,nrsamp,xx)

   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: nrsamp
   integer, intent(in) :: nrobs
   integer, intent(in) :: obspos(nrobs)
   integer, intent(inout) :: iprt
   real, intent(in) :: time
   real, intent(in) :: dx
   real, intent(in) :: ana(nx)
   real, intent(in) :: mem(nx,nrsamp)
   real, intent(in) :: ave(nx)
   real, intent(in) :: var(nx)
   real, intent(in) :: obs(nrobs)
   real, intent(in) :: obsvar
   integer reclA,i,m
   character(len=1) xx
   character(len=4) tag
   character(len=6) outtag
   logical ex

   iprt=iprt+1

   !print *,'dumpsol A'
   inquire(iolength=reclA)time,ana,ave,var
   open(10,file='solution.uf',form='unformatted',access='direct',recl=reclA)
      write(10,rec=iprt)time,ana,ave,sqrt(var+0.000001)
   close(10)


   inquire(file='Solution',exist=ex)
   if ( .not.ex ) call system('mkdir Solution')

   inquire(file='Members',exist=ex)
   if ( .not.ex ) call system('mkdir Members')

   !print *,'dumpsol B'
   outtag='      '
   write(outtag,'(f0.2)')time
   print *,'time tag',time,'--',outtag,'--'

   write(tag,'(i4.4)')iprt
   open(10,file='Solution/sol'//trim(outtag)//xx//'.dat')
      do i=1,nx
         write(10,'(5f10.4)')time,float(i-1)*dx,ana(i),ave(i),sqrt(var(i)+0.00001)
      enddo
   close(10)

   !print *,'dumpsol C'
   open(10,file='Solution/obs'//trim(outtag)//xx//'.dat')
      do m=1,nrobs
         write(10,'(5f10.4)')time,(obspos(m)-1)*dx,obs(m),-sqrt(obsvar),sqrt(obsvar)
      enddo
   close(10)

   !print *,'dumpsol D'
   open(10,file='Members/mem_'//tag//'.dat')
      do i=1,nx
         write(10,'(12f10.4)')time,float(i-1)*dx,mem(i,1:10)
      enddo
   close(10)
   !print *,'dumpsol E'

end subroutine dumpsol
end module m_dumpsol


