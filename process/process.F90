program post
use mod_dimensions

real, allocatable :: time(:)
real, allocatable :: ana(:,:)
real, allocatable :: ave(:,:)
real, allocatable :: std(:,:)
integer reclA,i,k,j,iprt,nrt
real, allocatable :: rms(:)

real :: dx=1.0

allocate (time(1))
allocate (ana(nx,1))
allocate (ave(nx,1))
allocate (std(nx,1))
inquire(iolength=reclA)time,ana,ave,std
deallocate (ana,ave,std)

open(10,file='solution.uf',form='unformatted',access='direct',recl=reclA)
   iprt=0
   do
      iprt=iprt+1
      read(10,rec=iprt,err=100)time(1)
      print *,'time',time(1)
   enddo
100 nrt=iprt-2
   deallocate(time)
   
   allocate (time(nrt))
   allocate (ana(nx,nrt))
   allocate (ave(nx,nrt))
   allocate (std(nx,nrt))
   do iprt=1,nrt
      read(10,rec=iprt)time(iprt),ana(:,iprt),ave(:,iprt),std(:,iprt)
   enddo
close(10)

open(10,file='solution.tec')
   write(10,*)'TITLE = "Rossby model"'
   write(10,*)'VARIABLES = "i" "j" "x" "Time" "ANA(1)" "AVE(1)" "STD(1)" "RES(1)"'
   write(10,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',nx,', J=',nrt,', K=1'
   write(10,'(20I5)')((i,i=1,nx),j=1,nrt)
   write(10,'(20I5)')((j,i=1,nx),j=1,nrt)

   write(10,'(10(1x,e12.5))')((dx*float(i-1),i=1,nx),j=1,nrt)
   write(10,'(10(1x,e12.5))')((time(j)      ,i=1,nx),j=1,nrt)

   write(10,'(10(1x,e12.5))')((ana(i,j)            ,i=1,nx),j=1,nrt)
   write(10,'(10(1x,e12.5))')((ave(i,j)            ,i=1,nx),j=1,nrt)
   write(10,'(10(1x,e12.5))')((std(i,j)            ,i=1,nx),j=1,nrt)
   write(10,'(10(1x,e12.5))')((ave(i,j)-ana(i,j)   ,i=1,nx),j=1,nrt)
close(10)

open(10,file='rms.dat')
   allocate (rms(nrt))
   do j=1,nrt
      rms(j)=sqrt(dot_product(ave(:,j)-ana(:,j),ave(:,j)-ana(:,j))/float(nx))
      write(10,'(i5,10g13.5)')j,time(j),rms(j)
   enddo
close(10)

open(10,file='rmstot.dat')
  write(10,'(3g13.5)')sum(rms(:))/float(nrt)
close(10)

end program

