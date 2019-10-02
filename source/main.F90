program main
   use mod_dimensions
   use m_multa
   use m_random
   use m_advect
   use m_analyt0
   use m_dumpsol
   use m_sample1D
   use m_ensmean
   use m_ensvar
   use m_comp_residual
   use m_set_random_seed2
   use m_give_random_seed
   use m_enkf
   use m_ensrank
   use m_measurements
   use m_random_normal
   use mod_shapiro
   implicit none

! Variables read from infile
   integer nrt                           ! Number of timesteps
   real u                                ! Advection spead
   integer iprtint                       ! Print time interval
   real rh                               ! Horizontal correlation
   real const                            ! mean of analytical solution
   integer startsol                      ! 0 zero, 1 analytical, 2 random
   real dx                               ! horizontal grid spacing
   real dt                               ! time step
   integer nrens                         ! ensemble size
   integer nre_ini,nre_sys,nre_obs       ! Number of ensembles used in sampling
   integer nrr_ini,nrr_sys,nrr_obs       ! Number of singular vectors used in improved sampling
   integer mode_analysis                 ! 1 standard, 2 fixed R
   logical samp_fix
   real inivar 
   real obsvar
   real sysvar
   integer iadv                          ! advection scheme
   integer nrobs                         ! Number of measurement per assimilation time
   logical mkobs                         ! Create or read measurements
   logical Rexact                        ! Use exact(true) or lowrank(false) R matrix
   real obsdt                            ! time between assimilation times
   logical :: lrandrot=.true.            ! random rotation in SQRT schemes
   logical :: lsymsqrt=.true.            ! Always use the symmetrical square root rather than one-sided

! inflation
   integer inflate                       ! 0--no inflation, 1--constant inflation, 2--adaptive inflation
   real infmult                          ! constant inflation or adjustment of adaptive inflation

! local analysis
   integer local
   real obs_radius,obs_truncation

   real truncation
   character(len=8) covmodel
   real rd                               ! Horizontal correlation of observation errors in Gaussian case
   character(len=9) mesopt


! other variables
   logical ex
   integer iprt                          ! print record counter
   integer iens,i,j,m,k,nn,iobs          ! Counters etc.


   real time                 
   real, allocatable :: mem(:,:)
   real, allocatable :: sample(:,:)

   real ana(nx)  ! analytical solution
   real ave(nx)  ! ensemble average
   real var(nx)  ! ensemble variance

   real,    allocatable :: enssing(:,:)
   real,    allocatable :: residual(:)
   real,    allocatable :: stddev(:)
   integer, allocatable :: obspos(:)     ! Position of data point
   real,    allocatable :: obs(:)        ! Observations
   integer nrobst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer :: iana=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shapiro filter variables (used with adaptive localization)
   logical lshapiro
   integer, parameter :: shdim=8
   real sh(shdim)
   integer :: ish=2
   real, allocatable :: x(:),y(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reading input data
   open(10,file='infile.in')
      read(10,*)nrt               ; print *,'nrt=         ',nrt
      read(10,*)u                 ; print *,'u=           ',u
      read(10,*)iprtint           ; print *,'iprtint=     ',iprtint
      read(10,*)rh                ; print *,'rh=          ',rh
      read(10,*)const             ; print *,'const=       ',const
      read(10,*)startsol          ; print *,'startsol=    ',startsol


      read(10,*)dx                ; print *,'dx=          ',dx
      read(10,*)dt                ; print *,'dt=          ',dt

      read(10,*)nrens             ; print *,'nrens=       ',nrens

      read(10,*)nre_sys,nrr_sys   ; print *,'nre_sys=     ',nre_sys,nrr_sys
      if (nre_sys > 1 .and. nrr_sys > nre_sys) stop 'nrr_sys > nre_sys)'

      read(10,*)nre_ini,nrr_ini   ; print *,'nre_ini=     ',nre_ini,nrr_ini
      if (nre_ini > 1 .and. nrr_ini > nre_ini) stop 'nrr_ini > nre_ini)'

      read(10,*)nre_obs,nrr_obs   ; print *,'nre_obs=     ',nre_obs,nrr_obs
      if (nre_obs > 1 .and. nrr_obs > nre_obs) stop 'nrr_obs > nre_obs)'

      read(10,'(1x,l1)')samp_fix  ; print *,'samp_fix=    ',samp_fix

      read(10,*)inivar            ; print *,'inivar=      ',inivar
      read(10,*)obsvar            ; print *,'obsvar=      ',obsvar
      read(10,*)sysvar            ; print *,'sysvar=      ',sysvar

      read(10,*)iadv              ; print *,'iadv=        ',iadv
      read(10,*)mode_analysis     ; print *,'mode_ana=    ',mode_analysis

      read(10,*)nrobs             ; print *,'nrobs=       ',nrobs
      read(10,*)obsdt             ; print *,'obsdt=       ',obsdt
      read(10,'(1x,l1)')mkobs     ; print *,'mkobs=       ',mkobs

      read(10,*)truncation        ; print *,'truncation=  ',truncation
      read(10,'(1x,a)')covmodel   ; print *,'covmodel=    ',trim(covmodel)
      read(10,*)rd                ; print *,'rd      =    ',rd
      read(10,'(1x,l1)')Rexact    ; print *,'Rexact=      ',Rexact
      read(10,'(1x,a)')mesopt     ; print *,'mesopt=      ',mesopt
      read(10,'(1x,l1)')lshapiro  ; print *,'lshapiro=    ',lshapiro
      read(10,'(1x,l1)')lrandrot  ; print *,'lrandrot=    ',lrandrot
      read(10,*)inflate,infmult   ; print *,'inflation=   ',inflate,infmult
      read(10,*)local,obs_radius,obs_truncation; print *,'localization=',local,obs_radius,obs_truncation
   close(10)

   call set_random_seed2

   call system('rm -f eigenvalues.dat')

   allocate (mem(nx,nrens))
   allocate (sample(nx,nrens))
   allocate(residual(0:nrt))
   allocate(stddev(0:nrt))

   iprt=0
   iobs=0
   time=0.0

   nrobst=nrt/obsdt
   print *,'nrobst= ',nrobst
   allocate(enssing(nrens,0:nrobst))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate and define measurement setup
   allocate(obs(nrobs))
   allocate(obspos(nrobs))

   select case(mesopt)
   case('uniformly')
!     Uniformly distributed measurements
      nn=nint(float(nx)/float(nrobs))
      obspos(1)=nint(real(nn)/2.0)
      do m=2,nrobs
         obspos(m)=min(obspos(m-1)+nn,nx)
      enddo
   case('clustered')
!     Measurement clustered on neighbouring grid points in second half of interval
      obspos(1)=nint(real(nx)/2.0)
      do m=2,nrobs
         obspos(m)=min(obspos(m-1)+2,nx)
      enddo
   case default
      stop 'invalid mesopt'
   end select

   open(10,file='obspos.dat')
      do m=1,nrobs
         write(10,'(i4,i4)')m,obspos(m)
      enddo
   close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The true solution is a perturbation around the value "const" where the 
! perturbation is a smooth pseudo random field drawn from  N(0,1,rh).
   call sample1D(ana,nx,1,1,1,dx,rh,.false.,.true.)
   ana=ana+const
   print *,'main: ana ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First guess solution
   call sample1D(ave,nx,1,1,1,dx,rh,.false.,.true.)
!   ave=ave+const  ! First guess is a random perturbation from N(0,1,rh) added to "const" (Gives residual=2)
!   ave=const      ! First guess is just equal to "const"
!   ave=ave+ana    ! First guess is a random perturbation from N(0,1,rh) added to the truth
   ave=(ave + ana-const)/sqrt(2.0) +  const   ! Modified 2019 to make consistent initial error
   print *,'main: fg ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of ensemble
   call sample1D(sample,nx,nrens,nre_ini,nrr_ini,dx,rh,samp_fix,.true.)
   do i=1,nrens
      mem(:,i)=ave(:) + sqrt(inivar)*sample(:,i)
   enddo

   print *,'main: ensemble ok'

   call ensmean(mem,ave,nx,nrens)
   call ensvar(mem,ave,var,nx,nrens)
   residual(0)=comp_residual(ave,ana,nx)
   stddev(0)=sqrt(sum(var(1:nx))/float(nx))
   obs=0.0
   call dumpsol(time,ana,ave,var,nx,iprt,dx,obs,obsvar,obspos,nrobs,mem,nrens)
   !call ensrank(mem,ave,nx,nrens,nint(time),enssing(:,0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time stepping
   print *,'main: start time stepping'
   do k=1,nrt   
      time=time+dt
      print '(a,i6,f10.2)','timestep ',k,time

! Advection
      if (iadv == 3) then
         do iens=1,nrens
            call analyt0(mem(:,iens),nx,u)
         enddo
      else
         do iens=1,nrens
            call advect(mem(:,iens),0.5*dt,dx,u,iadv)
            call advect(mem(:,iens),0.5*dt,dx,u,iadv)
         enddo
      endif
      call analyt0(ana,nx,u)

! System noise
      if (sysvar > 0.0) then
         call sample1D(sample,nx,nrens,nre_sys,nrr_sys,dx,rh,samp_fix,.true.)
         do i=1,nrens
            mem(:,i)=mem(:,i)+sqrt(2.0*sysvar*dt)*sample(:,i)
         enddo
      endif


! Assimilation step
      if (mod(time,obsdt).lt.0.05*dt) then
         iobs=iobs+1

         call ensmean(mem,ave,nx,nrens)
         call ensvar(mem,ave,var,nx,nrens)
         call measurements(ana,nx,obs,obspos,nrobs,obsvar,iobs,mkobs,time)
         
         call dumpsol(time,ana,ave,var,nx,iprt,dx,obs,obsvar,obspos,nrobs,mem,nrens)

         call enkf(mem,nx,nrens,obs,obsvar,obspos,nrobs,1,1,mode_analysis,&
                  &truncation,covmodel,dx,rh,Rexact,rd,lrandrot,lsymsqrt,&
                  &inflate,infmult,local,obs_radius,obs_truncation)

         if (lshapiro) then
            call shfact(ish,sh)
            allocate(x(nx),y(nx))
            do j=1,nrens
               x(:)=mem(:,j)
               call shfilt(ish,sh,nx,x,1,y,1,shdim)
               mem(:,j)=y(:)
            enddo
            deallocate(x,y)
         endif

         iana=iana+1

         call ensmean(mem,ave,nx,nrens)
         call ensvar(mem,ave,var,nx,nrens)
         call dumpsol(time,ana,ave,var,nx,iprt,dx,obs,obsvar,obspos,nrobs,mem,nrens)


         !call ensrank(mem,ave,nx,nrens,nint(time),enssing(:,iobs))
      endif

      call ensmean(mem,ave,nx,nrens)
      call ensvar(mem,ave,var,nx,nrens)
      residual(k)=comp_residual(ave,ana,nx)
      stddev(k)=sqrt(sum(var(1:nx))/float(nx))

   enddo

! Print residuals
   open(11,file='residual.dat')
      time=0.0
      do i=0,nrt
         write(11,'(i4,3f10.4)')i,time,sqrt( residual(i)/float(nx) ),stddev(i)
         time=time+dt
      enddo
   close(11)

   inquire(file='residual_tot.dat',exist=ex)
   if (ex) then
      open(11,file='residual_tot.dat',position='append')
   else 
      open(11,file='residual_tot.dat')
   endif
!   write(11,*)sqrt( sum(residual(0:nrt))/float(nx*(nrt+1)) )
   write(11,*)sqrt( sum(residual(1:nrt))/float(nx*(nrt)) )
!   write(11,*)sqrt( residual(nrt)/float(nx))
   print '(a,f12.5)','Total residual is',sqrt( sum(residual(0:nrt))/float(nx*(nrt+1)) )
   print '(a,f12.5)','Final residual is',sqrt( sum(residual(nrt-50:nrt))/float(nx*(50+1)) )
   close(11)

   open(10,file='enssing.tec')
   write(10,*)'TITLE = "Enssing"'
   write(10,*)'VARIABLES = "iens" "iobs" "enssing"'
   write(10,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',nrens,', J=',nrobst,', K=1'
   write(10,'(20I5)')((i                      ,i=1,nrens),j=1,nrobst)
   write(10,'(20I5)')((j                      ,i=1,nrens),j=1,nrobst)
   write(10,'(10(1x,e12.5))')((enssing(i,j)   ,i=1,nrens),j=1,nrobst)
   close(10)


end program main

