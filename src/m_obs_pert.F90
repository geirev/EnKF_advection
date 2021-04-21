module m_obs_pert
contains
subroutine obs_pert(E,nrens,nrobs,fixsamp,dx,rh,covmodel,obspos)
   use mod_dimensions
   use m_random
   use m_pseudo1D
   use m_fixsample1D
   implicit none
   integer, intent(in)    :: nrobs
   integer, intent(in)    :: nrens
   real,    intent(out)   :: E(nrobs,nrens)
   logical, intent(in)    :: fixsamp
   real,    intent(in)    :: dx
   real,    intent(in)    :: rh
   character(len=8),  intent(in)    :: covmodel
   integer,    intent(in)    :: obspos(nrobs)

   real, allocatable :: EEfield(:,:)
   integer iens,m

   if (trim(covmodel) == 'diagonal') then
      call random(E,nrobs*nrens)   

   elseif(trim(covmodel) == 'gaussian') then
      allocate(EEfield(nx,nrens))
      call pseudo1D(EEfield,nx,nrens,rh,dx,nx)

      do iens=1,nrens
      do m=1,nrobs
         E(m,iens)=EEfield(obspos(m),iens)
      enddo
      enddo
      deallocate(EEfield)
   else
      print *,'Problem with covmodel:+++',trim(covmodel),'+++'
      stop
   endif

   if (fixsamp) call fixsample1D(E,nrobs,nrens)

end subroutine obs_pert
end module m_obs_pert
