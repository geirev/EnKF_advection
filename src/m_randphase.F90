module m_randphase
contains
subroutine randphase(C,n1,nrens)
! This routine generates a real orthogonal random matrix.
! The algorithm is the one by
!    Francesco Mezzadri (2007), How to generate random matrices from the classical
!    compact groups, Notices of the AMS, Vol. 54, pp 592-604.
! 1. First a matrix with independent random normal numbers are simulated.
! 2. Then the QR decomposition is computed, and Q will then be a random orthogonal matrix.
! 3. The diagonal elements of R are extracted and we construct the diagonal matrix X(j,j)=R(j,j)/|R(j,j)|
! 4. An updated Q'=Q X is computed, and this is now a random orthogonal matrix with a Haar measure.

   implicit none
   integer, intent(in)  :: nrens
   integer, intent(in)  :: n1
   complex, intent(inout) :: C(n1/2+1,nrens)



   complex, dimension(nrens)  :: diagR
   complex sigma(nrens), work(10*nrens)
   real, parameter :: pi=3.14159253589
   integer i,j,ierr

   real phi(0:n1/2)
   real fampl(0:n1/2,2)

   
   do j=1,nrens
! Calculating uniformly distributed random numbers between 0 and 2 pi
      call random_number(phi)
      phi=2.0*pi*phi

! Generating complex numbers distributed uniformly on the unit circle in the complex plane
      do i=0,n1/2
         C(i+1,j)=cmplx(cos(phi(i)),sin(phi(i)))
      enddo
      fampl(0,2)=0.0
      print '(a,i5,2f12.2)','A1',j,C(1,j)
      C(1,j)=cmplx(abs(C(1,j)),0.0)
      print '(a,i5,2f12.2)','A2',j,C(1,j)
      print '(a,2i5,3f12.2)','B ',j,n1/2,C(10,j),abs(C(10,j))
   enddo


!$OMP CRITICAL
! QR factorization
   call zgeqrf(n1/2+1, nrens, C, n1/2+1, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: zgeqrf ierr=',ierr

   do i=1,nrens
      diagR(i)=C(i,i)/abs(C(i,i))
      print '(a,4f12.4)','R(i,i), diagR(i) ',C(i,i), diagR(i)
   enddo

! Construction of Q
   call zungqr(n1/2+1, nrens, nrens, C, n1/2+1, sigma, work, 10*nrens, ierr )
   if (ierr /= 0) print *, 'randrot: zungqr ierr=',ierr
!$OMP END CRITICAL

   do i=1,nrens
      C(:,i)=C(:,i)*diagR(i)
   enddo


   open(10,file='phase.dat')
      write(10,*)'TITLE = "Orthogonal random phases"'
      write(10,*)'VARIABLES = "i" "real(C)" "Img(C)" "abs(C)"'
      do j=1,100
      write(10,'(a,i3,a,i5,a)')' ZONE T="member(',j,')", F=POINT, I=',n1/2+1,' J=1 K=1'
      do i=1,n1/2+1
         write(10,'(i5,2f10.3,f10.3)')i,C(i,j),abs(C(i,j))
      enddo
      enddo
   close(10)


end subroutine randphase
end module m_randphase
