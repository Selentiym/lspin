!program test
!    use constants
!    use matrix
!    use lspinors
!    use laguerre
!    use integrate
!    implicit none
!
!    integer:: N, i,j, degree
!    real(WP), allocatable::x(:)
!    real(WP)::test1
!
!    N = 600
!    degree = 20
!    allocate(x(N))
!    do i=1,N
!        x(i) = i/10.0D+00 + 1.0D+00
!    end do
!    !call setLspinorGlobalParameters(N, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!    call setLspinorGlobalParameters(degree, 1.0D-10, 0.0D+00, N, x)
!    call setLspinorParameters(degree, 'U')
!    open(unit=18, file='test_deriv.dat', status='replace')
!    do i=2,N
!        j=i-1
!        test1 = (laguerrePoly(i, degree) - laguerrePoly(i-1, degree))/(x(i)-x(j))! / laguerrePolyDerivative(i, 3)
!        !if (abs(test1) > 5) then test1 = 0.0_WP; end if
!        write (18,*) x(i), laguerrePoly(i, degree), laguerrePolyDerivative(i, degree), test1
!    end do
!    close(18)
!    deallocate(x)
!end program test
