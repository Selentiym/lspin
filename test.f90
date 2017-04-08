!program test
!    use constants
!    use matrix
!    use lspinors
!    use laguerre
!    use integrate
!    implicit none
!
!    integer:: N, i,j, degree
!    character:: letter
!    real(WP), allocatable::x(:), lspin1(:),lspin2(:)
!    real(WP)::test1
!
!    N = 100
!    degree = 30
!    allocate(x(N))
!    allocate(lspin1(N))
!    allocate(lspin2(N))
!    do i=1,N
!        x(i) = i/10.0D+00 + 1.0D+00
!    end do
!    !call setLspinorGlobalParameters(N, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!    call setLspinorGlobalParameters(degree, 1.0D+00, 1D+00, N, x)
!    call setLspinorParameters(degree, 'U')
!    open(unit=18, file='test_deriv.dat', status='replace')
!    open(unit=19, file='test_lspinor_deriv.dat', status='replace')
!    !Проверяем производную полиномов Лагерра
!    do i=2,N
!        j=i-1
!        test1 = (laguerrePoly(i, degree) - laguerrePoly(i-1, degree))/(x(i)-x(j))! / laguerrePolyDerivative(i, 3)
!        !if (abs(test1) > 5) then test1 = 0.0_WP; end if
!        write (18,*) x(i), laguerrePoly(i, degree), laguerrePolyDerivative(i, degree), test1
!    end do
!    close(18)
!
!    degree = 10
!    letter = 'U'
!    !lspin1 = lspinorVector(degree, letter)
!
!    lspin1 = lspinorVectorWithFunc(ampl, degree, letter)
!
!    do i=2,N
!       j = i - 1
!       test1 = (lspin1(i) - lspin1(j))/(x(i)-x(j))
!       write (19,*) x(i), lspin1(i), lspinorDerivative(i) * ampl(x(i)), test1
!    end do
!    close(19)
!    deallocate(x)
!    deallocate(lspin1)
!    deallocate(lspin2)
!end program test
