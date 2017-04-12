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
!    lspin2 = lspinorDerivativeVector(degree, letter)
!
!    do i=2,N
!       j = i - 1
!       test1 = (lspin1(i) - lspin1(j))/(x(i)-x(j))
!       write (19,*) x(i), lspin1(i), lspin2(i) * ampl(x(i)), test1
!    end do
!    close(19)
!    deallocate(x)
!    deallocate(lspin1)
!    deallocate(lspin2)
!
!    call testLspinorOrthogonality()
!
!contains
!
!    subroutine testLaguerreOrthog()
!        real(WP), allocatable::poly1(:),poly2(:), gramTest(:,:), w1(:)
!        real(WP)::kappa, Z, summa
!        integer::N2, qL, i, k, di, dj
!
!        kappa = 0.1D-3
!        Z = 0.0_WP
!        !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
!        !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
!        ! 2*qLen - 1. Чтобы интегрировать произведения полиномов Лагерра (степень 2*N),
!        !нужно взять такую qLen:
!        N2 = 20
!        qL = N2+1
!        call setQuadratureParameters(qL, 2.0D+00 * calculateGamma(kappa,Z))
!        call setLspinorGlobalParameters(N2+1, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!        allocate (gramTest(0:N2,0:N2))
!        allocate (poly1(qL))
!        allocate (poly2(qL))
!
!        w1 = getIntegrateWeights()
!
!        gramTest(:,:) = 0.0D+00
!
!        do di=0,N2
!            do dj=0,N2
!                summa = 0.0D+00
!                do i=1,qL
!                    poly1(i) = laguerrePoly(i, di)
!                    poly2(i) = laguerrePoly(i, dj)
!                    summa = summa + poly1(i) * poly2(i) * w1(i)
!                end do
!                write(*, "(F20.7)" , advance="no") summa
!                gramTest(dj,di) = summa
!            end do
!            print *, ""
!        end do
!
!        print *, ""
!
!        summa = 0.0D+00
!        do i=1,qL
!            poly1(i) = laguerrePoly(i, N2)
!            poly2(i) = laguerrePoly(i, N2)
!            summa = summa + poly1(i)*poly2(i)*w1(i)
!        end do
!
!        call print_r8(gramTest, 0)
!        deallocate (gramTest, poly1,poly2)
!    end subroutine testLaguerreOrthog
!
!    subroutine testLspinorOrthogonality()
!        real(WP), allocatable::lspin1(:),lspin2(:), gramTest(:,:), w1(:)
!        real(WP)::kappa, Z, summa, saveCoeff
!        integer::N2, qL, i, k, di, dj
!
!        kappa = 0.1D-3
!        Z = 0.0_WP
!        !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
!        !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
!        ! 2*qLen - 1. Чтобы интегрировать произведения полиномов Лагерра (степень 2*N),
!        !нужно взять такую qLen:
!        N2 = 20
!        qL = N2+1
!        call setQuadratureParameters(qL, 2.0D+00 * calculateGamma(kappa,Z))
!        call setLspinorGlobalParameters(N2, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!
!        allocate(lspin1(qL))
!        allocate(lspin2(qL))
!        lspin1 = lspinorVector(10, 'L')
!        saveCoeff = -coeff1
!        lspin2 = lspinorVector(9, 'U')
!        summa = integrateOnGrid(lspin1, lspin2)
!
!        print *, summa, coeff2*saveCoeff
!
!        deallocate(lspin1)
!        deallocate(lspin2)
!
!    end subroutine testLspinorOrthogonality
!
!
!
!end program test
