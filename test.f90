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
!    real(WP)::test1, summa
!
!    N = 700
!    degree = 51
!    allocate(x(N))
!    allocate(lspin1(N))
!    allocate(lspin2(N))
!    do i=1,N
!        x(i) = i**(1.2)/10.0D+00
!    end do
!    !call setLspinorGlobalParameters(N, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!    call setLspinorGlobalParameters(degree,1.0D-00, 0.0D+00, N, x)
!    call setLspinorParameters(degree, 'U')
!    open(unit=18, file='test_deriv.dat', status='replace')
!
!    open(unit=19, file='test_lspinor_deriv.dat', status='replace')
!    open(unit=17, file='test_lag_deriv_vec.dat', status='replace')
!    lspin1 = lagVectorWithFunc(ampl, degree, 'U')
!    lspin2 = lagDerivativeVectorWithFunc(ampl,degree, 'U')
!    summa = 0.0D+00
!    !Проверяем производную полиномов Лагерра
!    do i=2,N
!        j=i-1
!        !test1 = (laguerrePoly(i, degree) - laguerrePoly(i-1, degree))/(x(i)-x(j))! / laguerrePolyDerivative(i, 3)
!        test1 = (lspin1(i) - lspin1(j))/(x(i)-x(j))! / laguerrePolyDerivative(i, 3)
!        !if (abs(test1) > 5) then test1 = 0.0_WP; end if
!
!        !write (18,*) x(i), laguerrePoly(i, degree), laguerrePolyDerivative(i, degree), test1
!        write (18,*) x(i), lspin1(i),lspin2(i), ampl(x(i))* &
!        calcExtendedDeriv(laguerrePoly, laguerrePolyDerivative, degree, degree, i) / sqrt(laguerrePolyNorm(degree)), test1
!
!        test1 = (lspin1(i) - lspin1(j))/(x(i)-x(j))
!        write (17,*) x(i), lspin1(i), lspin2(i), test1
!        summa = summa + lspin1(i) * lspin2(i) * (x(i)-x(j))
!    end do
!    !print *, "Lag integral=",summa
!    close(18)
!    close(17)
!    degree = 1
!    letter = 'U'
!    !lspin1 = lspinorVector(degree, letter)
!
!    lspin1 = lspinorVectorWithFunc(ampl, degree, letter)
!    lspin2 = lspinorDerivativeVectorWithFunc(ampl, degree, letter)
!    !lspin2 = lspinorDerivativeVector(degree, letter)
!
!    do i=2,N
!       j = i - 1
!       test1 = (lspin1(i) - lspin1(j))/(x(i)-x(j))
!       write (19,*) x(i), lspin1(i), lspin2(i), test1
!    end do
!    close(19)
!    deallocate(x)
!    deallocate(lspin1)
!    deallocate(lspin2)
!
!    !call testLaguerreDerivIntegralsAssymetry()
!    call testLspinorDerivIntegralsAssymetry()
!
!contains
!
!    subroutine testLaguerreDerivIntegralsAssymetry()
!        real(WP), allocatable::poly1(:),poly2(:), derivTest(:,:),derivTest2(:,:), w1(:), x1(:), printMatr1(:,:),printMatr2(:,:)
!        real(WP)::kappa, Z, summa, gammaRel1, test
!        integer::N2, qL, i, j, k, di, dj
!
!        kappa = 1.0D-00
!        Z = 0.0_WP
!        !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
!        !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
!        ! 2*qLen - 1. Чтобы интегрировать произведения полиномов Лагерра (степень 2*N),
!        !нужно взять такую qLen:
!        N2 = 5
!        !qL = N2*150
!        qL = 300
!        gammaRel1 = calculateGamma(kappa,Z)
!        !gammaRel1 = 0.0D+00
!        call setQuadratureParameters(qL, 2.0D+00 * gammaRel1)
!        call setLspinorGlobalParameters(N2+1, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!        allocate (derivTest(0:N2,0:N2))
!        allocate (derivTest2(0:N2,0:N2))
!        allocate (poly1(qL))
!        allocate (poly2(qL))
!        allocate (w1(qL))
!        allocate (x1(qL))
!        w1 = getIntegrateWeights()
!        x1 = getIntegrateGrid()
!
!!        call setWriteToFile("LookOnPolynoms.dat")
!!        allocate(printMatr1(0:N2,0:qL-1))
!!        allocate(printMatr2(0:N2,0:qL-1))
!!        do i=0,N2
!!            lspin1 = lagVectorWithFunc(ampl,i, 'U')
!!            lspin2 = lagDerivativeVector(i, 'U')
!!            printMatr1(i,0:qL-1) = lspin1
!!            printMatr2(i,0:qL-1) = lspin2
!!        end do
!!
!!        call print_r8(printMatr1, 0)
!!        print *, "test"
!!        call print_r8(printMatr2, 0)
!!        print *, "test"
!!        deallocate(printMatr1)
!!        deallocate(printMatr2)
!
!        do i=0,N2
!            do j=0,N2
!                summa = 0.0D+00
!                lspin1 = lagVectorWithFunc(unityFunc,i, 'U')
!                lspin2 = lagDerivativeVector(j, 'U')
!                do k=1,qL
!                    test = lspin1(k) * lspin2(k) * w1(k)
!                    !summa = summa + test
!!                    summa = summa + w1(k) * &
!!                    laguerrePoly(k, i) * ( laguerrePolyDerivative(k, j) &
!!                    + gammaRel1 * laguerrePoly(k, j) / x1(k) &
!!                    - 1/2 * laguerrePoly(k, j) )
!!                    test = laguerrePoly(k, i)
!!                    test = 0.0D+00
!!                    test = laguerrePolyDerivative(k, j)
!!                    test = 0.0D+00
!!                    test = laguerrePoly(k, j)
!!                    test = 0.0D+00
!!                    test = w1(k) * &
!!                    laguerrePoly(k, i) * ( laguerrePolyDerivative(k, j) &
!!                    + gammaRel1 * laguerrePoly(k, j) / x1(k) &
!!                    - 1/2 * laguerrePoly(k, j) )
!
!                    !summa = summa + power(k, i) * calcExtendedDeriv(power, powerDeriv,j,j, k) * w1(k)
!
!                    summa = summa + laguerrePoly(k, i) * &
!                    calcExtendedDeriv(laguerrePoly,laguerrePolyDerivative,j,j,k)/&
!                    sqrt(laguerrePolyNorm(j)*laguerrePolyNorm(i)) * w1(k)
!                    !summa = summa + x1(k) ** i *w1(k)
!                    !summa = summa + x1(k)**i * ((j + gammaRel1)*x1(k)**(j-1) -0.5_WP*x1(k)**j) * w1(k)
!                end do
!
!                derivTest(i,j) = summa
!                derivTest2(i,j) = integrateOnGrid(lspin1,lspin2)
!            end do
!        end do
!
!        print *, "derivative and lag poly matrix"
!        call print_r8(derivTest, 0)
!        print *, power(1,4)
!        call print_r8(derivTest2, 0)
!
!        deallocate (poly1)
!        deallocate (poly2)
!        deallocate (derivTest)
!        deallocate (derivTest2)
!        deallocate (w1)
!        deallocate (x1)
!
!    end subroutine testLaguerreDerivIntegralsAssymetry
!
!    subroutine testLspinorDerivIntegralsAssymetry()
!        real(WP), allocatable::lspin1(:),lspin2(:), derivTest(:,:),derivTest2(:,:), w1(:), x1(:), printMatr1(:,:),printMatr2(:,:)
!        real(WP)::kappa, Z, summa, gammaRel1, test
!        integer::N2, qL, i, j, k, di, dj
!
!        kappa = 1.0D-00
!        Z = 0.0_WP
!        !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
!        !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
!        ! 2*qLen - 1. Чтобы интегрировать произведения полиномов Лагерра (степень 2*N),
!        !нужно взять такую qLen:
!        N2 = 5
!        !qL = N2*150
!        qL = 300
!        gammaRel1 = calculateGamma(kappa,Z)
!        !gammaRel1 = 0.0D+00
!        call setQuadratureParameters(qL, 2.0D+00 * gammaRel1)
!        call setLspinorGlobalParameters(N2+1, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!        allocate (derivTest(0:N2,0:N2))
!        allocate (derivTest2(0:N2,0:N2))
!        allocate (lspin1(qL))
!        allocate (lspin2(qL))
!        allocate (w1(qL))
!        allocate (x1(qL))
!        w1 = getIntegrateWeights()
!        x1 = getIntegrateGrid()
!
!        do i=0,N2
!            do j=0,N2
!                lspin1 = lspinorVector(i, 'U')
!                lspin2 = lspinorDerivativeVector(j, 'L')
!
!                derivTest(i,j) = integrateOnGrid(lspin1,lspin2)
!                lspin1 = lspinorVector(i, 'L')
!                lspin2 = lspinorDerivativeVector(j, 'U')
!                derivTest2(i,j) = integrateOnGrid(lspin1,lspin2)
!!                do k=1,qL
!!                    test = lspin1(k) * lspin2(k) * w1(k)
!!                    !summa = summa + test
!!!                    summa = summa + w1(k) * &
!!!                    laguerrePoly(k, i) * ( laguerrePolyDerivative(k, j) &
!!!                    + gammaRel1 * laguerrePoly(k, j) / x1(k) &
!!!                    - 1/2 * laguerrePoly(k, j) )
!!!                    test = laguerrePoly(k, i)
!!!                    test = 0.0D+00
!!!                    test = laguerrePolyDerivative(k, j)
!!!                    test = 0.0D+00
!!!                    test = laguerrePoly(k, j)
!!!                    test = 0.0D+00
!!!                    test = w1(k) * &
!!!                    laguerrePoly(k, i) * ( laguerrePolyDerivative(k, j) &
!!!                    + gammaRel1 * laguerrePoly(k, j) / x1(k) &
!!!                    - 1/2 * laguerrePoly(k, j) )
!!
!!                    !summa = summa + power(k, i) * calcExtendedDeriv(power, powerDeriv,j,j, k) * w1(k)
!!
!!                    summa = summa + laguerrePoly(k, i) * &
!!                    calcExtendedDeriv(laguerrePoly,laguerrePolyDerivative,j,j,k)/&
!!                    sqrt(laguerrePolyNorm(j)*laguerrePolyNorm(i)) * w1(k)
!!                    !summa = summa + x1(k) ** i *w1(k)
!!                    !summa = summa + x1(k)**i * ((j + gammaRel1)*x1(k)**(j-1) -0.5_WP*x1(k)**j) * w1(k)
!!                end do
!
!
!            end do
!        end do
!
!
!        print *, "derivative and lspinor matrix UL"
!        call print_r8(derivTest, 0)
!        print *, "derivative and lspinor matrix LU"
!        call print_r8(derivTest2, 0)
!        print *, "have to be ULt = -LU"
!
!        deallocate (lspin1)
!        deallocate (lspin2)
!        deallocate (derivTest)
!        deallocate (derivTest2)
!        deallocate (w1)
!        deallocate (x1)
!
!    end subroutine testLspinorDerivIntegralsAssymetry
!
!    subroutine testLaguerreOrthog()
!        real(WP), allocatable::poly1(:),poly2(:), gramTest(:,:), w1(:)
!        real(WP)::kappa, Z, summa
!        integer::N2, qL, i, k, di, dj
!
!        kappa = 1.0D+00
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
!                    !summa = summa + poly1(i) * poly2(i) * w1(i)
!                end do
!                summa = integrateOnGrid(poly1, poly2)
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
