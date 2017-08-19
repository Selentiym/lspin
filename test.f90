!program test
!    use constants
!    use matrix
!    use lspinors
!!    use basis
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
!    degree = 5
!
!
!
!    allocate(x(N))
!    allocate(lspin1(N))
!    allocate(lspin2(N))
!    do i=1,N
!        x(i) = i/300.0D+00 + 1.0D-04
!    end do
!    !call setLspinorGlobalParameters(N, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!    call setLspinorGlobalParameters(degree,1.0D-00, 0.0D+00, N, x)
!    call setLspinorParameters(degree, 'U')
!    open(unit=18, file='test_deriv.dat', status='replace')
!
!    open(unit=19, file='test_lspinor_deriv.dat', status='replace')
!    open(unit=17, file='test_lag_deriv_vec.dat', status='replace')
!    do i=1,N
!        lspin1(i) = laguerrePoly(i, degree)
!        lspin2(i) = laguerrePolyDerivative(i, degree)
!    end do
!!    lspin1 = lagVectorWithFunc(ampl, degree, 'U')
!!    lspin2 = lagDerivativeVectorWithFunc(ampl,degree, 'U')
!!    summa = 0.0D+00
!!    !Проверяем производную полиномов Лагерра
!    do i=2,N
!        j=i-1
!        !test1 = (laguerrePoly(i, degree) - laguerrePoly(i-1, degree))/(x(i)-x(j))! / laguerrePolyDerivative(i, 3)
!        test1 = (lspin1(i) - lspin1(j))/(x(i)-x(j))! / laguerrePolyDerivative(i, 3)
!        !if (abs(test1) > 5) then test1 = 0.0_WP; end if
!        !write (18,*) x(i), laguerrePoly(i, degree), laguerrePolyDerivative(i, degree), test1
!!        write (18,*) x(i), lspin1(i),lspin2(i), ampl(x(i))* &
!!        calcExtendedDeriv(laguerrePoly, laguerrePolyDerivative, degree, degree, i) / sqrt(laguerrePolyNorm(degree)), test1
!        write (18, *) x(i), lspin1(i), lspin2(i), test1, 0.0
!!        test1 = (lspin1(i) - lspin1(j))/(x(i)-x(j))
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
!!    call testLaguerreDerivIntegralsAssymetry()
!    call drawSpinor(0, 'U')
!
!    print *,"end"
!contains
!
!
!    subroutine drawSpinor(iNr, iLetter)
!        integer::iNr, N, i
!        character::iLetter
!        real(WP), allocatable::x(:), lspin1(:),lspin2(:), lspin3(:)
!        real(WP) :: const, testKappa, test
!
!        N = 1000
!        testKappa = -1.0
!        allocate(x(N))
!        allocate(lspin1(N))
!        allocate(lspin2(N))
!        allocate(lspin3(N))
!        do i=1,N
!            x(i) = i/50.0D+00 + 1.0D-04
!        end do
!        !call setLspinorGlobalParameters(N, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!        call setLspinorGlobalParameters(iNr,testKappa, 0.0D+00, N, x)
!        open(unit=18, file='lspinor.dat', status='replace')
!        lspin1 = lspinorVector(iNr, 'L')
!        lspin2 = lspinorDerivativeVector(iNr, 'U')
!        lspin3 = lspinorVector(iNr, 'L')
!        const = 1.0
!        do i=1,N
!!            test = (x(i)/const)**2*exp(-x(i)/2/const)
!            test = (x(i)/const)*exp(-x(i)/2/const)
!            !write(18, *) x(i), lspin3(i)*ampl(x(i)), (lspin2(i) + testKappa/x(i)*lspin1(i))*ampl(x(i))
!            write(18, *) x(i), lspin1(i)*ampl(x(i)), test,0.0_WP!, lspin1(i)*ampl(x(i))/test
!        end do
!
!        close(18)
!        deallocate(x)
!        deallocate(lspin1)
!        deallocate(lspin2)
!        deallocate(lspin3)
!    end subroutine drawSpinor
!
!
!    subroutine testLaguerreDerivIntegralsAssymetry()
!        real(WP), allocatable::poly1(:),poly2(:), derivTest(:,:),derivTest2(:,:), w1(:), x1(:), printMatr1(:,:),printMatr2(:,:)
!        real(WP)::kappa, Z, summa, gammaRel1, test
!        integer::N2, qL, i, j, k, di, dj
!
!        kappa = -1.0D-00
!        Z = 0.0_WP
!        !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
!        !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
!        ! 2*qLen - 1. Чтобы интегрировать произведения полиномов Лагерра (степень 2*N),
!        !нужно взять такую qLen:
!        N2 = 5
!        !qL = N2*150
!        qL = 50
!        !qL = N2+1
!        gammaRel1 = calculateGamma(kappa,Z)
!        !gammaRel1 = 0.0D+00
!        !call setQuadratureParameters(qL, 2.0D+00 * gammaRel1, 1.0D+00)
!        call setQuadratureParameters(qL, 2.0D-00 * gammaRel1, 1.0D+00)
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
!        summa = 0.0D+00
!        do i=1,qL
!            poly1(i) = 1.0D+00
!            poly2(i) = 1.0D+00
!        end do
!        summa = integrateOnGrid(poly2, poly1)
!        print *, summa
!        return
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
!        call setQuadratureParameters(qL, 2.0D+00 * gammaRel1, 1.0D+00)
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
!        call setQuadratureParameters(qL, 2.0D+00 * calculateGamma(kappa,Z), 1.0D+00)
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
!        call setQuadratureParameters(qL, 2.0D+00 * calculateGamma(kappa,Z), 1.0D+00)
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
!    subroutine emulateMainProgram()
!
!        implicit none
!        !Тема диплома  L-spinor basis set for relativistic calculations of Highly Charged Ions
!        !На русском Использование базиса L-спиноров в релятивистских расчетах электронной структуры многозарядных ионов
!
!        !Длина базиса, по которому раскладывается радиальная часть ВФ
!        !Это число совпадает для верхней и нижней компонент
!        integer :: N
!        !Число узлов квадратурной формулы
!        integer :: qLen
!        !Хранят узлы и веса квадратурной формулы соответственно.
!        !Длина всегда N + 1, чтобы интегрировать перекрестные м. элементы точно
!        real(WP), allocatable::x(:),w(:)
!        !Счетчики
!        integer :: i, j, k
!        !Угловой параметр, от него многое зависит
!        real(WP)::kappa
!        !Зарядовое число. Непонятно, откуда оно взялось, если учесть, что
!        !потенциал пока что должен быть произвольным
!        real(WP)::Z
!        real(WP)::summa
!
!        !Тут начинаются переменные для матричных элементов
!        !Перекрестные скалярные произведения L-спиноров (пока что радиальных частей)
!        real(WP), allocatable::gram(:,:,:,:)
!        !Матричные элементы потенциала
!        real(WP), allocatable::V(:,:,:,:)
!        !Матричные элементы с 1/x
!        real(WP), allocatable::revLen(:,:,:,:)
!        !Матричные элементы с производной
!        real(WP), allocatable::deriv(:,:,:,:)
!        !Гамильтониан и правая часть
!        real(WP), allocatable::Ham(:,:), S(:,:), E(:)
!        !Задает тип базиса
!        character::letter
!        !Тестовые переменные
!        real(WP), allocatable::testV(:)
!        real(WP), allocatable::testM(:,:)
!
!
!        !print *, "Select basis type. P - laguerre Polynomials, L - L-spinors"
!        open(unit=18, file='program.ini')
!        read (18,*) letter
!        read(18, *) N
!        read(18, *) kappa
!        read(18, *) Z
!
!        open(unit=19,file="matrix_output.dat",status="replace")
!        write(19,*) ""
!        close(19)
!
!        !kappa = 1.0_WP
!        !Z = 1.0_WP
!        close(18)
!        !read *, letter
!        !Читаем длину базиса
!        !print *, "Specify basis length for each component"
!        !read *, N
!
!        !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
!        !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
!        ! 2*qLen - 1. Чтобы интегрировать произведения полиномов Лагерра (степень 2*N),
!        !нужно взять такую qLen:
!        qLen = 2*N + 2
!        !qLen = max(N+1, 3000)
!        call setQuadratureParameters(qLen, 2.0D+00 * calculateGamma(kappa, Z),0.0D+00)
!        call setLspinorGlobalParameters(N, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
!        !allocate(testV(qLen))
!        !Проблемы с 0 порядком
!        !testV = lspinorVector(0, 'U')
!        !print *, testV
!        !deallocate (testV)
!        !Вдруг понадобятся узлы квадратуры в глобальном контексте в дальнейшем
!        allocate (x(qLen))
!
!        allocate(gram(1:2,1:2,N,N))
!        allocate(V(1:2,1:2,N,N))
!        allocate(deriv(1:2,1:2,N,N))
!        allocate(revLen(1:2,1:2,N,N))
!        allocate(Ham(2*N,2*N))
!        allocate(S(2*N,2*N))
!        allocate(E(2*N))
!
!        if (letter == 'L') then
!            !Считаем все необходимые матричные элементы на L-spinor'ах
!            do i=1,N
!                do j=1,N
!                !i = k-1
!                !j = l-1
!                !Матрица неортогональности
!                gram(1,1,i,j) = integrateOnGrid(lspinorVector(i, 'U'),lspinorVector(j, 'U'))
!                !gram(1,2,i,j) = integrateOnGrid(lspinorVector(i, 'U'),lspinorVector(j, 'L'))
!               ! gram(2,1,i,j) = integrateOnGrid(lspinorVector(i, 'L'),lspinorVector(j, 'U'))
!                gram(2,2,i,j) = integrateOnGrid(lspinorVector(i, 'L'),lspinorVector(j, 'L'))
!                !Матричный элемент потенциала
!                V(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(VPot, j,'U'))
!                !V(1,2,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(VPot, j,'L'))
!                !V(2,1,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(VPot, j,'U'))
!                V(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(VPot, j,'L'))
!                !С 1/r
!                !revLen(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(oneToX, j,'U'))
!                revLen(1,2,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(oneToX, j,'L'))
!                revLen(2,1,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(oneToX, j,'U'))
!                !revLen(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(oneToX, j,'L'))
!                !Производные
!                deriv(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorDerivativeVector(j,'U'))
!                deriv(1,2,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorDerivativeVector(j,'L'))
!                deriv(2,1,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorDerivativeVector(j,'U'))
!                !deriv(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorDerivativeVector(j,'L'))
!                end do
!            end do
!        else
!            !Считаем все необходимые матричные элементы на полиномах Лагерра
!            do i=1,N
!                do j=1,N
!                !i = k-1
!                !j = l-1
!                !Матрица неортогональности
!                gram(1,1,i,j) = integrateOnGrid(lagVector(i, 'U'),lagVector(j, 'U'))
!                !gram(1,2,i,j) = integrateOnGrid(lagVector(i, 'U'),lagVector(j, 'L'))
!               ! gram(2,1,i,j) = integrateOnGrid(lagVector(i, 'L'),lagVector(j, 'U'))
!                gram(2,2,i,j) = integrateOnGrid(lagVector(i, 'L'),lagVector(j, 'L'))
!                !Матричный элемент потенциала
!                V(1,1,i,j) = integrateOnGrid(lagVector(i,'U'), lagVectorWithFunc(VPot, j,'U'))
!                !V(1,2,i,j) = integrateOnGrid(lagVector(i,'U'), lagVectorWithFunc(VPot, j,'L'))
!                !V(2,1,i,j) = integrateOnGrid(lagVector(i,'L'), lagVectorWithFunc(VPot, j,'U'))
!                V(2,2,i,j) = integrateOnGrid(lagVector(i,'L'), lagVectorWithFunc(VPot, j,'L'))
!                !С 1/r
!                !revLen(1,1,i,j) = integrateOnGrid(lagVector(i,'U'), lagVectorWithFunc(oneToX, j,'U'))
!                revLen(1,2,i,j) = integrateOnGrid(lagVector(i,'U'), lagVectorWithFunc(oneToX, j,'L'))
!                revLen(2,1,i,j) = integrateOnGrid(lagVector(i,'L'), lagVectorWithFunc(oneToX, j,'U'))
!                !revLen(2,2,i,j) = integrateOnGrid(lagVector(i,'L'), lagVectorWithFunc(oneToX, j,'L'))
!                !Производные
!                !deriv(1,1,i,j) = integrateOnGrid(lagVector(i,'U'), lagDerivativeVector(j,'U'))
!                deriv(1,2,i,j) = integrateOnGrid(lagVector(i,'U'), lagDerivativeVector(j,'L'))
!                deriv(2,1,i,j) = integrateOnGrid(lagVector(i,'L'), lagDerivativeVector(j,'U'))
!                !deriv(2,2,i,j) = integrateOnGrid(lagVector(i,'L'), lagDerivativeVector(j,'L'))
!                end do
!            end do
!        end if
!        !Костыль, не должно быть этого здесь
!    !    if (letter == 'L') then
!    !        do i=1,N
!    !            do j=1,N
!    !                deriv(1,2,i,j) = (deriv(1,2,i,j) - deriv(2,1,i,j))/2
!    !                deriv(2,1,j,i) = - deriv(1,2,i,j)
!    !            end do
!    !        end do
!    !    else
!    !        do i=1,N
!    !            do j=1,N
!    !                deriv(1,2,i,j) = (deriv(1,2,i,j) - deriv(1,2,j,i))/2
!    !                deriv(2,1,j,i) = deriv(1,2,i,j)
!    !            end do
!    !        end do
!    !    end if
!        !Ham(1:N,1:N) = c**2 * gram(1,1,:,:) + V(1,1,:,:)
!        Ham(1:N,1:N) = V(1,1,:,:)
!        Ham(1:N,N+1:2*N) = -c * deriv(1,2,:,:) + c * kappa * revLen(1,2,:,:)
!        Ham(N+1:2*N,1:N) = c * deriv(2,1,:,:) + c * kappa * revLen(2,1,:,:)
!        Ham(N+1:2*N,N+1:2*N) = - 2*c**2 * gram(2,2,:,:) + V(2,2,:,:)
!        S(:,:) = 0;
!        S(1:N,1:N) = gram(1,1,:,:)
!        S(N+1:2*N,N+1:2*N) = gram(2,2,:,:)
!        print *, "Derivatives, UL"
!        call print4dimAs2(1,2,deriv)
!        !print *, "Derivatives, UU"
!        !call print4dimAs2(1,1,deriv)
!        print *, "Derivatives, LU"
!        call print4dimAs2(2,1,deriv)
!        print *, "Hamiltonian"
!        call print_r8(Ham, 1)
!        print *,"S"
!        call print_r8(S, 1)
!
!        call general_eigvect_r8(Ham, S, E)
!
!        !call print_r8(Ham)
!
!        E(:) = E(:)
!
!        do i=2*N,N+1,-1
!            print *, "E=",E(i)," n=",Z/sqrt(-2*E(i))
!        end do
!
!        deallocate(Ham)
!        deallocate(S)
!        deallocate(E)
!
!        print *, "end"
!    !    print *, gram
!        deallocate(gram)
!        deallocate(V)
!        deallocate(deriv)
!        deallocate(revLen)
!    !    !Проверка ортогональности полиномов Лагерра
!    !    call setGlobalParameters(N - 1, 0.0_WP, 0.0_WP, qLen, x)
!    !    do k = 0, N -1
!    !        do j = 0,N - 1
!    !            summa=0.0_WP
!    !            do i=1,qLen
!    !                summa = summa + w(i)*laguerrePoly(i, k)*laguerrePoly(i, j)
!    !            end do
!    !            print *, "j=",j,"k=",k," rez=",summa
!    !        end do
!    !        print * , ""
!    !    end do
!        deallocate(x)
!
!
!
!    end subroutine emulateMainProgram
!
!end program test
