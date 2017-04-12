program compute
    use constants
    use matrix
    use lspinors
    use laguerre
    use integrate
    implicit none
    !Тема диплома  L-spinor basis set for relativistic calculations of Highly Charged Ions
    !На русском Использование базиса L-спиноров в релятивистских расчетах электронной структуры многозарядных ионов

    !Длина базиса, по которому раскладывается радиальная часть ВФ
    !Это число совпадает для верхней и нижней компонент
    integer :: N
    !Число узлов квадратурной формулы
    integer :: qLen
    !Хранят узлы и веса квадратурной формулы соответственно.
    !Длина всегда N + 1, чтобы интегрировать перекрестные м. элементы точно
    real(WP), allocatable::x(:),w(:)
    !Счетчики
    integer :: i, j, k
    !Угловой параметр, от него многое зависит
    real(WP)::kappa
    !Зарядовое число. Непонятно, откуда оно взялось, если учесть, что
    !потенциал пока что должен быть произвольным
    real(WP)::Z
    real(WP)::summa

    !Тут начинаются переменные для матричных элементов
    !Перекрестные скалярные произведения L-спиноров (пока что радиальных частей)
    real(WP), allocatable::gram(:,:,:,:)
    !Матричные элементы потенциала
    real(WP), allocatable::V(:,:,:,:)
    !Матричные элементы с 1/x
    real(WP), allocatable::revLen(:,:,:,:)
    !Матричные элементы с производной
    real(WP), allocatable::deriv(:,:,:,:)
    !Гамильтониан и правая часть
    real(WP), allocatable::Ham(:,:), S(:,:), E(:)
    !Тестовые переменные
    real(WP), allocatable::testV(:)
    real(WP), allocatable::testM(:,:)


    !Читаем длину базиса
    print *, "Specify basis length for each component"
    read *, N
    kappa = 1.0_WP
    Z = 1.0_WP
    !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
    !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
    ! 2*qLen - 1. Чтобы интегрировать произведения полиномов Лагерра (степень 2*N),
    !нужно взять такую qLen:
    qLen = N + 1
    call setQuadratureParameters(qLen, 2.0D+00 * calculateGamma(kappa, Z))
    call setLspinorGlobalParameters(N, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
    !allocate(testV(qLen))
    !Проблемы с 0 порядком
    !testV = lspinorVector(0, 'U')
    !print *, testV
    !deallocate (testV)
    !Вдруг понадобятся узлы квадратуры в глобальном контексте в дальнейшем
    allocate (x(qLen))

    allocate(gram(1:2,1:2,N,N))
    allocate(V(1:2,1:2,N,N))
    allocate(deriv(1:2,1:2,N,N))
    allocate(revLen(1:2,1:2,N,N))
    allocate(Ham(2*N,2*N))
    allocate(S(2*N,2*N))
    allocate(E(2*N))
    !Считаем по очереди все необходимые матричные элементы
    do i=1,N
        do j=1,N
        !i = k-1
        !j = l-1
        !Матрица неортогональности
        gram(1,1,i,j) = integrateOnGrid(lspinorVector(i, 'U'),lspinorVector(j, 'U'))
        !gram(1,2,i,j) = integrateOnGrid(lspinorVector(i, 'U'),lspinorVector(j, 'L'))
       ! gram(2,1,i,j) = integrateOnGrid(lspinorVector(i, 'L'),lspinorVector(j, 'U'))
        gram(2,2,i,j) = integrateOnGrid(lspinorVector(i, 'L'),lspinorVector(j, 'L'))
        !Матричный элемент потенциала
        V(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(VPot, j,'U'))
        !V(1,2,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(VPot, j,'L'))
        !V(2,1,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(VPot, j,'U'))
        V(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(VPot, j,'L'))
        !С 1/r
        !revLen(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(oneToX, j,'U'))
        revLen(1,2,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(oneToX, j,'L'))
        revLen(2,1,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(oneToX, j,'U'))
        !revLen(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(oneToX, j,'L'))
        !Производные
        !deriv(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorDerivativeVector(j,'U'))
        deriv(1,2,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorDerivativeVector(j,'L'))
        deriv(2,1,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorDerivativeVector(j,'U'))
        !deriv(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorDerivativeVector(j,'L'))
        end do
    end do

    Ham(1:N,1:N) = c**2 * gram(1,1,:,:) + V(1,1,:,:)
    Ham(1:N,N+1:2*N) = -c * deriv(1,2,:,:) + c * kappa * revLen(1,2,:,:)
    Ham(N+1:2*N,1:N) = c * deriv(2,1,:,:) + c * kappa * revLen(2,1,:,:)
    Ham(N+1:2*N,N+1:2*N) = - c**2 * gram(2,2,:,:) + V(2,2,:,:)
    S(:,:) = 0;
    S(1:N,1:N) = gram(1,1,:,:)
    S(N+1:2*N,N+1:2*N) = gram(2,2,:,:)

    call print_r8(S, 1)

    call general_eigvect_r8(Ham, S, E)

    !call print_r8(Ham)

    E(:) = E(:) - c**2

    do i=1,N
        print *, E(2*N - i)
    end do

    deallocate(Ham)
    deallocate(S)
    deallocate(E)

    print *, "end"
!    print *, gram
    deallocate(gram)
    deallocate(V)
    deallocate(deriv)
    deallocate(revLen)
!    !Проверка ортогональности полиномов Лагерра
!    call setGlobalParameters(N - 1, 0.0_WP, 0.0_WP, qLen, x)
!    do k = 0, N -1
!        do j = 0,N - 1
!            summa=0.0_WP
!            do i=1,qLen
!                summa = summa + w(i)*laguerrePoly(i, k)*laguerrePoly(i, j)
!            end do
!            print *, "j=",j,"k=",k," rez=",summa
!        end do
!        print * , ""
!    end do
    deallocate(x)
end program
