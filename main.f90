program compute
    use constants
    use matrix
    use lspinors
    use laguerre
    use integrate
    implicit none
    !Тема диплома  L-spinor basis set for relativistic calculations of Highly Charged Ions
    !На русском Использование базиса L-спиноров в релятивистских расчетах электронной структуры многозарядных ионов

!    integer, parameter::N=5
!    complex::A(N,N),lworkINIT(1),lrworkINIT(1),liworkINIT(1)
!    real::OUTp(N),rwork(3*N-2)
!    complex, allocatable :: work(:)
!    integer::i,j,lwork,liwork,lrwork
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
    !Гамильтониан и правая часть
    real(WP), allocatable::Ham(:,:), S(:,:)
    !Тестовые переменные
    real(WP), allocatable::testV(:)
    real(WP), allocatable::testM(:,:)


    !Задаем длину базиса, потом будет как-нибудь вводиться
    N = 10
    kappa = 1.0_WP
    Z = 1.0_WP
    !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
    !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
    ! 2*qLen - 1. Чтобы интегрировать произведения полиномов Лагерра (степень 2*N),
    !нужно взять такую qLen:
    qLen = N + 1
    call setQuadratureParameters(qLen)
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
    allocate(revLen(1:2,1:2,N,N))
    allocate(Ham(2*N,2*N))
    allocate(S(2*N,2*N))
    !Считаем по очереди все необходимые матричные элементы
    do i=1,N
        do j=1,N
        !Матрица неортогональности
        gram(1,1,i,j) = integrateOnGrid(lspinorVector(i, 'U'),lspinorVector(j, 'U'))
        gram(1,2,i,j) = integrateOnGrid(lspinorVector(i, 'U'),lspinorVector(j, 'L'))
        gram(2,1,i,j) = integrateOnGrid(lspinorVector(i, 'L'),lspinorVector(j, 'U'))
        gram(2,2,i,j) = integrateOnGrid(lspinorVector(i, 'L'),lspinorVector(j, 'L'))
        !Матричный элемент потенциала
        V(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(VPot, j,'U'))
        V(1,2,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(VPot, j,'L'))
        V(2,1,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(VPot, j,'U'))
        V(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(VPot, j,'L'))
        !С 1/r
        revLen(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(oneToX, j,'U'))
        revLen(1,2,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(oneToX, j,'L'))
        revLen(2,1,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(oneToX, j,'U'))
        revLen(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(oneToX, j,'L'))
        !Производные
        end do
    end do

    deallocate(Ham)
    deallocate(S)

    print *, "end"
!    print *, gram
    deallocate(gram)
    deallocate(V)
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
