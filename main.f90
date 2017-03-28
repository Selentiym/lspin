program compute
    use constants
    use matrix
    use lspinors
    use laguerre

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
    !Счетчики
    integer :: i, j, k
    !Хранят узлы и веса квадратурной формулы соответственно.
    !Длина всегда N + 1, чтобы интегрировать перекрестные м. элементы точно
    real(WP), allocatable::x(:),w(:)
    !Угловой параметр, от него многое зависит
    real(WP)::kappa

    real(WP)::summa

    !Тут начинаются переменные для матричных элементов
    !Перекрестные скалярные произведения L-спиноров (пока что радиальных частей)
    real(WP), allocatable::gram(:,:)
    !Матричные элементы потенциала
    real(WP), allocatable::V(:,:)


    !Задаем длину базиса, потом будет как-нибудь вводиться
    N = 5
    !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
    !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
    ! 2*qLen - 1. Чтобы интегрировать произведения полиномов Лагерра (степень 2*N),
    !нужно взять такую qLen:
    qLen = N + 1
    !allocate(x(1))
    !x(1) = 0
    !qLen = 1
    allocate (x(qLen), w(qLen))
    !Получаем узлы с весами
    call l_quadrature_rule(qLen, x, w)

    !Считаем по очереди все необходимые матричные элементы
!    do i=1,N
!        do j=1,N
!        gram(i,j) = 1
!        end do
!    end do

    call setGlobalParameters(N - 1, 0.0_WP, 0.0_WP, qLen, x)



    do k = 0, N -1
    do j = 0,N - 1
    summa=0.0_WP
    do i=1,qLen
        summa = summa + w(i)*laguerrePoly(i, k)*laguerrePoly(i, j)
    end do
    print *, "j=",j,"k=",k," rez=",summa
    end do
    print * , ""
    end do
    deallocate(x,w)

end program
