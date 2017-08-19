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
    !Номер первой базисной функции.
    integer :: offset
    !Число узлов квадратурной формулы
    integer :: qLen
    !Хранят узлы и веса квадратурной формулы соответственно.
    !Длина всегда N + 1, чтобы интегрировать перекрестные м. элементы точно
    real(WP), allocatable::x(:),w(:)
    !Счетчики
    integer :: i, j, k, l, il, jl, minQuant, maxLevel, quantN, quantL
    !Квантовые числа. Пока что главное число - kappa
    real(WP):: quantJ, kappa
    !Зарядовое число. Непонятно, откуда оно взялось, если учесть, что
    !потенциал пока что должен быть произвольным
    real(WP)::Z
    !Всякий мусор сюда можно кидать, а также какое-нибудь трапециевое интегрирование
    real(WP)::summa
    !Параметр растяжения по координате. Предполагается, что радиальная часть
    !волновой функции раскладывается по L(xScale * x), в статье было
    !xScale = 2*lambda, однако мне удобнее вынести соответствующие множители в
    real(WP)::xScale
    !Сдвиг степени квадратуры. Смотри модуль интегрирования за более подробными комментариями
    real(WP)::qOffset
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
    !Задает тип базиса и для вывода буквенное выражение числа L
    character::letter, orbitalLetter
    !Для подсчета радиусов и других средних
    real(WP), allocatable::f1(:,:), f2(:,:)
    !Тестовые переменные
    real(WP), allocatable::testV(:)
    real(WP), allocatable::testM(:,:)

!***************************************************************************
!program starts
!***************************************************************************

    !print *, "Select basis type. P - laguerre Polynomials, L - L-spinors"
    open(unit=18, file='program.ini')
    read (18,*) letter
    read(18, *) N
    read(18, *) kappa
    read(18, *) Z
    read(18, *) xScale
    read(18, *) qOffset
    read(18, *) offset
    close(18)

    if (xScale <= 0.0_WP) then
        xScale = real(Z,8)
    end if


    if (offset <= 0) then
        if (kappa > 0.0_WP) then
            offset = 1
        else
            offset = 0
        end if
    end if

    quantJ = abs(kappa) - 0.5_WP

    open(unit=19,file="matrix_output.dat",status="replace")
    write(19,*) ""
    close(19)

    !Вычисляем число узлов квадратуры и выделяем память на узлы и веса
    !алгебраическая точность квадратурной формулы Гусса-Лагерра равна
    ! 2*qLen - 1. Мы будем интегрировать произведения полиномов Лагерра (степень 2*N)
    !Точность берем с запасом
    qLen = 2*N + 2

    !Матрица производных не получается антисимметричной, так как мы считаем производную
    !именно от базисной функции, а там дифференцируется X**alpha, который должен уходить в
    !формулу квадратуры, а после дифференцирования получается alpha*X**alpha / X,
    !который как раз и дает проблему!
    !Но если ввести сдвиг квадратуры, см qOffset, то все ок
    call setQuadratureParameters(qLen, 2.0D+00 * calculateGamma(kappa, Z),qOffset)
    call setQuadratureArgumentScale(xScale)
    call setBasisScale(xScale)
    !Поскольку есть сдвиг базиса, то понадобится спинор более высокого порядка, чем N.
    !Еще одна единица добавлена на всякий случай
    call setLspinorGlobalParameters(N+offset+1, kappa, Z, getIntegrateGridLength(), getIntegrateGrid())
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

    if (letter == 'L') then
        !Считаем все необходимые матричные элементы на L-spinor'ах
        do k=1,N
            do l=1,N
            i = k+offset - 1
            j = l+offset - 1
            il = k+offset - 1
            jl = l+offset - 1
            !Матрица неортогональности
            gram(1,1,k,l) = integrateOnGrid(lspinorVector(i, 'U'),lspinorVector(j, 'U'))
            !gram(1,2,i,j) = integrateOnGrid(lspinorVector(i, 'U'),lspinorVector(j, 'L'))
           ! gram(2,1,i,j) = integrateOnGrid(lspinorVector(i, 'L'),lspinorVector(j, 'U'))
            gram(2,2,k,l) = integrateOnGrid(lspinorVector(il, 'L'),lspinorVector(jl, 'L'))
            !Матричный элемент потенциала
            V(1,1,k,l) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(VPot, j,'U'))
            !V(1,2,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(VPot, j,'L'))
            !V(2,1,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(VPot, j,'U'))
            V(2,2,k,l) = integrateOnGrid(lspinorVector(il,'L'), lspinorVectorWithFunc(VPot, jl,'L'))
            !С 1/r
            !revLen(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(oneToX, j,'U'))
            revLen(1,2,k,l) = integrateOnGrid(lspinorVector(i,'U'), lspinorVectorWithFunc(oneToX, jl,'L'))
            revLen(2,1,k,l) = integrateOnGrid(lspinorVector(il,'L'), lspinorVectorWithFunc(oneToX, j,'U'))
            !revLen(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorVectorWithFunc(oneToX, j,'L'))
            !Производные
            !deriv(1,1,i,j) = integrateOnGrid(lspinorVector(i,'U'), lspinorDerivativeVector(j,'U'))
            deriv(1,2,k,l) = integrateOnGrid(lspinorVector(i,'U'), lspinorDerivativeVector(jl,'L'))
            deriv(2,1,k,l) = integrateOnGrid(lspinorVector(il,'L'), lspinorDerivativeVector(j,'U'))
            !deriv(2,2,i,j) = integrateOnGrid(lspinorVector(i,'L'), lspinorDerivativeVector(j,'L'))
            end do
        end do
    else
        !Считаем все необходимые матричные элементы на полиномах Лагерра
        do k=1,N
            do l=1,N
            i = k+offset - 1
            j = l+offset - 1
            il = k+offset - 1
            jl = l+offset - 1
            !Матрица неортогональности
            gram(1,1,k,l) = integrateOnGrid(lagVector(i, 'U'),lagVector(j, 'U'))
            !gram(1,2,i,j) = integrateOnGrid(lagVector(i, 'U'),lagVector(j, 'L'))
           ! gram(2,1,i,j) = integrateOnGrid(lagVector(i, 'L'),lagVector(j, 'U'))
            gram(2,2,k,l) = integrateOnGrid(lagVector(il, 'L'),lagVector(jl, 'L'))
            !Матричный элемент потенциала
            V(1,1,k,l) = integrateOnGrid(lagVector(i,'U'), lagVectorWithFunc(VPot, j,'U'))
            !V(1,2,i,j) = integrateOnGrid(lagVector(i,'U'), lagVectorWithFunc(VPot, j,'L'))
            !V(2,1,i,j) = integrateOnGrid(lagVector(i,'L'), lagVectorWithFunc(VPot, j,'U'))
            V(2,2,k,l) = integrateOnGrid(lagVector(il,'L'), lagVectorWithFunc(VPot, jl,'L'))
            !С 1/r
            !revLen(1,1,i,j) = integrateOnGrid(lagVector(i,'U'), lagVectorWithFunc(oneToX, j,'U'))
            revLen(1,2,k,l) = integrateOnGrid(lagVector(i,'U'), lagVectorWithFunc(oneToX, jl,'L'))
            revLen(2,1,k,l) = integrateOnGrid(lagVector(il,'L'), lagVectorWithFunc(oneToX, j,'U'))
            !revLen(2,2,i,j) = integrateOnGrid(lagVector(i,'L'), lagVectorWithFunc(oneToX, j,'L'))
            !Производные
            !deriv(1,1,i,j) = integrateOnGrid(lagVector(i,'U'), lagDerivativeVector(j,'U'))
            deriv(1,2,k,l) = integrateOnGrid(lagVector(i,'U'), lagDerivativeVector(jl,'L'))
            deriv(2,1,k,l) = integrateOnGrid(lagVector(il,'L'), lagDerivativeVector(j,'U'))
            !deriv(2,2,i,j) = integrateOnGrid(lagVector(i,'L'), lagDerivativeVector(j,'L'))
            end do
        end do
    end if
    !Считаем матрицу гамильтониана и неортогональности, используя готовые
    !матричные элементы
    Ham(1:N,1:N) = V(1,1,:,:)
    Ham(1:N,N+1:2*N) = -c * xScale *deriv(1,2,:,:) + c * kappa * revLen(1,2,:,:)
    Ham(N+1:2*N,1:N) = +c * xScale * deriv(2,1,:,:) + c * kappa * revLen(2,1,:,:)
    Ham(N+1:2*N,N+1:2*N) = - 2*c**2 * gram(2,2,:,:) + V(2,2,:,:)

    S(:,:) = 0;
    S(1:N,1:N) = gram(1,1,:,:)
    S(N+1:2*N,N+1:2*N) = gram(2,2,:,:)
    !Если нужно, можно раскомментировать для печати матриц на экран
!    print *, "Derivatives, UL"
!    call print4dimAs2(1,2,deriv)
!    !print *, "Derivatives, UU"
!    !call print4dimAs2(1,1,deriv)
!    print *, "Derivatives, LU"
!    call print4dimAs2(2,1,deriv)
!    print *, "Hamiltonian"
!    call print_r8(Ham, 1)
!    print *,"S"
!    call print_r8(S, 1)

    !Тут происходит основная работа, используется lapack'овская процедура
    !для решения обобщенной задачи на СЧ и СВ. general_eigvect_r8 - моя
    !оболочка для сложной лапаковской процедуры, лежит в matrix.f90
    !В E лежат СЧ, в Ham лежат СВ
    call general_eigvect_r8(Ham, S, E)

    !call print_r8(Ham)

    !Переводим в знакомые нерелятивистские числа
    minQuant = floor(abs(kappa))
    if (kappa > 0) then
        minQuant = minQuant + 1
    end if
    quantL = minQuant - 1

    !Далее отображательный код, я его всюду закомментил
!    open(unit=19,file=letter//".rez",status="replace")
!    write(19,*) "kappa=",kappa," scale=",xScale, " Z=",Z, " basis Length=",N," basis type=",letter
!
!    SELECT CASE (quantL)
!    CASE (0)
!        orbitalLetter = 's'
!    CASE (1)
!        orbitalLetter = 'p'
!    CASE (2)
!        orbitalLetter = 'd'
!    CASE (3)
!        orbitalLetter = 'f'
!    CASE DEFAULT
!        orbitalLetter = 'U'
!    END SELECT
    !call print_r8(Ham, 1)
100 format (F20.13)
    maxLevel = 10
    summa = 0.0_WP
    !Тут считаем радиусы: получается не очень для Z=92 и kappa=+1.
    !Ниже будет более подробный пример расчета значений волновой функции
    do quantN = minQuant, maxLevel
        i = quantN - minQuant + N + 1
        !write (*, 100) rightEnergy(quantN, quantJ)
!        print *, quantN, ' ', E(i), rightEnergy(quantN, quantJ), E(i) - rightEnergy(quantN, quantJ),":"
        !глобальная переменная, используется в функции powerPreset
        powerCached = 1.0_WP
        allocate(f1(1:2,2*N))
        allocate(f2(1:2,2*N))
        f1 = calcEigFunc(Ham, i, offset)
        f2 = calcEigFuncWithFunc(Ham, i, powerPreset, offset)
        summa = integrateOnGrid(f1(1,:), f2(1,:)) + integrateOnGrid(f1(2,:), f2(2,:))
        print *, summa/(integrateOnGrid(f1(1,:), f1(1,:)) + integrateOnGrid(f1(2,:), f1(2,:))), rightRadius(quantN, quantL, quantJ)
        deallocate(f1)
        deallocate(f2)
!        write (*, 100) E(i)
    end do


    !Теперь хочу посчитать какой-нибудь интеграл на сетке или выгрузить базисные функции
    !Для этого нужно пересчитать lspinor'ы (или полиномы лаггера) в новых точках.
    !1 аргумент - порядок максимально необходимого полинома лагерра - остается неизменным
    !так как зависит от размера решенной ранее задачи. kappa, Z также неизменны,
    !меняем только сетку
    maxLevel = 20000 !используем переменную не по назначению, будет хранить длину сетки
    allocate(testV(maxLevel))
    do i=1,maxLevel
        testV(i)=i/1000.0_WP
    end do
    call setLspinorGlobalParameters(N+offset+1, kappa, Z, maxLevel, testV)
    !Теперь все вызовы типа calcEigFuncWithFunc будут использовать новую сетку
    allocate(f1(1:2,maxLevel)) !сюда будем сохранять волновые функции
    allocate(f2(1:2,maxLevel)) !сюда будем сохранять волновые функции

    !minQuant можно было бы заменить на любое другое число, соответствующее
    !вычисленной функции, а вот offset нужно оставить, тк оно потом
    !учитывается в сопоставлении коэффициента с базисной функцией.
    !Сдвиг на N, чтобы получить верхний спектр
    f1=calcEigFuncExp(Ham, 1, offset)
    !Нормируем результат
    f1(:,:)=f1(:,:)/sqrt(integrateOnGrid(f1(1,:), f1(1,:)) + integrateOnGrid(f1(2,:), f1(2,:)))
    f2=calcEigFuncExp(Ham, N+1, offset)
    f2(:,:)=f2(:,:)/sqrt(integrateOnGrid(f2(1,:), f2(1,:)) + integrateOnGrid(f2(2,:), f2(2,:)))
    summa = integrateOnGrid(f1(1,:), f2(1,:)) + integrateOnGrid(f1(2,:), f2(2,:))
    print *, summa

    open(unit=19,file="wavefunc.dat",status="replace")
    do i=1,maxLevel
        write(19, "(F10.4)" , advance="no") testV(i)
        write(19, "(F15.5)" , advance="no") f1(1,i)
        write(19, "(F15.5)") f1(2,i)
    end do
    close(19)
    deallocate(f1)
    deallocate(f2)
    deallocate(testV)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Ниже идет много комметариев: это код, который я использовал
!для генерации теховских табличек для вставки в диплом.
!Он бесполезен, но удалять рука не подниамется :)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    open(unit=19, file='')
!    maxLevel = 5
!    write (*, *) "\bt{h}"
!    write (*, "(AI2AI2AI2AI2.2AI2AI2A)") "\caption{\captionTable{",floor(kappa),&
!    "}{",floor(Z),"}{",N,"}}\label{tab:Z",floor(Z),"ka",floor(kappa),"N",N,"}"
!    write(*, *) "\begin{tabulary}{\tableWidthMine}{| C | C | C | C |} \hline"
!    !Вставляем данные
!    write (*, *) " Состояние & Расчет & Точная & Разность \\ \hline"
!    do quantN = minQuant, maxLevel
!        i = quantN - minQuant + N + 1
!        write (*, "(A)", advance='no') trim(prettyStateName(quantL, quantN, quantJ))
!        write (*, "(A)", advance='no') " & "
!        write (*, 100, advance='no') E(i)
!        write (*, "(A)", advance='no') " & "
!        write (*, 100, advance='no') rightEnergy(quantN, quantJ)
!        write (*, "(A)", advance='no') " & "
!        write (*, "(E8.1)", advance='no') E(i) - rightEnergy(quantN, quantJ)
!        write (*, "(A)") " \\ \hline"
!        !write (*, "(A4I1A1A2I1A4)", advance='no') " & $", quantN,orbitalLetter,"_{",floor(quantJ*2),"/2}$"
!    end do
!    write(* ,*) "\end{tabulary}"
!    write(* ,*) "\end{table}"


!    !пишем заголовок
!    do quantN = minQuant, maxLevel
!        i = quantN - minQuant + N + 1
!        write (*, "(A4I1A1A2I1A4)", advance='no') " & $", quantN,orbitalLetter,"_{",floor(quantJ*2),"/2}$"
!    end do
!    write (*, *) "\\ \hline"
!    !пишем энергии из программы
!    write (*, "(A)", advance='no') "\lspы"
!    do quantN = minQuant, maxLevel
!        i = quantN - minQuant + N + 1
!        write (*, "(A)", advance='no') " & "
!        write (*, 100, advance='no')  E(i)
!    end do
!    write (*, *) "\\ \hline"
!    !пишем точные энергии
!    write (*, "(A)", advance='no') "Точная энергия"
!    do quantN = minQuant, maxLevel
!        i = quantN - minQuant + N + 1
!        write (*, "(A)", advance='no') " & "
!        write (*, 100, advance='no') rightEnergy(quantN, quantJ)
!    end do
!    write (*, *) " \\ \hline"

!    quantN = minQuant + maxLevel
    !print *, quantJ, rightEnergy(1, 0.5_WP)
!    close(19)
!    do i=N + 1 + maxLevel, N+1, -1
!!        print *, "E=",E(i)," n=",Z/sqrt(-2*E(i))
!
!!        write(19, *)
!        write (*, 100) rightEnergy(quantN, quantJ)
!        !write (*, "(E8.1)")  E(i) - rightEnergy(quantN, quantJ)
!!        print *, "E=",E(i), "right Energy=", rightEnergy(quantN, quantJ), &
!!        ' relative diff=',(E(i) - rightEnergy(quantN, quantJ))/abs(E(i))
!        quantN = quantN - 1
!        !0.5_WP*(1/quantN**2-1/c**2/quantN**3( 1/(0.5_WP + quantJ) - 0.75_WP/quantN ))
!    end do
    close(19)

    deallocate(Ham)
    deallocate(S)
    deallocate(E)

!    print *, "end"
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
contains
    function prettyStateName(iL, iN, iJ)
        character(100)::prettyStateName
        character::orbitalLetter
        integer:: iL, iN
        real(WP):: iJ

        SELECT CASE (iL)
        CASE (0)
            orbitalLetter = 's'
        CASE (1)
            orbitalLetter = 'p'
        CASE (2)
            orbitalLetter = 'd'
        CASE (3)
            orbitalLetter = 'f'
        CASE DEFAULT
            orbitalLetter = 'U'
        END SELECT
        write (prettyStateName, "(A2I1A1A2I1A4)") "$ ", iN,orbitalLetter,"_{",floor(iJ*2),"/2}$"
    end function prettyStateName
end program
