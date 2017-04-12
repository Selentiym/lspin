module lspinors
    use constants
    use laguerre
    implicit none
    !Физический смысл переменных ниже условен, он имеет место
    !при релятивистском расчете одновалентных атомов и ионов
    real(WP)::Z                 !заряд ядра
    real(WP)::kappa             !собственное число (sigma,p)
    integer::Nr                !Главное радиальное квантовое число
    real(WP)::Norm              !Нормировка lspinor'a
    real(WP)::gammaRel          !Релятивистское квантовое число
    real(WP)::ULdiff            !Слагаемое, которое идет с разным знаком в зависимости от типа спинора
    real(WP)::coeff1,coeff2     !коэффициенты при полиномах лаггера в этом спиноре
    character::letter           !Буква, определяющая компоненту спинора
                                ! 'U' - верхняя, 'L' - нижняя
    real(WP), allocatable::lValues(:,:) !Хранит вычисленные однажды значения функции Лагерра
    integer::nDots, maxN    !Количество этих точек и максимальный порядок функции Лагерра
    real(WP), allocatable::dots(:)      !Точки, в которых нужны спиноры

    private::Z, kappa, Nr, Norm, gammaRel, ULdiff
!    private:: coeff1, coeff2
    private::nDots, maxN, dots, letter, lValues

contains
!    function lspinor(x) result(out)
!        implicit none
!        real(WP)::x, gammaRel !Релятивистское квантовое число
!        real(WP)::out
!        real(WP)::Norm !Нормировка lspinor'a
!        real(WP)::ULdiff !Слагаемое, которое идет с разным знаком в зависимости от типа спинора
!
!        gammaRel = sqrt(kappa**2 - (Z/c)**2)
!        Norm = sqrt(Nr**2 + 2*Nr*gammaRel + kappa**2)
!
!        ULdiff = ((Norm - kappa)/(Nr + 2.0*gammaRel)) * laguerrePoly(x, 2.0*gammaRel, Nr)
!
!        if (letter = 'L') then
!            ULdiff = -ULdiff
!        end if
!
!        out = Norm * (x**gammaRel) * exp(-x/2) * &
!        ( &
!            (-1.0_WP+deltaCronecker(Nr, 0.0))*laguerrePoly(x, 2.0*gammaRel, Nr-1.0) &
!            + ULdiff &
!        )
!
!        out = 2
!    end function lspinor

!    function testExt(ext,dot)
!        external ext
!        real(WP)::dot, ext, dot1, testExt
!        testExt = ext(dot)
!    end function testExt

    function oneToX(x)
        real(WP)::x, oneToX
        oneToX = 1/x
    end function oneToX

    function VPot(r)
        real(WP)::VPot,r
        !VPot = 0.0D+00
        VPot = -Z/r
    end function VPot

    function unityFunc(x)
        real(WP)::unityFunc,x
        unityFunc = 1.0_WP
    end function unityFunc

    function lspinorVector(iNr, iLetter) result(out)
        !На выходе столбец заданной в setGlobalParameters
        !длины nDots значений в точках dots
        !При фиксированном в setGlobalParameters значении
        !gammaRel lspinor задается iNr и iLetter
        integer::iNr
        character::iLetter
        real(WP)::out(nDots)
        !out = 0.0D+00
        out = lspinorVectorWithFunc(unityFunc, iNr, iLetter)
    end function lspinorVector

    function lspinorDerivativeVector(iNr, iLetter) result(out)
        !На выходе столбец заданной в setGlobalParameters
        !длины nDots значений в точках dots
        !При фиксированном в setGlobalParameters значении
        !gammaRel lspinor задается iNr и iLetter
        integer::iNr
        character::iLetter
        real(WP)::out(nDots)
        !out = 0.0D+00
        out = lspinorDerivativeVectorWithFunc(unityFunc, iNr, iLetter)
    end function lspinorDerivativeVector

    function lspinorVectorWithFunc(func, iNr, iLetter) result(out)
        !На выходе столбец заданной в setGlobalParameters
        !длины nDots значений в точках dots
        !При фиксированном в setGlobalParameters значении
        !gammaRel lspinor задается iNr и iLetter
        external func  !Для большей гибкости передаю функцию, на которую домножается Lspinor
        integer::iNr
        character::iLetter
        real(WP)::out(nDots), dot, func
        integer::k !Счетчик

        !Настраиваем нужный спинор
        call setLspinorParameters(iNr, iLetter)
!        print *, dots
        !Непосредственно считаем
        do k=1,nDots
            out(k) = lspinor(k) * func(dots(k))
        end do
    end function lspinorVectorWithFunc

    function lspinorDerivativeVectorWithFunc(func, iNr, iLetter) result(out)
        !На выходе столбец заданной в setGlobalParameters
        !длины nDots значений в точках dots
        !При фиксированном в setGlobalParameters значении
        !gammaRel lspinor задается iNr и iLetter
        external func  !Для большей гибкости передаю функцию, на которую домножается производная Lspinor
        integer::iNr
        character::iLetter
        real(WP)::out(nDots), dot, func
        integer::k !Счетчик

        !Настраиваем нужный спинор
        call setLspinorParameters(iNr, iLetter)
!        print *, dots
        !Непосредственно считаем
        do k=1,nDots
            out(k) = lspinorDerivative(k) * func(dots(k))
        end do
    end function lspinorDerivativeVectorWithFunc

    function lspinorDerivative(dotNum) result(out)
        integer::dotNum
        real(WP)::out
        out = gammaRel * lspinor(dotNum) / dots(dotNum) - &
        0.5_WP * lspinor(dotNum) + &
        (-coeff1 * laguerrePolyDerivative(dotNum, Nr - 1) + &
         coeff2 * laguerrePolyDerivative(dotNum, Nr))
    end function lspinorDerivative

    function lspinor(dotNum) result(out)
        implicit none
        real(WP)::out, test
        integer::dotNum
        !print *, coeff1
        !print *, coeff2
        test = gammaRel
        if (Nr .eq. -100) then
            out = 1.0D+00
        else
        !out = (dots(dotNum)**gammaRel) * &
        out =  & !x ** gammaRel тоже учитывается в интегралах
        !exp(-x/2) * & !Временно (или насовсем) убираем экпоненту, тк она при интегрировании учитывается
        ( &
            - coeff1 * laguerrePoly(dotNum, Nr-1) &
            + coeff2 * laguerrePoly(dotNum, Nr) &
        )

        !uncomment
        !out = laguerrePoly(dotNum, Nr)
        end if



    end function lspinor

    function laguerrePolyDerivative(dotNum, degree) result (out)
        implicit none
        real(WP)::out
        integer::dotNum, degree
        out = 1/dots(dotNum) * &
         ( &
            real(degree, WP) * laguerrePoly(dotNum, degree) &
             - (real(degree,WP) + 2.0_WP*gammaRel) * laguerrePoly(dotNum, degree - 1) &
         )
    end function laguerrePolyDerivative

    function laguerrePoly(dotNum, degree) result(out)
        implicit none
        real(WP)::out
        integer::dotNum, degree
        !lf_function
        if (degree < 0) then
            out = 0.0_WP
        else
            out = lValues(dotNum,degree)
        end if
    end function laguerrePoly

    function deltaCronecker(a,b) result (out)
        real(WP)::out
        integer::a,b
        if (a==b) then
            out = 1.0_WP
        else
            out = 0.0_WP
        end if
    end function deltaCronecker

    subroutine setLspinorGlobalParameters(iMaxN, iKappa, iZ, iNDots, iDots)
        real(WP),allocatable::debugLValues(:,:)
        !Они выделены, тк для них нужно запускать расчет функции Лагерра
        !Кроме того, для фиксированного kappa все функции Лагерра
        !в lspinor'ах одинаковы
        real(WP)::iKappa, iZ !for variable explanation look module comments
        !Максимальная степень полинома Лагерра, который понадобится
        ! и число точек, в которых нужны будут значения
        integer::iMaxN, iNDots
        !Сами точки, в которых считать
        real(WP)::iDots(iNDots)
        real(WP)::testGamma
        kappa = iKappa
        Z = iZ
        nDots = iNDots
        dots = iDots
        maxN = iMaxN
        maxN = maxN + 1
        gammaRel = calculateGamma(kappa, Z)
        testGamma = gammaRel
        !Переменная модуля, хранит значения функций Лагерра нужных порядков
        if (allocated(lValues)) deallocate (lValues)
        allocate(lValues(nDots,0:maxN))
        allocate(debugLValues(nDots,0:maxN))
        !Сохранили значения.

        call lf_function ( nDots, maxN, 2.0_WP * gammaRel, dots, debugLValues )
        lValues(:,:) = debugLValues(:,:)
        deallocate(debugLValues)
    end subroutine setLspinorGlobalParameters

    function calculateGamma(iKappa, iZ) result(out)
        real(WP)::out, iKappa, iZ
        out = sqrt(iKappa**2 - (iZ/c)**2)
    end function calculateGamma

    subroutine setLspinorParameters(iNr, iLetter)
        real(WP)::testGamma
        !Задаем параметры lspinor'а

        integer::iNr
        character::iLetter !for variable explanation look module comments
        !Сохраняем переданные параметры
        Nr = iNr
        letter = iLetter

        !На основе полученных параметров вычисляем коэффициенты при полиномах лагерра
        Norm = sqrt(Nr**2 + 2*Nr*gammaRel + kappa**2)

        testGamma = Norm

        ULdiff = ((Norm - kappa)/(Nr + 2.0*gammaRel))

        if (letter == 'L') then
            ULdiff = -ULdiff
        end if

        coeff1 = (-1.0_WP+deltaCronecker(Nr, 0)) * Norm
        testGamma = coeff1

        coeff2 = Norm * ULdiff
        testGamma = coeff2

    end subroutine setLspinorParameters

    function f(x)
        real(WP)::f, x
        f = x**4
    end function f

    function ampl(x)
        real(WP)::x, ampl
        ampl = exp(-x/2) * x ** (gammaRel)
    end function ampl

end module lspinors
