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

    private::kappa, gammaRel, Z, Nr, ULdiff, coeff1, coeff2, letter, lValues
    private::nDots, maxN, dots

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

    function lspinor(dotNum) result(out)
        implicit none
        real(WP)::out
        integer::dotNum
        out = (dots(dotNum)**gammaRel) * &
        !exp(-x/2) * & !Временно (или насовсем) убираем экпоненту, тк она при интегрировании учитывается
        ( &
            coeff1 * laguerrePoly(dotNum, Nr-1) &
            + coeff2 * laguerrePoly(dotNum, Nr) &
        )
    end function lspinor

    function laguerrePoly(dotNum, degree) result(out)
        implicit none
        real(WP)::out
        integer::dotNum, degree
        !lf_function
        out = lValues(dotNum,degree)
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

    subroutine setGlobalParameters(iMaxN, iKappa, iZ, iNDots, iDots)
        real(WP),allocatable::debugLValues(:,:)
        !Они выделены, тк для них нужно запускать расчет функции Лагерра
        real(WP)::iKappa, iZ !for variable explanation look module comments
        !Максимальная степень полинома Лагерра, который понадобится
        ! и число точек, в которых нужны будут значения
        integer::iMaxN, iNDots
        !Сами точки, в которых считать
        real(WP)::iDots(iNDots)
        kappa = iKappa
        Z = iZ
        nDots = iNDots
        dots = iDots
        maxN = iMaxN

        gammaRel = sqrt(kappa**2 - (Z/c)**2)
        !Переменная модуля, хранит значения функций Лагерра нужных порядков
        allocate(lValues(nDots,0:maxN))
        allocate(debugLValues(nDots,0:maxN))
        !Сохранили значения.

        call lf_function ( nDots, maxN, 2.0_WP * gammaRel, dots, debugLValues )
        lValues(:,:) = debugLValues(:,:)
        deallocate(debugLValues)
    end subroutine setGlobalParameters

    subroutine setParameters(iNr, iLetter)
        !Задаем параметры lspinor'а

        integer::iNr
        character::iLetter !for variable explanation look module comments
        !Сохраняем переданные параметры
        Nr = iNr
        letter = iLetter

        !На основе полученных параметров вычисляем коэффициенты при lspinor'ах

        Norm = sqrt(Nr**2 + 2*Nr*gammaRel + kappa**2)

        ULdiff = ((Norm - kappa)/(Nr + 2.0*gammaRel))

        if (letter == 'L') then
            ULdiff = -ULdiff
        end if

        coeff1 = (-1.0_WP+deltaCronecker(Nr, 0)) * Norm
        coeff2 = Norm * ULdiff


    end subroutine setParameters

    function f(x)
        real(WP)::f, x
        f = x**4
    end function f

end module lspinors
