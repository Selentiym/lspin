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

    function lspinor(x) result(out)
        implicit none
        real(WP)::x, out

        out = (x**gammaRel) * exp(-x/2) * &
        ( &
            coeff1 * laguerrePoly(x, 2.0*gammaRel, Nr-1) &
            + coeff2 * laguerrePoly(x, 2.0*gammaRel, Nr) &
        )

    end function lspinor

    function laguerrePoly(x, upperIndex, lowerIndex) result(out)
        implicit none
        real(WP)::x,upperIndex, out
        integer::lowerIndex
        out = 1
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

    subroutine setParameters(iNr, iKappa, iZ, iLetter)
        !Задаем параметры lspinor'а
        real(WP)::iKappa, iZ !for variable explanation look module comments
        integer::iNr
        character::iLetter !for variable explanation look module comments
        !Сохраняем переданные параметры
        Nr = iNr
        kappa = iKappa
        Z = iZ
        letter = iLetter

        !На основе полученных параметров вычисляем коэффициенты при lspinor'ах
        gammaRel = sqrt(kappa**2 - (Z/c)**2)
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
