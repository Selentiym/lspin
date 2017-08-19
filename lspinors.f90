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
    real(WP)::realNorm          !нормировка L-spinor'ов
    character::letter           !Буква, определяющая компоненту спинора
                                ! 'U' - верхняя, 'L' - нижняя
    real(WP), allocatable::lValues(:,:) !Хранит вычисленные однажды значения функции Лагерра
    real(WP), allocatable::lDerivValues(:,:)    !Хранит вычисленные однажды значения
                                                !производной функции Лагерра
    integer::nDots, maxN    !Количество этих точек и максимальный порядок функции Лагерра
    real(WP), allocatable::dots(:)      !Точки, в которых нужны спиноры

    real(WP)::coordScale = 1.0_WP   !Поскольку может понадобиться растянуть базисные функции,
                                    !а это делать неудобно из-за интегрирования, растягивать будем в другую сторону все
                                    ! все остальные посредством coordsScale
    logical ::addExponent = .false. !Добавлять ли в спинор экспоненту и степень
    real(WP)::powerCached           !Хранит показатель степени, которая будет выдана
                                    !функцией powerPreset
    private::Z, kappa, Nr, Norm, gammaRel, ULdiff, realNorm, coordScale
!    private:: coeff1, coeff2
    private::nDots, maxN, dots, letter, lValues, lDerivValues

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
        real(WP)::out(nDots), dot, func, arg
        integer::k !Счетчик

        !Настраиваем нужный спинор
        call setLspinorParameters(iNr, iLetter)
!        print *, dots
        !Непосредственно считаем
        if (addExponent) then
            do k=1,nDots
                arg = coordScale * dots(k)
                out(k) = lspinor(k) * func(arg) * exp(-arg*0.5_WP) * arg**gammaRel
!                print *, arg, out(k),k
            end do
        else
            do k=1,nDots
                out(k) = lspinor(k) * func(coordScale * dots(k))
            end do
        end if

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
            out(k) = lspinorDerivative(k) * func(coordScale * dots(k))
        end do
    end function lspinorDerivativeVectorWithFunc

    function lspinorDerivative(dotNum) result(out)
        integer::dotNum
        real(WP)::out

        out = gammaRel * lspinor(dotNum) / dots(dotNum) - &
        0.5_WP * lspinor(dotNum) + 1.0_WP * &
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
        test = 0.0D+00
        test = coeff1
        test = 0.0D+00
        test = coeff2


        !out = (dots(dotNum)**gammaRel) * &
        out =  1.0_WP * & !x ** gammaRel тоже учитывается в интегралах
        !exp(-x/2) * & !Временно (или насовсем) убираем экпоненту, тк она при интегрировании учитывается
        ( &
            - coeff1 * laguerrePoly(dotNum, Nr-1) &
            + coeff2 * laguerrePoly(dotNum, Nr) &
        )

        !out = laguerrePoly(dotNum, Nr)
    end function lspinor

    function lagVector(iNr, iLetter) result(out)
        !На выходе столбец заданной в setGlobalParameters
        !длины nDots значений в точках dots
        !При фиксированном в setGlobalParameters значении
        !gammaRel lspinor задается iNr и iLetter
        integer::iNr
        character::iLetter
        real(WP)::out(nDots)
        !out = 0.0D+00
        out = lagVectorWithFunc(unityFunc, iNr, iLetter)
    end function lagVector

    function lagDerivativeVector(iNr, iLetter) result(out)
        !На выходе столбец заданной в setGlobalParameters
        !длины nDots значений в точках dots
        !При фиксированном в setGlobalParameters значении
        !gammaRel lspinor задается iNr и iLetter
        integer::iNr
        character::iLetter
        real(WP)::out(nDots)
        !out = 0.0D+00
        out = lagDerivativeVectorWithFunc(unityFunc, iNr, iLetter)
    end function lagDerivativeVector

    function lagVectorWithFunc(func, iNr, iLetter) result(out)
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
            out(k) = laguerrePoly(k, iNr) * func(coordScale * dots(k))
        end do
        out(:) = out(:)/sqrt(laguerrePolyNorm(iNr))
    end function lagVectorWithFunc

    function lagDerivativeVectorWithFunc(func, iNr, iLetter) result(out)
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
            !Базисная функция - не просто полином Лагерра, а домноженный
            !на exp(-x/2)*x**gammaRel
            out(k) = func(coordScale * dots(k)) * &
            (&
                laguerrePolyDerivative(k, iNr) &
                -0.5_WP * laguerrePoly(k, iNr)  &
                + gammaRel/dots(k) * laguerrePoly(k, iNr) &
            ) / sqrt(laguerrePolyNorm(iNr))
        end do
    end function lagDerivativeVectorWithFunc

    function laguerrePolyDerivative(dotNum, degree) result (out)
        implicit none
        real(WP)::out
        integer::dotNum, degree
        out = lDerivValues(dotNum, degree)
!        out = 1/dots(dotNum) * &
!         ( &
!            real(degree, WP) * laguerrePoly(dotNum, degree) &
!             - (real(degree,WP) + 2.0_WP*gammaRel) * laguerrePoly(dotNum, degree - 1) &
!         )
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
        real(WP),allocatable::debugLValues(:,:), debugLDerivValues(:,:)
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
        if (allocated(lDerivValues)) deallocate (lDerivValues)
        allocate(lValues(nDots,0:maxN))
        allocate(lDerivValues(nDots,0:maxN))
        allocate(debugLValues(nDots,0:maxN))
        allocate(debugLDerivValues(nDots,0:maxN))
        !Сохранили значения.

        call lf_function_derivative ( nDots, maxN, 2.0_WP * gammaRel, dots, debugLValues, debugLDerivValues )
        lValues(:,:) = debugLValues(:,:)
        lDerivValues(:,:) = debugLDerivValues(:,:)
        deallocate(debugLValues)
        deallocate(debugLDerivValues)
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

!        if ((Nr == 0).and.(kappa > 0)) then
!            Norm = -Norm
!        end if

        ULdiff = ((Norm - kappa)/(Nr + 2.0*gammaRel))

!        if (ULdiff == 0.0_WP) then
!            ULdiff = 1.0_WP
!        end if

        if (letter == 'L') then
            ULdiff = -ULdiff
        end if

        coeff1 = (1.0_WP-deltaCronecker(Nr, 0)) * Norm
        testGamma = coeff1

        coeff2 = Norm * ULdiff
        testGamma = coeff2



        realNorm = sqrt(coeff1 ** 2 *laguerrePolyNorm(Nr - 1) + coeff2 ** 2 * laguerrePolyNorm(Nr))

        if (realNorm == 0.0D+00) then
            coeff1 = 0.0D+00
            coeff2 = 0.0D+00
        else
            coeff1 = coeff1 / realNorm
            coeff2 = coeff2 / realNorm
        end if

    end subroutine setLspinorParameters

    function laguerrePolyNorm(iNr) result(out)
        integer::iNr
        real(WP)::out
        if (iNr >= 0) then
            out = gamma(2.0D+00*gammaRel + real(iNr + 1,WP))/gamma(real(iNr + 1,WP))
        else
            out = 0.0D+00
        end if
    end function laguerrePolyNorm

    function f(x)
        real(WP)::f, x
        f = x**4
    end function f

    function ampl(x)
        real(WP)::x, ampl
        ampl = exp(-x/2) * x ** (gammaRel)
    end function ampl

    function calcEigFuncExp(matr, num, qOffset) result(funcRez)
        real(WP), allocatable :: matr(:,:)
        integer :: num, qOffset
        real(WP) :: funcRez(2,nDots)
        funcRez = calcEigFuncWithFuncAndExp(matr, num, unityFunc, qOffset)
    end function calcEigFuncExp

    function calcEigFuncWithFuncAndExp(matr, num, func, io) result(tempMatr)
        external func
        real(WP), allocatable :: matr(:,:)
        real(WP) :: funcRez(2,nDots), lsp(nDots), func, temp1, tempMatr(2,nDots)
        character :: tempLetter
        integer :: num, s, j, ind, N, io, offs
        funcRez(:,:) = 0.0_WP
        offs = io - 1
        !qOffset = -qOffset + 1
        print *, offs
        s = size(matr,1)
        N = floor(s/2.0_WP)
        addExponent = .true.
        do j=1,N
!            temp1 =
            ind = j
            tempLetter = 'U'

            lsp = lspinorVectorWithFunc(func, j + offs, tempLetter)
            funcRez(1,:) = funcRez(1,:) + lsp(:) * matr(ind,num)

            ind = j + N
            tempLetter = 'L'

            lsp = lspinorVectorWithFunc(func, j + offs, tempLetter)
            funcRez(2,:) = funcRez(2,:) + lsp(:) * matr(ind,num)
        end do
        addExponent = .false.
        tempMatr(:,:) = funcRez(:,:)
    end function calcEigFuncWithFuncAndExp


    function calcEigFunc(matr, num, qOffset) result(funcRez)
        real(WP), allocatable :: matr(:,:)
        integer :: num, qOffset
        real(WP) :: funcRez(2,nDots)
        funcRez = calcEigFuncWithFunc(matr, num, unityFunc, qOffset)
    end function calcEigFunc

    function calcEigFuncWithFunc(matr, num, func, io) result(tempMatr)
        external func
        real(WP), allocatable :: matr(:,:)
        real(WP) :: funcRez(2,nDots), lsp(nDots), func, temp1, tempMatr(2,nDots)
        character :: tempLetter
        integer :: num, s, j, ind, N, io, offs
        funcRez(:,:) = 0.0_WP
        offs = io - 1
        !qOffset = -qOffset + 1
        print *, offs
        s = size(matr,1)
        N = floor(s/2.0_WP)
        do j=1,N
!            temp1 =
            ind = j
            tempLetter = 'U'

            lsp = lspinorVectorWithFunc(func, j + offs, tempLetter)
            funcRez(1,:) = funcRez(1,:) + lsp(:) * matr(ind,num)

            ind = j + N
            tempLetter = 'L'

            lsp = lspinorVectorWithFunc(func, j + offs, tempLetter)
            funcRez(2,:) = funcRez(2,:) + lsp(:) * matr(ind,num)
        end do
        tempMatr(:,:) = funcRez(:,:)
    end function calcEigFuncWithFunc

    function calcExtendedDeriv(func, deriv, fParam,dParam, point) result (out)

        !external func
        !external deriv
        real(WP)::out, func, deriv, test123
        integer::point, fParam, dParam


        out = -0.5_WP * func(point, fParam) &
        + gammaRel * func(point, fParam) / dots(point) + deriv(point, dParam)

    end function calcExtendedDeriv

    function powerPreset(point) result(out)
        real(WP)::out, point
        out = point ** powerCached
    end function powerPreset

    function power(point, powerC) result(out)
        real(WP)::out
        integer::point, powerC
        out = dots(point) ** real(powerC,WP)
    end function power

    function powerDeriv(point, power) result(out)
        real(WP)::out,powerDerivCoeff
        integer::point, power

        powerDerivCoeff = real(power,WP)
        if (powerDerivCoeff == 0.0_WP) then
            out = 0.0_WP
        else
            out = powerDerivCoeff * dots(point) ** (powerDerivCoeff - 1.0_WP)
        end if
    end function powerDeriv


    subroutine lf_function_derivative ( m, n, alpha, x, cx, cd )

!*****************************************************************************80
!
!! LF_FUNCTION evaluates the Laguerre function Lf(n,alpha,x).
!
!  Recursion:
!
!    Lf(0,ALPHA,X) = 1
!    Lf(1,ALPHA,X) = 1+ALPHA-X
!
!    Lf(N,ALPHA,X) = (2*N-1+ALPHA-X)/N * Lf(N-1,ALPHA,X)
!                      - (N-1+ALPHA)/N * Lf(N-2,ALPHA,X)
!
!  Restrictions:
!
!    -1 < ALPHA
!
!  Special values:
!
!    Lf(N,0,X) = L(N,X).
!    Lf(N,ALPHA,X) = LM(N,ALPHA,X) for ALPHA integral.
!
!  Norm:
!
!    Integral ( 0 <= X < +oo ) exp ( - X ) * Lf(N,ALPHA,X)^2 dX
!    = Gamma ( N + ALPHA + 1 ) / N!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order function to compute.
!
!    Input, real ( kind = 8 ) ALPHA, the parameter.  -1 < ALPHA is required.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(1:M,0:N), the functions of
!    degrees 0 through N evaluated at the points X.

!    Output, real ( kind = 8 ) CD(1:M,0:N), the derivatives of functions of
!    degrees 0 through N evaluated at the points X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) cx(1:m,0:n)
  real ( kind = 8 ) cd(1:m,0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(1:m)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LF_FUNCTION - Fatal error!'
    write ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
    write ( *, '(a)' ) '  but ALPHA must be greater than -1.'
    stop
  end if

  if ( n < 0 ) then
    return
  end if

  cx(1:m,0) = 1.0D+00
  cd(1:m,0) = 0.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1:m,1) = 1.0D+00 + alpha - x(1:m)
  cd(1:m,1) = -1.0D+00

  do i = 2, n
    cx(1:m,i) = ( &
      ( real ( 2 * i - 1, kind = 8 ) + alpha - x(1:m) ) * cx(1:m,i-1)   &
    + ( real (   - i + 1, kind = 8 ) - alpha          ) * cx(1:m,i-2) ) &
      / real (     i,     kind = 8 )
    cd(1:m,i) = ( &
      ( real ( 2 * i - 1, kind = 8 ) + alpha - x(1:m) ) * cd(1:m,i-1) - cx(1:m,i-1)   &
    + ( real (   - i + 1, kind = 8 ) - alpha          ) * cd(1:m,i-2) ) &
      / real (     i,     kind = 8 )
!    cd(1:m,i) = ( &
!      ( real ( 2 * i - 1, kind = 8 ) + alpha - x(1:m) ) * cd(1:m,i-1)   &
!    + ( real (   - i + 1, kind = 8 ) - alpha          ) * cd(1:m,i-2)   &
!    - cx(1:m,i-1)*0 ) &
!      / real (     i,     kind = 8 )
  end do

  return
end subroutine lf_function_derivative

function rightRadius(inpN, inpL, inpJ)
    integer::inpN, inpL
    real(WP)::rightRadius, inpJ, en, eps
    en = rightEnergy(inpN, inpJ)
    eps = 1 + en/c**2
    rightRadius = - 0.5_WP/Z * (1.5_WP * (Z**2/eps) *(c**2 + eps)/(c**2 + 0.5_WP * eps) &
    - kappa * (1.0_WP + kappa * eps) )
    !print *, inpN, inpL, inpJ
    rightRadius = 0.5_WP/Z * (3*inpN**2 - inpL*(inpL+1))
end function rightRadius

function rightEnergy(inpN,j)
    integer :: inpN
    real(WP)::n,j, rightEnergy, nr, num, root, summTemp
    n = real(inpN,8)
    nr = n - gammaRel
    num = sqrt(n**2 - 2*nr*(abs(kappa) - gammaRel))
    !rightEnergy = -Z**2/(num*(num + nr + gammaRel))
    !rightEnergy = -0.5_WP*Z**2*(1/n**2 + 1/(c**2*n**3) * ( 1/(0.5_WP + j) - 0.75_WP/n ))
    !rightEnergy = -Z**2/(num*(num+nr+gammaRel))
    nr=n-abs(kappa)
    summTemp = nr+gammaRel
    root = sqrt(1.0_WP+(Z/c)**2/summTemp**2)
    rightEnergy = -Z**2/summTemp**2/(root*(root+1.0_WP))
end function rightEnergy

subroutine setBasisScale(iScale)
    real(WP) :: iScale
    coordScale = 1.0D+00/iScale
end subroutine setBasisScale

end module lspinors
