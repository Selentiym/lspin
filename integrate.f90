module integrate
    use laguerre
    use constants
    implicit none

    integer::qLen       !Длина квадратурной формулы, совпадает с порядком
                        !полинома Лагерра, корнями которого являются узлы
    real(WP), allocatable::x(:),w(:) !Узлы и веса квадратуры соответственно
    real(WP) :: offset  !Сдвигквадратура относительно базисных функций. Для правильного
                        !подсчета интегралов типа 1/x**param, param <= offset, param целый
                        !Квадратура будет считаться для exp(-x)*x**(alpha - offset)
    real(WP) :: coordScale = 1.0_WP
    private qLen, x, w, offset, coordScale

contains
    subroutine setQuadratureParameters(iQLen, iAlpha, iOffset)
        !Генерирует и сохраняет точки и веса для дальнейшего
        !интегрирования. Считаются интегралы типа
        !exp(-x)*x**iAlpha*f(x)
        integer::iQLen
        real(WP)::iAlpha, iOffset, temp
        !Сохраянем на будущее число точек
        qLen = iQLen
        !Выделяем память и вычисляем узлы и веса
        if (allocated(x)) deallocate(x)
        if (allocated(w)) deallocate(w)
        allocate(x(qLen),w(qLen))
        temp = iAlpha - iOffset
        !print *, temp
        call lf_quadrature_rule(qLen, temp, x, w)
        offset = iOffset
    end subroutine setQuadratureParameters

    function integrateOnGrid(f1,f2) result(out)
        !На входе должно быть два вектора длины qLen со
        !значениями функции в узлах квадратуры
        real(WP)::f1(qLen), f2(qLen)
        real(WP)::out
        integer::i
        out = 0.0D+00
        do i=1,qLen
            out = out + f1(i)*f2(i)*w(i)*x(i)**(offset)
        end do
        out = out / coordScale
    end function integrateOnGrid

    function getIntegrateGrid()
        !Возвращает точки, в которых необходимо знать значение
        !функции для интегрирования
        real(WP)::getIntegrateGrid(qLen)
        getIntegrateGrid = x
    end function getIntegrateGrid

    function getIntegrateGridLength() result(out)
        integer::out
        out = qLen
    end function getIntegrateGridLength

    function getIntegrateWeights() result (wOut)
        real(WP),allocatable::wOut(:)
        allocate(wOut(qLen))
        wOut(:) = w(:)
    end function getIntegrateWeights

    subroutine setQuadratureArgumentScale(iScale)
        real(WP)::iScale
        coordScale = iScale
    end subroutine setQuadratureArgumentScale

end module integrate
