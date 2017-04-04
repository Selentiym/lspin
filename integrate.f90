module integrate
    use laguerre
    use constants
    implicit none

    integer::qLen       !Длина квадратурной формулы, совпадает с порядком
                        !полинома Лагерра, корнями которого являются узлы
    real(WP), allocatable::x(:),w(:) !Узлы и веса квадратуры соответственно

    private qLen, x, w

contains
    subroutine setQuadratureParameters(iQLen)
        integer::iQLen
        !Сохраянем на будущее число точек
        qLen = iQLen
        !Выделяем память и вычисляем узлы и веса
        allocate(x(qLen),w(qLen))
        call l_quadrature_rule(qLen, x, w)
    end subroutine setQuadratureParameters

    function integrateOnGrid(f1,f2) result(out)
        !На входе должно быть два вектора длины qLen со
        !значениями функции в узлах квадратуры
        real(WP)::f1(qLen), f2(qLen)
        real(WP)::out
        integer::i
        do i=1,qLen
            out = out + f1(i)*f2(i)*w(i)
        end do
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

end module integrate
