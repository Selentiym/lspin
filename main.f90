!program compute
!    use constants
!    use matrix
!    use lspinors
!    use laguerre
!!    integer, parameter::N=5
!!    complex::A(N,N),lworkINIT(1),lrworkINIT(1),liworkINIT(1)
!!    real::OUTp(N),rwork(3*N-2)
!!    complex, allocatable :: work(:)
!!    integer::i,j,lwork,liwork,lrwork
!    integer :: N, i
!    real(WP), allocatable::x(:),w(:)
!    real(WP)::summa
!    N = 5
!    allocate (x(N), w(N))
!
!    call l_quadrature_rule(N, x, w)
!
!    summa=0.0_WP
!
!    do i=1,N
!        summa = summa + w(i)*f(x(i))
!    end do
!
!    deallocate(x,w)
!    print *, summa
!
!end program
