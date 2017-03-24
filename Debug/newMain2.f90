!program compute
!    use constants
!    use matrix
!    complex(8), allocatable::A(:,:)
!    real(8), allocatable::o(:)
!    integer::N
!    N=2
!    allocate(A(N,N))
!    allocate(o(N))
!!    do i=1,N
!!        do j=1,N
!!            A(i,j)=0
!!        end do
!!    end do
!!    do i=1,N
!!        A(i,i) = 2
!!        j = i-1
!!        if (j > 0) then
!!            A(j,i)=-1
!!        end if
!!        j=i+1
!!        if (j <= N) then
!!            A(j,i)=-1
!!        end if
!!
!!!        B(i)=0
!!    end do
!
!    A(1,1) = 0
!    A(1,2) = (0,-1)
!    A(2,1) = (0,1)
!    A(2,2) = 0
!
!    call eigvect_c8(A, o)
!
!    print *, o
!    do i=1,N
!        print *, A(i,:)
!    end do
!
!end program
