module matrix
    contains
    subroutine eigvect_c8 (A, OUTp)
        integer::N
        complex(8)::A(:,:)
        real(8), allocatable::OUTp(:),rwork(:)
        complex(8), allocatable :: work(:)
        integer::i,j,lwork

        N=size(A,1)

        allocate(rwork(3*N-2))

!        do i=1,N
!            print *, A(i,:)
!        end do

        lwork = (ilaenv(1,'ZHETRD','VU',N,N,N,N) + 1) * N
        allocate(work(lwork))
        call zheev('V','U',N,A,N,OUTp,work,lwork,rwork,info)

!        do i=1,N
!            print *, A(i,:)
!        end do
    end subroutine eigvect_c8

end module matrix
