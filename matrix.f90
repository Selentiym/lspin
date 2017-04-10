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

    subroutine general_eigvect_r8 (H, S, OUTp)
    real(8),allocatable::OUTp(:), work(:), iwork(:), H(:,:),S(:,:)
    integer:: N, itype, lda, ldb, lwork, liwork, info
    character:: jobz, uplo

    N=size(H,1)

    itype = 1 !Ac=\lambda Bc

    jobz = 'V' !Both: values and vectors are needed

    uplo = 'U' !Upper triangles are stored

    lda = N !all matrix has to be diagonalized
    ldb = N

    lwork = 1 + 6*N + 2*N**2
    allocate(work(lwork))

    liwork = 3 + 5 * N
    allocate(iwork(liwork))

    call dsygvd(itype, jobz, uplo, N, H, lda, S, ldb, OUTp, work, lwork, iwork, liwork, info)

    deallocate(iwork)
    deallocate(work)

    end subroutine general_eigvect_r8

end module matrix
