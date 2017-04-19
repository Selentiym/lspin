module matrix
    character::fileNameInMod
    contains
    subroutine eigvect_c8 (A, OUTp)
        integer::N
        complex(8)::A(:,:)
        real(8), allocatable::OUTp(:),rwork(:)
        complex(8), allocatable :: work(:)
        integer::i,j,lwork
        integer::toFile

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

    subroutine print_r8(A, firstIndex)
        real(8), allocatable::A(:,:)
        integer::i,k, firstIndex, upperBound1, upperBound2
        upperBound1 = size(A,1) - 1 + firstIndex
        upperBound2 = size(A,2) - 1 + firstIndex
        if (toFile == 1) then
            open(unit=17, file=fileNameInMod,Access = 'append',status="old")
            do i=firstIndex,upperBound1
                do j=firstIndex,upperBound2
                    write(17, "(F10.2)" , advance="no") A(i,j)
                end do
                write(17,*) ""
            end do
            write(17, *) ""
            close(17)
        else
            do i=firstIndex,upperBound1
                do j=firstIndex,upperBound2
                    write(*, "(F10.4)" , advance="no") A(i,j)
                end do
                write(*,*) ""
            end do
        end if
    end subroutine print_r8

    function from4to2dim(f,s,A) result(out)
        real(8), allocatable::A(:,:,:,:), out(:,:)
        integer::f,s,fSize,sSize
        fSize = size(A,3)
        sSize = size(A,4)
        allocate (out(fSize,sSize))
        out(:,:) = A(f,s,:,:)
    end function from4to2dim

    subroutine print4dimAs2(f,s,A)
        real(8), allocatable::A(:,:,:,:)
        integer::f,s
        call print_r8(from4to2dim(f, s, A), 1)
    end subroutine print4dimAs2

    subroutine setWriteToFile(iFileName)
        character::iFileName
        fileNameInMod = iFileName
        open(unit=17,file=fileNameInMod,status="replace")
        close(17)
        toFile = 1
    end subroutine setWriteToFile

end module matrix
