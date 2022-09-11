module io_routines
implicit none
    
contains
    subroutine write_mat_int(mat)
        integer :: i,j
        integer, dimension(:,:), intent(in) :: mat
        do i=1,size(mat,dim=1)
            do j=1,size(mat,dim=2)
                write (*,'(100(i5,2x))',advance='no') mat(i,j)
            end do
            print *, ' '
        end do
    end subroutine write_mat_int

    subroutine write_vec_int(vec)
        integer :: i
        integer, dimension(:), intent(in) :: vec
        do i=1,size(vec)
            write (*,'(100(i5,2x))',advance='no') vec(i)
        end do
    end subroutine write_vec_int
    
    subroutine write_mat_real(mat)
        integer :: i,j
        real, dimension(:,:), intent(in) :: mat
        do i=1,size(mat,dim=1)
            do j=1,size(mat,dim=2)
                write (*,'(100(f9.3,2x))',advance='no') mat(i,j)
            end do
            print *, ' '
        end do
    end subroutine write_mat_real

    subroutine write_vec_real(vec)
        integer :: i
        real, dimension(:), intent(in) :: vec
        do i=1,size(vec)
            write (*,'(100(f10.3,2x))',advance='no') vec(i)
        end do
    end subroutine write_vec_real
    
end module io_routines