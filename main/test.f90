program test

    use kinds, only: dp
    implicit none
    
    real(dp), dimension(:), allocatable :: x
    integer :: i,n,j
    
    read *, n
    allocate(x(n))
    
    call rlxdinit(1,rand(time()))
    
    do j=1, 2, 1
        call ranlxdf(x,n)
        do i=1, n, 1
            print *, x(i)
        end do
        print *, " "
    end do
    
end program test
