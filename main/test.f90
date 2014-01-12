program test

    use kinds, only: dp
    implicit none
    
    
    real(dp), dimension(:), allocatable :: x
    integer :: i,n
    
    read *, n
    
    allocate(x(0:n-1))
    call ranlxdf(x,n)
    
    do i=0, n-1, 1
        print *, x(i)
    end do
    
end program test
