program testrandom
    implicit none
    
    double precision, dimension(:), allocatable :: x
    integer :: i,n
    read *, n
    
    allocate(x(0:n-1))
    call randf(x,n)
    
    do i=0, n-1, 1
        print *, x(i)
    end do
    
end program testrandom
