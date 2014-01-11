module mod_vettori
    implicit none
    
contains
    
    subroutine vec_sum (x)
        integer, dimension(:) :: x
        integer :: i, j, N
        N = size(x)
        
        do i=1, N, 1
            x(i) = 0
            do j=1, i, 1
                x(i) = x(i) + j
            end do
        end do
    end subroutine
    
end module mod_vettori

program vettori
    use mod_vettori
    implicit none
    
    integer :: i, N
    integer, dimension(:), allocatable :: vec
    
    read*, N
    allocate(vec(N))
    do i=1, N, 1
        vec(i) = i
    end do
    
    call vec_sum(vec)
    
    do i=1, N, 1
        print "(i3)", "sum_{i=0}^", i, "i = ", vec(i)
    end do
    
    deallocate(vec)
end program vettori
