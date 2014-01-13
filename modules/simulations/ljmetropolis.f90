

module metropolis

    use kinds, only: dp
    use jackknife, only: JK, JK_cluster, JK_init !, JK_function
    implicit none

    private none
    public
    
    interface
        function ac (x,t)
            real(dp) :: ac
            integer :: i, D, t
            real(dp), dimension(:), allocatable :: x
        end function ac
        
        function weight (x)
            real(dp) :: weight
            real(dp), dimension(:), allocatable :: x
        end function weight
    end interface
    
    abstract interface
        subroutine metr (f, x, delta)
            procedure(weight), external, pointer :: P
            real(dp), dimension(:), allocatable :: state
            real(dp) :: delta, swap, x_new, acceptance
            integer :: i, D
            real(dp), dimension(2) :: u
        end subroutine metr
    end interface

contains

    function autocorrelation (x, t)
        real(dp) :: autocorrelation
        integer :: t
        real(dp), dimension(:), pointer :: x
        
        integer :: D = size(x)
        integer :: i
        real(dp), dimension(3) :: temp = (/ 0, 0, 0 /)
        
        do i=0, i<D-t
            temp(3) = temp(3) + x(i)*x(i+t)/real(dp)(D - t)
            temp(1) = temp(1) + x(i)/real(dp)(D - t)
            temp(2) = temp(2) + x(i)*x(i)/real(dp)(D - t)
        end do
    end function autocorrelation
    
    
    ! Metropolis algorithm specialized to "thermal" states, i.e., with
    ! probability densities proportional to exp{F}=exp{beta*f}
    subroutine metr_th (F,state,delta)
        procedure(weight), external, pointer :: F
        real(dp), dimension(:), pointer :: state
        real(dp) :: delta
        
        ! ausiliary variables
        integer :: i, D
        real(dp) :: swap, x_new, acceptance
        real(dp), dimension(2) :: u
        
        D = size(state)
        do i=0, D, 1
            ! ranlxd(u,2);
            ! AGGIUNGERE COMANDO PER INIZIALIZZARE RANDOM u
            x_new = state(i) + delta*(2*u(1) - 1.0)
            swap = state(i)
            state(i) = x_new
                acceptance = F(state)
            state(i) = swap
                acceptance = acceptance - F(state)
            acceptance = exp(acceptance)

            if(acceptance >= u(2))
                state(i) = x_new
        end do
    end subroutine metr_th
