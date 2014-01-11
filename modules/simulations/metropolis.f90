

module metropolis

    use jackknife, only: JK, JK_cluster, JK_init, JK_function
    implicit none

    private none
    public
    
    interface
        function ac (x,t)
            real :: ac
            integer :: i, D, t
            real, dimension(:), allocatable :: x
        end function ac
        
        function weight (x)
            real :: weight
            real, dimension(:), allocatable :: x
        end function weight
    end interface
    
    abstract interface
        subroutine metr (f, x, delta)
            procedure(weight), external, pointer :: P
            real, dimension(:), allocatable :: state
            real :: delta, swap, x_new, acceptance
            integer :: i, D
            real, dimension(2) :: u
        end subroutine metr
    end interface

contains

    function autocorrelation (x, t)
        real :: autocorrelation
        integer :: t
        real, dimension(:), pointer :: x
        
        integer :: D = size(x)
        integer :: i
        real, dimension(3) :: temp = (/ 0, 0, 0 /)
        
        do i=0, i<D-t
            temp(3) = temp(3) + x(i)*x(i+t)/real(D - t)
            temp(1) = temp(1) + x(i)/real(D - t)
            temp(2) = temp(2) + x(i)*x(i)/real(D - t)
        end do
    end function autocorrelation

    ! Most general metropolis algorithm
    subroutine metr_gen (P,state,delta)
        procedure(weight), external, pointer :: P
        real, dimension(:), pointer :: state
        real :: delta
        
        ! ausiliary variables
        integer :: i, D
        real :: swap, x_new, acceptance
        real, dimension(2) :: u
        
        D = size(state)
        do i=0, D, 1
            ! ranlxd(u,2);
            ! AGGIUNGERE COMANDO PER INIZIALIZZARE RANDOM u
            x_new = state(i) + delta*(2*u(1) - 1.0)
            swap = state(i)
            state(i) = x_new
                acceptance = P(state)
            state(i) = swap
                acceptance = acceptance / P(state)
		
            if(acceptance >= u(2))
                state(i) = x_new
        end do
    end subroutine metr_gen

    
    ! Metropolis algorithm specialized to "thermal" states, i.e., with
    ! probability densities proportional to exp{F}=exp{beta*f}
    subroutine metr_th (F,state,delta)
        procedure(weight), external, pointer :: F
        real, dimension(:), pointer :: state
        real :: delta
        
        ! ausiliary variables
        integer :: i, D
        real :: swap, x_new, acceptance
        real, dimension(2) :: u
        
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
