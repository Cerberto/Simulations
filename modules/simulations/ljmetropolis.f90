!-------------------------------------------------------------------------------
!
!   Module containing the definitions needed for the metropolis algorithm
!   applied to the Lennard-Jones fluid.
!
!-------------------------------------------------------------------------------

module ljmetr

    use kinds,      only: dp
    use ljmod
    implicit none

    public
    
    integer :: nth      ! thermalization sweeps
    integer :: nsw      ! effective sweeps
    integer :: skip     ! skipped sweeps (for uncorrelated data)
    integer :: ndat     ! effective size of the final sample
    integer :: nbin     ! number of sweeps per bin (for uncorrelated data)
    
    ! interface for the autocorrelation function
    interface
        function ac (x,t)
            use kinds, only: dp
            real(dp) :: ac
            integer :: i, D, t
            real(dp), dimension(:), allocatable :: x
        end function ac
    end interface

contains

    function autocorrelation (x, t)
        real(dp) :: autocorrelation
        integer :: t
        real(dp), dimension(:), pointer :: x

        integer :: i
        real(dp), dimension(3) :: temp = (/ 0, 0, 0 /)
        
        do i=1, ndat-t, 1
            temp(3) = temp(3) + x(i)*x(i+t)/(ndat - t)
            temp(1) = temp(1) + x(i)/(ndat - t)
            temp(2) = temp(2) + x(i)*x(i)/(ndat - t)
        end do
        autocorrelation = (temp(3) - temp(1)**2)/temp(2)
        
    end function autocorrelation
    

!
!   Routine implementing the Metropolis algorithm for a thermal distribution
!   in the canonical ensemble
!
    function thmetropolis_v (ptcls)
    
        real(dp) :: thmetropolis_v
        type(particle), dimension(:) :: ptcls
        integer :: i,j,k
        real(dp), dimension(:), allocatable :: u
        real(dp) :: t1, t2
        
        allocate(u(5))
        thmetropolis_v = 0
    
        do j=1, N
            call ranlxdf(u,5)
            k = int(N*u(5))
            k = modulo(k, N) + 1

            do i=1, 3, 1
                pstn_new(i) = ptcls(k)%pstn(i) + delta*(2*u(i)-1)
                t1 = pstn_new(i)/side
                call rintf(t1, t2)
                pstn_new(i) = pstn_new(i) - side*t2
            end do
    
            t1 = delta_interaction(ptcls,k)
            t2 = exp(-beta*t1)
            
            if(t2 >= u(4)) then
                ptcls(k)%pstn = pstn_new
                poten = poten + t1
                thmetropolis_v = thmetropolis_v + 1/N
            end if
        end do
        deallocate(u)
    end function thmetropolis_v


!
!   Routine implementing the Metropolis algorithm for a thermal distribution
!   in the isothermal-isobaric ensemble
!
    function thmetropolis_p (ptcls)
    
        real(dp) :: thmetropolis_p
        type(particle), dimension(:) :: ptcls
        integer :: i,j
        real(dp), dimension(:), allocatable :: u
        real(dp) :: t1, t2, side_t
        
        allocate(ptcls_new(N))
        thmetropolis_p = 0
        
        !
        !   Move randomly all the particles
        !
        allocate(u(3))
        do j=1, N
            call ranlxdf(u,3)

            do i=1, 3
                ptcls_new(j)%pstn(i) = ptcls(j)%pstn(i) + delta*(2*u(i)-1)
                t1 = ptcls_new(j)%pstn(i)/side
                call rintf(t1, t2)
                ptcls_new(j)%pstn(i) = ptcls_new(j)%pstn(i) - side*t2
            end do
        end do
        
        !
        !   Change the volume of the box and the particle positions accordingly
        !
        call ranlxdf(u,3)
        side_t = side*u(1)
        do j=1, N
            do i=1, 3
                ptcls_new(j)%pstn(i) = u(4)*ptcls_new(j)%pstn(i)
            end do
        end do
        
        t1 = total_interaction(ptcls_new)
        t2 = beta*(press*(side**3 - side_t**3) + poten - t1) + 3*N*log(u(1))
        t2 = exp(t2)
        
        if(t2 >= u(2)) then
            side = side_t
            do i=1, N
                ptcls(i)%pstn = ptcls_new(i)%pstn 
            end do
            poten = t1
            thmetropolis_p = thmetropolis_p + 1
        end if
        
        deallocate(ptcls_new)
        deallocate(u)
    end function thmetropolis_p
    
end module ljmetr
