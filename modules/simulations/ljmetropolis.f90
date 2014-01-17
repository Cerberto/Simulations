!-------------------------------------------------------------------------------
!
!   Module containing the definitions needed for the metropolis algorithm
!   applied to the Lennard-Jones fluid.
!
!   The 'thmetropolis' routine by default calculates also the kinetic energy via
!   the virial theorem. If you are not interested in this calculation, comment
!   the following lines
!       
!       real(dp), external :: delta_virial
!       kinen = kinen + delta_virial(ptcls,k)
!
!   in the 'thmetropolis' routine. To initialize the variable 'kinen' (defined
!   in th 'ljmod' module) one has to explicitly call the function 'ke_virial'
!
!-------------------------------------------------------------------------------

module ljmetr

    use kinds,      only: dp
    use ljmod
    use jackknife,  only: JK, JK_cluster, JK_init !, JK_function
    implicit none

    public
    
    integer :: tmax     ! maximum "time" for autocorrelation
    integer :: nth      ! thermalization sweeps
    integer :: nsw      ! effective sweeps
    integer :: skip     ! skipped sweeps (for uncorrelated data)
    integer :: ndat     ! effective size of the final sample
    
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
        
        do i=0, ndat-t, 1
            temp(3) = temp(3) + x(i)*x(i+t)/(ndat - t)
            temp(1) = temp(1) + x(i)/(ndat - t)
            temp(2) = temp(2) + x(i)*x(i)/(ndat - t)
        end do
    end function autocorrelation
    

!
!   Routine implementing the Metropolis algorithm for a thermal distribution.
!
    subroutine thmetropolis (ptcls)
        !use kinds, only: dp
        !use ljmod, only: particle, delta, pstn_new, poten, side
        !implicit none
    
        type(particle), dimension(:) :: ptcls
        integer :: i,k
        real(dp), dimension(:), allocatable :: u
        real(dp) :: t1, t2
        real(dp), external :: delta_interaction
        real(dp), external :: delta_virial
    
        allocate(u(5))
        call ranlxdf(u,5)
    
        ! select randomly the particle to move
            u(5) = N-(N-1)*u(5)
            k = int(t1)
        do i=1, 3, 1
            pstn_new(i) = ptcls(k)%pstn(i) + delta*(2*u(i)-1)
            t1 = pstn_new(i)/side
            call rintf(pstn_new(i), t2)
            pstn_new(i) = pstn_new(i) - side*t2
        end do
    
        t1 = delta_interaction(ptcls,k)
        t2 = exp(-t1)
            
        if(t2 >= u(4)) then
            kinen = kinen + delta_virial(ptcls,k)
            ptcls(k)%pstn = pstn_new
            poten = poten + t1
        end if
    
        deallocate(u)

    end subroutine thmetropolis

end module ljmetr
