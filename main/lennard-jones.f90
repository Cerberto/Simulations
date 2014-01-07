program LJ
    
    
    use jackknife, only: JK, JK_init, JK_cluster, JK_function
    implicit none
    
    type(particle), dimension(:), pointer :: ptcls
    

end program LJ

module LJ_global
    implicit none
    
    type :: particle
        real, dimension(3) :: x
    end type particle
    
    ! cube size
    real, parameter :: L = 100.0
    integer :: N
    
    interface
        subroutine init (p, s)
            real :: s
            type(particle), dimension(:), pointer :: p
        end subroutine init
    end interface

contains

    subroutine particle_init (ptcls, side)
        real :: side    ! cube side
        type(particle), dimension(:), pointer :: ptcls  ! set of particles
        real :: lspc    ! lattice spacing
        
    end subroutine particle_init
    
end module LJ_global
