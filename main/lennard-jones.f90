module LJ_global
    implicit none
    
    type :: particle
        real, dimension(3) :: x
    end type particle
    
    ! cube size
    real :: L
    integer :: N

contains

    subroutine particle_init (ptcls, side)
        real :: side    ! cube side
        type(particle), dimension(:) :: ptcls  ! set of particles
        real :: lspc    ! lattice spacing
        integer :: i,j,k
        real :: temp
        integer :: npsd     ! particles per side
        real :: x,y,z
        
        ! ATTENTION : may particles go outside the box?
        temp = size(ptcls)**(1/3.0)
        lspc = side/(temp+1)
        npsd = int(temp)
        
        
    end subroutine particle_init
    
end module LJ_global


program LJ

    ! use jackknife, only: JK, JK_init, JK_cluster, JK_function
    use LJ_global, only: particle, particle_init, L, N
    implicit none
    
    type(particle), dimension(:), allocatable :: ptcls
    
    print *, "Box side: "
    read *, L
    print *, "Number of particles: "
    read *, N
    

end program LJ
