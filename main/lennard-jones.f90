module LJ_global
    implicit none
    
    type :: particle
        real :: x, y, z
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
        
        ! ATTENTION : may particles go outside the box? check with numbers
        temp = size(ptcls)**(1/3.0)
        lspc = side/(temp+2)
        npsd = int(temp)+1
        
        ! initialization of particles inside a cubic box centered in 0:
        ! particles are distributed on a cubic lattice with spacing lspc
        x = -side + lspc*npsd/2
        y = x
        z = x
        do i=1, npsd, 1
            x = x + lspc
            do j=1, npsd, 1
                y = y + lspc
                do k=1, npsd, 1
                    z = z + lspc
                    ptcls(i+j+k-2)%x = x
                    ptcls(i+j+k-2)%y = y
                    ptcls(i+j+k-2)%z = z
                end do
            end do      
        end do
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
    
    allocate(ptcls(N))
    call particle_init(ptcls,L)

end program LJ
