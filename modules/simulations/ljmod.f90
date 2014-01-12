module ljmod
    implicit none
    public
    
    type :: particle
        real :: x, y, z
    end type particle
    
    real :: side        ! cube side
    integer :: N        ! number of particles
    real, parameter :: sigma = 1.0      ! length scale
    real, parameter :: eps = 1.0        ! energy scale

end module ljmod


subroutine particle_init (ptcls)
    use ljmod, only: side, N, particle
    implicit none
    
    type(particle), dimension(:) :: ptcls  ! set of particles
    real :: lspc    ! lattice spacing
    integer :: i,j,k
    real :: temp
    integer :: npsd     ! particles per side
    real :: x,y,z
    
    ! ATTENTION : may particles go outside the box? check with numbers
    temp = N**(1/3.0)
    lspc = side/(temp+2)
    npsd = int(temp)+1
    
    ! initialization of particles inside a cubic box centered in 0:
    ! particles are distributed on a cubic lattice with spacing lspc
    x = -side + lspc*npsd/2
    y = x
    z = x
    do i=0, npsd-1, 1
        x = x + lspc
        do j=0, npsd-1, 1
            y = y + lspc
            do k=0, npsd-1, 1
                z = z + lspc
                ptcls(i+j+k)%x = x
                ptcls(i+j+k)%y = y
                ptcls(i+j+k)%z = z
            end do
        end do      
    end do
end subroutine particle_init


function distance (ptcls, i, j)
    use ljmod, only: particle
    implicit none
    
    type(particle), dimension(:) :: ptcls
    real :: distance
    integer :: i,j
    distance = (ptcls(i)%x - ptcls(j)%x)**2
    distance = distance + (ptcls(i)%y - ptcls(j)%y)**2
    distance = distance + (ptcls(i)%z - ptcls(j)%z)**2
end function distance


function lj_potential (r)
    use ljmod, only: sigma, eps
    implicit none
    
    real :: lj_potential, r
    lj_potential = 4.0*eps*((sigma/r)**12 - (sigma/r)**6)
end function lj_potential


function interaction (ptcls)
    use ljmod, only: particle, N
    implicit none
    
    real :: interaction, temp
    type(particle), dimension(:) :: ptcls
    integer :: i,j
    
    interaction = 0
    do i=1, N-1, 1
        do j=0, i-1, 1
            temp = distance(ptcls,i,j)
            interaction = interaction + lj_potential(temp)
        end do
    end do

end function interaction
