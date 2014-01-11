module ljmod
    implicit none
    
    type :: particle
        real :: x, y, z
    end type particle
    
    real :: L           ! cube side
    integer :: N        ! number of particles
    real, parameter :: sigma = 1.0      ! length scale
    real, parameter :: eps = 1.0        ! energy scale

end module ljmod


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
    implicit none
    type(particle), dimension(:) :: ptcls
    real :: distance
    integer :: i,j
    distance = (ptcls(i)%x - ptcls(j)%x)**2
    distance = distance + (ptcls(i)%y - ptcls(j)%y)**2
    distance = distance + (ptcls(i)%z - ptcls(j)%z)**2
end function distance


function lj_potential (r) result(V)
    implicit none
    use ljmod, only: sigma, eps
    
    real :: V, r
    V = 4*eps*((sigma/r)**12 - (sigma/r)**6)
end function lj_potential


function potential_tot (ptcls) result(V)
    implicit none
    
    real :: V
    type(particle), dimension(:) :: ptcls
    integer :: i,j,npar
    
    npar = size(ptcls)
    V = 0
    do i=1, npar-1, 1
        do j=0, i-1, 1
            V = V + lj_potential(distance(ptcls,i,j))
        end do
    end do

end function interaction
