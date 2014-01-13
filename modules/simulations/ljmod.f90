!
!   Module containing global types and variables for the Lennard-Jones
!   Metropolis calculations
!
module ljmod
    use kinds, only: dp
    implicit none
    public
    
    type :: particle
        real(dp), dimension(3) :: pstn
        ! real(dp), dimension(3) :: vlct
    end type particle
    
    real(dp) :: side
    integer :: N
    real(dp), parameter :: sigma = 1.0
    real(dp), parameter :: eps = 1.0
    real(dp) :: delta, potential
    real(dp), dimension(3) :: pstn_new

!    interface
!        function func (x, a, b)
!            integer :: a,b
!            type(particle), dimension(:) :: x
!            real(dp) :: func
!        end function func
!    end interface

end module ljmod


!   
!   Distance between two points in 3D specified by vec1 and vec2
!
function distance (vec1, vec2)
    use kinds, only: dp
    use ljmod, only: particle 
    real(dp), dimension(3) :: vec1, vec2
    real(dp) :: distance
    integer :: k
    
    distance = 0
    do k=1, 3, 1
        distance = distance + (vec1(k) - vec2(k))**2
    end do

end function distance
    

!
!   Initialization of the particles on a cubic lattice inside a cubic box
!
subroutine particle_init (ptcls)
    use kinds, only: dp
    use ljmod, only: side, N, particle, potential
    implicit none
    
    type(particle), dimension(:) :: ptcls  ! set of particles
    real(dp) :: lspc    ! lattice spacing
    integer :: i,j,k
    real(dp) :: temp
    integer :: npsd     ! particles per side
    real(dp), dimension(3) :: x0 
    
    ! ATTENTION : may particles go outside the box? check with numbers
    temp = N**(1/3.0)
    lspc = side/(temp+2)
    npsd = int(temp)+1
    
    ! initialization of particles inside a cubic box centered in 0:
    ! particles are distributed on a cubic lattice with spacing lspc
    x0(1) = -side + lspc*npsd/2
    x0(2) = x0(1)
    x0(3) = x0(1)
    do i=0, npsd-1, 1
        x0(1) = x0(1) + lspc
        do j=0, npsd-1, 1
            x0(2) = x0(2) + lspc
            do k=0, npsd-1, 1
                x0(3) = x0(3) + lspc
                ptcls(i+j+k)%pstn = x0      ! SI PUÒ FARE?!
            end do
        end do      
    end do
    potential = interaction(ptcls)
    
end subroutine particle_init


!
!   Pair potential (Lennard-Jones type)
!
function lj_potential (r)
    use kinds, only: dp
    use ljmod, only: sigma, eps
    implicit none
    
    real(dp) :: lj_potential, r
    
    lj_potential = 4.0*eps*((sigma/r)**12.0 - (sigma/r)**6.0)
end function lj_potential


!
!   Total interaction energy
!
function interaction (ptcls)
    use kinds, only: dp
    use ljmod, only: particle, N !, distance
    implicit none
    
    real(dp) :: interaction, temp
    type(particle), dimension(:) :: ptcls
    integer :: i,j
    
    interaction = 0
    do i=1, N-1, 1
        do j=0, i-1, 1
            temp = distance(ptcls(i)%pstn, ptcls(j)%pstn)
            interaction = interaction + lj_potential(temp)
        end do
    end do

end function interaction


!
!   Interaction energy difference when displacing the k-th particle
!
function delta_interaction (ptcls, k) result(diff)
    use kinds, only: dp
    use ljmod, only: particle, N
    
    real(dp) :: diff = 0
    real(dp) :: t1, t2
    type(particle), dimension(:) :: ptcls
    integer :: k, i
    
    t2 = 0
    do i=0, N, 1
        if (i /= k) then
            t1 = distance(ptcls(i)%pstn, ptcls(k)%pstn)
            t1 = lj_potential(t1)
            t2 = distance(ptcls(i)%pstn, pstn_new)
            t2 = lj_potential(t2)
            
            diff = diff - t1 + t2
        end if
    end do

end function delta_interaction


!
!   Routine implementing the Metropolis algorithm for a thermal distribution.
!
subroutine thmetropolis (ptcls, k)
    use kinds, only: dp
    use ljmod, only: particle, delta, pstn_new, potential
    implicit none
    
    type(particle), dimension(:) :: ptcls
    integer :: i,k
    real(dp), dimension(:), allocatable :: u
    real(dp) :: t1, t2
    
    allocate(u(4))
    call ranlxdf(u,4)
    
    do i=1, 3, 1
        pstn_new(i) = ptcls(k)%pstn(i) + delta*(2*u(i)-1)
        pstn_new(i) = pstn_new(i) - side*rintf(pstn_new(i)/side)
    end do
    
    t1 = delta_interaction(ptcls,k)
    t2 = exp(-t1)
        
    if(t2 >= u(4)) then
        ptcls(k)%pstn = pstn_new        ! SI PUÒ FARE?!
        potential = potential + t1
    end if
    
    deallocate(u)

end subroutine thmetropolis
