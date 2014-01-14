!
!   Pair interaction
!
function pair_interaction (vec1, vec2)
    use kinds, only: dp
    use ljmod, only: sigma, eps, side
    
    real(dp) :: pair_interaction, r, t1, t2
    real(dp), dimension(3) :: vec1, vec2
    integer :: i
    
    r = 0
    do i=0, 3, 1
        t1 = (vec2(i) - vec1(i))*2/side
        call rintf(t1, t2)
        r = r + (vec1(i) - vec2(i) + side*t1)**2
    end do
    r = sqrt(r)
    pair_interaction = 4*eps*((sigma/r)**12 - (sigma/r)**6)
    
end function pair_interaction


!
!   Initialization of the particles on a cubic lattice inside a cubic box
!   and calculation of the initial value of the potentials
!
subroutine particle_init (ptcls)
    use kinds, only: dp
    use ljmod, only: side, N, particle, potential
    implicit none
    
    type(particle), dimension(:) :: ptcls  ! set of particles
    real(dp), external :: pair_interaction
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
    
    potential = 0
    do i=1, N-1, 1
        do j=0, i-1, 1
            potential = potential + &
                pair_interaction(ptcls(i)%pstn, ptcls(j)%pstn)
        end do
    end do
    
end subroutine particle_init


!
!   Interaction energy difference when displacing the k-th particle
!
function delta_interaction (ptcls, k)
    use kinds, only: dp
    use ljmod, only: particle, N, pstn_new
    implicit none
    
    real(dp) :: delta_interaction
    real(dp) :: t1, t2
    real(dp), external :: pair_interaction
    type(particle), dimension(N) :: ptcls
    integer :: k, i
    
    t2 = 0
    delta_interaction = 0
    do i=0, N, 1
        if (i /= k) then
            t1 = pair_interaction(ptcls(i)%pstn, ptcls(k)%pstn)
            t2 = pair_interaction(ptcls(i)%pstn, pstn_new)
            delta_interaction = delta_interaction - t1 + t2
        end if
    end do

end function delta_interaction


!
!   Routine implementing the Metropolis algorithm for a thermal distribution.
!
subroutine thmetropolis (ptcls, k)
    use kinds, only: dp
    use ljmod, only: particle, delta, pstn_new, potential, side
    implicit none
    
    type(particle), dimension(:) :: ptcls
    integer :: i,k
    real(dp), dimension(:), allocatable :: u
    real(dp) :: t1, t2
    real(dp), external :: delta_interaction
    
    allocate(u(4))
    call ranlxdf(u,4)
    
    do i=1, 3, 1
        pstn_new(i) = ptcls(k)%pstn(i) + delta*(2*u(i)-1)
        t1 = pstn_new(i)/side
        call rintf(pstn_new(i), t2)
        pstn_new(i) = pstn_new(i) - side*t2
    end do
    
    t1 = delta_interaction(ptcls,k)
    t2 = exp(-t1)
        
    if(t2 >= u(4)) then
        ptcls(k)%pstn = pstn_new
        potential = potential + t1
    end if
    
    deallocate(u)

end subroutine thmetropolis
