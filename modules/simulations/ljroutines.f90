
!
!   Initialization of the particles on a cubic lattice inside a cubic box
!   and calculation of the initial value of the potential energy
!
subroutine particle_init (ptcls)
    use kinds, only: dp
    use ljmod, only: side, N, particle, poten
    implicit none
    
    type(particle), dimension(:) :: ptcls  ! set of particles
    real(dp), external :: pair_interaction
    real(dp) :: lspc    ! lattice spacing
    integer :: i,j,k
    real(dp) :: temp
    integer :: npsd     ! particles per side
    real(dp), dimension(3) :: x0 

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
    
    poten = 0
    do i=1, N-1, 1
        do j=0, i-1, 1
            poten = poten + &
                pair_interaction(ptcls(i)%pstn, ptcls(j)%pstn)
        end do
    end do
    
end subroutine particle_init


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
!   Pair interaction
!
function pair_virial (vec1, vec2)
    use kinds, only: dp
    use ljmod, only: sigma, eps, side
    
    real(dp) :: pair_virial, r, t1, t2
    real(dp), dimension(3) :: vec1, vec2
    integer :: i
    
    r = 0
    do i=0, 3, 1
        t1 = (vec2(i) - vec1(i))*2/side
        call rintf(t1, t2)
        r = r + (vec1(i) - vec2(i) + side*t1)**2
    end do
    r = sqrt(r)
    pair_virial = 24*eps*(-2*(sigma/r)**12 + (sigma/r)**6)
    
end function pair_virial


!
!   Interaction energy difference when displacing the k-th particle
!
function delta_virial (ptcls, k)
    use kinds, only: dp
    use ljmod, only: particle, N, pstn_new
    implicit none
    
    real(dp) :: delta_virial
    real(dp) :: t1, t2
    real(dp), external :: pair_virial
    type(particle), dimension(N) :: ptcls
    integer :: k, i
    
    t2 = 0
    delta_virial = 0
    do i=0, N, 1
        if (i /= k) then
            t1 = pair_virial(ptcls(i)%pstn, ptcls(k)%pstn)
            t2 = pair_virial(ptcls(i)%pstn, pstn_new)
            delta_virial = delta_virial - t1 + t2
        end if
    end do

end function delta_virial


!
!   Total kinetic energy, using the virial theorem
!
function ke_virial (ptcls)
    use kinds, only: dp
    use ljmod, only: particle, N, sigma, eps, side
    implicit none
    
    real(dp) :: ke_virial, t1, t2, rsq
    real(dp), external :: pair_virial
    type(particle), dimension(N) :: ptcls
    integer :: i,j,k
    
    ke_virial = 0
    do i=1, N, 1
        do j=0, i-1, 1
            ke_virial = ke_virial + &
                pair_virial(ptcls(i)%pstn, ptcls(j)%pstn)
        end do
    end do
    
end function ke_virial
