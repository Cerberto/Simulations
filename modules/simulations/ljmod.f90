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
    end type particle
    
    type(particle), dimension(:), allocatable :: ptcls
    
    real(dp) :: side        ! side of the cubic box
    integer  :: N           ! number of particles
    real(dp) :: sigma       ! unit of length (0 of the LJ potential)
    real(dp) :: eps         ! unit of energy (depth of the LJ potential well)
    real(dp) :: delta       ! (twice the) maximum displacement in the metropolis
    real(dp) :: core        ! core radius (to avoid divergences) 
    
    ! potential energy and kinetic energy (virial)
    real(dp) :: poten, kinen
    
    ! auxiliary vector (position updating) for the metropolis
    real(dp), dimension(3) :: pstn_new
    
contains
    
    
    !
    !   Initialization of the particles on a cubic lattice inside a cubic box
    !   and calculation of the initial value of the potential energy
    !
    subroutine particle_init (ptcls)
        
        type(particle), dimension(:) :: ptcls  ! set of particles
        real(dp) :: lspc    ! lattice spacing
        integer  :: i,j,k
        real(dp) :: temp
        integer  :: npsd     ! particles per side
        real(dp), dimension(3) :: x0, x1
    
        temp = N**(1/3.0)
        npsd = int(temp)+1
        lspc = side/npsd
        print *, "side = ", side
        print *, "lspc = ", lspc
        print *, "npsd = ", npsd
        
        ! initialization of particles inside a cubic box centered in 0:
        ! particles are distributed on a cubic lattice with spacing lspc
        !counter = 0
        x0 = (/ -side/2 + lspc/2, -side/2 + lspc/2, -side/2 + lspc/2 /)
        do i=0, npsd-1, 1
            do j=0, npsd-1, 1
                do k=0, npsd-1, 1
                    x1 = (/ x0(1) + i*lspc, x0(2) + j*lspc, x0(3) + k*lspc /)
                    ptcls(i*npsd**2 + j*npsd + k + 1)%pstn = x1
                    if (i*npsd**2 + j*npsd + k + 1 == N) goto 10
                end do
            end do
        end do
    
10      poten = 0
        do i=2, N, 1
            do j=1, i-1, 1
                poten = poten + &
                  pair_interaction(ptcls(i)%pstn, ptcls(j)%pstn)
            end do
        end do
        
    end subroutine particle_init
    
    
    !
    !   Pair interaction
    !
    function pair_interaction (vec1, vec2)
        
        real(dp) :: pair_interaction, r, t1, t2
        real(dp), dimension(:) :: vec1, vec2
        integer :: i
        
        r = 0
        do i=1, 3, 1
            t1 = (vec2(i) - vec1(i))*2/side
            call rintf(t1, t2)
            r = r + (vec1(i) - vec2(i) + side*t1)**2
        end do
        r = sqrt(r)
        if (r>side/2) then
            pair_interaction = 0
            return
        else if (r<core) then
            r = core
        end if
        pair_interaction = 4*eps*((sigma/r)**12 - (sigma/r)**6)
        
    end function pair_interaction
    
    
    !
    !   Interaction energy difference when displacing the k-th particle
    !
    function delta_interaction (ptcls, k)
        
        real(dp) :: delta_interaction
        real(dp) :: t1, t2
        type(particle), dimension(:) :: ptcls
        integer :: k, i
        
        t2 = 0
        delta_interaction = 0
        do i=1, N, 1
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
        
        real(dp) :: pair_virial, r, t1, t2
        real(dp), dimension(:) :: vec1, vec2
        integer :: i
        
        r = 0
        do i=1, 3, 1
            t1 = (vec2(i) - vec1(i))*2/side
            call rintf(t1, t2)
            r = r + (vec1(i) - vec2(i) + side*t1)**2
        end do
        r = sqrt(r)
        if (r>side/2) then
            pair_virial = 0
            return
        else if (r<core) then
            r = core
        end if
        pair_virial = 12*eps*(-2*(sigma/r)**12 + (sigma/r)**6)
        
    end function pair_virial
    

    !
    !   Interaction energy difference when displacing the k-th particle
    !
    function delta_virial (ptcls, k)
    
        real(dp) :: delta_virial
        real(dp) :: t1, t2
        type(particle), dimension(:) :: ptcls
        integer :: k, i
    
        t2 = 0
        delta_virial = 0
        do i=1, N, 1
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
    
        real(dp) :: ke_virial, t1, t2, rsq
        type(particle), dimension(:) :: ptcls
        integer :: i,j,k
    
        ke_virial = 0
        do i=2, N, 1
            do j=1, i-1, 1
                ke_virial = ke_virial + &
                  pair_virial(ptcls(i)%pstn, ptcls(j)%pstn)
            end do
        end do
    
    end function ke_virial


end module ljmod
