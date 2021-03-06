
!-------------------------------------------------------------------------------
!
!   Module containing global types and variables for the Lennard-Jones
!   Metropolis calculations
!
!-------------------------------------------------------------------------------

module ljmod
    use kinds, only: dp
    implicit none
    
    public
    
    type :: particle
        real(dp), dimension(3) :: pstn
    end type particle
    
    type(particle), dimension(:), allocatable :: ptcls
    
    real(dp) :: side        ! side of the cubic box
    real(dp) :: rho         ! average density
    real(dp) :: lspc        ! lattice spacing
    integer  :: N           ! number of particles
    real(dp) :: sigma       ! unit of length (0 of the LJ potential)
    real(dp) :: eps         ! unit of energy (depth of the LJ potential well)
    real(dp) :: delta       ! (twice the) maximum displacement in the metropolis
    real(dp) :: beta        ! inverse temperature
    real(dp) :: press       ! pressure
    real(dp) :: rcutoff     ! cutoff radius
    real(dp) :: dv          ! maximum relative variation of side
    
    real(dp) :: poten       ! potential energy
    
    ! auxiliary vector (position updating) for the metropolis (canonical)
    real(dp), dimension(3) :: pstn_new
    
    ! auxiliary particle array for the metropolis (isoth-isob)
    type(particle), dimension(:), allocatable :: ptcls_new
    
contains
    
    
    !
    !   Initialization of the particles on a cubic lattice inside a cubic box
    !   and calculation of the initial value of the potential energy
    !
    subroutine particle_init (ptcls)
        
        type(particle), dimension(:) :: ptcls  ! set of particles
        integer  :: i,j,k,counter
        real(dp) :: temp
        integer  :: npsd     ! particles per side
        real(dp), dimension(3) :: x0, x1
    
        temp = N**(1/3.0)
        npsd = int(temp)+1
        lspc = side/npsd
        !write (6,*) "Initial side    = ", side
        !write (6,*) "Lattice spacing = ", lspc
        !write (6,*) "Cutoff radius   = ", rcutoff
        
        
        ! initialization of particles inside a cubic box centered in 0:
        ! particles are distributed on a cubic lattice with spacing lspc
        counter = 0
        x0 = (/ -side/2 + lspc/2, -side/2 + lspc/2, -side/2 + lspc/2 /)
        do i=0, npsd-1
            do j=0, npsd-1
                do k=0, npsd-1
                    counter = counter + 1
                    x1 = (/ x0(1) + i*lspc, x0(2) + j*lspc, x0(3) + k*lspc /)
                    ptcls(counter)%pstn = x1
                    if (counter == N) goto 10
                end do
            end do
        end do
        
10      poten = total_interaction(ptcls)
        
    end subroutine particle_init
    
    
    !
    !   Pair interaction
    !
    function pair_interaction (vec1, vec2)
        
        real(dp) :: pair_interaction, r2, t1, t2
        real(dp), dimension(:) :: vec1, vec2
        integer :: i
        
        pair_interaction = 0.d0
                
        r2 = 0.d0
        do i=1, 3
            t1 = (vec2(i) - vec1(i))/side
            call rintf(t1, t2)
            r2 = r2 + (vec1(i) - vec2(i) + side*t2)**2
        end do
        
        if (r2 < rcutoff**2) then
            pair_interaction = 4.d0*eps*(sigma**12/r2**6 - sigma**6/r2**3)
        end if
        
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
        do i=1, N
            if (i /= k) then
                t1 = pair_interaction(ptcls(i)%pstn, ptcls(k)%pstn)
                t2 = pair_interaction(ptcls(i)%pstn, pstn_new)
                delta_interaction = delta_interaction - t1 + t2
            end if
        end do
    end function delta_interaction


    !
    !   Total potential energy
    !
    function total_interaction (ptcls)
    
        real(dp) :: total_interaction
        type(particle), dimension(:) :: ptcls
        integer :: i,j
        total_interaction = 0
        do i=2, N
            do j=1, i-1
                total_interaction = total_interaction + &
                  pair_interaction(ptcls(i)%pstn, ptcls(j)%pstn)
            end do
        end do
    
    end function total_interaction

end module ljmod
