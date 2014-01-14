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
    
    real(dp) :: side        ! side of the cubic box
    integer :: N            ! number of particles
    real(dp) :: sigma       ! unit of length (0 of the LJ potential)
    real(dp) :: eps         ! unit of energy (depth of the LJ potential well)
    real(dp) :: delta       ! (twice the) maximum displacement in the metropolis 
    
    ! potential energy and kinetic energy (virial)
    real(dp) :: poten, kinen
    
    ! auxiliary vector (position updating) for the metropolis
    real(dp), dimension(3) :: pstn_new

end module ljmod
