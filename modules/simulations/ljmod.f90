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
    real(dp) :: sigma
    real(dp) :: eps
    real(dp) :: delta
    real(dp) :: potential
    real(dp), dimension(3) :: pstn_new

!    interface
!        function func (x, a, b)
!            integer :: a,b
!            type(particle), dimension(:) :: x
!            real(dp) :: func
!        end function func
!    end interface

end module ljmod

