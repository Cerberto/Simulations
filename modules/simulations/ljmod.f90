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

!    interface
!        function func (x, a, b)
!            integer :: a,b
!            type(particle), dimension(:) :: x
!            real(dp) :: func
!        end function func
!    end interface

end module ljmod


function distance (ptcls, i, j)
    use kinds, only: dp
    use ljmod, only: particle 
    type(particle), dimension(:) :: ptcls
    real(dp) :: distance
    integer :: i,j,k
    
    distance = 0
    do k=1, 3, 1
        distance = distance + (ptcls(i)%pstn(k) - ptcls(j)%pstn(k))**2
    end do

end function distance
    

subroutine particle_init (ptcls)
    use kinds, only: dp
    use ljmod, only: side, N, particle
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


function lj_potential (r)
    use kinds, only: dp
    use ljmod, only: sigma, eps
    implicit none
    
    real(dp) :: lj_potential, r
    
    lj_potential = 4.0*eps*((sigma/r)**12.0 - (sigma/r)**6.0)
end function lj_potential


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
            temp = distance(ptcls,i,j)
            interaction = interaction + lj_potential(temp)
        end do
    end do

end function interaction

function delta_interaction (ptcls, k) result(diff)
    use kinds, only: dp
    use ljmod, only: particle, N
    
    real(dp) :: diff
    type(particle), dimension(:) :: ptcls
    integer :: k
    

end function delta_interaction

function thmetropolis (ptcls, k)
    use kinds, only: dp
    use ljmod, only: particle, delta
    implicit none
    
    real(dp) :: thmetropolis = 0
    type(particle), dimension(:) :: ptcls
    integer :: i,k
    real(dp), dimension(:), allocatable :: displ, temp
    real(dp) :: t1, t2
    
    allocate(u(4))
    allocate(temp(3))
    call ranlxdf(u,4)
    
    do i=1, 3, 1
        temp(i) = ptcls(k)%pstn(i) + delta*(2*u(i)-1)
        temp(i) = temp(i) - side*rintf(temp(i)/side)
    end do
    
    t1 = delta_interaction(ptcls,k)
    t2 = exp(-t1)
        
    if(t2 >= u(4)) then
        ptcls(k)%pstn = temp        ! SI PUÒ FARE?!
        potential = potential + t1
    end if
    
    deallocate(u)
    deallocate(temp)

end function thmetropolis


double HOmetropolis (double* state, int state_dim)
{
	int i;
	double u[2];
	double x_new, temp1, temp2;
	double DS = 0;
	
	for(i=0; i<state_dim; i++)
	{
		ranlxd(u,2);
		x_new = state[i] + DELTA*(2*u[1] - 1);
		temp1 = deltaS(state,state_dim,x_new,i);
		temp2 = exp(-temp1);

		/* scelta nuova configurazione */
		if(temp2 >= u[0])
		{
			state[i] = x_new;
			DS += temp1;
		}
	}
	return DS;
}
