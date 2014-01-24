
!-------------------------------------------------------------------------------
!
!   Simulation of a Lennard-Jones fluid with the Metropolis algorithm.
!   The program computes the potential energy per particle and the specific heat
!   per particle.
!
!-------------------------------------------------------------------------------

program LJ
    
    use kinds,      only: dp
    use jackknife,  only: JK, JK_init, JK_cluster
    use ljmod,      only: N, side, sigma, eps, delta, beta, poten, rho, &
                          particle, ptcls, particle_init, total_interaction
    use ljmetr,     only: thmetropolis, nth, nsw, ndat, nbin
    implicit none
    
    
    type(JK) :: p_en, cv
    real(dp), dimension(:), allocatable :: p_en_array, cv_array

    real(dp) :: sum_p_en, sum_cv, acpt_rate
    integer :: sw       ! label for the 'sweeps' of the metropolis
    integer :: tmax     ! maximum time for the autocorrelation
    integer :: i, counter
    
    open (unit=8, file="output/particle_init.dat", status="replace", &
        action="write")
    open (unit=9, file="output/particle_therm.dat", status="replace", &
        action="write")
    open (unit=10, file="output/potential.dat", status="replace", &
        action="write")
    
    call rlxdinit(1,rand(time()))
    
    read *, N

    read *, nth
    read *, nsw
    read *, nbin
    read *, tmax

    read *, rho
    read *, delta
    read *, eps
    read *, sigma
    read *, beta
    
    acpt_rate = 0
    side    = (N/rho)**(1/3.0)
    
    allocate(ptcls(N))
    call particle_init(ptcls)
    
    !
    !   Print position of particles right after initialization
    !
    do i=1, N
        write (8,*), ptcls(i)%pstn
    end do
    call flush(8)
    close (unit=8)
    
    !
    !   Thermalization
    !
    acpt_rate = 0
    do sw=1, nth
        do i=1, N
            acpt_rate = acpt_rate + thmetropolis(ptcls)/(nth*N)
        end do
    end do
    print *, "Acceptance rate (in thermalization) :", acpt_rate

    
    !
    !   Print position of particles assumed thermalized
    !
    do i=1, N, 1
        write (9,*), ptcls(i)%pstn
    end do
    call flush (9)
    close (unit=9)
    
    ndat = nsw/nbin
    
    allocate(p_en_array(ndat))
    allocate(cv_array(ndat))
    counter = 1
    acpt_rate = 0
    sum_cv = 0
    sum_p_en = 0
    do sw=1, nsw, 1
        do i=1, N, 1
            acpt_rate = acpt_rate + thmetropolis(ptcls)/(nsw*N)
        end do
        
        sum_p_en = sum_p_en + poten/nbin
        sum_cv = sum_cv + poten**2
        if (mod(sw,nbin) == 0) then
            p_en_array(counter) = sum_p_en
            cv_array(counter) = sum_cv
            sum_p_en = 0
            sum_cv = 0
            counter = counter + 1
        end if
    end do
    write (6,*) 'Acceptance rate (when thermalized) :', acpt_rate
    call flush (10)
    close (10)
    
    !
    ! Compute mean and variance (of the mean) of the potential energy
    !
    call JK_init (p_en,ndat)
    p_en%vec = p_en_array
    call JK_cluster (p_en)
    
    do counter=1, ndat
        cv_array(counter) = (beta**2)*(cv_array(counter) - p_en%mean**2)
    end do
    
    !
    ! Compute mean and variance (of the mean) of the specific heat
    !
    call JK_init (cv,ndat)
    cv%vec = cv_array
    call JK_cluster (cv)
    
    print *, 'Energy per particle        :', p_en%mean/N, '+-', sqrt(p_en%var)/N
    print *, 'Specific heat per particle :', cv%mean/N, '+-', sqrt(cv%var)/N

end program LJ

!
! Super efficient debugging
!
subroutine pluto()
    write (6,*) "Pluto!"
    call flush (6)
end subroutine pluto
