
!-------------------------------------------------------------------------------
!
!   Simulation of a Lennard-Jones fluid with the Metropolis algorithm.
!   The program computes the potential energy per particle and the specific heat
!   per particle.
!
!   It is possible to choose to perform calculations at constant volume or
!   pressure. In the main routine change
!       'thmetropolis_v' into 'thmetropolis_p' or viceversa
!
!-------------------------------------------------------------------------------

program LJ
    
    use kinds,      only: dp
    use jackknife,  only: JK, JK_init, JK_cluster
    use ljmod,      only: N, side, sigma, eps, delta, beta, poten, rho, &
                          rcutoff, particle, ptcls, particle_init, &
                          total_interaction
    use ljmetr,     only: thmetropolis_v, thmetropolis_p, thmetropolis_p_alt, &
                          nth, nsw, ndat, nbin
    implicit none
    
    
    type(JK) :: p_en, cv

    real(dp) :: sum_p_en, sum_cv, acpt_rate, temperature, tempinit, deltatemp
    integer :: sw, tmax, i, counter
    
    ! Positions of particles initialized
    open (unit=8, file='lj_output/V_particle_init.dat', status='replace', &
        action='write')
        
    ! Positions of particles thermalized
    open (unit=9, file='lj_output/V_particle_therm.dat', status='replace', &
        action='write')
    
    ! Potential vs markovian time
    open (unit=10, file='lj_output/V_potential.dat', status='replace', &
        action='write')
        
    ! Energy and cv (with respective errors) vs temperature
    !open(unit=12, file='lj_output/V_X_vs_T.dat', status='replace', action='write')
    open(unit=12, file='lj_output/V_X_vs_T.dat', access='append', action='write')
    
    call rlxdinit(1,rand(time()))
    
    read *, N

    read *, nth
    read *, nsw
    read *, nbin
    read *, tmax
    
    ndat = nsw/nbin

    read *, rho
    read *, delta
    read *, eps
    read *, sigma
    ! read *, beta
    read *, tempinit
    read *, deltatemp
    
    allocate(ptcls(N))
    call JK_init (p_en,ndat)
    call JK_init (cv,ndat)
    
temperature = tempinit
side    = (N/rho)**(1/3.0)
rcutoff = side/2.d0
! DO LOOP ON TEMPERATURE
do while (temperature <= 8.d0)
    
    write (6,*) 'Number of particles     :', N
    write (6,*) 'Volume                  :', side**(3.d0)
    write (6,*) 'Temperature             :', temperature
    
    beta    = temperature**(-1.0)
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
        acpt_rate = acpt_rate + thmetropolis_v(ptcls)/nth
        !if (mod(sw,10)==0) then
        !    write (10,*) sw, poten
        !end if
    end do
    write (6,*) 'Acceptance rate (in thermalization) :', acpt_rate
    
    !
    !   Print position of particles assumed thermalized
    !
    do i=1, N, 1
        write (9,*), ptcls(i)%pstn
    end do
    
    counter = 1
    acpt_rate = 0
    sum_cv = 0
    sum_p_en = 0
    
    do sw=1, nsw
        acpt_rate = acpt_rate + thmetropolis_v(ptcls)/nsw
        
        sum_p_en = sum_p_en + poten/nbin
        sum_cv = sum_cv + poten**2/nbin
        if (mod(sw,nbin) == 0) then
        !    write (10,*) (sw+nth)/nbin, poten
            p_en%vec(counter) = sum_p_en
            sum_p_en = 0
            cv%vec(counter) = beta*beta*(sum_cv - sum_p_en*sum_p_en)
            sum_cv = 0
            counter = counter + 1
        end if
    end do
    write (6,*) 'Acceptance rate (when thermalized) :', acpt_rate
    
    !
    ! Compute mean and variance (of the mean) of potential energy
    !
    call JK_cluster (p_en)
    call JK_cluster (cv)
    
    !
    ! Compute mean and variance (of the mean) of the specific heat
    !
!    do counter=1, ndat
!        cv%vec(counter) = (beta**2)*(cv%vec(counter) - p_en%mean**2)
!    end do
!    call JK_cluster (cv)
    
    write (6,*) 'Energy / particle       :', p_en%mean/N, '+-', sqrt(p_en%var)/N
    write (6,*) 'Specific heat / particle:', cv%mean/N, '+-', sqrt(cv%var)/N
    write (6,*) ' '
    call flush (6)
    
    write (12,*) temperature, p_en%mean/N, sqrt(p_en%var)/N, &
                              cv%mean/N, sqrt(cv%var)/N

    temperature = temperature + deltatemp
end do

    close (8)
    close (9)
    close (10)
    close (12)

end program LJ


!
! Super efficient debugging
!
subroutine pluto()
    write (6,*) 'Pluto!'
    call flush (6)
end subroutine pluto
