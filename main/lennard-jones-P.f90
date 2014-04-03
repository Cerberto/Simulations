
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
    use ljmod,      only: N, side, sigma, eps, delta, beta, poten, rho, press, &
                          rcutoff, particle, ptcls, particle_init, &
                          total_interaction, dv
    use ljmetr,     only: thmetropolis_v, thmetropolis_p, thmetropolis_p_alt, &
                          nth, nsw, ndat, nbin
    implicit none
    
    
    type(JK) :: p_en, cv, vol
    real(dp), dimension(:), allocatable :: p_en_array, vol_array
    
    real(dp) :: sum_p_en, sum_cv, sum_vol, acpt_rate, deltainit, deltapress, &
                pressinit, ratio_pressdelta
    integer :: sw, tmax, i, j, counter
    character (len=100) :: output_filename
    integer :: output_channel
    
!    open (unit=8, file='lj_output/P_particle_init.dat', status='replace', &
!        action='write')
!    open (unit=9, file='lj_output/P_particle_therm.dat', status='replace', &
!        action='write')
!    open (unit=10, file='lj_output/P_potential.dat', status='replace', &
!        action='write')

    call rlxdinit(1,rand(time()))
    
    read *, N
    
    read *, nth
    read *, nsw
    read *, nbin
    read *, tmax
    
    ndat = nsw/nbin
    
    read *, rho
    read *, deltainit
    read *, eps
    read *, sigma
    read *, beta
    read *, pressinit
    read *, dv
    read *, deltapress
    
!    write (10,*) '# Potential during thermalization process '
!    write (10,*) '# p     =', press
!    write (10,*) '# delta =', deltainit
!    write (10,*) '# dv    =', dv
    
    allocate(ptcls(N))
    allocate(p_en_array(ndat))
    allocate(vol_array(ndat))
    call JK_init (p_en,ndat)
    call JK_init (cv,ndat)
    call JK_init (vol,ndat)
    
    ! 
    !   Do loop on pressure
    !
    side    = (N/rho)**(1/3.0)
    rcutoff = side/2.d0
    press   = pressinit
    delta   = deltainit
    ratio_pressdelta = deltainit/pressinit
    j=0
    do while (press < 0.51)

        delta = deltainit*(press/pressinit)**(0.33333)
        call particle_init(ptcls)
        
        !
        !   Print position of particles right after initialization
        !
        do i=1, N
            write (8,*), ptcls(i)%pstn
        end do
        
        !
        !   Open file where saving the potential as a function of the markovian
        !   time
        !
        write(output_filename,18) press
18      format('lj_output/P_potential_',f5.3)
        output_filename = trim(output_filename)
        output_channel  = 21 + 2*j
        open (unit=output_channel, file=output_filename, &
            status='replace', action='write')
        
        !
        !   Thermalization
        !
        acpt_rate = 0
        do sw=1, nth
            acpt_rate = acpt_rate + thmetropolis_p_alt(ptcls)/nth
            write (output_channel,*) sw, poten
        end do
        write (6,*) 'Acceptance rate (in thermalization) :', acpt_rate
        
        !
        !   Print position of particles assumed thermalized
        !
        do i=1, N, 1
            write (9,*), ptcls(i)%pstn
        end do
        
        counter   = 1
        acpt_rate = 0
        sum_p_en  = 0
        sum_cv    = 0
        sum_vol   = 0
        do sw=1, nsw
            acpt_rate = acpt_rate + thmetropolis_p_alt(ptcls)/nsw
            write (output_channel,*) sw+nth, poten
            
            sum_p_en = sum_p_en + poten/nbin
            sum_cv   = sum_cv + poten**2/nbin
            sum_vol  = sum_vol + side**3/nbin
            if (mod(sw,nbin) == 0) then
                p_en%vec(counter) = sum_p_en
                sum_p_en = 0
                cv%vec(counter) = beta*beta*(sum_cv - sum_p_en*sum_p_en)
                sum_cv = 0
                vol%vec(counter) = sum_vol
                sum_vol = 0
                counter = counter + 1
            end if
        end do
        write (6,*) 'Acceptance rate (when thermalized) :', acpt_rate
    
        call flush (output_channel)
        close (output_channel)
        
        !
        ! Compute mean and variance (of the mean) of energy, specific heat and volume
        !
        call JK_cluster (p_en)
        call JK_cluster (cv)
        call JK_cluster (vol)
        
        write (6,*) 'Number of particles     :', N
        write (6,*) 'Pressure                :', press
        write (6,*) 'Temperature             :', 1.d0/beta
        write (6,*) 'Energy / particle       :', p_en%mean/N, '+-', sqrt(p_en%var)/N
        write (6,*) 'Specific heat / particle:', cv%mean/N, '+-', sqrt(cv%var)/N
        write (6,*) 'Volume                  :', vol%mean, '+-', sqrt(vol%var)
        write (6,*) '<side>/2 - rcutoff      :', ((vol%mean)**(1.d0/3))/2 - rcutoff
        write (6,*) '<side> - rcutoff        :', (vol%mean)**(1.d0/3) - rcutoff
        write (6,*) ' '
        
        
        !
        ! Open (in append mode) a file for each value of the pressure where saving
        ! the average quantities
        !
        write(output_filename,19) press
19      format('lj_output/P_en-cv-vol_',f5.3)
        output_filename = trim(output_filename)
        output_channel  = 21 + 2*j + 1
        open(unit=output_channel, file=output_filename, access='append', &
            action='write')
        write (output_channel,*) press, p_en%mean/N, sqrt(p_en%var)/N, &
                                        cv%mean/N, sqrt(cv%var)/N, &
                                        vol%mean, sqrt(vol%var)
        
        call flush (output_channel)
        close (output_channel)
        
        j = j+1
        press = pressinit + real(j,dp)*deltapress
        call flush (6)
    end do
        
    close (8)
    close (9)
    close (10)
        
end program LJ

!
! Super efficient debugging
!
subroutine pluto()
    write (6,*) 'Pluto!'
    call flush (6)
end subroutine pluto
