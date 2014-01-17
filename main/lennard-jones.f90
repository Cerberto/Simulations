program LJ
    
    use kinds,      only: dp
    use jackknife,  only: JK, JK_init, JK_cluster ! , JK_function
    use ljmod,      only: particle, side, N, eps, sigma, delta, poten, kinen
    use ljmetr,     only: thmetropolis, autocorrelation, nth, nsw, ndat, skip
    implicit none
    
    type(particle), dimension(:), allocatable :: ptcls
    type(JK) :: k_en, p_en
    real(dp), external :: ke_virial
    integer :: sw       ! label for the 'sweeps' of the metropolis
    integer :: tmax     ! maximum time for the autocorrelation
    integer :: i, counter, check
    
    call rlxdinit(1,rand(time()))
    
    print *, "Number of particles: "
    read *, N
    
    allocate(ptcls(N))
    call particle_init(ptcls,side)

    nth     = 200
    nsw     = 1000000
    skip    = 100
    tmax    = 50
    
    side    = 100
    delta   = 0.001
    eps     = 1.0
    sigma   = 1.0
    side    = 1000.0
    
    !
    !   Thermalization
    !
    do sw=1, nth, 1
        do i=1, N, 1
            call thmetropolis(ptcls)
        end do
        !
        ! print somewhere the potential, to show the thermalization process:
        !
    end do
    
    kinen = ke_virial(ptcls)
    ndat = nsw/skip
    call JK_init(k_en, ndat)
    call JK_init(p_en, ndat)
    counter = 1
    do sw=1, nsw, 1
        do i=1, N, 1
            call thmetropolis(ptcls)
        end do
        
        if (mod(sw,skip) == 0) then
            k_en%vec(counter) = kinen
            p_en%vec(counter) = poten
            counter = counter + 1
        end if
    end do
    
    call JK_cluster(k_en)
    call JK_cluster(p_en)
    
    print *, 'Energy per particle        : ', k_en%mean
    print *, 'Specific heat per particle : ', k_en%var

end program LJ
