program LJ
    
    use kinds,      only: dp
    use jackknife,  only: JK, JK_init, JK_cluster
    use ljmod
    use ljmetr,     only: thmetropolis, autocorrelation, nth, nsw, ndat, skip
    implicit none
    
    type(particle), dimension(:), allocatable :: ptcls
    type(JK) :: k_en, p_en
    real(dp), dimension(:), allocatable :: k_en_array, p_en_array
!    real(dp), external :: ke_virial
    real(dp) :: sum1, sum2
    integer :: sw       ! label for the 'sweeps' of the metropolis
    integer :: tmax     ! maximum time for the autocorrelation
    integer :: i, counter, check
    
    call rlxdinit(1,rand(time()))
    
    print *, "Number of particles: "
    read *, N

    nth     = 200
    nsw     = 100000
    skip    = 100
    tmax    = 50
    
    side    = (2.0*N)**(1/3.0)
    delta   = 0.01
    eps     = 1.0
    sigma   = 1.0
    
    allocate(ptcls(N))
    call particle_init(ptcls)

    print *, poten
    !
    !   Thermalization
    !
    do sw=1, nth, 1
        do i=1, N, 1
            call thmetropolis(ptcls)
        end do
        print *, sw, poten
        !
        ! print somewhere the potential, to show the thermalization process:
        !
    end do
    
    kinen = ke_virial(ptcls)
    ndat = nsw/skip
    allocate(k_en_array(ndat))
    allocate(p_en_array(ndat))
    counter = 1
    do sw=1, nsw, 1
        do i=1, N, 1
            call thmetropolis(ptcls)
        end do
        
        if (mod(sw,skip) == 0) then
            k_en_array(counter) = kinen/N
            p_en_array(counter) = poten/N
            counter = counter + 1
        end if
    end do
    
    sum1 = 0
    sum2 = 0
    do counter=1, ndat, 1
        sum1 = sum1 + k_en_array(counter)
        sum2 = sum2 + (k_en_array(counter))**2
    end do
    sum1 = sum1/ndat
    sum2 = sum2/(ndat-1)
    
    !call JK_init(k_en, ndat)
    !call JK_init(p_en, ndat)
    
    !k_en%vec = k_en_array
    !p_en%vec = p_en_array
    
    !call JK_cluster(k_en)
    !call JK_cluster(p_en)
    
    print *, 'Energy per particle        : ', sum1
    print *, 'Specific heat per particle : ', sum2

end program LJ
