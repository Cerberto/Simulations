program LJ
    
    use kinds,      only: dp
    use jackknife,  only: JK, JK_init, JK_cluster
    use ljmod,      only: particle, ptcls, N, side, sigma, eps, delta, core, &
                          poten, particle_init
    use ljmetr,     only: thmetropolis, autocorrelation, nth, nsw, ndat, skip
    implicit none
    
    
    type(JK) :: p_en !, k_en
    real(dp), dimension(:), allocatable :: p_en_array !, k_en_array
!    real(dp), external :: ke_virial
    real(dp) :: sum1, sum2
    integer :: sw       ! label for the 'sweeps' of the metropolis
    integer :: tmax     ! maximum time for the autocorrelation
    integer :: i, counter, check
    
    call rlxdinit(1,rand(time()))
    
    print *, "Number of particles: "
    read *, N

    nth     = 2000
    nsw     = 100000
    skip    = 100
    tmax    = 50
    
    side    = (2.0*N)**(1/3.0)
    delta   = 0.01
    eps     = 1.0e-6
    sigma   = 1.0
    core    = 0
    
    allocate(ptcls(N))
    call particle_init(ptcls)
    
    open (unit=8, file="output/particle_init.dat", status="replace", &
        action="write")
    do i=1, N, 1
        write (unit=8, fmt=*), &
            ptcls(i)%pstn(1), ptcls(i)%pstn(2), ptcls(i)%pstn(3)
    end do
    call flush(8)
    close (unit=8, iostat=check, status="keep")
    print *, "File closure: ", check

    print *, poten
    !
    !   Thermalization
    !
    do sw=1, nth, 1
        do i=1, N, 1
            call thmetropolis(ptcls)
        end do
        
        if (mod(sw,200) .eq. 0) then
            print *, sw, poten
        end if
        !
        ! print somewhere the potential, to show the thermalization process:
        !
    end do
    
    ! kinen = ke_virial(ptcls)
    ndat = nsw/skip
    ! allocate(k_en_array(ndat))
    allocate(p_en_array(ndat))
    counter = 1
    do sw=1, nsw, 1
        do i=1, N, 1
            call thmetropolis(ptcls)
        end do
        
        if (mod(sw,skip) == 0) then
            ! k_en_array(counter) = kinen/N
            p_en_array(counter) = poten/N
            counter = counter + 1
        end if
    end do

    
    ! call JK_init(k_en, ndat)
    call JK_init(p_en, ndat)
    
    ! k_en%vec = k_en_array
    p_en%vec = p_en_array
    
    ! call JK_cluster(k_en)
    call JK_cluster(p_en)
    
    print *, 'Energy per particle        : ', p_en%mean

end program LJ
