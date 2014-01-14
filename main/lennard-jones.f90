program LJ
    
    use kinds,      only: dp
    use jackknife,  only: JK, JK_init, JK_cluster ! , JK_function
    use ljmod,      only: particle, side, N, eps, sigma, delta
    use ljmetr,     only: thmetropolis, autocorrelation, nth, nsw, ndat, skip
    implicit none
    
    type(particle), dimension(:), allocatable :: ptcls
    
    call rlxdinit(1,rand(time()))
    
    print *, "Number of particles: "
    read *, N
    
    allocate(ptcls(0:N-1))
    call particle_init(ptcls,side)

    delta = 0.001
    eps = 1.0
    sigma = 1.0
    side = 1000.0
    
    

end program LJ
