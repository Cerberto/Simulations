program LJ
    
    use kinds, only: dp
    use jackknife, only: JK, JK_init, JK_cluster ! , JK_function
    use ljmod, only: particle, side, N, eps, sigma, delta
    implicit none
    
    type(particle), dimension(:), allocatable :: ptcls
    
    rlxdinit(1,rand(time()))
    
    print *, "Box side: "
    read *, side
    print *, "Number of particles: "
    read *, N
    
    allocate(ptcls(0:N-1))
    call particle_init(ptcls,side)
    
    

end program LJ
