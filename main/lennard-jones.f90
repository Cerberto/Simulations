program LJ

    use jackknife, only: JK, JK_init, JK_cluster ! , JK_function
    use ljmod, only: particle, L, N, eps, sigma
    implicit none
    
    type(particle), dimension(:), allocatable :: ptcls
    
    print *, "Box side: "
    read *, L
    print *, "Number of particles: "
    read *, N
    
    allocate(ptcls(0:N-1))
    call particle_init(ptcls,L)

end program LJ
