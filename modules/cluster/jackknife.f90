
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!   Module containing definitions and routines for manipulating data using the
!   jackknife re-sampling method.
!
!   Routines are:
!
!   _ cluster_init -> initialization of jackknife cluster;
!
!   _ clusterJK -> assigns the mean and the variance of the mean.
!
!   _ functionJK -> returns a clusterJK for a secondary r.v.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module jackknife
    implicit none
    
    private none
    public

    ! definition of the jackknife cluster
    type :: JK
        real, dimension(:), allocatable :: vec  ! sample vector
        integer :: D                        ! dimension of the sample
        real :: mean, var                   ! mean and variance of the mean
    end type JK

    ! interface for a real function of a real variable
    interface
        function func (x)
            real :: x, func
        end function
    end interface
    
    ! interface for a JK cluster function
    abstract interface
        function jkf (f, X)
            type(JK) :: jkf, X
            procedure(func), pointer :: f
            real :: temp
            integer :: i
        end function jkf
    end interface    
    
contains

    ! Assignment of mean value and variance in a cluster structure
    subroutine JK_cluster (C)
        implicit none
        
        type(JK) :: C
        integer :: i
        integer :: D = C%D
        
        C%mean = 0
        do i=0, D-1, 1
            C%mean = C%mean + (C%vec(i))/(real(D))
        end do
	
        do i=0, D-1, 1
            C%vec(i) = C%mean + (C%mean - C%vec(i))/(real(D - 1))
        end do
	
        C%sigma = 0;
        do i=0, D-1, 1
            C%sigma = C%sigma + (C%vec(i) - C%mean)*(C%vec(i) - C%mean)
        end do
            C%sigma = C%sigma * (real(D - 1))/(real(D))
    end subroutine JK_cluster


    subroutine JK_init (C, D)
        implicit none
        type(JK) :: C
        integer :: D
        
        C%D = D
        C%mean = 0
        C%sigma = 0
        allocate(C%vec(0:D-1))
    end subroutine JK_init


    function JK_function (f, X) result(res)
        implicit none
        
        ! definition of arguments and result
        type(JK) :: X, res
        real, external, pointer :: f
        ! ausiliary variables
        integer :: i
        real :: temp = 0
        
        ! initialization of the cluster for the result
        call cluster_init(res,D);
        
        ! assignments
        res%mean = f(X%mean);
        do i=0, D-1, 1
            res%vec(i) = f(X%vec(i));
            temp = temp + (res%vec(i) - res%mean)*(res%vec(i) - res%mean);
        end do
        temp = temp * (real(D - 1)/real(D));
        res%sigma = temp
        
    end function JK_function

end module jackknife
