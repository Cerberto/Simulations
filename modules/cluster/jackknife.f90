
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
    use kinds, only: dp
    implicit none
    public

    ! definition of the jackknife cluster
    type JK
        real(dp), dimension(:), allocatable :: vec  ! sample vector
        integer :: D                        ! dimension of the sample
        real(dp) :: mean, var                   ! mean and variance of the mean
    end type JK

    ! interface for a real(dp) function of a real(dp) variable
!    interface
!        function func (x)
!            real(dp) :: x, func
!        end function
!    end interface
    
    ! interface for a JK cluster function
!    abstract interface
!        function jkf (f, X)
!            type (JK) :: jkf, X
!            procedure(func), pointer :: f
!            real(dp) :: temp
!            integer :: i
!        end function jkf
!    end interface    
    
contains

    ! Assignment of mean value and variance in a cluster structure
    subroutine JK_cluster (C)
        
        type(JK) :: C
        integer :: i,D
        
        D = C%D
        C%mean = 0
        do i=1, D, 1
            C%mean = C%mean + (C%vec(i))/(real(D, dp))
        end do
	
        do i=1, D, 1
            C%vec(i) = C%mean + (C%mean - C%vec(i))/(real(D-1, dp))
        end do
	
        C%var = 0
        do i=1, D, 1
            C%var = C%var + (C%vec(i) - C%mean)*(C%vec(i) - C%mean)
        end do
            C%var = C%var * (real(D-1, dp))/(real(D, dp))
    end subroutine JK_cluster


    subroutine JK_init (C, D)
        type(JK) :: C
        integer :: D
        
        C%D = D
        C%mean = 0
        C%var = 0
        allocate(C%vec(D))
    end subroutine JK_init


!    function JK_function (f, X) result(res)
        
        ! definition of arguments and result
!        type (JK) :: X, res
!        integer :: D
!        real(dp), external, pointer :: f
        ! ausiliary variables
!        integer :: i
!        real(dp) :: temp = 0
        
        ! initialization of the cluster for the result
!        D = X%D
!        call cluster_init(res,D);
        
        ! assignments
!        res%mean = f(X%mean);
!        do i=0, D-1, 1
!            res%vec(i) = f(X%vec(i));
!            temp = temp + (res%vec(i) - res%mean)*(res%vec(i) - res%mean);
!        end do
!        temp = temp * (real(dp)(D - 1)/real(dp)(D));
!        res%var = temp
        
!    end function JK_function

end module jackknife
