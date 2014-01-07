!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!		File "zeros.c"
!
! Files containing routines for finding zeros of functions
! 
! _ double Zbisection (double (*f)(double), double *v, double accuracy):
! 		uses bisection method with final uncertainty specified by "accuracy".
! 
! _ double Zsecant (double (*f)(double), double *v, double accuracy):
! 		uses secant method with final uncertainty specified by "accuracy".
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module find_zero
    use constants
    implicit none
    
    abstract interface
        function zero (f, v, acc)
            procedure(func), pointer :: f
            real, dimension(:), pointer :: v
            real :: acc
        end function zero
    end interface
    
    interface
        function func (x)
            real :: x, func
        end function func
    end interface
    
    
contains


    function Zbisection (f, v, acc)
        procedure(func), pointer :: f
        real, dimension(:), pointer :: v
        real :: acc , Zbisection
        
        integer :: i = 0
        real :: x, f0, f1, fx
        f0 = f(v(0))
        f1 = f(v(1))
        ! printf("\t\t x_min\t\t x_max\t\t f(x_min)\t f(x_max)\n");
        ! printf("Bisection %d :\t%e\t%e\t%e\t%e\n", i, v[0], v[1], f0, f1);
        print*, htab, htab, " x_min", htab, htab, " x_max", htab, htab, "f(x_min)", htab, "f(x_max)"
        print*, "Bisection ", i, " ", htab, v(0), htab, v(1), htab, f0, htab, f1
	
        if (f0==f1 || (f0*f1>0)) then
            ! error message and exit program
        endif
	
        do
            i = i+1
            x = (v(0)+v(1))/2.0
            fx = f(x)
		
            if((f0*fx)<0) then
            	v(1)=x
                f1 = fx
            else if((f0*fx)>0) then
                v(0)=x
                f0 = fx
            else
                exit
            endif
			
            if (i==200) then
                ! error message (no convergence) and exit program
            endif

            if (fabs(v(0)-v(1)) < accuracy) then ! OCCHIO A FABS
                exit
            endif
        end do        
        Zbisection = x
    end function Zbisection
    
    
    
    function Zsecant (f, v, acc)
        procedure(func), pointer :: f
        real, dimension(:), pointer :: v
        real :: acc , Zsecant
        
        integer :: i = 0
        real :: x, f0, f1, fx
        f0 = f(v(0))
        f1 = f(v(1))
        ! printf("\t\t x_min\t\t x_max\t\t f(x_min)\t f(x_max)\n");
        ! printf("Secant    %d :\t%e\t%e\t%e\t%e\n", i, v[0], v[1], f0, f1);
        ! ....
	
        if (f0==f1 || (f0*f1>0)) then
            ! error message and exit program
        endif
	
        do
            i = i+1
            x = v(0) + (v(0)-v(1))*f0/(f1-f0)
            fx = f(x)
		
            if((f0*fx)<0) then
            	v(1)=x
                f1 = fx
            else if((f0*fx)>0) then
                v(0)=x
                f0 = fx
            else
                exit
            endif
			
            if (i==200) then
                ! error message (no convergence) and exit program
            endif
            
            ! printf("Secant    %d :\t%e\t%e\t%e\t%e\n", i, v[0], v[1], f0, f1);

            if (fabs(v(0)-v(1)) < acc) then ! OCCHIO A FABS
                exit
            endif
        end do        
        Zsecant = x
    end function Zsecant
    
end module find_zero
