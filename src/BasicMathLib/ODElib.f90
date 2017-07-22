!i want to write some basic operations for ordinary differetial equations
!as the basic lib for Type<timedriver> or others
    

!nothing temporarily
module ODElib
use constants
implicit none

    
    private
    public::    ODEeuler
    public::    ODERK
    public::    absdydx
    
    
!-------------------------------------------------------------------    
    interface ODEeuler
        procedure:: ODEeuler_1step
        procedure:: ODEeuler_nstep
    end interface
    
!-------------------------------------------------------------------   
    interface ODERK
        procedure:: ODERK_1step
        procedure:: ODERK_nstep
    end interface
    
!-------------------------------------------------------------------
    abstract interface
        pure real(rp) function absdydx(x,y) result(dydx)
        import:: rp
        real(rp),intent(in)::   x,y 
        end function absdydx
    end interface

    
contains

    
    pure real(rp) function ODEeuler_1step(dydx,dx,x0,y0) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
        y = y0 + dydx(x0,y0)*dx
    end function ODEeuler_1step
    
    pure function ODEeuler_nstep(dydx,dx,x0,y0,n) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    integer(ip),intent(in)::    n
    real(rp),dimension(n)::     y
    integer(ip)::               i
        y(1) = y0
        do i=2,n
            y(i) = ODEeuler(dydx,dx,x0+(i-1)*dx,y(i-1))
        enddo
    end function ODEeuler_nstep
    
!------------------------------------------------------------------- 
contains


     pure real(rp) function ODERK_1step(dydx,dx,x0,y0) result(y)
     procedure(absdydx)::        dydx
     real(rp),intent(in)::       dx,x0,y0
     real(rp)::                  k1,k2,k3,k4
         k1 = dydx(x0,y0)
         k2 = dydx(x0 + 0.5d0*dx,y0 + 0.5d0*dx*k1)
         k3 = dydx(x0 + 0.5d0*dx,y0 + 0.5d0*dx*k2)
         k4 = dydx(x0 + dx,y0 + dx*k3)
         y = y0 + (1/6.d0)*dx*(k1 + 2*k2 + 2*k3 + k4)
     end function ODERK_1step

     pure real(rp) function ODERK_nstep(dydx,dx,x0,y0,n) result(y)
     procedure(absdydx)::        dydx
     real(rp),intent(in)::       dx,x0,y0
     integer(ip),intent(in)::    n
     real(rp),dimension(n)::     y
     integer(ip)::               i
         y(1)=y0
         do i=2,n
             y(i) = ODERK(dydx,dx,x0+(i-1)*dx,y(i-1))
         end do
     end function ODERK_nstep

end module ODElib
