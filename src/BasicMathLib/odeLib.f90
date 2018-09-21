!i want to write some basic operations for ordinary differetial equations
!as the basic lib for Type<timedriver> or others
    

!Here the odelib contains Euler solver and Runge-Kutta solver that uses the fourth order TVD format.
module odelib
use constants
implicit none

    
    private
    public::    odeEuler
    public::    odeRK4,odesRK4
    public::    odeRK2
    
    !--  
    interface odeEuler
        procedure:: odeEuler_1
        procedure:: odeEuler_n
    end interface odeEuler
    
    !--
    interface odeRK4
        procedure:: odeRK4_TVD_1
        procedure:: odeRK4_TVD_n
    end interface odeRK4
    
    !--
    interface odesRK4
        procedure:: odesRK4_g_TVD_1
    end interface odesRK4

    !--
    interface odeRK2
        procedure:: odeRK2_TVD_1
        procedure:: odeRK2_TVD_n
    end interface odeRK2
    
    
    
    
    !-------------------------------------------------------------------
    abstract interface
        pure real(rp) function absdydx(x,y) result(dydx)
        import:: rp
        real(rp),intent(in)::               x,y 
        end function absdydx
        
        pure real(rp) function absdydxSys(i,x,y) result(dydx)
        import:: rp,ip
        integer(ip),intent(in)::            i
        real(rp),intent(in)::               x,y
        end function absdydxSys
        
        pure real(rp) function absdydxSysg(i,x,y) result(dydx)
        import:: rp,ip
        integer(ip),intent(in)::            i
        real(rp),dimension(:),intent(in)::  x,y
        end function absdydxSysg
    end interface

    
    !------------------------------------------------------------------
    
contains

    !------------------------------------------------------------------
    pure real(rp) function odeEuler_1(dydx,dx,x0,y0) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
        y = y0 + dydx(x0,y0)*dx
    end function odeEuler_1
    
    pure function odeEuler_n(dydx,dx,x0,y0,n) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    integer(ip),intent(in)::    n
    real(rp),dimension(n)::     y
    integer(ip)::               i
        y(1)=odeEuler(dydx,dx,x0,y0)
        do i=2,n
            y(i) = odeEuler(dydx, dx, x0+(i-1)*dx, y(i-1))
        end do
    end function odeEuler_n
    
    !------------------------------------------------------------------- 
    pure real(rp) function odeRK4_TVD_1(dydx,dx,x0,y0) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    real(rp)::                  k1,k2,k3,k4
        k1 = dydx(x0,y0)
        k2 = dydx(x0 + 0.5_rp*dx, y0 + 0.5_rp*dx*k1)
        k3 = dydx(x0 + 0.5_rp*dx, y0 + 0.5_rp*dx*k2)
        k4 = dydx(x0 + dx, y0 + dx*k3)
        y = y0 + (1._rp/6._rp)*dx*(k1 + 2._rp*k2 + 2._rp*k3 + k4)
    end function odeRK4_TVD_1

    !--
    pure function odeRK4_TVD_n(dydx,dx,x0,y0,n) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    integer(ip),intent(in)::    n
    real(rp),dimension(n)::     y
    integer(ip)::               i
        y(1)=odeRK4(dydx,dx,x0,y0)
        do i=2,n
            y(i) = odeRK4(dydx, dx, x0+(i-1)*dx , y(i-1))
        end do
    end function odeRK4_TVD_n
    
    !--
    pure subroutine odesRK4_TVD_1(dydx,dx,x0,y0,y)
    procedure(absdydxSys)::             dydx
    real(rp),intent(in)::               dx
    real(rp),dimension(:),intent(in)::  x0,y0
    real(rp),dimension(:),intent(out):: y
    real(rp)::                          k1,k2,k3,k4
    integer(ip)::                       i
        do i=1,size(x0)
            y(i) = odeRK4(d, dx, x0(i), y0(i))
        enddo
    contains
        pure real(rp) function d(x,y)
        real(rp),intent(in)::   x,y
            d = dydx(i,x,y)
        end function d
    end subroutine odesRK4_TVD_1
    
    !--
    pure subroutine odesRK4_g_TVD_1(dydx,dx,x0,y0,y)
    procedure(absdydxSysg)::            dydx
    real(rp),intent(in)::               dx
    real(rp),dimension(:),intent(in)::  x0,y0
    real(rp),dimension(:),intent(out):: y
    real(rp)::                          k1,k2,k3,k4
    integer(ip)::                       i
        do i=1,size(x0)
            y(i) = odeRK4(d, dx, x0(i), y0(i))
        enddo
    contains
        pure real(rp) function d(x,y)
        real(rp),intent(in)::   x,y
            d = dydx(i,x0,y0)
        end function d
    end subroutine odesRK4_g_TVD_1
    

!---------------------------------------------------------------------
    pure real(rp) function odeRK2_TVD_1(dydx,dx,x0,y0) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    real(rp)::                  k1,k2
        k1 = dydx(x0, y0)
        k2 = dydx(x0+dx, y0+dx*k1)
        y = y0 + 0.5_rp*dx*(k1 + k2)
    end function odeRK2_TVD_1

    pure function odeRK2_TVD_n(dydx,dx,x0,y0,n) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    integer(ip),intent(in)::    n
    real(rp),dimension(n)::     y
    integer(ip)::               i
        y(1)=odeRK2(dydx,dx,x0,y0)
        do i=2,n
            y(i) = odeRK2(dydx, dx, x0+(i-1)*dx, y(i-1))
        end do
    end function odeRK2_TVD_n


end module odelib