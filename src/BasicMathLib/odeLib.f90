!i want to write some basic operations for ordinary differetial equations
!as the basic lib for Type<timedriver> or others
    

!Here the odelib contains Euler solver and Runge-Kutta solver that uses the fourth order Tvd format.
module odelib
use constants
implicit none

    
    private
    public::    odeEuler
    public::    odeRk4,odesRk4
    public::    odeRk2
    public::    odeRk3
    
    !--  
    interface odeEuler
        procedure:: odeEuler_1
        procedure:: odeEuler_n
    end interface odeEuler
    
    !--
    interface odeRk4
        procedure:: odeRk4Tvd_1
        procedure:: odeRk4Tvd_n
    end interface odeRk4
    
    !--
    interface odesRk4
        procedure:: odesRk4_gTvd_1
    end interface odesRk4

    !--
    interface odeRk2
        procedure:: odeRk2Tvd_1
        procedure:: odeRk2Tvd_n
    end interface odeRk2
    
    !--
    interface odeRk3
        procedure:: odeRk3Tvd_1
        procedure:: odeRk3Tvd_2
    end interface odeRk3
    
    
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
    pure real(rp) function odeRk4Tvd_1(dydx,dx,x0,y0) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    real(rp)::                  k1,k2,k3,k4
        k1 = dydx(x0,y0)
        k2 = dydx(x0 + 0.5_rp*dx, y0 + 0.5_rp*dx*k1)
        k3 = dydx(x0 + 0.5_rp*dx, y0 + 0.5_rp*dx*k2)
        k4 = dydx(x0 + dx, y0 + dx*k3)
        y = y0 + (1._rp/6._rp)*dx*(k1 + 2._rp*k2 + 2._rp*k3 + k4)
    end function odeRk4Tvd_1

    !--
    pure function odeRk4Tvd_n(dydx,dx,x0,y0,n) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    integer(ip),intent(in)::    n
    real(rp),dimension(n)::     y
    integer(ip)::               i
        y(1)=odeRk4(dydx,dx,x0,y0)
        do i=2,n
            y(i) = odeRk4(dydx, dx, x0+(i-1)*dx , y(i-1))
        end do
    end function odeRk4Tvd_n
    
    !--
    pure subroutine odesRk4Tvd_1(dydx,dx,x0,y0,y)
    procedure(absdydxSys)::             dydx
    real(rp),intent(in)::               dx
    real(rp),dimension(:),intent(in)::  x0,y0
    real(rp),dimension(:),intent(out):: y
    real(rp)::                          k1,k2,k3,k4
    integer(ip)::                       i
        do i=1,size(x0)
            y(i) = odeRk4(d, dx, x0(i), y0(i))
        enddo
    contains
        pure real(rp) function d(x,y)
        real(rp),intent(in)::   x,y
            d = dydx(i,x,y)
        end function d
    end subroutine odesRk4Tvd_1
    
    !--
    pure subroutine odesRk4_gTvd_1(dydx,dx,x0,y0,y)
    procedure(absdydxSysg)::            dydx
    real(rp),intent(in)::               dx
    real(rp),dimension(:),intent(in)::  x0,y0
    real(rp),dimension(:),intent(out):: y
    real(rp)::                          k1,k2,k3,k4
    integer(ip)::                       i
        do i=1,size(x0)
            y(i) = odeRk4(d, dx, x0(i), y0(i))
        enddo
    contains
        pure real(rp) function d(x,y)
        real(rp),intent(in)::   x,y
            d = dydx(i,x0,y0)
        end function d
    end subroutine odesRk4_gTvd_1
    

!---------------------------------------------------------------------
    pure real(rp) function odeRk2Tvd_1(dydx,dx,x0,y0) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    real(rp)::                  k1,k2
        k1 = dydx(x0, y0)
        k2 = dydx(x0+dx, y0+dx*k1)
        y = y0 + 0.5_rp*dx*(k1 + k2)
    end function odeRk2Tvd_1

    pure function odeRk2Tvd_n(dydx,dx,x0,y0,n) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    integer(ip),intent(in)::    n
    real(rp),dimension(n)::     y
    integer(ip)::               i
        y(1)=odeRk2(dydx,dx,x0,y0)
        do i=2,n
            y(i) = odeRk2(dydx, dx, x0+(i-1)*dx, y(i-1))
        end do
    end function odeRk2Tvd_n
    
!---------------------------------------------------------------------
    !<Efficient implementation of weighted ENO schemes>    
    !same as 3order Rk in time driver
    !dydx = rhs(y)
    subroutine odeRk3Tvd_1(dydx, dx, y0, y, rhs)
    real(rp),intent(in)::                   dx
    real(rp),dimension(:),intent(in)::      y0
    real(rp),dimension(:),intent(out)::     y, rhs
    abstract interface
        subroutine dydx(y, rhs)
        import rp
        real(rp),dimension(:),intent(in)::  y
        real(rp),dimension(:),intent(out):: rhs
        end subroutine dydx
    end interface
    
        call dydx(y0, rhs)
        y = y0 + dx*rhs
        call dydx(y, rhs)
        y = 0.25_rp*(3._rp*y0 + y + dx*rhs)
        call dydx(y, rhs)
        y = 1._rp/3._rp*(y0 + 2._rp*y + 2._rp*dx*rhs)
    
    end subroutine odeRk3Tvd_1
    
    !--
    subroutine odeRk3Tvd_2(dydx, dx, y0, y, rhs)
    real(rp),intent(in)::                   dx
    real(rp),dimension(:,:),intent(in)::    y0
    real(rp),dimension(:,:),intent(out)::   y, rhs
    abstract interface
        subroutine dydx(y, rhs)
        import rp
        real(rp),dimension(:,:),intent(in)::    y
        real(rp),dimension(:,:),intent(out)::   rhs
        end subroutine dydx
    end interface
    
        call dydx(y0, rhs)
        y = y0 + dx*rhs
        call dydx(y, rhs)
        y = 0.25_rp*(3._rp*y0 + y + dx*rhs)
        call dydx(y, rhs)
        y = 1._rp/3._rp*(y0 + 2._rp*y + 2._rp*dx*rhs)
    
    end subroutine odeRk3Tvd_2

end module odelib