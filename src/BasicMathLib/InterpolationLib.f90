!reconstruction method is mainly derived from conservation equation
!reconstruction is a special interpolation based on averaged information rather than point information
!refer to <High Order Weighted Essentially Nonoscillatory Schemes for Convection Dominated Problems>
module interpolationLib
use constants
use arrayOpsLib
use vector_
implicit none

    private
    !interpolation | LgrIntp = Lagrange interpolation
    public:: LgrIntp,cntLgrIntp,exLgrIntp
    !reconstruction |
    public:: exRctr3
    public:: muscl2, minmod, vanLeer, vanAlbada, superbee
    public:: muscl2c, muscl2c_sp, muscl2c_spstd
    public:: weno5, weno5z
    
    !----------------------------------------------------
    !enum of limiter |abondon doubleminmod due to its asymmetry
    integer(ip),parameter::         minmod = 347980
    integer(ip),parameter::         vanLeer = 4871730
    integer(ip),parameter::         vanAlbada = 14292870
    integer(ip),parameter::         superbee = 471004821
    
    
    !---------------------------------------------------
    interface LgrIntp
        procedure:: LgrIntp1d
    end interface LgrIntp
    
    interface exLgrIntp
        !(x--1--s1--1--s2)
        procedure:: exLgrIntp2_scalar
        procedure:: exLgrIntp2_vector
        procedure:: exLgrIntp2_array
        !(x--1--s1--1--s2--1--s3)
        procedure:: exLgrIntp3_scalar
        procedure:: exLgrIntp3_vector
        procedure:: exLgrIntp3_array
    end interface exLgrIntp
    
    interface cntLgrIntp
        !(l2--1--l1--1--x--1--r1--1--r2)
        procedure:: cntLgrIntp2_scalar
        procedure:: cntLgrIntp2_vector
        procedure:: cntLgrIntp2_array
    end interface cntLgrIntp
    
    
    !-----------------------
    !reconstruction
    interface muscl2
        !(p(1)--1--p(2)--1/2--x--1/2--p(3))
        procedure muscl2_scalar
        procedure muscl2_array
        procedure muscl2_arrayarray
    end interface muscl2
    
    interface muscl2c
        !(p(1)--1--p(2)--1/2--x--1/2--p(3))
        procedure muscl2c_scalar
        procedure muscl2c_array
        procedure muscl2c_arrayarray
    end interface muscl2c
    
    interface muscl2c_sp
        procedure muscl2c_SpOrnoCoef
        procedure muscl2c_SpOrnoCoefArray
    end interface muscl2c_sp
    
    interface muscl2c_spstd
        procedure muscl2c_SpOrnoStd
    end interface muscl2c_spstd
        
    !(p(1)--1--p(2)--1--p(3)--1/2--x--1/2--p(4)--1--p(5))
    interface weno5
        procedure:: weno5_scalar
        procedure:: weno5_array
        procedure:: weno5_arrayarray
    end interface weno5
    
    interface weno5z    
        procedure:: weno5z_scalar
        procedure:: weno5z_array
        procedure:: weno5z_arrayarray
    end interface weno5z
    
contains
    
    pure real(rp) function LgrIntp1d(x,vals,valsx) result(v)
    real(rp),intent(in)::               x
    real(rp),dimension(:),intent(in)::  vals,valsx
    integer(ip)::                       i,j,n
    real(rp)::                          cu,cl
        n = min(size(vals),size(valsx))
        v = zero
        do j=1,n
            cu = 1._rp
            cl = 1._rp
            do i=1,j-1
                cu = cu*(x-valsx(i))
                cl = cl*(valsx(j) - valsx(i))
            enddo
            do i=j+1,n
                cu = cu*(x-valsx(i))
                cl = cl*(valsx(j) - valsx(i))
            enddo
            v = v + vals(j)*cu/cl
        enddo
    end function LgrIntp1d

    
    !-----------------------------------------------------------
    pure real(rp) function exLgrIntp2_scalar(s1,s2) result(x)
    real(rp),intent(in):: s1,s2
        x = 2._rp*s1 - s2
    end function exLgrIntp2_scalar
    
    pure type(vector) function exLgrIntp2_vector(s1,s2) result(x)
    type(vector),intent(in):: s1,s2
        x = 2._rp*s1 - s2
    end function exLgrIntp2_vector
    
    pure function exLgrIntp2_array(s1,s2) result(x)
    real(rp),dimension(:),intent(in)::  s1,s2
    real(rp),dimension(size(s1))::      x
        x = 2._rp*s1 - s2
    end function exLgrIntp2_array
    
    pure real(rp) function exLgrIntp3_scalar(s1,s2,s3) result(x)
    real(rp),intent(in):: s1,s2,s3
        x = (11._rp*s1 - 7._rp*s2 + 2._rp*s3)/6._rp
    end function exLgrIntp3_scalar
    
    pure type(vector) function exLgrIntp3_vector(s1,s2,s3) result(x)
    type(vector),intent(in):: s1,s2,s3
        x = (11._rp*s1 - 7._rp*s2 + 2._rp*s3)/6._rp
    end function exLgrIntp3_vector
    
    pure function exLgrIntp3_array(s1,s2,s3) result(x)
    real(rp),dimension(:),intent(in)::  s1,s2,s3
    real(rp),dimension(size(s1))::      x
        x = (11._rp*s1 - 7._rp*s2 + 2._rp*s3)/6._rp
    end function exLgrIntp3_array
    
    !-------
    pure real(rp) function cntLgrIntp2_scalar(l2,l1,r1,r2) result(x)
    real(rp),intent(in):: l2,l1,r1,r2
        x = (4._rp*(l1 + r1) - l2 - r2)/6._rp
    end function cntLgrIntp2_scalar
    
    pure type(vector) function cntLgrIntp2_vector(l2,l1,r1,r2) result(x)
    type(vector),intent(in):: l2,l1,r1,r2
        x = (4._rp*( l1 + r1) - l2 - r2 )/6._rp
    end function cntLgrIntp2_vector
    
    pure function cntLgrIntp2_array(l2,l1,r1,r2) result(x)
    real(rp),dimension(:),intent(in)::  l2,l1,r1,r2
    real(rp),dimension(size(l2))::      x
        x = (4._rp*(l1 + r1) - l2 - r2)/6._rp
    end function cntLgrIntp2_array
    

    
    !refer to <Accurate, efficient and monotonic numerical methods for multi-dimensional 
    !compressible flows Part II: Multi-dimensional limiting process> [formula (6a) (7a)]
    !original expression
    !r(2) = delta(1)/(delta(2) + GlobalEps)
    !x = f(2) + 0.25_rp*((1._rp-kappa)*Phi(r(1),ilimiter) + &
    !                    (1._rp+kappa)*Phi(r(2),ilimiter)*r(1)) * delta(1)
    pure real(rp) function muscl2_scalar(f,ilimiter) result(x)
    real(rp),dimension(:),intent(in)::  f
    integer(ip),intent(in)::            ilimiter
    real(rp)::                          delta(2),r
    !real(rp),parameter::                kappa = 1._rp/3._rp
    
        delta(1) = f(2) - f(1)
        delta(2) = f(3) - f(2)
        if(delta(1)==zero) delta(1) = GlobalEps
        r = delta(2)/delta(1)
        !if limiter Phi is symmetric: Phi(r) = r*Phi(1/r), the orginal expression reduce to
        x = f(2) + 0.5_rp*Phi(r, ilimiter)*delta(1)

    contains
    
        pure real(rp) function Phi(r,i)
        real(rp),intent(in)::       r
        integer(ip),intent(in)::    i
            select case(i)
            !symmetric limiter
            case(minmod)
                Phi = max(zero, min(r, 1._rp))
            case(vanLeer)
                Phi = (r + abs(r))/(1._rp + abs(r))
            case(vanAlbada)
                Phi = (r**2 + r)/(1._rp + r**2)
            case(superbee)
                Phi = max(zero, min(2._rp*r, 1._rp), min(r, 2._rp))
            !asymmetric limiter
            !case(doubleMinmod)
            !    Phi = max(zero,min(2._rp*r, 1._rp, (1._rp+r)/2._rp))
            case default
                call disableProgram
            end select
        end function Phi
        
    end function muscl2_scalar
    
    !--
    pure function muscl2_array(f,ilimiter) result(x)
    real(rp),dimension(:,:),intent(in)::    f
    integer(ip),intent(in)::                ilimiter
    real(rp),dimension(size(f,1))::         x
    integer(ip)::                           i
        forall(i=1:size(f,1)) x(i) = muscl2(f(i,:), ilimiter)
    end function muscl2_array
    
    !--
    pure function muscl2_arrayarray(f,ilimiter) result(x)
    real(rp),dimension(:,:,:),intent(in)::      f
    integer(ip),intent(in)::                    ilimiter
    real(rp),dimension(size(f,1),size(f,2))::   x
    integer(ip)::                               i
        forall(i=1:size(f,1)) x(i,:) = muscl2(f(i,:,:), ilimiter)
    end function  muscl2_arrayarray
    
    
    !---------------------------------------------------------------
    !refer to <Semi-discrete central-upwind scheme for hyperbolic conservation laws 
    !           and Hamilton-Jacobi equations>
    !implementation below leads to nan in some situation, abondon!
    !leave here as a reference to compare with muscl2
    !delta(1) = (f(2) - f(1))
    !delta(2) = (f(3) - f(2))
    !delta(3) = (f(3) - f(1))/2._rp
    !r = delta(2)/(delta(1) + GlobalEps)
    !beta = delta(3)/(delta(1) + GlobalEps)
    !x = f(2) + 0.5_rp*max(zero, min(theta, theta*r, beta))*delta(1)
    pure real(rp) function muscl2c_scalar(f,theta) result(x)
    real(rp),dimension(:),intent(in)::  f
    real(rp),intent(in)::               theta
    real(rp),dimension(3)::             delta
        delta(1) = (f(2) - f(1))*theta
        delta(2) = (f(3) - f(2))*theta
        delta(3) = (f(3) - f(1))/2._rp
        if(delta(1)*delta(2)<0._rp) then
            x = f(2)
        else
            x = f(2) + 0.5_rp*delta(minloc(abs(delta), 1))
        endif
    end function  muscl2c_scalar
    !--
    pure function muscl2c_array(f,theta) result(x)
    real(rp),dimension(:,:),intent(in)::    f
    real(rp),intent(in)::                   theta
    real(rp),dimension(size(f,1))::         x
    integer(ip)::                           i
        forall(i=1:size(f,1)) x(i) = muscl2c(f(i,:), theta)
    end function  muscl2c_array
    !--
    pure function muscl2c_arrayarray(f,theta) result(x)
    real(rp),dimension(:,:,:),intent(in)::      f
    real(rp),intent(in)::                       theta
    real(rp),dimension(size(f,1),size(f,2))::   x
    integer(ip)::                               i
        forall(i=1:size(f,1)) x(i,:) = muscl2c(f(i,:,:), theta)
    end function  muscl2c_arrayarray
    
    !--Orno => orthonormal
    pure function muscl2c_SpOrnoCoef(f, theta, beta) result(x)
    real(rp),dimension(:,:),intent(in)::        f
    real(rp),intent(in)::                       theta, beta
    real(rp),dimension(size(f,1))::             x
    real(rp),dimension(size(f,1),3)::           delta
    integer(ip)::                               loc
    real(rp)::                                  nor,eps
        delta(:,1) = (f(:,2) - f(:,1))*theta
        delta(:,2) = (f(:,3) - f(:,2))*theta
        delta(:,3) = (f(:,3) - f(:,1))/2._rp
        !functional approah 1. check direction; 2. check variation
        nor = norm(delta(:,1))*norm(delta(:,2))
        eps = nor*globaleps
        if((delta(:,1).ip.delta(:,2)) + eps < nor*cos(beta)) then
            x(:) = f(:,2)
        else
            loc = minloc([norm2(delta(:,1)), norm2(delta(:,2)), norm2(delta(:,3))], 1)
            x(:) = f(:,2) + 0.5_rp*delta(:,loc)
        endif
        
    end function muscl2c_SpOrnoCoef
    !--
    pure function muscl2c_SpOrnoCoefArray(f,theta,beta) result(x)
    real(rp),dimension(:,:,:),intent(in)::      f
    real(rp),intent(in)::                       theta, beta
    real(rp),dimension(size(f,1),size(f,2))::   x
    integer(ip)::                               i
        forall(i=1:size(f,2)) x(:,i) = muscl2c_sp(f(:,i,:), theta, beta)
    end function  muscl2c_SpOrnoCoefArray
    
    !Orno => orthonormal
    pure function muscl2c_SpOrnoStd(f,theta) result(x)
    real(rp),dimension(:,:),intent(in)::        f
    real(rp),intent(in)::                       theta
    real(rp),dimension(size(f,1))::             x
    integer(ip)::                               n
        n = size(f,1)
        x(1) = muscl2c(f(1,:), theta)
        x(2:n) = muscl2c_sp(f(2:n,:), theta, hfpi)
    end function muscl2c_SpOrnoStd

    !-------------
    !reconstruction
    !e.g. for exRctr3
    !p(x) = a0 + a1*x + a2*x**2
    !int_x0^{x0+dx} p(x) = f1*dx
    !int_{x0+dx}^{x0+2*dx} p(x) = f2*dx
    !int_{x0+2*dx}^{x0+3*dx} p(x) = f3*dx
    !ASK for [a,b,c] to satisfy
    !p(x0+3*dx) = x = sum([a,b,c]*f)
    !that is [a = 1._rp/3._rp, b = -7._rp/6._rp, c = 11._rp/6._rp]
    !--------------------------------------------------
    !|-f1--f2--f3-x|
    pure real(rp) function exRctr3(f) result(x)
    real(rp),dimension(:),intent(in)::  f
    real(rp),parameter::    a = 1._rp/3._rp, b = -7._rp/6._rp, c = 11._rp/6._rp
        x = sum([a,b,c]*f)
    end function exRctr3
    
    !|-f1--f2-x-f3-|
    pure real(rp) function RsRctr3(f) result(x)
    real(rp),dimension(:),intent(in)::  f
    real(rp),parameter::    a = -1._rp/6._rp, b = 5._rp/6._rp, c = 1._rp/3._rp
        x = sum([a,b,c]*f)
    end function RsRctr3
    
    !|-f1-x-f2--f3-|
    pure real(rp) function LsRctr3(f) result(x)
    real(rp),dimension(:),intent(in)::  f
        x = RsRctr3(f(3:1:-1))
    end function LsRctr3
    
    !---------
    !refer to <High Order Weighted Essentially Nonoscillatory Schemes for Convection Dominated Problems>
    pure function weno5_scalar(f) result(x)
    real(rp),dimension(:),intent(in)::f
    real(rp)::              x,s0,s1,s2,w0,w1,w2,w
    real(rp),parameter::    eps = 1.e-40_rp
    integer(ip),parameter:: p = 2
    real(rp),parameter::    alpha1 = 0.25_rp, alpha2 = 13._rp/12._rp
    
        !nonlinear indicator of smoothness
        s0 = alpha1*(f(1) - 4._rp*f(2) + 3._rp*f(3))**2 &
           + alpha2*(f(1) - 2._rp*f(2) + f(3))**2
                                
        s1 = alpha1*(f(2) - f(4))**2 &
           + alpha2*(f(2) - 2._rp*f(3) + f(4))**2

        s2 = alpha1*(3._rp*f(3) - 4._rp*f(4) + f(5))**2 &
           + alpha2*(f(3) - 2._rp*f(4) + f(5))**2

        !the coefficients for the 5th order Rctr from 3th order Rctr
        w0 = 0.1_rp/(s0+eps)**p
        w1 = 0.6_rp/(s1+eps)**p
        w2 = 0.3_rp/(s2+eps)**p
        
        w = w0 + w1 + w2
        w0 = w0/w
        w1 = w1/w
        w2 = w2/w

        x = w0*exRctr3(f(1:3)) + w1*RsRctr3(f(2:4)) + w2*LsRctr3(f(3:5))

    end function weno5_scalar
    
    !--dim(f,1) is field dimension and dim(f,2) is the node index
    pure function weno5_array(f) result(x)
    real(rp),dimension(:,:),intent(in)::f
    real(rp),dimension(size(f,1))::     x
    integer(ip)::                       i
        forall(i=1:size(f,1)) x(i) = weno5(f(i,:))
    end function weno5_array
    
    !--dim(f,1&2) is field dimension and dim(f,3) is the node index
    pure function weno5_arrayarray(f) result(x)
    real(rp),dimension(:,:,:),intent(in)::      f
    real(rp),dimension(size(f,1),size(f,2))::   x
    integer(ip)::                               i
        forall(i=1:size(f,1)) x(i,:) = weno5(f(i,:,:))
    end function weno5_arrayarray
    
    !refer to <An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws>
    !refer to <jcp - An improved WENO-Z scheme>
    pure function weno5z_scalar(f) result(x)
    real(rp),dimension(:),intent(in)::	f
    real(rp)::							x,s0,s1,s2,w0,w1,w2,w
    real(rp),parameter::				eps = 1.e-20_rp
    integer(ip),parameter::				p = 2
    real(rp),parameter::				alpha1 = 0.25_rp, alpha2 = 13._rp/12._rp
    
        !nonlinear indicator of smoothness
        s0 = alpha1*(f(1) - 4._rp*f(2) + 3._rp*f(3))**2 &
           + alpha2*(f(1) - 2._rp*f(2) + f(3))**2
                                
        s1 = alpha1*(f(2) - f(4))**2 &
           + alpha2*(f(2) - 2._rp*f(3) + f(4))**2

        s2 = alpha1*(3._rp*f(3) - 4._rp*f(4) + f(5))**2 &
           + alpha2*(f(3) - 2._rp*f(4) + f(5))**2

        !the coefficients for the 5th order Rctr from 3th order Rctr
        w0 = 0.1_rp*(1._rp + ((s0-s2)/(s0+eps))**p)
        w1 = 0.6_rp*(1._rp + ((s0-s2)/(s1+eps))**p)
        w2 = 0.3_rp*(1._rp + ((s0-s2)/(s2+eps))**p)
        
        w = w0 + w1 + w2
        w0 = w0/w
        w1 = w1/w
        w2 = w2/w

        x = w0*exRctr3(f(1:3)) + w1*RsRctr3(f(2:4)) + w2*LsRctr3(f(3:5))

    end function weno5z_scalar
    
    !--dim(f,1) is field dimension and dim(f,2) is the node index
    pure function weno5z_array(f) result(x)
    real(rp),dimension(:,:),intent(in)::f
    real(rp),dimension(size(f,1))::     x
    integer(ip)::                       i
        forall(i=1:size(f,1)) x(i) = weno5z(f(i,:))
    end function weno5z_array
    
    !--dim(f,1&2) is field dimension and dim(f,3) is the node index
    pure function weno5z_arrayarray(f) result(x)
    real(rp),dimension(:,:,:),intent(in)::      f
    real(rp),dimension(size(f,1),size(f,2))::   x
    integer(ip)::                               i
        forall(i=1:size(f,1)) x(i,:) = weno5z(f(i,:,:))
    end function weno5z_arrayarray

end module interpolationLib