!numerical integeration is also denoted as cubature
!quadrature is more applied to one dimensional integral
!and cubature is referred to over more than one dimension
module IntegrationLib
use constants
use arrayOpsLib
use stringOpsLib
use SpecialFunctionLib
use fftWrapperLib
use laWrapperLib
use odelib
implicit none

    private
    
    public:: integrate
    !--
    public:: GaussCev   !based on the mkl symtridiag, more memory, more accurate
    public:: GaussTom   !based on the ACM TOM, less memory, less accurate
    !--
    public:: QuadratureRule
    public:: GaussLegendre
    public:: GaussHermite
    public:: ClenshawCurtis
    public:: ClenshawCurtisNative
    !--
    public:: SparseGridSize
    public:: SparseGrid
    !--
    public:: integrateTrapezoid_f1
    public:: integrateAdaptiveSimpson_f1
    
    
!------------------------------------------------
    interface integrate
        procedure:: integrateAdaptiveSimpson_f1
        procedure:: integrateQuadratureRule
    end interface integrate
    !--

!--------------------------------------------------
    interface GaussHermite
        procedure:: GaussHermitePhys
    end interface
    
contains
    
    !-------------------------------------------------------------------------
    pure function integrateTrapezoid_f1(func,lo,up,epsilon) result(integral)
    procedure(absf1)::                      func
    real(rp),intent(in)::                   lo,up
    real(rp),intent(in),optional::          epsilon
    integer(ip)::                           n
    real(rp)::                              integral0,integral,epsrun,eps
        
        epsrun  = 1.d100
        n       = 101
        eps     = merge(epsilon,GlobalEps,present(epsilon))
        
        integral = trapInt()
        do while(epsrun>max(integral*eps,eps).and.n<1e9)
            integral0   = integral
            n           = n*1.618_rp      !golden
            integral    = trapInt()
            epsrun      = abs(integral-integral0)
        enddo
        
    contains
    
        pure function trapInt()
        real(rp)::      trapInt,x,y0,y1,dx
        integer(ip)::   i
        
            trapInt = zero
            dx  = (up-lo)/dfloat(n)
            x   = lo
            y1  = func(lo)
            
            do i=1,n
                x   = lo + i*dx
                y0  = y1
                y1  = func(x)
                trapInt = trapInt + 0.5_rp*(y0+y1)*dx
            enddo
            
        end function trapInt
        
    end function integrateTrapezoid_f1
    
    
    !---------refer to wiki(adaptive simpson's method) and 
    !https: //www.zhihu.com/question/30141678
    pure real(rp) function integrateAdaptiveSimpson_f1(f,lo,up,eps,maxLoop) result(integral)
    procedure(absf1)::                  f
    real(rp),intent(in)::               lo,up
    real(rp),intent(in),optional::      eps
    integer(ip),intent(in),optional::   maxLoop
        integral = integrateAdaptiveSimpson_f1_rcs(f,lo,f(lo),up,f(up), &
                                    simpson_f1(f,lo,f(lo),up,f(up)),    &
                                    merge(eps,Globaleps,present(eps)),  &
                                    merge(maxLoop,15_ip,present(maxLoop)))
    end function integrateAdaptiveSimpson_f1
    
    !--read flo and fup can reduce the computing, notice the <fmiddle>
    pure recursive real(rp) function integrateAdaptiveSimpson_f1_rcs(f,lo,flo,up,fup,s,eps,max) result(integral)
    procedure(absf1)::          f
    real(rp),intent(in)::       lo,flo,up,fup,s,eps
    integer(ip),intent(in)::    max
    real(rp)::                  mi,fmi,l,r
        mi  = (lo+up)/2._rp
        fmi = f(mi)
        l   = simpson_f1(f,lo,flo,mi,fmi)
        r   = simpson_f1(f,mi,fmi,up,fup)
        if(max==0.or.abs(l+r-s)<=15._rp*eps) then   !15._rp comes from error analysis
            integral = l+r+(l+r-s)/15._rp
        else
            integral = integrateAdaptiveSimpson_f1_rcs(f,lo,flo,mi,fmi,l,eps/2._rp,max-1)+ &
                        integrateAdaptiveSimpson_f1_rcs(f,mi,fmi,up,fup,r,eps/2._rp,max-1)
        endif
    end function integrateAdaptiveSimpson_f1_rcs

    !--
    pure function simpson_f1(f,lo,flo,up,fup) result(r)
    procedure(absf1)::          f
    real(rp),intent(in)::       lo,up,flo,fup
    real(rp)::                  r,mi
        mi = lo + (up-lo)/2._rp
        r = (flo + 4._rp*f(mi) + fup)*(up-lo)/6._rp
    end function simpson_f1
    
    
!---------------------------------------------------------------------------------
!quadrature rule related
!---------------------------------------------------------------------------------
    !refer to <Remark on algorithm 726: ORTHPOL---A package of routines for 
    !           generating orthogonal polynomials and Gauss-type quadrature rules>
    !P_{n+1} = (x-alpha_n)*P_n - beta_n*P_{n-1}
    pure subroutine GaussCev(alpha,beta,quadx,quadw)
    real(rp),dimension(:),intent(in)::  alpha,beta
    real(rp),dimension(:),intent(out):: quadx,quadw
    integer(ip)::                       n
    real(rp),dimension(size(quadx),size(quadx))::ev   
        
        n = size(quadx)
        quadx = alpha
        quadw(1:n-1) = sqrt(beta(2:n))  !temp for edge
        call eigenSymTriDiagonal(quadx,quadw,ev)
        quadw = beta(1)*ev(1,:)**2
    
    end subroutine GaussCev
    
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/toms726/toms726.html
    pure subroutine GaussTom(alpha,beta,quadx,quadw)
    real(rp),dimension(:),intent(in)::  alpha,beta
    real(rp),dimension(:),intent(out):: quadx,quadw
    integer(ip)::                       i,j,k,l,m,m2,mml
    real(rp)::                          b,c,f,g,p,r,s
    real(rp),dimension(size(quadx))::   e
    real(rp),parameter::                eps = 10._rp*GlobalEps
    integer(ip)::                       n
    
        n = size(quadx)
        
        quadx(1) = alpha(1)
        quadw(1) = beta(1)
        if(n==1) return

        quadw(1) = 1._rp
        e(n) = 0._rp
        !diag = alpha(1:n), edge = sqrt(beta(2:n))
        do k = 2, n
            quadx(k) = alpha(k)
            e(k-1) = sqrt(beta(k))
            quadw(k) = 0._rp
        end do

        !QR factor
        do l = 1, n
            j = 0
            do
                do m2 = l, n
                    m = m2
                    if(m==n) exit
                    if(abs(e(m)) <= eps*(abs(quadx(m)) + abs(quadx(m+1)))) exit
                end do

                p = quadx(l)
                if(m==l) exit
                if(30 <= j) call disableprogram
                j = j + 1
                g = (quadx(l+1) - p)/(2._rp * e(l))
                r = sqrt(g*g + 1._rp)
                g = quadx(m) - p + e(l)/(g + sign(r,g))
                s = 1._rp
                c = 1._rp
                p = 0._rp
                mml = m - l

                do i = m-1, 1, -1
                    f = s * e(i)
                    b = c * e(i)
                    if (abs(f) < abs(g)) then
                        s = f / g
                        r = sqrt(s*s + 1._rp)
                        e(i+1) = g * r
                        c = 1._rp / r
                        s = s * c
                    else
                        c = g / f
                        r = sqrt(c*c + 1._rp)
                        e(i+1) = f * r
                        s = 1._rp / r
                        c = c * s
                    end if
                    g = quadx(i+1) - p
                    r = (quadx(i) - g)*s + 2._rp*c*b
                    p = s * r
                    quadx(i+1) = g + p
                    g = c * r - b
                    f = quadw(i+1)
                    quadw(i+1) = s * quadw(i) + c * f
                    quadw(i) = c * quadw(i) - s * f
                end do
                
                quadx(l) = quadx(l) - p
                e(l) = g
                e(m) = 0._rp
            end do
        end do

        !sorting
        do i=1,n-1
            k = minloc(quadx(i:n),1)
            if(k/=1) then
                call swap(quadx(i),quadx(k+i-1))
                call swap(quadw(i),quadw(k+i-1))
            endif
        enddo
        quadw(1:n) = beta(1)*quadw(1:n)**2
    
    end subroutine GaussTom
    
    
    !-------------------------------------------------
    real(rp) function integrateQuadratureRule(f,ruleType,Order,normalCoef) result(r)
    procedure(absf1)::              f
    character(*),intent(in)::       ruleType
    integer(ip),intent(in)::        Order
    real(rp),optional,intent(in)::  normalCoef
    real(rp),dimension(ishft(Order,-1)+1):: x,w
    
        call QuadratureRule(ruleType,x,w)
        r = merge(normalCoef, 1._rp, present(normalCoef))*sum(f(x)*w)
        
    end function integrateQuadratureRule
    
    
    !---------------------------------------------------
    subroutine quadratureRule(rule,quadx,quadw)
    character(*),intent(in)::               rule
    real(rp),dimension(:),intent(out)::     quadx,quadw
    character(len(rule))::                  r
    
        r = rule;   call lowerstring(r)
        
        select case(adjustl(r))
        case('gausslegendre','gl','legendre')
            call GaussLegendre(quadx,quadw)
        case('gausshermite','gh','hermite')
            call GaussHermite(quadx,quadw)
        case('gausshermite_prob','gh_prob','hermite_prob')
            call GaussHermite(quadx,quadw,'prob')
        case('clenshawcurtis','cc')
            call ClenshawCurtis(quadx,quadw)
        case default
            stop 'error: IntegrationLib/quadratureRule get an unrecoganized rule type'
        end select
        
    end subroutine quadratureRule
    
    
    !refer to https://github.com/chebfun/chebfun/blob/34f92d12ca51003f5c0033bbeb4ff57ac9c84e78/legpts.m
    !maybe a better choice https://github.com/Pazus/Legendre-Gauss-Quadrature/blob/master/legzo.m
    !calculate \[\int_{-1}^{1} f(x) dx\],(refer to wiki, \pi(x)=1 rather than \pi(x)=\frac{1}{2})
    !asymptotic method unavailiable temporarily
    !n points reach 2n-1 order accuracy | N+1 points reach 2N+1 order accuracy
    pure subroutine GaussLegendre(quadx,quadw)
    real(rp),dimension(:),intent(out)::         quadx,quadw
    integer(ip)::                               n
    real(rp)::                                  a,b,c,d
        !--
        n = size(quadx)
        !--
        if(n<6) then
        
            select case(n)
            case(1)
                quadx = 0._rp
                quadw = 2._rp
            case(2)
                c = 1._rp/sqrt(3._rp)
                !--
                quadx = [-c , c]
                quadw = [1._rp , 1._rp]
            case(3)
                c = sqrt(3._rp/5._rp)
                !--
                quadx = [-c , 0._rp , c]
                quadw = [5._rp/9._rp , 8._rp/9._rp , 5._rp/9._rp]
            case(4)
                a = 2._rp/7._rp * sqrt(6._rp/5._rp)
                b = sqrt(3._rp/7._rp - a)
                c = sqrt(3._rp/7._rp + a)
                d = sqrt(30._rp)/36._rp
                !--
                quadx = [-c , -b , b , c]            
                quadw = [0.5_rp-d , 0.5_rp+d , 0.5_rp+d , 0.5_rp-d]
            case(5)
                a = 2._rp * sqrt(10._rp/7._rp)
                b = 1._rp/3._rp * sqrt(5._rp - a)
                c = 1._rp/3._rp * sqrt(5._rp + a)
                d = 13._rp * sqrt(70._rp)
                !--
                quadx = [-c , -b , 0._rp , b , c]
                quadw = [ (322._rp-d)/900._rp , (322._rp+d)/900._rp , 128._rp/225._rp,  &
                            (322._rp+d)/900._rp , (322._rp-d)/900._rp ]
            end select
            
        elseif(n<1e6) then  !60,80,100: all ok. temperarialy relax the restriction
        
            call iterativeMethod(n)
            
        else
        
            call asymptoticMethod(n)
            
        endif
    
    contains
    
        !use newton's method to find the roots of Legendre Polynomials
        pure subroutine iterativeMethod(n)
        integer(ip),intent(in)::            n
        integer(ip)::                       k,iter
        real(rp),dimension(ishft(n+1,-1)):: x0,dx,pm1,pm2,ppm1,ppm2,p,pp

            !Asymptotic formula (Tricomi), only for positive x, as the initialization
            x0 = dfloat( [ ishft(n+1,-1) : 1 : -1 ] )
            x0 = pi * (4._rp*x0 - 1._rp) / (4._rp*n + 2._rp)
            x0 = (1._rp - (n-1._rp)/(8._rp*n**3) - &
                1._rp/(384._rp*n**4)*(39._rp-28._rp/sin(x0)**2) ) * cos(x0)
                
            !--
            !initialize
            pm1 = x0
            pm2 = 1._rp
            ppm1 = 1._rp
            ppm2 = 0._rp
            dx = 100._rp
            iter = 0
            !---
            do while(norm(dx,-1)>GlobalEps.and.iter<10)
                iter = iter + 1
                !
                do k=1,n-1
                    p = ((2._rp*k + 1._rp)*pm1*x0 - k*pm2)/(k+1._rp)
                    pm2 = pm1
                    pm1 = p
                    pp = ((2._rp*k + 1._rp)*(Pm2 + x0*PPm1) - k*PPm2 ) / (k+1._rp)
                    ppm2 = ppm1
                    ppm1 = pp
                enddo
            
                dx = - p/pp
                x0 = x0 + dx
                
                !--
                pm1 = x0
                pm2 = 1._rp
                ppm1 = 1._rp
                ppm2 = 0._rp
            enddo
            
            do k=1,n-1
                p = ((2._rp*k + 1._rp)*pm1*x0 - k*pm2)/(k+1._rp)
                pm2 = pm1
                pm1 = p
                pp = ((2._rp*k + 1._rp)*(Pm2 + x0*PPm1) - k*PPm2 ) / (k+1._rp)
                ppm2 = ppm1
                ppm1 = pp
            enddo
            
            quadx = [-x0(ishft(n+1,-1):1+ibits(n,0,1):-1),x0]
            !deriavtives
            quadw = [ pp(ishft(n+1,-1):1+ibits(n,0,1):-1),pp]
            !--refer to wiki(GausssQuadrature)
            quadw = 2._rp /( (1._rp-quadx**2) * quadw**2 )
            
        end subroutine iterativeMethod
        
        !--
        pure subroutine asymptoticMethod(n)
        integer(ip),intent(in)::        n
        integer(ip)::                   m
        real(rp),dimension(ishft(n+1,-1)):: jk
            
            !give up, too many, waiting for some day
            call disableProgram
            
            m = ishft(n+1,-1)
            jk = besseljn_roots(0,m)
            
        
        end subroutine asymptoticMethod
        
    end subroutine GaussLegendre
    
    !------------------------------------------------------------------------------------
    !refer to https://github.com/chebfun/chebfun/blob/development/hermpts.m
    !calculate \[\int_{-inf}^{inf} f(x) exp(-x^2) dx\],(refer to wiki)
    !hermite asymptotic method unavailiable temporarily
    !GaussHermite returns hermite points and gauss-hermite quadrature weights.
    !by default these are roots of the 'physicist'-type hermite polynomials, which are orthogonal with respect to the weight exp(-x.^2).\
    !if(var=='prob'), it calculate \[\int_{-inf}^{inf} f(x) exp(-x^2/2) dx\]
    !------------------------------------------------------------------------------------
    pure subroutine GaussHermitePhys(quadx,quadw,var)
    real(rp),dimension(:),intent(out)::     quadx,quadw
    character(*),intent(in),optional::      var
    integer(ip)::                           n
    
        n=size(quadx)
        if(n==0) then
            quadx = 0._rp
            quadw = 0._rp
        elseif(n==1) then
            quadx = 0._rp 
            quadw = spi 
        else if(n<21) then
            !call GolubWelsch_eigenvalueMethod(n)
            call GlaserLiuRokhlinAlgorithm(n)
        else if(n<200) then
            call hermpts_rec(n)
        else
            !normally,we won't use this subroutine which is not accomplished
            call hermpts_asy(n)
        end if
        
        !Normalise so that sum(w) = spi
        quadw = (spi/sum(quadw))*quadw
        
        !--for probability weight exp(-x^2/2)
        if(present(var)) then
            if(var=='prob') then
                quadx = quadx * sqrt(2._rp)
                quadw = quadw * sqrt(2._rp)
            else
                call disableprogram
            endif
        endif
            
    contains
        
        !-------------------------------------------------------
        pure subroutine GlaserLiuRokhlinAlgorithm(n)
        integer(ip),intent(in)::            n
        real(rp),dimension(n)::             ders
        
            !'GLR' uses Glaser-Liu-Rokhlin fast algorithm which is much faster for large N
            call alg0_Herm(n,quadx,ders)
            quadw = (2._rp*exp(-quadx**2)/ders**2)  ! Quadrature weights
            
        end subroutine GlaserLiuRokhlinAlgorithm
        
        !-------------------------------------------------------
        pure subroutine alg0_Herm(n,x,ders)
        integer(ip),intent(in)::            n
        real(rp),dimension(:),intent(out):: x,ders
        
            ! Compute coefficients of H_m(0), H_m'(0), m = 0,..,N
            call first_root(n,x,ders)       !find first root
            call alg1_Herm(n,x,ders)        !compute roots and derivatives
        
        end subroutine alg0_Herm
        
        !-------------------------------------------------------
        pure subroutine alg1_Herm(n,roots,ders)
        integer(ip),intent(in)::            n
        real(rp),dimension(:),intent(out):: roots,ders
        integer(ip)::                       s,N1,j,l,i,k
        integer(ip),parameter::             m = 30  !number of terms in Taylor expansion 
        real(rp)::                          x,h,m1,c1,c2,c3,step,x0,x1,x11,y0,dx 
        real(rp),dimension(m)::             z,z1,z2
        real(rp),dimension(10)::            x2
        real(rp),dimension(m+1)::           hh1,u,up,hh
        
            s=ibits(n,0,1)
            N1=ishft(n,-1)
            !initialise
            hh1= 1._rp
            u  = 0._rp
            up = 0._rp
        
            do j=N1+1,n-1
                !previous root
                x = roots(j)
                !initial approx
                x0=  pi/2._rp
                x11=-pi/2._rp
                y0= x
                dx=(x11-x0)/10._rp
                x2=odeRK2(dydx,dx,x0,y0,10)
                x1=x2(10)
                h=x1-x
                !scaling
                m1 = 1._rp/h            
                !recurrence relation for Hermite polynomials
                c1 = -(2._rp*n+1-x**2)/(m1**2) 
                c2 = 2._rp*x/(m1**3)
                c3 = 1._rp/(m1**4)
                u(1) = 0._rp
                u(2) = ders(j)/m1 
                u(3) = 0.5_rp*c1*u(1)
                u(4) = (c1*u(2)+c2*u(1))/6._rp
                up(1) = u(2)
                up(2) = 2._rp*u(3)*m1
                up(3) = 3._rp*u(4)*m1
                up(m+1) = 0._rp
                do k = 2,m-2
                    u(k+3) = (c1*u(k+1) + c2*u(k) + c3*u(k-1))/((k+1)*(k+2))
                    up(k+2) = (k+2)*u(k+3)*m1
                end do
                ! flip for more accuracy in inner product calculation
                u = [(u(i) ,i = m+1,1,-1)]
                up =[(up(i),i = m+1,1,-1)]
                !Newton iteration
                hh = hh1
                hh(m+1) = m1
                step = maxrp   
                l = 0
                z = 0._rp
                do while ((abs(step)>GlobalEps).and.(l<10))
                    l = l + 1              
                    step = (u.ip.hh)/(up.ip.hh)
                    h = h - step
                    !powers of h (This is the fastest way!)
                    z1 = m1*h+z
                    z2 = cumprod(z1)
                    hh = [m1,(z2(i),i=1,m)]           
                    !flip for more accuracy in inner product calculation
                    hh(1:m+1) = hh(m+1:1:-1)
                end do
                roots(j+1) = x + h
                !ders(j+1) =  dot_product(up,hh)
                ders(j+1) = (up.ip.hh)
            end do   
            !nodes are symmetric
            roots(1:N1+s) = -roots(n:N1+1:-1)
            ders(1:N1+s) = ders(n:N1+1:-1)
            
        end subroutine alg1_Herm
                
        !-------------------------------------------------------
        pure subroutine first_root(n,x,ders)
        integer(ip),intent(in)::            n
        real(rp),dimension(:),intent(out):: x,ders
        real(rp)::                          Hm2, Hm1, Hpm2, Hpm1,H,Hp
        integer(ip)::                       k                    
        
            Hm2 = 0._rp
            Hm1 = pi**(-0.25_rp)
            Hpm2 = 0._rp
            Hpm1 = 0._rp
            do k = 0,n-1
                H = -sqrt(1._rp*k/(k+1))*Hm2
                Hp = sqrt(2._rp/(k+1))*Hm1 - sqrt(1._rp*k/(k+1))*Hpm2
                Hm2 = Hm1
                Hm1 = H
                Hpm2 = Hpm1
                Hpm1 = Hp
            end do
            x=0._rp
            ders=0._rp
            if (ibits(n,0,1)==1) then
                !zero is a root:
                x(ishft((n-1),-1)) = 0._rp  !x((n-1)/2) = 0 
                ders(ishft((n+1),-1)) = Hp  !ders((n+1)/2) = Hp        
            else
                !find first root
                call alg2_Herm(H,n,x(ishft(n,-1)+1),ders(ishft(n,-1)+1))
            end if
        
        end subroutine first_root
        
        !-------------------------------------------------------
        pure subroutine alg2_Herm(H,n,x1,d1)
        real(rp),intent(in)::               H
        integer(ip),intent(in)::            n
        real(rp),intent(out)::              x1,d1
        integer(ip)::                       k,i,l
        integer(ip),parameter::             m = 30  !number of terms in Taylor expansion
        real(rp)::                          m1,c,step,dx,x0,x11,y0
        real(rp),dimension(m)::             z,z1,z2
        real(rp),dimension(10)::            x2
        real(rp),dimension(m+1)::           u,up,x1k
            
            ! find the first root (note H_n'(0) = 0)
            ! advance ODE via Runge-Kutta for initial approx
            !x1 = odeRK2(0._rp,-pi/2,0._rp,n)  
            !   x1 = odeRK2(0._rp,-pi/2,0._rp,n)
            !   odeRK2_TVD_1step(dydx,dx,x0,y0)
            x0 = 0._rp
            x11 = -pi/2._rp
            y0 = 0._rp
            dx = (x11 - x0)/10._rp
            x2 = odeRK2(dydx, dx, x0, y0, 10)
            x1 = x2(10)
            !scaling
            m1 = 1._rp/x1
            !initialise
            u = 0._rp
            up = 0._rp
            
            !recurrence relation for Legendre polynomials
            u(1) = H
            u(3) = -0.5_rp*(2._rp*n+1)*u(1)/m1**2   !attention
            up(1) = 0._rp 
            up(2) = 2._rp*u(3)*m1
            do k = 2,m-2,2
                u(k+3) = (-(2._rp*n+1)*u(k+1)/m1**2 + u(k-1)/m1**4)/((k+1)*(k+2))
                up(k+2) = (k+2)*u(k+3)*m1
            end do
            
            !flip for more accuracy in inner product calculation
            u = [(u(i) ,i = m+1,1,-1)]
            up =[(up(i),i = m+1,1,-1)]
            z = 0._rp                                                  
            z1 = m1*x1+z
            z2 = cumprod(z1)                                          
            x1k = [m1,(z2(i),i=1,m)]
            step = maxrp                                                
            l = 0
            
            !Newton iteration
            do while ((abs(step)>GlobalEps).and.(l<10))
                l = l + 1
                step = (u.ip.x1k)/(up.ip.x1k)
                x1 = x1 - step
                !powers of h (This is the fastest way!)
                z1 = m1*x1+z
                z2 = cumprod(z1)
                x1k = [1._rp,(z2(i),i=1,m)]
                x1k = [(x1k(i),i=m+1,1,-1)]
            end do
            
            !Update derivative
            d1 = (up.ip.x1k)
            
        end subroutine alg2_Herm
        
        !-------------------------------------------------------
        pure real(rp) function dydx(x0,y0) result(y)
        real(rp),intent(in)::       x0,y0
            y=-1._rp/(sqrt(2._rp*n+1-y0**2) - 0.5_rp*y0*sin(2._rp*x0)/(2._rp*n+1-y0**2))
        end function dydx
                
        !-------------------------------------------------------
        pure subroutine hermpts_rec(n)
        integer(ip),intent(in)::            n
        real(rp),dimension(:),allocatable:: x1,w1,v1    
        integer(ip)::                       m
        
            m=ishft(n,-1)
            if(ibits(n,0,1)==1) then
                allocate(x1(m+1),w1(m+1),v1(m+1))
            else
                allocate(x1(m),w1(m),v1(m))
            end if  
            
            call Hermpts_recParameter(n,m,x1,w1,v1)
            
            if (ibits(n,0,1)==1) then   !fold out
                quadx = [-x1(m+1:1:-1),x1(2:m+1)]
                quadw = [w1(m+1:1:-1),w1(2:m+1)]
                quadw = (spi/sum(quadw))*quadw
            else
                quadx = [-x1(m:1:-1),x1(:)]
                quadw = [w1(m:1:-1),w1(:)]
                quadw = (spi/sum(quadw))*quadw
            end if
        
        end subroutine hermpts_rec
        
        !-------------------------------------------------------
        pure subroutine Hermpts_recParameter(n,m,x1,w1,v1)
        integer(ip),intent(in)::            n
        integer(ip),intent(inout)::         m
        real(rp),dimension(:),intent(out):: x1,w1,v1    
        real(rp),dimension(:),allocatable:: x0,val,dval,dx
        integer(ip)::                       kk,i,size_para  
        
            size_para=size(x1)
            allocate(x0(size_para),val(size_para),dval(size_para),dx(size_para))
            
            call HermiteInitialGuesses(m,n,x0)
            x0 = x0*sqrt(2._rp)
            do kk = 1,10
                call hermpoly_rec(n,x0,val,dval)  
                dx = val/dval
                do i = 1,size(dx)
                    if(isnan(dx(i))) dx(i) = 0._rp
                end do
                x0 = x0 - dx
                if(norm(dx,-1)<sqrt(GlobalEps)) exit
            end do
        
            x1 = x0/sqrt(2._rp)
            w1 = (exp(-x1**2)/dval**2)      !quadrature weights
            v1 = exp(-x1**2/2._rp)/dval     !Barycentric weights
            
        end subroutine Hermpts_recParameter
        
        !-------------------------------------------------------
        pure subroutine HermiteInitialGuesses(m,n,x)
        integer(ip),intent(in)::            m,n
        real(rp),dimension(:),intent(out):: x
        integer::                           i,k
        real(rp)::                          a,nu,p,tin
        real(rp),dimension(m)::             airyrts,x_init,x_init_airy,Tnk0,rhs,val,dval,dTnk0,tnk,x_init_sin 
        real(rp),dimension(10)::            airyrts_exact
        
            ! hermiteintitialguesses(n), Initial guesses for Hermite zeros.
            if(ibits(n,0,1)==1) then
                a =  0.5_rp
            else
                a = -0.5_rp
            end if  
            nu = 4*m + 2*a + 2   
            do i =1,m
                airyrts(i) = -T(3._rp/8._rp*pi*(4._rp*i-1._rp))
            end do
        
            ! Exact Airy roots.
            airyrts_exact = [-2.338107410459762_rp,-4.087949444130970_rp,-5.520559828095555_rp,&    
                             -6.786708090071765_rp,-7.944133587120863_rp,-9.022650853340979_rp,&
                            -10.040174341558084_rp,-11.008524303733260_rp,-11.936015563236262_rp,&
                            -12.828776752865757_rp]
            
            airyrts(1:10) = airyrts_exact(1:10)     ! correct first 10.
        
            x_init = real(sqrt(cmplx(nu + 2**(2._rp/3._rp)*airyrts*nu**(1._rp/3._rp) + &
                1._rp/5._rp*2**(4._rp/3._rp)*airyrts**2*nu**(-1._rp/3._rp) + &
                (11._rp/35._rp-a**2-12._rp/175._rp*airyrts**3)/nu + &
                (16._rp/1575._rp*airyrts+92._rp/7875._rp*airyrts**4)*2**(2._rp/3._rp)*nu**(-5._rp/3._rp) - &
                (15152._rp/3031875._rp*airyrts**5+1088._rp/121275._rp*airyrts**2)*2**(1._rp/3._rp)*nu**(-7._rp/3._rp),kind=rp)),kind=rp)  
            x_init_airy(1:m) = x_init(m:1:-1)

            ! Tricomi initial guesses. Equation (2.1) in [1]. Originally in [2].
            ! These initial guesses are good near x = 0 . Note: zeros of besselj(+/-.5,x)
            ! are integer and half-integer multiples of pi.
            ! x_init_bess =  bess/sqrt(nu).*sqrt((1+ (bess.^2+2*(a^2-1))/3/nu^2) );
            Tnk0 = [(pi/2,i=1,m)]
            nu = 4*m+2*a+2
            rhs = [((4*m-4*i+3)/nu*pi,i=1,m)]
            do k = 1,7
                val = Tnk0 - sin(Tnk0) - rhs 
                dval = 1 - cos(Tnk0)
                dTnk0 = val/dval 
                Tnk0 = Tnk0 - dTnk0 
            end do
      
            tnk = cos(Tnk0/2._rp)**2
            x_init_sin = sqrt(nu*tnk - (5._rp/(4._rp*(1._rp-tnk)**2) - 1._rp/(1._rp-tnk)-1._rp+3._rp*a**2)/3._rp/nu)
            !Patch together
            p = 0.4985_rp+GlobalEps
            x_init = [x_init_sin(1:floor(p*n)),x_init_airy(ceiling(p*n):m)]                   
            if (ibits(n,0,1)==1) then
                x = [0._rp,x_init]                                                             
                x = x(1:m+1)
            else
                x = x_init(1:m)
            end if
            
        end subroutine HermiteInitialGuesses
        
        !-------------------------------------------------------
        pure real(rp) function T(x) result(y)
        real(rp),intent(in)::               x
        
            y = x**(2._rp/3._rp)*(1._rp + 5._rp/48._rp*x**(-2) - 5._rp/36._rp*x**(-4) + &
                                (77125._rp/82944._rp)*x**(-6) - 108056875._rp/6967296._rp*x**(-8) + &
                                162375596875._rp/334430208._rp*x**(-10))
        
        end function T
        
        !-------------------------------------------------------
        pure subroutine hermpoly_rec(n,x0,val,dval)
        integer(ip),intent(in)::            n
        real(rp),dimension(:),intent(in)::  x0
        real(rp),dimension(:),intent(out):: val,dval
        integer(ip)::                       k
        real(rp),dimension(size(x0))::      Hold,H,Hnew      
            !HERMPOLY_rec evaluation of scaled Hermite poly using recurrence
            !evaluate:
            Hold = exp(-x0**2/4._rp) 
            H = x0*exp(-x0**2/4._rp)
            do k = 1,n-1
                Hnew = (x0*H/sqrt((k+1)*1._rp) - Hold/sqrt((1+1._rp/k)*1._rp))
                Hold = H 
                H = Hnew
            end do
            !evaluate derivative:
            val = Hnew
            dval = (-x0*Hnew + n**(1._rp/2._rp)*Hold)
        end subroutine hermpoly_rec
        
        !--------------------------------------------------------------------------
        pure subroutine hermpts_asy(n)
        integer(ip),intent(in)::n
            !HERMPTS_ASY, fast algorithm for computing Gauss-Hermite nodes and weights 
            !using Newton's method with polynomial evaluation via asymptotic expansions.
            !x = Gauss-Hermite nodes, w = quad weights, v = bary weights.
            call disableProgram    
        end subroutine hermpts_asy
        
    end subroutine GaussHermitePhys
    
    !----------------------------------------------------------------------------------------------------------------------------------
    pure subroutine GaussHermiteProb(quadx,quadw)
    real(rp),dimension(:),intent(out):: quadx,quadw
   
        call GaussHermitePhys(quadx,quadw)
        quadx = quadx*sqrt(2._rp)
        quadw = quadw*sqrt(2._rp)
       
    end subroutine GaussHermiteProb
    
    !-----------------------------------------------------------------------------------------------
    pure subroutine GaussHermitePhys_GW(quadx,quadw)
    real(rp),dimension(:),intent(out):: quadx,quadw
    integer(ip)::                       n
    
        n=size(quadx)
        call GolubWelsch_eigenvalueMethod(n)
        
    contains
    
        pure subroutine GolubWelsch_eigenvalueMethod(n)
        integer(ip),intent(in)::            n
        integer(ip)::                       i
        integer(ip),dimension(n)::          index_
        integer(ip),dimension(ishft(n,-1))::ii
        real(rp),dimension(n)::             Value1,D_G
        real(rp),dimension(n,n)::           Z_real
        real(rp),dimension(ishft(n,-1))::   x_half,w_half,v_half
        
            call GolubWelschParameter(n,D_G,Z_real,index_)
            quadx = D_G
            quadw = spi*Z_real(1,index_)**2
            !Enforce symmetry
            ii = [(i,i=1,ishft(n,-1))]                                      
            x_half = quadx(ii)
            w_half = quadw(ii)       
            if (ibits(n,0,1)==1) then                                          
                quadx = [x_half(:),0._rp,-x_half(ishft(n,-1):1:-1)]
                quadw = [w_half(:),spi-sum(2._rp*w_half(:)),w_half(ishft(n,-1):1:-1)]
            else
                quadx = [x_half(:),-x_half(ishft(n,-1):1:-1)]
                quadw = [w_half(:),w_half(ishft(n,-1):1:-1)]
            end if    
            
        end subroutine GolubWelsch_eigenvalueMethod
        
        !--------------------------------------------------------------------------
        pure subroutine GolubWelschParameter(n,D_G,Z_real,index_)
        integer(ip),intent(in)::                n
        real(rp),dimension(:),intent(out)::     D_G
        real(rp),dimension(:,:),intent(out)::   Z_real
        integer(ip),dimension(:),intent(out)::  index_
        integer(ip)::                           i
        real(rp),dimension(n-1)::               beta,D
        real(rp),dimension(n)::                 Value1
        real(rp),dimension(n,n)::               T
        complex(rp),dimension(n,n)::            Z
        
            beta = [(sqrt(0.5_rp*i),i=1,n-1)]
            T = diag(beta,1)+diag(beta,-1)                  
            Value1 = [(T(i,i),i=1,n)]
            D = [(T(i,i+1),i=1,n-1)]
            call eigenTriDiagonal(Value1,D,Z)
            Z_real = real(Z,kind=rp)
            D_G = Value1
            call sort(D_G,index_)
            
        end subroutine GolubWelschParameter
        
    end subroutine GaussHermitePhys_GW
    
    !--------------------------------------------------------------------------------------------------
    !https://github.com/chebfun/chebfun/blob/34f92d12ca51003f5c0033bbeb4ff57ac9c84e78/%40chebtech2/chebpts.m
    !https://github.com/chebfun/chebfun/blob/34f92d12ca51003f5c0033bbeb4ff57ac9c84e78/%40chebtech2/quadwts.m
    !also refer to <lwReference>/<Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules>
    !calculate $\int_{-1}^{1} f(x) dx$
    !n points reach n-1 order accuracy | N+1 points reach N orderaccuracy
    subroutine ClenshawCurtis(quadx,quadw)
    real(rp),dimension(:),intent(out)::         quadx,quadw
    integer(ip)::                               i,n,nm1,nhc
    complex(rp),dimension(size(quadx)-1)::      c
    
        n = size(quadx)
        nm1 = n - 1
        nhc = ishft(size(quadx)+1,-1)   !n_half_ceil
        
        !--
        if(n==1) then
            quadx(1) = 0._rp
            quadw(1) = 2._rp
            return
        endif
        
        !--
        forall(i=0:nm1) quadx(i+1) = - cos( i * pi / real(nm1,kind=rp) )
        
        !--Exact integrals of $\int_{-1}^{1} T_k (x) dx (k \in even number of [0,n-1])$
        !Warning...
        !this array constructor leads to stack overflow, if necessary, modify it.
        c(1:nhc) = cmplx( 2._rp/[1._rp , 1._rp - real([2:n-1:2]**2,kind=rp)] , kind=rp )
        c(1:nm1) = [c(1:nhc),c(ishft(n,-1):2:-1)]
        
        call ifft(c)
        
        quadw(1:nm1) = real(c,kind=rp)
        quadw(1) = quadw(1)/2._rp
        quadw(n) = quadw(1)
    
    end subroutine ClenshawCurtis
    
    !a quite slower method, as a reference
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/f_src.html %sparse_grid_cc
    pure subroutine ClenshawCurtisNative(quadx,quadw)
    real(rp),dimension(:),intent(out)::         quadx,quadw
    integer(ip)::                               n,nm1,i,j
    real(rp)::                                  b,theta(size(quadx))
          
        n = size(quadx)
        nm1 = n-1
        forall(i=1:n) theta(i) = (i-1) * pi / nm1
        
        !--
        if(n==1) then
            quadx(1) = 0._rp
            quadw(1) = 2._rp
            return
        endif

        !--
        forall(i=1:n) quadx(i) = - cos(theta(i))
        
        !--
        do i = 1 , n
            quadw(i) = 1._rp
            do j = 1 , nm1 / 2
                b = merge(1._rp , 2._rp , 2*j==nm1)
                quadw(i) = quadw(i) - b * cos(2._rp * j * theta(i)) / (4._rp * j**2 - 1._rp)
            end do
        end do

        quadw = quadw / nm1
        quadw(2:nm1) = 2._rp * quadw(2:nm1)
    
    end subroutine ClenshawCurtisNative
    
    
    
    !----------------------------------------------------
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/sandia_sparse/sandia_sparse.f90
    pure integer(ip) function SparseGridSize(dim,level,quadratureRule) result(n)
    integer(ip),intent(in)::        dim,level
    character(*),intent(in)::       quadratureRule
    character(len(quadratureRule))::rule
    
        rule = quadratureRule
        call lowerstring(rule)
        
        if(rule == 'cc' .or. rule == 'clenshawcurtis') then
            n = cc(dim,level)           !close set| cc[-1,1]
        elseif(rule == 'gh' .or. rule == 'gausshermite' .or. &
            rule == 'gl' .or. rule == 'gausslegendre') then
            n = openWeakNest(dim,level) !open set | legendre(-1,1)
        else
            call disableprogram
        endif
    
    contains
    
        pure integer(ip) function cc(dim,level)
        integer(ip),intent(in)::        dim,level
        integer(ip)::                   i,j,h,t
        integer(ip),dimension(0:level)::newPoint1d
        integer(ip),dimension(dim)::    level1d
        logical(lp)::                   more
        
            if(level<0) then
                cc = 0
            elseif(level==0) then
                cc = 1
            else
                !due to the np for close is | np = merge(1 , 2**level + 1 , level==0) | 1,3,5,9,17,...
                !so for different level | newpoint(i|i>1) = 2**(i-1) | 1,2,2,4,8...
                newPoint1d(0) = 1; newPoint1d(1) = 2
                j = 1
                do i=2,level
                    j = j*2
                    newPoint1d(i) = j
                enddo
                
                cc = 0
                !traverse indices[level1d] norm(level1d,1) <= level | the smolyak rule
                do i=0,level
                    more = .false.; h = 0; t = 0
                    !traverse indices[level1d] norm(level1d,1) = i
                    do; call compositionNext(i,dim,level1d,more,h,t)
                        cc = cc + product(newPoint1d(level1d))
                        if(.not.more) exit
                    enddo
                enddo
            endif
            
        end function cc
        
        !--
        pure integer(ip) function openWeakNest(dim,level) result(sz)
        integer(ip),intent(in)::        dim,level
        integer(ip)::                   i,j,h,t,levelstart
        integer(ip),dimension(dim)::    level1d,newPoint1d
        logical(lp)::                   more
        
            if(level==0) then
                sz = 1; return
            endif
        
            if(dim==1) then
                levelstart = level
                sz = 1
            else
                levelstart = 0
                sz = 0
            endif
            
            do i=levelstart,level
                more = .false.
                h = 0; t = 0
                do; call compositionNext(i,dim,level1d,more,h,t)!level at different dim
                    newPoint1d = NpAtLevelOpen(level1d)
                    do j=1,dim
                        if(newPoint1d(j)>1) newPoint1d(j) = newPoint1d(j) - 1
                        !due to weakly nest, onlyt 1 point nest
                    enddo
                    sz = sz  + product(newPoint1d)
                    if(.not.more) exit
                enddo 
            enddo
        end function openWeakNest
        
    end function sparseGridSize
    
    
    !-----------------------------
    !cfn -> close fully nest    |cc
    !ofn -> open fully nest     |fejer1(f1), fejer2(f2), gauss petterson(gp)
    !onn -> open none nest      |gauss laguree(gla)
    !own -> open weakly nest    |gl,gh
    !np -> number of points | npSg -> number of points for sparse grid | npTensor -> number of points |
    !level -> in each dimension, a level determines the resolustion in that dimension, the points increases with level up
    subroutine sparseGrid(level,QuadratureRule,cubx,cubw)
    integer(ip),intent(in)::                level
    character(*),intent(in)::               QuadratureRule
    real(rp),dimension(:),intent(out)::     cubw
    real(rp),dimension(:,:),intent(out)::   cubx
    integer(ip)::                           dimSparseGrid,npSparseGrid
    character(len(QuadratureRule))::        rule

        !transfer to local string
        rule = QuadratureRule
        call lowerstring(rule)
        if(rule=='clenshawcurtis')  rule = 'cc'
        if(rule=='gausshermite')    rule = 'gh'
        if(rule=='gausslegendre')   rule = 'gl'
        !--------------------------------
        
        dimSparseGrid = size(cubx, 1)
        npSparseGrid= size(cubx, 2)
        
        if(rule == 'cc') then
            !--close fully nested
            call ClosefullyNest(dimSparseGrid,level,npSparseGrid, cubx,cubw)
        elseif(rule == 'gh'.or. rule == 'gl') then
            !--open weakly nested, only nest the middle point
            call OpenWeaklyNest(dimSparseGrid,level,npSparseGrid, rule(1:2), cubx,cubw)
        else
            call disableprogram
        endif
    
    contains
    
        !--
        subroutine CloseFullyNest(dim,maxlvl,npSg,x,w)
        integer(ip),intent(in)::                dim,maxlvl,npSg
        real(rp),dimension(:,:),intent(out)::   x
        real(rp),dimension(:),intent(out)::     w
        integer(ip)::                           npMaxlvl,ipt,id,i
        integer(ip),dimension(dim,npSg)::       gi,gb   !gi->gridindex:  gb->gridbase
            
            !gi map one-dimensional index to sparse multi-dimensional index with index coordinate like (0,0,...0)
            !notice index from 0 to np-1
            call levelIndexCfn(dim,maxlvl,npSg, gi,gb)
            
            npMaxlvl = NpAtLevelClose(maxlvl)
            
            !calculate x
            do ipt = 1,npSg
                do id = 1,dim
                    i = gi(id,ipt) + 1  !from range[0:2**maxlvl] to range [1:2**maxlvl+1]or[1:npMaxlvl]
                    if(npMaxlvl==1) then  ! one point only
                        x(id,ipt) = 0._rp
                    elseif(2*(npMaxlvl-i)==npMaxlvl-1) then !middle point
                        x(id,ipt) = 0._rp
                    else
                        ! pi -> 0, the same as | - cos(i * pi / real(nm1,kind=rp)) | i = 0:nm1
                        x(id,ipt) = cos(real(npMaxlvl - i, kind = rp)*pi &
                                    / real(npMaxlvl - 1, kind = rp))
                    endif
                enddo
            enddo
            
            !calculate w
            call WeightCfn(dim,maxlvl,npSg,gi,w)
        
        end subroutine CloseFullyNest
        !--
        subroutine levelIndexCfn(dim,maxlvl,npSg, gridIndex,gridBase)
        integer(ip),intent(in)::                dim,maxlvl,npSg
        integer(ip),dimension(:,:),intent(out)::gridIndex,gridBase
        integer(ip),dimension(dim)::            lvl_d,np_d
        integer(ip),dimension(:,:),allocatable::giTensor
        integer(ip),dimension(:),allocatable::  glvl
        integer(ip)::                           j,ilvl,n,npTensor,h,t,npmlv,lv4id,s,lv
        logical(lp)::                           more
        
            n = 0
            !traverse indices[lvl_d] norm(lvl_d,1) <= maxlvl | the smolyak rule
            do ilvl=0,maxlvl
                more = .false.; h = 0; t = 0
                
                !traverse indices[lvl_d] norm(lvl_d,1) = ilvl
                do; call compositionNext(ilvl,dim,lvl_d,more,h,t)
                    np_d = NpAtLevelClose(lvl_d)
                    npTensor = product(np_d)  !tensor point number under this lvl_d(1:dim)
                    allocate(giTensor(dim,npTensor),glvl(npTensor))
                    
                    !giTensor map one-dimensional index to tensor multi-dimensional index with index coordinate like (0,0,...0)
                    call TensorGridIndexCfn(dim,np_d, giTensor)
                    
                    !--scale to reflect the level according to multiply index a power of 2
                    !--the indices [giTensor] give the coordinate of point index at maxlvl meaning
                    call TensorGridScaleClose(dim,maxlvl,lvl_d, giTensor)

                    !--Determine the first level of appearance of each of the points, and calculate the norm_1 of the level of point
                    !--glvl(ipoint) = norm(lvl(gi(:,ipoint) , 1) | for sparse grid, we just need the point [glvl(ipoint)<=maxlvl]
                    call NormOneLevelClose(dim,maxlvl,gitensor, glvl)
                    
                    do j=1,npTensor
                        if(glvl(j) == ilvl) then
                            n = n + 1
                            gridBase(:,n) = np_d
                            gridIndex(:,n) = giTensor(:,j)
                        endif
                    enddo
                    
                    deallocate(giTensor,glvl)
                    if(.not.more) exit
                enddo
            enddo
        end subroutine levelIndexCfn
        !--
        subroutine TensorGridIndexCfn(dim,np_d,giTensor)
        integer(ip),intent(in)::                dim
        integer(ip),dimension(:),intent(in)::   np_d        !np in each dimension
        integer(ip),dimension(:,:),intent(out)::giTensor
        integer(ip)::                           i
        logical(lp)::                           more
        integer(ip),dimension(dim)::            a
            more = .false.
            !traverse the tensor grid by first dimension
            do i=1,size(giTensor,dim=2) ! = product(np_d) = npTensor
                !cycle to (0,0,...0) -> (np_d(1)-1,np_d(2)-1,...,np_d(d)-1) | total npTensor
                call colexNext(dim,np_d,a,more)
                giTensor(:,i) = a
            enddo
        end subroutine TensorGridIndexCfn
        !--
        subroutine TensorGridScaleClose(dim,maxlvl,lvl_d,gi)
        integer(ip),intent(in)::                    dim,maxlvl
        integer(ip),dimension(:),intent(in)::       lvl_d
        integer(ip),dimension(:,:),intent(inout)::  gi
        integer(ip)::                               i
            do i=1,dim
                if(lvl_d(i)==0) then 
                    !if lvl==0, the middle point only, and it located in the (npMaxlvl-1)/2 |
                    !maxlvl=2 npMaxlvl = 5 | 0,1,2,3,4 | the middle point located at 2 = (5-1)/2
                    gi(i,:) = ishft(NpAtLevelClose(maxlvl)-1,-1)
                else
                    !lvl_d(i) = 1 | 0,1,2 | maxlvl=2 | 0,1,2,3,4 | || So we have the relation below due to the close points
                    !0{lvl=2} = 0{lvl=1}*2**(2-1)| 2{lvl=2} = 1{lvl=1}*2**(2-1) | 4{lvl=2} = 2{lvl=1}*2**(2-1)
                    gi(i,:) = gi(i,:) * 2**(maxlvl-lvl_d(i))
                endif
            enddo
        end subroutine TensorGridScaleClose
        !--
        subroutine NormOneLevelClose(dim,maxlvl,giTensor,glvl)
        integer(ip),intent(in)::                dim,maxlvl
        integer(ip),dimension(:,:),intent(in):: giTensor
        integer(ip),dimension(:),intent(out)::  glvl
        integer(ip)::                           i,j,s,npMaxlvl,lvl
            if(maxlvl<=0) then
                glvl = 0
            else
                npMaxlvl = NpAtLevelClose(maxlvl)
                !traverse all points at grid
                do j=1,size(giTensor,dim=2)
                    glvl(j) = 0
                    do i=1,dim
                        !giTensor(i,j) range [0:2**maxlvl], and npMaxlvl=2**maxlvl + 1
                        s = modulo(giTensor(i,j),npMaxlvl)
                        if(s==0) then !only giTensor(i,j) = 0, first appearance in lvl = 1
                            lvl = 0
                        else
                            lvl = maxlvl
                            do while(ibits(s,0,1)==0) !if s is odd, then first appearance in this lvl
                                s = ishft(s,-1)
                                lvl = lvl - 1
                            enddo
                        endif
                        if(lvl==0) then !if s=0 or s=2**maxlvl, transfer [lvl=0] to [lvl=1]
                            lvl = 1
                        elseif(lvl==1) then !if s = 2**(maxlvl-1), transfer [lvl=1] to [lvl=0]
                            lvl = 0
                        endif
                        glvl(j) = glvl(j) + lvl !norm one for non-negative integer
                    enddo
                enddo
            endif
        end subroutine NormOneLevelClose
        !--
        subroutine WeightCfn(dim,maxlvl,npSg,gi,w)
        integer(ip),intent(in)::                dim,maxlvl,npSg
        integer(ip),dimension(:,:),intent(in):: gi
        real(rp),dimension(:),intent(out)::     w
        integer(ip)::                           ilvl,h,t,npTensor,i,j
        real(rp)::                              coef
        logical(lp)::                           more
        integer(ip),dimension(dim)::            np_d,lvl_d
        integer(ip),dimension(:,:),allocatable::giTensor
        real(rp),dimension(:),allocatable::     gw
            
            if(maxlvl==0) then
                w = 2._rp ** dim
                return
            endif
            
            w = 0._rp
            !strange index
            do ilvl=max(0,maxlvl+1-dim),maxlvl
            
                more = .false.; h=0; t=0
                
                do; call compositionNext(ilvl,dim,lvl_d,more,h,t) !sum(lvl_d) = ilvl and do recycle
                    np_d = NpAtLevelClose(lvl_d)
                    npTensor = product(np_d)  !tensor point number under this lvl_d(1:dim)
                    allocate(giTensor(dim,npTensor),gw(npTensor))
                    
                    !giTensor map one-dimensional index to tensor multi-dimensional index with index coordinate like (0,0,...0)
                    call TensorGridIndexCfn(dim,np_d, giTensor)
                    
                    !--
                    call TensorWeight(dim,np_d,'cc', gw)
                    
                    !--scale to reflect the level according to multiply index a power of 2
                    call TensorGridScaleClose(dim,maxlvl,lvl_d, giTensor)
                    
                    !--
                    coef = merge(1._rp,-1._rp,ibits(maxlvl-ilvl,0,1)==0) * combination(dim-1,maxlvl-ilvl)
                    do j=1,npTensor
                        do i=1,npSg
                            if(all(giTensor(:,j)==gi(:,i))) w(i) = w(i) + coef*gw(j)
                        enddo
                    enddo
                    
                    deallocate(giTensor,gw)
                    if(.not.more) exit
                enddo
            enddo
        
        end subroutine WeightCfn
        !--
        subroutine OpenWeaklyNest(dim,maxlvl,npSg,rule,x,w)
        integer(ip),intent(in)::                dim,maxlvl,npSg
        character(*),intent(in)::               rule
        real(rp),dimension(:,:),intent(out)::   x
        real(rp),dimension(:),intent(out)::     w
        integer(ip),dimension(dim)::            lvl_d,np_d,gbase
        integer(ip),dimension(:,:),allocatable::giTensor
        integer(ip),dimension(:),allocatable::  glvl
        real(rp),dimension(:,:),allocatable::   gx
        real(rp),dimension(:),allocatable::     gw
        integer(ip)::                           minlvl,ilvl,h,t,npTensor,i,j,ipSg
        logical(lp)::                           more
        real(rp)::                              coef
            
            w = 0._rp; ipSg = 0
            minlvl = max(0 , maxlvl + 1 - dim)
            
            !--
            do ilvl = merge(minlvl,0,dim==1),maxlvl
                more=.false.; h=0; t=0
                
                !sum(lvl_d) = ilvl and do recycle, range [0:ilvl]
                do; call compositionNext(ilvl,dim,lvl_d,more,h,t)
                
                    np_d = NpAtLevelOpen(lvl_d)
                    gbase = pidb2(np_d-1)
                    npTensor = product(np_d)
                    allocate(giTensor(dim,npTensor), gx(dim,npTensor), gw(npTensor), glvl(npTensor))
                    
                    call TensorWeight(dim,np_d,rule, gw,gx)
                    coef = merge(1._rp,-1._rp,ibits(maxlvl-ilvl,0,1)==0) * combination(dim-1,maxlvl-ilvl)
                    
                    !index from(-M,M); M = (np_d(d)-1)/2
                    call TensorGridIndexOwn(dim,np_d, giTensor)
                    
                    !determine each first level(norm 1) the point appearance, glvl(i) = norm(lvl(giTensor(i)),1)
                    call TensorLevelIndexOwn(dim,maxlvl,ilvl,lvl_d,gbase,giTensor, glvl)
                    
                    do i=1,npTensor
                        if(glvl(i)==ilvl) then
                            ipSg = ipSg + 1
                            x(:,ipSg) = gx(:,i)
                            if(minlvl<=ilvl) w(ipSg) = coef*gw(i)
                        elseif(minlvl<=ilvl) then
                            do j=1,ipSg
                                if(all(x(:,j)==gx(:,i))) exit
                            enddo                            
                            w(j) = w(j) + coef*gw(i)
                        endif
                    enddo

                    deallocate(giTensor, gx, gw, glvl)                    
                    if(.not.more) exit
                enddo
            enddo
            
        end subroutine OpenWeaklyNest
        
        !--
        subroutine TensorGridIndexOwn(dim,np_d,giTensor)
        integer(ip),intent(in)::                dim
        integer(ip),dimension(:),intent(in)::   np_d        !np in each dimension
        integer(ip),dimension(:,:),intent(out)::giTensor
        integer(ip)::                           i
        logical(lp)::                           more
        integer(ip),dimension(dim)::            a
            more = .false.
            !traverse the tensor grid by first dimension
            do i=1,size(giTensor,dim=2) ! = product(np_d) = npTensor
                !cycle to (0,0,...0) -> (np_d(1)-1,np_d(2)-1,...,np_d(d)-1) | total npTensor
                call colexNext(dim,np_d,a,more)
                !from (0:np_d(d)-1) to (-M,M) where M = (np_d(d) -1)/2, np_d(i) is always odd
                giTensor(:,i) = a - ishft(np_d - 1, -1)
            enddo
        end subroutine TensorGridIndexOwn
        !--
        subroutine TensorLevelIndexOwn(dim,maxlvl,ilvl,lvl_d,gbase,giTensor, glvl)
        integer(ip),intent(in)::                    dim,maxlvl,ilvl
        integer(ip),dimension(:),intent(in)::       lvl_d,gbase
        integer(ip),dimension(:,:),intent(in)::     giTensor
        integer(ip),dimension(:),intent(out)::      glvl
        integer(ip)::                               minlvl,i,j

            minlvl = merge(maxlvl,0,dim==1)
            do i = 1,size(giTensor,dim=2)
                glvl(i) = max(ilvl,minlvl)
                do j=1,dim
                    !if giTenor==0 the lvl should decrease due to index=0 point come from lvl=0
                    if(giTensor(j,i)==0) glvl(i) = max(glvl(i) - gbase(j), minlvl)
                enddo
            enddo
        
        end subroutine TensorLevelIndexOwn
        
    end subroutine sparseGrid
 
    !----------------------------------------------------------------------
    !    For the Fejer Type 1, Fejer Type 2, and Gauss Patterson rules, the point
    !    growth is nested.  If we have np points on a particular LEVEL, the next
    !    level includes all these old points, plus np+1 new points, formed in the
    !    gaps between successive pairs of old points plus an extra point at each
    !    end.
    !
    !    Level      np = New + Old
    !
    !    0          1   =  1  +  0
    !    1          3   =  2  +  1
    !    2          7   =  4  +  3
    !    3         15   =  8  +  7
    !    4         31   = 16  + 15
    !    5         63   = 32  + 31
    !
    !    If we use a series of Gauss Legendre rules, then there is almost no
    !    nesting, except that the central point is shared.  If we insist on
    !    producing a comparable series of such points, then the "nesting" behavior
    !    is as follows:
    !
    !    Level      np = New + Old
    !
    !    0          1   =  1  +  0
    !    1          3   =  2  +  1
    !    2          7   =  6  +  1
    !    3         15   = 14  +  1
    !    4         31   = 30  +  1
    !    5         63   = 62  +  1
    elemental integer(ip) function NpAtLevelOpen(level) result(np)
    integer(ip),intent(in)::    level ![0:]
        np = 2**(level + 1) - 1
    end function NpAtLevelOpen
    !-----------------------------------------------------------------------------
    !    For the Clenshaw Curtis rules, the point growth
    !    is nested.  If we have ORDER points on a particular LEVEL, the next
    !    level includes all these old points, plus ORDER-1 new points, formed
    !    in the gaps between successive pairs of old points.
    !
    !    Level      Order = New + Old
    !
    !    0          1   =  1  +  0
    !    1          3   =  2  +  1
    !    2          5   =  2  +  3
    !    3          9   =  4  +  5
    !    4         17   =  8  +  9
    !    5         33   = 16  + 17
    elemental integer(ip) function NpAtLevelClose(level) result(np)
    integer(ip),intent(in)::    level ![0:]
        np = merge(1 , 2**level + 1 , level==0)
    end function NpAtLevelClose
    
    !-----------------------------------------------------------------------------
    subroutine TensorWeight(dim,np_d,rule, gw,gx)
    integer(ip),intent(in)::                dim
    integer(ip),dimension(:),intent(in)::   np_d
    character(*),intent(in)::               rule
    real(rp),dimension(:),intent(out)::     gw
    real(rp),dimension(:,:),intent(out),optional::  gx
    real(rp),dimension(:),allocatable::     x,w
    integer(ip)::                           id
    integer(ip)::                           contig,skip,rep,start,j,k
        
        gw = 1._rp
        contig = 1; skip = 1; rep = size(gw)
        do id = 1,dim
            allocate(x(np_d(id)),w(np_d(id)))
            
            if(rule=='cc') then
                call ClenshawCurtis(x,w)
            elseif(rule=='gh') then
                call GaussHermite(x,w)
            elseif(rule=='gl') then
                call GaussLegendre(x,w)
            else
                call disableprogram
            endif
            
            !prod(id,np_d(id),w,dim,npTensor,gw) | size(w) = np_id(id)
            rep = rep / np_d(id)
            skip = skip * np_d(id)
            do j=1,np_d(id)
                start = 1 + (j-1)*contig
                do k=1,rep
                    !gw(np1,np2,np3...npn), dim=k:
                    !w(j) multiple all points(np1*np2..np(k-1)), for(j=1,npk), repeats times (np(k+1)*,...*npn)
                    !skip (np1*np2..npk)
                    gw(start:start+contig-1) = gw(start:start+contig-1) * w(j)
                    if(present(gx)) gx(id,start:start+contig-1)  = x(j)
                    start = start + skip
                enddo
            enddo
            contig = contig * np_d(id)
            
            deallocate(x,w)
        enddo
        
    end subroutine TensorWeight
    
end module IntegrationLib