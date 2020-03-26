!some algorithm for special polynomial construction
module polynomialLib
use constants
use arrayOpsLib
use stringOpsLib
use SpecialFunctionLib
implicit none
    
    private
    public:: mChebAlgo, Stieltjes
	public:: recurCoef
    
contains

    !use modified Chebyshev algorithm to generate the Orthogonal Polynomials which is generally represented
    !by recurrence form: p_{k+1}(x) = (x - alpha_{k})*p_{k}(x) - beta_{k}*p_{k-1}(x) | p_{-1}=0, p_{0}=1
    !input---
    !a(2n-1),b(2n-1): the recursive coef for kernel polynomial h
    !mom(2n): the moment \int{h_i(x), dP(x)} | if a=0 and b=0 -> h_i = x^i which is original chebyshev algo
    !output---
    !alpha,beta: the desired recusive coef
    !s: the normalizer s(i) = \int{p_i(x)**2, dP(x)}
    !ierr:
    !0, normal return,
    !1, if  abs(fnu(0)) is less than the machine zero,
    !2, if  n  is out of range,
    !-k if S(K), k = 0,1,2,...,n-1, is about to underflow
    !+k if S(K) is about to overflow.
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/toms726/toms726.html
    !refer to <On Generating Orthogonal Polynomials>
    pure subroutine mChebAlgo(a, b, mom, alpha, beta, s)
    real(rp),dimension(:),intent(in)::  a,b,mom
    real(rp),dimension(:),intent(out):: alpha,beta,s
    integer(ip)::                       n,k,l,lk
    real(rp),dimension(2*size(alpha)):: s0,s1,s2
    
        n = size(alpha)

        alpha(1) = a(1) + mom(2)/mom(1)
        beta(1) = mom(1)
        if(n==1) return

        s(1) = mom(1)
        s0(1:2*n) = 0._rp
        s1(1:2*n) = mom(1:2*n)

        do k=2,n
            lk=2*n - k + 1
            do l=k, 2*n - k + 1
                s2(l) = s1(l+1) - (alpha(k-1) - a(l))*s1(l) - beta(k-1)*s0(l) + b(l)*s1(l-1)
                if (l==k) s(k) = s2(k)
            end do

            if(abs(s(k))<10._rp*tinrp) call disableprogram
            if(0.1_rp*maxrp < abs(s(k))) call disableprogram

            alpha(k) = a(k) + s2(k+1)/s2(k) - s1(k)/s1(k-1)
            beta(k) = s2(k)/s1(k-1)

            s0(k:lk) = s1(k:lk)
            s1(k:lk) = s2(k:lk)
        end do
        
    end subroutine mChebAlgo

    
    !coef for common-used orthogonal polynomial
	!refer to Tom726
    !p_{k+1}(x) = (x - alpha_{k})*p_{k}(x) - beta_{k}*p_{k-1}(x) | p_{-1}=0, p_{0}=1
    pure subroutine recurCoef(polytype, a, b, par)
    character(*),intent(in)::                   polytype
    real(rp),dimension(0:),intent(out)::        a, b
	real(rp),dimension(:),optional,intent(in):: par
	character(len(polytype))::				    pt
    integer(ip)::                               n, k
	real(rp)::                                  alpha, beta, sab, alpha2, beta2
	
	    pt = lowerStringFc(polytype)
		n = ubound(a, 1)
		
        a = 0._rp
        if(pt=='legendre') then
		
            b(0) = 2._rp
            do k=1,n
                b(k) = 1._rp/(4._rp - 1._rp/k**2)
            enddo
			
        elseif(pt=='chebyshev1') then
		
            b(0) = 4._rp*pi
            if(n==0) return
			
            b(1) = 0.5_rp
            b(2:n) = 0.25_rp
			
		!physical hermite polynomial, weight e^{-x^2}
        elseif(pt=='hermite') then
		
		    b(0) = srpi
			do k=1,n
			    b(k) = 0.5_rp*k
			enddo
		
		elseif(pt=='jacobi') then
		
		    alpha = par(1)
			beta = par(2)
			
			sab = alpha + beta
			a(0) = (beta - alpha)/(sab + 2._rp)
		    b(0) = 2._rp**(sab + 1._rp)*betafc(alpha + 1._rp, beta + 1._rp)
            if(n == 0) return
            
            alpha2 = alpha**2
            beta2 = beta**2
            a(1) = (beta2 - alpha2)/((sab + 2._rp)*(sab + 4._rp))
            b(1) = 4._rp*(alpha + 1._rp)*(beta + 1._rp)/((sab + 3._rp)*(sab + 2._rp)**2)
            
            do k=2,n

                a(k) = 0.25_rp*(beta2 - alpha2) &
                    /(k**2*(1._rp + 0.5_rp*sab/k)*(1._rp + 0.5_rp*(sab + 2._rp)/k))
            
                b(k) = 0.25_rp*(1._rp + alpha/k)*(1._rp + beta/k)*(1._rp + sab/k)&
                    /((1._rp + 0.5_rp*(sab + 1._rp)/k)&
					    *(1._rp + 0.5_rp*(sab - 1._rp)/k)&
						*(1._rp + 0.5_rp*sab/k)**2)
            
            end do
		
		else
		
            call disableprogram
			
        endif
    
    end subroutine recurCoef
    
    !--
    !input
    !n is the number of polynomial, highest x^{n-1}
    !ncap is the number of quadrature point
    !x is the quadrature point
    !w is the quadrature weight
    !output
    !alpha, beta, ierr, what want
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/toms726/toms726.html
    !refer to <Multi-Element Generalized Polynomial Chaos for Arbitrary Probability Measures>
    pure subroutine Stieltjes(x, w, alpha, beta)
    real(rp),dimension(:),intent(in)::  x,w
    real(rp),dimension(:),intent(out):: alpha,beta
    integer(ip)::                       n,ncap,k,m
    real(rp)::                          sum0,sum1,sum2,t
    real(rp),dimension(size(x))::       p0,p1,p2
    
        ncap = size(x)
        n = size(alpha)
        if (n<=0 .or. ncap<n) call disableprogram

        sum0 = sum(w(1:ncap))
        sum1 = w(1:ncap).ip.x(1:ncap)
        alpha(1) = sum1/sum0
        beta(1) = sum0

        if (n == 1) return

        p1(1:ncap) = 0._rp
        p2(1:ncap) = 1._rp
        do k = 1, n-1
            sum1 = 0._rp
            sum2 = 0._rp
            do m = 1, ncap
                if(w(m)/= 0._rp) then
                    p0(m) = p1(m)
                    p1(m) = p2(m)
                    p2(m) = (x(m) - alpha(k))*p1(m) - beta(k)*p0(m)

                    if(0.1_rp*maxrp<abs(p2(m)) .or. 0.1_rp*maxrp<abs(sum2)) call disableprogram

                    t = w(m)*p2(m)*p2(m)
                    sum1 = sum1 + t
                    sum2 = sum2 + t*x(m)
                end if
            end do

            if(abs(sum1)<10._rp*tinrp) call disableprogram

            alpha(k+1) = sum2 / sum1
            beta(k+1) = sum1 / sum0
            sum0 = sum1
        end do

    end subroutine Stieltjes
    
end module polynomialLib