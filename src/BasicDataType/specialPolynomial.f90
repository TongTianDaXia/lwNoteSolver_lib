module specialPolynomial_
use constants
use IntegrationLib
use SpecialFunctionLib
use polynomialLib
use polynomial_
implicit none

	private
	
	!--based on the recursive definition P_{i+1} = (x-alpha_i)*P_i - beta_i*P_{i-1}
    public:: orthnormalPolynomial
    public:: orthnormalPolynomialSet
	
	!--Legendre, use recursive method rather than explict expression to 
    !   avoid failure due to large n for binomialCoef/factorial method
    public:: LegendrePolynomial
    public:: LegendrePolynomialSet
    !the following polynomial is normalized for 
    !\int_{-1}^{1} p_n p_m dx = \delta_{nm} | tip: the weight function is 1 here
    !for probability density, the density function is {1/2}
    public:: normalLegendrePolynomial
    public:: normalLegendrePolynomialSet
    
    !--Hermite, use recursive method rather than explict expression
    public:: HermitePolynomial
    public:: HermitePolynomialSet
    !the following polynomial is normalized for
    !\int_{-\infty}^{\infty} p_n p_m e^{-x^2} dx = \delta_{nm}
    !for probability density, the density function is {1/sqrt(2\pi)e^{-x^2/2}}
    public:: normalHermitePolynomial
    public:: normalHermitePolynomialSet
    
	!--
	public:: jacobiPolynomial
	public:: jacobiPolynomialSet
	!weight (1-x)**alpha(1+x)**beta
	public:: normalJacobiPolynomial
    public:: normalJacobiPolynomialSet
	
	
    !--chebyshev
    public:: ChebyshevPolynomialT
    public:: ChebyshevPolynomialTset
    public:: ChebyshevPolynomialT_Clenshaw
	

contains

    
	!--alpha(0:n-1),beta(0:n-1):-> generate poly(n)
    !!--based on the recursive definition P_{i+1} = (x-alpha_i)*P_i - beta_i*P_{i-1}
    pure type(polynomial) function orthnormalPolynomial(alpha,beta) result(poly)
    real(rp),dimension(0:),intent(in)::         alpha,beta
    type(polynomial)::                          x,tmp1,tmp2
    integer(ip)::                               i,n
    real(rp),dimension(size(alpha))::           qx,qw
    real(rp)::                                  normalizer
    
        n = size(alpha)
        
        x = [0._rp,1._rp]
        poly = [1._rp]
        poly = (x - alpha(0))*poly
        if(n==1) return
        
        tmp2 = [1._rp]
        do i=2,n
            tmp1 = poly
            poly = (x-alpha(i-1))*tmp1 - beta(i-1)*tmp2
            tmp2 = tmp1
        enddo
        
        call GaussGw(alpha, beta, qx, qw)
        
        normalizer = sum(poly%funcval(qx)**2*qw)
        poly = poly/sqrt(normalizer)
    
    end function orthnormalPolynomial
	
    !--
    pure function orthnormalPolynomialSet(alpha, beta) result(poly)
    real(rp),dimension(0:),intent(in)::         alpha,beta
    type(polynomial),dimension(0:size(alpha)):: poly
    type(polynomial)::                          x
    integer(ip)::                               i,n
    real(rp),dimension(size(alpha))::           qx,qw
    real(rp)::                                  normalizer
    
        n = size(alpha)
        
        x = [0._rp,1._rp]
        poly(0) = [1._rp]        
        poly(1) = (x - alpha(0))*poly(0)
        if(n==1) return
        
        do i=2,n
            poly(i) = (x-alpha(i-1))*poly(i-1) - beta(i-1)*poly(i-2)
        enddo
        
        call GaussGw(alpha, beta, qx, qw)
        
        do i=0,n
            normalizer = sum(poly(i)%funcval(qx)**2*qw)
            poly(i) = poly(i)/sqrt(normalizer)
        enddo
    
    end function orthnormalPolynomialSet
	

    !--
    !n P_n = (2n-1) x P_{n-1} - (n-1) P_{n-2}
    elemental type(polynomial) function LegendrePolynomial(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial)::                  tm2,tm1,x
    integer(ip)::                       i
	
        if(n<=0) then
            poly = [1._rp]
        elseif(n==1) then
            poly = [0._rp, 1._rp]
        else
            x   = [0._rp, 1._rp]
            tm2 = [1._rp]
            tm1 = x
            do i = 2 , n
                poly = (2._rp*i - 1._rp)*x*tm1 - (i - 1._rp)*tm2
                poly = poly / real(i,kind=rp)
                tm2 = tm1
                tm1 = poly
            enddo
        endif
		
    end function LegendrePolynomial
    
    !--
    pure function LegendrePolynomialSet(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial),dimension(0:n)::   poly
    type(polynomial)::                  x
    integer(ip)::                       i
	
        if(n<=0) then
            poly(0) = [1._rp]
        elseif(n==1) then
            poly(0) = [1._rp]
            poly(1) = [0._rp, 1._rp]
        else
            x = [0._rp,1._rp]
            poly(0) = [1._rp]
            poly(1) = [0._rp, 1._rp]
            do i = 2, n
                poly(i) = (2._rp*i - 1._rp)*x*poly(i-1) - (i - 1._rp)*poly(i-2)
                poly(i) = poly(i)/real(i, rp)
            enddo
        endif
		
    end function LegendrePolynomialSet
    
    !--
    elemental type(polynomial) function normalLegendrePolynomial(n) result(poly)
    integer(ip),intent(in)::                n
	
        poly = sqrt(real(2*n + 1, rp)/2._rp)*LegendrePolynomial(n)
		
    end function normalLegendrePolynomial
    
    !--
    pure function normalLegendrePolynomialSet(n) result(poly)
    integer(ip),intent(in)::                n
    type(polynomial),dimension(0:n)::       poly
    integer(ip)::                           i
	
        poly = LegendrePolynomialSet(n)
        do i=0,n
            poly(i) = sqrt(real(2*i + 1, rp)/2.d0) * poly(i)
        enddo
		
    end function normalLegendrePolynomialSet
    

    !wiki(chebyshev polynomial) T(n) for first kind and U(n) for second kind
    elemental type(polynomial) function ChebyshevPolynomialT(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial)::                  tm2,tm1,x
    integer(ip)::                       i
	
        if(n<=0) then
            poly = [1._rp]
        elseif(n==1) then
            poly = [0._rp, 1._rp]
        else
            x = [0._rp, 1._rp]
            tm2 = [1._rp]
            tm1 = x
            do i = 2, n
                poly = 2._rp*x*tm1 - tm2
                tm2 = tm1
                tm1 = poly
            enddo
        endif
		
    end function ChebyshevPolynomialT
    
    !--
    pure function ChebyshevPolynomialTset(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial),dimension(0:n)::   poly
    type(polynomial)::                  x
    integer(ip)::                       i
	
        if(n<=0) then
            poly(0) = [1._rp]
        elseif(n==1) then
            poly(0) = [1._rp]
            poly(1) = [0._rp, 1._rp]
        else
            x       = [0._rp, 1._rp]
            poly(0) = [1._rp]
            poly(1) = [0._rp, 1._rp]
            do i = 2 , n
                poly(i) = 2._rp*x*poly(i-1) - poly(i-2)
            enddo
        endif
		
    end function ChebyshevPolynomialTset
    
    !--compute sum_0^n(c*ChebPoly), sum = c(0)*T(0)%funcval(x)+......+c(n)*T(n)%funcval(x)
    !use Clenshaw algorithm, see wiki
    !https://github.com/chebfun/chebfun/blob/development/%40chebtech/clenshaw.m
    pure real(rp) function ChebyshevPolynomialT_Clenshaw(x,c) result(s)
    real(rp),intent(in)::               x
    real(rp),dimension(0:),intent(in):: c
    real(rp)::                          bk1,bk2,b,x2
    integer(ip)::                       n,k
	
        bk1 = 0._rp
        bk2 = bk1
        x2 = 2._rp*x   !double
        n = ubound(c,dim=1)
        do k=n,2,-2
            bk2 = c(k) + x2*bk1 - bk2
            bk1 = c(k-1) + x2*bk2 - bk1
        enddo
        if(ibits(n,0,1)==1) then
            b = bk1
            bk1 = c(1) + x2*bk1 - bk2
            bk2 = b
        endif
        s = c(0) + x * bk1 - bk2
		
    end function ChebyshevPolynomialT_Clenshaw
    
    
    !---------------------
    !H_n = 2x H_{n-1} - 2(n-1) H_{n-2}
    !for w(x) = e^{-x^2} | <Hi,Hi>=sqrt(pi) 2^n n!
    elemental type(polynomial) function HermitePhysPolynomial(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial)::                  tm2,tm1,x
    integer(ip)::                       i
	
        if(n<=0) then
            poly = [1._rp]
        elseif(n==1) then
            poly = [0._rp, 2._rp]
        else
            x = [0._rp, 1._rp]
            tm2 = [1._rp]
            tm1 = 2._rp*x
            do i = 2, n
                poly = 2._rp*x*tm1 - 2._rp*(i - 1._rp)*tm2
                tm2 = tm1
                tm1 = poly
            enddo
        endif
		
    end function HermitePhysPolynomial
    
    !--
    pure function HermitePhysPolynomialSet(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial),dimension(0:n)::   poly
    type(polynomial)::                  x
    integer(ip)::                       i
	
        if(n<=0) then
            poly(0) = [1._rp]
        elseif(n==1) then
            poly(0) = [1._rp]
            poly(1) = [0._rp, 2._rp]
        else
            x = [0._rp, 1._rp]
            poly(0) = [1._rp]
            poly(1) = [0._rp,2._rp]
            do i = 2, n
                poly(i) = 2._rp*x*poly(i-1) - 2._rp*(i - 1._rp)*poly(i-2)
            enddo
        endif
		
    end function HermitePhysPolynomialSet
    
    !--
    elemental type(polynomial) function normalHermitePhysPolynomial(n) result(poly)
    integer(ip),intent(in)::                n
	
        poly = (1._rp/sqrt(srpi*2**n*factorial(n)))*HermitePhysPolynomial(n)
		
    end function normalHermitePhysPolynomial
    
    !--
    pure function normalHermitePhysPolynomialSet(n) result(poly)
    integer(ip),intent(in)::                n
    type(polynomial),dimension(0:n)::       poly
    integer(ip)::                           i
	
        poly = HermitePhysPolynomialSet(n)
        do i=0,n
            poly(i) = (1._rp/sqrt(srpi*2**i*factorial(i)))*poly(i)
        enddo
		
    end function normalHermitePhysPolynomialSet
    
    
    !---------------------
    !H_n = x H_{n-1} - (n-1) H_{n-2}
    !for w(x) = e^{-x^2/2} | <Hi,Hi>=sqrt(2*pi) n!
    elemental type(polynomial) function HermiteProbPolynomial(n) result(poly)
    integer(ip),intent(in)::                n
    type(polynomial)::                      tm2, tm1, x
    integer(ip)::                           i
	
        if(n<=0) then
            poly = [1._rp]
        elseif(n==1) then
            poly = [0._rp,1._rp]
        else
            x   = [0._rp,1._rp]
            tm2 = [1._rp]
            tm1 = x
            do i = 2 , n
                poly = x * tm1 - (i - 1._rp) * tm2
                tm2 = tm1
                tm1 = poly
            enddo
        endif
		
    end function HermiteProbPolynomial
    
    !--
    pure function HermiteProbPolynomialSet(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial),dimension(0:n)::   poly
    type(polynomial)::                  x
    integer(ip)::                       i
	
        if(n<=0) then
            poly(0) = [1._rp]
        elseif(n==1) then
            poly(0) = [1._rp]
            poly(1) = [0._rp, 1._rp]
        else
            x   = [0._rp, 1._rp]
            poly(0) = [1._rp]
            poly(1) = [0._rp, 1._rp]
            do i = 2 , n
                poly(i) = x*poly(i-1) - (i - 1._rp)*poly(i-2)
            enddo
        endif
		
    end function HermiteProbPolynomialSet
    
    !--
    elemental type(polynomial) function normalHermiteProbPolynomial(n) result(poly)
    integer(ip),intent(in)::                n
	
        poly = (1._rp/sqrt(sqrt(2._rp)*srpi*factorial(n)))*HermiteProbPolynomial(n)
		
    end function normalHermiteProbPolynomial
    
    !--
    pure function normalHermiteProbPolynomialSet(n) result(poly)
    integer(ip),intent(in)::                n
    type(polynomial),dimension(0:n)::       poly
    integer(ip)::                           i
	
        poly = HermiteProbPolynomialSet(n)
        do i=0,n
            poly(i) = (1._rp/sqrt(sqrt(2._rp)*srpi*factorial(i)))*poly(i)
        enddo
		
    end function normalHermiteProbPolynomialSet
    
	
    !--
    elemental type(polynomial) function hermitePolynomial(n,var) result(p)
    integer(ip),intent(in)::            n
    character(*),optional,intent(in)::  var
    
        if(present(var)) then
            if(var=='phys') then
                p = HermitePhysPolynomial(n)
            elseif(var=='prob') then
                p = HermiteProbPolynomial(n)
            else
                call disableprogram
            endif
        else
            p = HermitePhysPolynomial(n)
        endif
    
    end function hermitePolynomial
    
    !--
    pure function HermitePolynomialSet(n,var) result(p)
    integer(ip),intent(in)::            n
    character(*),optional,intent(in)::  var
    type(polynomial),dimension(0:n)::   p
    
        if(present(var)) then
            if(var=='phys') then
                p = HermitePhysPolynomialSet(n)
            elseif(var=='prob') then
                p = HermiteProbPolynomialSet(n)
            else
                call disableprogram
            endif
        else
            p = HermitePhysPolynomialSet(n)
        endif
    
    end function HermitePolynomialSet
    
    !--
    elemental type(polynomial) function normalHermitePolynomial(n,var) result(p)
    integer(ip),intent(in)::            n
    character(*),optional,intent(in)::  var
    
        if(present(var)) then
            if(var=='phys') then
                p = normalHermitePhysPolynomial(n)
            elseif(var=='prob') then
                p = normalHermiteProbPolynomial(n)
            else
                call disableprogram
            endif
        else
            p = normalHermitePhysPolynomial(n)
        endif
    
    end function normalHermitePolynomial
    
    !--
    pure function normalHermitePolynomialSet(n,var) result(p)
    integer(ip),intent(in)::            n
    character(*),optional,intent(in)::  var
    type(polynomial),dimension(0:n)::   p
    
        if(present(var)) then
            if(var=='phys') then
                p = normalHermitePhysPolynomialSet(n)
            elseif(var=='prob') then
                p = normalHermiteProbPolynomialSet(n)
            else
                call disableprogram
            endif
        else
            p = normalHermitePhysPolynomialSet(n)
        endif
    
    end function normalHermitePolynomialSet
	
	
	!--
	!weight (1-x)^\alpha(1+x)^\beta
	!refer to 
	!https://people.sc.fsu.edu/~jburkardt/f_src/jacobi_polynomial/jacobi_polynomial.f90
	type(polynomial) function jacobiPolynomial(n, alpha, beta) result(poly)
	integer(ip),intent(in)::        n
	real(rp),intent(in)::           alpha, beta
	type(polynomial)::              x, tmp1, tmp2
	real(rp)::                      tmp
	integer(ip)::                   i
	
	    if(alpha <= -1._rp)&
		stop 'error: specialPolynomial/jacobiPolynomial get error alpha'
		
	    if(beta <= -1._rp)&
		stop 'error: specialPolynomial/jacobiPolynomial get error beta'
		
		if(n<=0) then
            poly = [1._rp]
        elseif(n==1) then
            poly = [alpha - beta, 2._rp + alpha + beta]/2._rp
        else
            x   = [0._rp, 1._rp]
            tmp2 = [1._rp]
            tmp1 = [alpha - beta, 2._rp + alpha + beta]/2._rp
            do i = 2, n
			    tmp  = 2._rp*i + alpha + beta
                poly = (tmp - 1._rp)*(alpha**2 - beta**2 + tmp*(tmp - 2._rp)*x)*tmp1
				poly = poly - 2._rp*(i - 1._rp + alpha)*(i - 1._rp + beta)*tmp*tmp2
				poly = poly/(2._rp*i*(i + alpha + beta)*(tmp - 2._rp))
                tmp2 = tmp1
                tmp1 = poly
            enddo
        endif
	
	end function jacobiPolynomial
	
	!--
	function jacobiPolynomialSet(n, alpha, beta) result(poly)
	integer(ip),intent(in)::            n
	real(rp),intent(in)::               alpha, beta
	type(polynomial),dimension(0:n)::   poly
	type(polynomial)::                  x
	real(rp)::                          tmp
	integer(ip)::                       i
	
	    if(alpha <= -1._rp)&
		stop 'error: specialPolynomial/jacobiPolynomial get error alpha'
		
	    if(beta <= -1._rp)&
		stop 'error: specialPolynomial/jacobiPolynomial get error beta'
		
		if(n<=0) then
            poly(0) = [1._rp]
        elseif(n==1) then
		    poly(0) = [1._rp]
            poly(1) = [alpha - beta, 2._rp + alpha + beta]/2._rp
        else
            x = [0._rp, 1._rp]
			poly(0) = [1._rp]
            poly(1) = [alpha - beta, 2._rp + alpha + beta]/2._rp
            do i = 2, n
			    tmp  = 2._rp*i + alpha + beta
                poly(i) = (tmp - 1._rp)*(alpha**2 - beta**2 + tmp*(tmp - 2._rp)*x)*poly(i-1)
				poly(i) = poly(i) - 2._rp*(i - 1._rp + alpha)*(i - 1._rp + beta)*tmp*poly(i-2)
				poly(i) = poly(i)/(2._rp*i*(i + alpha + beta)*(tmp - 2._rp))
            enddo
        endif
	
	end function jacobiPolynomialSet
	
	!--
	type(polynomial) function normalJacobiPolynomial(n, alpha, beta) result(poly)
	integer(ip),intent(in)::        n
	real(rp),intent(in)::           alpha, beta
	real(rp)::                      normer
	 
	    !it's not beta function
        normer = 2._rp**(alpha + beta + 1._rp)/(2._rp*n + alpha + beta + 1._rp)
		normer = normer*gamma(n + alpha + 1._rp)*gamma(n + beta + 1._rp)
		normer = normer/factorial(n)/gamma(n + alpha + beta + 1._rp)
		
		poly = jacobiPolynomial(n, alpha, beta)/sqrt(normer)
	    
	end function normalJacobiPolynomial
	
	!--
	function normalJacobiPolynomialSet(n, alpha, beta) result(poly)
	integer(ip),intent(in)::            n
	real(rp),intent(in)::               alpha, beta
	type(polynomial),dimension(0:n)::   poly
	real(rp)::                          normer
	integer(ip)::                       i
	
	    poly = jacobiPolynomialSet(n, alpha, beta)
		
		do i=0,n
		    normer = 2._rp**(alpha + beta + 1._rp)/(2._rp*i + alpha + beta + 1._rp)
		    normer = normer*gamma(i + alpha + 1._rp)*gamma(i + beta + 1._rp)
		    normer = normer/factorial(i)/gamma(i + alpha + beta + 1._rp)
		    poly(i) = poly(i)/sqrt(normer)
		enddo
	    
	end function normalJacobiPolynomialSet
	
end module specialPolynomial_