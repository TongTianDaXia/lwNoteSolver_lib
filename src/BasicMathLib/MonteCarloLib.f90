!refer to <2015 - Giles, Michael B. - Multilevel Monte Carlo methods>
module MonteCarloLib
use constants
use arrayOpsLib
use SamplingLib
use laWrapperLib
use clock_
implicit none


    interface gmlmc
        procedure:: geoMultiLevelMonteCarlo
    end interface gmlmc

    
    !--
    !for sum(1:4) = sum(Y^[1:4]), sum(5:6) = sum(P_fine^[1:2]) | Y = P_fine - P_coarse
    abstract interface
        pure function mlmcfn(lvl,N) result(s)
        import:: ip,rp
        integer(ip),intent(in)::    lvl,N
        real(rp),dimension(6)::     s
        end function mlmcfn
    end interface
    
contains

    !alpha -> weak error ~ O(2^{-alpha*l})
    !beta -> variance ~ O(2^{-beta*l})
    !gamma -> sample cost ~ O(2^{gamma*l})
    !sumsl(1,:) = sum(Y)
    !sumsl(2,:) = sum(Y^2)
    !eps -> rms error
    pure subroutine geoMultiLevelMonteCarlo(fn,maxlvl,eps,nl0,alpha0,beta0,gamma,r,nl,lvl1)
    procedure(mlmcfn)::                         fn
    integer(ip),intent(in)::                    maxlvl
    integer(ip),dimension(:),intent(in)::       nl0         !the number of samples on lvl 0,1,2
    real(rp),intent(in)::                       eps,alpha0,beta0,gamma
    real(rp),intent(out)::                      r
    integer(ip),intent(out)::                   lvl1,nl(maxlvl)
    integer(ip)::                               i,l,l1,lvl
    integer(ip),dimension(maxlvl)::             dnl,nsmp    !
    real(rp),dimension(maxlvl)::                ml,vl,cl,x  !mean and variance
    real(rp)::                                  alpha,beta,rem,sums(6),sumsl(6,maxlvl)
    integer(ip),dimension(:,:),allocatable::    idx
    real(rp),dimension(:,:),allocatable::       A
    
        alpha = max(0._rp,alpha0)
        beta = max(0._rp,beta0)
        
        lvl = 2; lvl1 = 3
        dnl = nl0(:)
        nl = 0
        sumsl = 0._rp
        
        do while(sum(dnl(1:lvl1))>0 .and. lvl1<maxlvl)
            
            do l = 0,lvl
                l1 = l+1
                if(dnl(l1)>0) then  !increase sampling number
                    sums = fn(l, dnl(l1))
                    nl(l1) = nl(l1) + dnl(l1)
                    sumsl(:,l1) = sumsl(:,l1) + sums(:)
                endif
            enddo
            
            !evaluate
            ml(1:lvl1) = sumsl(1,1:lvl1) / real(nl(1:lvl1),kind=rp)
            vl(1:lvl1) = max(0._rp , sumsl(2,1:lvl1)/real(nl(1:lvl1),kind=rp) - ml(1:lvl1)**2)
            
            !fix
            do l = 3,lvl1
                ml(l) = max(ml(l) , 0.5_rp*ml(l-1)/2._rp**alpha)    !extra
                vl(l) = max(vl(l) , 0.5_rp*vl(l-1)/2._rp**beta)     !extra
            enddo
            
            !linear regressian estimate alpha and beta
            !define the complete matrix for solving Ax=b
            if(alpha0<=0._rp) then
                call repmat(real([1:lvl],kind=rp), 1, 2, A)
                call repmat(trans([1:0:-1]), lvl, 1, idx)
                A = A**idx
                x(1:lvl) = log2(ml(2:lvl1))
                call solveGeneralLES(A,x(1:lvl))
                alpha = max(0.5_rp, -x(1))
            endif
            
            if(beta0<=0._rp) then
                call repmat(real([1:lvl],kind=rp), 1, 2, A)
                call repmat(trans([1:0:-1]), lvl, 1, idx)
                A = A**idx
                x(1:lvl) = log2(vl(2:lvl1))
                call solveGeneralLES(A,x(1:lvl))
                beta = max(0.5_rp, -x(1))
            endif
            
            !--
            cl(1:lvl1) = 2._rp**(gamma*[0:lvl])
            nsmp(1:lvl1) = ceiling(2._rp*vl(1:lvl1)/cl(1:lvl1) * sum(sqrt(vl(1:lvl1)*cl(1:lvl1))) / eps**2)
            dnl(1:lvl1) = max(0, nsmp(1:lvl1) - nl(1:lvl1))
            
            !--estimate remaning error and decide whether a new level is requreed
            !--see section 3.1| E(P-P_L) = sum_{\ell=L+1}^{\infty}E(P_\ell - P_{\ell-1}) = E(P_L-P{L-1})/(2^\alpha-1)
            if(all(dnl < 0.01_rp*nl)==.true.) then
                rem = maxval(ml(lvl1-2:lvl1) * 2._rp**(alpha*[-2:0])) / (2._rp**alpha - 1._rp)
                if(rem>eps/sqrt(2._rp) .and. lvl1<maxlvl) then
                    lvl = lvl1
                    lvl1 = lvl + 1
                    vl(lvl1) = vl(lvl) / 2**beta
                    nl(lvl1) = 0
                    sumsl(:,lvl1) = 0._rp
                    
                    cl(1:lvl1) = 2._rp**(gamma*[0:lvl])
                    nsmp(1:lvl1) = ceiling(2._rp*vl(1:lvl1)/cl(1:lvl1) * sum(sqrt(vl(1:lvl1)*cl(1:lvl1))) / eps**2)
                    dnl(1:lvl1) = max(0, nsmp(1:lvl1) - nl(1:lvl1))
                endif
            endif
        enddo
        
        r = sum(sumsl(1,1:lvl1)/nl(1:lvl1))
        
    end subroutine geoMultiLevelMonteCarlo
    
    !cf -> cost factor ~ 2**gamma
    !kurt -> kurtosis
    !N: sample number for convergence test
    !lvl: level number for convergence test
    !dN0: initial delta number for lvl 0,1,2
    subroutine geoMultiLevelMonteCarloTest(fn,maxlvl,cf,N,lvl,dn0,epsAr)
    procedure(mlmcfn)::                     fn
    integer(ip),dimension(:),intent(in)::   dN0
    integer(ip),intent(in)::                maxlvl,N,lvl
    real(rp),dimension(:),intent(in)::      epsAr
    real(rp),intent(in)::                   cf
    integer(ip)::                           nl(maxlvl),lvlc
    real(rp),dimension(lvl+1)::             del1,del2,var1,var2,cost,kurt,check
    integer(ip)::                           l,l1,l2,lvl1,i,maxl
    real(rp)::                              sums(6),eps,alpha,beta,gamma,poly(0:1),ef
    type(wclock)::                          clock
        
        lvl1 = lvl + 1
        
        del1 = 0._rp
        del2 = 0._rp
        var1 = 0._rp
        var2 = 0._rp
        
        !maybe helpful
        do l = 0, lvl
        
            l1 = l + 1
            
            call clock%clear()
            call clock%start()
            sums = fn(l, N)
            call clock%pauseReco()
            cost(l1) = clock%wtimeCost()
            
            sums = sums / real(N,kind=rp)
            
            !see http://mathworld.wolfram.com/SampleVarianceDistribution.html
            !and http://mathworld.wolfram.com/Kurtosis.html
            kurt(l1) = (sums(4) - 4*sums(3)*sums(1) + 6*sums(2)*sums(1)**2 &
                        - 3*sums(1)**4)/(sums(2)-sums(1)**2)**2
            
            del1(l1) = sums(1)
            del2(l1) = sums(5)
            var1(l1) = sums(2) - sums(1)**2
            var2(l1) = sums(6) - sums(5)**2
            var2(l1) = max(var2(l1), 1.e-12_rp) !fix for var=0
            
            if(l==0) then
                check(l1) = 0._rp
            else
                !see section 3.2
                check(l1) = abs(del1(l1)+del2(l)-del2(l1))&
                            /(3*(sqrt(var1(l1)) + sqrt(var2(l)) + sqrt(var2(l+1))))
                check(l1) = check(l1) / sqrt(real(N,kind=rp))
            endif
            
        enddo
        
        !estimate alpha,beta,gamma
        poly = polyfit(real([0:lvl],kind=rp), log2(abs(del1(1:lvl+1))), 1)
        alpha = -poly(1)
        
        poly = polyfit(real([0:lvl],kind=rp), log2(abs(var1(1:lvl+1))), 1)
        beta = -poly(1)
        
        gamma = log2(cost(lvl+1)/cost(lvl))
        
        
        !-check info
        if(maxval(check)>1._rp) then
            print*, 'warning: geoMLMC test maximum consistency error:'
            print*, maxval(check)
            print*, 'identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied'
        endif
        
        if(kurt(lvl+1)>100._rp) then
            print*, 'warning: geoMLMC test kurtosis on finest level:'
            print*, kurt(lvl+1)
            print*, 'indicates MLMC correction dominated by a few rare paths'
        endif
        
        
        maxl = 0
        do i=1,size(EpsAr)
            eps = epsAr(i)
            gamma = log2(cf)
            call geoMultiLevelMonteCarlo(fn,maxlvl,eps,dn0,alpha,beta,gamma,ef,nl,lvlc)
            maxl = max(lvlc-1, maxl)
        enddo
    
    end subroutine geoMultiLevelMonteCarloTest
    
    
end module MonteCarloLib