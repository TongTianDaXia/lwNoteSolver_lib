!This lib is used to sample from a prescribed distribution
    
module samplinglib
use constants
implicit none

    private
    
    public::    sampleUniform
    public::    sampleAcceptRejection
    public::    samplePullin    !normal variates with given sample mean and variance
    public::    sampleNormal    !normal variates with box-muller scheme
    public::    sampleBoxMuller
    !--
    public::    randomNum
    public::    randomSeed
    

    !----
    integer::   flagRandom = 0  !1: intrinsic | else: bird
    
    !----
    interface sampleUniform
        procedure:: sampleUniform_default
        procedure:: sampleUniform_range
    end interface sampleUniform
    
    !the acceptance-rejection method
    interface sampleAcceptRejection     
        procedure:: sampleAR_discrete
        procedure:: sampleAR_continuous
    end interface 
    
    !--
    interface randomSeed
        procedure:: randomSeed_input
        procedure:: randomSeed_null
    end interface randomSeed
    
    !--
    interface randomNum
        procedure:: random_intrinsic
        procedure:: random_bird
        procedure:: random_bird_thread
    end interface randomNum

!-------------------------------------------------
contains
!-------------------------------------------------

    subroutine randomSeed_input(seed)
    integer(ip),intent(in)::    seed
    integer(ip)::               i,n
        call random_seed()
        call random_seed(size = n)
        call random_seed(put = [(seed + i*37, i=1,n)])
        flagRandom = 1
    end subroutine randomSeed_input
    
    !--
    subroutine randomSeed_null
        call random_seed()
        flagRandom = 1
    end subroutine randomSeed_null
    
    !----
    real(rp) function random_intrinsic() result(rf)
        call random_number(rf)
    end function random_intrinsic
    
    !--
    !a random number generator, the part of Bird's original program of DSMC
    real(rp) function random_bird(idum) result(rf)
    integer(ip),intent(in)::    idum
    integer(ip)::               ma(55), inext, inextp, iff
    save::                      ma, inext, inextp, iff
    data                        iff/0/
        rf = randomNum(idum, ma, inext, inextp, iff)
    end function random_bird
    
    !--
    real(rp) function random_bird_thread(idum, ma, inext, inextp, iff) result(rf)
    integer(ip),intent(in)::    idum
    integer(ip),intent(inout):: ma(55), inext, inextp, iff
    integer(ip),parameter::     mbig = 1000000000, mseed = 161803398, mz = 0
    real(rp),parameter::        fac = 1.e-9_rp
    integer(ip)::               mj, mk, i, k, ii

        !when idum<0, reset psudo random number
        if(idum<0 .or. iff==0) then
            iff = 1
            mj = mseed - iabs(idum)
            mj = mod(mj, mbig)
            ma(55) = mj
            mk = 1
            do i=1,54
                ii = mod(21*i, 55)
                ma(ii) = mk
                mk = mj - mk
                if(mk<mz) mk = mk + mbig
                mj = ma(ii)
            end do
            do k=1,4
                do i=1,55
                    ma(i) = ma(i) - ma(1 + mod(i+30, 55))
                    if(ma(i)<mz) ma(i) = ma(i) + mbig
                end do
            end do
            inext = 0
            inextp = 31
        end if
        
        !kernel
        do; inext = inext + 1
            if(inext==56) inext = 1
            inextp = inextp + 1
            if(inextp==56) inextp = 1
            
            mj = ma(inext) - ma(inextp)
            if(mj<mz) mj = mj + mbig
            ma(inext) = mj
            rf = mj*fac
            if(rf>1.e-8_rp .and. rf<0.99999999_rp) return
        end do
        
    end function random_bird_thread

    !-------------------------------------------------
    !default
    real(rp) function sampleUniform_default() result(s)
        if(flagRandom==1) then
            s = randomNum()
        else
            s = randomNum(0)
        endif
    end function sampleUniform_default
    
    !range(lo,up)
    real(rp) function sampleUniform_range(lo,up) result(s)
    real(rp),intent(in)::   lo,up
        s = lo + (up - lo)*sampleUniform()
    end function sampleUniform_range
    
    !-------------------------------------------------
    !--Accept-Rejection method, refer to
    !https://zhuanlan.zhihu.com/p/25610149
    integer(ip) function sampleAR_discrete(f1, fmax, xmin, xmax) result(x)   
    real(rp),intent(in)::       fmax
    integer(ip),intent(in)::    xmin,xmax
    interface
        !the probability density function that is needed to input
        pure real(rp) function f1(i)
            import:: ip, rp
            integer(ip),intent(in)::    i
        end function f1
    end interface
        do; x = xmin + int(sampleUniform()*dfloat(xmax + 1 - xmin), ip)
            if(f1(x)/fmax > sampleUniform()) exit
        end do
    end function sampleAR_discrete
    
    !--
    real(rp) function sampleAR_continuous(f1,fmax,xmin,xmax) result(x)
    procedure(absf1)::          f1
    real(rp),intent(in)::       fmax,xmin,xmax
        do; x = xmin + sampleUniform()*(xmax-xmin)
            if(f1(x)/fmax > sampleUniform()) exit
        end do
    end function sampleAR_continuous
    
!-------------------------------------------------
    !generate samples from Gaussian/normal distribution
    function samplePullin(n,um,e) result(u)
    real(rp),intent(in)::       e,um            !e is the variance and um is the mean value of u
    integer(ip),intent(in)::    n               !the number of variates 
    real(rp),dimension(n)::     u
    real(rp)::                  mu,r,b,te,ee    
    integer(ip)::               neta,j,m,i
    real(rp),pointer,dimension(:)::    eprime,v,t
        if(n>2) then
            if(ibits(n - 1, 0, 1)==0) then
                mu = 0._rp
            else if (ibits(n - 1, 0, 1)==1) then
                mu = 0.5_rp
            end if
            r = 0.5_rp*(n - 1) - mu
            neta = int(r + 2._rp*mu, ip)
            allocate (eprime(neta), v(n + 1), t(neta - 1))
            eprime = e        
            v = 0._rp
            t = 0._rp
            do j=1,neta-1
                t(j) = sampleUniform()**(1._rp/(neta - j - mu))
            end do
            do j=1,neta
                do m=1,j-1
                    eprime(j) = eprime(j)*t(m)
                end do
                if(j==neta) exit
                eprime(j) = eprime(j)*(1._rp - t(j))
            end do
            do m=1,int(r)
                b = 2._rp*pi*sampleUniform()
                v(2*m) = sqrt(2._rp*eprime(m))*cos(b)
                v(2*m + 1) = sqrt(2._rp*eprime(m))*sin(b)
            end do
            if(ibits(n - 1, 0, 1)==1) then
                if(sampleUniform()<0.5) then
                    v(n) = sqrt(2._rp*eprime(neta))
                else
                    v(n) = -sqrt(2._rp*eprime(neta))
                end if
            end if
            u(1) = um - sqrt(dfloat(n) - 1._rp)*v(2)/sqrt(1._rp*dfloat(n))
            do i=2,n
                u(i) = u(i - 1) + (sqrt(dfloat(n) + 2._rp - dfloat(i))*v(i) - sqrt(dfloat(n) &
                    - 1._rp*dfloat(i))*v(i + 1))/sqrt(dfloat(n) + 1._rp - dfloat(i))
            end do
            deallocate (eprime,v,t)
        else
            if(sampleUniform()<0.5_rp) then
                ee = sqrt(e)
            else
                ee = -sqrt(e)
            end if
            u(1) = um + ee
            u(2) = 2._rp*um - u(1)
        end if
    end function samplePullin
    
    !--
    real(rp) function sampleNormal()
    real(rp)::                  u1,u2
        u1 = sampleUniform()
        u2 = sampleUniform()
        sampleNormal = sqrt(-2._rp*log(u1))*cos(2._rp*pi*u2)
    end function sampleNormal
    
    !--
    function sampleBoxMuller() result(s)
    real(rp)::                  u1,u2
    real(rp),dimension(2)::     s
        u1 = sampleUniform()
        u2 = sampleUniform()
        s(1) = sqrt(-2._rp*log(u1))*cos(2._rp*pi*u2)
        s(2) = sqrt(-2._rp*log(u1))*sin(2._rp*pi*u2)
    end function sampleBoxMuller
    
end module samplinglib
