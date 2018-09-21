module SpecialFunctionLib
use constants
use arrayOpsLib
implicit none

    private
    public:: factorial
    public:: binomialCoef, combination
    public:: besseljn_roots
    
    
!-------------------------------------------
    interface factorial
        procedure::  factorial_n
        procedure::  factorial_kn
    end interface factorial

    !--binomial coef is combinational number
    interface combination
        procedure::  BinomialCoef_int
        procedure::  BinomialCoef_general
    end interface combination
    
    interface BinomialCoef
        procedure::  BinomialCoef_int
        procedure::  BinomialCoef_general
    end interface BinomialCoef
    
    
contains

    !due to factorial increase very rapidly, so give a bigest integer type
    pure integer(8) function factorial_n(n) result(factorial)
    integer(ip),intent(in)::    n
    integer(ip)::               i
        factorial = 1
        if(n<=20) then
            do i=1,n
                factorial = factorial*i
            enddo
        else
            !integer(8) limit of factoral(20)
            call disableProgram
            !fast factorial, see wiki(Stirling's approximation, second order)
            factorial = nint( sqrt(2._rp*pi*n) * (n/e)**n )
        endif
    end function factorial_n
    
    pure integer(8) function factorial_kn(k,n) result(factorial)
    integer(ip),intent(in)::    k,n
    integer(ip)::               i
        factorial = 1
        do i=k,n
            factorial = factorial*i
        enddo
    end function factorial_kn
    
!--------------------------------------------------
    !refer to wiki <binomial coe>
    !----------------------avoid overflow and roundoff of integer
    !coef = factorial(n-k+1,n) / factorial(k)
    pure real(rp) function BinomialCoef_int(n,k) result(coef)
    integer(ip),intent(in)::    n,k
    integer(ip)::               mn,mx,i
        mn = min(k,n-k)
        if(mn<0) then
            coef = 0._rp
        elseif(mn==0) then
            coef = 1._rp
        else
            mx = max(k,n-k)
            coef = real(mx+1,kind=rp)
            do i=2,mn
                coef = coef * real(mx+i,kind=rp) / real(i,kind=rp)
            enddo
        endif
    end function BinomialCoef_int
    !----------------------
    pure real(rp) function BinomialCoef_general(a,k) result(coef)
    real(rp),intent(in)::       a
    integer(ip),intent(in)::    k
    integer(ip)::               i
    real(rp)::                  up
        up = 1._rp
        do i =1,k
            up = up*(a-i+1)
        enddo
        coef = up / dfloat( factorial(k) )
    end function BinomialCoef_general
    
    
!---------------------------------------------------
    !https://github.com/chebfun/chebfun/blob/34f92d12ca51003f5c0033bbeb4ff57ac9c84e78/besselroots.m
    pure function besseljn_roots(jn,n) result(rts)
    integer(ip),intent(in)::    jn,n
    real(rp),dimension(n)::     rts,b
    real(rp)::                  mu,a1,a3,a5,a7,a9,a11,a13
    integer(ip)::               i
    real(rp),dimension(20),parameter:: j0_roots20 = [           &
                                        2.4048255576957728_8,   &
                                        5.5200781102863106_8,   &
                                        8.6537279129110122_8,   &
                                        11.791534439014281_8,   &
                                        14.930917708487785_8,   &
                                        18.071063967910922_8,   &
                                        21.211636629879258_8,   &
                                        24.352471530749302_8,   &
                                        27.493479132040254_8,   &
                                        30.634606468431975_8,   &
                                        33.775820213573568_8,   &
                                        36.917098353664044_8,   &
                                        40.058425764628239_8,   &
                                        43.199791713176730_8,   &
                                        46.341188371661814_8,   &
                                        49.482609897397817_8,   &
                                        52.624051841114996_8,   &
                                        55.765510755019979_8,   &
                                        58.906983926080942_8,   &
                                        62.048469190227170_8]

! McMahon's expansion. This expansion gives very accurate approximation 
! for the sth zero (s >= 7) in the whole region V >=- 1, and moderate
! approximation in other cases.
        mu = 4._8*jn**2
        a1 = 1._8 / 8._8
        a3 = (7._8*mu-31._8)/384._8
        a5 = 4._8*(3779._8+mu*(-982._8+83._8*mu))/61440._8 !Evaluate via Horner's method
        a7 = 6._8*(-6277237._8+mu*(1585743._8+mu*(-153855._8+6949._8*mu))) / 20643840._8
        a9 = 144._8*(2092163573._8+mu*(-512062548._8+mu*(48010494._8 +  &
            mu*(-2479316._8+70197._8*mu)))) / 11890851840._8
        a11 = 720._8*(-8249725736393._8+mu*(1982611456181._8+mu*(-179289628602._8 +  &
            mu*(8903961290._8 + mu*(-287149133._8+5592657._8*mu))))) / 10463949619200._8
        a13 = 576._8*(423748443625564327._8 + mu*(-100847472093088506._8+mu*(8929489333108377._8 + &
            mu*(-426353946885548._8+mu*(13172003634537._8+mu*(-291245357370._8 + &
            mu*4148944183._8)))))) / 13059009124761600._8
        
        b = (4._8*[(i,i=1,n)] + 2._8*jn - 1._8 ) * pi * 0.25_8
        do i=1,n
            rts(i) = b(i) - (mu-1._8)*&
                polyval([0._rp,a1,0._rp,a3,0._rp,a5,0._rp,a7,0._rp,a9,0._rp,a11,0._rp,a13], 1._rp/b(i))
        enddo
        
        !correction
        if(jn == 0) rts(1:min(n,20)) = j0_roots20(1:min(n,20))
        !we negelect the correction for jn<1.or.jn>5
    !---------------------------------------------------------------------------------------------------------------------
    !real(rp),dimension(180),parameter:: c = [ &
    !   2.883975316228_8,  8.263194332307_8, 11.493871452173_8, 14.689036505931_8, 17.866882871378_8, 21.034784308088_8, &
    !   0.767665211539_8,  4.209200330779_8,  4.317988625384_8,  4.387437455306_8,  4.435717974422_8,  4.471319438161_8, &
    !  -0.086538804759_8, -0.164644722483_8, -0.130667664397_8, -0.109469595763_8, -0.094492317231_8, -0.083234240394_8, &
    !   0.020433979038_8,  0.039764618826_8,  0.023009510531_8,  0.015359574754_8,  0.011070071951_8,  0.008388073020_8, &
    !  -0.006103761347_8, -0.011799527177_8, -0.004987164201_8, -0.002655024938_8, -0.001598668225_8, -0.001042443435_8, &
    !   0.002046841322_8,  0.003893555229_8,  0.001204453026_8,  0.000511852711_8,  0.000257620149_8,  0.000144611721_8, &
    !  -0.000734476579_8, -0.001369989689_8, -0.000310786051_8, -0.000105522473_8, -0.000044416219_8, -0.000021469973_8, &
    !   0.000275336751_8,  0.000503054700_8,  0.000083834770_8,  0.000022761626_8,  0.000008016197_8,  0.000003337753_8, &
    !  -0.000106375704_8, -0.000190381770_8, -0.000023343325_8, -0.000005071979_8, -0.000001495224_8, -0.000000536428_8, &
    !   0.000042003336_8,  0.000073681222_8,  0.000006655551_8,  0.000001158094_8,  0.000000285903_8,  0.000000088402_8, &
    !  -0.000016858623_8, -0.000029010830_8, -0.000001932603_8, -0.000000269480_8, -0.000000055734_8, -0.000000014856_8, &
    !   0.000006852440_8,  0.000011579131_8,  0.000000569367_8,  0.000000063657_8,  0.000000011033_8,  0.000000002536_8, &
    !  -0.000002813300_8, -0.000004672877_8, -0.000000169722_8, -0.000000015222_8, -0.000000002212_8, -0.000000000438_8, &
    !   0.000001164419_8,  0.000001903082_8,  0.000000051084_8,  0.000000003677_8,  0.000000000448_8,  0.000000000077_8, &
    !  -0.000000485189_8, -0.000000781030_8, -0.000000015501_8, -0.000000000896_8, -0.000000000092_8, -0.000000000014_8, &
    !   0.000000203309_8,  0.000000322648_8,  0.000000004736_8,  0.000000000220_8,  0.000000000019_8,  0.000000000002_8, &
    !  -0.000000085602_8, -0.000000134047_8, -0.000000001456_8, -0.000000000054_8, -0.000000000004_8,              0._8, &
    !   0.000000036192_8,  0.000000055969_8,  0.000000000450_8,  0.000000000013_8,              0._8,              0._8, &
    !  -0.000000015357_8, -0.000000023472_8, -0.000000000140_8, -0.000000000003_8,              0._8,              0._8, &
    !   0.000000006537_8,  0.000000009882_8,  0.000000000043_8,  0.000000000001_8,              0._8,              0._8, &
    !  -0.000000002791_8, -0.000000004175_8, -0.000000000014_8,              0._8,              0._8,              0._8, &
    !   0.000000001194_8,  0.000000001770_8,  0.000000000004_8,              0._8,              0._8,              0._8, &
    !  -0.000000000512_8, -0.000000000752_8,              0._8,              0._8,              0._8,              0._8, &
    !   0.000000000220_8,  0.000000000321_8,              0._8,              0._8,              0._8,              0._8, &
    !  -0.000000000095_8, -0.000000000137_8,              0._8,              0._8,              0._8,              0._8, &
    !   0.000000000041_8,  0.000000000059_8,              0._8,              0._8,              0._8,              0._8, &
    !  -0.000000000018_8, -0.000000000025_8,              0._8,              0._8,              0._8,              0._8, &
    !   0.000000000008_8,  0.000000000011_8,              0._8,              0._8,              0._8,              0._8, &
    !  -0.000000000003_8, -0.000000000005_8,              0._8,              0._8,              0._8,              0._8, &
    !   0.000000000001_8,  0.000000000002_8,              0._8,              0._8,              0._8,              0._8]
    end function besseljn_roots

end module SpecialFunctionLib