!refer to https://gcc.gnu.org/onlinedocs/gfortran/STRUCTURE-and-RECORD.html
!period(.) is an old component access, and use percent(%) for
module constants
use,intrinsic:: ieee_arithmetic
implicit none
    
    !----------------------Global Accuracy control list----------------------
    integer,parameter ::            isp = selected_int_kind(9)
    integer,parameter ::            idp = selected_int_kind(13)
    integer,parameter ::            rsp = selected_real_kind(p=6,r=37)
    integer,parameter ::            rdp = selected_real_kind(p=15,r=307)
    integer,parameter::             ip  = isp
    integer,parameter::             rp  = rdp
    integer,parameter::             lp  = ip
    integer,parameter::             cl  = 32
    
    
    !----------------------character control list----------------------
    !this is standard(like *) output of different types, as a reference here
    character(10),parameter::       fmt_rdp  = '(E23.15E3)'
    character(9),parameter::        fmt_rsp  = '(E13.6E2)'
    character(5),parameter::        fmt_idp  = '(I20)'
    character(5),parameter::        fmt_isp  = '(I11)'
    integer(ip),parameter::         lowercase_a = ichar('a')
    integer(ip),parameter::         lowercase_z = ichar('z')
    integer(ip),parameter::         uppercase_a = ichar('A')
    integer(ip),parameter::         uppercase_z = ichar('Z')

    
    !----------------------physical and mathmatic constants-----------------------
    real(rp),parameter::            zero    = 0._rp
    real(rp),parameter::            pi      = 4._rp*atan(1._rp)
    real(rp),parameter::            spi     = sqrt(pi)
    real(rp),parameter::            pi2     = pi**2
    real(rp),parameter::            e       = exp(1._rp)
    real(rp),parameter::            k_b     = 1.38067852e-23_rp !boltzmann constant
    real(rp),parameter::            R_c     = 8.3144598_rp      !gas constant   ( J * [K^-1] * [mol^-1] )
    real(rp),parameter::            R_air   = 287.058_rp        !specific gas constant for ideal gas ( J * [kg^-1] * [mol^-1] )
    real(rp),parameter::            gm_diatomic = 1.4_rp        !specific gas
    real(rp),parameter::            P_atm   = 101325._rp        !pressure at one atmosphere
    
    
    !----------------------numeric constant-----------------------------------------
    integer(isp),parameter::        Minisp  = -2147483648_isp
    integer(idp),parameter::        Minidp  = -9223372036854775808_idp
    real(rdp),parameter::           minrdp  = -huge(1._rdp)
    real(rsp),parameter::           minrsp  = -huge(1._rsp)
    !--
    integer(ip),parameter::         maxip   = huge(1_ip)
    integer(ip),parameter::         minip   = maxip + 1_ip
    real(rp),parameter::            maxrp   = huge(1._rp)
    real(rp),parameter::            minrp   = - maxrp
    real(rp),parameter::            nanrp   = transfer(-1_rp,0._rp)
    real(rp),parameter::            infrp   = maxrp * (1._rp + epsilon(1._rp))
    real(rp),parameter::            tinrp   = tiny(1._rp)
    !--
    real(rp),parameter::            GlobalEps = epsilon(1._rp) * 10._rp
    

    
!----------------------------------------------------------------------------------------
    !a special ==, because of its unlimited polymorphism, don't override the intrinsic ==
    interface operator(.lweq.)
        procedure:: anyiseq
    end interface
    
    
    !------
    interface swap
        procedure:: swapInt
        procedure:: swapReal
        procedure:: swapChar
        procedure:: swapLogi
    end interface
    
    
    !here we offer two methods for disabling program [disableProgram]&[disableNumber]
    !and correspondingly offer a inquire function to check 
    !if the number has beed disabled[disableNumber]
    !Tips: don't use module procedure refer to
    !https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/721674
    interface disableProgram
        procedure:: disableProgram_
    end interface disableProgram

    !--
    interface disableNumber
        procedure::  rsp_nan
        procedure::  rdp_nan
        procedure::  isp_nan
        procedure::  idp_nan
    end interface disableNumber
    
    !--
    interface disableStat
        procedure::  disableStat_rsp
        procedure::  disableStat_rdp
        procedure::  disableStat_isp
        procedure::  disableStat_idp
    end interface disableStat
    
    !this is a basic abstract procedure, as a sample and a common procedureType here
    abstract interface
        elemental real(rp) function absf1(x) result(y)
        import:: rp
        real(rp),intent(in)::   x
        end function absf1
        
        elemental real(rp) function absf2(x,y) result(z)
        import:: rp
        real(rp),intent(in)::   x,y
        end function absf2
    end interface
    
!===================
contains
!===================

    !--compare two scalar byte by byte
    elemental logical(lp) function anyiseq(lhs,rhs) result(r)
    class(*),intent(in)::               lhs,rhs
    integer(1),dimension(sizeof(lhs)):: lb
    integer(1),dimension(sizeof(rhs)):: rb
    integer(ip)::                       i
        r   = .false.
        if(.not.same_type_as(lhs,rhs))  return
        if(.not.size(lb)==size(rb))     return
        lb  = transfer(lhs,mold=lb)
        rb  = transfer(rhs,mold=rb)
        r   = all(lb==rb)
    end function anyiseq
    
    !--positive integer mod 2
    elemental integer(ip) function pimod2(i)
    integer(ip),intent(in)::    i
        pimod2 = ibits(i,0,1)
    end function pimod2
    
    !positive integer divided by 2
    elemental integer(ip) function pidb2(i)
    integer(ip),intent(in)::    i
        pidb2 = ishft(i,-1)
    end function pidb2
    
    !--
    elemental real(rp) function log2(a)
    real(rp),intent(in)::       a
        log2 = log(a) / log(2._rp)
    end function log2
    
    
    !-------------------------------
    pure subroutine disableProgram_
    integer(ip),dimension(:),allocatable:: n
        n(1) = 0
    end subroutine disableProgram_

    !--
    elemental subroutine rsp_nan(r)
    real(rsp),intent(out)::     r
        r = nanrp
    end subroutine rsp_nan

    elemental subroutine rdp_nan(r)
    real(rdp),intent(out)::     r
        r = nanrp
    end subroutine rdp_nan
    
    elemental subroutine isp_nan(i)
    integer(isp),intent(out)::  i
        i = minisp
    end subroutine isp_nan
    
    elemental subroutine idp_nan(i)
    integer(idp),intent(out)::  i
        i = minidp
    end subroutine idp_nan
    
    !--
    elemental logical(lp) function disableStat_rsp(r) result(l)
    real(rsp),intent(in)::      r
        l = isnan(minrsp)
    end function disableStat_rsp
    
    elemental logical(lp) function disableStat_rdp(r) result(l)
    real(rdp),intent(in)::      r
        l = isnan(minrdp)
    end function disableStat_rdp
    
    elemental logical(lp) function disableStat_isp(i) result(l)
    integer(isp),intent(in)::   i
        l = i == minisp
    end function disableStat_isp
    
    elemental logical(lp) function disableStat_idp(i) result(l)
    integer(idp),intent(in)::   i
        l = i == minidp
    end function disableStat_idp
    
    !-----------
    elemental subroutine swapInt(a,b)
    integer(ip),intent(inout):: a,b
    integer(ip)::               t
        t = a; a = b; b = t
    end subroutine swapInt
    
    elemental subroutine swapReal(a,b)
    real(rp),intent(inout)::    a,b
    real(rp)::                  t
        t = a; a = b; b = t
    end subroutine swapReal
    
    elemental subroutine swapChar(a,b)
    character(*),intent(inout)::        a,b
    character(len=max(len(a),len(b))):: t
        t = a; a = b; b = t
    end subroutine swapChar
    
    elemental subroutine swapLogi(a,b)
    logical(lp),intent(inout):: a,b
    logical(lp)::               t
        t = a; a = b; b = t
    end subroutine swapLogi
    
    !--
    subroutine writefunc(func,lo,up,np,filename)
    procedure(absf1)::                  func
    real(rp),intent(in)::               lo,up
    integer(ip),intent(in)::            np
    character(*),optional,intent(in)::  filename
    character(cl)::                     fn
    integer(ip)::                       i
    real(rp)::                          x,dx
    
        if(present(filename)) then
            fn = filename
        else
            fn = 'func'
        endif
        dx = (up-lo)/(np-1)
        
        open(unit=99, file=trim(fn)//'.dat')
        write(99,*) 'variables=x,'//trim(fn)
        do i=1,np
            x = lo + (i-1)*dx
            write(99,*) x, func(x)
        enddo
        close(99)
    
    end subroutine writefunc
    
end module constants