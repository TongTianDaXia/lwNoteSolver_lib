module stringOpsLib
use constants
implicit none

    private
    !--
    public:: readKeyVal
    public:: countsubstring, lowerString, upperString, bondString
    public:: i2c
    
    
    !--------------------------------------------------
    interface readKeyVal
        procedure::  readKeyVal_integer
        procedure::  readKeyVal_integer1rank
        procedure::  readKeyVal_integer2rank
        procedure::  readKeyVal_real
        procedure::  readKeyVal_char
    end interface readKeyVal
    
    
!=======================================
contains
!=======================================

    !--
    elemental subroutine lowerString(str)
    character(*),intent(inout)::str
    integer::                   i,n
        do i=1,len(str)
            n = ichar(str(i:i))
            if(n>=uppercase_a.and.n<=uppercase_z) then
                str(i:i) = char(ichar(str(i:i)) - uppercase_a + lowercase_a)
            endif
        enddo
    end subroutine lowerString
    
    !--
    elemental subroutine upperString(str)
    character(*),intent(inout)::str
    integer::                   i,n
        do i=1,len(str)
            n = ichar(str(i:i))
            if(n>=lowercase_a.and.n<=lowercase_z) then
                str(i:i) = char(ichar(str(i:i)) - lowercase_a + uppercase_a)
            endif
        enddo
    end subroutine upperString
    
    !--
    elemental function countsubstring(string,substring) result(n)
    character(*),intent(in)::   string,substring
    integer(ip)::               n,loc,ls,i
    
        n   = 0
        ls  = len(substring)
        if( ls <= 0 .or. ls > len(string)) return

        loc = 1
        do while(.true.)
            i = index(string(loc:),substring)
            if(i > 0) then
                n   = n + 1
                loc = loc + i - 1 + ls
            else
                exit
            endif
        enddo
        
    end function countsubstring
    
    !--
    pure function bondString(s1,s2)
    character(*),intent(in)::       s1,s2
    character(len(s1)+len(s2))::    bondString
        bondString = s1//s2
    end function bondString
    
    !===================================================
    pure subroutine readKeyVal_char(string,key,val)
    character(*),intent(in)::   string,key
    character(*),intent(out)::  val
    integer(ip)::               st
    
        st = index(string,key)
        if(st == 0 ) then 
            val = 'null'; return
        endif
        st = st + len(key)
        
        !check '='
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            val = 'null'; return
        endif
        st = st + 1

        read(string(st:),*) val
        
    end subroutine readKeyVal_char

    !--    
    pure subroutine readKeyVal_real(string,key,val)
    character(*),intent(in)::   string,key
    real(rp),intent(out)::      val
    integer(ip)::               st
    
        st = index(string,key)
        if(st == 0 ) then 
            call disableNumber(val); return
        endif
        st = st + len(key)
        
        !check '='
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            call disableNumber(val); return
        endif
        st = st + 1
        
        read(string(st:),*) val
        
    end subroutine readKeyVal_real
    
    !--
    pure subroutine readKeyVal_integer(string,key,val)
    character(*),intent(in)::   string,key
    integer(ip),intent(out)::   val
    integer(ip)::               st
    
        st = index(string,key)
        if(st == 0 ) then 
            call disableNumber(val); return
        endif
        st = st + len(key)
        
        !check '='
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            call disableNumber(val); return
        endif
        st = st + 1
        
        read(string(st:),*) val
        
    end subroutine readKeyVal_integer
    
    !--read | Key = (n,n,n)
    pure subroutine readKeyVal_integer1rank(string,key,val)
    character(*),intent(in)::   string,key
    integer(ip),dimension(:),intent(out)::   val
    integer(ip)::               st,ed
    
        st = index(string,key)
        if(st == 0 ) then 
            call disableNumber(val); return
        endif
        st = st + len(key)
        
        !check '='
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            call disableNumber(val); return
        endif
        st = st + 1

        !check '('
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='(') then
            call disableNumber(val); return
        end if
        st = st + 1
        
        read(string(st:),*) val
        
    end subroutine readKeyVal_integer1rank
    
    !--read | Key = (n,n,n)(n,n,n) | rank(3,2)
    pure subroutine readKeyVal_integer2rank(string,key,val)
    character(*),intent(in)::   string,key
    integer(ip),dimension(:,:),intent(out)::   val
    integer(ip)::               st,n,i
    
        st = index(string,key)
        if(st == 0 ) then 
            call disableNumber(val); return
        endif
        
        st = st + len(key)
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            call disableNumber(val); return
        endif
        st = st + 1
        
        n  = size(val,dim=2)
        do i=1,n
            st = verify(string(st:),' ') + st - 1
            if(string(st:st)/='(') then
                call disableNumber(val); return
            end if
            st = st + 1
            read(string(st:),*) val(:,i)
            st = index(string(st:),')') + st
        enddo
        
    end subroutine readKeyVal_integer2rank
    
    
    !--
    character(20) function i2c(i)
    integer(ip),intent(in)::    i
    integer(ip)::               absi, maxExp, pos, p, k
    character(5)::              fmt
        
        !Refer to max number of integer. it can't denot number above 1e**(maxExp+1)
        if(ip==4) then
            maxExp = 9
        elseif(ip==8) then
            maxExp = 18
        else
            stop 'error: stringOpsLib/i2c function reject integer precision'
        endif
        
        absi = abs(i)
        pos = 1
        do k=1,maxExp
            pos = pos*10
            if(absi<pos) then
                p = k
                if(i<0) p = p + 1
                if(p<10) then
                    write(fmt, '("(i",i1,")")') p
                else
                    write(fmt, '("(i",i2,")")') p
                endif
                exit
            endif
        enddo
        
        write(i2c,trim(fmt)) i
    
    end function i2c

end module stringOpsLib