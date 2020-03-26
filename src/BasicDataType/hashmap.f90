!refer to:
!https://en.wikipedia.org/wiki/List_of_hash_functions
!src\OpenFOAM\primitives\hashes\Hasher // 
!https://github.com/vkandy/jenkins-hash-java/blob/master/src/JenkinsHash.java
!https://www.cnblogs.com/napoleon_liu/archive/2010/12/29/1920839.html
!https://github.com/nncarlson/petaca/issues/2
!https://www.cnblogs.com/napoleon_liu/archive/2010/12/29/1920839.html
!
!some tisps from fortran on bit: no unsigned type in fortran
!b'1' -> binary constants
!o'7' -> octal constants
!z'f' -> Hexadecimal Constants
!n#number -> different interpreted number
!=> 16#number == z'number'
!--------------------------------------------------------------------
!some bit operations
!iand   ==  &
!ior    ==  |
!ieor   ==  ^
!ishft  ==  <<
!--
!func<same_type_as>
module hashmap_
use constants
use arrayOpsLib
implicit none

    private
    public:: hashmap
    public:: hashfunc
    
    !----------------------------------------------------------------------------------------------------
    integer(isp),parameter::    MapSizeRepository(28) = [               &
            17,         37,         79,         163,        331,        &
            673,        1361,       2729,       5471,       10949,      &
            21911,      43853,      87719,      175447,     350899,     &
            701819,     1403641,    2807303,    5614657,    11229331,   &
            22458671,   44917381,   89834777,   179669557,  359339171,  &
            718678369,  1437356741, 2147483647 ]

!--------------------------------------------------------------------
    type:: hashEntry
        private
        !key is stored by encoded Byte Array, don't care about its type
        integer(1),dimension(:),&
        allocatable::               key_
        class(*),allocatable::      val_
        !--
        class(hashEntry),pointer::  next_   => null()
    contains
    
        procedure::                 depth   => depth_Entry
        procedure::                 delet   => delet_Entry
        procedure::                 clear   => clear_Entry
        procedure::                 set     => set_Entry
        procedure::                 get     => get_Entry
        procedure::                 iseq    => iseq_Entry
        
    end type hashEntry
    
    
!--------------------------------------------------------------------
    type:: hashmap
        private
        class(hashEntry),dimension(:),allocatable:: buckets_
        class(hashEntry),pointer::                  iterEntry_ => null()
    contains
    
        generic::                   init    => init_n
        procedure,private::         init_n
        
        !-----------
        procedure::                 selectMapSize
        procedure::                 mapsize
        procedure::                 bid
        procedure::                 maxcol
        procedure::                 ncol
        procedure::                 nkey
        
        !--
        generic::                   set => set_si
        procedure::                 set_si
        !--
        generic::                   get => get_si,get_iter
        procedure::                 get_si

        !----
        procedure::                 get_iter
        procedure::                 iterEntry
        procedure::                 startIterEntry
        
        !----
        final::                     f_hashmap
    end type hashmap

    
!--------------------------------------------------------------------
    interface hashfunc
        procedure:: javahash
        procedure:: b3hs_hash_key_jenkins
    end interface hashfunc
!--------------------------------------------------------------------
    
    
    
contains
    
!-------------------------------------------------------------------------
!space for hash entry
!-------------------------------------------------------------------------
    integer(ip) function depth_Entry(this) result(depth)
    class(hashEntry),target,intent(in)::    this
    class(hashEntry),pointer::              ptr
        depth = 0
        if(.not.allocated(this%key_)) return
        depth = 1;  ptr => this
        do while(associated(ptr%next_))
            depth = depth + 1
            ptr => ptr%next_
        enddo
    end function depth_Entry
    
    !--
    pure recursive subroutine set_Entry(this,key,val)
    class(hashEntry),intent(inout)::        this
    integer(1),dimension(:),intent(in)::    key
    class(*),intent(in)::                   val
        if(.not.allocated(this%key_)) then
            allocate(this%key_,source=key)
            allocate(this%val_,source=val)
        else
            if(all(this%key_==key)) then
                deallocate(this%val_)
                allocate(this%val_,source=val)
            else
                if(.not.associated(this%next_)) allocate(this%next_)
                call this%next_%set(key,val)
            endif
        endif
    end subroutine set_Entry
    
    !--
    recursive function get_Entry(this,key) result(val)
    class(hashEntry),target,intent(in)::    this
    integer(1),dimension(:),intent(in)::    key
    class(*),pointer::                      val
        if(.not.allocated(this%key_)) then
            val => null()
        elseif(all(this%key_==key)) then
            val => this%val_
        elseif(associated(this%next_)) then
            val => this%next_%get(key)
        endif
    end function get_Entry
    
    !--not allocated key, no next
    pure recursive subroutine delet_Entry(this,key,previous)
    class(hashEntry),intent(inout)::            this
    integer(1),dimension(:),intent(in)::        key
    class(hashEntry),optional,intent(inout)::   previous
    class(hashEntry),pointer::                  tp
        !means the start Entry in Buckets
        if(.not.allocated(this%key_)) return
        if(all(this%key_ == key)) then  !need to delet
            if(present(previous)) then  !means it's the son list
                tp => previous%next_
                previous%next_ => this%next_
                deallocate(tp)
            else! means start Entry in Buckets
                if(associated(this%next_)) then
                    tp => this%next_
                    deallocate(this%key_,this%val_)
                    allocate(this%val_,source=tp%val_)
                    this%key_   = tp%key_
                    this%next_  => tp%next_
                    deallocate(tp)
                else
                    deallocate(this%key_,this%val_)
                endif
            endif
        else
            if(associated(this%next_)) call this%next_%delet(key,this)
        endif
    end subroutine delet_Entry
    
    !--
    pure recursive subroutine clear_Entry(this)
    class(hashEntry),intent(inout)::    this
        if(associated(this%next_)) then
            call this%next_%clear
            deallocate(this%next_)
            this%next_ => null()
        endif
        deallocate(this%key_,this%val_)
    end subroutine clear_Entry
    
    !--
    pure logical(lp) function iseq_Entry(this,key,val) result(is)
    class(hashEntry),intent(in)::           this
    integer(1),dimension(:),intent(in)::    key
    class(*),intent(in)::                   val
        is = .false.
        if(.not.allocated(this%key_))   return
        if(size(key)/=size(this%key_))  return
        if(.not.all(key==this%key_))    return
        is = this%val_ .lweq. val
    end function iseq_Entry

    
    
!------------------------------------------------------------------------
!space for class hashmap
!--------------------------------------------------------------------------
    pure subroutine init_n(this,n)
    class(hashmap),intent(out)::    this
    integer(ip),intent(in)::        n
        allocate(this%buckets_( this%selectMapSize(n)))
    end subroutine init_n

    !---------------------------------
    pure integer(ip) function mapsize(this)
    class(hashmap),intent(in)::     this
        mapsize = size(this%buckets_)    
    end function mapsize
    
    pure integer(ip) function selectMapSize(this,n)
    class(hashmap),intent(in)::     this
    integer(ip),intent(in)::        n
    integer(ip)::                   i
        do i=1,28
            if(mapsizeRepository(i)>n) exit
        enddo
        selectMapSize = mapsizeRepository(i)
    end function selectMapSize
    
    !-----------------------------------
    pure integer(ip) function bid(this,key)
    class(hashmap),intent(in)::             this
    integer(1),dimension(:),intent(in)::    key
        bid = modulo(hashfunc(key),this%mapsize()) + 1
    end function bid
    
    integer(ip) function maxcol(this) result(mc)
    class(hashmap),intent(in)::             this
    integer(ip)::                           i,d
        mc = 0
        do i=1,this%mapsize()
            d = this%buckets_(i)%depth()
            if(d>max(mc,2)) mc = d - 1
        enddo
    end function maxcol
    
    integer(ip) function ncol(this) result(nc)
    class(hashmap),intent(in)::             this
    integer(ip)::                           i,d
        nc = 0
        do i=1,this%mapsize()
            d = this%buckets_(i)%depth()
            if(d<2) cycle
            nc = nc + d - 1
        enddo
    end function ncol
    
    !--
    integer(ip) function nkey(this) result(nk)
    class(hashmap),intent(in)::             this
    integer(ip)::                           i
        nk = 0
        do i=1,this%mapsize()
            nk = nk + this%buckets_(i)%depth()
        enddo
    end function nkey
    
    !--
    pure subroutine set_si(this,key,val)
    class(hashmap),intent(inout)::          this
    integer(ip),intent(in)::                key
    class(*),intent(in)::                   val
    integer(1),dimension(sizeof(key))::     tkey
        tkey = transfer(key,mold=tkey)
        call this%buckets_( this%bid(tkey) )%set(tkey,val)
    end subroutine set_si
    
    !--
    function get_si(this,key)
    class(hashmap),target,intent(in)::      this
    integer(ip),intent(in)::                key
    integer(1),dimension(sizeof(key))::     tkey
    class(*),pointer::                      get_si
        tkey = transfer(key,mold=tkey)
        get_si => this%buckets_( this%bid(tkey) )%get(tkey)
    end function get_si
    
    !-------------------------------------------
    function get_iter(this)
    class(hashmap),target,intent(in)::      this
    class(*),pointer::                      get_iter
        if(associated(this%iterEntry_)) then
            get_iter => this%iterEntry_%val_
        else
            get_iter => null()
        endif
    end function get_iter
    
    !--
    pure subroutine startIterEntry(this)
    class(hashmap),target,intent(inout)::   this
    integer(ip)::                           i
        do i=1,this%mapsize()
            if(allocated(this%buckets_(i)%key_)) then
                this%iterEntry_ => this%buckets_(i)
                return
            endif
        enddo
        this%iterEntry_ => null()
    end subroutine startIterEntry
    
    !--
    elemental subroutine IterEntry(this)
    class(hashmap),target,intent(inout)::   this
    integer(ip)::                           i,ic
        if(associated(this%iterEntry_%next_)) then
            this%iterEntry_ => this%iterEntry_%next_
            return
        else
            ic = this%bid(this%iterEntry_%key_)
            do i=ic+1,this%mapsize()
                if(allocated(this%buckets_(i)%key_)) then
                    this%iterEntry_ => this%buckets_(i)
                    return
                endif
            enddo
        endif
        this%iterEntry_ => null()
    end subroutine iterEntry
    
    
    
    !--------------------
    pure subroutine f_hashmap(this)
    type(hashmap),intent(inout)::   this
    integer(ip)::                   i
        do i=1,this%mapsize()
            call this%buckets_(i)%clear()
        enddo
    end subroutine f_hashmap
    !----------------------
    
    
    
    
    
    
    
    
!--------------------------------------------------------------
!space for different hash function
!---------------------------------------------------------------
    !https://www.zhihu.com/question/51784530
    pure integer(isp) function javahash(key)
    integer(1),dimension(:),intent(in)::    key
        javahash = bkdrhash(key)
        javahash = ieor(javahash,ieor(ishft(javahash,-20),ishft(javahash,-12)))
        javahash = ieor(ieor(javahash,ishft(javahash,-7)),ishft(javahash,-4))
    end function javahash
    
    !https://www.zhihu.com/question/20507188
    pure integer(isp) function bkdrhash(key)
    integer(1),dimension(:),intent(in)::key
    integer(isp),parameter::            s = 131
    integer(isp)::                      n,i
        n = size(key); bkdrhash = 0_isp
        do i=1,n
            bkdrhash = bkdrhash * s + key(i)
        enddo
        bkdrhash = iand(bkdrhash,z'7fffffff')
    end function bkdrhash
    
!---------------------------------------------------------------
    !http://computer-programming-forum.com/49-fortran/0596e59d0fa2e5e4.htm
    pure function b3hs_hash_key_jenkins (key, range) result (code) 
    character(*), intent(in)    :: key 
    integer(isp), intent(in)    :: range 
    integer(isp)                :: code 
    integer(isp)                :: len_key 
    integer(isp)                :: a 
    integer(isp)                :: b 
    integer(isp)                :: c 
    integer(isp)                :: k 
! hash the key into a code, using the algorithm 
! described by bob jenkins at: 
! http://burtleburtle.net/bob/hash/doobs.html 
! 
! note that range should be a power of 2, and 
! that the 32-bit algorithm is used
        len_key = len_trim(key)
        a = z'9e3779b9'
        b = a 
        c = z'12345678'
        k = 1 
        char_loop : do while(len_key>=12) 
! pack the key into 32 bits 
            a = a + ichar(key(k+0:k+0))  + ishft(ichar(key(k+1:k+1)), 8) + & 
            &       ishft(ichar(key(k+2:k+2)), 16) + ishft(ichar(key(k+3:k+3)), 24) 
            b = b + ichar(key(k+4:k+4))  + ishft(ichar(key(k+5:k+5)), 8) + & 
            &       ishft(ichar(key(k+6:k+6)), 16) + ishft(ichar(key(k+7:k+7)), 24) 
            c = c + ichar(key(k+8:k+8))  + ishft(ichar(key(k+9:k+9)), 8) + & 
            &       ishft(ichar(key(k+10:k+10)), 16) + ishft(ichar(key(k+11:k+11)), 24) 
! mix it up 
            call b3hs_hash_key_jenkins_mix_(a,b,c)
            k = k + 12 
            len_key = len_key - 12 
        end do char_loop 
        c = c + len_key
        
! process remaining bits
        select case(len_key)
        case(11) 
            c = c + ishft(ichar(key(k+10:k+10)), 24) + ishft(ichar(key(k+9:k+9)), 16) + & 
            &       ishft(ichar(key(k+8:k+8)), 8) 
            b = b + ishft(ichar(key(k+7:k+7)), 24) + ishft(ichar(key(k+6:k+6)), 16) + & 
            &       ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(10) 
            c = c + ishft(ichar(key(k+9:k+9)), 16) + ishft(ichar(key(k+8:k+8)), 8) 
            b = b + ishft(ichar(key(k+7:k+7)), 24) + ishft(ichar(key(k+6:k+6)), 16) + & 
            &       ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(9) 
            c = c + ishft(ichar(key(k+8:k+8)), 8) 
            b = b + ishft(ichar(key(k+7:k+7)), 24) + ishft(ichar(key(k+6:k+6)), 16) + & 
            &       ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(8) 
            b = b + ishft(ichar(key(k+7:k+7)), 24) + ishft(ichar(key(k+6:k+6)), 16) + & 
            &       ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(7) 
            b = b + ishft(ichar(key(k+6:k+6)), 16) + ishft(ichar(key(k+5:k+5)), 8) + & 
            &       ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(6) 
            b = b + ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(5) 
            b = b + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(4) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(3) 
            a = a + ishft(ichar(key(k+2:k+2)), 16) + ishft(ichar(key(k+1:k+1)), 8) + & 
            &       ichar(key(k:k)) 
        case(2) 
            a = a + ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(1) 
            a = a + ichar(key(k:k))
        end select 
        call b3hs_hash_key_jenkins_mix_(a,b,c) 
        code = iand(c, range - 1) + 1

    contains
        pure subroutine b3hs_hash_key_jenkins_mix_(a,b,c)
        integer(isp),intent(inout)::    a,b,c 
            a = ieor(a - b - c, ishft(c, -13))
            b = ieor(b - c - a, ishft(a, 8)) 
            c = ieor(c - a - b, ishft(b, -13)) 
            a = ieor(a - b - c, ishft(c, -12)) 
            b = ieor(b - c - a, ishft(a, 16)) 
            c = ieor(c - a - b, ishft(b, -5)) 
            a = ieor(a - b - c, ishft(c, -3)) 
            b = ieor(b - c - a, ishft(a, 10)) 
            c = ieor(c - a - b, ishft(b, -15)) 
        end subroutine b3hs_hash_key_jenkins_mix_
    end function b3hs_hash_key_jenkins 

end module hashmap_