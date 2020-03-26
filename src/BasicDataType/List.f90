module list_
use constants
implicit none

    private
    
    public:: list
    
    !--
    type::  list
        
        private
        integer(ip),allocatable,dimension(:)::  ls
        
    contains
    
        !constructor
        generic::           init    =>  init_n,     &
                                        init_n_val, &
                                        init_vals
        procedure::         init_n
        procedure::         init_n_val
        procedure::         init_vals
        
        !--
        generic::           init_movealloc  => init_movealloc_vals, &
                                                init_movealloc_list
        procedure::         init_movealloc_vals
        procedure::         init_movealloc_list
        
        !--
        generic::           sval    =>  sval_val,   &
                                        sval_nval,  &
                                        sval_vals,  &
                                        sval_sectionvals
        procedure::         sval_val
        procedure::         sval_nval
        procedure::         sval_vals
        procedure::         sval_sectionvals
        
        !--
        generic::           ptr => ptr_0,ptr_section
        procedure,private:: ptr_0
        procedure,private:: ptr_section
        
        !--
        generic::           val => val_i,val_0,val_section
        procedure::         val_i
        procedure::         val_0
        procedure::         val_section
        
        !--
        procedure::         lsize
        procedure::         lsum
        procedure::         lmax
        procedure::         adjust
        procedure::         swap => swap_list
        !--
        procedure::         prints
        
    end type list
    
!
contains

    !-------------
    subroutine init_n(this,n)
    class(list),intent(out)::               this
    integer(ip)::                           n
        allocate(this%ls(n))
        this%ls = 0
    end subroutine init_n
	
    !--
    subroutine init_n_val(this,n,val)
    class(list),intent(out)::               this
    integer(ip)::                           n,val
        allocate(this%ls(n))
        this%ls = val
    end subroutine init_n_val
	
    !--
    subroutine init_vals(this,vals)
    class(list),intent(out)::               this
    integer(ip),dimension(:),intent(in)::   vals
        allocate(this%ls,source=vals)
    end subroutine init_vals
	
	!---------------
    subroutine init_movealloc_vals(this,vals)
    class(list),intent(out)::               this
    integer(ip),dimension(:),&
    allocatable,intent(inout)::             vals
	
        call move_alloc(vals,this%ls)
		
    end subroutine init_movealloc_vals
	
    !--
    subroutine init_movealloc_list(this,ls)
    class(list),intent(out)::               this
    class(list),intent(inout)::             ls

        call move_alloc(ls%ls,this%ls)
		
    end subroutine init_movealloc_list
    
    !----------------
    elemental subroutine sval_val(this,val)
    class(list),intent(inout)::             this
    integer(ip),intent(in)::                val
	
        this%ls(:) = val
		
    end subroutine sval_val
    
	!---
    elemental subroutine sval_nval(this,n,val)
    class(list),intent(inout)::             this
    integer(ip),intent(in)::                n,val
	
        this%ls(n) = val
		
    end subroutine sval_nval
    
	!--
    pure subroutine sval_vals(this,vals)
    class(list),intent(inout)::             this
    integer(ip),dimension(:),intent(in)::   vals
        this%ls(:) = vals
    end subroutine sval_vals
    
	!--
    pure subroutine sval_sectionvals(this,s,e,vals)
    class(list),intent(inout)::             this
    integer(ip),intent(in)::                s,e
    integer(ip),dimension(:),intent(in)::   vals
        this%ls(s:e) = vals
    end subroutine sval_sectionvals
    
    !--------------
    function ptr_0(this)
    class(list),target,intent(in)::         this
    integer(ip),dimension(:),pointer::      ptr_0
        ptr_0 => this%ls
    end function ptr_0
    
    function ptr_section(this,s,e)
    class(list),target,intent(in)::         this
    integer(ip),intent(in)::                s,e
    integer(ip),dimension(:),pointer::      ptr_section
        ptr_section => this%ls(s:e)
    end function ptr_section
    
    !---
    elemental integer(ip) function val_i(this,n)
    class(list),intent(in)::                this
    integer(ip),intent(in)::                n
        val_i = this%ls(n)
    end function val_i
    
    !---
    pure function val_0(this) result(vals)
    class(list),intent(in)::                this
    integer(ip),dimension(:),allocatable::  vals
	
        allocate(vals, source=this%ls)
		
    end function val_0
    
    !---
    pure function val_section(this,s,e) result(vals)
    class(list),intent(in)::        this
    integer(ip),intent(in)::        s,e
    integer(ip),dimension(e-s+1)::  vals
        vals = this%ls(s:e)
    end function val_section
    
    !----
    elemental integer(ip) function lsize(this)
    class(list),intent(in)::        this
        lsize = size(this%ls)
    end function lsize
    
    !---
    elemental integer(ip) function lsum(this)
    class(list),intent(in)::        this
        lsum = sum(this%ls)
    end function lsum
    
    !--
    elemental integer(ip) function lmax(this)
    class(list),intent(in)::        this
        lmax = maxval(this%ls)
    end function lmax
    
    !--
    elemental subroutine adjust(this,n)
    class(list),intent(inout)::     this
    integer(ip),intent(in)::        n
    integer(ip),dimension(n)::      ls
    integer(ip)::                   cap
	
        ls = 0
        cap = min(this%lsize(), n)
        ls(1:cap) = this%ls(1:cap)
        deallocate(this%ls)
		allocate(this%ls, source = ls)
		
    end subroutine adjust
    
    !--
    pure subroutine swap_list(this,i,j)
    class(list),intent(inout)::     this
    integer(ip),intent(in)::        i,j
	
        call swap(this%ls(i), this%ls(j))
		
    end subroutine swap_list
    
    !--
    subroutine prints(this)
    class(list),intent(in)::    this
	
        print*, '--------------------------'
        print*, 'listprint Starts'
        print*, 'listSize:',size(this%ls)
        print*, 'listContents:'
        print*, this%ls
        print*, 'listprint Ends'
		
    end subroutine prints
    
end module list_