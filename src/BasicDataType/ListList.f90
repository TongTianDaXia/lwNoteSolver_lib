module listlist_
use constants
use list_
implicit none

    private
    
    public:: listlist
    
    !--
    type::  listlist
    
        private
        type(list),dimension(:),allocatable::	ls
        
    contains
	
        !init
        generic::       init    =>  init_n,     &
                                    init_size,  &
                                    init_sizeval
        procedure::     init_n
        procedure::     init_size
        procedure::     init_sizeval
        
        !--
        procedure::     init_movealloc
        
        !--
        generic::       sval    => sval_val,   &
                                    sval_vals
        procedure::     sval_val
        procedure::     sval_vals
        
        !---
        procedure::     val => val_ij
        procedure::     val_ij
        
        !--
        generic::       ptr => ptr_i
        procedure::     ptr_i
        procedure::     list => ptr_i
        
        !--
        generic::       lsize   => nlist,   &
                                    lsize_i
        procedure::     nlist
        procedure::     lsize_i
        !--
        generic::       lmax    => element_max_i
        procedure::     element_max_i
        
        !--
        procedure::     adjust
		
        !--
        procedure::     multiInverse
        
        !--
        procedure::     prints
        
    end type listlist
    

contains

    subroutine init_n(this,lsize,nlist)
    class(listlist),intent(out)::   this
    integer(ip),intent(in)::        lsize,nlist
    integer(ip)::                   i
	
        allocate(this%ls(nlist))
        do i = 1,nlist
            call this%ls(i)%init(lsize)
        enddo
		
    end subroutine init_n
	
    !--
    subroutine init_size(this,sizes)
    class(listlist),intent(out)::			this
    integer(ip),dimension(:),intent(in)::	sizes
    integer(ip)::							i
	
        allocate(this%ls(size(sizes)))
        do i = 1,size(sizes)
            call this%ls(i)%init(sizes(i))
        enddo
		
    end subroutine init_size
	
    !--
    subroutine init_sizeVal(this,sizes,vals)
    class(listlist),intent(out)::			this
    integer(ip),dimension(:),intent(in)::	sizes,vals
    integer(ip)::							i,disp
	
        allocate(this%ls(size(sizes)))
        disp = 0
        do i = 1,size(sizes)
            call this%ls(i)%init(vals(disp + 1 : disp+sizes(i)))
            disp = disp + sizes(i)
        enddo
		
    end subroutine init_sizeVal
    
    !--
    subroutine init_movealloc(this,lsls)
    class(listlist),intent(out)::   this
    class(listlist),intent(inout):: lsls
    integer(ip)::                   i,n
	
        n = lsls%lsize()
        allocate(this%ls(n))
        do i=1,n
            call this%ls(i)%init_movealloc(lsls%ls(i))
        enddo
		
    end subroutine init_movealloc

    !--
    elemental subroutine sval_val(this, i, j, val)
    class(listlist),intent(inout):: this
    integer(ip),intent(in)::        i, j, val
	
        call this%ls(j)%sval(i, val)
		
    end subroutine sval_val
	
	!--
    pure subroutine sval_vals(this,n,vals)
    class(listlist),intent(inout)::      this
    integer(ip),intent(in)::             n
    integer(ip),dimension(:),intent(in)::vals
	
        call this%ls(n)%sval(vals)
		
    end subroutine sval_vals
    
    !--
    elemental function val_ij(this,i,j)
    class(listlist),intent(in)::    this
    integer(ip),intent(in)::        i,j
    integer(ip)::                   val_ij
	
        val_ij = this%ls(j)%val(i)
		
    end function val_ij
    
    !--
    function ptr_i(this,i) result(ls)
    class(listlist),intent(in),target::	this
    integer(ip),intent(in)::			i
    type(list),pointer::				ls
	
        ls => this%ls(i)
		
    end function ptr_i
    
    !--
    elemental function nlist(this) result(n)
    class(listlist),intent(in)::		this
    integer(ip)::						n
	
        n = size(this%ls)
		
    end function nlist
    
    !--
    elemental function lsize_i(this,i) result(n)
    class(listlist),intent(in)::        this
    integer(ip),intent(in)::            i
    integer(ip)::                       n
	
        n =  this%ls(i)%lsize()
		
    end function lsize_i
    
    !--
    pure integer(ip) function element_max_i(this) result(m)
    class(listlist),intent(in)::            this
    integer(ip)::                           i
	
        m = minisp
        do i=1,this%lsize()
            m = max(m,this%ls(i)%lmax())
        enddo
		
    end function element_max_i
    
    !--
    pure subroutine adjust(this,sizes)
    class(listlist),intent(inout)::         this
    integer(ip),dimension(:),intent(in)::   sizes
	
        call this%ls%adjust(sizes)
		
    end subroutine adjust
    
    !--
    subroutine multiInverse(this,inv_n,inv)
    class(listlist),intent(in)::    this
    integer(ip),intent(in)::        inv_n
    type(listlist),intent(out)::    inv
    integer(ip)::                   i,j,n,k
    integer(ip),dimension(inv_n)::  inv_sizes
    
        n = this%lsize()
        inv_sizes = 0
        
        do j=1,n
            do i=1,this%lsize(j)
                k = this%val(i,j)
                inv_sizes(k) = inv_sizes(k) + 1
            enddo
        enddo
        
        call inv%init(inv_sizes)
        inv_sizes = 0
        
        do j=1,n
            do i=1,this%ls(j)%lsize()
                k = this%val(i,j)
                inv_sizes(k) = inv_sizes(k) + 1
                call inv%sval(inv_sizes(k), k, j)
            enddo
        enddo
        
    end subroutine multiInverse
    
    
    !--
    subroutine prints(this)
    class(listlist),intent(in)::    this
    integer(ip)                     i
	
        print*, 'listlistprint Starts'
        print*, 'nlist:',size(this%ls)
        print*, '-------------------------------'
        do i=1,this%nlist()
            print*, '!--'
            print*, 'listnumber:',i
            call this%ls(i)%prints
        enddo
        print*, '-------------------------------'
        print*, 'listlistprint Ends'
		
    end subroutine prints
    
end module listlist_